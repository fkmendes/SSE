package SSE;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.app.BeastMCMC;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.*;

@Description("Cladogenetic State change Speciation and Extinction (ClaSSE) model")
@Citation(value="Goldberg EE, Igic B (2012) Tempo and mode in plant breeding system evolution. Evolution 16(12):3701-3709",
year=2012, firstAuthorSurname="Goldberg", DOI="10.1111/j.1558-5646.2012.01730.x")
public class HiddenStateDependentSpeciationExtinctionProcess extends Distribution {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<TraitStash> traitStashInput = new Input<>("traitStash", "TraitStash object containing the observed character state for each species.", Validate.REQUIRED);
	final public Input<HiddenInstantaneousRateMatrix> hirmInput = new Input<>("hiddenInstantaneousRateMatrix", "HiddenInstantaneousRateMatrix object containing anagenenetic rates for both observed and hidden states.", Validate.REQUIRED);
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies at root.", Validate.REQUIRED);
	final public Input<Boolean> cladoFlagInput = new Input<>("incorporateCladogenesis", "Whether or not to incorporate cladogenetic events.", Validate.REQUIRED);

	final public Input<Boolean> useThreadsInput = new Input<>("useThreads", "calculated the distributions in parallel using threads (default false)", false);
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

	// input
	private Tree tree;
	private TraitStash traitStash;
	private HiddenInstantaneousRateMatrix q;
	private CladogeneticSpeciationRateStash cladoStash;
	private Double[] lambda;
	private Double[] mu;
	private Double[] pi; // root eq freqs
	private int numStates; // number of observed states
	private int numHiddenStates;
	private int totalNumStates;
	private double rate;
	private boolean incorporateCladogenesis;
	
	/* Original version: sliced branches into chunks, aimed at fixed step size ODE solvers */ 
	// private double dt; // time slice size (ctor populates)
	// private double rootAge;
	// private int numTimeSlices;
	
	// members used for lk computation
	// private double[][] nodePartialScaledLksPreOde;
	private double[][] nodePartialScaledLksPostOde;
	private double[] scalingConstants;
	
	// cache used for lk computation
	// private double[][] storedNodePartialScaledLksPreOde;
	private double[][] storedNodePartialScaledLksPostOde;
	private double[] storedScalingConstants;
	
	double finalLogLk;
	double finalLk;
	int rootIdx;
	
    boolean useThreads;
	int nrOfThreads;
	
	/**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths;
    protected double[] storedBranchLengths;
	
    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;
    
	@Override
	public void initAndValidate() {		
		super.initAndValidate();
		tree = treeInput.get();
		traitStash = traitStashInput.get();
		q = hirmInput.get();
		mu = muInput.get().getValues();
		pi = piInput.get().getValues();
		numStates = q.getNumStates();
		numHiddenStates = q.getNumHiddenStates();
		totalNumStates = numStates + numHiddenStates;
		rate = 1.0;
		incorporateCladogenesis = cladoFlagInput.get();
		
		if (incorporateCladogenesis) {
			cladoStash = cladoStashInput.get();
		}
		else { lambda = lambdaInput.get().getValues(); }
		
		// Original version V0: fixed-step-size ode-related stuff
		// numTimeSlices = 1;
		// rootAge = tree.getRoot().getHeight();
		// dt = rootAge / ((double) (numTimeSlices * 50));
		
		// likelihood-related
		// nodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; // tips have initialization lks, internal nodes (and root) just after merge
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][totalNumStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		scalingConstants = new double[tree.getNodeCount()]; // equivalent to diversitree's lq (but not in log-scale), these are used as denominators during the likelihood computation
		
		Arrays.fill(scalingConstants, 1.0);		
		finalLk = 0.0;
		finalLogLk = 0.0;
		
		// cache-related
		// storedNodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; 
		storedNodePartialScaledLksPostOde = new double[tree.getNodeCount()][totalNumStates*2];  
		storedScalingConstants = new double[tree.getNodeCount()]; 
		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		hasDirt = Tree.IS_FILTHY;

        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
		nrOfThreads = useThreads ? BeastMCMC.m_nThreads : 1;
		if (useThreads && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		}
		if (useThreads) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}
	}
	
/* Original constructor before interfacing with BEAST 2 */
//	public StateDependentSpeciationExtinctionProcess(TreeParser tree, double[] lambda, double[] mu, double[] pi, int numStates,
//			TraitStash traitStash, CladogeneticSpeciationRateStash clado_stash, InstantaneousRateMatrix q, double rate,
//			boolean incorporateCladogenesis) {
//		this.lambda = lambda;
//		this.mu = mu;
//		this.pi = pi;
//		this.numStates = numStates;
//		this.traitStash = traitStash;
//		this.cladoStash = clado_stash;
//		Q = q;
//		this.rate = rate;
//		this.incorporateCladogenesis = incorporateCladogenesis;
//		
//		numTimeSlices = 1;
//		double rootAge = tree.getRoot().getHeight();
//		dt = rootAge / ((double) (numTimeSlices * 50));
//		// System.out.println("Root age (height): " + Double.toString(rootAge));
//		// System.out.println("dt: " + Double.toString(dt));
//		
//		// initializing members for lk computation
//		nodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; // tips have initialization lks, internal nodes (and root) just after merge
//		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
//		
//		scalingConstants = new double[tree.getNodeCount()];
//		Arrays.fill(scalingConstants, 1.0);
//		
//		finalLk = 0.0;
//		finalLogLk = 0.0;
//	}
	
	@Override
	public double calculateLogP() {
		muInput.get().getValues(mu); // every time this method is called, we need to update mu; instead of creating a new vector and assigning that to mu, we can just copy the contents of muInput into it (note that getValues is called with an argument here)
		piInput.get().getValues(pi); // same as above

		if (!incorporateCladogenesis) { lambdaInput.get().getValues(lambda); }

		// when a biogeographical parameter changes, tree is filthy, we can use threads
		// otherwise, caching allows us to only recompute part of the tree likelihood, and thread overhead not worth it
		if (hasDirt == Tree.IS_FILTHY && useThreads) {
			computeNodeLkUsingThreads();
		} else {
			computeNodeLk(tree.getRoot(), true);
		}
		logP = finalLogLk;
		return logP;
	}

	/* Original version: fixed-step size ODE solvers could have very small (~0) chunks to integrate over, creating problems */
	// The smallest time slice that we will attempt to numerically integrate.
	// if the lattice of dt's falls closer than this to a node, that sliver of time is ignored in the integration.
	// private static final double VERY_SMALL_TIME_SLIVER = 1e-15;

    class ThreadRunner implements Runnable {

        ThreadRunner() {
        }

        @Override
		public void run() {
        	Integer nodeIdx = null;
            try {            
            	//int k = 0;
            	while (true) {
            		nodeIdx = next();
            		//k++;
            		if (nodeIdx == null) {
                        countDown.countDown();
                        //System.err.println("Done " + k);
            			return;
            		}
            		Node node = tree.getNode(nodeIdx);
            		computeNodeLk(node, false);
            		done[nodeIdx] = true;            		
            	}
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " +
                		nodeIdx == null ? "something resulting in a null" : 
                			tree.getNode(nodeIdx).getID());
                e.printStackTrace();
                System.exit(1);
            }
        }

    } // CoreRunnable

    CountDownLatch countDown;
    ExecutorService exec;
    PriorityQueue<Integer> queue;
    boolean [] done;

    /**
     * comparator used by priority queue*
     */
    class TupleComparator implements Comparator<Integer> {
        @Override
		public int compare(final Integer o1, final Integer o2) {
        	Node node1 = tree.getNode(o1);
        	Node node2 = tree.getNode(o2);
        	
        	// leafs go first
        	if (node1.isLeaf() && !node2.isLeaf()) {
        		return -1;
        	}
        	if (!node1.isLeaf() && node2.isLeaf()) {
        		return 1;
        	}
        	
        	// both nodes are leafs OR both nodes are internal
            if (node1.getHeight() < node2.getHeight()) {
                return -1;
            }
            if (node1.getHeight() == node2.getHeight()) {
                return 0;
            }
            return 1;
        }
    }
    
    Integer next() {
    	synchronized(this) {
    		Integer i = queue.peek();
    		if (i == null) {
    			return null;
    		}
    		Node node = tree.getNode(i);
    		if (node.isLeaf()) {
        		return queue.poll();    			
    		}
    		if (!done[node.getChild(0).getNr()] || !done[node.getChild(1).getNr()]) {
    			return null;
    		}
    		return queue.poll();
    	}
    }

    private void computeNodeLkUsingThreads() {
        try {
        	// set up queue
        	done = new boolean[tree.getNodeCount()];
        	queue = new PriorityQueue<>(tree.getNodeCount(), new TupleComparator());
        	for (int i = 0; i < tree.getNodeCount(); i++) {
        		queue.add(i);
        	}
        	
            countDown = new CountDownLatch(nrOfThreads);
            // kick off the threads
            for (int i = 0; i < nrOfThreads; i++) {
            	ThreadRunner coreRunnable = new ThreadRunner();
            	exec.execute(coreRunnable);
            }
            countDown.await();
        } catch (RejectedExecutionException | InterruptedException e) {
            useThreads = false;
            Log.err.println("Stop using threads: " + e.getMessage());
            computeNodeLk(tree.getRoot(), true); // in case of error, call regular non-thread computeNodeLk
        }
    }
	
	private int computeNodeLk(Node node, boolean recurse) {
		int nodeIdx = node.getNr();
		
		// cache-related stuff
        int update = (node.isDirty() | hasDirt);
        if (branchLengths[node.getNr()] != node.getLength()) {
        	update |= Tree.IS_DIRTY;
        }
        
        // may not be necessary if there is no (relaxed) clock
        branchLengths[node.getNr()] = node.getLength();

        double [] nodePartial = this.nodePartialScaledLksPostOde[nodeIdx];
        
        if (node.isLeaf()) {
        	if (update != Tree.IS_CLEAN) {
        		// nodePartialScaledLksPreOde[nodeIdx] = traitStash.getSpLks(node.getID());
        		// nodePartialScaledLksPostOde[nodeIdx] = traitStash.getSpLks(node.getID()).clone();
        		System.arraycopy(traitStash.getSpLks(node.getID()), 0, nodePartial, 0, nodePartial.length);
        	}
        	
			// System.out.println("Leaf " + node.getID() + " has node idx: " + node.getNr());
			// System.out.println("Leaf " + node.getID() + " scaled pre-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_pre_ode[nodeIdx]));
			// System.out.println("Leaf " + node.getID() + " scaled post-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_post_ode[nodeIdx]));
		}
		
		// this block initializes node_partial_lks for internal nodes
		else {
			Node left = node.getChild(0);
			Node right = node.getChild(1);

			int leftIdx = left.getNr();
			int rightIdx = right.getNr();

			// recursion over here
			if (recurse) {
				// note that computeNodeLk returns update (int)
				update |= computeNodeLk(left, recurse); // if either left child dirty
				update |= computeNodeLk(right, recurse); // or right child dirty, I (parent) am dirty
			}
			
            if (update != Tree.IS_CLEAN) {
				// System.out.println("Recurring back to internal node " + Integer.toString(node.getNr()));
				
            	// note that children's partial lks are returned already scaled by computeNodeLk (differently from original version V0)
				final double[] leftLks = nodePartialScaledLksPostOde[leftIdx]; // at this point, left_lks has not yet been updated (scaled) for its parent merging
				final double[] rightLks = nodePartialScaledLksPostOde[rightIdx];
				
				HashMap<int[], Double> eventMap = new HashMap<int[], Double>();
				Double[] speciationRates = new Double[totalNumStates];
				if (incorporateCladogenesis) {
					eventMap = cladoStash.getEventMap();
					// System.out.println("Event map inside computeNodeLk");
					// System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
				}
				
				else {
					speciationRates = lambda;
				}
				
				// merge descendant lks
				for (int i = 0; i < numStates; ++i) {
					// E's
					
					/* Original version (V0) dealt w/ scaling here at merge */
					// node_partial_unscaled_lks[nodeIdx][i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
					// nodePartialScaledLksPreOde[nodeIdx][i] = nodePartialScaledLksPostOde[nodeIdx][i]; // same for pre-ODE
					nodePartial[i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
							
					// merging with cladogenetic component
					if (incorporateCladogenesis) {
						double likeSum = 0.0;
						
						for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
							int[] states = entry.getKey();
							int j = states[1] - 1;
							int k = states[2] - 1;
							double speciationRate = entry.getValue();
							
							// if parent state (correcting for offset) is the state of the i-th node
							if ((states[0]-1) == i) {
								// System.out.println("Currently looking at state = " + Integer.toString(i + 1));
								// System.out.println("States are: " + Arrays.toString(states));
								// System.out.println("Lambda_ijk is: " + Double.toString(speciation_rate));
								
								// D's pre-merging steps
								
								/*
								 * RevBayes version of Equation A3 (with scaling)
								 */
								// double likelihoods = leftLks[numStates + j] / scalingConstants[leftIdx] * 
								// 		 rightLks[numStates + k] / scalingConstants[rightIdx];
								// likeSum += likelihoods * speciationRate;
								
								// System.out.println("Left lk array pre-scaling: " + Arrays.toString(leftLks));
								// System.out.println("Right lk array pre-scaling: " + Arrays.toString(rightLks));
								// System.out.println("Left j: " + Double.toString(leftLks[numStates + j]) + "\nRight k: " + Double.toString(rightLks[numStates + k]));
								// System.out.println("Left k: " + Double.toString(leftLks[numStates + k]) + "\nRight j: " + Double.toString(rightLks[numStates + j]));
								// System.out.println("Left scale factor: " + Double.toString(scaling_constants[leftIdx]) + "\nRight scale factor: " + Double.toString(scaling_constants[rightIdx]));
								// System.out.println("Scaled left j: " + Double.toString(leftLks[numStates + j]/scalingConstants[leftIdx]) + "\nScaled right k: " + Double.toString(right_lks[numStates + k]/scaling_constants[rightIdx]));
								// System.out.println("Scaled left k: " + Double.toString(leftLks[numStates + k]/scalingConstants[leftIdx]) + "\nScaled right j: " + Double.toString(right_lks[numStates + j]/scaling_constants[rightIdx]));
	
								/*
								 * Equation A3 (original version V0 took care of scaling here -- not anymore)
								 */
								// double dnjDmk = (leftLks[numStates + j] / scalingConstants[leftIdx]) *
								//     (rightLks[numStates + k] / scalingConstants[rightIdx]);
								// double dnkDmj = (leftLks[numStates + k] / scalingConstants[leftIdx]) *
								//     (rightLks[numStates + j] / scalingConstants[rightIdx]);
								// double dnjDmkPlusDnkDmj = dnjDmk + dnkDmj;
								double dnjDmk = (leftLks[numStates + j] ) * (rightLks[numStates + k]);
								double dnkDmj = (leftLks[numStates + k] ) *	(rightLks[numStates + j]);
								double dnjDmkPlusDnkDmj = dnjDmk + dnkDmj;
								likeSum += 0.5 * dnjDmkPlusDnkDmj * speciationRate;
								
								// System.out.println("Root age (height): " + Double.toString(rootAge));
								// System.out.println("Current scaled likelihood sum: " + Double.toString(likeSum));
							}
						}
	
						/* Original version V0 scaled here */
						// nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] / scalingConstants[leftIdx];
						// nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / scalingConstants[rightIdx];
						
						// finalizing merging for state Di
						nodePartial[numStates + i] = likeSum;

						// preOde vector was not being used at all (even in original version, it was just to match diversitree's)
						// nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
						// so preOde for internal nodes contains post-merging, but before ODE of parent kicks in
					}
					
					// merging w/o cladogenetic component
					else {		
						// D's pre-merging steps (scaling left and right children)
						
						/* Original version V0 had scaling being done here */
						// nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] /
						//	   dLeftScalingConstant; // now we update left child lk with scaled value
						// nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / 
						//     dRightScalingConstant; // now we update right child lk with scaled value
						nodePartial[numStates + i] = leftLks[numStates + i] *
								rightLks[numStates + i];
						nodePartial[numStates + i] *= speciationRates[i];
						
						// preOde vector was not being used at all (even in original version, it was just to match diversitree's)
						// keeping track of likelihoods right before ODE, at all nodes (so at internal nodes, it's post scaling and merging)
						// nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
					}
				}
							
				// System.out.println("Normalized lks left: " + Arrays.toString(nodePartialScaledLksPostOde[leftIdx]));
				// System.out.println("Normalized lks right: " + Arrays.toString(nodePartialScaledLksPostOde[rightIdx]));
				// System.out.println("(Merged) Normalized lks from internal node: " + Arrays.toString(nodePartialScaledLksPostOde[nodeIdx]));
            }
		}
		
		// numerical integration is carried out for all branches starting at this node (if this node is not the root), up to its parent
		if (!node.isRoot()) {
            if (update != Tree.IS_CLEAN) {
				// we are going from present (begin) to past (end)
				double beginAge = node.getHeight();
				double endAge = node.getParent().getHeight();
				
				numericallyIntegrateProcess(nodePartial, beginAge, endAge);

//				// System.out.println("Initial conditions: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
//				int currentDt = 0; // counter used to multiply dt
//				while ((currentDt * dt) + beginAge < endAge) {
//					double currentDtStart = (currentDt * dt) + beginAge;
//					double currentDtEnd = ((currentDt + 1) * dt) + beginAge;
//					
//					if (currentDtEnd > endAge) {
//						currentDtEnd = endAge;
//					}
//	
//					double timeslice = currentDtEnd - currentDtStart;
//					if (timeslice >= VERY_SMALL_TIME_SLIVER) {
//						numericallyIntegrateProcess(nodePartial, currentDtStart, currentDtEnd);
//					} else {
//						// DO NOTHING BECAUSE TOO LITTLE TIME HAS PAST AND nodePartialScaledLksPostOde[nodeIdx] will be unaffected
//					}
//					
//		            currentDt++;
//				}
				
				// scaling is done not at time of merging as original V0 version, but prior to returning when recurring
				double dScalingConstant = sum(nodePartial, numStates, nodePartial.length, -1, false); // -1 means don't ignore any item
				scalingConstants[nodeIdx] = dScalingConstant;
				for (int i = 0; i < numStates; i++) {
					nodePartial[numStates + i] /= dScalingConstant;
				}
				
            }
		}
		
		// if we reach root, no more numerical integration, children have already been joined above,
		// now multiply by prior, populate final_prob
		else {
			int rootIdx = nodeIdx;
			
			// for (int i = 0; i < nodePartialScaledLksPreOde.length; ++i) {
			//     System.out.println("Pre-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPreOde[i]));
			// }
			// for (int i = 0; i < nodePartialScaledLksPostOde.length; ++i) {
			//     System.out.println("Post-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPostOde[i]));
			// }
			
			double prob = 0.0;
			for (int i = 0; i < numStates; ++i) {
				prob += pi[numStates + i] * nodePartial[numStates + i];
			}
			
			boolean takeLog = true;
			finalLogLk = Math.log(prob) + sum(scalingConstants, 0, scalingConstants.length, rootIdx, takeLog);
			
			// System.out.println("Root node state = " + Double.toString(nodeIdx));
			// System.out.println("(Sum over states) Pi * lk of state = " + Double.toString(prob));
			// System.out.println("Normalizing constants = " + Arrays.toString(scalingConstants));
			// System.out.println("Lk: " + Double.toString(finalLk));
			// System.out.println("LnLk: " + Double.toString(finalLogLk));
		}
		
		return update; // this is the reason why computeNodeLk isn't void() as in the original V0 version (we need update to carry out caching)
	}
	
	private void numericallyIntegrateProcess(double[] likelihoods, double beginAge, double endAge) {
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-6, 1.0e-6);
		HiddenSSEODE ode = new HiddenSSEODE(mu, q, rate, incorporateCladogenesis);
		
		if (incorporateCladogenesis) {
			HashMap<int[], Double> eventMap = cladoStash.getEventMap();
			// System.out.println("Event map inside ODE");
			// System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
			
			ode.setEventMap(eventMap);
			dp853.integrate(ode, beginAge, likelihoods, endAge, likelihoods);
			// System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		}
		
		else {
			ode.setSpeciationRates(lambda);
			dp853.integrate(ode, beginAge, likelihoods, endAge, likelihoods);
			// System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		}
	}

	// helper
	public static double sum(double[] arr, int fromIdx, int toIdx, int idxToIgnore, boolean takeLog) {
		double result = 0.0;
		for (int i = fromIdx; i < toIdx; ++i) {
			if (i != idxToIgnore) {
				if (takeLog) {
					result += Math.log(arr[i]);
				}
				else { result += arr[i]; }
			}
		}
		// System.out.println("Lk array being summed to get normalizing constant: " + Arrays.toString(arr));
		// System.out.println("Sum is: " + Double.toString(result));
		return result;
	}
	
	public static double prod(double[] arr, int fromIdx, int toIdx, int idxToIgnore) {
		double result = 1.0;
		for (int i = fromIdx; i < toIdx; ++i) {
			if (i != idxToIgnore) {
				result *= arr[i];
			}
		}
		return result;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
    @Override
	protected boolean requiresRecalculation() {
//    	if (true) {
//			hasDirt = Tree.IS_FILTHY;
//			return true;
//    	}
        hasDirt = Tree.IS_CLEAN;
		if ((hirmInput.get().isDirtyCalculation()) ||
			(lambdaInput.get() != null && lambdaInput.get().somethingIsDirty()) || 
			(muInput.get().somethingIsDirty()) || 
			(piInput.get().somethingIsDirty()) ||
			(cladoStash != null && cladoStash.isDirtyCalculation())) {
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
        return treeInput.get().somethingIsDirty();
        // return true;
	}
    
    @Override
    public void store() {
    	System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
    	// for (int i = 0; i < nodePartialScaledLksPreOde.length; i++) {
    	//  	System.arraycopy(nodePartialScaledLksPreOde[i], 0, storedNodePartialScaledLksPreOde[i], 0, nodePartialScaledLksPreOde[i].length);    	     		
      	// }
    	
    	// before operator changes something, we store partial likelihoods of
    	// E's and D's
    	for (int i = 0; i < nodePartialScaledLksPostOde.length; i++) {
    	 	System.arraycopy(nodePartialScaledLksPostOde[i], 0, storedNodePartialScaledLksPostOde[i], 0, nodePartialScaledLksPostOde[i].length);    	     		
    	}
    	// and also store scaling constants, so that if operator changes yield worse posterior, we can go back to pre-operation
    	System.arraycopy(scalingConstants, 0, storedScalingConstants, 0, scalingConstants.length);

    	super.store();
    }
    
    @Override
    public void restore() {
    	double [] tmp = storedBranchLengths;
    	storedBranchLengths = branchLengths;
    	branchLengths = tmp;
    	
    	tmp = storedScalingConstants;
    	storedScalingConstants = scalingConstants;
    	scalingConstants = tmp;
    	
    	// double [][] tmp2 = storedNodePartialScaledLksPreOde;
    	// storedNodePartialScaledLksPreOde = nodePartialScaledLksPreOde;
    	// nodePartialScaledLksPreOde = tmp2;
    	
    	double [][] tmp2 = storedNodePartialScaledLksPostOde;
    	storedNodePartialScaledLksPostOde = nodePartialScaledLksPostOde;
    	nodePartialScaledLksPostOde = tmp2;
    	
    	super.restore();
    }
}
