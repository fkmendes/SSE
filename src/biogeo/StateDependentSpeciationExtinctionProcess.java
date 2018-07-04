package biogeo;

import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.evolution.tree.*;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcess {

	// input
	private double[] lambda;
	private double[] mu;
	private double[] pi; // root eq freqs
	private int numStates;
	private TraitStash traitStash;
	private CladogeneticSpeciationRateStash cladoStash;
	private InstantaneousRateMatrix Q;
	private double rate;
	private double dt; // time slice size (ctor populates)
	private int numTimeSlices;
	
	private boolean incorporateCladogenesis;
	
	// members used for lk computation
	private double[][] nodePartialScaledLksPreOde;
	private double[][] nodePartialScaledLksPostOde;
	private double[] scalingConstants;
	double finalLogLk;
	double finalLk;
	int rootIdx;
	
	public StateDependentSpeciationExtinctionProcess(TreeParser tree, double[] lambda, double[] mu, double[] pi, int numStates,
			TraitStash traitStash, CladogeneticSpeciationRateStash clado_stash, InstantaneousRateMatrix q, double rate,
			boolean incorporateCladogenesis) {
		this.lambda = lambda;
		this.mu = mu;
		this.pi = pi;
		this.numStates = numStates;
		this.traitStash = traitStash;
		this.cladoStash = clado_stash;
		Q = q;
		this.rate = rate;
		this.incorporateCladogenesis = incorporateCladogenesis;
		
		numTimeSlices = 1;
		double rootAge = tree.getRoot().getHeight();
		dt = rootAge / ((double) (numTimeSlices * 50));
		// System.out.println("Root age (height): " + Double.toString(rootAge));
		// System.out.println("dt: " + Double.toString(dt));
		
		// initializing members for lk computation
		nodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; // tips have initialization lks, internal nodes (and root) just after merge
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		
		scalingConstants = new double[tree.getNodeCount()];
		Arrays.fill(scalingConstants, 1.0);
		
		finalLk = 0.0;
		finalLogLk = 0.0;
	}

	public double getLogLk() {
		return finalLogLk;
	}
	
	public void computeNodeLk(Node node, int nodeIdx) {		
		if (node.isLeaf() == true) {
			nodePartialScaledLksPreOde[nodeIdx] = traitStash.getSpLks(node.getID());
			nodePartialScaledLksPostOde[nodeIdx] = traitStash.getSpLks(node.getID()).clone();
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
			computeNodeLk(left, leftIdx);
			computeNodeLk(right, rightIdx);
			
			// System.out.println("Recurring back to internal node " + Integer.toString(node.getNr()));
			
			double[] leftLks = nodePartialScaledLksPostOde[leftIdx].clone(); // at this point, left_lks has not yet been updated (scaled) for its parent merging
			double[] rightLks = nodePartialScaledLksPostOde[rightIdx].clone();
			
			double dLeftScalingConstant = sum(leftLks, numStates, leftLks.length);
			scalingConstants[leftIdx] = dLeftScalingConstant;
			
			double dRightScalingConstant = sum(rightLks, numStates, rightLks.length);
			scalingConstants[rightIdx] = dRightScalingConstant;
			
			HashMap<int[], Double> eventMap = new HashMap<int[], Double>();
			double[] speciationRates = new double[numStates];
			if (incorporateCladogenesis == true) {
				eventMap = cladoStash.getEventMap();
			}
			
			else {
				speciationRates = lambda;
			}
			
			// merge descendant lks
			for (int i = 0; i < numStates; ++i) {
				// E's
				// node_partial_unscaled_lks[nodeIdx][i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
				nodePartialScaledLksPostOde[nodeIdx][i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
				nodePartialScaledLksPreOde[nodeIdx][i] = nodePartialScaledLksPostOde[nodeIdx][i]; // same for pre-ODE
						
				// merging with cladogenetic component
				if (incorporateCladogenesis == true) {
					double likeSum = 0.0;
					
					for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
						int[] states = entry.getKey();
						int j = states[1] - 1;
						int k = states[2] - 1;
						double speciationRate = entry.getValue();
						
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
							 * Equation A3
							 */
							double dnjDmk = (leftLks[numStates + j] / scalingConstants[leftIdx]) *
									(rightLks[numStates + k] / scalingConstants[rightIdx]);
							double dnkDmj = (leftLks[numStates + k] / scalingConstants[leftIdx]) *
									(rightLks[numStates + j] / scalingConstants[rightIdx]);
							double dnjDmkPlusDnkDmj = dnjDmk + dnkDmj;
							likeSum += 0.5 * dnjDmkPlusDnkDmj * speciationRate;
							
							// System.out.println("Root age (height): " + Double.toString(rootAge));
							// System.out.println("Current scaled likelihood sum: " + Double.toString(likeSum));
						}
					}

					nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] / scalingConstants[leftIdx];
					nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / scalingConstants[rightIdx];
					
					// finalizing merging for state Di
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] = likeSum;
					nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
					// so pre_ode for internal nodes contains post-merging, but before ODE of 
				}
				
				// merging w/o cladogenetic component
				else {		
					// D's pre-merging steps (scaling left and right children)
					
					// scaling
					nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] /
							dLeftScalingConstant; // now we update left child lk with scaled value
					nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / 
							dRightScalingConstant; // now we update right child lk with scaled value
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[leftIdx][numStates + i] *
							nodePartialScaledLksPostOde[rightIdx][numStates + i];
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] *= speciationRates[i];
					
					// keeping track of likelihoods right before ODE, at all nodes (so at internal nodes, it's post scaling and merging)
					nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
				}
			}
						
			// System.out.println("Normalized lks left: " + Arrays.toString(nodePartialScaledLksPostOde[leftIdx]));
			// System.out.println("Normalized lks right: " + Arrays.toString(nodePartialScaledLksPostOde[rightIdx]));
			// System.out.println("(Merged) Normalized lks from internal node: " + Arrays.toString(nodePartialScaledLksPostOde[nodeIdx]));
		}
		
		// numerical integration is carried out for all branches starting at this node, up to its parent
		// but if root, then no more numerical integration
		if (node.isRoot() == false) {
			// we are going from present (begin) to past (end)
			double beginAge = node.getHeight();
			double endAge = node.getParent().getHeight();
			
			// System.out.println("Initial conditions: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
			int currentDt = 0; // counter used to multiply dt
			while ((currentDt * dt) + beginAge < endAge) {
				double currentDtStart = (currentDt * dt) + beginAge;
				double currentDtEnd = ((currentDt + 1) * dt) + beginAge;
				
				if (currentDtEnd > endAge) {
					currentDtEnd = endAge;
				}
				
				numericallyIntegrateProcess(nodePartialScaledLksPostOde[nodeIdx], currentDtStart, currentDtEnd);
	
	            currentDt++;
			}
		}
		
		// if we reach root, no more numerical integration, children have already been joined above,
		// now multiply by prior, populate final_prob
		else {
			int rootIdx = nodeIdx;
			
			for (int i = 0; i < nodePartialScaledLksPreOde.length; ++i) {
				System.out.println("Pre-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPreOde[i]));
			}
			for (int i = 0; i < nodePartialScaledLksPostOde.length; ++i) {
				System.out.println("Post-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPostOde[i]));
			}
			
			double prob = 0.0;
			for (int i = 0; i < numStates; ++i) {
				prob += pi[numStates + i] * nodePartialScaledLksPostOde[nodeIdx][numStates + i];
			}
			
			finalLk = prob * prod(scalingConstants, 0, scalingConstants.length, rootIdx); // de-scaling if scaled
			finalLogLk = Math.log(finalLk);
			
			// System.out.println("Root node state = " + Double.toString(nodeIdx));
			// System.out.println("(Sum over states) Pi * lk of state = " + Double.toString(prob));
			// System.out.println("Normalizing constants = " + Arrays.toString(scalingConstants));
			// System.out.println("Lk: " + Double.toString(finalLk));
			System.out.println("LnLk: " + Double.toString(finalLogLk));
		}
	}
	
	public void numericallyIntegrateProcess(double[] likelihoods, double beginAge, double endAge) {
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		SSEODE ode = new SSEODE(mu, Q, rate, incorporateCladogenesis);
		
		if (incorporateCladogenesis == true) {
			HashMap<int[], Double> event_map = cladoStash.getEventMap();			
			ode.setEventMap(event_map);
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
	public static double sum(double[] arr, int fromIdx, int toIdx) {
		double result = 0;
		for (int i = fromIdx; i < toIdx; ++i) {
			result += arr[i];
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
}
