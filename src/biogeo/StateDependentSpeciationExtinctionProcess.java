package biogeo;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.evolution.tree.*;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcess {

	// input
	private TreeParser tree; 
	private double[] lambda;
	private double[] mu;
	private double[] pi; // root eq freqs
	private int num_states;
	private TraitStash trait_stash;
	private CladogeneticSpeciationRateStash clado_stash;
	private InstantaneousRateMatrix Q;
	private double rate;
	private double dt; // time slice size (ctor populates)
	private int num_time_slices;
	
	private boolean incorporate_cladogenesis;
	
	// members used for lk computation
	private double[][] node_partial_normalized_lks_pre_ode;
	private double[][] node_partial_normalized_lks_post_ode;
	private boolean[][] node_is_normalized;
	private double[] normalizing_constants;
	private double[] normalizing_test;
	double final_log_prob;
	double final_prob;
	int root_idx;
	
	public StateDependentSpeciationExtinctionProcess(TreeParser tree, double[] lambda, double[] mu, double[] pi, int num_states,
			TraitStash trait_stash, CladogeneticSpeciationRateStash clado_stash, InstantaneousRateMatrix q, double rate,
			boolean incorporate_cladogenesis) {
		this.tree = tree;
		this.lambda = lambda;
		this.mu = mu;
		this.pi = pi;
		this.num_states = num_states;
		this.trait_stash = trait_stash;
		this.clado_stash = clado_stash;
		Q = q;
		this.rate = rate;
		this.incorporate_cladogenesis = incorporate_cladogenesis;
		
		num_time_slices = 1;
		double root_age = tree.getRoot().getHeight();
		System.out.println("Root age (height): " + Double.toString(root_age));
		dt = root_age / ((double) (num_time_slices * 50));
		System.out.println("dt: " + Double.toString(dt));
		
		// initializing members for lk computation
		node_partial_normalized_lks_pre_ode = new double[tree.getNodeCount()][num_states*2]; // tips have initialization lks, internal nodes (and root) just after merge
		node_partial_normalized_lks_post_ode = new double[tree.getNodeCount()][num_states*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		
		normalizing_test = new double[tree.getNodeCount()];
		normalizing_constants = new double[tree.getNodeCount()];
		Arrays.fill(normalizing_constants, 1.0);
		node_is_normalized = new boolean[tree.getNodeCount()][num_states*2]; // tips and internal nodes have true/false for normalization (when incorporating cladogenesis)
		
		final_prob = 0.0;
		final_log_prob = 0.0;
	}

	public double getProb() {
		return final_prob;
	}
	
	public void computeNodeLk(Node node, int node_idx) {		
		if (node.isLeaf() == true) {
			node_partial_normalized_lks_pre_ode[node_idx] = trait_stash.getSpLks(node.getID());
			node_partial_normalized_lks_post_ode[node_idx] = trait_stash.getSpLks(node.getID()).clone();
			System.out.println("Leaf " + node.getID() + " has node idx: " + node.getNr());
			System.out.println("Leaf " + node.getID() + " normalized pre-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_pre_ode[node_idx]));
			System.out.println("Leaf " + node.getID() + " normalized post-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
		}
		
		// this block initializes node_partial_lks for internal nodes
		else {
			Node left = node.getChild(0);
			Node right = node.getChild(1);

			int left_idx = left.getNr();
			int right_idx = right.getNr();

			// recursion over here
			computeNodeLk(left, left_idx);
			computeNodeLk(right, right_idx);
			
			System.out.println("Recurring back to internal node " + Integer.toString(node.getNr()));
			
			// double[] left_lks = node_partial_lks[left_idx];
			double[] left_lks = node_partial_normalized_lks_post_ode[left_idx].clone(); // at this point, left_lks has not yet been updated (scaled) for its parent merging
			
			// double[] right_lks = node_partial_lks[right_idx];
			double[] right_lks = node_partial_normalized_lks_post_ode[right_idx].clone(); // at this point (pre-merging), right_lks has not yet been scaled
			
			double D_left_normalizing_constant = sum(left_lks, num_states, left_lks.length);
			normalizing_constants[left_idx] = D_left_normalizing_constant;
			// (((H,C),G),O): does not match that of internal nodes when ClaSSE
			// hc-gorilla in diversitree is 0.2225018, but here my sum is larger = 0.22279
			// System.out.println("Left normalizing constant:" + Double.toString(D_left_normalizing_constant));
			
			double D_right_normalizing_constant = sum(right_lks, num_states, right_lks.length);
			normalizing_constants[right_idx] = D_right_normalizing_constant;
			// System.out.println("Right normalizing constant:" + Double.toString(D_right_normalizing_constant));
			
			// resetting matrix of flags 
			for (boolean[] row: node_is_normalized) {
				Arrays.fill(row, false);	
			}
//			System.out.println("Reset node_is_normalized");
//			for (boolean[] row: node_is_normalized) {
//				System.out.println(Arrays.toString(row));
//			}
			
			HashMap<int[], Double> event_map = new HashMap<int[], Double>();
			double[] speciation_rates = new double[num_states];
			if (incorporate_cladogenesis == true) {
				event_map = clado_stash.getEventMap();
			}
			
			else {
				speciation_rates = lambda;
				
//				D_left_normalizing_constant = sum(left_lks, num_states, left_lks.length);
//				normalizing_constants[left_idx] = D_left_normalizing_constant;
				// System.out.println("Left normalizing constant:" + Double.toString(D_left_normalizing_constant));
				
//				D_right_normalizing_constant = sum(right_lks, num_states, right_lks.length);
//				normalizing_constants[right_idx] = D_right_normalizing_constant;
				// System.out.println("Right normalizing constant:" + Double.toString(D_right_normalizing_constant));	
			}
			
			// merge descendant lks
			for (int i = 0; i < num_states; ++i) {
				// E's
				// node_partial_unnormalized_lks[node_idx][i] = (left_lks[i] + right_lks[i])/2; // filling out Es using avg of children
				node_partial_normalized_lks_post_ode[node_idx][i] = (left_lks[i] + right_lks[i])/2; // filling out Es using avg of children
				node_partial_normalized_lks_pre_ode[node_idx][i] = node_partial_normalized_lks_post_ode[node_idx][i]; // same for pre-ODE
						
				// merging with cladogenetic component
				if (incorporate_cladogenesis == true) {
					double like_sum = 0.0;
					
					for (HashMap.Entry<int[], Double> entry: event_map.entrySet()) {
						int[] states = entry.getKey();
						int j = states[1] - 1;
						int k = states[2] - 1;
						double speciation_rate = entry.getValue();
						
						if ((states[0]-1) == i) {
							 System.out.println("Currently looking at state = " + Integer.toString(i + 1));
							 System.out.println("States are: " + Arrays.toString(states));
							 System.out.println("Lambda_ijk is: " + Double.toString(speciation_rate));
							
							// D's pre-merging steps
							
							/*
							 * RevBayes version of Equation A3
							 */
//							double likelihoods = left_lks[num_states + j] / normalizing_constants[left_idx] * 
//									right_lks[num_states + k] / normalizing_constants[right_idx];
							// double likelihoods = left_lks[num_states + j] * right_lks[num_states + k];
//							like_sum += likelihoods * speciation_rate;
							
							System.out.println("Left lk array pre-scaling: " + Arrays.toString(left_lks));
							System.out.println("Right lk array pre-scaling: " + Arrays.toString(right_lks));
							System.out.println("Left j: " + Double.toString(left_lks[num_states + j]) + "\nRight k: " + Double.toString(right_lks[num_states + k]));
							System.out.println("Left k: " + Double.toString(left_lks[num_states + k]) + "\nRight j: " + Double.toString(right_lks[num_states + j]));
							System.out.println("Left scale factor: " + Double.toString(normalizing_constants[left_idx]) + "\nRight scale factor: " + Double.toString(normalizing_constants[right_idx]));
							System.out.println("Normalized left j: " + Double.toString(left_lks[num_states + j]/normalizing_constants[left_idx]) + "\nNormalized right k: " + Double.toString(right_lks[num_states + k]/normalizing_constants[right_idx]));
							System.out.println("Normalized left k: " + Double.toString(left_lks[num_states + k]/normalizing_constants[left_idx]) + "\nNormalized right j: " + Double.toString(right_lks[num_states + j]/normalizing_constants[right_idx]));

							/*
							 * Equation A3
							 */
							double dnj_dmk = (left_lks[num_states + j] / normalizing_constants[left_idx]) *
									(right_lks[num_states + k] / normalizing_constants[right_idx]);
							System.out.println("Dnj * Dmk " + Double.toString(dnj_dmk));
							double dnk_dmj = (left_lks[num_states + k] / normalizing_constants[left_idx]) *
									(right_lks[num_states + j] / normalizing_constants[right_idx]);
							System.out.println("Dnk * Dmj " + Double.toString(dnk_dmj));
							double dnj_dmk_plus_dnk_dmj = dnj_dmk + dnk_dmj;
							like_sum += 0.5 * dnj_dmk_plus_dnk_dmj * speciation_rate;
							System.out.println("Current normalized likelihood sum: " + Double.toString(like_sum));
						}
					}

					node_partial_normalized_lks_post_ode[left_idx][num_states + i] = left_lks[num_states + i] / normalizing_constants[left_idx];
					node_partial_normalized_lks_post_ode[right_idx][num_states + i] = right_lks[num_states + i] / normalizing_constants[right_idx];
					
					// finalizing merging for state Di
					node_partial_normalized_lks_post_ode[node_idx][num_states + i] = like_sum;
					node_partial_normalized_lks_pre_ode[node_idx][num_states + i] = node_partial_normalized_lks_post_ode[node_idx][num_states + i]; // this is a double, so no deep copy necessary
					// so pre_ode for internal nodes contains post-merging, but before ODE of 
				}
				
				// merging w/o cladogenetic component
				else {		
					// D's pre-merging steps (normalizing left and right children)
					
					// scaling
					node_partial_normalized_lks_post_ode[left_idx][num_states + i] = left_lks[num_states + i] /
							D_left_normalizing_constant; // now we update left child lk with normalized value
					node_partial_normalized_lks_post_ode[right_idx][num_states + i] = right_lks[num_states + i] / 
							D_right_normalizing_constant; // now we update right child lk with normalized value
					node_partial_normalized_lks_post_ode[node_idx][num_states + i] = node_partial_normalized_lks_post_ode[left_idx][num_states + i] *
							node_partial_normalized_lks_post_ode[right_idx][num_states + i];
					node_partial_normalized_lks_post_ode[node_idx][num_states + i] *= speciation_rates[i];
					
					// keeping track of likelihoods right before ODE, at all nodes (so at internal nodes, it's post scaling and merging)
					node_partial_normalized_lks_pre_ode[node_idx][num_states + i] = node_partial_normalized_lks_post_ode[node_idx][num_states + i]; // this is a double, so no deep copy necessary
				}
			}
						
			System.out.println("Normalized lks left: " + Arrays.toString(node_partial_normalized_lks_post_ode[left_idx]));
			System.out.println("Normalized lks right: " + Arrays.toString(node_partial_normalized_lks_post_ode[right_idx]));
			System.out.println("(Merged) Normalized lks from internal node: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
		}
		
		// numerical integration is carried out for all branches starting at this node, up to its parent
		// but if root, then no more numerical integration
		if (node.isRoot() == false) {
			// we are going from present (begin) to past (end)
			double begin_age = node.getHeight();
			double end_age = node.getParent().getHeight();
			
			System.out.println("Initial conditions: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
			int current_dt = 0; // counter used to multiply dt
			while ((current_dt * dt) + begin_age < end_age) {
				double current_dt_start = (current_dt * dt) + begin_age;
				double current_dt_end = ((current_dt + 1) * dt) + begin_age;
				
				if (current_dt_end > end_age) {
					current_dt_end = end_age;
				}
				
				numericallyIntegrateProcess(node_partial_normalized_lks_post_ode[node_idx], current_dt_start, current_dt_end);
	
	            current_dt++;
			}
		}
		
		// if we reach root, no more numerical integration, children have already been joined above,
		// now multiply by prior, populate final_prob
		else {
			int root_idx = node_idx;
			System.out.println("Root node state = " + Double.toString(node_idx));
			System.out.println("Normalizing constants for nodes = " + Arrays.toString(normalizing_constants));
			
//			for (int i = 0; i < node_partial_normalized_lks_pre_ode.length; ++i) {
//				System.out.println("Unscaled Post-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(node_partial_unnormalized_lks[i]));
//			}
			for (int i = 0; i < node_partial_normalized_lks_pre_ode.length; ++i) {
				System.out.println("Pre-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(node_partial_normalized_lks_pre_ode[i]));
			}
			for (int i = 0; i < node_partial_normalized_lks_post_ode.length; ++i) {
				System.out.println("Post-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(node_partial_normalized_lks_post_ode[i]));
			}
			
			double prob = 0.0;
//			for (int i = 0; i < num_states*2; ++i) {
			for (int i = 0; i < num_states; ++i) {
				// System.out.println("Pi for state " + i + ": " + Double.toString(pi[num_states + i]));
				// de-scaling it requires multiplying by all normalizing constants 
				prob += pi[num_states + i] * node_partial_normalized_lks_post_ode[node_idx][num_states + i];
			}
						
			System.out.println("(Sum over states) Pi * lk of state = " + Double.toString(prob));
			System.out.println("Normalizing constants = " + Arrays.toString(normalizing_constants));
			System.out.println("Normalizing constants for internal nodes = " + Arrays.toString(normalizing_test));

//			final_prob = prob; // no need to de-scale if not scaled
			final_prob = prob * prod(normalizing_constants, 0, normalizing_constants.length, root_idx); // de-scaling if scaled
			System.out.println("Lk: " + Double.toString(final_prob));
			final_log_prob = Math.log(final_prob);
			System.out.println("LnLk: " + Double.toString(final_log_prob));
		}
	}
	
	public void numericallyIntegrateProcess(double[] likelihoods, double begin_age, double end_age) {
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		SSEODE ode = new SSEODE(mu, Q, rate, incorporate_cladogenesis);
		
		if (incorporate_cladogenesis == true) {
			HashMap<int[], Double> event_map = clado_stash.getEventMap();			
			ode.setEventMap(event_map);
			dp853.integrate(ode, begin_age, likelihoods, end_age, likelihoods);
			System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		}
		
		else {
			ode.setSpeciationRates(lambda);
			dp853.integrate(ode, begin_age, likelihoods, end_age, likelihoods);
			System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		}
	}
	
	public void checkAndNormalize(double[] likelihoods, double normalizing_factor,
			int node_idx, int some_state, int num_states) {
		if (this.node_is_normalized[node_idx][num_states + some_state] == false) {
			this.node_is_normalized[node_idx][num_states + some_state] = true;
			this.node_partial_normalized_lks_post_ode[node_idx][num_states + some_state] = likelihoods[num_states + some_state] /
					normalizing_factor;
		}
//		else { System.out.println("Already normalized state: " + Integer.toString(some_state+1)); }
	}
	
	// helper
	public static double sum(double[] arr, int from_idx, int to_idx) {
		double result = 0;
		for (int i = from_idx; i < to_idx; ++i) {
			result += arr[i];
		}
		System.out.println("Lk array being summed to get normalizing constant: " + Arrays.toString(arr));
		System.out.println("Sum is: " + Double.toString(result));
		return result;
	}
	
	public static double prod(double[] arr, int from_idx, int to_idx, int idx_to_ignore) {
		double result = 1.0;
		for (int i = from_idx; i < to_idx; ++i) {
			if (i != idx_to_ignore) {
				result *= arr[i];
			}
		}
		return result;
	}
}
