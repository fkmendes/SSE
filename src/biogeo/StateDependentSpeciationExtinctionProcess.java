package biogeo;

import java.util.Arrays;
import java.util.HashMap;

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
	double final_prob;
	
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
		final_prob = 0.0;
	}

	public void computeNodeLk(Node node, int node_idx) {		
		if (node.isLeaf() == true) {
			node_partial_normalized_lks_pre_ode[node_idx] = trait_stash.getSpLks(node.getID());
			node_partial_normalized_lks_post_ode[node_idx] = trait_stash.getSpLks(node.getID()).clone();
			System.out.println("Leaf " + node.getID() + " lks are: " + Arrays.toString(node_partial_normalized_lks_pre_ode[node_idx]));
			System.out.println("Leaf " + node.getID() + " normalized lks are: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
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
			double[] left_lks = node_partial_normalized_lks_post_ode[left_idx]; // at this point (pre-merging), left_lks has not been updated to its normalized value
			double D_left_normalizing_constant = sum(left_lks, num_states, left_lks.length);
			
			// double[] right_lks = node_partial_lks[right_idx];
			double[] right_lks = node_partial_normalized_lks_post_ode[right_idx]; // at this point (pre-merging), right_lks has not been updated to its normalized value
			double D_right_normalizing_constant = sum(right_lks, num_states, right_lks.length);
			
			double[] speciation_rates = new double[num_states];
			if (incorporate_cladogenesis == true) {
				
			}
			
			else {
				speciation_rates = lambda;
			}
			
			// merge descendant lks
			for (int i = 0; i < num_states; ++i) {
				node_partial_normalized_lks_post_ode[node_idx][i] = (left_lks[i] + right_lks[i])/2; // filling out Es using avg of children
				node_partial_normalized_lks_pre_ode[node_idx][i] = node_partial_normalized_lks_post_ode[node_idx][i];
				
				if (incorporate_cladogenesis == true) {
					
				}
				
				else {
					node_partial_normalized_lks_post_ode[left_idx][num_states + i] = left_lks[num_states + i] /
							D_left_normalizing_constant; // now we update left child lk with normalized value
					node_partial_normalized_lks_post_ode[right_idx][num_states + i] = right_lks[num_states + i] / 
							D_right_normalizing_constant; // now we update right child lk with normalized value
					
					node_partial_normalized_lks_post_ode[node_idx][num_states + i] = node_partial_normalized_lks_post_ode[left_idx][num_states + i] *
							node_partial_normalized_lks_post_ode[right_idx][num_states + i];
					node_partial_normalized_lks_post_ode[node_idx][num_states + i] *= speciation_rates[i];
					node_partial_normalized_lks_pre_ode[node_idx][num_states + i] = node_partial_normalized_lks_post_ode[node_idx][num_states + i];
//					node_partial_normalized_lks[node_idx][num_states + i] = left_lks[num_states + i] * right_lks[num_states + i];
//					node_partial_normalized_lks[node_idx][num_states + i] *= speciation_rates[i];
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
			
			// System.out.println("Initial conditions: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
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
		
		// if we reach root, no more numerical integration, we just join children, multiply by prior, return log
		else {
//			System.out.println("Final lks at root: "+ Arrays.toString(node_partial_lks[node_idx]));
			for (int i = 0; i < node_partial_normalized_lks_pre_ode.length; ++i) {
					System.out.println("Pre-ODE lks for state i = " + Integer.toString(i) + ": " + Arrays.toString(node_partial_normalized_lks_pre_ode[i]));
			}
			for (int i = 0; i < node_partial_normalized_lks_post_ode.length; ++i) {
				System.out.println("Post-ODE lks for state i = " + Integer.toString(i) + ": " + Arrays.toString(node_partial_normalized_lks_post_ode[i]));
			}
			
			double prob = 0.0; 
			for (int i = 0; i < num_states; ++i) {
				// System.out.println("Pi for state " + i + ": " + Double.toString(pi[num_states + i]));
				prob += pi[num_states + i] * node_partial_normalized_lks_post_ode[node_idx][num_states + i];
			}
							
			final_prob = Math.log(prob);
		}
	}

	
	
//	public double computeRootLk() {
//		Node[] nodes = tree.getNodesAsArray();
//		
//		Node root = tree.getRoot();
//		Node left = root.getChild(0);
//		Node right = root.getChild(1);
//		
//		int root_idx = root.getNr();
//		int left_idx = left.getNr();
//		int right_idx = right.getNr();
//		
//		double[] root_lks = node_partial_lks[root_idx];
//		double[] left_lks = node_partial_lks[left_idx];
//		double[] right_lks = node_partial_lks[right_idx];
//		
//		double[] speciation_rates;
//		if (incorporate_cladogenesis == true) {
//			
//		}
//		
//		else {
//			speciation_rates = lambda;
//		}
//		
//		// merge descendant lks
//		for (int i = 0; i < num_states; ++i) {
//			root_lks = left_lks;
//			
//			if (incorporate_cladogenesis == true) {
//				
//			}
//			
//			else {
//				root_lks[num_states + i] = left_lks[num_states + i] * right_lks[num_states + i];
//			}
//		}
//		
//		double prob = 0.0; 
//		for (int i = 0; i < num_states; ++i) {
//			prob += pi[num_states + i] * root_lks[num_states + i];
//		}
//		
//		return Math.log(prob);
//	}
	
	public void numericallyIntegrateProcess(double[] likelihoods, double begin_age, double end_age) {
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		SSEODE ode = new SSEODE(mu, Q, rate, incorporate_cladogenesis);
		ode.setSpeciationRates(lambda);
		dp853.integrate(ode, begin_age, likelihoods, end_age, likelihoods);
		// System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
	}
	
	// helper
	public static double sum(double[] arr, int from_idx, int to_idx) {
		double result = 0;
		for (int i = from_idx; i < to_idx; ++i) {
			result += arr[i];
		}
		return result;
	}
}
