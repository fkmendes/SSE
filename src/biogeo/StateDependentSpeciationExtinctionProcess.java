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
	private double num_time_slices;
	
	private boolean incorporate_cladogenesis;
	
	// members used for lk computation
	private double[][] node_partial_lks;
	double final_prob;
	
	public StateDependentSpeciationExtinctionProcess(TreeParser tree, double[] lambda, double[] mu, int num_states,
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
		
		num_time_slices = 500.0;
		double root_age = tree.getRoot().getHeight();
		dt = root_age / num_time_slices * 50;
		
		// initializing members for lk computation
//		System.out.println("Number of nodes is " + String.valueOf(tree.getNodeCount()));
//		System.out.println("Number of states is " + String.valueOf(num_states * 2));
		node_partial_lks = new double[tree.getNodeCount()][num_states*2];
		final_prob = 0.0;
	}

	public void computeNodeLk(Node node, int node_idx) {		
		if (node.isLeaf() == true) {
			node_partial_lks[node_idx] = trait_stash.getSpLks(node.getID());
			System.out.println("Leaf " + node.getID() + " lks are: " + Arrays.toString(node_partial_lks[node_idx]));
		}
		
		// this block initializes node_partial_lks for internal nodes
		else {
			Node left = node.getChild(0);
			Node right = node.getChild(1);

			int left_idx = left.getNr();
			int right_idx = right.getNr();

			// recursion over here
			System.out.println("Got here.");
			computeNodeLk(left, left_idx);
			computeNodeLk(right, right_idx);
			
			System.out.println("Recurring back to internal node " + Integer.toString(node.getNr()));
			
			double[] left_lks = node_partial_lks[left_idx];
			double[] right_lks = node_partial_lks[right_idx];

			double[] speciation_rates;
			if (incorporate_cladogenesis == true) {
				
			}
			
			else {
				speciation_rates = lambda;
			}
			
			// merge descendant lks
			for (int i = 0; i < num_states; ++i) {
				node_partial_lks[node_idx][i] = left_lks[i]; // filling out Es from one of children
				
				if (incorporate_cladogenesis == true) {
					
				}
				
				else {
					node_partial_lks[node_idx][num_states + i] = left_lks[num_states + i] * right_lks[num_states + i];
				}
			}
		}
		
		// numerical integration is carried out for all branches starting at this node, up to its parent
		// but if root, then no more numerical integration
		if (node.isRoot() == false) {
			// we are going from present (begin) to past (end)
			double begin_age = node.getHeight();
			double end_age = node.getParent().getHeight();
			
			System.out.println("Initial conditions: " + Arrays.toString(node_partial_lks[node_idx]));
			int current_dt = 0; // counter used to multiply dt
			while ((current_dt * dt) + begin_age < end_age) {
				double current_dt_start = (current_dt * dt) + begin_age;
				double current_dt_end = ((current_dt + 1) * dt) + begin_age;
				
				if (current_dt_end > end_age) {
					current_dt_end = end_age;
				}
				
				numericallyIntegrateProcess(node_partial_lks[node_idx], current_dt_start, current_dt_end);
	
	            current_dt++;
			}
		}
		
		// if we reach root, no more numerical integration, we just join children, multiply by prior, return log
		else {
			if (node.isRoot() == true) {
				double prob = 0.0; 
				for (int i = 0; i < num_states; ++i) {
					prob += pi[num_states + i] * node_partial_lks[node_idx][num_states + i];
				}
							
				final_prob = Math.log(prob);
			}
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
		System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		// HashMap<String[], Double> event_map = this.cladogenesis_matrix.getEventMap();
		// ode.setEventMap(event_map);
	}
}
