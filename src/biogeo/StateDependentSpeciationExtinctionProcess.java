package biogeo;

import java.util.HashMap;
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
		double[][] node_partial_lks = new double[tree.getNodeCount()][num_states*2];
	}

	public void computeNodeLk(Node node, int node_idx) {		
		if (node.isLeaf() == true) {
			node_partial_lks[node_idx] = trait_stash.getSpLks(node.getID());
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
		
		// defining starting and ending points for numerical integrator
		// we are going from present (begin) to past (end)
		double begin_age = node.getHeight();
		double end_age = node.getParent().getHeight();
		
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
	
//	public void NumericallyIntegrateProcess(double[] likelihoods, double begin_age, double end_age) {
//		SSEODE ode = new SSEODE(this.mu, this.Q, this.rate, false);
//		HashMap<String[], Double> event_map = this.cladogenesis_matrix.getEventMap();
//		ode.setEventMap(event_map);
//	}
}
