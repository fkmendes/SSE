package biogeo;

import java.util.HashMap;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcess {

	private TreeParser tree_parser; 
	private double[] lambda;
	private double[] mu;
	private CladogeneticSpeciationRateStash cladogenesis_matrix;
	private InstantaneousRateMatrix Q;
	private double rate;
	
	public StateDependentSpeciationExtinctionProcess(TreeParser my_tree, double[] lambda, double[] mu,
			CladogeneticSpeciationRateStash cladogenesis_matrix, InstantaneousRateMatrix q, double rate) {
		super();
		this.tree_parser = my_tree;
		this.lambda = lambda;
		this.mu = mu;
		this.cladogenesis_matrix = cladogenesis_matrix;
		Q = q;
		this.rate = rate;
	}

//	public void NumericallyIntegrateProcess(double[] likelihoods, double begin_age, double end_age) {
//		SSEODE ode = new SSEODE(this.mu, this.Q, this.rate, false);
//		HashMap<String[], Double> event_map = this.cladogenesis_matrix.getEventMap();
//		ode.setEventMap(event_map);
//	}
}
