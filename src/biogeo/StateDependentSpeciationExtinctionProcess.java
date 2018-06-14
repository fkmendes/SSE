package biogeo;

import java.util.HashMap;

public class StateDependentSpeciationExtinctionProcess {

	private double[] lambda;
	private double[] mu;
	private CladogeneticSpeciationRateStash cladogenesis_matrix;
	private InstantaneousRateMatrix Q;
	private double rate;
	
	public void NumericallyIntegrateProcess(double[] likelihoods, double begin_age, double end_age) {
		SSEODE ode = new SSEODE(this.mu, this.Q, this.rate, false);
		HashMap<String[], Double> event_map = this.cladogenesis_matrix.getEventMap();
		ode.setEventMap(event_map);
	}
}
