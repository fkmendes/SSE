package biogeo;

import java.util.HashMap;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SSEODE implements FirstOrderDifferentialEquations {

	private int num_states;
	private double[] mu;
	private double[] lambda;
	private InstantaneousRateMatrix Q;
	private HashMap<String[], Double> event_map;
	private double rate;
	
	/*
	 * Constructor
	 * Speciation rates and event map are set independently so more or less general models can use this class
	 */
	public SSEODE(double[] mu, InstantaneousRateMatrix Q, double rate) {
		this.mu = mu;
		this.Q = Q;
		this.rate = rate;
		num_states = Q.getNumStates();
	}

	// setters and getters
	public void setSpeciationRates(double[] speciation_rates) {
		lambda = speciation_rates;
	}
	
	public void setEventMap(HashMap<String[], Double> event_map) {
		this.event_map = event_map;
	}
	
	// for integrator
	public int getDimension() {
		return num_states;
	}

	// for integrator (this is where we specify the diff eqn) 
	public void computeDerivatives(double t, double[] x, double[] dxdt) {
		
		double lambda_sum = 0.0;
		
		// for each event, grab respective sp rate (lambda) and keep adding	
		for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {
			double ith_lambda = entry.getValue();
			lambda_sum += ith_lambda;
		}
	}
}
