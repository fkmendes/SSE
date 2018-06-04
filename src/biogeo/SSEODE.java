package biogeo;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SSEODE implements FirstOrderDifferentialEquations {

	private int num_states;
	private double[] mu;
	private InstantaneousRateMatrix Q;
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

	// for integrator
	public int getDimension() {
		return num_states;
	}

	// for integrator (this is where we specify the diff eqn) 
	public void computeDerivatives(double t, double[] x, double[] dxdt) {
		
	}
}
