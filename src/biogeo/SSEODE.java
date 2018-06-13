package biogeo;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

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
		
		// I haven't personally checked this, but Will seems to have
		// noticed this behavior: every new x can have negative probs
		// or probs greater than 1 coming from the ODE stepper, due to
		// rounding error; so fixing this here 
		double[] safe_x = x;
		for (int i = 0; i < num_states * 2; ++i) {
			safe_x[i] =  x[i] < 0.0 ? 0.0 : x[i];
			safe_x[i] =  x[i] > 1.0 ? 1.0 : x[i];
		}
		
		// iterating over states (ranges)
		for (int i = 0; i < num_states; ++i) {
			
			/*
			 * Step 1: get sum of lambdas to be used in A2 and A1
			 */
			double lambda_sum = 0.0;
			
			// for each event, grab respective sp rate (lambda) and keep adding	
			for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {
				double ith_lambda = entry.getValue();
				lambda_sum += ith_lambda;
			}

			/*
			 * Step 2: equation A2 
			 */

			// extinction
			dxdt[i] = mu[i];
			
			// no event (what's inside parentheses)
			double no_event_rate = mu[i] + lambda_sum;
	        for (int j = 0; j < num_states; ++j) {
	            if (i != j) {
	                no_event_rate += Q.getCell(i, j, rate);
	            }
	        }
	        dxdt[i] -= no_event_rate * safe_x[i];

	        // speciation
	        for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {	
				String[] states = entry.getKey();
				int j = Integer.parseInt(states[1]);
				int k = Integer.parseInt(states[2]);
	        	double this_lambda = entry.getValue();
	        	
	        	// if parent state (range) is the same as current (ith) state
	        	if (i == Integer.parseInt(states[0])) {
            		dxdt[i] += this_lambda * safe_x[j] * safe_x[k];
            	}
	        }
	        
//	        // for MuSSE later
//	        else {
//	        	dxdt[i] += lambda[i] * x[i] * x[i];
//	        }

			// anagenetic change
			for (int j = 0; j < num_states; ++j) {
				if (i != j) {
					dxdt[i] += Q.getCell(i, j, rate) * safe_x[j];
				}
			}
			
			/*
			 * Step 3: equation A1 
			 */

			// no event
            dxdt[i + num_states] = -no_event_rate;

            // speciation
            for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {
            	String[] states = entry.getKey();
				int j = Integer.parseInt(states[1]);
				int k = Integer.parseInt(states[2]);
            	double this_lambda = entry.getValue();
            	
            	// if parent state (range) is the same as current (ith) state
            	if (i == Integer.parseInt(states[0])) {
            		double dnj_times_ek = safe_x[j + num_states] * safe_x[k]; // D_Nj * E_k
            		double dnk_times_ej = safe_x[k + num_states] * safe_x[j]; // D_Nj * E_k
            		dxdt[i + num_states] = this_lambda * (dnj_times_ek + dnk_times_ej);
            	}
            }
            
            // anagenetic change
            for (int j = 0; j < num_states; ++j) {
            	if (i != j) {
            		dxdt[i + num_states] += Q.getCell(i, j, rate) * safe_x[j + num_states];
            	}
            }
		}
		// end of for iteration over states (ranges)
		
	}
}