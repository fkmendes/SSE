package biogeo;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SSEODE implements FirstOrderDifferentialEquations {

	private int num_states; // ctor populates
	private double rate; // ctor populates
	private double[] mu; // ctor arg
	private double[] lambda; // ctor arg
	private InstantaneousRateMatrix Q; // ctor arg
	private boolean incorporate_cladogenesis; // ctor arg
	private HashMap<String[], Double> event_map; // setter 
	
	/*
	 * Constructor
	 * Speciation rates and event map are set independently so more or less general models can use this class
	 */
	public SSEODE(double[] mu, InstantaneousRateMatrix Q, double rate, boolean incorporate_cladogenesis) {
		this.mu = mu;
		this.Q = Q;
		this.rate = rate;
		this.incorporate_cladogenesis = incorporate_cladogenesis;
		num_states = Q.getNumStates();
		System.out.println(num_states);
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
		return num_states * 2; // num_states for E and for D
	}

	// for integrator (this is where we specify the diff eqn) 
	public void computeDerivatives(double t, double[] x, double[] dxdt) {
		
		// I haven't personally checked this, but Will seems to have
		// noticed this behavior: every new x can have negative probs
		// or probs greater than 1 coming from the ODE stepper, due to
		// rounding error; so fixing this here 
		double[] safe_x = x;
		for (int i = 0; i < num_states * 2; ++i) {
			safe_x[i] = (x[i] < 0.0 ? 0.0 : x[i]);
			safe_x[i] = (x[i] > 1.0 ? 1.0 : x[i]);
		}
		
		// System.out.println(Arrays.toString(safe_x));
		
		// iterating over states (ranges)
		for (int i = 0; i < num_states; ++i) {
			
			/*
			 * Step 1: get sum of lambdas to be used in A2 and A1
			 */
			double lambda_sum = 0.0;
			
			if (incorporate_cladogenesis == true) {
				
				// for each event, grab respective sp rate (lambda) and keep adding	
				for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {
					double ith_lambda = entry.getValue();
					lambda_sum += ith_lambda;
				}
			}
			else {
				lambda_sum = lambda[i];
			}

			/*
			 * Step 2: equation A2 (getting E's, first half of dxdt)
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
	        if (incorporate_cladogenesis == true) {
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
	        }
	        else {
	        	dxdt[i] += lambda[i] * safe_x[i] * safe_x[i];
	        }
	        
			// anagenetic change
			for (int j = 0; j < num_states; ++j) {
				if (i != j) {
					dxdt[i] += Q.getCell(i, j, rate) * safe_x[j];
				}
			}
			
			/*
			 * Step 3: equation A1 (getting D's, second half of dxdt)
			 */

			// no event
//			System.out.println("DNi before: " + Double.toString(dxdt[i + num_states]));
            dxdt[i + num_states] = (-no_event_rate) * safe_x[i + num_states];
//          System.out.println("DNi after: " + Double.toString(dxdt[i + num_states]));

            // speciation
            if (incorporate_cladogenesis == true) {
	            for (HashMap.Entry<String[], Double> entry : event_map.entrySet()) {
	            	String[] states = entry.getKey();
					int j = Integer.parseInt(states[1]);
					int k = Integer.parseInt(states[2]);
	            	double this_lambda = entry.getValue();
	            	
	            	// if parent state (range) is the same as current (ith) state
	            	if (i == Integer.parseInt(states[0])) {
	            		double dnj_times_ek = safe_x[j + num_states] * safe_x[k]; // D_Nj * E_k
	            		double dnk_times_ej = safe_x[k + num_states] * safe_x[j]; // D_Nj * E_k
	            		dxdt[i + num_states] += this_lambda * (dnj_times_ek + dnk_times_ej);
	            	}
	            }
            }
            else {
            	dxdt[i + num_states] += 2.0 * lambda[i] * safe_x[i] * safe_x[i + num_states];
//            	System.out.println("speciation: " + Double.toString(2.0 * lambda[i] * safe_x[i] + safe_x[i + num_states]));
//            	System.out.println("speciation (Ei): " + Double.toString(safe_x[i]));
            }
            
            // anagenetic change
            for (int j = 0; j < num_states; ++j) {
            	if (i != j) {
            		dxdt[i + num_states] += Q.getCell(i, j, rate) * safe_x[j + num_states];
//            		System.out.println("anagen: " + Double.toString(Q.getCell(i, j, rate) * safe_x[j + num_states]));
            	}
            }
		}
		// end of for iteration over states (ranges)
		
	}
}