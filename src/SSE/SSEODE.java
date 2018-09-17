package biogeo;

import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SSEODE implements FirstOrderDifferentialEquations {

	private int numStates; // ctor populates
	private double rate; // ctor populates
	private Double[] mu; // ctor arg
	private Double[] lambda; // ctor arg
	private InstantaneousRateMatrix q; // ctor arg
	private boolean incorporateCladogenesis; // ctor arg
	private HashMap<int[], Double> eventMap; // setter 
	
	/*
	 * Constructor
	 * Speciation rates and event map are set independently so more or less general models can use this class
	 */
	public SSEODE(Double[] mu, InstantaneousRateMatrix q, double rate, boolean incorporateCladogenesis) {
		this.mu = mu;
		this.q = q;
		this.rate = rate;
		this.incorporateCladogenesis = incorporateCladogenesis;
		numStates = q.getNumStates();
		// System.out.println("SSEODE: Self-initialized " + Integer.toString(num_states) + " states.");
	}

	// setters and getters
	public void setSpeciationRates(Double[] speciationRates) {
		lambda = speciationRates;
	}
	
	public void setEventMap(HashMap<int[], Double> eventMap) {
		this.eventMap = eventMap;
	}
	
	// for integrator
	public int getDimension() {
		return numStates * 2; // num_states for E and for D
	}

	// for integrator (this is where we specify the diff eqn) 
	public void computeDerivatives(double t, double[] x, double[] dxdt) {
		
		// I haven't personally checked this, but Will seems to have
		// noticed this behavior: every new x can have negative probs
		// or probs greater than 1 coming from the ODE stepper, due to
		// rounding error; so fixing this here 
		double[] safeX = x;
		for (int i = 0; i < numStates * 2; ++i) {
			safeX[i] = (x[i] < 0.0 ? 0.0 : x[i]);
			safeX[i] = (x[i] > 1.0 ? 1.0 : x[i]);
		}
		
		// System.out.println(Arrays.toString(safe_x));
		
		// iterating over states (ranges)
		for (int i = 0; i < numStates; ++i) {
			
			/*
			 * Step 1: get sum of lambdas to be used in A2 and A1
			 */
			double lambda_sum = 0.0;
			
			if (incorporateCladogenesis) {
				
				// for each event, grab respective sp rate (lambda) and keep adding	
				for (HashMap.Entry<int[], Double> entry : eventMap.entrySet()) {
					int[] states = entry.getKey();
					double this_lambda = entry.getValue();
					
					if (i == (states[0]-1)) {
						// System.out.println("Matched " + Double.toString(this_lambda));
						lambda_sum += this_lambda;
					}
					// else { System.out.println("Did not match " + Double.toString(this_lambda)); }
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
	        for (int j = 0; j < numStates; ++j) {
	            if (i != j) {
	                no_event_rate += q.getCell(i, j, rate);
	            }
	        }
	        dxdt[i] -= no_event_rate * safeX[i];

	        // speciation
	        if (incorporateCladogenesis) {
		        for (HashMap.Entry<int[], Double> entry : eventMap.entrySet()) {	
					int[] states = entry.getKey();
					int j = states[1]-1;
					int k = states[2]-1;
		        	double thisLambda = entry.getValue();
		        	
		        	// if parent state (range) is the same as current (ith) state
		        	if (i == (states[0]-1)) {
	            		dxdt[i] += thisLambda * safeX[j] * safeX[k];
	            	}
		        }
	        }
	        
	        else {
	        	dxdt[i] += lambda[i] * safeX[i] * safeX[i];
	        }
	        
			// anagenetic change
			for (int j = 0; j < numStates; ++j) {
				if (i != j) {
					dxdt[i] += q.getCell(i, j, rate) * safeX[j];
				}
			}
			
			/*
			 * Step 3: equation A1 (getting D's, second half of dxdt)
			 */

			// no event
//			System.out.println("DNi before: " + Double.toString(dxdt[i + num_states]));
            dxdt[i + numStates] = (-no_event_rate) * safeX[i + numStates];
//          System.out.println("DNi after: " + Double.toString(dxdt[i + num_states]));

            // speciation
            if (incorporateCladogenesis) {
	            for (HashMap.Entry<int[], Double> entry : eventMap.entrySet()) {
	            	int[] states = entry.getKey();
					int j = states[1]-1;
					int k = states[2]-1;
	            	double thisLambda = entry.getValue();
	            	
	            	// if parent state (range) is the same as current (ith) state
	            	if (i == (states[0]-1)) {
	            		double dnjTimesEk = safeX[j + numStates] * safeX[k]; // D_Nj * E_k
	            		double dnkTimesEj = safeX[k + numStates] * safeX[j]; // D_Nj * E_k
	            		dxdt[i + numStates] += thisLambda * (dnjTimesEk + dnkTimesEj);
	            	}
	            }
            }
            
            else {
            	dxdt[i + numStates] += 2.0 * lambda[i] * safeX[i] * safeX[i + numStates];
//            	System.out.println("speciation: " + Double.toString(2.0 * lambda[i] * safe_x[i] + safe_x[i + num_states]));
//            	System.out.println("speciation (Ei): " + Double.toString(safe_x[i]));
            }
            
            // anagenetic change
            for (int j = 0; j < numStates; ++j) {
            	if (i != j) {
            		dxdt[i + numStates] += q.getCell(i, j, rate) * safeX[j + numStates];
//            		System.out.println("anagen: " + Double.toString(Q.getCell(i, j, rate) * safe_x[j + num_states]));
            	}
            }
		}
		// end of for iteration over states (ranges)
		
	}
}