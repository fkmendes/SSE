package SSE;

import java.util.HashMap;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SSEODE implements FirstOrderDifferentialEquations {

	protected int numStates; // ctor populates
	protected double rate; // ctor populates
	protected Double[] mu; // ctor arg
	protected Double[] lambda; // ctor arg
	private InstantaneousRateMatrix q; // ctor arg

	protected boolean incorporateCladogenesis; // ctor arg
	protected HashMap<int[], Double> eventMap; // setter 
	protected boolean backwardTime;
	protected boolean extinctionOnly;

	/*
	 * Constructor
	 * Speciation rates and event map are set independently so more or less general models can use this class
	 */
	public SSEODE(Double[] mu, double rate, boolean incorporateCladogenesis) {
		this.mu = mu;
		this.rate = rate;
		this.incorporateCladogenesis = incorporateCladogenesis;
	}

	public SSEODE(Double[] mu, double rate, boolean incorporateCladogenesis, boolean backwardTime, boolean extinctionOnly) {
		this(mu, rate, incorporateCladogenesis);
		this.backwardTime = backwardTime;
		this.extinctionOnly = extinctionOnly;
	}

	public SSEODE(Double[] mu, InstantaneousRateMatrix q, double rate, boolean incorporateCladogenesis) {
		this(mu, rate, incorporateCladogenesis);
		this.q = q;
		numStates = q.getNumObsStates();
	}

	public SSEODE(Double[] mu, InstantaneousRateMatrix q, double rate, boolean incorporateCladogenesis, boolean backwardTime, boolean extinctionOnly) {
        this(mu, q, rate, incorporateCladogenesis);
		this.backwardTime = backwardTime;
		this.extinctionOnly = extinctionOnly;
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
		return numStates * 2; // numStates for E and for D
	}

	// for integrator (this is where we specify the diff eqn)
	protected double getQCell(int from, int to, double rate) {
		return q.getCell(from, to, rate);
	}
	
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

		// System.out.println(Arrays.toString(safeX));

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
					double thisLambda = entry.getValue();

					if (i == (states[0]-1)) {
						// System.out.println("Matched " + Double.toString(thisLambda));
						lambda_sum += thisLambda;
					}
					// else { System.out.println("Did not match " + Double.toString(thisLambda)); }
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
	            	no_event_rate += this.getQCell(i, j, rate);
//	                no_event_rate += q.getCell(i, j, rate);
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
					dxdt[i] += this.getQCell(i, j, rate) * safeX[j];
//					dxdt[i] += q.getCell(i, j, rate) * safeX[j];
				}
			}

            if (!backwardTime) {
                dxdt[i] *= -1;
            }

			/*
			 * Step 3: equation A1 (getting D's, second half of dxdt)
			 */
			if (extinctionOnly) {
				continue;
			}

			// no event
//			System.out.println("DNi before: " + Double.toString(dxdt[i + numStates]));
            dxdt[i + numStates] = (-no_event_rate) * safeX[i + numStates];
//          System.out.println("DNi after: " + Double.toString(dxdt[i + numStates]));

            // speciation
            if (incorporateCladogenesis) {
	            for (HashMap.Entry<int[], Double> entry : eventMap.entrySet()) {
	            	int[] states = entry.getKey();
                    // Ancestor, Left, Right
                    int a = states[0]-1;
					int l = states[1]-1;
					int r = states[2]-1;
	            	double thisLambda = entry.getValue();

                    if (backwardTime) {
	            	    // if parent state (range) is the same as current (ith) state
	            	    if (i == a) {
	            		    double dnj_times_ek = safeX[l + numStates] * safeX[r]; // D_Nj * E_k
	            		    double dnk_times_ej = safeX[r + numStates] * safeX[l]; // D_Nj * E_k
	            		    dxdt[i + numStates] += thisLambda * (dnj_times_ek + dnk_times_ej);
	            	    }
                    }
                    else {
                        if (i == l) {
	            		    dxdt[i + numStates] += thisLambda * safeX[a + numStates] * safeX[r];
                        }
                        if (i == r) {
	            		    dxdt[i + numStates] += thisLambda * safeX[a + numStates] * safeX[l];
                        }
                    }

                }
            }

            else {
            	dxdt[i + numStates] += 2.0 * lambda[i] * safeX[i] * safeX[i + numStates];
//            	System.out.println("speciation: " + Double.toString(2.0 * lambda[i] * safeX[i] + safeX[i + numStates]));
//            	System.out.println("speciation (Ei): " + Double.toString(safeX[i]));
            }

            // anagenetic change
            for (int j = 0; j < numStates; ++j) {
            	if (i != j) {
                    if (backwardTime) {
            		    dxdt[i + numStates] += this.getQCell(i, j, rate) * safeX[j + numStates];
//            		    System.out.println("anagen: " + Double.toString(q.getCell(i, j, rate) * safeX[j + numStates]));
                    } else {
            		    dxdt[i + numStates] += this.getQCell(j, i, rate) * safeX[j + numStates];
                    }
            	}
            }
		}
		// end of for iteration over states (ranges)

	}
}
