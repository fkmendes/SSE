package biogeo;

import java.util.Arrays;
import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

public class SSEODETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		String[] sp_names = new String[] { "Human" };
		int num_states = 2; // BiSSE
		double[] lambda = new double[] {0.222222222, 0.222222222};
		double[] mu = new double[] {0.0, 0.0}; // pure birth
		
		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
		Q.setCell(0, 0, 1.0);
		Q.setCell(1, 1, 1.0);
		Q.setCell(0, 1, 0.0);
		Q.setCell(1, 0, 0.0);
		
		/* testing other values */
		// double[] mu = new double[] {0.1, 0.1}; // testing other values
		// Q.setCell(0, 0, 0.9); 
		// Q.setCell(1, 1, 0.9);
		// Q.setCell(0, 1, 0.1);
		// Q.setCell(1, 0, 0.1);
		
		// Q.printMatrix(); // checking Q
		StateStash state_stash = new StateStash(num_states, sp_names);
		state_stash.setSpStateValue("Human", "1");
		state_stash.printLksMap();
		double[] y = state_stash.getSpLk("Human"); // e0t = 0.0, e1t = 0.0, d0t = 0.0, d1t = 1.0
		// double[] y = new double[] { 0.0, 0.0, 0.0, 1.0 };
		
		// initializing integrator and ODE
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		
		SSEODE ode = new SSEODE(mu, Q, 1.0, false); // incorporate_cladogenesis = false (BiSSE, or MuSSE)
		ode.setSpeciationRates(lambda);
		
		// running
		System.out.println("Initial conditions: " + Arrays.toString(y));
		for (double i = 0.0; i < 1.0; i += 0.02) {
			double end = i + 0.02;
			dp853.integrate(ode, i, y, end, y);
			System.out.println("Conditions at time " + end + ": " + Arrays.toString(y));
		}
	}

}
