package drivers.SSE;

import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.nonstiff.*;

import SSE.TestODE;

public class TestODETestDriver {

	public static void main(String[] args) {
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
				
		FirstOrderDifferentialEquations ode = new TestODE();
		
		double[] y = new double[] { 1.0 }; // initial condition of y
		
		System.out.println("Initial conditions: " + y[0]);
		for (double i = 0.0; i < 2.5; i += 0.5) {
			double end = i + 0.5;
			dp853.integrate(ode, i, y, end, y);
			System.out.println("Conditions at time " + end + ": " + y[0]);
		}
	}
}
