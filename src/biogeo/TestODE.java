package biogeo;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class TestODE implements FirstOrderDifferentialEquations {

	private double C;
	
	public TestODE() {
		this.C = 1;
	}
	
	public int getDimension() {
		return 1;
	}
	
	// y is is an array containing our variables (will be updated with the values of yDot with every step)
	// yDot is an array containing the derivatives for our variables 
	public void computeDerivatives(double t, double[] y, double[] yDot) {
		yDot[0] = this.C * -y[0]; // should be equivalent to the sol'n for this diff eq: C + exp(-t)
	}
}