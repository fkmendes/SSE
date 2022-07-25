import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math3.distribution.NormalDistribution;

@Description("Mosse tip likelihood using regression model of tip rate given traits and parameters beta and epsilon")
public class MosseTipLikelihood extends CalculationNode {

	final public Input<RealParameter> beta0Input = new Input<>("beta0", "beta coefficient of first trait", Input.Validate.REQUIRED);
	final public Input<RealParameter> beta1Input = new Input<>("beta1", "beta coefficient of second trait", Input.Validate.REQUIRED);
	final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "error term of regression model", Input.Validate.REQUIRED);

	private RealParameter beta0;
	private RealParameter beta1;
	private RealParameter epsilon;

	public MosseTipLikelihood() {

	}

	/**
	 * Returns P(tip rate | beta0, beta1, epsilon, trait0, trait1) ~
	 * Gaussian(mean = beta0 * trait0 + beta1 * trait1 + epsilon, sd = epsilon)
	 * within the tip rate interval (a,b)
	 * @param a start value of tip rate interval
	 * @param b end value of tip rate interval
	 * @param trait0 value of first trait
	 * @param trait1 value of second trait
	 * @return
	 */
	public double getTipLikelihood(double a, double b, double trait0, double trait1) {
		double mean = beta0.getValue() * trait0 + beta1.getValue() * trait1 + epsilon.getValue();
		double sd = epsilon.getValue();
		// Gaussian distribution
		NormalDistribution normalDist = new NormalDistribution(mean, sd);
		double startProb = normalDist.cumulativeProbability(a);
		double endProb = normalDist.cumulativeProbability(b);
		return endProb - startProb;
	}

	@Override
	public void initAndValidate() {
		beta0 = beta0Input.get();
		beta1 = beta1Input.get();
		epsilon = epsilonInput.get();
	}
}