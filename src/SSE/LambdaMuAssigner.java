package SSE;

import java.util.regex.Pattern;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class LambdaMuAssigner extends BEASTObject {

	final public Input<Integer> TotalNstatesInput = new Input<>("totalNumberOfStates", "How many states (observed + hidden) or geographical ranges can affect speciation and extinction.");
	final public Input<String> LambdasToStatesAssignerStringInput = new Input<>("lambdasToStates", "Comma-separated integer string, one lambda per state.");
	final public Input<String> MusToStatesAssignerStringInput = new Input<>("musToStates", "Comma-separated integer string, one mu per state.");
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	
	private int totalNumberOfStates;
	private String lambdaToStatesString;
	private String muToStatesString;
	private int[] lambdaAssignments;
	private int[] muAssignments;
	private Double[] distinctLambdas;
	private Double[] distinctMus;
	private Double[] lambda;
	private Double[] mu;
		
	@Override
	public void initAndValidate() {
		totalNumberOfStates = TotalNstatesInput.get();
		Pattern comma = Pattern.compile(",");
		lambdaToStatesString = LambdasToStatesAssignerStringInput.get();
		lambdaAssignments = comma.splitAsStream(lambdaToStatesString).mapToInt(Integer::parseInt).toArray(); // convert comma-separated string into array of ints
		
		muToStatesString = MusToStatesAssignerStringInput.get();
		muAssignments = comma.splitAsStream(muToStatesString).mapToInt(Integer::parseInt).toArray();
		// muAssignments = comma.splitAsStream(muToStatesString).map(Double::parseDouble).toArray(Double[]::new); // if Double[]
		
		updateLambdasAndMus();
	}

	public void updateLambdasAndMus() {
		lambdaInput.get().getValues(distinctLambdas);
		muInput.get().getValues(distinctMus);
		
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int lambdaAssignmentIdx = lambdaAssignments[i];
			int muAssignmentIdx = muAssignments[i];
			lambda[i] = distinctLambdas[lambdaAssignmentIdx];
			mu[i] = distinctMus[muAssignmentIdx];
		}
	}
	
	// getters
	public Double[] getLambdas() {
		return lambda;
	}
	
	public Double[] getMus() {
		return mu;
	}
}
