package SSE;

import java.util.regex.Pattern;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class LambdaMuAssigner extends CalculationNode {

	final public Input<Integer> TotalNstatesInput = new Input<>("totalNumberOfStates", "How many states (observed + hidden) or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> NDistinctMusInput = new Input<>("nDistinctMus", "How many distinct mu values.");
	final public Input<String> MusToStatesAssignerStringInput = new Input<>("musToStates", "Comma-separated integer string, one mu per state.");
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<String> LambdasToStatesAssignerStringInput = new Input<>("lambdasToStates", "Comma-separated integer string, one lambda per state.", Validate.XOR, cladoStashInput);
	final public Input<Integer> NDistinctLambdasInput = new Input<>("nDistinctLambdas", "How many distinct lambda values.", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	
	private int totalNumberOfStates;
	private int numberOfDistinctLambdas;
	private int numberOfDistinctMus;
	private String lambdaToStatesString;
	private String muToStatesString;
	private int[] lambdaAssignments;
	private int[] muAssignments;
	private Double[] distinctLambdas;
	private Double[] distinctMus;
	private Double[] lambda;
	private Double[] mu;
	private boolean assignerDirty = true;
		
	@Override
	public void initAndValidate() {
		populateAssigner();
	}

	public void populateAssigner() {
		Pattern comma = Pattern.compile(",");
		totalNumberOfStates = TotalNstatesInput.get();
		
		// no incorporate cladogenesis
		if (cladoStashInput.get() == null) {
			lambda = new Double[totalNumberOfStates];
			numberOfDistinctLambdas = NDistinctLambdasInput.get();
			distinctLambdas = new Double[numberOfDistinctLambdas];
			lambdaToStatesString = LambdasToStatesAssignerStringInput.get();
			lambdaAssignments = comma.splitAsStream(lambdaToStatesString).mapToInt(Integer::parseInt).toArray(); // convert comma-separated string into array of ints
			updateLambdasNoClado();
		}
				
		mu = new Double[totalNumberOfStates];
		numberOfDistinctMus = NDistinctMusInput.get();
		distinctMus = new Double[numberOfDistinctMus];
		muToStatesString = MusToStatesAssignerStringInput.get();
		muAssignments = comma.splitAsStream(muToStatesString).mapToInt(Integer::parseInt).toArray();
//		muAssignments = comma.splitAsStream(muToStatesString).map(Double::parseDouble).toArray(Double[]::new); // if Double[]
		updateMus();
	}
	
	public void updateLambdasNoClado() {
		// if lambdas were operated on, things are dirty, we need to update
		lambdaInput.get().getValues(distinctLambdas); // not creating new object, just writing on it
					
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int lambdaAssignmentIdx = lambdaAssignments[i];
			lambda[i] = distinctLambdas[lambdaAssignmentIdx];
		}
			
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	public void updateMus() {
		muInput.get().getValues(distinctMus); // not creating new object, just writing on it
					
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int muAssignmentIdx = muAssignments[i];
			
			mu[i] = distinctMus[muAssignmentIdx];
		}
			
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	// getters
	public Double[] getLambdas() {
		// if lambdas were operated on, things are dirty, we need to update
		if (assignerDirty) {
			updateLambdasNoClado();
		}
		return lambda;
	}
	
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (assignerDirty) {
			updateMus();
		}
		return mu;
	}
	
	protected boolean requiresRecalculation() {
		assignerDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateAssigner();
		super.restore();
	}
}
