package SSE;

import java.util.Arrays;
import java.util.regex.Pattern;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class LambdaMuAssigner extends CalculationNode {

	final public Input<Integer> totalNstatesInput = new Input<>("totalNumberOfStates", "How many states (observed + hidden) or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> nDistinctMusInput = new Input<>("nDistinctMus", "How many distinct mu values.");
	final public Input<String> musToStatesAssignerStringInput = new Input<>("musToStates", "Comma-separated integer string, one mu per state.");
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<String> lambdasToStatesAssignerStringInput = new Input<>("lambdasToStates", "Comma-separated integer string, one lambda per state.", Validate.XOR, cladoStashInput);
	final public Input<Integer> nDistinctLambdasInput = new Input<>("nDistinctLambdas", "How many distinct lambda values.", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	
	private int totalNumberOfStates;
	private int numberOfDistinctLambdas;
	private int numberOfDistinctMus;
	private String lambdaToStatesString;
	private String muToStatesString;
	private int[] lambdaAssignments;
	private int[] muAssignments;
	private Double[] lambdasContent;
	private Double[] musContent;
	private Double[] lambda;
	private Double[] mu;
	private boolean assignerDirty = true;
	Pattern comma;
		
	@Override
	public void initAndValidate() {
		comma = Pattern.compile(",");
		populateAssigner();
	}
	
	/*
	 * called by initAndValidate and getters when assigner is dirty 
	 */
	public void populateAssigner() {
		lambdaToStatesString = lambdasToStatesAssignerStringInput.get();
		lambdaAssignments = comma.splitAsStream(lambdaToStatesString).mapToInt(Integer::parseInt).toArray(); // comma-separated string comes from xml (done only once)
		muToStatesString = musToStatesAssignerStringInput.get();
		muAssignments = comma.splitAsStream(muToStatesString).mapToInt(Integer::parseInt).toArray();
//		muAssignments = comma.splitAsStream(muToStatesString).map(Double::parseDouble).toArray(Double[]::new); // if Double[]
		lambdasContent = lambdaInput.get().getValues();
		musContent = muInput.get().getValues();
		totalNumberOfStates = totalNstatesInput.get();
		numberOfDistinctLambdas = nDistinctLambdasInput.get();
		numberOfDistinctMus = nDistinctMusInput.get();
		
		// updating lambdas (no incorporate cladogenesis support for now)
		if (cladoStashInput.get() == null) {
			lambda = new Double[totalNumberOfStates];
			updateLambdasNoClado(lambdasContent);
		}
		
		// updating mus
		mu = new Double[totalNumberOfStates];
		updateMus(musContent);
	}
	
	/*
	 * called by masquerade ball when applying mask (and updating assigner)
	 */
	public void populateAssigner(int totalNStates, int nDistinctLambdas, int nDistinctMus, Double[] aLambdaContent, Double[] aMuContent, int[] aLambdaAssignment, int[] aMuAssignment) {
		totalNumberOfStates = totalNStates;
		numberOfDistinctLambdas = nDistinctLambdas;
		numberOfDistinctMus = nDistinctMus;
		lambdaAssignments = aLambdaAssignment;
		muAssignments = aMuAssignment;
		
		// updating lambdas (no incorporate cladogenesis support for now)
		if (cladoStashInput.get() == null) {
			lambda = new Double[totalNumberOfStates];
			updateLambdasNoClado(aLambdaContent);
		}
		
		// updating mus
		mu = new Double[totalNumberOfStates];
		updateMus(aMuContent);
	}
	
	public void updateLambdasNoClado(Double[] aLambdaContent) {
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int lambdaAssignmentIdx = lambdaAssignments[i];
			lambda[i] = aLambdaContent[lambdaAssignmentIdx];
		}
			
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	public void updateMus(Double[] aMuContent) {					
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int muAssignmentIdx = muAssignments[i];
			mu[i] = aMuContent[muAssignmentIdx];
		}
			
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	// getters
	// returns lambda real parameters
	public RealParameter getLambdasRealParameter() {
		return lambdaInput.get();
	}
	
	// returns only lambdas that matter (after applying mask, for example)
	public Double[] getLambdas() {
		// if lambdas were operated on, things are dirty, we need to update
		if (assignerDirty) {
			populateAssigner();
		}
		
		System.out.println("Getting lambdas from LambdaMuAssigner: " + Arrays.toString(lambda));
		
		return lambda;
	}
	
	// returns what was obtained from Input (always the same number of lambdas)
	public Double[] getLambdaContent() {
		return lambdasContent;
	}
	
	// returns mu real parameters
	public RealParameter getMusRealParameter() {
		return muInput.get();
	}
		
	// returns only mus that matter (after applying mask, for example)
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (assignerDirty) {
			populateAssigner();
		}
		
		System.out.println("Getting mus from LambdaMuAssigner: " + Arrays.toString(mu));
		
		return mu;
	}
	
	// returns what was obtained from Input (always the same number of mus)
	public Double[] getMuContent() {
		return musContent;
	}
	
	// setters
	// for CID
//	public void setLambdas(int nDistinctLambdas, int[] aLambdaToStatesArray) {
//		numberOfDistinctLambdas = nDistinctLambdas;
//		lambdaAssignments = aLambdaToStatesArray;
//		updateLambdasNoClado();
//	}
//	
//	public void setMus(int nDistinctMus, int[] aMuToStatesArray) {
//		numberOfDistinctMus = nDistinctMus;
//		muAssignments = aMuToStatesArray;
//		updateMus();
//	}
//	
//	// for non-CID
//	public void setLambdas(int[] aLambdaToStatesArray) {
//		lambdaAssignments = aLambdaToStatesArray;
//		updateLambdasNoClado();
//	}
//	
//	public void setMus(int[] aMuToStatesArray) {
//		muAssignments = aMuToStatesArray;
//		updateMus();
//	}
	
//	public void setNDistinctLambdas(int nDistinctLambdas) {
//		numberOfDistinctLambdas = nDistinctLambdas;
//	}
//	
//	public void setNDistinctMus(int nDistinctMus) {
//		numberOfDistinctMus = nDistinctMus;
//	}
//	
//	public void setTotalNumberOfStates(int aTotalNumberOfStates) {
//		totalNumberOfStates = aTotalNumberOfStates;
//		System.out.println("Just set total number of states of lambdaMuAssigner to " + totalNumberOfStates);
//	}
	
	// getters
	public int getNDistinctLambdas() {
		return numberOfDistinctLambdas;
	}
	
	public int getNDistinctMus() {
		return numberOfDistinctMus;
	}
	
	// mcmc
	protected boolean requiresRecalculation() {
		assignerDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateAssigner();
		super.restore();
	}
}
