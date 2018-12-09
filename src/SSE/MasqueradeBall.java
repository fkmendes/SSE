package SSE;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

public class MasqueradeBall extends CalculationNode {
	final public Input<IntegerParameter> statesMaskFlagInput = new Input<>("stateMask", "Series of integers that collectively determine what observed states have a hidden state.");
	final public Input<IntegerParameter> cidMaskFlagInput = new Input<>("cidMask", "Integers that turns CID on (1) and off(0).");
	final public Input<HiddenTraitStash> hiddenTraitStashInput = new Input<>("hiddenTraitStash", "TraitStash object containing the observed character state for each species.", Validate.REQUIRED);
	final public Input<HiddenInstantaneousRateMatrix> hirmInput = new Input<>("hiddenInstantaneousRateMatrix", "HiddenInstantaneousRateMatrix object containing anagenenetic rates for both observed and hidden states.", Validate.REQUIRED);
	final public Input<LambdaMuAssigner> lambdaMuAssignerInput = new Input<>("lambdaMuAssigner", "LambdaMuAssigner object that assigns distinct parameters to each state.", Validate.REQUIRED);

	private Integer[] statesMaskPart;
	private Integer[] cidMaskPart;
	private Integer[] mask;
	private boolean masqueradeBallDirty = true;

	private int numberOfStates; // number of observed states
	private int numberOfHiddenStatesInMask; // number of observed states
	private int totalNumberOfStates;
	
	private HiddenTraitStash hiddenTraitStash;
	
	private Double[] matrixContent;
	private int[] hiddenToObsAssignment;
	private int[] obsStatesToSymmetrify;
	private Set<Integer> hiddenStateIdxToIgnore;
	private Map<Integer, int[]> realParameterToQCell;
	private HiddenInstantaneousRateMatrix hirm;
	
	private int nDistinctLambdas;
	private int nDistinctMus;
	private int[] lambdaAssignmentArray;
	private int[] muAssignmentArray;
	private Double[] lambdaContent;
	private Double[] muContent;
	private LambdaMuAssigner lambdaMuAssigner;
	
	@Override
	public void initAndValidate() {
		hiddenTraitStash = hiddenTraitStashInput.get(); // initialized once, then updated every time apply mask is called
		hirm = hirmInput.get(); // this will be done again every time apply mask is called
		numberOfStates = hirm.getNumObsStates(); // numberOfStates never changes, so done only once (same below)
		hiddenToObsAssignment = new int[numberOfStates];
		obsStatesToSymmetrify = new int[numberOfStates];
		statesMaskPart = statesMaskFlagInput.get().getValues();
		cidMaskPart = cidMaskFlagInput.get().getValues();
		
		applyMask(statesMaskPart, cidMaskPart);
	}
	
	public void applyMask(Integer[] aStatesMaskPart, Integer[] aCIDFlag) {
		
		// mask
		mask = ArrayUtils.addAll(aStatesMaskPart, aCIDFlag);
//		System.out.println("My mask is = " + Arrays.toString(mask));
				
		/*
		 * preparing inputs to set associate hidden and observed states;
		 * < (mask.length-1) because we do not care about CID indicator at this point, which would be =(mask.length-1)
		 */
		int hiddenStateIdx = 0;
		hiddenStateIdxToIgnore = new HashSet<Integer>();
		for (int i = 0; i < mask.length-1; ++i) {
			double maskItem = mask[i];
			
			// no hidden state for this particular observed state
			if (maskItem == 0) {
				obsStatesToSymmetrify[i] = -1;
				// System.out.println("Not symmetrifying state " + i);
				hiddenToObsAssignment[i] = -1;
				hiddenStateIdxToIgnore.add(numberOfStates + i); // used to update both hirm and lambdaMuAssigner
			}
								
			// add hidden state for this particular observed state, but transition from and to are symmetrical
			if (maskItem == 1) {
				obsStatesToSymmetrify[i] = i;
				// System.out.println("Symmetrifying state " + i);
				hiddenToObsAssignment[i] = hiddenStateIdx;
				hiddenStateIdx++;
			};
								
			// add hidden state for this particular observed state, with different transition rates
			if (maskItem == 2) {
				// System.out.println("Not symmetrifying state (with mask value 2.0) " + i);
				obsStatesToSymmetrify[i] = -1;
				hiddenToObsAssignment[i] = hiddenStateIdx;
				hiddenStateIdx++;				
			}
		}
										
		// System.out.println("Hidden state indexes to ignore: " + Arrays.toString(hiddenStateIdxToIgnore.toArray()));
		// System.out.println("Hidden to Obs map: " + Arrays.toString(hiddenToObsAssignment));
		
		numberOfHiddenStatesInMask = 0;
		for (int i=0; i<numberOfStates; ++i) {
			if (hiddenToObsAssignment[i] != -1 ) {
				numberOfHiddenStatesInMask += 1;
			}
		}
		totalNumberOfStates = numberOfStates + numberOfHiddenStatesInMask;
		
		// transition matrix stuff
		hirm.setHiddenObsStateAssignment(hiddenToObsAssignment); // updating HiddenObservedStateMapper inside hirm	

		updateHIRM();
		for (int obsStateIdx: obsStatesToSymmetrify) {
			if (obsStateIdx != -1) {
				hirm.symmetrifyAcrossDiagonal(obsStateIdx);
			}
		}
		// hirm.printMatrix();
		
		// lambda and mu stuff
		lambdaMuAssigner = lambdaMuAssignerInput.get();
		lambdaAssignmentArray = new int[totalNumberOfStates];
		muAssignmentArray = new int[totalNumberOfStates];
		updateLambdaMuAssigner();
		
		// trait stash stuff
		hiddenTraitStash.setHiddenStateAssignment(hiddenToObsAssignment);
		updateTraitStash();
		
		masqueradeBallDirty = false;
	} // end apply mask
	
	// apply mask steps
	private void updateTraitStash() {
		hiddenTraitStash.populateSpLksMap(totalNumberOfStates, numberOfStates, numberOfHiddenStatesInMask);
		// System.out.println("New trait stash below");
		hiddenTraitStash.printLksMap();
	}
	
	private void updateHIRM() {
//		matrixContent = hirm.getMatrixContent();
		matrixContent = hirm.getQRealParameter().getValues();
		int newMatrixContentLength = (int) (Math.pow(numberOfStates, 2) - numberOfStates + // top-left
				numberOfHiddenStatesInMask*2 + // upper-right and bottom-left
				(Math.pow(numberOfHiddenStatesInMask, 2) - numberOfHiddenStatesInMask)); // bottom-right
		Double[] newMatrixContent = new Double[newMatrixContentLength];
		realParameterToQCell = hirm.getRealParameterToQCellMap();	
		
//		System.out.println("Number of observed states = " + numberOfStates);
//		System.out.println("Number of hidden states = " + numberOfHiddenStatesInMask);
//		System.out.println("Total number of states = " + totalNumberOfStates);
//		System.out.println("Size of matrix content = " + matrixContent.length);
//		System.out.println("Size of new matrix content = " + newMatrixContentLength);
		
		int j = 0;
		for (int i=0; i<matrixContent.length; ++i) {
			int[] qCell = realParameterToQCell.get(i); // for each of the RealParameters, get the cell they occupy in Q matrix
			
			// now see if either row or column of that cell has a hidden state that is inactive; if active, include that (ith) RealParameter in the newMatrixContent
			if ( !(hiddenStateIdxToIgnore.contains(qCell[0]) || hiddenStateIdxToIgnore.contains(qCell[1])) ) {  
				newMatrixContent[j] = matrixContent[i];
				j++;
			}
		}
		
		// System.out.println("matrixContent: " + Arrays.toString(matrixContent));
		// System.out.println("newMatrixContent: " + Arrays.toString(newMatrixContent));
		
		hirm.populateIRM(true, true, -1, numberOfStates, numberOfHiddenStatesInMask, newMatrixContent);
	}
	
	private void updateLambdaMuAssigner() {
		lambdaContent = lambdaMuAssigner.getLambdaContent();
		muContent = lambdaMuAssigner.getMuContent();
		Double[] newLambdaContent = new Double[totalNumberOfStates];
		Double[] newMuContent = new Double[totalNumberOfStates];
		
		int j = 0;
		// there are up to a total of 2 * number of observed states (i.e., all observed states have hidden states)
		for (int i=0; i<2*numberOfStates; ++i) {
			if (!hiddenStateIdxToIgnore.contains(i)) {
				newLambdaContent[j] = lambdaContent[i];
				newMuContent[j] = muContent[i];
				++j;
			}
		}
		
		// dealing with lambda and mus
		Integer cidIndicator = mask[mask.length-1];
		// System.out.println("CID indicator=" + cidIndicator);
		// CID
		if (cidIndicator == 1) {
			for (int i=0; i<totalNumberOfStates; ++i) {
				if (i < numberOfStates) {
					lambdaAssignmentArray[i] = 0;
					muAssignmentArray[i] = 0;
				}
				else {
//					lambdaAssignmentArray[i] = 1;
//					muAssignmentArray[i] = 1;
					lambdaAssignmentArray[i] = numberOfStates; // we give the first hidden state lambda to all hidden states
					muAssignmentArray[i] = numberOfStates;
				}
			}

			nDistinctLambdas = 2;
			nDistinctMus = nDistinctLambdas;
		}
		// not CID
		else {
			for (int i=0; i<totalNumberOfStates; ++i) {
				lambdaAssignmentArray[i] = i;
				muAssignmentArray[i] = i;
			}
			
			nDistinctLambdas = numberOfStates + numberOfHiddenStatesInMask;
			nDistinctMus = nDistinctLambdas;
		}

		 lambdaMuAssigner.populateAssigner(totalNumberOfStates, nDistinctLambdas, nDistinctMus, newLambdaContent, newMuContent, lambdaAssignmentArray, muAssignmentArray);
	}
	
	// getters
	public int getNumStates() {
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return totalNumberOfStates;
	}
	
	public int getNumHiddenStatesInMask() {
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return numberOfHiddenStatesInMask;
	}
	
	public Double[] getPis() {		
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return lambdaMuAssigner.getPis();
	}
	
	public Double[] getLambdas() {		
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return lambdaMuAssigner.getLambdas();
	}
	
	public RealParameter getLambdasRealParameter() {
		return lambdaMuAssigner.getLambdasRealParameter();
	}
	
	public CladogeneticSpeciationRateStash getCladoStash() {
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return lambdaMuAssigner.getCladoStash();
	}
	
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return lambdaMuAssigner.getMus();
	}
	
	public RealParameter getMusRealParameter() {
		return lambdaMuAssigner.getMusRealParameter();
	}
	
	public RealParameter getPisRealParameter() {
		return lambdaMuAssigner.getPisRealParameter();
	}
	
	public HiddenInstantaneousRateMatrix getHIRM() {
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return hirm;
	}
	
	public RealParameter getQRealParameter() {
		return hirm.getQRealParameter();
	}
	
	public HiddenTraitStash getHTS() {
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return hiddenTraitStash;
	}
	
	// for testing (actual likelihood calls getHIRM)
	public Double[][] getQs() {
		// if qs were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			statesMaskFlagInput.get().getValues(statesMaskPart);
			cidMaskFlagInput.get().getValues();
			applyMask(statesMaskPart, cidMaskPart);
		}
		return hirm.getQ();
	}
	
	// setters
	// for testing
	public void setMask(Integer[] aStatesMaskPart, Integer[] aCIDMaskPart) {
		applyMask(aStatesMaskPart, aCIDMaskPart);
	}
	
	protected boolean requiresRecalculation() {
		masqueradeBallDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		statesMaskFlagInput.get().getValues(statesMaskPart);
		cidMaskPart = cidMaskFlagInput.get().getValues();
		applyMask(statesMaskPart, cidMaskPart);
		super.restore();
	}
}
