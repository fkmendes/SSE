package SSE;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class MasqueradeBall extends CalculationNode {
	final public Input<RealParameter> modelMaskInput = new Input<>("modelMask", "Series of integers that collectively determine what (sub)model to use.");
	final public Input<HiddenInstantaneousRateMatrix> hirmInput = new Input<>("hiddenInstantaneousRateMatrix", "HiddenInstantaneousRateMatrix object containing anagenenetic rates for both observed and hidden states.", Validate.REQUIRED);
	final public Input<LambdaMuAssigner> lambdaMuAssignerInput = new Input<>("lambdaMuAssigner", "LambdaMuAssigner object that assigns distinct parameters to each state.", Validate.REQUIRED);

	private Double[] mask;
	private boolean masqueradeBallDirty = true;

	private int numberOfStates; // number of observed states
	private int numberOfHiddenStatesInMask; // number of observed states
	private int totalNumberOfStates;
	
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
		hirm = hirmInput.get(); // this will be done again every time apply mask is called
		numberOfStates = hirm.getNumStates(); // numberOfStates never changes, so done only once (same below)
		hiddenToObsAssignment = new int[numberOfStates];
		obsStatesToSymmetrify = new int[numberOfStates];
		mask = new Double[numberOfStates + 1];
		
		applyMask();
	}
	
	public void applyMask() {
		
		// mask
		modelMaskInput.get().getValues(mask);
				
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
				hiddenToObsAssignment[i] = -1;
				hiddenStateIdxToIgnore.add(numberOfStates + i); // used to update both hirm and lambdaMuAssigner
			}
								
			// add hidden state for this particular observed state, but transition from and to are symmetrical
			if (maskItem == 1) {
				obsStatesToSymmetrify[i] = i;
				// hirm.symmetrifyAcrossDiagonal(obsStateToSymmetrifyIdx);
				hiddenToObsAssignment[i] = hiddenStateIdx;
				hiddenStateIdx++;
			};
								
			// add hidden state for this particular observed state, with different transition rates
			if (maskItem == 2) {
				hiddenToObsAssignment[i] = hiddenStateIdx;
				hiddenStateIdx++;				
			}
		}
										
		System.out.println("Hidden state indexes to ignore: " + Arrays.toString(hiddenStateIdxToIgnore.toArray()));
		System.out.println("Hidden to Obs map: " + Arrays.toString(hiddenToObsAssignment));
		
		numberOfHiddenStatesInMask = 0;
		for (int i=0; i<numberOfStates; ++i) {
			if (hiddenToObsAssignment[i] != -1 ) {
				numberOfHiddenStatesInMask += 1;
			}
		}
		
		// transition matrix stuff
		hirm.setHiddenObsStateAssignment(hiddenToObsAssignment); // updating HiddenObservedStateMapper inside hirm	

		updateHIRM();
		for (int obsStateIdx: obsStatesToSymmetrify) {
			hirm.symmetrifyAcrossDiagonal(obsStateIdx);
		}
		hirm.printMatrix();
		
		// lambda and mu stuff
		totalNumberOfStates = numberOfStates + numberOfHiddenStatesInMask;
		lambdaMuAssigner = lambdaMuAssignerInput.get();
		lambdaAssignmentArray = new int[totalNumberOfStates];
		muAssignmentArray = new int[totalNumberOfStates];
		updateLambdaMuAssigner();
		
		masqueradeBallDirty = false;
	}
	
	// apply mask steps
	private void updateHIRM() {
		matrixContent = hirm.getMatrixContent();
		int newMatrixContentLength = (int) (Math.pow(numberOfStates, 2) - numberOfStates + // top-left
				numberOfHiddenStatesInMask*2 + // upper-right and bottom-left
				(Math.pow(numberOfHiddenStatesInMask, 2) - numberOfHiddenStatesInMask)); // bottom-right
		Double[] newMatrixContent = new Double[newMatrixContentLength];
		realParameterToQCell = hirm.getRealParameterToQCellMap();	
		
		System.out.println("Number of observed states = " + numberOfStates);
		System.out.println("Number of hidden states = " + numberOfHiddenStatesInMask);
		System.out.println("Total number of states = " + totalNumberOfStates);
		System.out.println("Size of matrix content = " + matrixContent.length);
		System.out.println("Size of new matrix content = " + newMatrixContentLength);
		
		int j = 0;
		for (int i=0; i<matrixContent.length; ++i) {
			int[] qCell = realParameterToQCell.get(i); // for each of the RealParameters, get the cell they occupy in Q matrix
			// System.out.println("qCell = " + Arrays.toString(qCell));
			
			// now see if either row or column of that cell has a hidden state that is inactive; if active, include that (ith) RealParameter in the newMatrixContent
			if ( !(hiddenStateIdxToIgnore.contains(qCell[0]) || hiddenStateIdxToIgnore.contains(qCell[1])) ) {  
				newMatrixContent[j] = matrixContent[i];
				// System.out.println("Putting " + i + "th RealParameter of transition rates into newMatrixContent. Its qCell was " + Arrays.toString(qCell));
				j++;
			}
		}
		
		System.out.println(Arrays.toString(matrixContent));
		System.out.println(Arrays.toString(newMatrixContent));
		
		hirm.populateIRM(true, true, 0, numberOfStates, numberOfHiddenStatesInMask, newMatrixContent);
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
		Double cidIndicator = mask[mask.length-1];
		// System.out.println("CID indicator=" + cidIndicator);
		// CID
		if (cidIndicator == 1.0) {
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
	public Double[] getPis() {		
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return lambdaMuAssigner.getPis();
	}
	
	public Double[] getLambdas() {		
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return lambdaMuAssigner.getLambdas();
	}
	
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return lambdaMuAssigner.getMus();
	}
	
	public Double[][] getQs() {
		// if qs were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return hirm.getQ();
	}
	
	protected boolean requiresRecalculation() {
		masqueradeBallDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		applyMask();
		super.restore();
	}
}
