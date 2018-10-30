package SSE;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.google.common.base.Joiner;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class HiddenInstantaneousRateMatrix extends CalculationNode {

	final public Input<Integer> NstatesInput = new Input<>("numberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> NHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> FlatQmatrixInput = new Input<>("flatQMatrix", "Array (matrix whose rows were pasted) containing the instantaneous transition rate between character states.");
	final public Input<Boolean> DisallowDoubleTransitionsInput = new Input<>("disallowDoubleTransitions", "Whether or not to set double transition parameters to zero.", Validate.REQUIRED);
	final public Input<HiddenObservedStateMapper> HiddenObservedStateMapperInput = new Input<>("hiddenObsStateMapper", "Maps hidden states onto observed states and vice-versa.", Validate.REQUIRED);
	
	private int numberOfStates; // number of observed states (never changes)
	private int numberOfHiddenStates;
	private int totalNumberOfStates;
	private boolean disallowDoubleTransitions;
	private HiddenObservedStateMapper hiddenObsStateMapper;
	private Map<Integer, List<Integer>> doubleTransitionIRMCellsMap = new HashMap<>();
	private boolean ignoreDiagonal = true; // right now, always true (we never query Qij for i=j in SSEODE)
	private boolean irmDirty = true;
	private Double[][] q;
	
	@Override
	public void initAndValidate() {
		numberOfStates = NstatesInput.get();
			
		disallowDoubleTransitions = DisallowDoubleTransitionsInput.get();
		if (disallowDoubleTransitions) {
			hiddenObsStateMapper = HiddenObservedStateMapperInput.get();
		}

		populateIRM(ignoreDiagonal, disallowDoubleTransitions);
				
		// Q.initByName("minordimension", numberOfStates);
	}
		
	public void populateIRM(boolean ignoreDiagonal, boolean disallowDoubleTransitions) {
		Double[] matrixContent = FlatQmatrixInput.get().getValues();
		numberOfHiddenStates = NHiddenStatesInput.get(); // when rjMCMC implemented, this can change at different steps
		totalNumberOfStates = numberOfStates + numberOfHiddenStates; // this will also change as a result
		q = new Double[totalNumberOfStates][totalNumberOfStates]; // and so we need to vary the size of q
		int numberOfInputElements = FlatQmatrixInput.get().getDimension();
		int q = 0;
		int diagEntry = 0;
		
		// Making sure the right number of transition parameters was provided
		if (ignoreDiagonal) {
			if (!disallowDoubleTransitions) {
				String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements, but not double transitions), but was given ";
				if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates)) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values.");
				} else if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates)) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
			
			else if (numberOfHiddenStates > 0 && disallowDoubleTransitions){
				String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements and also double transitions), but was given ";
				if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates) - (2 * (numberOfStates * numberOfHiddenStates - numberOfHiddenStates))) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values."); // subtract diagonal and then double transitions (twice b/c it's top-right and bottom-left)
				} else if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates) - (2 * (numberOfStates * numberOfHiddenStates - numberOfHiddenStates))) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
		}
				
		// grabbing cells with double transitions (populating doubleTransitionIRMCellsMap), which will be ignored when populating IRM and remain 0's
		if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
			findDoubleTransitionIRMCells(numberOfHiddenStates, totalNumberOfStates);
			// printDoubleTransitionIRMCellsMap(); // for debugging
		}
		
		List<Integer> hiddenIdxList = new ArrayList<Integer>(); // will store list of hidden idx associated to observed idx (and be rewritten in a loop)
		for (int i=0; i<totalNumberOfStates; ++i) {
			for (int j=0; j<totalNumberOfStates; ++j) {
				if (ignoreDiagonal && j == diagEntry) { 
					setCell(i, j, 0.0);
					continue;
				} // ignore diagonal element (make it zero)
				
				if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
					hiddenIdxList = doubleTransitionIRMCellsMap.get(i);
					
					if (hiddenIdxList != null) {
						if (hiddenIdxList.contains(j)) {
							setCell(i, j, 0.0);
							continue;
						}
					}
				} // ignore double transitions (make it zero)

				setCell(i, j, matrixContent[q]);
				q += 1;
			}
				 
			diagEntry += 1; // index of element we should ignore
		}
		
		irmDirty = false; // after re-population of IRM, things are clean
	}
	
	/* Later I should make sure this is made just once for each k when rjMCMC is implemented, or even when any transition rate is operated on */
	public void findDoubleTransitionIRMCells(int numberOfHiddenStates, int totalNumberOfStates) {
		for (int i=0; i<totalNumberOfStates; ++i) {
			// i-th row

			// recording double transition cells from top-right of IRM
			if (i < numberOfStates) {
				for (int j=numberOfStates; j<totalNumberOfStates; ++j) {
					// j-th col

					int h = hiddenObsStateMapper.getHiddenFromObs(i);
					
					// this is checking the diagonal of the hidden states sub-matrix (off-diagonal are 0!)
					if ((j - numberOfStates) != h) {
						// j is hidden state idx
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
						
						// System.out.println("Top-right i=" + i + " j=" + j);
					}
				}
			}
			
			// recording double transition cells from bottom-left of IRM
			if (i >= numberOfStates) {
				int h = i-numberOfStates;
				Collection<Integer> obsIdxColl = hiddenObsStateMapper.getObsCollFromHidden(h);
				
				for (int j=0; j < numberOfStates; ++j) {
					if (!obsIdxColl.contains(j)) {			
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
					}

					// System.out.println("Bottom-left i=" + i + " j=" +j);
				}
			}
		}
	}
	
	public void setCell(int from, int to, double prob) {
		// Q.setMatrixValue(from, to, prob);
        this.q[from][to] = prob;
	}

	// getters
	public int getNumStates() {
		return numberOfStates;
	}
	
	public int getNumHiddenStates() {
		return numberOfHiddenStates;
	}
	
	public int getTotalNumStates () {
		return totalNumberOfStates;
	}
	
	public double getCell(int from, int to, double rate) {
		// return Q.getMatrixValue(from, to) * rate;
		if (irmDirty) {
			populateIRM(ignoreDiagonal, disallowDoubleTransitions); // only re-populate IRM if some transition rate was operated on
		}
		return q[from][to] * rate;
	}
	
	// helper
	public void printMatrix() {
		for (int i = 0; i < totalNumberOfStates; ++i) {
			for (int j = 0; j < totalNumberOfStates; ++j) {
				// System.out.print(Double.toString(Q.getMatrixValue(i, j)) + " ");
				System.out.print(Double.toString(q[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	public void printDoubleTransitionIRMCellsMap() {
		System.out.println("IRM cells that contain double transitions and that should get 0s.");
		Iterator<Entry<Integer, List<Integer>>> iter = doubleTransitionIRMCellsMap.entrySet().iterator();
		while (iter.hasNext()) {
			Entry<Integer, List<Integer>> entry = iter.next();
			String k = Integer.toString(entry.getKey());
			String v = Joiner.on(",").join(entry.getValue());
			System.out.println("i=" + k + " j(s)=" + v);
		}
	}
	
	protected boolean requiresRecalculation() {
		irmDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateIRM(ignoreDiagonal, disallowDoubleTransitions);
		super.restore();
	}
}