package SSE;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class HiddenInstantaneousRateMatrix extends CalculationNode {

	final public Input<Integer> NstatesInput = new Input<>("numberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> NHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> FlatQmatrixInput = new Input<>("flatQMatrix", "Array (matrix whose rows were pasted) containing the instantaneous transition rate between character states.");
	final public Input<Boolean> DisallowDoubleTransitionsInput = new Input<>("disallowDoubleTransitions", "Whether or not to set double transition parameters to zero.", Validate.REQUIRED);
	
	private int numberOfStates;
	private int numberOfHiddenStates;
	private int totalNumberOfStates;
	private boolean disallowDoubleTransitions;
	private boolean ignoreDiagonal = true; // right now, always true (we never query Qij for i=j in SSEODE)
	private boolean irmDirty = true;
	private Double[][] q;
	
	@Override
	public void initAndValidate() {
		numberOfStates = NstatesInput.get();
		int numberOfInputElements = FlatQmatrixInput.get().getDimension();
		disallowDoubleTransitions = DisallowDoubleTransitionsInput.get();
		q = new Double[numberOfStates][numberOfStates];
		
		// Making sure the right number of transition parameters was provided
		if (ignoreDiagonal) {
			String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements), but was given ";
			
			if (((numberOfStates + numberOfHiddenStates) * (numberOfStates + numberOfHiddenStates) - (numberOfStates + numberOfHiddenStates)) > numberOfInputElements) {
				throw new RuntimeException(errorMsg + "too few values.");
			} else if (((numberOfStates + numberOfHiddenStates) * (numberOfStates + numberOfHiddenStates) - (numberOfStates + numberOfHiddenStates)) < numberOfInputElements) {
				throw new RuntimeException(errorMsg + "too many values.");
			}
		}

		populateIRM(ignoreDiagonal);
				
		// Q.initByName("minordimension", numberOfStates);
	}
		
	public void populateIRM(boolean ignoreDiagonal) {
		Double[] matrixContent = FlatQmatrixInput.get().getValues();
		numberOfHiddenStates = NHiddenStatesInput.get(); // when rjMCMC implemented, this can change at different steps
		totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		int q = 0;
		int diagEntry = 0;
		
		for (int i=0; i<totalNumberOfStates; ++i) {
			for (int j=0; j<totalNumberOfStates; ++j) {
				if (ignoreDiagonal && j == diagEntry) { 
					setCell(i, j, 0.0);
					continue;
				} // ignore diagonal element (make it zero)
				
				setCell(i, j, matrixContent[q]);
				q += 1;
			}
				 
			diagEntry += 1; // index of element we should ignore
		}
		
		// setting double transitions to 0 
		if (disallowDoubleTransitions) {
			
		}
		
		irmDirty = false; // after re-population of IRM, things are clean
	}
	
	public void setCell(int from, int to, double prob) {
		// Q.setMatrixValue(from, to, prob);
        this.q[from][to] = prob;
	}

	// getters
	public int getNumStates() {
		return numberOfStates;
	}
	
	public double getCell(int from, int to, double rate) {
		// return Q.getMatrixValue(from, to) * rate;
		if (irmDirty) {
			populateIRM(ignoreDiagonal); // only re-populate IRM if some transition rate was operated on
		}
		return q[from][to] * rate;
	}
	
	// helper
	public void printMatrix() {
		for (int i = 0; i < numberOfStates; ++i) {
			for (int j = 0; j < numberOfStates; ++j) {
				// System.out.print(Double.toString(Q.getMatrixValue(i, j)) + " ");
				System.out.print(Double.toString(q[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	protected boolean requiresRecalculation() {
		irmDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateIRM(ignoreDiagonal);
		super.restore();
	}
}