package biogeo;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public class InstantaneousRateMatrix extends CalculationNode {

	final public Input<Integer> NstatesInput = new Input<>("NumberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> FlatQmatrixInput = new Input<>("FlatQMatrix", "Array (matrix whose rows were pasted) containing the instantaneous transition rate between character states.");
	
	private int numberOfStates;
	private RealParameter Q;
	
	@Override
	public void initAndValidate() {
		numberOfStates = NstatesInput.get();
		Q = FlatQmatrixInput.get();
	}
	
//	// ctor
//	public InstantaneousRateMatrix(int num_states) {
//		this.num_states = num_states;
//		mat = new double[num_states][num_states];
//	}
	
	public void setCell(int from, int to, double prob) {
		Q.setMatrixValue(from, to, prob);
//		q[from][to] = prob;
	}

	// getters
	public int getNumStates() {
		return numberOfStates; // do I need 'this'?
	}
	
	public double getCell(int from, int to, double rate) {
		return Q.getMatrixValue(from, to) * rate;
//		return q[from][to] * rate;
	}
	
	// helper
	public void printMatrix() {
		for (int i = 0; i < numberOfStates; ++i) {
			for (int j = 0; j < numberOfStates; ++j) {
				System.out.print(Double.toString(Q.getMatrixValue(i, j)) + " ");
			}
			System.out.println();
		}
	}
}
