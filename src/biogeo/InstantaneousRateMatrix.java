package biogeo;

public class InstantaneousRateMatrix {

	double rate;
	int num_states;
	double[][] mat;
	
	// ctor
	public InstantaneousRateMatrix(int num_states) {
		this.num_states = num_states;
		mat = new double[num_states][num_states];
	}

	// setters
	public void setRate(double rate) {
		this.rate = rate;
	}
	
	public void setCell(int from, int to, double prob) {
		mat[from][to] = prob;
	}

	// getters
	public int getNumStates() {
		return num_states; // do I need 'this'?
	}
	
	public double getCell(int from, int to) {
		return mat[from][to];
	}
	
	// other
	public void printMatrix() {
		for (int i = 0; i < num_states; ++i) {
			for (int j = 0; j < num_states; ++j) {
				System.out.print(mat[i][j] + " ");
			}
			System.out.println();
		}
	}
}
