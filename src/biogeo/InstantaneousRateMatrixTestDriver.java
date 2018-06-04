package biogeo;

public class InstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		int num_states = 4;
		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
		Q.setCell(0, 0, 1.0);
		Q.setCell(3, 3, 0.5);
		double q44 = Q.getCell(3, 3);
		Q.printMatrix();
		System.out.println(q44);
	}
}
