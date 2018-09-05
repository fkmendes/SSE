package biogeo;

public class InstantaneousRateMatrixOldTestDriver {

	public static void main(String[] args) {
		InstantaneousRateMatrixOld irm = new InstantaneousRateMatrixOld();
		String flatQMatrixString = "1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0";
		irm.initByName("NumberOfStates", 3, "FlatQMatrix", flatQMatrixString);
		irm.printMatrix();
		
//		int num_states = 4;
//		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
//		Q.setCell(0, 0, 1.0);
//		Q.setCell(3, 3, 0.5);
//		double q44 = Q.getCell(3, 3, 1);
//		Q.printMatrix();
//		System.out.println(q44);
	}
}
