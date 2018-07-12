package biogeo;

public class InstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("NumberOfStates", 2, "FlatQMatrix", "0.0 0.9 0.0 0.9");
	    irm.setCell(1, 1, 0.8);
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
