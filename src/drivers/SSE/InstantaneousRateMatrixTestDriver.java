package drivers.SSE;

import SSE.InstantaneousRateMatrix;

public class InstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		String flatQMatrixString = "1.0 1.0 1.0 2.0 2.0 3.0";
		irm.initByName("numberOfStates", 3, "flatQMatrix", flatQMatrixString);
		irm.printMatrix();
		System.out.println(irm.getCell(2, 1, 1.0));
	}
}
