package SSE;

public class InstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		String flatQMatrixString = "1.0 1.0 1.0 2.0 2.0 3.0";
		irm.initByName("NumberOfStates", 3, "FlatQMatrix", flatQMatrixString);
		irm.printMatrix();
		System.out.println(irm.getCell(2, 1, 1.0));
	}
}