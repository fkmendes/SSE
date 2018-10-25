package SSE;

public class HiddenInstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		HiddenInstantaneousRateMatrix irm = new HiddenInstantaneousRateMatrix();
		// String flatQMatrixString = "0.1 0.2 0.3 0.4 0.0 0.0 1.0 1.2 1.3 0.0 1.5 0.0 2.0 2.1 0.0 2.3 0.0 0.0 2.6 3.0 3.1 3.2 0.0 3.6 4.0 0.0 0.0 0.0 4.5 4.6 0.0 5.1 0.0 0.0 5.4 5.6 0.0 0.0 6.2 6.3 6.4 6.5"; // no diagonals, with 0.0s in place of double transitions
		String flatQMatrixString = "0.1 0.2 0.3 0.4 1.0 1.2 1.3 1.5 2.0 2.1 2.3 2.6 3.0 3.1 3.2 3.6 4.0 4.5 4.6 5.1 5.4 5.6 6.2 6.3 6.4 6.5"; // no diagonals and no double transitions!
		
		String hiddenStatesString = "0,1,2,2"; // observed state 3 and 4 will both transition into hidden state 3
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
				
		irm.initByName("numberOfStates", 4, "numberOfHiddenStates", 3, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", true, "hiddenObsStateMapper", stateMapper);
		irm.printMatrix();
//		System.out.println(irm.getCell(2, 1, 1.0));
	}
}
