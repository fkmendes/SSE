package drivers;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;

public class HiddenInstantaneousRateMatrixTestDriver {

	public static void main(String[] args) {
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		// String flatQMatrixString = "0.1 0.2 0.3 0.4 0.0 0.0 1.0 1.2 1.3 0.0 1.5 0.0 2.0 2.1 0.0 2.3 0.0 0.0 2.6 3.0 3.1 3.2 0.0 3.6 4.0 0.0 0.0 0.0 4.5 4.6 0.0 5.1 0.0 0.0 5.4 5.6 0.0 0.0 6.2 6.3 6.4 6.5"; // no diagonals, with 0.0s in place of double transitions
		String flatQMatrixString = "0.1 0.2 0.3 0.4 1.0 1.2 1.3 1.5 2.0 2.1 2.3 2.6 3.0 3.1 3.2 3.6 4.0 4.5 4.6 5.1 5.4 5.6 6.2 6.3 6.4 6.5"; // no diagonals and no double transitions!
		
		String hiddenStatesString = "0,1,2,2"; // observed state 3 and 4 will both transition into hidden state 3 (but all observed states can still transition into hidden states)
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
				
		System.out.println("If number of hidden states > 0, matrix is expanded:");
		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 3, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", true, "symmetrifyAcrossDiagonal", -1, "hiddenObsStateMapper", stateMapper);
		hirm.printMatrix();
		
		System.out.println("If number of hidden states is zero, it's just plain ClaSSE:");
		String flatQMatrixString2 = "0.2 0.3 0.4 1.0 1.2 1.3 2.0 2.1 2.3 3.0 3.1 3.2";
		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 0, "flatQMatrix", flatQMatrixString2, "disallowDoubleTransitions", false, "symmetrifyAcrossDiagonal", -1, "hiddenObsStateMapper", stateMapper);
		hirm.printMatrix();
		
		HiddenInstantaneousRateMatrix hirm2 = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString3 = "0.1 0.2 0.3 1.0 1.2 1.3 1.4 2.0 2.1 2.3 3.0 3.1 3.2 3.5 4.1 4.5 5.3 5.4"; // no double transitions and no q's for observed states that are not linked to hidden states
		
		String hiddenStatesString2 = "-1,0,-1,1"; // observed state 2 and 4 will transition into hidden states 1 and 2 (here two observed states don't transition into hidden states)
		HiddenObservedStateMapper stateMapper2 = new HiddenObservedStateMapper();
		stateMapper2.initByName("hiddenStates", hiddenStatesString2);
		
		System.out.println("If number of hidden states is smaller than number of observed states, some rows and cols are zero on top of diagonal and double transitions.");
		hirm2.initByName("numberOfStates", 4, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString3, "disallowDoubleTransitions", true, "symmetrifyAcrossDiagonal", -1, "hiddenObsStateMapper", stateMapper2);
		hirm2.printMatrix();
		
		String flatQMatrixString4 = "0.1 0.2 0.3 1.0 1.2 1.3 1.5 2.0 2.1 2.3 3.0 3.1 3.2 4.1"; // no double transitions and no q's for observed states that are not linked to hidden states
		int[] aHiddenStateAssignment = new int[] { -1, 0, -1, -1 };
		stateMapper2.setHiddenStateStrings(aHiddenStateAssignment);
		hirm2.initByName("numberOfStates", 4, "numberOfHiddenStates", 1, "flatQMatrix", flatQMatrixString4, "disallowDoubleTransitions", true, "symmetrifyAcrossDiagonal", -1, "hiddenObsStateMapper", stateMapper2);
		hirm2.printMatrix();
		
//		System.out.println("Now symmetrifying diagonals for observed state 2 (hidden state 1). Note that element row=4 col=2 (bottom-left) matches that of row=2 col=4 (top-right)");
//		hirm2.initByName("numberOfStates", 4, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString3, "disallowDoubleTransitions", true, "symmetrifyAcrossDiagonal", 1, "hiddenObsStateMapper", stateMapper2);
//		hirm2.printMatrix();
		
		
	}
}
