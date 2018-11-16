package SSE;

import java.util.Arrays;

import beast.core.parameter.RealParameter;

public class MasqueradeBallTestDriver {

	public static void main(String[] args) {
		String hiddenStatesString = "0,1,2,3";
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
				
		String lambdasToStatesString = "0,1,2,3,4,5,6,7";
		Double lambda1 = 0.1; // 0A
		Double lambda2 = 0.15; // 1A
		Double lambda3 = 0.2; // 2A
		Double lambda4 = 0.1; // 3A
		Double lambda5 = 0.2; // 0B
		Double lambda6 = 0.3; // 1B
		Double lambda7 = 0.4; // 2B
		Double lambda8 = 0.2; // 3B
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8 };
		RealParameter lambda = new RealParameter(lambdas);
		
		String musToStatesString = "0,1,2,3,4,5,6,7";
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double mu5 = 0.06;
		Double mu6 = 0.09;
		Double mu7 = 0.12;
		Double mu8 = 0.06;
		Double[] mus = { mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8 };
		RealParameter mu = new RealParameter(mus);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 8, "nDistinctLambdas", 8, "nDistinctMus", 8, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu);
		System.out.println("Lambdas: " + Arrays.toString(lambdaMuAssigner.getLambdas()));
		System.out.println("Mus: " + Arrays.toString(lambdaMuAssigner.getMus()));
		
		boolean disallowDoubleTransitions = true; // not used
		int symmetrifyAcrossDiagonal = -1;
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.1 0.2 0.3 0.4 1.0 1.2 1.3 1.5 2.0 2.1 2.3 2.6 3.0 3.1 3.2 3.7 4.0 4.5 4.6 4.7 5.1 5.4 5.6 5.7 6.2 6.4 6.5 6.7 7.3 7.4 7.5 7.6";
		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper); // MuSSE
		
		Double[] mask1Array = { 0.0, 0.0, 0.0, 0.0, 0.0 }; // first four states, then CID/Not-CID	
		RealParameter mask1 = new RealParameter(mask1Array);
		Double[] mask2Array = { 1.0, 0.0, 0.0, 0.0, 0.0 };
		RealParameter mask2 = new RealParameter(mask2Array);
		Double[] mask3Array = { 1.0, 1.0, 0.0, 0.0, 0.0 };
		RealParameter mask3 = new RealParameter(mask3Array);
		Double[] mask4Array = { 1.0, 1.0, 1.0, 0.0, 0.0 };
		RealParameter mask4 = new RealParameter(mask4Array);
		Double[] mask5Array = { 1.0, 1.0, 1.0, 1.0, 0.0 };
		RealParameter mask5 = new RealParameter(mask5Array);
		Double[] mask6Array = { 1.0, 1.0, 1.0, 1.0, 1.0 };
		RealParameter mask6 = new RealParameter(mask6Array);
		
		Double[] mask7Array = { 0.0, 2.0, 0.0, 0.0, 0.0 };
		RealParameter mask7 = new RealParameter(mask7Array);
		Double[] mask8Array = { 1.0, 2.0, 0.0, 0.0, 0.0 };
		RealParameter mask8 = new RealParameter(mask8Array);
			
		Double[] mask9Array = { 0.0, 2.0, 0.0, 1.0, 1.0 };
		RealParameter mask9 = new RealParameter(mask9Array);
		Double[] mask10Array = { 0.0, 2.0, 0.0, 2.0, 1.0 };
		RealParameter mask10 = new RealParameter(mask10Array);
		
		Double[] mask11Array = { 0.0, 0.0, 0.0, 0.0, 1.0 };
		RealParameter mask11 = new RealParameter(mask11Array);
		
		Double[] mask12Array = { 2.0, 2.0, 2.0, 2.0, 1.0 };
		RealParameter mask12 = new RealParameter(mask12Array);
		
		MasqueradeBall maskBall = new MasqueradeBall();
		maskBall.initByName("modelMask", mask9, "hiddenInstantaneousRateMatrix", hirm, "lambdaMuAssigner", lambdaMuAssigner);
		System.out.println(Arrays.toString(maskBall.getLambdas()));
		System.out.println(Arrays.toString(maskBall.getMus()));
		
//		maskBall.initByName("modelMask", mask6, "hiddenInstantaneousRateMatrix", hirm, "lambdaMuAssigner", lambdaMuAssigner);
//		System.out.println(Arrays.toString(maskBall.getLambdas()));
	}

}
