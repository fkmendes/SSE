package test;

import java.util.ArrayList;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.LambdaMuAssigner;
import SSE.MasqueradeBall;
import beast.core.parameter.RealParameter;

public class MasqueradeBallTest {
	
	private MasqueradeBall maskBall;

	HiddenObservedStateMapper stateMapper;
	String hiddenStatesString;
	String flatQMatrixString;
	String lambdasToStatesString;
	String musToStatesString;
	RealParameter lambda;
	RealParameter mu;
	boolean disallowDoubleTransitions;
	int symmetrifyAcrossDiagonal;
	
	RealParameter mask1;
	RealParameter mask2;
	RealParameter mask3;
	RealParameter mask4;
	RealParameter mask5;
	RealParameter mask6;
	ArrayList<RealParameter> masks;
	
	public interface Instance {
        Double[] getLambdas();

        Double[] getMus();

        Double[][] getQ();
    }
	
	Instance test1 = new Instance() {
		@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, .2, .3, .4 }, 
				{ 1.1, 0.0, 1.3, 1.4 }, 
				{ 2.1, 2.2, 0.0, 2.4 }, 
				{ 3.1, 3.2, 3.3, 0.0 } 
			};
		}
    };
    
    Instance test2 = new Instance() {
		@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .1 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .03 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, .2, .3, .4, .5 }, 
				{ 1.1, 0.0, 1.3, 1.4, 0.0 }, 
				{ 2.1, 2.2, 0.0, 2.4, 0.0 }, 
				{ 3.1, 3.2, 3.3, 0.0, 0.0 },
				{ .5, 0.0, 0.0, 0.0, 0.0 }
			};
		}
    };
    
    Instance test3 = new Instance() {
		@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .1, .15 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .03, .045 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, .2, .3, .4, .5, 0.0 }, 
				{ 1.1, 0.0, 1.3, 1.4, 0.0, 1.6 }, 
				{ 2.1, 2.2, 0.0, 2.4, 0.0, 0.0 }, 
				{ 3.1, 3.2, 3.3, 0.0, 0.0, 0.0 },
				{ .5, 0.0, 0.0, 0.0, 0.0, 4.6 },
				{ 0.0, 5.2, 0.0, 0.0, 5.5, 0.0 },
			};
		}
    };
    
	@Before
	public void setUp() throws Exception {
				
		hiddenStatesString = "0,1,2,3";
		stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
				
		lambdasToStatesString = "0,1,2,3,4,5,6,7";
		Double lambda1 = 0.1; // 0A
		Double lambda2 = 0.15; // 1A
		Double lambda3 = 0.2; // 2A
		Double lambda4 = 0.1; // 3A
		Double lambda5 = 0.1; // 0B
		Double lambda6 = 0.15; // 1B
		Double lambda7 = 0.2; // 2B
		Double lambda8 = 0.1; // 3B
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8 };
		lambda = new RealParameter(lambdas);
		
		musToStatesString = "0,1,2,3,4,5,6,7";
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double mu5 = 0.03;
		Double mu6 = 0.045;
		Double mu7 = 0.06;
		Double mu8 = 0.03;
		Double[] mus = { mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8 };
		mu = new RealParameter(mus);
		
		disallowDoubleTransitions = true; // not used
		symmetrifyAcrossDiagonal = -1;
		flatQMatrixString = "0.2 0.3 0.4 0.5 1.1 1.3 1.4 1.6 2.1 2.2 2.4 2.7 3.1 3.2 3.3 3.8 4.1 4.6 4.7 4.8 5.2 5.5 5.7 5.8 6.3 6.5 6.6 6.8 7.8 7.5 7.6 7.7";		
		
		Double[] mask1Array = { 0.0, 0.0, 0.0, 0.0, 0.0 }; // first four states, then CID/Not-CID	
		mask1 = new RealParameter(mask1Array);
		
		Double[] mask2Array = { 1.0, 0.0, 0.0, 0.0, 0.0 };
		mask2 = new RealParameter(mask2Array);
		
		Double[] mask3Array = { 1.0, 1.0, 0.0, 0.0, 0.0 };
		mask3 = new RealParameter(mask3Array);
		
		Double[] mask4Array = { 1.0, 1.0, 1.0, 0.0, 0.0 };
		mask4 = new RealParameter(mask4Array);
		
		Double[] mask5Array = { 1.0, 1.0, 1.0, 1.0, 0.0 };
		mask5 = new RealParameter(mask5Array);
		
		Double[] mask6Array = { 1.0, 1.0, 1.0, 1.0, 1.0 };
		mask5 = new RealParameter(mask6Array);
		
		masks = new ArrayList<RealParameter>();
		masks.add(mask1);
		masks.add(mask2);
		masks.add(mask3);
		masks.add(mask4);
		masks.add(mask5);
		masks.add(mask6);
	}

	Instance[] all = { test1, test2, test3 };
	
	@Test
	public void testMaskeradeBall() {
		
		int i = 0;
		for (Instance test : all) {
			LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
			lambdaMuAssigner.initByName("totalNumberOfStates", 8, "nDistinctLambdas", 8, "nDistinctMus", 8, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu);
			
			HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
			hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper);
			
			maskBall = new MasqueradeBall();
			maskBall.initByName("modelMask", masks.get(i), "hiddenInstantaneousRateMatrix", hirm, "lambdaMuAssigner", lambdaMuAssigner);
	
			Assert.assertArrayEquals(test.getLambdas(), maskBall.getLambdas());
			Assert.assertArrayEquals(test.getMus(), maskBall.getMus());
			Assert.assertArrayEquals(test.getQ(), maskBall.getQs());

			++i;
		}
	}
	
}
