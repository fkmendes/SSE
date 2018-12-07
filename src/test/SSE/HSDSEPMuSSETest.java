package test.SSE;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenStateDependentSpeciationExtinctionProcess;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

public class HSDSEPMuSSETest {
	final static double EPSILON = 1E-4;
	private double negLnl;

	@Before
	public void setUp() throws Exception {
		int numberOfStates = 4;
		int numberOfHiddenStates = 0;
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames = new String[] { "sp4", "sp6", "sp10", "sp11", "sp12", "sp15", "sp16", "sp17", "sp18", "sp20", "sp21", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper(); // don't need to put anything in it (see later how to make this optional, etc.)
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp4=1,sp6=4,sp10=2,sp11=3,sp12=1,sp15=1,sp16=4,sp17=1,sp18=4,sp20=2,sp21=3,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=4,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		hiddenTraitStash.printLksMap();
		
		String lambdasToStatesString = "0,1,2,3"; // first lambda to first state, second lambda to second state, and so no...
		Double lambda1 = 0.1;
		Double lambda2 = 0.15;
		Double lambda3 = 0.2;
		Double lambda4 = 0.1;
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4 };
		RealParameter lambda = new RealParameter(lambdas);
		
		String musToStatesString = "0,1,2,3";
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double[] mus = { mu1, mu2, mu3, mu4 };
		RealParameter mu = new RealParameter(mus);
		
		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		RealParameter pi = new RealParameter(pis);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 4, "nDistinctLambdas", 4, "nDistinctMus", 4, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu, "pi", pi);
		System.out.println("Lambdas: " + Arrays.toString(lambdaMuAssigner.getLambdas()));
		System.out.println("Mus: " + Arrays.toString(lambdaMuAssigner.getMus()));

		boolean disallowDoubleTransitions = false; // not used
		int symmetrifyAcrossDiagonal = -1;
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.05 0.05 0.0 0.05 0.0 0.05 0.05 0.0 0.05 0.0 0.05 0.05"; // passing the 0.0s of "double transitions" explicitly (even though there are no hidden states)
		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 0, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper); // MuSSE
		
		System.out.println("Qs:");
		hirm.printMatrix();
		
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        HiddenStateDependentSpeciationExtinctionProcess hsdsep = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash,
        		"lambdaMuAssigner", lambdaMuAssigner,
        		"hiddenInstantaneousRateMatrix", hirm,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	negLnl = hsdsep.calculateLogP();
    	System.out.println(negLnl); // -122.88014232305267
	}

	@Test
	public void againstDiversitreeMuSSE() {
		Assert.assertEquals(-122.8801, negLnl, EPSILON); 
	}

	@Test
	public void againstMyMuHiSSE() {
		Assert.assertEquals(-122.88014179920914, negLnl, EPSILON);
	}

}
