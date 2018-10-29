package SSE;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

public class HiddenStateDependentSpeciationExtinctionProcessHiSSETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 2;
		int numberOfHiddenStates = 2;
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames1 = new String[] { "sp4", "sp6", "sp10", "sp11", "sp12", "sp15", "sp16", "sp17", "sp18", "sp20", "sp21", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" }; // test 1
		String[] spNames2 = new String[] { "sp6", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp17", "sp18", "sp19", "sp21", "sp22", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42", "sp43", "sp44", "sp45", "sp46", "sp47", "sp48", "sp49", "sp50", "sp51", "sp52", "sp53", "sp54", "sp55", "sp56", "sp57", "sp58", "sp59", "sp60" }; // test 2
		List<Taxon> taxaList1 = Taxon.createTaxonList(Arrays.asList(spNames1));
		List<Taxon> taxaList2 = Taxon.createTaxonList(Arrays.asList(spNames2));
		TaxonSet taxonSet1 = new TaxonSet(taxaList1);
		TaxonSet taxonSet2 = new TaxonSet(taxaList2);
		
		HiddenTraitStash hiddenTraitStash1 = new HiddenTraitStash();
		HiddenTraitStash hiddenTraitStash2 = new HiddenTraitStash();
		hiddenTraitStash1.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet1, "value", "sp4=1,sp6=2,sp10=2,sp11=1,sp12=1,sp15=1,sp16=2,sp17=1,sp18=2,sp20=2,sp21=1,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=2,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		hiddenTraitStash2.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet2, "value", "sp6=1,sp10=1,sp11=1,sp12=1,sp14=2,sp15=1,sp16=2,sp17=2,sp18=1,sp19=1,sp21=1,sp22=1,sp23=1,sp24=1,sp25=1,sp26=2,sp27=2,sp28=1,sp29=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=1,sp36=1,sp37=2,sp38=2,sp39=1,sp40=1,sp41=1,sp42=2,sp43=2,sp44=2,sp45=2,sp46=2,sp47=2,sp48=2,sp49=1,sp50=1,sp51=2,sp52=2,sp53=1,sp54=1,sp55=2,sp56=2,sp57=1,sp58=1,sp59=1,sp60=1");
		
		hiddenTraitStash1.printLksMap();
		hiddenTraitStash2.printLksMap();
		
		Double lambda1 = 0.2860762;
		Double lambda2 = 0.1313947;
		Double lambda3 = 0.01910449;
		Double lambda4 = 2.116117;
		Double[] lambdas1 = { lambda1, lambda2, lambda3, lambda4 }; // 0A, 1A, 0B, 1B
		System.out.println("Lambdas test 1: " + Arrays.toString(lambdas1));
		
		lambda1 = 6.440928e-10;
		lambda2 = 0.1128975;
		lambda3 = 0.1308483;
		lambda4 = 85.03492;
		Double[] lambdas2 = { lambda1, lambda2, lambda3, lambda4 }; // 0A, 1A, 0B, 1B
		System.out.println("Lambdas test 2: " + Arrays.toString(lambdas2));
		
		Double mu1 = 5.89647E-10;
		Double mu2 = 0.0007677269;
		Double mu3 = 3.937728E-11;
		Double mu4 = 4.361642E-09;
		Double[] mus1 = { mu1, mu2, mu3, mu4 };		
		System.out.println("Mus test 1: " + Arrays.toString(mus1));
		
		mu1 = 1.417061E-09;
		mu2 = 0.04597049;
		mu3 = 0.03687651;
		mu4 = 1.266694E-06;
		Double[] mus2 = { mu1, mu2, mu3, mu4 };
		System.out.println("Mus test 2: " + Arrays.toString(mus2));
		
		RealParameter muTest1 = new RealParameter(mus1);
		RealParameter muTest2 = new RealParameter(mus2);
		muTest1.initByName("minordimension", 1);
		muTest2.initByName("minordimension", 1);
		
		RealParameter lambdaTest1 = new RealParameter(lambdas1);
		RealParameter lambdaTest2 = new RealParameter(lambdas2);

		HiddenInstantaneousRateMatrix hirm1 = new HiddenInstantaneousRateMatrix();
		HiddenInstantaneousRateMatrix hirm2 = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString1 = "2.078263E-09 0.8410903 3.315091E-02 2.061154E-09 2.061154E-09 2.061154E-09 2.383310 0.4754780"; // test 1
		String flatQMatrixString2 = "2.11773E-09 2.061154E-09 6.432635E-03 2.479683E-09 2.355085E-09 1.565876E-02 1.00000E+02 1.000000E+02"; // test 2
		String hiddenStatesString = "0,1"; // observed state 1 and 2 will transition to hidden states 1 and 2, respectively
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
		boolean disallowDoubleTransitions = true;
		hirm1.initByName("numberOfStates", 2, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString1, "disallowDoubleTransitions", disallowDoubleTransitions, "hiddenObsStateMapper", stateMapper); // HiSSE
		hirm1.printMatrix();
		hirm2.initByName("numberOfStates", 2, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString2, "disallowDoubleTransitions", disallowDoubleTransitions, "hiddenObsStateMapper", stateMapper); // HiSSE
		hirm2.printMatrix();
		
		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr1 = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;"; // test 1
		String treeStr2 = "(((sp15:10.27880339,(sp57:0.4327353378,sp58:0.4327353378):9.846068053):21.30935137,((((sp49:1.322566942,sp50:1.322566942):6.531246386,(((((sp42:1.618558172,sp43:1.618558172):1.249323508,sp37:2.86788168):0.4105311845,sp36:3.278412865):1.110829025,sp28:4.38924189):2.453996398,((sp53:0.6765630317,sp54:0.6765630317):5.834067793,sp21:6.510630824):0.3326074635):1.01057504):6.546385565,sp12:14.40019889):3.891878236,((((sp18:8.595427361,((sp19:6.988162304,((sp39:1.941330272,(sp59:0.4256083779,sp60:0.4256083779):1.515721894):1.374985348,sp35:3.31631562):3.671846684):1.028692949,(sp24:5.527011086,(sp25:5.478875203,(sp40:1.898502308,sp41:1.898502308):3.580372894):0.04813588287):2.489844168):0.5785721075):0.8605508177,((sp47:1.324188282,sp48:1.324188282):1.210143714,sp38:2.534331996):6.921646183):1.848794077,(sp22:6.144323416,sp23:6.144323416):5.160448839):4.752352041,sp10:16.0571243):2.234952832):13.29607763):8.9940146,(sp6:33.80408947,(((sp29:4.271294196,sp30:4.271294196):3.963360008,(sp46:1.515605972,(sp51:0.6842469553,sp52:0.6842469553):0.8313590168):6.719048232):21.69107479,((((sp44:1.517683119,sp45:1.517683119):13.83340518,((sp33:3.451233406,sp34:3.451233406):7.318030201,sp14:10.76926361):4.581824694):2.3268441,((sp31:3.988873926,sp32:3.988873926):13.39833,(sp26:5.46221229,sp27:5.46221229):11.92499164):0.2907284735):12.10203097,((sp16:9.676541191,sp17:9.676541191):11.55054389,(sp11:16.00734921,(sp55:0.6152478573,sp56:0.6152478573):15.39210136):5.219735869):8.552878292):0.1457656227):3.878360468):6.778079891):0.0;"; // test 2
        TreeParser myTree1 = new TreeParser(treeStr1, false, false, true, 0); // true b/c species are labelled, offset=0
        TreeParser myTree2 = new TreeParser(treeStr2, false, false, true, 0);
        
        boolean incorporateCladogenesis = false;
                
        HiddenStateDependentSpeciationExtinctionProcess hsdsep1 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep1.initByName(
        		"tree", myTree1,
        		"hiddenTraitStash", hiddenTraitStash1,
        		"hiddenInstantaneousRateMatrix", hirm1,
        		"lambda", lambdaTest1,
        		"mu", muTest1,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        HiddenStateDependentSpeciationExtinctionProcess hsdsep2 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep2.initByName(
        		"tree", myTree2,
        		"hiddenTraitStash", hiddenTraitStash2,
        		"hiddenInstantaneousRateMatrix", hirm2,
        		"lambda", lambdaTest2,
        		"mu", muTest2,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
        System.out.println(hsdsep1.calculateLogP());
        
    	System.out.println(hsdsep2.calculateLogP());
	}
}
