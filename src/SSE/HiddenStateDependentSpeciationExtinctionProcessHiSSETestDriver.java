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
		String[] spNames = new String[] { "sp1", "sp2", "sp3", "sp4", "sp5" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "value", "sp1=1,sp2=1,sp3=2,sp4=1,sp5=2");
		hiddenTraitStash.printLksMap();
		
		Double birthRate = 0.2;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate, deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] lambdas = new Double[totalNumberOfStates];
		Arrays.fill(lambdas, birthRate);
		RealParameter lambda = new RealParameter(lambdas);

		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1";
		String hiddenStatesString = "0,1"; // observed state 1 and 2 will transition to hidden states 1 and 2, respectively
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
		hirm.initByName("numberOfStates", 2, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", true, "hiddenObsStateMapper", stateMapper);
		hirm.printMatrix();
		
		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,(sp4:1.0,sp5:1.0):3.0)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
//        HiddenStateDependentSpeciationExtinctionProcess hsdsep = new HiddenStateDependentSpeciationExtinctionProcess();
//        hsdsep.initByName(
//        		"tree", myTree,
//        		"traitStash", traitStash,
//        		"instantaneousRateMatrix", irm,
//        		"lambda", lambda,
//        		"mu", mu,
//        		"pi", pi,
//        		"incorporateCladogenesis", incorporateCladogenesis,
//              "disallowDoubleTransitions", true
//        		);
//    	
//    	System.out.println(hsdsep.calculateLogP());
	}
}
