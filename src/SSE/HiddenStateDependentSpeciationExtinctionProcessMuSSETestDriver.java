package SSE;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

public class HiddenStateDependentSpeciationExtinctionProcessMuSSETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 4;
		int numberOfHiddenStates = 0;
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames = new String[] { "sp4", "sp6", "sp10", "sp11", "sp12", "sp15", "sp16", "sp17", "sp18", "sp20", "sp21", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "value", "sp4=1,sp6=4,sp10=2,sp11=3,sp12=1,sp15=1,sp16=4,sp17=1,sp18=4,sp20=2,sp21=3,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=4,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		hiddenTraitStash.printLksMap();
		
		Double lambda1 = 0.1;
		Double lambda2 = 0.15;
		Double lambda3 = 0.2;
		Double lambda4 = 0.1;
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4 };
		System.out.println("Lambdas: " + Arrays.toString(lambdas));
		
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double[] mus = { mu1, mu2, mu3, mu4 };
		
		System.out.println("Mus: " + Arrays.toString(mus));
		
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		RealParameter lambda = new RealParameter(lambdas);

		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.05 0.05 0.0 0.05 0.0 0.05 0.05 0.0 0.05 0.0 0.05 0.05"; // MuSSE
		String hiddenStatesString = "0,1"; // observed state 1 and 2 will transition to hidden states 1 and 2, respectively
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
		boolean disallowDoubleTransitions = false;

		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 0, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "hiddenObsStateMapper", stateMapper); // MuSSE
		hirm.printMatrix();
		
		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25String treeStr = "(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,(sp4:1.0,sp5:1.0):2.0)0.0;";
		
		
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        HiddenStateDependentSpeciationExtinctionProcess hsdsep = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash,
        		"hiddenInstantaneousRateMatrix", hirm,
        		"lambda", lambda,
        		"mu", mu,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	System.out.println(hsdsep.calculateLogP()); // -122.88014232305367
	}
}
