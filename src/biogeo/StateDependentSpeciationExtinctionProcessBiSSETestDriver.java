package biogeo;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcessBiSSETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numStates = 2; // BiSSE
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash(numStates);
		traitStash.initByName("taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2");
		traitStash.printLksMap();
		
		Double birthRate = 0.222222222;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] lambdas = new Double[numStates];
		Arrays.fill(lambdas, birthRate);
		RealParameter lambda = new RealParameter(lambdas);
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("NumberOfStates", numStates, "FlatQMatrix", "0.9 0.1 0.9 0.1");
		irm.printMatrix();
		
		Double[] piEs = new Double[numStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numStates];
		Arrays.fill(piDs, (1.0/numStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        StateDependentSpeciationExtinctionProcess sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"TreeParser", myTree,
        		"TraitStash", traitStash,
        		"InstantaneousRateMatrix", irm,
        		"Lambda", lambda,
        		"Mu", mu,
        		"Pi", pi,
        		"IncorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	sdsep.computeLk();
	}
}
