package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import biogeo.InstantaneousRateMatrix;
import biogeo.StateDependentSpeciationExtinctionProcess;
import biogeo.TraitStash;

public class SDSEPBiSSETest {
	final static double EPSILON = 1e-10;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("NumberOfStates", numberOfStates, "taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2");
		traitStash.printLksMap();
				
		// initializing birth-death parameters
		Double birthRate = 0.222222222;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] lambdas = new Double[numberOfStates];
		Arrays.fill(lambdas, birthRate);
		RealParameter lambda = new RealParameter(lambdas);
	
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("NumberOfStates", numberOfStates, "FlatQMatrix", "0.1 0.1");
		irm.printMatrix();
		
		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
        
        sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"TreeParser", myTree,
        		"TraitStash", traitStash,
        		"InstantaneousRateMatrix", irm,
        		"Lambda", lambda,
        		"Mu", mu,
        		"Pi", pi,
        		"IncorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	System.out.println(sdsep.calculateLogP());
	}

	@Test
	public void test() {
		Assert.assertEquals(-5.588460032653, sdsep.calculateLogP(), EPSILON);
	}

}
