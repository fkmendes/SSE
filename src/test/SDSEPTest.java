package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import biogeo.CladoTriplet;
import biogeo.CladogeneticSpeciationRateStash;
import biogeo.InstantaneousRateMatrix;
import biogeo.StateDependentSpeciationExtinctionProcess;
import biogeo.TraitStash;
import biogeo.CladoTriplet.speciationType;

public class SDSEPTest {
	final static double EPSILON = 1e-10;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
		int numStates = 4; // state = 1 (index 0) is 'null state'
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash(numStates);
		traitStash.initByName("taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2,Orang=3");
		traitStash.printLksMap();
				
		// initializing birth-death parameters
		double sympProb = 1.0; // DEC-like
		double subsympProb = 1.0 / 6.0;
		double vicProb = 1.0 / 6.0;
		double jProb = 0.0; // no jump dispersal
		
		double birthRate = 0.32222224;
		double deathRate = 0.1; // DEC-like
				
		Double[] mus = { deathRate, deathRate, deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] sSpeciationRate = {sympProb * birthRate};
		Double[] ssSpeciationRate = {subsympProb * birthRate};
		Double[] vSpeciationRate = {vicProb * birthRate};
		Double[] jSpeciationRate = {jProb * birthRate}; // 0.0
		RealParameter sympatricSpeciationRate = new RealParameter(sSpeciationRate);
		RealParameter subSympatricSpeciationRate = new RealParameter(ssSpeciationRate);
		RealParameter vicariantSpeciationRate = new RealParameter(vSpeciationRate);
		RealParameter jumpSpeciationRate = new RealParameter(jSpeciationRate);
		
		CladoTriplet nullTriplet = new CladoTriplet();
		nullTriplet.initByName("ParentState", 1,
				"LeftChildState", 1,
				"RightChildState", 1,
				"SpeciationType", speciationType.SYMPATRY);
		
		CladoTriplet sTriplet1 = new CladoTriplet();
		sTriplet1.initByName("ParentState", 2,
				"LeftChildState", 2,
				"RightChildState", 2,
				"SpeciationType", speciationType.SYMPATRY);
		
		CladoTriplet sTriplet2 = new CladoTriplet();
		sTriplet2.initByName("ParentState", 3,
				"LeftChildState", 3,
				"RightChildState", 3,
				"SpeciationType", speciationType.SYMPATRY);
		
		CladoTriplet jTriplet1 = new CladoTriplet();
		jTriplet1.initByName("ParentState", 2,
				"LeftChildState", 2,
				"RightChildState", 3,
				"SpeciationType", speciationType.JUMPDISPERSAL);
		
		CladoTriplet jTriplet2 = new CladoTriplet();
		jTriplet2.initByName("ParentState", 3,
				"LeftChildState", 2,
				"RightChildState", 3,
				"SpeciationType", speciationType.JUMPDISPERSAL);
		
		CladoTriplet vTriplet1 = new CladoTriplet();
		vTriplet1.initByName("ParentState", 4,
				"LeftChildState", 2,
				"RightChildState", 3,
				"SpeciationType", speciationType.VICARIANCE);
		
		CladoTriplet ssTriplet1 = new CladoTriplet();
		ssTriplet1.initByName("ParentState", 4,
				"LeftChildState", 2,
				"RightChildState", 4,
				"SpeciationType", speciationType.SUBSYMPATRY);
		
		CladoTriplet ssTriplet2 = new CladoTriplet();
		ssTriplet2.initByName("ParentState", 4,
				"LeftChildState", 3,
				"RightChildState", 4,
				"SpeciationType", speciationType.SUBSYMPATRY);
				
		List<CladoTriplet> cladoTripletList = new ArrayList<CladoTriplet>();
		Collections.addAll(cladoTripletList, nullTriplet, sTriplet1, sTriplet2, jTriplet1, jTriplet2, vTriplet1, ssTriplet1, ssTriplet2);
		
		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("CladoTriplets", cladoTripletList,
				"SympatricRate", sympatricSpeciationRate,
				"SubsympatricRate", subSympatricSpeciationRate,
				"VicariantRate", vicariantSpeciationRate,
				"JumpRate", jumpSpeciationRate);
		csrt.printEventMap();
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		String FlatQMatrixString = "0.0 0.0 0.0 0.0 0.01 0.0 0.0 0.01 0.01 0.0 0.0 0.01 0.0 0.01 0.01 0.00";
		                           
		irm.initByName("NumberOfStates", numStates, "FlatQMatrix", FlatQMatrixString);
		irm.printMatrix();
		
		Double[] piEs = new Double[numStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numStates];
		Arrays.fill(piDs, 1.0/(numStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		
        boolean incorporateCladogenesis = true;
        
        StateDependentSpeciationExtinctionProcess sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"TreeParser", myTree,
        		"TraitStash", traitStash,
        		"InstantaneousRateMatrix", irm,
        		"CladogeneticStash", csrt,
        		"Mu", mu,
        		"Pi", pi,
        		"IncorporateCladogenesis", incorporateCladogenesis
        		);
        
        System.out.println(sdsep.calculateLogP());
	}

	@Test
	public void test() {
		Assert.assertEquals(-10.59346884351, sdsep.calculateLogP(), EPSILON);
	}

}
