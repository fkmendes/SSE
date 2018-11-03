//
/*
Compare posterior probabilities at all internal nodes for drawJoint and drawStochastic
 */
package src.test;

import SSE.CladoTriplet;
import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import SSE.CladogeneticSpeciationRateStash;
import SSE.CladoTriplet.speciationType;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.*;

public class SDSEPSampleCLaSSEComparisonTest {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	private StateDependentSpeciationExtinctionProcess sdsep;

	public void runExperiment(String treeStr, String spAttr, String[] spNames, String expName,
							  int numTimeSlices) throws Exception {
		// initializing states
		int numberOfStates = 4; // CLaSSE
		int numSpecies = spNames.length;
		int numNodes = numSpecies * 2 - 1;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", spAttr);
		traitStash.printLksMap();

		// initializing birth-death parameters
		double sympProb = 1.0; // DEC-like
		double subsympProb = 1.0 / 6.0;
		double vicProb = 1.0 / 6.0;
		double jProb = 0.0; // no jump dispersal

		double birthRate = 0.32222224;
		double deathRate = 0.1; // DEC-like

		Double[] mus = {deathRate, deathRate, deathRate, deathRate };
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
		nullTriplet.initByName("parentState", 1,
				"leftChildState", 1,
				"rightChildState", 1,
				"speciationType", speciationType.SYMPATRY);

		CladoTriplet sTriplet1 = new CladoTriplet();
		sTriplet1.initByName("parentState", 2,
				"leftChildState", 2,
				"rightChildState", 2,
				"speciationType", speciationType.SYMPATRY);

		CladoTriplet sTriplet2 = new CladoTriplet();
		sTriplet2.initByName("parentState", 3,
				"leftChildState", 3,
				"rightChildState", 3,
				"speciationType", speciationType.SYMPATRY);

		CladoTriplet jTriplet1 = new CladoTriplet();
		jTriplet1.initByName("parentState", 2,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.JUMPDISPERSAL);

		CladoTriplet jTriplet2 = new CladoTriplet();
		jTriplet2.initByName("parentState", 3,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.JUMPDISPERSAL);

		CladoTriplet vTriplet1 = new CladoTriplet();
		vTriplet1.initByName("parentState", 4,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.VICARIANCE);

		CladoTriplet ssTriplet1 = new CladoTriplet();
		ssTriplet1.initByName("parentState", 4,
				"leftChildState", 2,
				"rightChildState", 4,
				"speciationType", speciationType.SUBSYMPATRY);

		CladoTriplet ssTriplet2 = new CladoTriplet();
		ssTriplet2.initByName("parentState", 4,
				"leftChildState", 3,
				"rightChildState", 4,
				"speciationType", speciationType.SUBSYMPATRY);

		List<CladoTriplet> cladoTripletList = new ArrayList<CladoTriplet>();
		Collections.addAll(cladoTripletList, nullTriplet, sTriplet1, sTriplet2, jTriplet1, jTriplet2, vTriplet1, ssTriplet1, ssTriplet2);

		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("cladoTriplets", cladoTripletList,
				"sympatricRate", sympatricSpeciationRate,
				"subsympatricRate", subSympatricSpeciationRate,
				"vicariantRate", vicariantSpeciationRate,
				"jumpRate", jumpSpeciationRate);
		csrt.printEventMap();

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		String FlatQMatrixString = "0.0 0.0 0.0 0.01 0.0 0.01 0.01 0.0 0.01 0.0 0.01 0.01";
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", FlatQMatrixString);
		irm.printMatrix();

		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		// fixed equilibrium frequency
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);

		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		// tree is given

		boolean incorporateCladogenesis = true;

		sdsep = new StateDependentSpeciationExtinctionProcess();
		sdsep.initByName(
				"tree", myTree,
				"traitStash", traitStash,
				"instantaneousRateMatrix", irm,
				"cladogeneticStash", csrt,
				"mu", mu,
				"pi", pi,
				"incorporateCladogenesis", incorporateCladogenesis
		);

		// Sample many times with drawStochasticChar and calculate the posterior
		sdsep.setNumTimeSlices(numTimeSlices);
		double[][] posteriorStoc = sdsep.sampleAndSummarizeCLaSSE(numTrials, false);
		double[][] posteriorJoint = sdsep.sampleAndSummarizeCLaSSE(numTrials, true);

		for (int i = 0; i < numNodes; i++) {
			System.out.println("Stoch for node " + i + Arrays.toString(posteriorStoc[i]));
			System.out.println("Joint for node " + i + Arrays.toString(posteriorJoint[i]));
		}

		TestHelper.compareNestedArr(posteriorJoint, posteriorStoc);
	}

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void test() throws Exception {
		String treeStr = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
		String spAttr = "Human=3,Chimp=1,Gorilla=4,Orang=3";

		runExperiment(treeStr, spAttr, spNames, "CLaSSE", 50);
	}

}
