/*
Under CLaSSE, we sample the nodes of a fixed tree many times and get a marginal posterior distribution on the states for each node
Then, for each internal node, we get a 95% credible set - the smallest set of states such that it contains 95% of the probability
For each node, we see if the true state lies in the 95% credible set
We report how frequent the truth lies in the credible set (I call this accuracy tho it may not be the best word in context)
Do this for many trees with Joint and Stoc. Our accuracy better be 95%. (or whatever threshold we set)
 */
package test;

import SSE.CladogeneticSpeciationRateStash;
import SSE.CladoTriplet.speciationType;
import SSE.CladoTriplet;
import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Array;
import java.util.*;

public class SDSEPSampleCLaSSECredibleSetTest {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	private StateDependentSpeciationExtinctionProcess sdsep;

	public double runExperiment(String treeStr, String spAttr, String[] spNames,
								int numTimeSlices, String[] divLbls, String[] divStates,
								boolean useJoint, double credibleThreshold) throws Exception {
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

		// Sample many times, compute posterior probabilities and construct credible set
		sdsep.setNumTimeSlices(numTimeSlices);
		double[][] posterior = sdsep.sampleAndSummarizeCLaSSE(numTrials, useJoint);
		posterior = TestHelper.trimTipsCLaSSE(posterior);
		System.out.println("Node 4 posterior: " + Arrays.toString(posterior[3]));
		Set<Integer>[] credibleSets = TestHelper.constructCredibleSets(posterior, credibleThreshold);

		// Compare credible set with truth
		HashMap<String, Double> divMap = TestHelper.getDivMap(divLbls, divStates);
		String[] idxLabelMapper = sdsep.getNodeIndexNameMapper();
		double accuracy = TestHelper.computeProportionTruthInCredibleSet(divMap, credibleSets, idxLabelMapper);
		System.out.println("Proportion of internal nodes in their credible sets with joint: " + accuracy);

		return accuracy;
	}

	@Before
	public void setUp() throws Exception {
	}

	public double runWrapper(String testName, double credibleThreshold) throws Exception {
		String treeSuffix = ".tree";
		String spAttrSuffix = "-beast_str.txt";
		String tipSuffix = "-tips.csv";
		String nodeSuffix = "-node_truth.csv";

		String treeStr = TestHelper.readFirstLine(testName + treeSuffix);
		String spAttr = TestHelper.readFirstLine(testName + spAttrSuffix);
		String tips = TestHelper.readFirstLine(testName + tipSuffix);
		ArrayList<String> divData = TestHelper.readFile(testName + nodeSuffix);

		String[] spNames = tips.split(",");
		String[] divLbls = divData.get(0).split(",");
		String[] divStates = divData.get(1).split(",");

		return runExperiment(treeStr, spAttr, spNames, 50, divLbls, divStates, true, credibleThreshold);
	}

	@Test
	public void test() throws Exception {
	    String baseTestName = "data/test/test"; // dir + exp name
	    int numTrees = 10;
		double credibleThreshold = 0.7;

		double totalAcc = 0.0;
		for (int i = 1; i <= numTrees; i++) {
			totalAcc += runWrapper(baseTestName + i, credibleThreshold);
		}
		double avgAcc = totalAcc / numTrees;
		System.out.println(avgAcc);

		Assert.assertTrue(credibleThreshold < avgAcc + EPSILON);


		// This test was to debug the system
//		String treeStr = "(((Human:1.0,Chimp:1.0)nd1:1.0,Gorilla:2.0)nd2:1.0,Orang:3.0)nd3;";
//		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
//		String spAttr = "Human=3,Chimp=1,Gorilla=2,Orang=4";
//		// truth data is fake right now
//		String[] divLbls = new String[] {"nd1", "nd2", "nd3"};
//		String[] divStates = new String[] {"1","0","3"};
//
//		double avgAcc = runExperiment(treeStr, spAttr, spNames, 100, divLbls, divStates, false, credibleThreshold);
//
//		Assert.assertEquals(credibleThreshold, avgAcc, EPSILON);
	}

}
