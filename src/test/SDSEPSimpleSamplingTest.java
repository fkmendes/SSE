package src.test;

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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class SDSEPSimpleSamplingTest {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	private StateDependentSpeciationExtinctionProcess sdsep;

	public void runExperiment(String treeStr, String spAttr, String[] spNames, String expName) throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE
		int numSpecies = spNames.length;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", spAttr);
		traitStash.printLksMap();

		// initializing birth-death parameters
		Double[] mus = { 0.01, 0.01 };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols

		Double[] lambdas = {0.08, 0.08};
		RealParameter lambda = new RealParameter(lambdas);

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.01 0.01");
		irm.printMatrix();

		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		// fixed equilibrium frequency
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);

		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		// tree is given

		boolean incorporateCladogenesis = false;

		sdsep = new StateDependentSpeciationExtinctionProcess();
		sdsep.initByName(
				"tree", myTree,
				"traitStash", traitStash,
				"instantaneousRateMatrix", irm,
				"lambda", lambda,
				"mu", mu,
				"pi", pi,
				"incorporateCladogenesis", incorporateCladogenesis
		);

		// Sample many times with drawStochasticChar and calculate the posterior
		sdsep.setNumTimeSlices(500);
		double[] posteriorStoc = sdsep.sampleAndSummarize(numTrials, false);
		double[] posteriorJoint = sdsep.sampleAndSummarize(numTrials, true);
//		TestHelper.compareArr(posteriorStoc, posteriorJoint);

		double[] posterior = posteriorStoc;
//		double[] posterior = posteriorJoint;

		// Write only the ancestral states to csv
		posterior = TestHelper.trimTips(posterior, numSpecies);
		String dir = "/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/scm";
		String fileName = expName + ".csv";
		TestHelper.writeToCSV(dir, fileName, posterior, sdsep);

//		String[] divLbls = {"nd2","nd3","nd5","nd4"};
//		String[] divLks = {"0.87785842562904", "0.829891574906887", "0.00242853993720819", "0.995986779987027"};
//		String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
//		TestHelper.compareDiv(divLbls, divLks, indexNameMapper, posterior);
	}

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void test() throws Exception {
		String treeStr = "(((sp5:1.753921149,sp6:1.753921149)nd5:10.54206596,sp2:12.2959871)nd3:5.60266132,(sp3:6.029064844,sp4:6.029064844)nd4:11.86958358)nd2;";
		String spAttr = "sp2=1,sp3=1,sp4=1,sp5=2,sp6=2";
		String[] spNames = new String[] {"sp2","sp3","sp4","sp5","sp6"};
		runExperiment(treeStr, spAttr, spNames, "beast_small");
		Assert.assertEquals(-17.97332, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-17.97332, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver

		treeStr = "(((((((sp23:0.1790943213,sp24:0.1790943213)nd24:1.092201813,(sp25:0.07202419731,sp26:0.07202419731)nd25:1.199271937)nd13:11.36119223,sp4:12.63248837)nd11:1.446039834,(((sp9:10.42249592,sp10:10.42249592)nd15:0.404458011,sp6:10.82695394)nd14:0.008741084797,sp5:10.83569502)nd12:3.24283318)nd8:5.258830888,(sp7:10.79794068,sp8:10.79794068)nd9:8.539418405)nd5:10.54206596,((((sp13:6.968607434,(sp21:0.8635028776,sp22:0.8635028776)nd22:6.105104556)nd19:0.7298402864,(sp15:7.029426078,sp12:7.029426078)nd20:0.6690216423)nd18:2.616431804,(sp19:1.17454796,sp20:1.17454796)nd23:9.140331564)nd10:7.268558415,sp3:17.58343794)nd6:12.2959871)nd3:5.60266132,((sp17:2.68800466,sp18:2.68800466)nd7:20.92449812,sp2:23.61250278)nd4:11.86958358)nd2;";
		spAttr = "sp2=1,sp3=2,sp4=2,sp5=2,sp6=2,sp7=2,sp8=2,sp9=2,sp10=2,sp12=1,sp13=1,sp15=1,sp17=1,sp18=1,sp19=1,sp20=1,sp21=1,sp22=1,sp23=2,sp24=2,sp25=2,sp26=2";
		spNames = new String[] {"sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp12","sp13","sp15","sp17","sp18","sp19","sp20","sp21","sp22","sp23","sp24","sp25","sp26"};
		runExperiment(treeStr, spAttr, spNames, "beast_large");
		Assert.assertEquals(-81.12881, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-81.12881, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
	}

}
