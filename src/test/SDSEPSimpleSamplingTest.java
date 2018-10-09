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
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		String treeStr = "(((sp5:0.1753921149,sp6:0.1753921149)nd5:1.054206596,sp2:1.22959871)nd3:0.560266132,(sp3:0.6029064844,sp4:0.6029064844)nd4:1.186958358)nd2;";
		String spAttr = "sp2=1,sp3=1,sp4=1,sp5=2,sp6=2";
		// initializing states
		int numberOfStates = 2; // BiSSE

		String[] spNames = new String[] {"sp2","sp3","sp4","sp5","sp6"};
		int numSpecies = spNames.length;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", spAttr);
		traitStash.printLksMap();

		// initializing birth-death parameters
		Double[] mus = { 0.1, 0.1 };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols

		Double[] lambdas = {0.8, 0.8};
		RealParameter lambda = new RealParameter(lambdas);

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.1 0.1");
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
		double[] posteriorStoc = sdsep.sampleAndSummarize(10000, false);
//		double[] posteriorJoint = sdsep.sampleAndSummarize(10000, true);
//		compareArr(posteriorStoc, posteriorJoint);

		double[] posterior = posteriorStoc;
//		double[] posterior = posteriorJoint;

//		 Write only the ancestral states to csv
		posterior = TestHelper.trimTips(posterior, numSpecies);
		TestHelper.writeToCSV("beast.csv", posterior, sdsep);

		String[] divLbls = {"nd2","nd3","nd5","nd4"};
		String[] divLks = {"0.87785842562904", "0.829891574906887", "0.00242853993720819", "0.995986779987027"};
		String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
		TestHelper.compareDiv(divLbls, divLks, indexNameMapper, posterior);
	}

	@Test
	public void test() {
		Assert.assertEquals(-8.311203, sdsep.calculateLogP(), 1e-1); // Used in original version with fixed-step size ODE solver
//		sdsep.setSampleCharacterHistory(true);
//		Assert.assertEquals(-8.311203, sdsep.calculateLogP(), 1e-1); // Used in original version with fixed-step size ODE solver
	}

}
