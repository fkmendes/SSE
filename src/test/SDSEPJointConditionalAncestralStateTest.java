package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

public class SDSEPJointConditionalAncestralStateTest {
	final static double EPSILON = 1e-10;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE
		String[] spNames = new String[] { "sp1", "sp2", "sp3" };
		int numSpecies = spNames.length;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "sp1=1,sp2=1,sp3=1");
		traitStash.printLksMap();

		// initializing birth-death parameters
		Double[] mus = { 0.001, 0.1 };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols

		Double[] lambdas = {0.2, 0.4};
		RealParameter lambda = new RealParameter(lambdas);

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.1 0.4");
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

		String treeStr = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
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

		// Run the sampling many times
        int numTrials = 10000;
		int[][] samples = new int[numTrials][2 * numSpecies - 1];
		for (int i = 0; i < numTrials; i++) {
			sdsep.drawJointConditionalAncestralStates();
			int[] drawnAncestralEnd = sdsep.endStates;
			System.arraycopy(drawnAncestralEnd, 0, samples[i], 0, 2 * numSpecies - 1);
		}
		System.out.println(Arrays.toString(samples[0]));
		System.out.println(Arrays.toString(samples[8]));

		// Calculate the posterior probabilities by counting the frequency the node is in state one
		double[] posterior = new double[2 * numSpecies - 1];
		int numStateOne;
		for (int nIdx = 0; nIdx < 2 * numSpecies - 1; nIdx++) {
			numStateOne = 0;
			for (int nTrial = 0; nTrial < numTrials; nTrial++) {
				if (samples[nTrial][nIdx] == 1) {
					numStateOne++;
				}
			}
			System.out.println(numStateOne);
			posterior[nIdx] = 1.0 * numStateOne / numTrials;
        }
        System.out.println("Summary of sampled states: " + Arrays.toString(posterior));

		// Write only the ancestral states to csv
		BufferedWriter br = new BufferedWriter(new FileWriter("beast.csv"));
		StringBuilder sb = new StringBuilder();

		for (int i = numSpecies; i < 2 * numSpecies - 1; i++) {
			double element = posterior[i];
			sb.append(Double.toString(element));
			sb.append(",");
		}

		br.write(sb.toString());
		br.close();
	}

	@Test
	public void test() {
		// Assert.assertEquals(-5.588460032653, sdsep.calculateLogP(), EPSILON); // Used in original version with fixed-step size ODE solver
	}

}
