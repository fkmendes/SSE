package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

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
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2");
		traitStash.printLksMap();
				
		// initializing birth-death parameters
		Double birthRate = 0.222222222;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols
		
		Double[] lambdas = new Double[numberOfStates];
		Arrays.fill(lambdas, birthRate);
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
    	
    	System.out.println(sdsep.calculateLogP());
	}

	private void assert2DArrayEquals(double[][] arr1, double[][] arr2) {
		Assert.assertEquals(arr1.length, arr2.length);

		for (int i = 0; i < arr1.length; i++) {
			for (int j = 0; j < arr1[i].length; j++) {
				Assert.assertEquals(arr1[i][j], arr2[i][j], 1e-7);
			}
		}
	}

	private double[][] deep2DArrayCopy(double[][] arr) {
	    double[][] ret = new double[arr.length][arr[0].length];

	    for (int i = 0; i < arr.length; i++) {
	        System.arraycopy(arr[i], 0, ret[i], 0, arr[i].length);
        }
        return ret;
	}

	@Test
	public void test() {
		// Assert.assertEquals(-5.588460032653, sdsep.calculateLogP(), EPSILON); // Used in original version with fixed-step size ODE solver
		Assert.assertEquals(-5.5884600307, sdsep.calculateLogP(), EPSILON);
		System.out.println("Passed likelihood test!");
		double[][] nodePartialsNoCharHist = deep2DArrayCopy(sdsep.getNodePartialScaledLksPostOde());

		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-5.5884600307, sdsep.calculateLogP(), 1e-6);
		System.out.println("Passed likelihood with sample character history test!");
		double[][] nodePartialsWCharHist = sdsep.getNodePartialScaledLksPostOde();

		assert2DArrayEquals(nodePartialsNoCharHist, nodePartialsWCharHist);
		System.out.println("Passed matching likelihood test!");
	}

}
