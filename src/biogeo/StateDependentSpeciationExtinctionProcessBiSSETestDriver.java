package biogeo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.*;
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
		
		double rate = 1.0;
		double[] lambda = new double[] {0.222222222, 0.222222222};
		double[] mu = new double[] {0.1, 0.1};
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("NumberOfStates", numStates, "FlatQMatrix", "0.9 0.1 0.9 0.1");
		irm.printMatrix();
		
		double[] piEs = new double[numStates];
		double[] piDs = new double[numStates];
		Arrays.fill(piDs, 1.0/(double)numStates);
		double[] pi = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5

		System.out.println("Pi is: " + Arrays.toString(pi));
		
		String treeStr = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        Node root = myTree.getRoot();
        // Node[] nodes = my_tree.getNodesAsArray();
        
        boolean incorporateCladogenesis = false;
        
        CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash(); // not using it
        
    	StateDependentSpeciationExtinctionProcess sdsep = 
		new StateDependentSpeciationExtinctionProcess(myTree, lambda, mu, pi, numStates,
				traitStash, csrt, irm, rate, incorporateCladogenesis);
    	
    	sdsep.computeNodeLk(root, root.getNr());
	}
}
