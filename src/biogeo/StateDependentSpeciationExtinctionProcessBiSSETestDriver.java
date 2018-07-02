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
		int num_states = 2; // BiSSE
		String[] sp_names = new String[] { "Human", "Chimp", "Gorilla" };
		List<Taxon> taxa_list = Taxon.createTaxonList(Arrays.asList(sp_names));
		TaxonSet taxon_set = new TaxonSet(taxa_list);
		TraitStash trait_stash = new TraitStash(num_states);
		trait_stash.initByName("taxa", taxon_set, "value", "Human=2,Chimp=2,Gorilla=2");
		trait_stash.printLksMap();
		
//		List<TraitStash> trait_stashes = new ArrayList<TraitStash>();
//		for (int h = 1; h<=2; ++h) {
//			for (int c = 1; c<=2; ++c) {
//				for (int g = 1; g<=2; ++g) {
//					String trait_str = "Human=" + Integer.toString(h) +
//							",Chimp=" + Integer.toString(c) +
//							",Gorilla=" + Integer.toString(g);
//					TraitStash this_stash = new TraitStash(num_states);
//					this_stash.initByName("taxa", taxon_set, "value", trait_str);
//					trait_stashes.add(this_stash);
//				}
//			}
//		}
		
		double rate = 1.0;
		double[] lambda = new double[] {0.222222222, 0.222222222};
		double[] mu = new double[] {0.1, 0.1};
		
		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
		Q.setCell(0, 0, 0.9); // q11
		Q.setCell(1, 1, 0.9); // q22
		Q.setCell(0, 1, 0.1); // q12
		Q.setCell(1, 0, 0.1); // q21
		Q.printMatrix();
		
		// not using it
		int[][] cladogenetic_events = {{2, 2, 1},{1, 1, 1}};
		CladogeneticSpeciationRateStash clado_stash = new CladogeneticSpeciationRateStash(cladogenetic_events, lambda);
		
		double[] pi_es = new double[num_states];
		double[] pi_ds = new double[num_states];
		Arrays.fill(pi_ds, 1.0/(double)num_states);
		double[] pi = ArrayUtils.addAll(pi_es, pi_ds); // 0.0, 0.0, 0.5, 0.5
		
//		double[] pi = new double[num_states*2];
//		Arrays.fill(pi, 1.0/((double)num_states*2)); // 0.25, 0.25, 0.25, 0.25
		
		System.out.println("Pi is: " + Arrays.toString(pi));
		
		String tree_str = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
        TreeParser my_tree = new TreeParser(tree_str, false, false, true, 0); // true b/c species are labelled, offset=0
        Node root = my_tree.getRoot();
        // Node[] nodes = my_tree.getNodesAsArray();
        
        boolean incorporate_cladogenesis = false;
        
    	StateDependentSpeciationExtinctionProcess sdsep = 
		new StateDependentSpeciationExtinctionProcess(my_tree, lambda, mu, pi, num_states,
				trait_stash, clado_stash, Q, rate, incorporate_cladogenesis);
    	
    	sdsep.computeNodeLk(root, root.getNr());
        
//        double[] res = new double[trait_stashes.size()];
//        double total_likelihood = 0.0;
//        for (int i = 0; i < trait_stashes.size(); ++i) {
//        	System.out.println(trait_stashes.get(i).getStringValue("Human"));
//        	System.out.println(trait_stashes.get(i).getStringValue("Chimp"));
//        	System.out.println(trait_stashes.get(i).getStringValue("Gorilla"));
//        	StateDependentSpeciationExtinctionProcess sdsep = 
//            		new StateDependentSpeciationExtinctionProcess(my_tree, lambda, mu, pi, num_states,
//            				trait_stashes.get(i), clado_stash, Q, rate, incorporate_cladogenesis);
//            
//            sdsep.computeNodeLk(root, root.getNr());
//            res[i] = sdsep.getProb();     
//            total_likelihood += res[i];
//
//        }
//        System.out.println(Arrays.toString(res));
//        System.out.println(total_likelihood);
	}
}
