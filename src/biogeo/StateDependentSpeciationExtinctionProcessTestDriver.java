package biogeo;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.*;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcessTestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		String[] sp_names = new String[] { "Human", "Chimp", "Gorilla" };
		List<Taxon> taxa_list = Taxon.createTaxonList(Arrays.asList(sp_names));
		TaxonSet taxon_set = new TaxonSet(taxa_list);
		TraitStash trait_stash = new TraitStash();
		trait_stash.initByName("taxa", taxon_set, "value", "Human=1,Chimp=2,Gorilla=2");
		trait_stash.printLksMap();
		
		int num_states = 2; // BiSSE
		double rate = 1.0;
		double[] lambda = new double[] {0.222222222, 0.222222222};
		double[] mu = new double[] {0.0, 0.0}; // pure birth
		String[][] cladogenetic_events = {{"12","12","1"},{"1","1","1"}};
		
		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
		Q.setCell(0, 0, 1.0);
		Q.setCell(1, 1, 1.0);
		Q.setCell(0, 1, 0.0);
		Q.setCell(1, 0, 0.0);
		
		CladogeneticSpeciationRateStash clado_stash = new CladogeneticSpeciationRateStash(cladogenetic_events, lambda);
		
		double[] pi_es = new double[num_states];
		double[] pi_ds = new double[num_states];
		Arrays.fill(pi_ds, 1.0/(double)num_states);
		double[] pi = ArrayUtils.addAll(pi_es, pi_ds);
		System.out.println(Arrays.toString(pi));
		
//		StateStash state_stash = new StateStash(num_states, sp_names);
//		state_stash.setSpStateValue("Human", "1");
//		state_stash.setSpStateValue("Chimp", "1");
//		state_stash.setSpStateValue("Gorilla", "1");
//		state_stash.printLksMap();
		
		String tree_str = "((Human:1.0,Chimp:1.0):1.0,Gorilla:1.0)0.0;";
        TreeParser my_tree = new TreeParser(tree_str, false, false, true, 0); // true b/c species are labelled, offset=0
        Node root = my_tree.getRoot();
        // Node[] nodes = my_tree.getNodesAsArray();
        
        boolean incorporate_cladogenesis = false;
        
        StateDependentSpeciationExtinctionProcess sdsep = 
        		new StateDependentSpeciationExtinctionProcess(my_tree, lambda, mu, num_states,
        				trait_stash, clado_stash, Q, rate, incorporate_cladogenesis);
        
        sdsep.computeNodeLk(root, root.getNr());
	}
}
