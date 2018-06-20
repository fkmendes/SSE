package biogeo;

import beast.evolution.tree.*;
import beast.util.TreeParser;

public class StateDependentSpeciationExtinctionProcessTestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		String[] sp_names = new String[] { "Human", "Chimp", "Gorilla" };
		int num_states = 2; // BiSSE
		double[] lambda = new double[] {0.222222222, 0.222222222};
		double[] mu = new double[] {0.0, 0.0}; // pure birth
				
		InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
		Q.setCell(0, 0, 1.0);
		Q.setCell(1, 1, 1.0);
		Q.setCell(0, 1, 0.0);
		Q.setCell(1, 0, 0.0);
		
		StateStash state_stash = new StateStash(num_states, sp_names);
		state_stash.setSpStateValue("Human", "1");
		state_stash.setSpStateValue("Chimp", "1");
		state_stash.setSpStateValue("Gorilla", "1");
		state_stash.printLksMap();
		
		String tree_str = "((Human:1.0,Chimp:1.0):1.0,Gorilla:1.0)0.0;";
        TreeParser my_tree = new TreeParser(tree_str, false, false, true, 0); // true b/c species are labelled, offset=0
        Node root = my_tree.getRoot();
        int num_nodes = my_tree.getNodeCount();
        System.out.println(num_nodes);
        Node[] nodes = my_tree.getNodesAsArray();
	}
}
