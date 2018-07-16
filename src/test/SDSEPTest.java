package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import biogeo.CladogeneticSpeciationRateStash;
import biogeo.InstantaneousRateMatrix;
import biogeo.StateDependentSpeciationExtinctionProcess;
import biogeo.TraitStash;

public class SDSEPTest {
	final static double EPSILON = 1e-10;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
				int num_states = 4; // state = 1 (index 0) is 'null state'
				String[] sp_names = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
				List<Taxon> taxa_list = Taxon.createTaxonList(Arrays.asList(sp_names));
				TaxonSet taxon_set = new TaxonSet(taxa_list);
				TraitStash trait_stash = new TraitStash(num_states);
				trait_stash.initByName("taxa", taxon_set, "value", "Human=2,Chimp=2,Gorilla=2,Orang=3");
				trait_stash.printLksMap();
						
				// initializing birth-death parameters
				double rate = 1.0;
				double birth_rate = 0.32222224;
				double death_rate = 0.1; // DEC-like	
				double symp_prob = 1.0; // DEC-like
				double j_prob = 0.0; // no jump dispersal	
				double[] mu = new double[] {death_rate, death_rate, death_rate, death_rate};
				double[] lambdas_clado_stash = new double[] {
						birth_rate, 
						symp_prob*birth_rate, 
						symp_prob*birth_rate,
						j_prob*birth_rate, 
						j_prob*birth_rate,
						(double)1/6*birth_rate, 
						(double)1/6*birth_rate, 
						(double)1/6*birth_rate}; // the "rotated" version of these ones are unidentifiable and set to 0 (see Goldberg and Igic 2012)
				int[][] cladogenetic_events = {
						{1,1,1},
						{2,2,2},
						{3,3,3},
						{2,2,3}, 
						{3,2,3},
						{4,2,3},
						{4,2,4},
						{4,3,4}};
				CladogeneticSpeciationRateStash clado_stash = new CladogeneticSpeciationRateStash(cladogenetic_events, lambdas_clado_stash);
				clado_stash.printEventMap();
				double[] lambda = new double[num_states];
				Arrays.fill(lambda, birth_rate);

				// initializing trait parameters
				InstantaneousRateMatrix Q = new InstantaneousRateMatrix(num_states);
				Q.setCell(0, 1, 0.00); // q12 
				Q.setCell(0, 2, 0.00); // q13 
				Q.setCell(0, 3, 0.00); // q14 
				Q.setCell(1, 0, 0.01); // q21
				Q.setCell(1, 2, 0.00); // q23
				Q.setCell(1, 3, 0.01); // q24
				Q.setCell(2, 0, 0.01); // q31
				Q.setCell(2, 1, 0.00); // q32
				Q.setCell(2, 3, 0.01); // q34
				Q.setCell(3, 0, 0.00); // q41
				Q.setCell(3, 1, 0.01); // q42
				Q.setCell(3, 2, 0.01); // q43
						
				// initializing prior
				double[] pi_es = new double[num_states];
				double[] pi_ds = new double[num_states];
				Arrays.fill(pi_ds, 1.0/((double)(num_states)));
				double[] pi = ArrayUtils.addAll(pi_es, pi_ds); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
						
				// initializing tree
				String tree_str = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
				TreeParser my_tree = new TreeParser(tree_str, false, false, true, 0); // true b/c species are labelled, offset=0
				Node root = my_tree.getRoot();
				        
				// we do want cladogenesis
				boolean incorporate_cladogenesis = true;
				        
				// finally initialize process and compute lk
				sdsep = new StateDependentSpeciationExtinctionProcess(my_tree, lambda, mu, pi, num_states,
						trait_stash, clado_stash, Q, rate, incorporate_cladogenesis);
				sdsep.computeNodeLk(root, root.getNr());
	}

	@Test
	public void test() {
		Assert.assertEquals(-10.59346884351, sdsep.getLogLk(), EPSILON);
	}

}