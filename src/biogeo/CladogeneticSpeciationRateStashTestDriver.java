package biogeo;

public class CladogeneticSpeciationRateStashTestDriver {

	public static void main(String[] args) {
		String[][] cladogenetic_events = {{"12","12","1"},{"1","1","1"}};
		double[] lambda = {0.1, 0.2};
		CladogeneticSpeciationRateStash clado_stash = new CladogeneticSpeciationRateStash(cladogenetic_events, lambda);
		clado_stash.printEventMap();		
	}
}

