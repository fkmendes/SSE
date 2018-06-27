package biogeo;

import java.util.HashMap;

public class CladogeneticSpeciationRateStash {
	
	private int[][] cladogenetic_events;
	private double[] speciation_rates;
	private HashMap<int[], Double> event_map = new HashMap<int[], Double>(); // ctor populates this
	
	// ctor (populates event_map)
	public CladogeneticSpeciationRateStash(int[][] cladogenetic_events, double[] speciation_rates) {
		this.cladogenetic_events = cladogenetic_events;
		this.speciation_rates = speciation_rates;
		int num_events = cladogenetic_events.length;
		
		if (num_events != this.speciation_rates.length) {
			System.out.println("The number of cladogenetic events did not match the number of speciation rates. Exiting...");
			throw new RuntimeException();
		}
		
		for (int i = 0; i < num_events; ++i) {
			if (this.cladogenetic_events[i].length != 3) {
				System.out.println("We need cladogenetic events as parent-daughter1-daughter2. Found something different. Exiting...");
			}
			
			int[] key = new int[3];
			key[0] = this.cladogenetic_events[i][0];
			key[1] = this.cladogenetic_events[i][1];
			key[2] = this.cladogenetic_events[i][2];
			this.event_map.put(key, this.speciation_rates[i]);
		}
	}
	
	// setters and getters
	HashMap<int[], Double> getEventMap() {
		return event_map;
	}
	
	// helper
	public void printEventMap() {
		System.out.println(new PrettyPrintHashMap<int[], Double>(event_map));
	}
}
