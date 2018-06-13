package biogeo;

import java.util.HashMap;

public class CladogeneticSpeciationRateStash {
	
	String[][] cladogenetic_events;
	double[] speciation_rates;
	HashMap<String[], Double> event_map = new HashMap<String[], Double>(); // ctor populates this
	
	// ctor (populates event_map)
	public CladogeneticSpeciationRateStash(String[][] cladogenetic_events, double[] speciation_rates) {
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
			
			String[] key = new String[3];
			key[0] = this.cladogenetic_events[i][0];
			key[1] = this.cladogenetic_events[i][1];
			key[2] = this.cladogenetic_events[i][2];
			this.event_map.put(key, this.speciation_rates[i]);
		}
	}
	
	// setters and getters
	HashMap<String[], Double> getEventMap() {
		return event_map;
	}
	
	// helper
	public void printEventMap() {
		System.out.println(new PrettyPrintHashMap<String[], Double>(event_map));
	}
}
