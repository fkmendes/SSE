package biogeo;

import java.util.Arrays;
import java.util.HashMap;

public class StateStash {

	int num_states;
	String[] sp_names;
	private HashMap<String, String> spname_state_map = new HashMap<String, String>(); // set by user
	private HashMap<String, double[]> spname_lks_map = new HashMap<String, double[]>(); // initialized by ctor
	
	public StateStash(int num_states, String[] sp_names) {
		this.num_states = num_states;
		this.sp_names = sp_names;
		
		for (int i = 0; i < sp_names.length; ++i) {
			double[] lks = new double[num_states * 2]; // Es and Ds, so num_states * 2
			spname_lks_map.put(sp_names[i], lks); // initializing the state lks of each species to 0.0
		}
	}

	// setters and getters
	public void setSpStateValue(String spname, String state) {
		spname_state_map.put(spname, state);
		int sp_state_idx = num_states + Integer.parseInt(state); // we add num_states because of Es
		spname_lks_map.get(spname)[sp_state_idx] = 1.0;
	}
	
	public HashMap<String, double[]> getSpLkMap() {
		return spname_lks_map;
	}
	
	public double[] getSpLk(String spname) {
		return spname_lks_map.get(spname);
	}
	
	public HashMap<String, String> getSpStateMap() {
		return spname_state_map;
	}
	
	public String getSpState(String spname) {
		return spname_state_map.get(spname);
	}
	
	// helper
	public void printLksMap() {
		for (HashMap.Entry<String, double[]> entry : spname_lks_map.entrySet()) {	
			String spname = entry.getKey();
			double[] lks = entry.getValue();
			System.out.println(spname + " initial lks = " + Arrays.toString(lks));
		}
	}
}

