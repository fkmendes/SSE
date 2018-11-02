package SSE;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.*;


import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.BEASTObject;
import beast.core.Input;

public class HiddenObservedStateMapper extends BEASTObject {

	final public Input<String> hiddenStatesInFocalOrderInput = new Input<>("hiddenStates", "Comma-separated integer strings corresponding to hidden states, one per focal states (in the order of focal states).");
	private Map<Integer, Integer> obs2HiddenMap = new HashMap<>(); // makeMaps populates this
	private Multimap<Integer, Integer> hidden2ObsMap = HashMultimap.create();
	// private HashMap<Integer, Integer[]> hidden2ObsMap = new HashMap<>(); // makeMaps populates this
	
	@Override
	public void initAndValidate() {
		makeMaps();
	}

	public void makeMaps() {
		String hiddenStateString = hiddenStatesInFocalOrderInput.get();
		String[] hiddenStateStrings = hiddenStateString.split(",");
	
		for (int i=0; i<hiddenStateStrings.length; i++) {
			obs2HiddenMap.put(i, Integer.parseInt(hiddenStateStrings[i]));
		}
		
		printObs2HiddenMap(); // for debugging
		
		for (Entry<Integer, Integer> entry: obs2HiddenMap.entrySet()) {
			hidden2ObsMap.put(entry.getValue(), entry.getKey());
		}
		
		printHidden2ObsMap(); // for debugging

	}
	
	// helper
	public int getHiddenFromObs(int obsIdx) {
		return obs2HiddenMap.get(obsIdx);
	}
	
	public Collection<Integer> getObsCollFromHidden(int hiddenIdx) {
		return hidden2ObsMap.get(hiddenIdx);
	}
	
	public void printObs2HiddenMap() {
		for (Integer obs: obs2HiddenMap.keySet()) {
			String obsString = Integer.toString(obs);
			String hiddenString = obs2HiddenMap.get(obs).toString();
			System.out.println("Obs=" + obsString + " hidden=" + hiddenString);
		}
	}
	
	public void printHidden2ObsMap() {
		for (Entry<Integer, Collection<Integer>> entry: hidden2ObsMap.asMap().entrySet()) {
			System.out.println("hidden=" + entry.getKey() + " obs=" + entry.getValue());
		}
	}
}
