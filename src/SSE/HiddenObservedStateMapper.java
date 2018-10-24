package SSE;

import java.util.HashMap;

import beast.core.BEASTObject;
import beast.core.Input;

public class HiddenObservedStateMapper extends BEASTObject {

	final public Input<String> hiddenStatesInFocalOrderInput = new Input<>("hiddenStates", "Comma-separated integer strings corresponding to hidden states, one per focal states (in the order of focal states).");
	private HashMap<Integer, Integer> obs2HiddenMap = new HashMap<>(); // makeMaps populates this
	private HashMap<Integer, Integer[]> hidden2ObsMap = new HashMap<>(); // makeMaps populates this
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
	}

	public void makeMaps() {
		String hiddenStateString = hiddenStatesInFocalOrderInput.get();
		String[] hiddenStateStrings = hiddenStateString.split(",");
	
		for (int i=0; i<hiddenStateStrings.length; i++) {
			obs2HiddenMap.put(i, Integer.parseInt(hiddenStateStrings[i]));
		}
		
		printObs2HiddenMap();
		
		// TODO Create hidden2ObsMap, then pass this class instance to the HiddenIRM class so that it can disallow double transitions
	}
	
	// helper
	public void printObs2HiddenMap() {
		System.out.println("OBS HIDDEN");
		for (Integer obs: obs2HiddenMap.keySet()) {
			// PRINT MAP
			String obsString = Integer.toString(obs);
			String hiddenString = obs2HiddenMap.get(obs).toString();
			System.out.println(obsString + " " + hiddenString);
		}
	}
}
