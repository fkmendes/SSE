package SSE;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.tree.*;

public class HiddenTraitStash extends TraitSet {

	final public Input<Integer> nStatesInput = new Input<>("numberOfStates", "How many observed states or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> nHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<Boolean> passingHidden = new Input<>("passingHidden", "Whether both and hidden states are being passed", Validate.REQUIRED);
	
	// if Human=1,Chimp=1,Gorilla=1
	// inheriting taxonValues, map, and values variables
	// protected String[] taxonValues; // state values as str, e.g., ["1", "1", "1"]
	// Map<String, Integer> map; // "spname":sp index, e.g. {"Human":0, "Chimp":1, "Gorilla":2}
	// double[] values; // state values as doubles, e.g., [1.0, 1.0, 1.0]
	private HashMap<String, double[]> spNameLksMap = new HashMap<String, double[]>(); // initialized by ctor
	private int numberOfStates;
	private int numberOfHiddenStates;
	private int totalNumberOfStates;
	private String[] taxonHiddenValues;
	protected Map<String, Integer> hiddenMap;
	
	public HiddenTraitStash() {
		traitNameInput.setRule(Input.Validate.FORBIDDEN);
		dateTimeFormatInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	public void initAndValidate() {
		map = new HashMap<>();
		numberOfStates = nStatesInput.get();
		numberOfHiddenStates = nHiddenStatesInput.get();
		totalNumberOfStates = numberOfStates + numberOfHiddenStates;
        List<String> labels = taxaInput.get().asStringList();
        String[] traits = traitsInput.get().split(","); // ["Human=1", "Chimp=1", "Gorilla=1"]
        boolean hiddenWasPassed = passingHidden.get();
        taxonValues = new String[labels.size()]; // one observed state string per taxon
        taxonHiddenValues = new String[labels.size()];
        values = new double[labels.size()];
        
        for (String trait : traits) {
            trait = trait.replaceAll("\\s+", " ");
            String[] strs = trait.split("=");
            
            if (strs.length != 2) {
                throw new IllegalArgumentException("Could not parse trait: " + trait);
            }
                        
            String taxonID = normalize(strs[0]); // taxon name
            int taxonNr = labels.indexOf(taxonID); // taxon idx
            
            if (taxonNr < 0) {
                throw new IllegalArgumentException("Trait (" + taxonID + ") is not a known taxon. Spelling error perhaps?");
            }
            
            if (!hiddenWasPassed) {
            	taxonValues[taxonNr] = normalize(strs[1]);
            } // only observed passed
            else {
            	String[] obsAndHiddenStrs = strs[1].split("|");            	
            	taxonValues[taxonNr] = normalize(obsAndHiddenStrs[0]);
            	taxonHiddenValues[taxonNr] = normalize(obsAndHiddenStrs[2]);
            } // "observed,hidden" was passed (only used for validation)
            
            // converting state to double
            try {
            	values[taxonNr] = convertValueToDouble(taxonValues[taxonNr]);
            } catch (IllegalArgumentException ex) {
                Log.err.println("Failed to parse trait '" + taxonValues[taxonNr] + "'.");
                System.exit(1);
            }
            
            map.put(taxonID, taxonNr); // {taxon name: taxon idx}
        }
        
        // num_states = numOfValues(taxonValues);
        
        // find extremes
        minValue = values[0];
        maxValue = values[0];
        for (double value : values) {
            minValue = Math.min(minValue, value);
            maxValue = Math.max(maxValue, value);
        }
        
        // minValue should be 1
        if (minValue == 0) {
        	Log.warning.println("WARNING: Trait state values should start at 1, not 0. Exiting...");
        	System.exit(1);
        }
        
        // if all species share the same state, we assume binary state
        // developers: this is useful for testing likelihood computations
        // if (numberOfStates == 1) {
        //    Log.warning.println("WARNING: All species had same trait state. Assuming binary trait.");
        //    numberOfStates = 2;
        // }
        
        // checking all distinct states are consecutive order (1,2,4 is not allowed)
         if ((double) numberOfStates < maxValue) {
             Log.warning.println("WARNING: The larger observed state identifier (number) is greater than the total count of observed states (i.e., you specified 3 states, but then one species is =4). Exiting...");
             System.exit(1);
         }
        
        // checking all species were specified some trait value
        for (int i = 0; i < labels.size(); i++) {
            if (taxonValues[i] == null) {
            	Log.warning.println("WARNING: No trait specified for " + labels.get(i) + ". Exiting...");
            	System.exit(1);
            }
        }
        
        // populating spNameLksMap
        for (Entry<String, Integer> entry : map.entrySet()) {
        	double[] lks = new double[totalNumberOfStates*2];
			String spName = entry.getKey();
			int spIdx = entry.getValue();
			spNameLksMap.put(spName, lks);
			// System.out.println(taxonValues[sp_idx]);
			spNameLksMap.get(spName)[totalNumberOfStates - 1 + Integer.parseInt(taxonValues[spIdx])] = 1.0;
			
			if (!hiddenWasPassed) {
				spNameLksMap.get(spName)[totalNumberOfStates + numberOfStates - 1 + Integer.parseInt(taxonValues[spIdx])] = 1.0;
//				for (int h=0; h<numberOfHiddenStates; ++h) {
//					spNameLksMap.get(spName)[totalNumberOfStates + numberOfStates + h] = 1.0; // skip E's, then skip observed D's
//				}
			} // all hidden states are initialized to 1 if hidden states not passed
			else {
				spNameLksMap.get(spName)[totalNumberOfStates + numberOfStates - 1 + Integer.parseInt(taxonHiddenValues[spIdx])] = 1.0;
			} // hidden states are passed and initialized like observed states (only used for validation)
        }
	}
	
	// getters and setters
	public double[] getSpLks(String spname) {
		return spNameLksMap.get(spname);
	}
	
	// helper
	public int numOfValues(String[] arr) {
	    return (int) Arrays.stream(arr).distinct().count();
	}
	
	public void printLksMap() {
		for (HashMap.Entry<String, double[]> entry : spNameLksMap.entrySet()) {	
			String spname = entry.getKey();
			double[] lks = entry.getValue();
			System.out.println(spname + " initial lks = " + Arrays.toString(lks));
		}
	}
}
