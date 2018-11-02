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
	final public Input<HiddenObservedStateMapper> HiddenObservedStateMapperInput = new Input<>("hiddenObsStateMapper", "Maps hidden states onto observed states and vice-versa.", Validate.REQUIRED);
	
	// if Human=1,Chimp=1,Gorilla=1
	// inheriting taxonValues, map, and values variables
	// protected String[] taxonValues; // state values as str, e.g., ["1", "1", "1"]
	// Map<String, Integer> map; // "spname":sp index, e.g. {"Human":0, "Chimp":1, "Gorilla":2}
	// double[] values; // state values as doubles, e.g., [1.0, 1.0, 1.0]
	private HashMap<String, double[]> spNameLksMap = new HashMap<String, double[]>(); // initialized by ctor
	private int numberOfStates;
	private int numberOfHiddenStates;
	private int totalNumberOfStates;
	protected Map<String, Integer> hiddenMap;
	private HiddenObservedStateMapper hiddenObsStateMapper;
	
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
        String[] traits = traitsInput.get().split(","); // ["sp1=1", "sp2=1", "sp3=1"]
        taxonValues = new String[labels.size()]; // one observed state string per taxon
        values = new double[labels.size()];
        
        if (numberOfHiddenStates > 0) { hiddenObsStateMapper = HiddenObservedStateMapperInput.get(); }
        
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
            
            taxonValues[taxonNr] = normalize(strs[1]);
            
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
        int hiddenIdx;
        for (Entry<String, Integer> entry : map.entrySet()) {
        	double[] lks = new double[totalNumberOfStates*2];
			String spName = entry.getKey();
			int spIdx = entry.getValue();
			spNameLksMap.put(spName, lks);
			// System.out.println(taxonValues[sp_idx]);
			spNameLksMap.get(spName)[totalNumberOfStates - 1 + Integer.parseInt(taxonValues[spIdx])] = 1.0;
			
			if (numberOfHiddenStates > 0) {
				hiddenIdx = hiddenObsStateMapper.getHiddenFromObs(Integer.parseInt(taxonValues[spIdx])-1);
				
				// only if hiddenIdx is matched with some observed state
				if (hiddenIdx != -1) {
					spNameLksMap.get(spName)[totalNumberOfStates + numberOfStates + hiddenIdx] = 1.0;
				}
			} // hidden states match observed states according to hiddenObsStateMapper
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
