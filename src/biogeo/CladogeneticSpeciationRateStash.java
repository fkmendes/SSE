package biogeo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import biogeo.CladoTriplet.speciationType;

public class CladogeneticSpeciationRateStash extends CalculationNode {
	
	final public Input<List<CladoTriplet>> cladogeneticEventsInput = new Input<>("CladoTriplets", "List of CladoTriplet objects that contain states of a parent and its children.", new ArrayList<>());
	final public Input<RealParameter> sympatricRateInput = new Input<>("SympatricRate", "Speciation rate for sympatric event.");
	final public Input<RealParameter> subSympatricRateInput = new Input<>("SubsympatricRate", "Speciation rate for subsympatric event.");
	final public Input<RealParameter> vicariantRateInput = new Input<>("VicariantRate", "Speciation rate for vicariant event.");
	final public Input<RealParameter> jumpRateInput = new Input<>("JumpRate", "Speciation rate for jump dispersal event.");
	
//	private int[][] cladogeneticEvents;
	private double sympatricRate;
	private double subSympatricRate;
	private double vicariantRate;
	private double jumpRate;
	private HashMap<int[], Double> eventMap = new HashMap<>(); // ctor populates this
	private boolean stashClean = false;
	
	@Override
	public void initAndValidate() {
		populateStash();
	}

	public void populateStash() {

		//Parameter sympatricRateParam = sympatricRateInput.get();

		sympatricRate = sympatricRateInput.get().getValue();
		subSympatricRate = subSympatricRateInput.get().getValue();
		vicariantRate = vicariantRateInput.get().getValue();
		jumpRate = jumpRateInput.get().getValue();

		List<CladoTriplet> cladoTripletList = cladogeneticEventsInput.get();

		for (CladoTriplet cladoTriplet : cladoTripletList) {
			int[] cladogeneticEvent = cladoTriplet.getCladogeneticEvent();
			speciationType speciationEvent = cladoTriplet.getSpeciationEvent();

			switch (speciationEvent) {
				case SYMPATRY:
					eventMap.put(cladogeneticEvent, sympatricRate);
					break;
				case SUBSYMPATRY:
					eventMap.put(cladogeneticEvent, subSympatricRate);
					break;
				case VICARIANCE:
					eventMap.put(cladogeneticEvent, vicariantRate);
					break;
				case JUMPDISPERSAL:
					eventMap.put(cladogeneticEvent, jumpRate);
					break;
			}
		}
		stashClean = true;
	}

	protected boolean requiresRecalculation() {
		stashClean = !(sympatricRateInput.get().somethingIsDirty() || subSympatricRateInput.get().somethingIsDirty() ||
				vicariantRateInput.get().somethingIsDirty() || jumpRateInput.get().somethingIsDirty());
		return stashClean;
	}


//	// ctor (populates event_map)
//	public CladogeneticSpeciationRateStash(int[][] cladogenetic_events, double[] speciation_rates) {
//		this.cladogeneticEvents = cladogenetic_events;
//		this.speciationRates = speciation_rates;
//		int num_events = cladogenetic_events.length;
//		
//		if (num_events != this.speciationRates.length) {
//			System.out.println("The number of cladogenetic events did not match the number of speciation rates. Exiting...");
//			throw new RuntimeException();
//		}
//		
//		for (int i = 0; i < num_events; ++i) {
//			if (this.cladogeneticEvents[i].length != 3) {
//				System.out.println("We need cladogenetic events as parent-daughter1-daughter2. Found something different. Exiting...");
//			}
//			
//			int[] key = new int[3];
//			key[0] = this.cladogeneticEvents[i][0];
//			key[1] = this.cladogeneticEvents[i][1];
//			key[2] = this.cladogeneticEvents[i][2];
//			this.eventMap.put(key, this.speciationRates[i]);
//		}
//	}
	
	// setters and getters
	HashMap<int[], Double> getEventMap() {
		if (!stashClean) {
			populateStash();
		}
		return eventMap;
	}
	
	// helper
	public void printEventMap() {
		System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
	}
}
