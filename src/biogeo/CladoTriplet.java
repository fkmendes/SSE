package biogeo;

import beast.core.BEASTObject;
import beast.core.Input;

public class CladoTriplet extends BEASTObject {
		
	final public Input<Integer> parentStateInput = new Input<>("ParentState", "Character state of parent.");
	final public Input<Integer> leftChildInput = new Input<>("LeftChildState", "Character state of left child.");
	final public Input<Integer> rightChildInput = new Input<>("RightChildState", "Character state of right child.");	
	final public Input<speciationType> speciationTypeInput = new Input<>("SpeciationType", "Speciation type (sympatric, subsympatric, vicariance).");
	
	public enum speciationType { SYMPATRY, SUBSYMPATRY, VICARIANCE };
	private int[] cladogeneticEvent;
	private speciationType speciationEvent;
	
	@Override
	public void initAndValidate() {
		cladogeneticEvent = new int[3];
		cladogeneticEvent[0] = parentStateInput.get();
		cladogeneticEvent[1] = leftChildInput.get();
		cladogeneticEvent[2] = rightChildInput.get();
		speciationEvent = speciationTypeInput.get();
	}
	
	public int[] getCladogeneticEvent() {
		return cladogeneticEvent;
	}
	
	public speciationType getSpeciationEvent() {
		return speciationEvent;
	}
}
