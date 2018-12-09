package SSE;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.Input.Validate;

public class MaskOperatorDist extends Distribution {

	final public Input<IntegerParameter> stateMaskInput = new Input<>("stateMask", "State part of mask.", Validate.REQUIRED);
	final public Input<IntegerParameter> cidMaskInput = new Input<>("cidMask", "CID part of mask.", Validate.REQUIRED);
	
	IntegerParameter stateMask;
	IntegerParameter cidMask;
	
	@Override
	public void initAndValidate() {
		stateMask = stateMaskInput.get();
		cidMask = cidMaskInput.get();
	}
	
	@Override
	public double calculateLogP() {
		return 1.0;
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
