package SSE;

import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;

public class CIDMaskOperator extends Operator {
    final public Input<IntegerParameter> cidMaskInput =
            new Input<>("mask", "integer (1=CID, 0=Not-CID) to operate on.", Validate.REQUIRED);
    
	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final IntegerParameter mask = cidMaskInput.get();
		final double value = mask.getValue(0);
		
		// if at this position is 0 or 2, we can only (and will!) go to 1
		if (value == 0) {
			mask.setValue(0, 1);
		} else {
			mask.setValue(0, 0);
		}
		
		return 0.0;
	}
}
