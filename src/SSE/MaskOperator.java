package SSE;

import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;

public class MaskOperator extends Operator {
    final public Input<IntegerParameter> stateMaskInput =
            new Input<>("mask", "the array of doubles to operate a random walk on.", Validate.REQUIRED);
    
	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final IntegerParameter mask = stateMaskInput.get();
		final int pos = Randomizer.nextInt(mask.getDimension());
		final double value = mask.getValue(pos);
		
		// no matter if CID or hidden state flag
		if (value == 0.0) {
			mask.setValue(pos, 1); // if value = 0.0 always -> 1.0
		}
		
		// CID flag
		if (pos == (mask.getDimension()-1)) {
			if (value == 1.0) { mask.setValue(pos, 0); }
			return 0.0;
		}
		
		// hidden states flag
		else {
			if (value == 2.0) {
				mask.setValue(pos, 1); // 2.0 always -> 1.0
				return Math.log(0.5);
			}
			
			else {
				if (value == 1.0) {
					if (Math.random() < 0.5) {
						mask.setValue(pos, 0); // half chance of going 0.0
					}
					else { mask.setValue(pos, 2); } // half chance of going 2.0
				
					return Math.log(2.0); // if value == 1.0
				}
				
				return Math.log(0.5); // if value == 0.0
			}
		}
	}
}
