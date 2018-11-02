package SSE;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class MasqueradeBall extends CalculationNode {
	final public Input<RealParameter> modelMaskInput = new Input<>("modelMask", "Series of integers that collectively determine what (sub)model to use.");
	final public Input<HiddenInstantaneousRateMatrix> hirmInput = new Input<>("hiddenInstantaneousRateMatrix", "HiddenInstantaneousRateMatrix object containing anagenenetic rates for both observed and hidden states.", Validate.REQUIRED);
	final public Input<LambdaMuAssigner> lambdaMuAssignerInput = new Input<>("lambdaMuAssigner", "LambdaMuAssigner object that assigns distinct parameters to each state.", Validate.REQUIRED);

	private Double[] mask;
	private Integer[] hiddenStatesOnOrOff;
	private boolean masqueradeBallDirty = true;

	private int numberOfStates; // number of observed states
	private Double[] lambda;
	private Double[] mu;
	private Double[][] q;
	private HiddenInstantaneousRateMatrix hirm;
	
//	private boolean disallowDoubleTransitions;
//	private int symmetrifyAcrossDiagonal;
//	private boolean ignoreDiagonal; // right now, always true (we never query Qij for i=j in SSEODE)
//	private String lambdaToStatesString;
//	private String muToStatesString;

	@Override
	public void initAndValidate() {
		mask = new Double[numberOfStates];
		hiddenStatesOnOrOff = new Integer[numberOfStates];
		hirm = hirmInput.get();

		applyMask();
	}
	
	public void applyMask() {
		modelMaskInput.get().getValues(mask);
		int hiddenStateCount = 0;
		
		for (int i = 0; i < mask.length; ++i) {
			double maskItem = mask[i];
			
			// no hidden state for this particular observed state
			if (maskItem == 0) {
				hiddenStatesOnOrOff[i] = -1;
			}
			
			// add hidden state for this particular observed state, but transition from and to are symmetrical
			if (maskItem == 1) {
				int obsStateToSymmetrifyIdx = i;
				hirm.symmetrifyAcrossDiagonal(obsStateToSymmetrifyIdx);
				hiddenStatesOnOrOff[i] = hiddenStateCount;
				hiddenStateCount++;
			};
			
			// add hidden state for this particular observed state, with different transition rates
			if (maskItem == 2) {
				hiddenStatesOnOrOff[i] = hiddenStateCount;
				hiddenStateCount++;
			}
		}
		
		masqueradeBallDirty = false;
	}
	
	// getters
	public Double[] getLambdas() {
		// if lambdas were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return lambda;
	}
	
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return mu;
	}
	
	public Double[][] getQs() {
		// if qs were operated on, things are dirty, we need to update
		if (masqueradeBallDirty) {
			applyMask();
		}
		return q;
	}
	
	protected boolean requiresRecalculation() {
		masqueradeBallDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		applyMask();
		super.restore();
	}
}
