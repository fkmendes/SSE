package SSE;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class LogisticFunction extends Macro2QuantLinkFn {

    final public Input<RealParameter> curveMaxBaseInput = new Input<>("curveMaxBase", "Curve maximum base value.", Input.Validate.REQUIRED);
    final public Input<RealParameter> added2CurveMaxInput = new Input<>("added2CurveMax", "How much is added to curve max base (after subtracting curve max base from this).", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmoidMidpointInput = new Input<>("sigmoidMidpoint", "Midpoint of sigmoid curve.", Input.Validate.REQUIRED);
    final public Input<RealParameter> logisticGrowthRateInput = new Input<>("logisticGrowthRate", "Growth rate of logistic curve.", Input.Validate.REQUIRED);

    private double y0, y1, x0, r, curveMax;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public void refreshParams() {
        x = xInput.get().getValues();
        y0 = curveMaxBaseInput.get().getValue();
        y1 = added2CurveMaxInput.get().getValue();
        r = logisticGrowthRateInput.get().getValue();
        x0 = sigmoidMidpointInput.get().getValue();
    }

    @Override
    public Double[] getMacroParams() {
        refreshParams();

        for (int i=0; i<x.length; i++) {
            curveMax = y1 - y0;
            y[i] = y0 + curveMax / (1.0 + Math.exp(r * (x0 - x[i])));
        }

        return y;
    }
}
