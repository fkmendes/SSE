package SSE;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;

import java.util.List;
import java.util.Random;

public class LogisticFunction extends Distribution {

    final public Input<RealParameter> xInput = new Input<>("x", "Independent variable of logistic function.", Input.Validate.REQUIRED);
    final public Input<RealParameter> curveMaxBaseInput = new Input<>("curveMaxBase", "Curve maximum base value.", Input.Validate.REQUIRED);
    final public Input<RealParameter> added2CurveMaxInput = new Input<>("added2CurveMax", "How much is added to curve max base (after subtracting curve max base from this).", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmoidMidpointInput = new Input<>("sigmoidMidpoint", "Midpoint of sigmoid curve.", Input.Validate.REQUIRED);
    final public Input<RealParameter> logisticGrowthRateInput = new Input<>("logisticGrowthRate", "Growth rate of logistic curve.", Input.Validate.REQUIRED);

    private double y0, y1, x0, r, curveMax;
    private Double[] x, y;

    @Override
    public void initAndValidate() {
        refreshParams();
        y = new Double[x.length];
    }

    public void refreshParams() {
        x = xInput.get().getValues();
        y0 = curveMaxBaseInput.get().getValue();
        y1 = added2CurveMaxInput.get().getValue();
        r = logisticGrowthRateInput.get().getValue();
        x0 = sigmoidMidpointInput.get().getValue();
    }

    public Double[] getY() {
        refreshParams();

        for (int i=0; i<x.length; i++) {
            curveMax = y1 - y0;
            y[i] = y0 + curveMax / (1.0 + Math.exp(r * (x0 - x[i])));
        }

        return y;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
