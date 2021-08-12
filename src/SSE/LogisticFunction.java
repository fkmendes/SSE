package SSE;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public class LogisticFunction extends BEASTObject implements LinkFn {

    final public Input<RealParameter> curveMaxBaseInput = new Input<>("curveMaxBase", "Curve maximum base value.", Input.Validate.REQUIRED);
    final public Input<RealParameter> added2CurveMaxInput = new Input<>("added2CurveMax", "How much is added to curve max base (after subtracting curve max base from this).", Input.Validate.REQUIRED);
    final public Input<RealParameter> sigmoidMidpointInput = new Input<>("sigmoidMidpoint", "Midpoint of sigmoid curve.", Input.Validate.REQUIRED);
    final public Input<RealParameter> logisticGrowthRateInput = new Input<>("logisticGrowthRate", "Growth rate of logistic curve.", Input.Validate.REQUIRED);

    private double y0, y1, x0, r, curveMax;
    private static final String LINKFUNCTION = "logistic";

    @Override
    public void initAndValidate() {
        y0 = curveMaxBaseInput.get().getValue();
        y1 = added2CurveMaxInput.get().getValue();
        r = logisticGrowthRateInput.get().getValue();
        x0 = sigmoidMidpointInput.get().getValue();
    }

    @Override
    public boolean refreshParams() {

        boolean refreshedSomething = false;

        if (curveMaxBaseInput.get().somethingIsDirty()) {
            y0 = curveMaxBaseInput.get().getValue();
            refreshedSomething = true;
        }

        if (added2CurveMaxInput.get().somethingIsDirty()) {
            y1 = added2CurveMaxInput.get().getValue();
            refreshedSomething = true;
        }

        if (logisticGrowthRateInput.get().somethingIsDirty()) {
            r = logisticGrowthRateInput.get().getValue();
            refreshedSomething = true;
        }

        if (sigmoidMidpointInput.get().somethingIsDirty()) {
            x0 = sigmoidMidpointInput.get().getValue();
            refreshedSomething = true;
        }

        return refreshedSomething;
    }

    @Override
    public double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        boolean refreshedSomething = false;
        if (!ignoreRefresh) refreshedSomething = refreshParams(); // if something changed in deterministic function parameters, we need to repopulate macroevol arrays

        /*
         * if either we don't care about refreshing, or we do and something was refreshed,
         * we repopulate y
         */
        if (ignoreRefresh || refreshedSomething) {
            if (x.length != y.length) throw new RuntimeException("Sizes of x (qu trait) and y (macroevol param) differ. Exiting...");

            for (int i = 0; i < x.length; i++) {
                curveMax = y1 - y0;
                y[i] = y0 + curveMax / (1.0 + Math.exp(r * (x0 - x[i])));
            }
        }

        return y;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

}
