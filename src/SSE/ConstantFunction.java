package SSE;

import beast.core.Input;
import beast.core.parameter.RealParameter;

/*
 * Applies the same y (macroevol param) value to all x (qu trait) bins.
 * If you put a prior on the y value, it's the same as assuming y is
 * distributed according to that prior, and independent of x
 */
public class ConstantFunction extends Quant2MacroLinkFn {

    final public Input<RealParameter> yValueInput = new Input<>("yV", "Constant value of dependent variable (quantitative trait).", Input.Validate.REQUIRED);

    double yValue;
    private static final String LINKFUNCTION = "constant";

    @Override
    public void initAndValidate() {
        yValue = yValueInput.get().getValue();
    }

    @Override
    protected void refreshParams() { }

    @Override
    public double[] getMacroParams(double[] x, double[] y) {

        if (x.length != y.length) throw new RuntimeException("Sizes of x (qu trait) and y (macroevol param) differ. Exiting...");

        for (int i=0; i<x.length; i++) {
            y[i] = yValue;
        }

        return y;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

}
