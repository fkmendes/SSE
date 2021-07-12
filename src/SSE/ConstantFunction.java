package SSE;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class ConstantFunction extends Macro2QuantLinkFn {

    final public Input<RealParameter> yValueInput = new Input<>("yV", "Constant value of dependent variable (quantitative trait).", Input.Validate.REQUIRED);

    double yValue;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    protected void refreshParams() {
        x = xInput.get().getValues();
        yValue = yValueInput.get().getValue();
    }

    @Override
    public Double[] getMacroParams() {

        for (int i=0; i<x.length; i++) {
            y[i] = yValue;
        }

        return y;
    }
}
