package SSE;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public abstract class Macro2QuantLinkFn extends BEASTObject  {

    final public Input<RealParameter> xInput = new Input<>("x", "Independent variable (quantitative trait) of logistic function.", Input.Validate.REQUIRED);

    protected Double[] x, y; // y is the macroevolutionary parameter

    @Override
    public void initAndValidate() {
        refreshParams();
        y = new Double[x.length];
    }

    protected abstract void refreshParams();

    protected abstract Double[] getMacroParams();

}
