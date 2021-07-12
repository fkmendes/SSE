package SSE;

import beast.core.Input;
import beast.core.parameter.RealParameter;

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
    public Double[] getMacroParams(Double[] xHi) {

        if (yHi == null) throw new RuntimeException("Quantitative trait ruler has not been initialized. Exiting...");

        for (int i=0; i<xHi.length; i++) {
            yHi[i] = yValue;
        }

        return yHi;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

    @Override
    protected double getXMid() {
        return 0;
    }
}
