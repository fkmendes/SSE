package SSE;

import beast.core.BEASTObject;

public abstract class Quant2MacroLinkFn extends BEASTObject {

    // x rulers and y array values, low and high resolution
    // y is the macroevolutionary parameter
    protected Double[] xLo, yLo, xHi, yHi;

    @Override
    public void initAndValidate() {}

    protected abstract void refreshParams();

    protected abstract Double[] getMacroParams(Double[] x);

    public void setYSize(int ySize) {
        yHi = new Double[ySize];
    }

    protected abstract String getLinkFnName();

    protected abstract double getXMid();
}
