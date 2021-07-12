package SSE;

import beast.core.BEASTObject;

public abstract class quant2MacroLinkFn extends BEASTObject {

    protected Double[] xLo, yLo, xHi, yHi; // y is the macroevolutionary parameter

    @Override
    public void initAndValidate() {}

    protected abstract void refreshParams();

    protected abstract Double[] getMacroParams();

    public void setYSizeAndX(int ySize, Double[] x) {
        xHi = x;
        yHi = new Double[ySize];
    }

}
