package SSE;

import beast.core.BEASTObject;

public abstract class Quant2MacroLinkFn extends BEASTObject {

    @Override
    public void initAndValidate() {}

    protected abstract void refreshParams();

    // x = qu trait, y = macroevol parameter
    protected abstract double[] getMacroParams(double[] x, double[] y);

    protected abstract String getLinkFnName();

}
