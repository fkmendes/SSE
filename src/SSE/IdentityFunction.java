package SSE;

public class IdentityFunction extends Quant2MacroLinkFn {

    private static final String LINKFUNCTION = "identity";

    @Override
    protected void refreshParams() {}

    @Override
    public Double[] getMacroParams(Double[] x) {

        if (yHi == null) throw new RuntimeException("Quantitative trait ruler has not been initialized. Exiting...");

        return x;
    }

    public String getLinkFnName() {
        return LINKFUNCTION;
    }

    @Override
    protected double getXMid() {
        return 0;
    }
}
