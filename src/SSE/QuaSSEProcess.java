package SSE;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Specifies a quantitative trait(s) state-dependent speciation and" +
        "extinction birth-death process.")
public abstract class QuaSSEProcess extends Distribution {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    // final public Input<RealParameter> quTraitsInput = new Input<>("quTraits", "Quantitative trait values observed at tips", Input.Validate.REQUIRED);
    final public Input<RealParameter> dtInput = new Input<>("dt", "Length of time interval over which integration is carried out.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> nXbinsInput = new Input<>("nX", "Total number of quantitative trait bins after discretization at low resolution.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dXBinInput = new Input<>("dX", "Width of quantitative trait bins.", Input.Validate.REQUIRED);
    final public Input<RealParameter> xMidInput = new Input<>("xMid", "Midpoint to center the quantitative trait space.", Input.Validate.REQUIRED);
    final public Input<RealParameter> driftInput = new Input<>("drift", "Drift term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> flankWidthScalerInput = new Input<>("flankWidthScaler", "Multiplier of normal standard deviation when determining number of flanking bins.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> highLowRatioInput = new Input<>("hiLoRatio", "Scale nX by this when at high resolution.", Input.Validate.REQUIRED);

    protected Tree tree;
    protected RealParameter quTraits;

    // state for dimensioning things and setting up resolution of integration
    protected double dt, tc;
    protected double dXbin, flankWidthScaler, xMinLo, xMinHi, xMid;
    protected int nXbins, nXbinsHi, nUsefulXbinsHi, nUsefulXbinsLo, nXbinsLo, hiLoRatio;
    protected int[] nUsefulXbins, nLeftNRightFlanksHi, nLeftNRightFlanksLo;
    protected double[] xLo, xHi; // x rulers

    // quantitative trait evolution
    protected double changeInXNormalMean; // (=diversitree's drift)
    protected double changeInXNormalSd; // (=diversitree's diffusion)

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        // quTraits = quTraitsInput.get();

        nLeftNRightFlanksLo = new int[2];
        nLeftNRightFlanksHi = new int[2];
        nUsefulXbins = new int[2];

        dt = dtInput.get().getValue();
        dXbin = dXBinInput.get().getValue();
        nXbins = nXbinsInput.get().getValue();
        xMid = xMidInput.get().getValue();
        flankWidthScaler = flankWidthScalerInput.get().getValue();
        hiLoRatio = highLowRatioInput.get().getValue();
        changeInXNormalMean = driftInput.get().getValue() * dt;
        changeInXNormalSd = Math.sqrt(diffusionInput.get().getValue() * dt);

        prepareDimensionsInPlace(); // in parent class
        prepareXRulers(); // in parent class
    }

    /*
     *
     */
    protected void prepareDimensionsInPlace() {

        // left flank bins
        nLeftNRightFlanksLo[0] = (int)(Math.ceil(-(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin));
        nLeftNRightFlanksHi[0] = hiLoRatio * nLeftNRightFlanksLo[0];

        // right flank bins
        nLeftNRightFlanksLo[1] = (int)(Math.ceil(-(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin));
        nLeftNRightFlanksHi[1] = hiLoRatio * nLeftNRightFlanksLo[1];

        nXbinsLo = nXbins;
        nXbinsHi = hiLoRatio * nXbinsLo;

        nUsefulXbinsLo = nXbinsLo - (nLeftNRightFlanksLo[0] + 1 + nLeftNRightFlanksLo[1]);
        nUsefulXbinsHi = nXbinsHi - (nLeftNRightFlanksHi[0] + 1 + nLeftNRightFlanksHi[1]);
        nUsefulXbins[0] = nUsefulXbinsLo;
        nUsefulXbins[1] = nUsefulXbinsHi;
    }

    /*
     *
     */
    protected void prepareXRulers() {

        xMinLo = xMid - dXbin * Math.ceil((nUsefulXbinsLo - 1.0) / 2.0);
        xMinHi = xMinLo - dXbin * (1.0 - 1.0 / hiLoRatio);

        // preparing x rulers
        xLo = new double[nUsefulXbinsLo];
        xLo[0] = xMinLo;
        for (int i = 1; i < nUsefulXbinsLo; i++) {
            xLo[i] = xLo[i-1] + dXbin;
        }

        xHi = new double[nUsefulXbinsHi];
        xHi[0] = xMinHi;
        for (int i = 1; i< nUsefulXbinsHi; i++) {
            xHi[i] = xHi[i-1] + dXbin / hiLoRatio;
        }
    }

    /*
     *
     */
    protected abstract void populateMacroevolParams(boolean ignoreRefresh);

    /*
     *
     */
    protected abstract void initializeEsDs(int nNodes, int nDimensionsFFT, int nXbinsHi);

    /*
     *
     */
    protected abstract void populateTipsEsDs(int nDimensionsFFT, int nXbins, boolean ignoreRefresh);

    /*
     *
     */
    protected abstract void pruneTree();

    /*
     *
     */
    protected abstract void processBranch(Node node);

    /*
     *
     */
    protected abstract void processInternalNode();

    /*
     *
     */
    protected abstract void processRootNode();

    /*
     *
     */
    protected abstract void doIntegrate(double[][] esDsAtNode, double startTime, boolean isFirstDt, boolean lowRes);

    /*
     *
     */
    protected abstract void propagateT(double[][] esDsAtNode, boolean lowRes);

    /*
     *
     */
    protected abstract void propagateX(boolean lowRes);

    /*
     *
     */
    protected abstract void convolve();

    public int getnXbins(boolean lowRes) {
        if (lowRes) return nXbins;
        else return nXbinsHi;
    }

    public int getNUsefulXbins(boolean lowRes) {
        if (lowRes) return nUsefulXbinsLo;
        else return nUsefulXbinsHi;
    }

    public int getNLeftFlanks(boolean lowRes) {
        if (lowRes) return nLeftNRightFlanksLo[0];
        else return nLeftNRightFlanksHi[0];
    }

    public int getNRightFlanks(boolean lowRes) {
        if (lowRes) return nLeftNRightFlanksLo[1];
        else return nLeftNRightFlanksHi[1];
    }

    public double getXMinLo() {
        return xMinLo;
    }

    public double getXMinHi() {
        return xMinHi;
    }

    public double[] getX(boolean lowRes) {
        if (lowRes) return xLo;
        else return xHi;
    }

}
