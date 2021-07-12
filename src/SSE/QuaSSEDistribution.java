package SSE;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends Distribution implements QuaSSEProcess {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    // final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each quantitative trait state bin.", Input.Validate.REQUIRED);
    // final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each quantitative trait state bin.", Input.Validate.REQUIRED);
    final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies at root.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dtInput = new Input<>("dt", "Length of time interval over which integration is carried out.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> nXbinsInput = new Input<>("nX", "Total number of quantitative trait bins after discretization.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dXBinInput = new Input<>("dX", "Width of quantitative trait bins.", Input.Validate.REQUIRED);
    final public Input<RealParameter> driftInput = new Input<>("drift", "Drift term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> flankWidthScalerInput = new Input<>("flankWidthScaler", "Multiplier of normal standard deviation when determining number of flanking bins.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> highLowRatioInput = new Input<>("hiLoRatio", "Scale nX by this when at high resolution.", Input.Validate.REQUIRED);
    final public Input<Quant2MacroLinkFn> q2mLambdaInput = new Input<>("q2mLambda", "Function converting quantitative trait into lambda parameter.", Input.Validate.REQUIRED);
    final public Input<Quant2MacroLinkFn> q2mMuInput = new Input<>("q2mMu", "Function converting quantitative trait into mu parameter.", Input.Validate.REQUIRED);

    /*
     * changeInXNormalMean (=diversitree's drift)
     * changeInXNormalVar (=diversitree's diffusion)
     */
    private double dt, tc, dXbin, changeInXNormalMean, changeInXNormalSd;
    private double flankWidthScaler, nUsefulBinsHi, nUsefulBinsLo, xMinLo, xMinHi;
    private int nXbins, nXbinsHi, nXbinsLo, hiLoRatio;
    private double[] nLeftNRightFlanksHi, nLeftNRightFlanksLo;
    private Double[] xHi, lambdaHi, muHi; // xHi is the x ruler at high res
    private Double[] xLo, lambdaLo, muLo; // xLo is the x ruler at low res
    private double[][] esDs, scratch;
    private Quant2MacroLinkFn q2mLambda, q2mMu;

    @Override
    public void initAndValidate() {

        dt = dtInput.get().getValue();
        dXbin = dXBinInput.get().getValue();
        nXbins = nXbinsInput.get().getValue();
        flankWidthScaler = flankWidthScalerInput.get().getValue();
        hiLoRatio = highLowRatioInput.get().getValue();
        changeInXNormalMean = driftInput.get().getValue() * dt;
        changeInXNormalSd = Math.sqrt(diffusionInput.get().getValue()) * dt;
        q2mLambda = q2mLambdaInput.get();
        q2mMu = q2mMuInput.get();

        prepareDimensionsInPlace();
        populateXLambdaMu();
    }

    private void prepareDimensionsInPlace() {

        // left flank bins
        nLeftNRightFlanksLo[0] = Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        nLeftNRightFlanksHi[0] = hiLoRatio * Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        // right flank bins
        nLeftNRightFlanksLo[1] = Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        nLeftNRightFlanksHi[1] = hiLoRatio * Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;

        nXbinsLo = nXbins;
        nXbinsHi = hiLoRatio * nXbinsLo;

        nUsefulBinsLo = nXbinsLo - (nLeftNRightFlanksLo[0] + 1 + nLeftNRightFlanksLo[1]);
        nUsefulBinsHi = nXbinsHi - (nLeftNRightFlanksHi[0] + 1 + nLeftNRightFlanksHi[1]);

        if (q2mLambda.getLinkFnName() == "logistic") {
            xMinLo = q2mLambda.getXMid() - dXbin * Math.ceil((nUsefulBinsLo - 1.0) / 2.0);
        }

        xMinHi = xMinLo - dXbin * (1.0 - 1.0 / hiLoRatio);

        // preparing x rulers
        xLo[0] = xMinLo;
        for (int i=1; i<nUsefulBinsLo; i++) {
            xLo[i] = xLo[-1] + dXbin;
        }

        xHi[0] = xMinHi;
        for (int i=1; i<nUsefulBinsLo; i++) {
            xHi[i] = xHi[-1] + dXbin / hiLoRatio;
        }
    }

    private void populateXLambdaMu() {
        lambdaLo = q2mLambda.getMacroParams(xLo);
        lambdaHi = q2mLambda.getMacroParams(xHi);
    }

    @Override
    public void initializeTips() {

    }

    @Override
    public void pruneTree() {

    }

    @Override
    public void processBranch() {

    }

    @Override
    public void processInternalNode() {

    }

    @Override
    public void processRootNode() {

    }

    @Override
    public void doIntegrate() {

    }

    @Override
    public void propagateT() {

    }

    @Override
    public void propagateX() {

    }

    @Override
    public void convolve() {

    }

    /*
     * Getters and setters
     */

    public int getNXbinsHi() {
        return nXbinsHi;
    }

    public int getNXbinsLo() {
        return nXbinsLo;
    }

    public double getNUsefulBinsHi() {
        return nUsefulBinsHi;
    }

    public double getNUsefulBinsLo() {
        return nUsefulBinsLo;
    }

    public double getNLeftFlankHi() {
        return nLeftNRightFlanksHi[0];
    }

    public double getNRightFlankHi() {
        return nLeftNRightFlanksHi[1];
    }

    public double getNLeftFlankLo() {
        return nLeftNRightFlanksLo[0];
    }

    public double getNRightFlankLo() {
        return nLeftNRightFlanksLo[1];
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
