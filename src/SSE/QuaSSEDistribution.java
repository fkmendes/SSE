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
    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each quantitative trait state bin.", Input.Validate.REQUIRED);
    final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each quantitative trait state bin.", Input.Validate.REQUIRED);
    final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies at root.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dtInput = new Input<>("dt", "Length of time interval over which integration is carried out.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> nXbinsInput = new Input<>("nX", "Total number of quantitative trait bins after discretization.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dXBinInput = new Input<>("dX", "Width of quantitative trait bins.", Input.Validate.REQUIRED);
    final public Input<RealParameter> driftInput = new Input<>("drift", "Drift term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> flankWidthScalerInput = new Input<>("flankWidthScaler", "Multiplier of normal standard deviation when determining number of flanking bins.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> highResScalerInput = new Input<>("hiScaler", "Scale nX by this when at high resolution.", Input.Validate.REQUIRED);


    /*
     * changeInXNormalMean (=diversitree's drift)
     * changeInXNormalVar (=diversitree's diffusion)
     */
    private double dt, tc, dXbin, changeInXNormalMean, changeInXNormalSd;
    private double flankWidthScaler, nUsefulBinsHi, nUsefulBinsLo;
    private int nXbins, nXbinsHi, nXbinsLo, highScaler;
    private double[] nLeftNRightFlanksHi, nLeftNRightFlanksLo;
    private double[] xHi, lambdaHi, muHi;
    private double[] xLo, lambdaLo, muLo;
    private double[][] esDs, scratch;

    @Override
    public void initAndValidate() {

        dt = dtInput.get().getValue();
        dXbin = dXBinInput.get().getValue();
        nXbins = nXbinsInput.get().getValue();
        flankWidthScaler = flankWidthScalerInput.get().getValue();
        highScaler = highResScalerInput.get().getValue();
        changeInXNormalMean = driftInput.get().getValue() * dt;
        changeInXNormalSd = Math.sqrt(diffusionInput.get().getValue()) * dt;

        prepareDimensionsInPlace();
        populateXLambdaMu();
    }

    private void prepareDimensionsInPlace() {

        // left flank bins
        nLeftNRightFlanksLo[0] = Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        nLeftNRightFlanksHi[0] = highScaler * Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        // right flank bins
        nLeftNRightFlanksLo[1] = Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;
        nLeftNRightFlanksHi[1] = highScaler * Math.ceil(changeInXNormalMean - flankWidthScaler * changeInXNormalSd) / dXbin;

        nXbinsLo = nXbins;
        nXbinsHi = highScaler * nXbinsLo;

        nUsefulBinsLo = nXbinsLo - (nLeftNRightFlanksLo[0] + 1 + nLeftNRightFlanksLo[1]);
        nUsefulBinsHi = nXbinsHi - (nLeftNRightFlanksHi[0] + 1 + nLeftNRightFlanksHi[1]);
    }

    private void populateXLambdaMu() {
        ;
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
