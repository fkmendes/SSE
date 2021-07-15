package SSE;

import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;

import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends QuaSSEProcess {

    final public Input<Quant2MacroLinkFn> q2mLambdaInput = new Input<>("q2mLambda", "Function converting quantitative trait into lambda parameter.", Input.Validate.REQUIRED);
    final public Input<Quant2MacroLinkFn> q2mMuInput = new Input<>("q2mMu", "Function converting quantitative trait into mu parameter.", Input.Validate.REQUIRED);

    private double[] lambdaLo, muLo, lambdaHi, muHi; // macroevol parameters
    private Quant2MacroLinkFn q2mLambda, q2mMu;

    // state that matters for calculateLogP
    private double[][] esDs, scratch;

    @Override
    public void initAndValidate() {

        super.initAndValidate(); // read in all dimension-related stuff

        q2mLambda = q2mLambdaInput.get();
        lambdaLo = new double[nUsefulXbinsLo];
        lambdaHi = new double[nUsefulXbinsHi];

        q2mMu = q2mMuInput.get();
        muLo = new double[nUsefulXbinsLo];
        muHi = new double[nUsefulXbinsHi];

        populateMacroevolParams();

        int nDimensionsFFT = 2;
        initializeEsDs(nDimensionsFFT, nXbins);
    }

    @Override
    public void populateMacroevolParams() {
        lambdaLo = q2mLambda.getMacroParams(xLo, lambdaLo);
        lambdaHi = q2mLambda.getMacroParams(xHi,lambdaHi);
        muLo = q2mMu.getMacroParams(muLo, muLo);
        muHi = q2mMu.getMacroParams(muHi, muHi);
    }

    @Override
    public void initializeEsDs(int nDimensionsFFT, int nXbins) {
        // do stuff with nDimensionsFFT
    }

    @Override
    public void initializeTips() {

    }

    @Override
    public void pruneTree() {

    }

    @Override
    public void processBranch(Node node) {

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

    public double[] getLambda(boolean lowRes) {
        if (lowRes) return lambdaLo;
        else return lambdaHi;
    }

    public double[] getMu(boolean lowRes) {
        if (lowRes) return muLo;
        else return muHi;
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
