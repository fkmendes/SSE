package SSE;

import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;

import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends QuaSSEProcess {

    final public Input<Quant2MacroLinkFn> q2mLambdaInput = new Input<>("q2mLambda", "Function converting quantitative trait into lambda parameter.", Input.Validate.REQUIRED);
    final public Input<Quant2MacroLinkFn> q2mMuInput = new Input<>("q2mMu", "Function converting quantitative trait into mu parameter.", Input.Validate.REQUIRED);

    private int nDimensions = 2; // 1 for E, 1 for D
    protected double[] birthRatesLo, deathRatesLo, birthRatesHi, deathRatesHi; // macroevol parameters
    private Quant2MacroLinkFn q2mLambda, q2mMu;

    // state that matters for calculateLogP
    private double[][][] esDs; // first dimension are nodes, second is Es and Ds, third is each E (or D) along X ruler
    private double[][] scratch;

    @Override
    public void initAndValidate() {

        super.initAndValidate(); // read in all dimension-related stuff

        q2mLambda = q2mLambdaInput.get();
        birthRatesLo = new double[nUsefulXbinsLo];
        birthRatesHi = new double[nUsefulXbinsHi];

        q2mMu = q2mMuInput.get();
        deathRatesLo = new double[nUsefulXbinsLo];
        deathRatesHi = new double[nUsefulXbinsHi];

        populateMacroevolParams(true);

        int nDimensionsFFT = 2;
        initializeEsDs(tree.getNodeCount(), nDimensionsFFT, nXbins);

        populateTipsEsDs(nDimensionsFFT, nUsefulXbinsHi);
    }

    @Override
    public void populateMacroevolParams(boolean ignoreRefresh) {
        birthRatesLo = q2mLambda.getMacroParams(xLo, birthRatesLo, ignoreRefresh);
        birthRatesHi = q2mLambda.getMacroParams(xHi, birthRatesHi, ignoreRefresh);
        deathRatesLo = q2mMu.getMacroParams(deathRatesLo, deathRatesLo, ignoreRefresh);
        deathRatesHi = q2mMu.getMacroParams(deathRatesHi, deathRatesHi, ignoreRefresh);
    }

    @Override
    public void initializeEsDs(int nNodes, int nDimensionsFFT, int nXbins) {
        esDs = new double[nNodes][nDimensionsFFT][nXbins];
        // do stuff with nDimensionsFFT
    }

    @Override
    public void populateTipsEsDs(int nDimensionsFFT, int nUsefulXbins) {
        for (Node tip: tree.getExternalNodes()) {
            String tipName = tip.getID();
            int nodeIdx = tip.getNr();

            for (int i=0; i < nDimensionsFFT; i++) {
                for (int j = 0; j < nUsefulXbins; j++) {
                    if (j == 0) esDs[nodeIdx][i][j] = 0.0; // E's = 1.0 - sampling.f
                    else {
                        esDs[nodeIdx][i][j] = 0.0; // TODO: dnorm()
                    }

                    System.out.println(tipName + " = " + quTraits.getValue(tipName));
                }
            }
        }
    }

    @Override
    public void pruneTree() {

    }

    @Override
    public void processBranch(Node node) {

        int nodeIdx = node.getNr();
        double[][] esDsAtNode = esDs[nodeIdx];

        boolean isFirstDt = true;
        boolean lowRes = false;

        double startTime = node.getHeight(); // we're going backwards in time, toward the root

        while ((startTime + dt) <= node.getParent().getHeight()) {
            if (startTime > tc) lowRes = true;

            doIntegrate(esDsAtNode, startTime, isFirstDt, lowRes);

            isFirstDt = false;
            startTime += dt;
        }

    }

    @Override
    public void processInternalNode() {

    }

    @Override
    public void processRootNode() {

    }

    @Override
    public void doIntegrate(double[][] esDsAtNode, double startTime, boolean isFirstDt, boolean lowRes) {

        // integrate over birth and death events (low or high resolution inside)
        propagateT(esDsAtNode, lowRes);

        // integrate over diffusion of substitution rate
        propagateX(lowRes);
    }

    @Override
    public void propagateT(double[][] esDsAtNode, boolean lowRes) {
        // grab scratch, dt and nDimensions from QuaSSEDistribution state
        if (lowRes) SSEUtils.propagateEandDinTQuaSSE(esDsAtNode, scratch, birthRatesLo, deathRatesLo, dt, nUsefulXbinsLo, nDimensions);
        else SSEUtils.propagateEandDinTQuaSSE(esDsAtNode, scratch, birthRatesHi, deathRatesHi, dt, nUsefulXbinsHi, nDimensions);
    }

    @Override
    public void propagateX(boolean lowRes) {

    }

    @Override
    public void convolve() {

    }

    @Override
    public double calculateLogP() {

        populateMacroevolParams(false);

        return 0.0;
    }


    /*
     * Getters and setters
     */

    public double[] getLambda(boolean lowRes) {
        if (lowRes) return birthRatesLo;
        else return birthRatesHi;
    }

    public double[] getMu(boolean lowRes) {
        if (lowRes) return deathRatesLo;
        else return deathRatesHi;
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
