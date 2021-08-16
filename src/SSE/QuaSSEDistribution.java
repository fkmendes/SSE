package SSE;

import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends QuaSSEProcess {

    final public Input<LinkFn> q2mLambdaInput = new Input<>("q2mLambda", "Function converting quantitative trait into lambda parameter.", Input.Validate.REQUIRED);
    final public Input<LinkFn> q2mMuInput = new Input<>("q2mMu", "Function converting quantitative trait into mu parameter.", Input.Validate.REQUIRED);
    final public Input<LinkFn> q2dInput = new Input<>("q2d", "Function converting quantitative trait into initial D values.", Input.Validate.REQUIRED);

    // dimension state
    private int nDimensions = 2;
    private int nDimensionsE = 1;
    private int nDimensionsD = 1;

    // macroevo state (for propagate in t)
    protected double[] birthRatesLo, deathRatesLo, birthRatesHi, deathRatesHi; // macroevol parameters
    private LinkFn q2mLambda, q2mMu, q2d;

    // state that matters directly for calculateLogP
    private double[][][] esDsLo, esDsHi; // first dimension are nodes, second is Es and Ds, third is each E (or D) along X ruler
    private double[][][] scratch;

    @Override
    public void initAndValidate() {

        super.initAndValidate(); // read in all dimension-related stuff
        int nNodes = tree.getNodeCount();

        q2mLambda = q2mLambdaInput.get();
        birthRatesLo = new double[nUsefulXbinsLo];
        birthRatesHi = new double[nUsefulXbinsHi];

        q2mMu = q2mMuInput.get();
        deathRatesLo = new double[nUsefulXbinsLo];
        deathRatesHi = new double[nUsefulXbinsHi];

        q2d = q2dInput.get();

        scratch = new double[nNodes][nDimensions][2 * nXbinsHi]; // for real and complex part after FFT

        populateMacroevolParams(true);

        int nDimensionsFFT = nDimensions;
        initializeEsDs(nNodes, nDimensionsFFT, nXbinsLo, nXbinsHi);

        populateTipsEsDs(nDimensionsFFT, nXbinsHi, true);
    }

    @Override
    public void populateMacroevolParams(boolean ignoreRefresh) {
        birthRatesLo = q2mLambda.getY(xLo, birthRatesLo, ignoreRefresh);
        birthRatesHi = q2mLambda.getY(xHi, birthRatesHi, ignoreRefresh);
        deathRatesLo = q2mMu.getY(deathRatesLo, deathRatesLo, ignoreRefresh);
        deathRatesHi = q2mMu.getY(deathRatesHi, deathRatesHi, ignoreRefresh);
    }

    @Override
    public void initializeEsDs(int nNodes, int nDimensionsFFT, int nXbinsLo, int nXbinsHi) {
        esDsLo = new double[nNodes][nDimensionsFFT][2 * nXbinsLo]; // 2 * for real and complex part after FFT
        esDsHi = new double[nNodes][nDimensionsFFT][2 * nXbinsHi];
        // do stuff with nDimensionsFFT
    }

    @Override
    public void populateTipsEsDs(int nDimensionsFFT, int nXbinsHi, boolean ignoreRefresh) {

        for (Node tip: tree.getExternalNodes()) {
            String tipName = tip.getID();
            int nodeIdx = tip.getNr();

            for (int i=0; i < nDimensionsFFT; i++) {
                // E's
                if (i == 0) {
                    // high res
                    for (int j=0; j < nXbinsHi; j++) esDsHi[nodeIdx][i][j] = 0.0; // E's = 1.0 - sampling.f

                    // low res (just for debugging)
                    for (int j=0; j < nXbinsLo; j++) esDsLo[nodeIdx][i][j] = 0.0; // E's = 1.0 - sampling.f
                }

                // D's
                /*
                 * TODO: if tips are not contemporanous (at present moment), but instead fossil, then we need to initialize esDsLo if t > tc
                 */
                else {
                    // high res
                    esDsHi[nodeIdx][i] = q2d.getY(xHi, esDsHi[nodeIdx][i], nLeftNRightFlanksHi, tipName, ignoreRefresh);

                    // low res (just for debugging)
                    esDsLo[nodeIdx][i] = q2d.getY(xLo, esDsLo[nodeIdx][i], nLeftNRightFlanksLo, tipName, ignoreRefresh);
                }
            }
        }
    }

    @Override
    public void pruneTree() {

    }

    @Override
    public void processBranch(Node node) {

        boolean lowRes = false;
        int nodeIdx = node.getNr();

        double startTime = node.getHeight(); // we're going backwards in time, toward the root
        if (startTime > tc) lowRes = true; // entire branch might already be > tc

        double[][] esDsAtNode;
        if (lowRes) esDsAtNode = esDsLo[nodeIdx];
        else esDsAtNode = esDsHi[nodeIdx];

        double[][] scratchAtNode = scratch[nodeIdx];

        boolean isFirstDt = true;
        while ((startTime + dt) <= node.getParent().getHeight()) {
            if (startTime > tc) lowRes = true;

            System.out.println("startTime + dt = " + (startTime + dt));
            System.out.println("node.getParent().getHeight() = " + node.getParent().getHeight());
            System.out.println("calling doIntegrateInPlace");
            doIntegrateInPlace(esDsAtNode, scratchAtNode, startTime, isFirstDt, lowRes);

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
    public void populatefY(boolean ignoreRefresh, boolean doFFT) {
        super.populatefY(ignoreRefresh, doFFT);
    }

    @Override
    public void doIntegrateInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double startTime, boolean isFirstDt, boolean lowRes) {

        // integrate over birth and death events (low or high resolution inside)
        propagateTInPlace(esDsAtNode, scratchAtNode, lowRes);

        // make normal kernel and FFTs it
        populatefY(true, true);

        // integrate over diffusion of substitution rate
        // propagateXInPlace(esDsAtNode, lowRes);
        propagateXInPlace(esDsAtNode, scratchAtNode, lowRes);

        // return esDsAtNode;
    }

    @Override
    public void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, boolean lowRes) {
        // grab scratch, dt and nDimensions from QuaSSEDistribution state
        if (lowRes) SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, dt, nUsefulXbinsLo, nDimensionsD);
        else SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, dt, nUsefulXbinsHi, nDimensionsD);
    }

    @Override
    public void propagateXInPlace(double[][] esDsAtNode, double[][] scratchAtNode, boolean lowRes) {

        // TODO: check if this should not be moved into SSEUtils
        for (int ithDim=0; ithDim < nDimensions; ithDim++) {
            for (int i=0; i < esDsAtNode[ithDim].length; i++) {
                scratchAtNode[ithDim][i] = esDsAtNode[ithDim][i];
            }
        }

        // System.out.println("esDsAtNode[1] = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratch[1] = " + Arrays.toString(scratch[1]));

        // grab dt and nDimensions from state
        if (lowRes) SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYLo, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], nDimensionsE, nDimensionsD, fftForEandDLo);
        else SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYHi, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], nDimensionsE, nDimensionsD, fftForEandDHi);
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

    public double[] getfY(boolean lowRes) {
        if (lowRes) return fYLo;
        else return fYHi;
    }

    public double[][][] getEsDs(boolean lowRes) {
        if (lowRes) return esDsLo;
        else return esDsHi;
    }

    public double[][] getEsDsAtNode(int nodeIdx, boolean lowRes) {
        if (lowRes) return esDsLo[nodeIdx];
        return esDsHi[nodeIdx];
    }

    public int getNUsefulTraitBins(boolean lowRes) {
        if (lowRes) return nUsefulXbinsLo;
        else return nUsefulXbinsHi;
    }

    public void setEsDsAtNodeElementAtDim(int nodeIdx, int dim, int eleIdx, double val, boolean lowRes) {
        if (lowRes) esDsLo[nodeIdx][dim][eleIdx] = val;
        else esDsHi[nodeIdx][dim][eleIdx] = val;
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
