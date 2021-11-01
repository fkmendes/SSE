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
    private double[][][] scratchLo, scratchHi;
    private double[] logNormalizationFactors;

    @Override
    public void initAndValidate() {

        super.initAndValidate(); // read in all dimension-related stuff, populates fYLo and fYHi

        int nNodes = tree.getNodeCount();

        logNormalizationFactors = new double[nNodes];

        q2mLambda = q2mLambdaInput.get();
        birthRatesLo = new double[nUsefulXbinsLo];
        birthRatesHi = new double[nUsefulXbinsHi];

        q2mMu = q2mMuInput.get();
        deathRatesLo = new double[nUsefulXbinsLo];
        deathRatesHi = new double[nUsefulXbinsHi];

        q2d = q2dInput.get();

        scratchLo = new double[nNodes][nDimensions][2 * nXbinsLo]; // for real and complex part after FFT
        scratchHi = new double[nNodes][nDimensions][2 * nXbinsHi]; // for real and complex part after FFT

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
    public void startRecursionAtRootNode(Node rootNode) {
        processInternalNode(rootNode); // start recursion

        // do stuff after recursion is done
        // TODO: stuff
    }

    @Override
    public void processInternalNode(Node aNode) {
        // recur if internal node or sampled ancestor
        if (!aNode.isLeaf()) {
            for (Node childNode: aNode.getChildren()) {
                processInternalNode(childNode); // recur
            }
        }

        // after recursion, do its own branch
        integrateBranch(aNode);
    }

    @Override
    public void integrateBranch(Node aNode) {
        boolean lowRes = false;
        int nodeIdx = aNode.getNr();

        // dealing with branch lengths
        double startTime = aNode.getHeight(); // we're going backwards in time, toward the root
        double endTime = aNode.getParent().getHeight();
        double branchLength2Integrate = endTime - startTime;

        // dealing with dt
        double dt, nIntervals, nonIntegratedDt;
        if (dynamicallyAdjustDt) {
            nIntervals = Math.ceil(branchLength2Integrate / dtMax);
            dt = branchLength2Integrate / nIntervals; // dynamically adjusted dt
        } else {
            nIntervals = Math.floor(branchLength2Integrate / dtMax);
            dt = dtMax;
            nonIntegratedDt = branchLength2Integrate % dtMax;
        }

        /*
         * Sorting out E's and D's prior to integration
         * (1) Grabbing esDs and scratch depending on resolution
         * (2) Calculating normalization factor
         * (3) Storing normalization factor to unnormalize at the very end of pruning
         * (4) Normalize DsThis is done only once if node is already at low, rather than
         * at every loop of the while block below
         */
        double[][] esDsAtNode, scratchAtNode;

        // option 1: entire branch might already be > tc (all low-res)
        if (startTime >= tc) {
            lowRes = true;
            esDsAtNode = esDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];
            // debug
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));
        }
        // option 2: entire branch < tc (all high-res)
        // option 3: tc happens inside branch (tip-end part in high-res, root-end part in low-res)
        else {
            esDsAtNode = esDsHi[nodeIdx];
            scratchAtNode = scratchHi[nodeIdx];
            // debug
            // System.out.println("Unnormalized D's at high-res = " + Arrays.toString(esDsAtNode[1]));
        }

        /*
         * Normalizing D's to avoid underflow
         * This is done prior to integration, at the right resolution already
         * (since we have already gone through options 1 -- 3 above)
         */
        double normalizationFactorFromDs = 0.0;
        if (lowRes) normalizationFactorFromDs = calculateNormalizationFactor(esDsAtNode[1], dXbin);
        else normalizationFactorFromDs = calculateNormalizationFactor(esDsAtNode[1], (dXbin / hiLoRatio));
        logNormalizationFactors[nodeIdx] = normalizationFactorFromDs; // storing normalization factors for when returning log-lik we can de-normalize it

        for (int i=0; i<esDsAtNode[1].length; i++) {
            esDsAtNode[1][i] /= normalizationFactorFromDs;
            // there will be an additional step of normalization when in option 3
            // as a result of integrating the tip-end part in high-res
        }

        // debug
        // System.out.println("Normalized D's at either low or high-res = " + Arrays.toString(esDsAtNode[1]));

        // integrating!
        double integrationTStart = startTime;
        for (int i=0; i<nIntervals; i++) {
            // tc happens inside branch (tip-end part in high-res, root-end part in low-res)
            if ((integrationTStart + dt) >= tc) {
                // integration of this tip-end in high res
                double dtAtHi = tc - integrationTStart;
                doIntegrateInPlace(esDsAtNode, scratchAtNode, dtAtHi, false);

                // one more normalization for this high-res stretch
                normalizationFactorFromDs = calculateNormalizationFactor(esDsAtNode[1], (dXbin / hiLoRatio));
                for (int j=0; j<esDsAtNode[1].length; j++) esDsAtNode[1][j] /= normalizationFactorFromDs;
                logNormalizationFactors[nodeIdx] = normalizationFactorFromDs; // storing normalization factors for when returning log-lik we can de-normalize it

                // integration of the root-end in low res
                esDsAtNode = esDsLo[nodeIdx];
                scratchAtNode = scratchLo[nodeIdx];
                double dtAtLo = integrationTStart + dt - tc;
                doIntegrateInPlace(esDsAtNode, scratchAtNode, dtAtLo, true);
            }
            // otherwise
            else {
                doIntegrateInPlace(esDsAtNode, scratchAtNode, dt, lowRes);
                integrationTStart += dt;
            }
        }
    }

    /*
     * Math-y methods start below
     */
    @Override
    public void populatefY(boolean ignoreRefresh, boolean doFFT) {
        super.populatefY(ignoreRefresh, doFFT);
    }

    private double calculateNormalizationFactor(double[] dS, double binSize) {
        double normalizationFactorFromDs = 0.0;
        for (double d: dS) normalizationFactorFromDs += d;
        return (normalizationFactorFromDs *= binSize);
    }

    @Override
    public void populatePriorProbAtRoot(String rootPriorType) {
        if (rootPriorType == FLAT) {
            ;
        }
        else if (rootPriorType == OBS) {
            ;
        }
        else if (givenPriorProbsAtRoot != null) {
            // GIVEN!
        }
        else { throw new RuntimeException(); }
    }

    @Override
    public void doIntegrateInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes) {

        // debugging
        System.out.println("esDsAtNode before propagate in t = " + Arrays.toString(esDsAtNode[1]));

        // integrate over birth and death events (low or high resolution inside)
        propagateTInPlace(esDsAtNode, scratchAtNode, dt, lowRes);

        // debugging
        // System.out.println("esDsAtNode after propagate in t and before x = " + Arrays.toString(esDsAtNode[1]));

        // make normal kernel and FFTs it
        populatefY(true, true);

        // integrate over diffusion of quantitative trait
        // propagateXInPlace(esDsAtNode, lowRes);
        propagateXInPlace(esDsAtNode, scratchAtNode, lowRes);

        // debugging
        // System.out.println("esDsAtNode after propagate in t and x = " + Arrays.toString(esDsAtNode[1]));
    }

    @Override
    public void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes) {
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
    public double getLogPFromRelevantObjects() {
        double denormalizationLogLik = 0.0;
        for (double d: logNormalizationFactors) denormalizationLogLik += d;


        return 0.0;
    }

    @Override
    public double calculateLogP() {

        // refreshing parameters
        populateMacroevolParams(false);

        // refreshing tree
        tree = treeInput.get();

        // start recursion for likelihood calculation
        startRecursionAtRootNode(tree.getRoot());

        // put it all together
        double myLogP = getLogPFromRelevantObjects();

        return myLogP;
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

    public double[][][] getScratch(boolean lowRes) {
        if (lowRes) return scratchLo;
        else return scratchHi;
    }

    public double[][] getEsDsAtNode(int nodeIdx, boolean lowRes) {
        if (lowRes) return esDsLo[nodeIdx];
        return esDsHi[nodeIdx];
    }

    public double[][] getScratchAtNode(int nodeIdx, boolean lowRes) {
        if (lowRes) return scratchLo[nodeIdx];
        else return scratchHi[nodeIdx];
    }

    public int getNUsefulTraitBins(boolean lowRes) {
        if (lowRes) return nUsefulXbinsLo;
        else return nUsefulXbinsHi;
    }

    public void setEsDsAtNodeElementAtDim(int nodeIdx, int dimIdx, int eleIdx, double val, boolean lowRes) {
        if (lowRes) esDsLo[nodeIdx][dimIdx][eleIdx] = val;
        else esDsHi[nodeIdx][dimIdx][eleIdx] = val;
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
