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
        int rootIdx = rootNode.getNr();
        double rootHeight = rootNode.getHeight();

        // start recursion
        processInternalNode(rootNode);

        // we're done pruning, let's deal with the prior probs at root now
        double[][] esDsAtRoot;
        double dxAtRightRes;
        int nXBinsRightRes, nUsefulXBinsRightRes;
        if (rootHeight > tc) {
            esDsAtRoot = esDsLo[rootIdx];
            dxAtRightRes = dXbin;
            nUsefulXBinsRightRes = nUsefulXbinsLo;
            nXBinsRightRes = nXbinsLo;
            got2LowRes = true;
        } else {
            esDsAtRoot = esDsHi[rootIdx];
            dxAtRightRes = dXbin / hiLoRatio;
            nUsefulXBinsRightRes = nUsefulXbinsHi;
            nXBinsRightRes = nXbinsHi;
            got2LowRes = false;
        }

        if (priorProbsAtRoot == null) priorProbsAtRoot = new double[nUsefulXBinsRightRes]; // if root prior probs not given, initialize it!

        populatePriorProbAtRoot(esDsAtRoot[1], dxAtRightRes, nUsefulXBinsRightRes, rootPriorType); // ok, prior probs are set!
    }

    @Override
    public void processInternalNode(Node aNode) {

        // debugging
        System.out.println("\nDoing node " + aNode.getID());

        // recur if internal node or sampled ancestor
        if (!aNode.isLeaf()) {
            for (Node childNode: aNode.getChildren()) {
                processInternalNode(childNode); // recur
            }

            // after recursion, prepare initial conditions for integrating this node
            mergeChildrenNodes(aNode);
        }

        // unless we're at the root, now we integrate this branch
        if (!aNode.isRoot()) processBranch(aNode, true);
    }

    @Override
    public void processBranch(Node aNode, boolean normalize) {
        boolean lowRes = false;
        int nodeIdx = aNode.getNr();

        // dealing with branch lengths
        double startTime = aNode.getHeight(); // we're going backwards in time, toward the root
        double endTime = aNode.getParent().getHeight();
        double branchLength2Integrate = endTime - startTime;

        /*
         * Sorting out E's and D's prior to integration
         * (1) Grabbing esDs and scratch depending on resolution
         * (2) Calculating normalization factor
         * (3) Storing normalization factor to unnormalize at the very end of pruning
         * (4) Normalize DsThis is done only once if node is already at low, rather than
         * at every loop of the while block below
         */
        double[][] esDsAtNode, scratchAtNode;

//        // option 0: (rare) branch starts at tc (all low-res)
//        if (startTime == tc) {
//            for (int ithDim=0; ithDim<nDimensions; ithDim++)
//                SSEUtils.hiToLoTransferInPlace(esDsHi[nodeIdx][ithDim], esDsLo[nodeIdx][ithDim], hiLoIdxs4Transfer);
//            esDsAtNode = esDsLo[nodeIdx];
//            scratchAtNode = scratchLo[nodeIdx];
//        }
        // option 1: entire branch might already be > tc (all low-res)
        // note that here the high to low res transfer should have already happened
        if (startTime > tc) {
            esDsAtNode = esDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];

            // debugging
            System.out.println("OPTION 1: E's size should be 2 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor
            if (normalize) logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin); // normalize D's and returns factor, which we record

            // now integrate whole branch
            integrateLength(esDsAtNode, scratchAtNode, branchLength2Integrate, dynamicallyAdjustDt, dtMax, true);
        }
        // option 2: entire branch < tc (all high-res)
        else if ((startTime + branchLength2Integrate) < tc) {
            esDsAtNode = esDsHi[nodeIdx];
            scratchAtNode = scratchHi[nodeIdx];

            // debugging
            System.out.println("OPTION 2: E's size should be 4 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor
            if (normalize) logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin/hiLoRatio); // normalize D's and returns factor, which we record

            // now integrate whole branch
            integrateLength(esDsAtNode, scratchAtNode, branchLength2Integrate, dynamicallyAdjustDt, dtMax, false);
        }
        // option 3: tc happens inside branch (tip-end part in high-res, root-end part in low-res)
        else {
            esDsAtNode = esDsHi[nodeIdx];
            scratchAtNode = scratchHi[nodeIdx];

            // debugging
            System.out.println("OPTION 3: E's size should be 2 * 4 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at high-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor
            if (normalize) logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin/hiLoRatio); // normalize D's and returns factor, which we record

            // high res part
            double lenHi = branchLength2Integrate - tc;

            // debugging
            System.out.println("high-res part, lenHi = " + lenHi);

            integrateLength(esDsAtNode, scratchAtNode, lenHi, dynamicallyAdjustDt, dtMax, false);

            // normalize and record normalization factor
            logNormalizationFactors[nodeIdx] += calculateNormalizationFactor(esDsAtNode[1], dXbin/hiLoRatio);

            // transferring high-res esDs to low-res esDs
            for (int ithDim=0; ithDim<nDimensions; ithDim++)
                SSEUtils.hiToLoTransferInPlace(esDsHi[nodeIdx][ithDim], esDsLo[nodeIdx][ithDim], hiLoIdxs4Transfer);
            esDsAtNode = esDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];

            // low res part
            double lenLo = branchLength2Integrate - lenHi;

            // debugging
            System.out.println("low-res part, lenLo = " + lenLo);

            integrateLength(esDsAtNode, scratchAtNode, lenLo, dynamicallyAdjustDt, dtMax, true);
        }
    }

    @Override
    public void integrateLength(double[][] esDsAtNode, double[][] scratchAtNode, double aLength, boolean dynamicallyAdjust, double maxDt, boolean lowRes) {
        // dealing with dt
        double dt, nIntervals, nonIntegratedDt;
        if (dynamicallyAdjust) {
            nIntervals = Math.ceil(aLength / maxDt);
            dt = aLength / nIntervals; // dynamically adjusted dt
        } else {
            nIntervals = Math.floor(aLength / maxDt);
            dt = maxDt;
            nonIntegratedDt = aLength % maxDt;
        }

        // integrating!
        for (int i=0; i<nIntervals; i++) {
            doIntegrateInPlace(esDsAtNode, scratchAtNode, dt, lowRes);
        }
    }

    @Override
    public double normalizeDs(double[] dsAtNode, double dxAtRightRes) {
        double logNormalizationFactorFromDs;
        logNormalizationFactorFromDs = calculateNormalizationFactor(dsAtNode, dxAtRightRes);

        for (int i = 0; i < dsAtNode.length; i++) {
            dsAtNode[i] /= logNormalizationFactorFromDs;
            // there will be an additional step of normalization when in option 3
            // as a result of integrating the tip-end part in high-res
        }

        return logNormalizationFactorFromDs;
    }

    @Override
    public void mergeChildrenNodes(Node aNode) {
        int nodeIdx = aNode.getNr();
        double nodeHeight = aNode.getHeight();

        List<Node> childrenNode = aNode.getChildren();
        Node leftChild = childrenNode.get(0);
        int leftChildIdx = leftChild.getNr();
        Node rightChild = childrenNode.get(1);
        int rightChildIdx = rightChild.getNr();

        double[][][] esDs;
        double[] birthRatesAtRightRes;
        int nUsefulXbinAtRightRes;
        if (nodeHeight < tc) {
            esDs = esDsHi;
            nUsefulXbinAtRightRes = nUsefulXbinsHi;
            birthRatesAtRightRes = birthRatesHi;
        } else {
            esDs = esDsLo;
            nUsefulXbinAtRightRes = nUsefulXbinsLo;
            birthRatesAtRightRes = birthRatesLo;
        }

        // debugging
        System.out.println("\nMerging children");
        System.out.println("E's before merging = " + Arrays.toString(esDs[leftChildIdx][0]));
        System.out.println("D's left before merging = " + Arrays.toString(esDs[leftChildIdx][1]));
        System.out.println("D's right before merging = " + Arrays.toString(esDs[rightChildIdx][1]));
        System.out.println("Birth rates = " + Arrays.toString(birthRatesAtRightRes));

        // proper merging
        for (int i=0; i<(esDs[nodeIdx][0].length/2); ++i) {
            if (i < nUsefulXbinAtRightRes) {
                esDs[nodeIdx][0][i] = esDs[leftChildIdx][0][i]; // E's
                esDs[nodeIdx][1][i] = esDs[leftChildIdx][1][i] * esDs[rightChildIdx][1][i] * birthRatesAtRightRes[i]; // D's
            } else {
                esDs[nodeIdx][0][i] = 0.0; // rest is 0's
                esDs[nodeIdx][1][i] = 0.0;
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

        normalizationFactorFromDs *= binSize;
        return normalizationFactorFromDs;
    }

    @Override
    public void populatePriorProbAtRoot(double[] dsAtRoot, double dxAtRightRes, int nXbinAtRightRes, String aRootPriorType) {
        // making sure things have the right dimensions (dividing dsAtRoot by two because it carries real and imaginary parts)
        if (providedPriorAtRoot && got2LowRes && ((priorProbsAtRoot.length / hiLoRatio) != nUsefulXbinsLo))
            throw new RuntimeException("ERROR: You need to specify " + nUsefulXbinsHi + " prior probability values, but you specified " + priorProbsAtRoot.length + ". Exiting...");
        else if (providedPriorAtRoot && !got2LowRes && (priorProbsAtRoot.length != nUsefulXbinsHi)) {
            System.out.println("nUsefulXbinsHi = " + nUsefulXbinsHi);
            System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);
            throw new RuntimeException("ERROR: The number of prior probabilities did not match the number of discrete quantitative bins at high resolution. Exiting...");
        }

        System.out.println("rootPriorType = " + aRootPriorType);
        System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);

        // now we populate prior probs
        if (aRootPriorType.equals(FLAT)) {
            for (int i=0; i<priorProbsAtRoot.length; ++i) {
                priorProbsAtRoot[i] = 1.0 / ((nXbinAtRightRes - 1) * dxAtRightRes);
            }
        }
        else if (aRootPriorType.equals(OBS)) {
            double sumOfDsAtRoot = 0.0;
            for (int i=0; i<priorProbsAtRoot.length; ++i) {
                sumOfDsAtRoot += dsAtRoot[i];
            }

            for (int i=0; i<priorProbsAtRoot.length; ++i) {
                priorProbsAtRoot[i] = dsAtRoot[i] / (sumOfDsAtRoot * dxAtRightRes);
            }
        }
        else { throw new RuntimeException("ERROR: You specified an invalid prior probability distribution for the root. Exiting..."); }
    }

    @Override
    public void doIntegrateInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes) {

        // debugging
        System.out.println("esAtNode before propagate in t = " + Arrays.toString(esDsAtNode[0]));
        System.out.println("dsAtNode before propagate in t = " + Arrays.toString(esDsAtNode[1]));

        // integrate over birth and death events (low or high resolution inside)
        propagateTInPlace(esDsAtNode, scratchAtNode, dt, lowRes);

        // debugging
        System.out.println("esAtNode after propagate in t and before x = " + Arrays.toString(esDsAtNode[0]));
        System.out.println("dsAtNode after propagate in t and before x = " + Arrays.toString(esDsAtNode[1]));

        // make normal kernel and FFTs it
        populatefY(true, true);

        // integrate over diffusion of quantitative trait
        // debugging
        // System.out.println("D's length before propagating in X = " + esDsAtNode[1].length);
        propagateXInPlace(esDsAtNode, scratchAtNode, lowRes);

        // debugging
        System.out.println("esAtNode after propagate in t and x = " + Arrays.toString(esDsAtNode[0]));
        System.out.println("dsAtNode after propagate in t and x = " + Arrays.toString(esDsAtNode[1]));
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

        // debugging
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
    public double getLogPFromRelevantObjects(double[][] esDsAtRoot, double sumOfLogNormalizationFactors, double[] birthRates, double dXAtRightRes) {

        // debugging
        System.out.println("\nProcessing quantitites at root");

        double[] esAtRoot = esDsAtRoot[0];
        double[] dsAtRoot = esDsAtRoot[1];

        // debugging
        System.out.println("priorProbsAtRoot = " + Arrays.toString(priorProbsAtRoot));
        System.out.println("esAtRoot = " + Arrays.toString(esAtRoot));
        System.out.println("dsAtRoot = " + Arrays.toString(dsAtRoot));
        System.out.println("esAtRoot.length = " + dsAtRoot.length);
        System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);
        System.out.println("birthRates = " + Arrays.toString(birthRates));
        System.out.println("birthRates.length = " + birthRates.length);

        // conditioning on survival of lineages
        double denomSumForConditioning = 0.0;
        // priorProbsAtRoot.length should be = nUsefulXBins at the right resolution
        for (int i=0; i<priorProbsAtRoot.length; ++i)
            denomSumForConditioning += (priorProbsAtRoot[i] * birthRates[i] * Math.pow(1 - esAtRoot[i], 2));

        // priorProbsAtRoot.length should be = nUsefulXBins at the right resolution
        for (int i=0; i<priorProbsAtRoot.length; ++i) {
            dsAtRoot[i] = dsAtRoot[i] / denomSumForConditioning * dXAtRightRes;
        }

        // debugging
        System.out.println("dsAtRoot after conditioning on survival = " + Arrays.toString(dsAtRoot));

        // incorporating prior probabilities at root
        double prodSumPriorDs = 0.0;
        for (int i=0; i<priorProbsAtRoot.length; ++i) {
            prodSumPriorDs += (priorProbsAtRoot[i] * dsAtRoot[i]);
        }

        double thisLogLik = Math.log(prodSumPriorDs * dXAtRightRes) + sumOfLogNormalizationFactors; // denormalizing by adding sumOfLogNormalizationFactors

        return thisLogLik;
    }

    @Override
    public double calculateLogP() {

        // refreshing parameters
        populateMacroevolParams(false);

        // refreshing tree and root
        tree = treeInput.get();
        Node rootNode = tree.getRoot();
        int rootIdx = rootNode.getNr();

        // start recursion for likelihood calculation
        startRecursionAtRootNode(rootNode);

        // put it all together
        double sumOfLogNormalizationFactors = 0.0;
        for (double d: logNormalizationFactors) sumOfLogNormalizationFactors += d;
        double[][] esDsAtRootAtRightRes;
        double[] birthRatesAtRightRes;
        double dxAtRightRes;
        if (got2LowRes) {
            esDsAtRootAtRightRes = esDsLo[rootIdx]; // low res
            birthRatesAtRightRes = birthRatesLo;
            dxAtRightRes = dtMax;
        }
        else {
            esDsAtRootAtRightRes = esDsHi[rootIdx]; // high res
            birthRatesAtRightRes = birthRatesHi;
            dxAtRightRes = dtMax / hiLoRatio;
        }

        double myLogP = getLogPFromRelevantObjects(esDsAtRootAtRightRes, sumOfLogNormalizationFactors, birthRatesAtRightRes, dxAtRightRes);

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

    // for debugging
    public void setEsDsAtNodeElementAtDim(int nodeIdx, int dimIdx, int eleIdx, double val, boolean lowRes) {
        if (lowRes) esDsLo[nodeIdx][dimIdx][eleIdx] = val;
        else esDsHi[nodeIdx][dimIdx][eleIdx] = val;
    }

    // for debugging
    public void setGot2LowRes() {
        got2LowRes = true;
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
