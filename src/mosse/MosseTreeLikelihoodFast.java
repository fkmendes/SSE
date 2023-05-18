package mosse;

import SSE.LinkFn;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Kylie Chen
 */
public class MosseTreeLikelihoodFast extends TreeLikelihood {


    protected MosseLikelihoodCore mosseLikelihoodCore;
    // trait data
    final public Input<List<TraitSet>> traitListInput = new Input<>("traits", "list of traits", new ArrayList<>());
    // tip model, species diversification model and trait model
    final public Input<MosseTipLikelihood> tipModelInput = new Input<>("tipModel", "model of tip probabilities", Input.Validate.REQUIRED);
    final public Input<Distribution> treeModelInput = new Input<>("treeModel", "species diversification model", Input.Validate.REQUIRED);
    // substitution rate parameters

    // TODO: remove redundant parameter in substitution rate inputs
    final public Input<RealParameter> startSubsRateInput = new Input<>("startSubsRate", "lower range for substitution rate", Input.Validate.REQUIRED);
    final public Input<RealParameter> endSubsRateInput = new Input<>("endSubsRate", "upper range for substitution rate", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> numRateBinsInput = new Input<>("numRateBins", "number of bins for substitution rate", Input.Validate.REQUIRED);
    final public Input<LinkFn> lambdaFuncInput = new Input<>("lambdaFunc", "function for birth rate lambda", Input.Validate.REQUIRED);
    final public Input<LinkFn> muFuncInput = new Input<>("muFunc", "function for death rate mu", Input.Validate.REQUIRED);

    final public Input<LinkFn> rootFuncInput = new Input<>("rootFunc", "function for root", Input.Validate.OPTIONAL);

    final public Input<IntegerParameter> rootOptionInput = new Input<>("rootOption", "option for root calculation", Input.Validate.OPTIONAL);

    final public Input<IntegerParameter> resolutionOptionInput = new Input<>("resolution", "resolution scale factor", new IntegerParameter("1"), Input.Validate.OPTIONAL);


//    final public int SUBST_NUM_STATES = 4; // for testing

    // root options
    final public int ROOT_FLAT = 1;
    final public int ROOT_OBS = 2;
    final public int ROOT_EQUI = 3;
    final public int ROOT_GIVEN = 4;

    protected int resolution;

    protected int rootOption;

    protected LinkFn rootFunc;

    protected LinkFn lambdaFunc;
    protected LinkFn muFunc;
    protected List<TraitSet> traits;
    protected MosseTipLikelihood tipModel;
    protected MosseDistribution treeModel;
    protected double startSubsRate;
    protected double endSubsRate;
    protected int numRateBins;

    protected double[] lambdas;

    protected double[] mus;

    double[] flatTransitionMatrices;

    private boolean updateTips = true;

    private boolean updateSiteModel = true;

    @Override
    public void initAndValidate() {
        traits = traitListInput.get();
        tipModel = tipModelInput.get();
        treeModel = (MosseDistribution) treeModelInput.get();

        startSubsRate = startSubsRateInput.get().getValue();
        endSubsRate = endSubsRateInput.get().getValue();
        numRateBins = numRateBinsInput.get().getValue();

        lambdaFunc = lambdaFuncInput.get();
        muFunc = muFuncInput.get();

        if (rootOptionInput.get() != null) {
            rootOption = rootOptionInput.get().getValue();
        }
        if (rootFuncInput.get() != null) {
            rootFunc = rootFuncInput.get();
        }
        if (resolutionOptionInput.get() != null) {
            resolution = resolutionOptionInput.get().getValue();
        }

        // input checking
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            String message = String.format(
                    "The number of leaves in tree (%d) does not match the number of sequences (%d).",
                    treeInput.get().getLeafNodeCount(),
                    dataInput.get().getTaxonCount());
            throw new IllegalArgumentException(message);
        } else if (numRateBins <= 0) {
            throw new IllegalArgumentException("numRateBins input must be a positive integer");
        } else if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        } else if (branchRateModelInput.get() != null) {
            System.err.println("Ignoring clock model " + branchRateModelInput.get().getID());
        }

        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        // remove requirement for clock model
        branchRateModelInput.setRule(Input.Validate.OPTIONAL);
        branchRateModel = null;
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();

        // set likelihood core number of states and number of rate bins
        int padLeft = treeModel.padLeftInput.get().getValue();
        int padRight = treeModel.padRightInput.get().getValue();
        mosseLikelihoodCore = new MosseLikelihoodCore(stateCount, numRateBins, padLeft, padRight);

        String className = getClass().getSimpleName();
        Alignment alignment = dataInput.get();

        // logging likelihood class
        Log.info.println(className + "(" + getID() + ") uses " + mosseLikelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));

        setProportionInvariant(m_siteModel.getProportionInvariant());
        m_siteModel.setPropInvariantIsCategory(false);
        if (getProportionInvariant() > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        // set up likelihood core and initialize partials
        this.initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount * numRateBins];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }

    /**
     * set leaf partials using tip GLM likelihood model *
     */
    @Override
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            double[] traitValues = getTraits(node);
            double[] partials = new double[patternCount * (states + 1) * numRateBins];
            int k = 0;
            int taxonIndex = data.getTaxonIndex(node.getID());
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                double[] tipLikelihoods = tipModel.getTipLikelihoods(traitValues, numRateBins, startSubsRate, endSubsRate);
                int stateCount = data.getPattern(taxonIndex, patternIndex);
                boolean[] stateSet = data.getStateSet(stateCount);
                // E initial values are zero
                for (int i = 0; i < numRateBins; i++) {
                    partials[k++] = 0.0;
                }
                // D initial values
                for (int state = 0; state < states; state++) {
                    if (stateSet[state]) {
                        // set likelihoods for nucleotide in data
                        for (int i = 0; i < numRateBins; i++) {
                            partials[k++] = tipLikelihoods[i];
                        }
                    } else {
                        // otherwise set likelihoods to zero
                        for (int i = 0; i < numRateBins; i++) {
                            partials[k++] = 0.0;
                        }
                    }
                }
            }
            mosseLikelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    /**
     * set leaf states (not applicable for this class, use setPartials instead)
     */
    @Override
    protected void setStates(Node node, int patternCount) {
        throw new UnsupportedOperationException();
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        mosseLikelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true,
                m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        // set up tip partials
        setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());

        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            mosseLikelihoodCore.createNodePartials(extNodeCount + i);
        }
    }


    /**
     *
     * @param nx number of bins for substitution rate
     * @param dx distance between xs
     * @param r for resolution scale factor
     * @param QMatrix Q instantaneous rate matrix
     * @param result root node result matrix of D and E values
     * @param conditionSurv whether to condition on survival
     * @return log probability for root
     */
    private double makeRootFuncMosse(int nx, double dx, int r, double[] QMatrix, double[] result, boolean conditionSurv) {
        double dxScaled = dx * r;
        int ntypes = 4;
        double[][] vals = new double[nx][ntypes+1];
        int count = 0;
        for (int j = 0; j < ntypes + 1; j++) { // nucleotide types columns
            for (int i = 0; i < nx; i++) { // nx rows
                vals[i][j] = result[count];
                count++;
            }
        }

        double[][] dRoot = getDValues(vals); // get root D values in last column
        double[] eRoot = null;

        double[] x = getSubstitutionRates(nx);
        // root options
        double[][] rootP = getRootProb(dRoot, x, nx, QMatrix, rootOption, rootFunc, ntypes);

        if (conditionSurv) {
            eRoot = getColumn(vals, 1); // get root E values as a column
            // apply function on dRoot
            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; i < nx; j++) {
                    double lambdaX = lambdas[j];
                    // element-wise division of d column
                    dRoot[i][j] = dRoot[i][j] / (lambdaX * Math.pow((1 - eRoot[j]), 2));
                }
            }
        }

        double[][] rootProduct = getProduct(rootP, dRoot);
        double lq = getSum(dRoot); // lq value is sum(D vector) for numerical underflow

        double logProb = Math.log(getSum(rootProduct) * dxScaled) + lq; // log for numerical underflow

        return logProb;
    }

    private double getColumnSum(double[][] dRoot, int column) {
        double sum = 0.0;
        int nrow = dRoot.length;
        for (int i = 0; i < nrow; i++) {
            sum = sum + dRoot[i][column];
        }
        return sum;
    }

    private double[][] getDValues(double[][] vals) {
        // all columns except first column
        int ntypes = 4;
        int nrow = vals.length;
        int ncol = vals[0].length;
        double[][] dValues = new double[nrow][ntypes];
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                dValues[i][j] = vals[i][j + 1];
            }
        }
        return dValues;
    }

    private double[][] getRootProb(double[][] dRoot, double[] x, int nx, double[] QOrig, int rootOption, LinkFn rootFunc, int ntypes) {
        double dx = x[1] - x[0];
        int numSubstBins = dRoot[0].length;
        double[][] p = new double[dRoot.length][numSubstBins];

        if (rootOption == ROOT_FLAT) {
            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; j < numSubstBins; j++) {
                    p[i][j] = 1 / ((nx - 1) * ntypes * dx);
                }
            }
        } else if (rootOption == ROOT_OBS)  {
            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; j < numSubstBins; j++) {
                    p[i][j] = dRoot[i][j] / (getSum(dRoot) * dx);
                }
            }
        } else {
            double[] rootI = substitutionModel.getFrequencies(); // equilibrium freqs
            if (rootOption == ROOT_EQUI) {
                for (int i = 0; i < ntypes; i++) {
                    for (int j = 0; j < numSubstBins; j++) {
                        p[i][j] = rootI[i] * dRoot[i][j] / (getColumnSum(dRoot, i) * dx); // mapply
                    }
                }
            } else if (rootOption == ROOT_GIVEN){
                double[] y = new double[x.length];
                y = rootFunc.getY(x, y, true);
                for (int i = 0; i < ntypes; i++) {
                    for (int j = 0; j < numSubstBins; j++) {
                        p[i][j] = rootI[i] * y[j]; // mapply
                    }
                }
            }
        }

        return p;
    }

    private double[][] getProduct(double[][] array1, double[][] array2) {
        assert(array1.length == array2.length);
        assert(array1[0].length == array2[0].length);
        double[][] res = new double[array1.length][array1[0].length];
        for (int i = 0; i < res.length; i++) {
            for (int j = 0; j < res[0].length; j++) {
                res[i][j] = array1[i][j] * array2[i][j];
            }
        }
        return res;
    }

    private double getSum(double[] array) {
        double res = 0.0;
        for (double i: array) {
            res = res + i;
        }
        return res;
    }

    private double getSum(double[][] array) {
        double res = 0.0;
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                res = res + array[i][j];
            }
        }
        return res;
    }

    private double[] getColumn(double[][] values, int colIndex) {
        int rows = values.length;
        double[] column = new double[rows];
        if (colIndex == -1) {
            colIndex = values[0].length - 1;
        }
        for (int i = 0; i < rows; i++) {
            column[i] = values[i][colIndex];
        }
        return column;
    }

    public double[] getFlatTransitionMatrices() {
        return flatTransitionMatrices;
    }

    private double[] getSubstitutionRates(int numEntries) {
        double[] res = new double[numEntries];
        res[0] = startSubsRate;
        double interval = (endSubsRate - startSubsRate) / numEntries;
        for(int i = 1; i < numEntries; i++) {
            res[i] = res[i - 1] + interval;
        }
        return res;
    }

    private double[] getTraits(Node node) {
        String taxonName = node.getID();
        double[] traitValues = new double[traits.size()];
        for (int i = 0; i < traits.size(); i++) {
            double traitValue = traits.get(i).getValue(taxonName);
            traitValues[i] = traitValue;
        }
        return traitValues;
    }

    private double[] matrixPow(double[] matrixData, double scale) {
        DoubleMatrix matrix = new DoubleMatrix(matrixData);
        DoubleMatrix result = MatrixFunctions.pow(scale, matrix);
        return result.toArray();
    }

    public double calculateLogP() {
        final TreeInterface tree = treeInput.get();

        // cached version of likelihood
        try {
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN) {
                calcLogP();
            }
        } catch (Exception e) {
            return Double.NEGATIVE_INFINITY;
        }

        return 0.0;
    }

    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }
    }

    /**
     * traverse tree with optimized caching
     * @param node tree node
     * @return update flag
     */
    @Override
    protected int traverse(final Node node) {
        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchTime = node.getLength();

        // udpate lambda and mu
        double[] lambda = new double[1024];
        double[] mu = new double[1024];

        if (node.isLeaf() && (update != Tree.IS_CLEAN) && updateTips) {
            // update tips from GLM if node is a leaf
            updateTips = false;
        }

        if (node.isLeaf() && (update != Tree.IS_CLEAN) && updateSiteModel) {
            // update site transition matrices
            updateSiteModel = false;
        }

        if (!node.isLeaf()) {
            // if node is internal, update the partial likelihoods
            final Node child1 = node.getLeft();
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // if either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                mosseLikelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    mosseLikelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    // update each branch
                    double deltaT = 0.001;
                    double rate = 1.0;
                    int numStates = dataInput.get().getDataType().getStateCount(); // substitutionModel.getStateCount();
                    double[] transitionMatrix = new double[numStates * numStates];
                    double[] vars = new double[1024]; // vars from children nodes
                    double[] result = new double[vars.length];
                    substitutionModel.getTransitionProbabilities(node, 0, deltaT, rate, transitionMatrix);
                    treeModel.calculateBranchLogP(branchTime, vars, lambda, mu, transitionMatrix, result);
                    // update internal node
                    mosseLikelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                }

                if (node.isRoot()) {
                    // root of the beast tree
                    // calculate the pattern likelihoods
                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    mosseLikelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }
                    mosseLikelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }
            }
        }

        return update;
    }

    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            updateSiteModel = true;
            return true;
        }
        if (tipModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            updateTips = true;
            return true;
            // update leaf partials from tip model
        }

        return treeInput.get().somethingIsDirty();
    }


}
