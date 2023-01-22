package mosse;

import SSE.ConstantLinkFn;
import SSE.LinkFn;
import SSE.LogisticFunction;
import beast.core.Description;
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

import java.lang.UnsupportedOperationException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * @author Kylie Chen
 */

@Description("Mosse likelihood class calculates the probability of sequence and trait data on a tree")
public class MosseTreeLikelihood extends TreeLikelihood {

    protected MosseLikelihoodCore mosseLikelihoodCore;
    // trait data
    final public Input<List<TraitSet>> traitListInput = new Input<>("traits", "list of traits", new ArrayList<>());
    // tip model, species diversification model and trait model
    final public Input<MosseTipLikelihood> tipModelInput = new Input<>("tipModel", "model of tip probabilities", Input.Validate.REQUIRED);
    final public Input<Distribution> treeModelInput = new Input<>("treeModel", "species diversification model", Input.Validate.REQUIRED);
    // substitution rate parameters
    final public Input<RealParameter> startSubsRateInput = new Input<>("startSubsRate", "lower range for substitution rate", Input.Validate.REQUIRED);
    final public Input<RealParameter> endSubsRateInput = new Input<>("endSubsRate", "upper range for substitution rate", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> numRateBinsInput = new Input<>("numRateBins", "number of bins for substitution rate", Input.Validate.REQUIRED);
    final public Input<LinkFn> lambdaFuncInput = new Input<>("lambdaFunc", "function for birth rate lambda", Input.Validate.REQUIRED);
    final public Input<LinkFn> muFuncInput = new Input<>("muFunc", "function for death rate mu", Input.Validate.REQUIRED);

//    final public int SUBST_NUM_STATES = 4; // for testing

    protected LinkFn lambdaFunc;
    protected LinkFn muFunc;
    protected List<TraitSet> traits;
    protected MosseTipLikelihood tipModel;
    protected MosseDistribution treeModel;
    protected double startSubsRate;
    protected double endSubsRate;
    protected int numRateBins;

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
     * calculate log P without caching
     * @return log P
     */
    public double calculateLogPFull() {
        final TreeInterface tree = treeInput.get();

        traverseFull(tree.getRoot());

        return 0.0; // update using root partials
    }

    private void traverseFull(Node node) {
        int numPlan = 5; // dimensions
        double deltaT = 0.001; // dt
        double rate = 1.0; // dx
        int padLeft = treeModel.padLeftInput.get().getValue();
        int padRight = treeModel.padRightInput.get().getValue();
        int numEntries = numRateBins - padLeft - padRight - 1; // num non-zero entries (length of lambda and mu)
        int numStates = dataInput.get().getDataType().getStateCount();
        int numPattern = dataInput.get().getPatternCount();

        double[] transitionMatrix = new double[numStates * numStates];
        double[][] transitionMatrices = new double[numEntries][transitionMatrix.length];
        // P(0) = exp(dx * Q * dt)
        substitutionModel.getTransitionProbabilities(node, 0, deltaT, rate, transitionMatrix);
        transitionMatrices[0] = transitionMatrix;
        // update transitionMatrices
        for (int i = 1; i < numEntries; i++) {
            double[] prevMatrix = transitionMatrices[i - 1];
            transitionMatrices[i] = matrixPow(prevMatrix, 2);
        }
        flatTransitionMatrices = Arrays.stream(transitionMatrices)
                .flatMapToDouble(Arrays::stream)
                .toArray();

        // get lambdas and mus
        double[] x = getSubstitutionRates(numEntries); // substitution rates
        double[] lambda = new double[numEntries];
        double[] mu = new double[numEntries];
        lambda = lambdaFunc.getY(x, lambda, true);
        mu = muFunc.getY(x, mu, true);

        if (node.isLeaf()) {
            // leaf node
            // leaf partials size = nrPatterns * (nrStates + 1) * numBins
            // columns = (4x D's for each nucleotide, 1x E), rows = bins for substitution rate
            setPartials(node, numPattern); // all site patterns

        } else if (!node.isRoot()) {
            // internal node
            traverseFull(node.getLeft());
            traverseFull(node.getRight());

            double branchTimeLeft = node.getHeight() - node.getLeft().getHeight();
            double branchTimeRight = node.getHeight() - node.getRight().getHeight();

            double[] patternPartialsLeft = new double[numPattern * numPlan * numRateBins];
            double[] patternPartialsRight = new double[numPattern * numPlan * numRateBins];
            double[] partialsAllPatterns = new double[numPattern * numPlan * numRateBins];

            // get child node partials all patterns
            mosseLikelihoodCore.getNodePartials(node.getLeft().getNr(), patternPartialsLeft);
            mosseLikelihoodCore.getNodePartials(node.getRight().getNr(), patternPartialsRight);
            int k = 0;
            for (int pattern = 0; pattern < numPattern; pattern++) {
                // partial for single pattern
                int partialSize = numPlan * numRateBins;
                int startPos = pattern * partialSize;
                double[] partialsLeft =  new double[numPlan * numRateBins];
                System.arraycopy(patternPartialsLeft, startPos, partialsLeft, 0, partialSize);
                double[] partialsRight = new double[numPlan * numRateBins];
                System.arraycopy(patternPartialsRight, startPos, partialsRight, 0, partialSize);
                double[] partialsCombined = new double[partialsLeft.length];

                // propagate each child branch
                treeModel.calculateBranchLogP(branchTimeLeft, partialsLeft, lambda, mu, flatTransitionMatrices, partialsLeft);
                treeModel.calculateBranchLogP(branchTimeRight, partialsRight, lambda, mu, flatTransitionMatrices, partialsRight);

                // assumes t less than tc threshold
                for (int i = 0; i < numRateBins; i++) {
                    double lambdaX = lambda[i]; // birth rate at substitution rate x
                    for (int j = 0; j < numPlan; j++) {
                        int index = i * numPlan + j; // TODO: add unit test for index
                        if (j == 0) {
                            // E is topology independent
                            partialsCombined[index] = partialsLeft[index]; // for testing
                            partialsAllPatterns[k] = partialsLeft[index];
                        } else {
                            // D_left * D_right * lambda(x)
                            partialsCombined[index] = partialsLeft[index] * partialsRight[index] * lambdaX; // for testing
                            partialsAllPatterns[k] = partialsLeft[index] * partialsRight[index] * lambdaX;
                        }
                        k++;
                    }
                }
            }
            // set internal node partials for all patterns
            mosseLikelihoodCore.setNodePartials(node.getNr(), partialsAllPatterns);

        } else {
            // root node
            traverseFull(node.getLeft());
            traverseFull(node.getRight());

            // propagate child branches of root


            // do calculation at root
        }
    }

    public double[] getFlatTransitionMatrices() {
        return flatTransitionMatrices;
    }

    private double[] getSubstitutionRates(int numEntries) {
        double[] res = new double[numEntries];
        // test case values
        // startSubsRate = 0.0164
        // endSubsRate = 0.3932
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

//        try {
//            if (traverse(tree.getRoot()) != Tree.IS_CLEAN) {
//                calcLogP();
//            }
//        }
//        catch (ArithmeticException e) {
//            return Double.NEGATIVE_INFINITY;
//        }

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
