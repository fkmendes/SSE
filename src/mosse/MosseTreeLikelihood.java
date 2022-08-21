package mosse;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("Mosse likelihood class calculates the probability of sequence and trait data on a tree")
public class MosseTreeLikelihood extends TreeLikelihood {

    // trait data
    final public Input<List<TraitSet>> traitListInput = new Input<>("traits", "list of traits", new ArrayList<>());
    // tip model, species diversification model and trait model
    final public Input<MosseTipLikelihood> tipModelInput = new Input<>("tipModel", "model of tip probabilities", Input.Validate.REQUIRED);
    final public Input<Distribution> treeModelInput = new Input<>("treeModel", "species diversification model", Input.Validate.REQUIRED);
    final public Input<Distribution> traitModelInput = new Input<>("traitModel", "model of trait evolution (e.g., brownian motion)", Input.Validate.REQUIRED);
    // substitution rate parameters
    final public Input<RealParameter> startSubsRateInput = new Input<>("startSubsRate", "lower range for substitution rate", Input.Validate.REQUIRED);
    final public Input<RealParameter> endSubsRateInput = new Input<>("endSubsRate", "upper range for substitution rate", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> numRateBinsInput = new Input<>("numRateBins", "number of bins for substitution rate", Input.Validate.REQUIRED);

    protected List<TraitSet> traits;
    protected MosseTipLikelihood tipModel;
    protected Distribution treeModel;
    protected Distribution traitModel;
    protected double startSubsRate;
    protected double endSubsRate;
    protected int numRateBins;

    @Override
    public void initAndValidate() {
        traits = traitListInput.get();
        tipModel = tipModelInput.get();
        treeModel = treeModelInput.get();
        traitModel = traitModelInput.get();

        startSubsRate = startSubsRateInput.get().getValue();
        endSubsRate = endSubsRateInput.get().getValue();
        numRateBins = numRateBinsInput.get().getValue();

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
        }

        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();

        // set likelihood core number of states and number of rate bins
        likelihoodCore = new MosseLikelihoodCore(stateCount, numRateBins);

        String className = getClass().getSimpleName();
        Alignment alignment = dataInput.get();

        // logging likelihood class
        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
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
     * set leaf partials using tip likelihood model *
     */
    @Override
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            String taxonName = node.getID();
            // get traits
            double[] traitValues = new double[traits.size()];
            for (int i = 0; i < traits.size(); i++) {
                double traitValue = traits.get(i).getValue(taxonName);
                traitValues[i] = traitValue;
            }
            double[] partials = new double[patternCount * states * numRateBins];
            int k = 0;
            int taxonIndex = data.getTaxonIndex(node.getID());
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                double[] tipLikelihoods = tipModel.getTipLikelihoods(traitValues, numRateBins, startSubsRate, endSubsRate);
                int stateCount = data.getPattern(taxonIndex, patternIndex);
                    boolean[] stateSet = data.getStateSet(stateCount);
                    for (int state = 0; state < states; state++) {
                        if (stateSet[state] == true) {
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
            likelihoodCore.setNodePartials(node.getNr(), partials);

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
        throw new NotImplementedException();
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
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
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    public double calculateLogP() {
        return 0.0;
    }

    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (tipModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
            // update leaf partials from tip model
        }

        return treeInput.get().somethingIsDirty();
    }
}
