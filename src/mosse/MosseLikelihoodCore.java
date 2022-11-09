package mosse;

import beast.core.Description;
import beast.evolution.likelihood.BeerLikelihoodCore;


/**
 * @author Kylie Chen
 */

@Description("Mosse likelihood core calculation class")
public class MosseLikelihoodCore extends BeerLikelihoodCore {

    protected int numRateBins;
    protected int padLeft;
    protected int padRight;

    protected int lambdaSize;

    public MosseLikelihoodCore(int nrOfStates, int numRateBins, int padLeft, int padRight) {
        super(nrOfStates);
        this.numRateBins = numRateBins;
        this.padLeft = padLeft;
        this.padRight = padRight;
    }

    /**
     * initializes partial likelihood arrays.
     *
     * @param nodeCount           the number of nodes in the tree
     * @param patternCount        the number of patterns
     * @param matrixCount         the number of matrices (i.e., number of categories for gamma rate heterogeneity)
     * @param integrateCategories whether sites are being integrated over all matrices
     * @param useAmbiguities      whether to use ambiguious characters
     */
    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities) {
        this.nrOfNodes = nodeCount;
        this.nrOfPatterns = patternCount;
        this.nrOfMatrices = matrixCount;
        this.integrateCategories = integrateCategories;

        if (integrateCategories) {
            partialsSize = patternCount * (nrOfStates + 1) * matrixCount * numRateBins;
        } else {
            partialsSize = patternCount * (nrOfStates + 1) * numRateBins;
        }

        partials = new double[2][nodeCount][];
        lambdaSize = patternCount * numRateBins - padLeft - padRight - 1;

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        states = new int[nodeCount][];

        for (int i = 0; i < nodeCount; i++) {
            partials[0][i] = null;
            partials[1][i] = null;
            states[i] = null;
        }

        matrixSize = nrOfStates * nrOfStates;

        matrices = new double[2][nodeCount][matrixCount * matrixSize];
    }

    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {
        if (this.partials[0][nodeIndex] == null) {
            this.partials[0][nodeIndex] = new double[partialsSize];
            this.partials[1][nodeIndex] = new double[partialsSize];
        }
        if (partials.length < partialsSize) {
            // rate heterogeneity categories
            int k = 0;
            for (int i = 0; i < nrOfMatrices; i++) {
                System.arraycopy(partials, 0, this.partials[0][nodeIndex], k, partials.length); // update to current cache index
                k += partials.length;
            }
        } else {
            System.arraycopy(partials, 0, this.partials[0][nodeIndex], 0, partials.length); // update to current cache index
        }
    }

}
