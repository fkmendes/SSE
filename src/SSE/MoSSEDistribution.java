package SSE;

import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

import java.util.List;
import java.util.Random;

public class MoSSEDistribution extends QuaSSEDistribution {

    final public Input<Alignment> alignInput = new Input<>("seqAlign", "Molecular alignment.", Input.Validate.REQUIRED);
    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "Molecular site model.", Input.Validate.REQUIRED);

    private int nDimensions = 5;
    private Alignment algn;
    private SiteModel.Base siteModel;
    private SubstitutionModel substModel;
    private MoSSELikelihoodCore likelihoodCore;

    private int nSitePat;
    private double[] instRateMat; // flattened out matrix

    // state that matters for calculateLogP
    private double[][][] esDs, scratch; // first dimension are all nodes, second dimension are E's and D's, third are the bins


    @Override
    public void initAndValidate() {

        // super will (1) populate macroevol parameters, ...
        super.initAndValidate();

        String what2check = "algn";
        algn = alignInput.get();
        areInputsOK(what2check); // checks if tree and MSA are OK

        what2check = "sitemodel";
        areInputsOK(what2check);
        siteModel = (SiteModel.Base) siteModelInput.get();

        nSitePat = algn.getPatternCount();

        int dtRes = 1;
        instRateMat = new double[4*4];

        int nDimensionsFFT = 5; // E, and 4 D's (one for each nuc)
        initializeEsDs(nSitePat, nDimensionsFFT, this.nXbins);
    }

    // needs extra dimension of nSitePat
    public void initializeEsDs(int nSitePat, int nDimensionsFFT, int nXbins) {
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
        super.processBranch(node);
    }

    @Override
    public void processInternalNode() {
    }

    @Override
    public void processRootNode() {

    }

    @Override
    public void doIntegrate(double[][] esDsAtNode, double startTime, boolean isFirstDt, boolean lowRes) {

        // integrate over birth and death events (either at low or high resolution inside)
        propagateT(esDsAtNode, lowRes);

        /*
         * For the step below, I am working on the MoSSELikelihoodCore
         * class.
         *
         * This class will never deal with partials at internal nodes as
         * LikelihoodCore does, because internal nodes are dealt with
         * by the MoSSE class, which just does the D * D product.
         *
         * Instead, this class will update a single partial (per bin) along
         * a branch, at every dt. P(dt,r) will be constant, and the
         * partial update function should take the partial as input and
         * just update it in place.
         *
         * Then I'd need two separate functions for updating partials:
         * (1) For the very first dt of a leaf (where missing data is dealt with)
         * (2) For all other dts
         */

        // MoSSE: integrate over substitution events
        propagateSubst(esDsAtNode, startTime, dt, isFirstDt, lowRes);

        // integrate over diffusion of substitution rate
        propagateX(lowRes);
    }

    @Override
    public void propagateT(double[][] esDsAtNode, boolean lowRes) {
        super.propagateT(esDsAtNode, lowRes);
    }

    @Override
    public void propagateX(boolean lowRes) {
        super.propagateX(lowRes);
    }

    public void propagateSubst(double[][] esDsAtNode, double startTime, double aDt, boolean isFirstDt, boolean lowRes) {

        double endTime = startTime + aDt;
        double[] xRuler;
        if (lowRes) xRuler = xLo;
        else xRuler = xHi;

            for (int i=0; i<nXbins; i++) {
                double xValue = xRuler[i];
                double[] esDsAtNodeAtBin = esDsAtNode[i];
                substModel.getTransitionProbabilities(null, startTime, endTime, xValue, instRateMat); // P(dt,r), where r = x = subst rate
                likelihoodCore.updateMatrices(instRateMat); // set Q matrix inside likelihoodCore

                if (isFirstDt) likelihoodCore.updatePartialInPlaceFirstDt(esDsAtNodeAtBin, aDt); // first dt deals with missing data
                else likelihoodCore.updatePartialInPlace(esDsAtNodeAtBin, aDt);
            }

    }

    @Override
    public void convolve() {

    }

    @Override
    public double calculateLogP() {

        populateMacroevolParams(false);

        // iterating over all unique site patterns (multiply by their frequencies)
//        for (int i=0; i<nSitePat; ++i) {
//            // process(root);
//        }

        return 0.0;
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

    private void areInputsOK(String what2Check) {

        // alignment should have same #taxa as tree, and only works with nucleotide alignments
        if (what2Check == "algn") {
            if (alignInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences. Exiting...");
            if (alignInput.get().getMaxStateCount() != 4) throw new IllegalArgumentException("MoSSE only works with multiple nucleotide alignments. Exiting...");
        }

        // site model should have specific type
        if (what2Check == "sitemodel" && (!(siteModelInput.get() instanceof SiteModel.Base))) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base. Exiting...");
        }
    }
}
