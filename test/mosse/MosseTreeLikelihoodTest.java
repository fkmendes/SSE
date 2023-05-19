package mosse;

import SSE.ConstantLinkFn;
import SSE.LogisticFunction;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GTR;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

import org.junit.Test;
import test.beast.evolution.substmodel.GTRTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNotNull;


/**
 * @author Kylie Chen
 */

public class MosseTreeLikelihoodTest {

    private static double DELTA = 1e-7;

    /**
     * returns an Alignment of nucleotide sequences
     * @param numLeaves number of taxa (or leaves in tree)
     * @param sequences nucleotide sequences for each taxa
     * @return an Alignment of nucleotide sequences
     */
    public Alignment getAlignment(int numLeaves, String[] sequences) {
        List<Sequence> seqList = new ArrayList<>();

        for (int i = 0; i < numLeaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, sequences[i]));
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");

        return alignment;
    }

    /**
     * tests initialization of MosseTreeLikelihood does not throw any errors
     * using a simple two taxa tree (t0: 0.5, t1: 0.5) with tips "A" and "C"
     */
    @Test
    public void testMosseLikelihood() {
        int numLeaves = 2;
        String[] sequences = {"A", "C"};
        Alignment alignment = getAlignment(numLeaves, sequences);
        String newick = "(t0: 0.5, t1: 0.5);";
        Tree tree = new Tree(newick);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", "1.0",
                "gammaCategoryCount", 1,
                "substModel", JC);

        Double[] betasArray = {0.1, 0.2};
        double meanSubst = 0.01; // mean substitution rate
        double epsilon = 0.01;

        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta", new RealParameter(betasArray),
                "subst", Double.toString(meanSubst),
                "epsilon", Double.toString(epsilon),
                "logscale", "false"
        );
        tipModel.initAndValidate();

        TaxonSet taxonSet = new TaxonSet(alignment);
        int numTraits = 2;
        // trait 0
        TraitSet trait0 = new TraitSet();
        String trait0Values = "t0=0.3, t1=0.3";
        trait0.initByName(
                "traitname", "trait0",
                "taxa", new TaxonSet(alignment),
                "value", trait0Values);
        // trait 1
        TraitSet trait1 = new TraitSet();
        String trait1Values = "t0=0.15, t1=0.1";
        trait1.initByName(
                "traitname", "trait1",
                "taxa", new TaxonSet(alignment),
                "value", trait1Values);
        List<TraitSet> traitsList = new ArrayList<>(numTraits);
        traitsList.add(trait0);
        traitsList.add(trait1);

        // lambda and mu functions
        // logistic
        Double[] x0 = new Double[] { 0.0 };
        Double[] y1 = new Double[] { 0.2 };
        Double[] y0 = new Double[] { 0.1 };
        Double[] r = new Double[] { 2.5 };
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);
        LogisticFunction logFunc = new LogisticFunction();
        logFunc.initByName( "curveYBaseValue",
                y0rp, "curveMaxY", y1rp,
                "sigmoidMidpoint", x0rp,
                "logisticGrowthRate", rrp);
        // constant
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValueRP = new RealParameter(yValue);
        ConstantLinkFn constFunc = new ConstantLinkFn();
        constFunc.initByName("yV", yValueRP);

        double startSubsRate = 0.01;
        int numBins = 1024;

        MosseTreeLikelihood likelihood = new MosseTreeLikelihood();
        likelihood.initByName(
                "data", alignment,
                "tree", tree,
                "siteModel", siteModel,
                "tipModel", tipModel,
                "treeModel", new MosseDistribution(),
                "traits", traitsList,
                "startSubsRate", Double.toString(startSubsRate),
                "numRateBins", Integer.toString(numBins),
                "lambdaFunc", logFunc,
                "muFunc", constFunc
                );

        // using observed root
        double result = likelihood.calculateLogP();

        assert !Double.isNaN(result);
        assert !Double.isInfinite(result);

        System.out.println("testMosseLikelihood logP = " + result);
    }

    // test matrix exponentiation from test-mosse2.R
    //   Q <- matrix(c(-0.6,0.2,0.1,0.1,
    //                  0.1,-0.8,0.3,0.4,
    //                  0.2,0.2,-0.6,0.2,
    //                  0.3,0.4,0.2,-0.7),4,4)
    // expm(Q)
    @Test
    public void testMosseExponentiatedMatrix() {
        RealParameter f = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25});
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);

        Double[] qMatrix = new Double[]{
                -0.6, 0.1, 0.2, 0.3,
                0.2, -0.8, 0.2, 0.4,
                0.1, 0.3, -0.6, 0.2,
                0.1, 0.4, 0.2, -0.7};

        RealParameter customRates = new RealParameter(qMatrix);
        CustomSubstitutionModel substModel = new CustomSubstitutionModel();
        substModel.initByName("frequencies", freqs, "customRates", customRates);

        double startTime = 1;
        double endTime = 0;
        double rate = 1;

        int len = substModel.getStateCount();
        double[] transitionProbMatrix = new double[len * len];

        // testing transition probability exp(Q * t)
        substModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, transitionProbMatrix, false);

        double[] expectedTransitionMatrix = {
                0.57265280, 0.1018599, 0.1376678, 0.1878195,
                0.12128431, 0.5143445, 0.1376678, 0.2267034,
                0.08036085, 0.1816673, 0.5869967, 0.1509751,
                0.08240038, 0.2185117, 0.1376678, 0.5614202};

        assertArrayEquals(transitionProbMatrix, expectedTransitionMatrix, DELTA);
    }


    public void testMosseLikelihoodRoot() {
        // test root node treatments
    }
}
