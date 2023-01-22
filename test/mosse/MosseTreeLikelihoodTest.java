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
    public void testMosseLikelihoodInit() {
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
        double epsilon = 0.01;

        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta", new RealParameter(betasArray),
                "epsilon", Double.toString(epsilon)
        );
        tipModel.initAndValidate();

        TaxonSet taxonSet = new TaxonSet(alignment);
        int numTraits = 2;
        // trait 0
        TraitSet trait0 = new TraitSet();
        String trait0Values = "t0=1.0, t1=10.0";
        trait0.initByName(
                "traitname", "trait0",
                "taxa", new TaxonSet(alignment),
                "value", trait0Values);
        // trait 1
        TraitSet trait1 = new TraitSet();
        String trait1Values = "t0=15.0, t1=20.0";
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

        double startSubsRate = 1E-10;
        double endSubsRate = 1E-8;
        int numBins = 100;

        MosseTreeLikelihood likelihood = new MosseTreeLikelihood();
        likelihood.initByName(
                "data", alignment,
                "tree", tree,
                "siteModel", siteModel,
                "tipModel", tipModel,
                "treeModel", new MosseDistribution(),
                "traits", traitsList,
                "startSubsRate", Double.toString(startSubsRate),
                "endSubsRate", Double.toString(endSubsRate),
                "numRateBins", Integer.toString(numBins),
                "lambdaFunc", logFunc,
                "muFunc", constFunc
                );
        likelihood.calculateLogP();
    }

    public void testMosseExponentiatedMatrix() {
        // test case from test-mosse2.R
        //   Q <- t(matrix(c(-0.5,0.2,0.1,0.1,
        //                  0.1,-0.8,0.3,0.4,
        //                  0.2,0.2,-0.6,0.2,
        //                  0.3,0.4,0.2,-0.7),4,4))
        // expm(Q)
        double[] expectedTransitionMatrix = {
                0.6319373, 0.1273384, 0.0840574, 0.08617252,
                0.1062949, 0.5146194, 0.1818267, 0.21867370,
                0.1442823, 0.1381027, 0.5872501, 0.13792546,
                0.1970986, 0.2273229, 0.1513365, 0.56178771};
    }


    public void testMosseLikelihoodRootNode() {
        // test root node treatments
    }

    @Test
    public void testMosseLikelihoodOnTree() {
        int numLeaves = 2;
        String[] sequences = {"A", "C"};
        Alignment alignment = getAlignment(numLeaves, sequences);
        String newick = "(t0: 0.5, t1: 0.5);";
        Tree tree = new Tree(newick);

        RealParameter f = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25});
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);

        Double[] relativeRates = new Double[]{
                -0.5, 0.2, 0.1, 0.1,
                0.1, -0.8, 0.3, 0.4,
                0.2, 0.2, -0.6, 0.2,
                0.3, 0.4, 0.2, -0.7};
        RealParameter customRates = new RealParameter(relativeRates);
        RealParameter rates = new RealParameter(new Double[]{0.2, 0.1, 0.1, 0.1, 0.3, 0.4, 0.2, 0.2, 0.2, 0.3, 0.4, 0.2});
        CustomSubstitutionModel substModel = new CustomSubstitutionModel();
        substModel.initByName("frequencies", freqs, "rates", rates, "customRates", customRates);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", "1.0",
                "gammaCategoryCount", 1,
                "substModel", substModel);

        double startTime = 1;
        double endTime = 0;
        double rate = 1;

        int len = substModel.getStateCount();
        double[] transitionProbMatrix = new double[len*len];
        // testing transition probability exp(Q * t)
        substModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, transitionProbMatrix, false);
        System.out.println("P = exp(Q * t)");
        int count = 1;
        for (double x: transitionProbMatrix) {
            System.out.print(x + " ");
            if (count % 4 == 0) {
                System.out.println();
            }
            count++;
        }
        // test case from test-mosse2.R
        //   Q <- t(matrix(c(-0.5,0.2,0.1,0.1,
        //                  0.1,-0.8,0.3,0.4,
        //                  0.2,0.2,-0.6,0.2,
        //                  0.3,0.4,0.2,-0.7),4,4))
        // expm(Q)
        double[] expectedTransitionMatrix = {
                0.6319373, 0.1273384, 0.0840574, 0.08617252,
                0.1062949, 0.5146194, 0.1818267, 0.21867370,
                0.1442823, 0.1381027, 0.5872501, 0.13792546,
                0.1970986, 0.2273229, 0.1513365, 0.56178771};

        assertArrayEquals(transitionProbMatrix, expectedTransitionMatrix, DELTA);

        Double[] betasArray = {0.1, 0.2};
        double epsilon = 0.01;

        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta", new RealParameter(betasArray),
                "epsilon", Double.toString(epsilon)
        );
        tipModel.initAndValidate();

        TaxonSet taxonSet = new TaxonSet(alignment);
        int numTraits = 2;
        // trait 0
        TraitSet trait0 = new TraitSet();
        String trait0Values = "t0=1.0, t1=10.0";
        trait0.initByName(
                "traitname", "trait0",
                "taxa", new TaxonSet(alignment),
                "value", trait0Values);
        // trait 1
        TraitSet trait1 = new TraitSet();
        String trait1Values = "t0=15.0, t1=20.0";
        trait1.initByName(
                "traitname", "trait1",
                "taxa", new TaxonSet(alignment),
                "value", trait1Values);
        List<TraitSet> traitsList = new ArrayList<>(numTraits);
        traitsList.add(trait0);
        traitsList.add(trait1);

        // lambda and mu functions
        // logistic function
        Double[] y0 = new Double[] { 0.0 };
        Double[] y1 = new Double[] { 0.1 };
        Double[] x0 = new Double[] { 0.1 }; // xmid
        Double[] r = new Double[] { 0.01 };

        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);
        LogisticFunction logFunc = new LogisticFunction();
        logFunc.initByName( "curveYBaseValue", y0rp,
                "curveMaxY", y1rp,
                "sigmoidMidpoint", x0rp,
                "logisticGrowthRate", rrp);
        // constant function
        Double[] yValue = new Double[] { 0.01 }; // const value
        RealParameter yValueRP = new RealParameter(yValue);
        ConstantLinkFn constFunc = new ConstantLinkFn();
        constFunc.initByName("yV", yValueRP);

        double startSubsRate = 0.0164; // TODO: double check value from test case test-mosse2.R (log scale?)
        double endSubsRate = 0.3932;
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
                "endSubsRate", Double.toString(endSubsRate),
                "numRateBins", Integer.toString(numBins),
                "lambdaFunc", logFunc,
                "muFunc", constFunc
        );
        likelihood.initAndValidate();
        likelihood.calculateLogPFull();
    }
}
