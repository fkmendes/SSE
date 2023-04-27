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
        double meanSubst = 0.01; // mean substitution rate
        double epsilon = 0.01;

        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta", new RealParameter(betasArray),
                "subst", Double.toString(meanSubst),
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


    public void testMosseLikelihoodRootNode() {
        // test root node treatments
    }

    private Alignment getAlignmentLarge() {
        return null;
    }

    public void testMosseLikelihoodOnTree() {
        String newick = "(sp2:13.77320255,(sp1:12.76688384,((((sp12:1.170387028,sp13:1.170387028)nd16:0.9837720325,sp9:2.154159061)nd11:5.451401092,((sp5:4.311645343,(sp14:0.8910055279,sp15:0.8910055279)nd14:3.420639815)nd9:2.536663776,((sp16:0.3011866125,sp17:0.3011866125)nd12:4.264383667,(sp6:3.95083843,sp7:3.95083843)nd13:0.6147318498)nd10:2.282738839)nd8:0.7572510339)nd5:2.554739141,((sp10:2.059478202,sp11:2.059478202)nd15:0.4198789018,sp8:2.479357104)nd6:7.68094219)nd4:2.60658455)nd3:1.006318707)nd1;";
        Tree tree = new Tree(newick);

        RealParameter f = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25}); // update pis
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

        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", "1.0",
                "gammaCategoryCount", 1,
                "substModel", substModel);

        Double[] betasArray = {0.1, 0.2}; // dummy GLM variables
        double meanSubst = 0.01;
        double epsilon = 0.01;

        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta", new RealParameter(betasArray),
                "subst", Double.toString(meanSubst),
                "epsilon", Double.toString(epsilon)
        );
        tipModel.initAndValidate();

        // nucleotides
        // types <- c(1, 1, 1, 3, 4, 1, 2, 1, 2, 3, 1, 3, 3, 4, 1)
        // names(types)
        // [1] "sp1"  "sp2"  "sp5"  "sp6"  "sp7"  "sp8"  "sp9"  "sp10" "sp11" "sp12" "sp13" "sp14"
        //[13] "sp15" "sp16" "sp17"
        Alignment alignment = getAlignmentLarge();
        //int numLeaves = 2;
        //String[] sequences = {"A", "C"};
        //Alignment alignment = getAlignment(numLeaves, sequences);
        TaxonSet taxonSet = new TaxonSet(alignment);
        int numTraits = 1;
        // traits
        TraitSet trait0 = new TraitSet();
        String trait0Values = "t0=1.0, t1=10.0";
        // taxa names
        // "sp1"  "sp2"  "sp5"  "sp6"  "sp7"  "sp8"  "sp9"  "sp10" "sp11" "sp12" "sp13" "sp14"
        //[13] "sp15" "sp16" "sp17"
        // trait values
        // sp1=0.010707127,
        // sp2=0.009700834,
        // sp5=0.010450199,
        // sp6=0.009184059,
        // sp7=0.010065188,
        // sp8=0.010411688,
        // sp9=0.010025265,
        // sp10=0.010276946,
        // sp11=0.009451115,
        // sp12=0.010506482,
        // sp13=0.009989017,
        // sp14=0.009375481,
        // sp15=0.008813563,
        // sp16=0.009812874,
        // sp17=0.010246622
        trait0.initByName(
                "traitname", "trait0",
                "taxa", new TaxonSet(alignment),
                "value", trait0Values);
        List<TraitSet> traitsList = new ArrayList<>(numTraits);
        traitsList.add(trait0);

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
