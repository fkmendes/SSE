package mosse;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Kylie Chen
 */

public class MosseTreeLikelihoodTest {

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
                "numRateBins", Integer.toString(numBins)
                );
        likelihood.calculateLogP();
    }

    public void testMosseExponentiatedMatrix() {

    }


    public void testMosseLikelihoodRootNode() {

    }

    public void testMosseLikelihoodOnTree() {

    }
}
