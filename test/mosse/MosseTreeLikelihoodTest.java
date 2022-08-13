package mosse;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.treeparser.NewickParser;
import com.sun.xml.internal.bind.v2.model.core.EnumLeafInfo;
import org.junit.Test;
import test.beast.BEASTTestCase;
import test.beast.evolution.tree.TraitSetTest;

import java.util.ArrayList;
import java.util.List;

public class MosseTreeLikelihoodTest {

    private Alignment getAlignment(int numLeaves, String[] sequences) {
        List<Sequence> seqList = new ArrayList<Sequence>();

        for (int i = 0; i < numLeaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, sequences[i]));
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");

        return alignment;
    }

    @Test
    public void testMosseLikelihoodInit() throws Exception {
        int numLeaves = 2;
        String[] sequences = {"A", "C"};
        Alignment alignment = getAlignment(numLeaves, sequences); // BEASTTestCase.getAlignment();
        String newick = "(t0: 0.5, t1: 0.5);";
        Tree tree = new Tree(newick);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName(
                "mutationRate", "1.0",
                "gammaCategoryCount", 1,
                "substModel", JC);

        double beta0 = 0.1;
        double beta1 = 0.2;
        double epsilon = 0.01;
        MosseTipLikelihood tipModel = new MosseTipLikelihood();
        tipModel.initByName(
                "beta0", Double.toString(beta0),
                "beta1", Double.toString(beta1),
                "epsilon", Double.toString(epsilon)
        );
        tipModel.initAndValidate();


        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet traitData = new TraitSet();
        traitData.initByName(
                "traitname", "traits",
                "taxa", taxonSet,
                "value", "t0=1.0, t1=10.0");

        double startSubsRate = 1E-10;
        double endSubsRate = 1E-8;
        int numBins = 100;

        MosseTreeLikelihood likelihood = new MosseTreeLikelihood();
        likelihood.initByName(
                "data", alignment,
                "tree", tree,
                "siteModel", siteModel,
                "tipModel", tipModel,
                "traits", traitData,
                "startSubsRate", Double.toString(startSubsRate),
                "endSubsRate", Double.toString(endSubsRate),
                "numRateBins", Integer.toString(numBins)
                );
        likelihood.calculateLogP();
    }
}
