package src.test;

import org.junit.Assert;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import SSE.StateDependentSpeciationExtinctionProcess;

public class TestHelper {
    final static double EPSILON = 1e-5;

    // Get number of tips from number of total nodes (internal nodes and tips)
    public static int numNodesToNumTips(int numNodes) {
        return (int) (1.0 * (numNodes + 1) / 2);
    }

    // Get just the internal nodes indices (deep copy)
    public static double[] trimTips(double[] posterior) {
        int numSpecies = numNodesToNumTips(posterior.length);
        double[] newPosterior = new double[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }

    // Get just the internal nodes indices (deep copy)
    public static int[] trimTipsInt(int[] posterior) {
        int numSpecies = numNodesToNumTips(posterior.length);
        int[] newPosterior = new int[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }

    /**
     * First, store all of the names/labels of the nodes in order. Then, store all of the values of arr.
     * ex.
     *      arr = [1,2,3,4]
     *
     *      foo.csv
     *      sp1, sp2, sp3, sp5
     *      1, 2, 3, 4
     *
     */
    public static void writeToCSV(String dir, String name, double[] arr, SSE.StateDependentSpeciationExtinctionProcess sdsep) throws Exception {
        BufferedWriter br = new BufferedWriter(new FileWriter(new File(dir, name)));
        StringBuilder sb = new StringBuilder();

        String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
        for (int i = 0; i < indexNameMapper.length; i++) {
            String id = indexNameMapper[i];
            sb.append(id);
            if (i != indexNameMapper.length - 1) {
                sb.append(",");
            }
        }

        sb.append("\n");

        for (int i = 0; i < arr.length; i++) {
            double element = arr[i];
            sb.append(Double.toString(element));
            if (i != arr.length - 1) {
                sb.append(",");
            }
        }

        br.write(sb.toString());
        br.close();
    }

    /**
     * helps us work with the diversitree data in the format we receive it
     * @param divLbls: an array of the labels that diversitree gives the nodes e.g. ["sp2", "sp5", "sp7"]
     * @param divLks: an array of the likelihoods or true states obtained through diversitree for the nodes
     *              e.g. [0.4, 0.6, 0.8]
     *              e.g. [1, 2, 2]
     * @return a mapping from the labels to the likelihood or state
     */
    public static HashMap<String, Double> getDivMap(String[] divLbls, String[] divLks) {
        HashMap<String, Double> divData = new HashMap<String, Double>();
        for (int i = 0; i < divLbls.length; i++) {
            divData.put(divLbls[i], Double.valueOf(divLks[i]));
        }
        return divData;
    }

    /**
     * @param divMap: diversitree mapping from labels to true states to {0, 1}
     * @param idxLabelMapper: our mapping from node number to node label/ID/name
     *                      array where idxLabelMapper[i] is the label of node i (internal nodes indexed from 0)
     *                      e.g. ["sp4", "sp7"] means node 0 is sp4, node 1 is sp7
     * @param parsimony: most likely states for all nodes separately. parsimony[i] is the most likely state for node i
     *                 index by 1
     * @return accuracy measure. How good parsimony does compared to diversitree truth
     */
    public static double compareDivTruth(HashMap<String, Double> divMap, String[] idxLabelMapper, int[] parsimony) {
        int numNodes = idxLabelMapper.length;
        int numCorrect = 0;
        for (int i = 0; i < numNodes; i++) {
            String lbl = idxLabelMapper[i];
            int expectedState = parsimony[i] - 1;
            int trueState = divMap.get(lbl).intValue();
            if (expectedState == trueState) {
                numCorrect++;
            }
        }
        return 1.0 * numCorrect / numNodes;
    }

    /**
     *  Compares diversitree asr.marginal posterior with our method's posterior. Assert will end process if not same
     * @param divMap See above
     * @param idxLabelMapper See above
     * @param post discrete marginal posterior distribution for all nodes. post[i] is the posterior probability that
     *             node i is in state 1
     */
    public static void compareDivPosterior(HashMap<String, Double> divMap, String[] idxLabelMapper, double[] post) {
        /*
        idxLabelMapper - mapping from nodeIndex to node label, indexing for internal nodes starting at 0
         */
        String lbl;
        double postBeast, postDiv;
        if (idxLabelMapper.length != post.length) {
            post = trimTips(post);
        }
        for (int i = 0; i < idxLabelMapper.length; i++) {
            lbl = idxLabelMapper[i];
            postBeast = post[i];
            postDiv = divMap.get(lbl);
//			System.out.println("" + postBeast + ", " + postDiv);
            Assert.assertEquals(postDiv, postBeast, 1e-1);
        }
    }

    // Verifies that the arrays have the same values for the first arr1.length elements
    public static void compareArr(double[] arr1, double[] arr2) {
        for (int i = 0; i < arr1.length; i++) {
            Assert.assertEquals(arr1[i], arr2[i], EPSILON);
        }
    }

    /**
     * prepare data by trimming tips and then write to csv
     */
    public static void prepareAndWriteToCSV(double[] posteriorWTips, String expName, StateDependentSpeciationExtinctionProcess sdsep) throws Exception {
        int numSpecies = (int) ((posteriorWTips.length + 1) / 2);
        double[] posterior = TestHelper.trimTips(posteriorWTips);
        String dir = "/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/scm/beast";
        String fileName = expName + ".csv";
        TestHelper.writeToCSV(dir, fileName, posterior, sdsep);
    }

    /**
     * for each element/node, determines the state in which the node is most likely to be in
     * This only works for 2 states
     * TODO Extend to CLaSSE
     * @param arr arr[i] is the likelihood that node i is in state 1
     * @return parsimony array. ret[i] is most likely state that node i is in
     */
    public static int[] parsimony(double[] arr) {
        int[] ret = new int[arr.length];
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > 0.5) {  // since the value is likelihood of state 0
                ret[i] = 1;
            } else {
                ret[i] = 2;
            }
        }
        return ret;
    }
}
