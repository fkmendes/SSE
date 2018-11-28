package test;

import org.junit.Assert;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Arrays;

import SSE.StateDependentSpeciationExtinctionProcess;

public class TestHelper {
    final static double EPSILON = 1e-5;

    public static int argmax(double[] arr) {
        int index = 0;
        double max = arr[0];
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > max) {
                max = arr[i];
                index = i;
            }
        }
        return index;
    }

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

    /**
     *
     * @param divLabelMapper: See above divMap
     * @param beastLabelMapper: See above idxLabelMapper
     * @return mapping from diversitree indices to beast indices
     */
    public static Integer[] getIndexMapping(String[] divLabelMapper, String[] beastLabelMapper) {
        Integer[] divBeastMapper = new Integer[divLabelMapper.length];

        for (int i = 0; i < divLabelMapper.length; i++) {
            String divLabel = divLabelMapper[i];
            int beastLabelIdx = Arrays.asList(beastLabelMapper).indexOf(divLabel);
            divBeastMapper[i] = beastLabelIdx;
        }
        return divBeastMapper;
    }


    // Verifies that the arrays have the same values for the first arr1.length elements
    public static void compareArr(double[] arr1, double[] arr2) {
        for (int i = 0; i < arr1.length; i++) {
            Assert.assertEquals(arr1[i], arr2[i], 1e-1);
        }
    }

    public static void compareNestedArr(double[][] arr1, double[][] arr2) {
        for (int i = 0; i < arr1.length; i++) {
            compareArr(arr1[i], arr2[i]);
        }
    }

    /**
     * prepare data by trimming tips and then write to csv
     */
    public static void prepareAndWriteToCSV(double[] posteriorWTips, String expName, StateDependentSpeciationExtinctionProcess sdsep) throws Exception {
        int numSpecies = (int) ((posteriorWTips.length + 1) / 2);
        double[] posterior = TestHelper.trimTips(posteriorWTips);
        String dir = "/Users/jeff/Documents/Research/Phylogenetics/SSE/validation";
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

    /**
     * TODO Should I make it a deep copy
     * @param posteriors posterior probabilities for all tree nodes for all states
     *                   size numTotalNodes x numStates
     * @return posterior probabilities for all internal nodes for all states
     */
    public static double[][] trimTipsCLaSSE(double[][] posteriors) {
        int numTotalNodes = posteriors.length;
        int numTips = numNodesToNumTips(numTotalNodes);
        int numInternalNodes = numTotalNodes - numTips;
        int numStates = posteriors[0].length;

        double[][] ret = new double[numInternalNodes][numStates];

        for (int i = 0; i < numInternalNodes; i++) {
            ret[i] = posteriors[numTips + i];
        }

        return ret;
    }

    /**
     * TODO Should I put this in the CLaSSECredibleSetTest instead of here? might do this for BiSSE too, so put here
     * Note: this method will destory probs. Maybe make a deep copy
     * @param probs probs[i][j] is posterior probability node i is in state j
     * @return a credible set for each node
     */
    public static Set<Integer>[] constructCredibleSets(double[][] probs, double credibleThreshold) {
        int numIntNodes = probs.length;
        Set<Integer>[] credibleSets = new HashSet[numIntNodes];

        double credibility;
        for (int i = 0; i < numIntNodes; i++) {
            credibility = 0.0;
            credibleSets[i] = new HashSet<>();
            while (credibility < credibleThreshold) {
                int likelyState = argmax(probs[i]) + 1;
                credibleSets[i].add(likelyState);
                credibility += probs[i][likelyState - 1];
                probs[i][likelyState - 1] = -1;  // so that we won't pick it again
            }
        }
        return credibleSets;
    }

    /**
     * TODO Should I put this in the CLaSSECredibleSetTest instead of here?
     * @param divMap Map from node label/ID/name to the true state that the node is in.
     *               in diversitree, states are indexed by 0
     * @param credibleSets a set for every node in which
     *                     for every node, the sum of the posteriors of those states adds up to be >= credible threshold
     * @param idxLabelMapper idxLabelMapper[nodeIdx] provides the node label/ID/name
     * @return accuracy measure. How often is the ground truth in our credible set
     */
    public static double computeProportionTruthInCredibleSet(HashMap<String, Double> divMap, Set<Integer>[] credibleSets, String[] idxLabelMapper) {
        int numIntNodes = idxLabelMapper.length;

        int numInside = 0;
        for (int i = 0; i < numIntNodes; i++) {
            String lbl = idxLabelMapper[i];
            Set<Integer> credibleSet = credibleSets[i];
            if (credibleSet.isEmpty()) {
                Assert.fail("The credible set is empty!");
            }
            int trueState = divMap.get(lbl).intValue();
            if (credibleSet.contains(trueState)) {
                numInside++;
            }
        }

        return 1.0 * numInside / numIntNodes;
    }

    public static String readFirstLine(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        try {
            line = br.readLine();
        } finally {
            br.close();
        }
        return line;
    }

    public static ArrayList<String> readFile(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        ArrayList<String> strArr = new ArrayList<>();
        try {
            String line = br.readLine();

            while (line != null) {
                strArr.add(line);
                line = br.readLine();
            }
        } finally {
            br.close();
        }
        return strArr;
    }
}
