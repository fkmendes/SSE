package src.test;

import org.junit.Assert;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import SSE.StateDependentSpeciationExtinctionProcess;

public class TestHelper {
    final static double EPSILON = 1e-5;
    public static int numNodesToNumTips(int numNodes) {
        return (int) (1.0 * (numNodes + 1) / 2);
    }

    public static double[] trimTips(double[] posterior) {
        // Remove entries of the tips, deep copy
        int numSpecies = numNodesToNumTips(posterior.length);
        double[] newPosterior = new double[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }


    public static int[] trimTipsInt(int[] posterior) {
        // Remove entries of the tips
        int numSpecies = numNodesToNumTips(posterior.length);
        int[] newPosterior = new int[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }

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

    public static HashMap<String, Double> getDivMap(String[] divLbls, String[] divLks) {
        HashMap<String, Double> divData = new HashMap<String, Double>();
        for (int i = 0; i < divLbls.length; i++) {
            divData.put(divLbls[i], Double.valueOf(divLks[i]));
        }
        return divData;
    }

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
        double accuracy = 1.0 * numCorrect / numNodes;
        return accuracy;
    }

    public static void compareDivPosterior(HashMap<String, Double> divMap, String[] idxLabelMapper, double[] post) {
        String lbl;
        double postBeast, postDiv;
        for (int i = 0; i < idxLabelMapper.length; i++) {
            lbl = idxLabelMapper[i];
            postBeast = post[i];
            postDiv = divMap.get(lbl);
//			System.out.println("" + postBeast + ", " + postDiv);
            Assert.assertEquals(postDiv, postBeast, 1e-1);
        }
    }

    public static void compareArr(double[] arr1, double[] arr2) {
        for (int i = 0; i < arr1.length; i++) {
            Assert.assertEquals(arr1[i], arr2[i], EPSILON);
        }
    }

    public static void prepareAndWriteToCSV(double[] posteriorWTips, String expName, StateDependentSpeciationExtinctionProcess sdsep) throws Exception {
        int numSpecies = (int) ((posteriorWTips.length + 1) / 2);
        double[] posterior = TestHelper.trimTips(posteriorWTips);
        String dir = "/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/scm/beast";
        String fileName = expName + ".csv";
        TestHelper.writeToCSV(dir, fileName, posterior, sdsep);
    }

    public static int[] parsimony(double[] arr) {
        // for each element, check the state given the max element entry
        // Picks the state in which the node is most likely to be in
        // This only works for 2 states TODO Extend to CLaSSE
        int[] ret = new int[arr.length];
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > 0.5) {
                ret[i] = 2;
            } else {
                ret[i] = 1;
            }
        }
        return ret;
    }
}
