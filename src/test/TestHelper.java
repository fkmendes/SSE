package src.test;

import org.junit.Assert;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;

public class TestHelper {
    final static double EPSILON = 1e-5;

    public static double[] trimTips(double[] posterior, int numSpecies) {
        // Remove entries of the tips
        double[] newPosterior = new double[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }

    public static int[] trimTipsInt(int[] posterior, int numSpecies) {
        // Remove entries of the tips
        int[] newPosterior = new int[posterior.length - numSpecies];
        System.arraycopy(posterior, numSpecies, newPosterior, 0, newPosterior.length);

        return newPosterior;
    }

    public static void writeToCSV(String name, double[] arr, SSE.StateDependentSpeciationExtinctionProcess sdsep) throws Exception {
        BufferedWriter br = new BufferedWriter(new FileWriter(name));
        StringBuilder sb = new StringBuilder();

        String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
        for (String id: indexNameMapper) {
            sb.append(id);
            sb.append(",");
        }

        sb.append("\n");

        for (int i = 0; i < arr.length; i++) {
            double element = arr[i];
            sb.append(Double.toString(element));
            sb.append(",");
        }

        br.write(sb.toString());
        br.close();
    }

    public static void compareDiv(String[] divLbls, String[] divLks, String[] idxLabelMapper, double[] post) {
        HashMap<String, Double> divData = new HashMap<String, Double>();
        for (int i = 0; i < divLbls.length; i++) {
            divData.put(divLbls[i], Double.valueOf(divLks[i]));
        }

        String lbl;
        double postBeast, postDiv;
        for (int i = 0; i < idxLabelMapper.length; i++) {
            lbl = idxLabelMapper[i];
            postBeast = post[i];
            postDiv = divData.get(lbl);
//			System.out.println("" + postBeast + ", " + postDiv);
            Assert.assertEquals(postDiv, postBeast, 1e-1);
        }
    }

    public static void compareArr(double[] arr1, double[] arr2) {
        for (int i = 0; i < arr1.length; i++) {
            Assert.assertEquals(arr1[i], arr2[i], EPSILON);
        }
    }
}
