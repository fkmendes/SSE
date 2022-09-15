package mosse;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class MosseDistributionTest {

    private static double DELTA = 1e-7;

    private double[] readArray(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = reader.readLine();
        List<Double> list = new ArrayList<Double>();
        while (line != null) {
            String[] elements = line.trim().split(", ");
            for (String i: elements) {
                if (i != "") {
                    list.add(Double.parseDouble(i.trim()));
                }
            }
            line = reader.readLine();
        }
        reader.close();
        double[] arr = new double[list.size()];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }

    @Test
    public void testIntegrateBranchOne() throws IOException {
        MosseDistribution d = new MosseDistribution();

        int[] nd = {5}; // number of dimensions for each fft plan
        int nx = 1024; // number of bins
        double dx_input = 0.0001;
        int r = 4;
        double dx = dx_input * r;

        // reading in test case input
        double[] vars = readArray("testcase/integrate_vars.txt");
        double[] lambda = readArray("testcase/integrate_lambda.txt");
        double[] mu = readArray("testcase/integrate_mu.txt");
        double[] Q = readArray("testcase/integrate_Q.txt");

        assertEquals(5120, vars.length);
        assertEquals(943, lambda.length);
        assertEquals(943, mu.length);
        assertEquals(15088, Q.length);

        // diffusion parameters
        double drift = 0.0;
        double diffusion = 0.001;

        // time step parameters
        double dt_max = 0.01;
        double len = 2; // total branch time
        int nt = (int) Math.ceil(len / dt_max); // number of time steps (200 for low res)

        int pad_left = 40;
        int pad_right = pad_left;

        int ncol = vars.length / nx;
        double[][] ans = new double[nx][ncol];
        double[] result = d.doIntegration(
                nx, dx, nd, 0,
                vars, lambda, mu, drift, diffusion, Q, nt, dt_max,
                pad_left, pad_right);
        double logP = d.calculateLogP(result, nx, ncol, dx, ans);

        double expected_ans_r = -0.2653163;
        assertEquals(expected_ans_r, logP, DELTA);
    }
}
