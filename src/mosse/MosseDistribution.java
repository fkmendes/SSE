package mosse;

import beast.core.Description;
import beast.core.State;
import beast.evolution.tree.TreeDistribution;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.assertEquals;

@Description("Mosse tree model")
public class MosseDistribution extends TreeDistribution {

    final static double DELTA = 1E-8;

    static {
        System.loadLibrary("test");
    }

    private native long makeMosseFFT(int nx, double dx, int[] array_nd, int flags);

    private native void mosseFinalize(long obj_ptr);

    private native double[] doIntegrateMosse(long obj_ptr, double[] vars, double[] lambda, double[] mu,
            double drift, double diffusion, double[] Q, int nt, double dt, int pad_left, int pad_right);

    public static void main(String[] args) throws IOException {
        MosseDistribution d = new MosseDistribution();

        int[] nd = {5}; // number of dimensions for each fft plan
        int nx = 1024; // number of bins
        double dx_input = 0.0001;
        int r = 4;
        double dx = dx_input * r;

        // reading in test case input
        double[] vars = d.readArray("testcase/integrate_vars.txt");
        double[] lambda = d.readArray("testcase/integrate_lambda.txt");
        double[] mu = d.readArray("testcase/integrate_mu.txt");
        double[] Q = d.readArray("testcase/integrate_Q.txt");

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
        System.out.println("nt: " + nt);

        int pad_left = 40;
        int pad_right = pad_left;

        double[][] ans = new double[nx][vars.length / nx];
        double logP = d.doIntegration(
                nx, dx, nd, 0,
                vars, lambda, mu, drift, diffusion, Q, nt, dt_max,
                pad_left, pad_right, ans);

        // print results
        System.out.println("logP: " + logP);

        double expected_ans_r = -0.2653163;
        assertEquals(expected_ans_r, logP, DELTA);
    }

    private double[][] eProbs;
    private double[][] dProbs;

    public void initAndValidate() {

    }

    public double doIntegration(int nx, double dx, int[] nd, int flags,
                              double[] vars, double[] lambda, double[] mu,
                              double drift, double diffusion,
                              double[] Q, int nt, double dt_max,
                              int pad_left, int pad_right, double[][] ans) {
        // make mosse fft object pointer
        long ptr = makeMosseFFT(nx, dx, nd, flags);
        // integrate using C propagate x and propagate t
        double[] result = doIntegrateMosse(ptr, vars, lambda, mu, drift, diffusion, Q, nt, dt_max, pad_left, pad_right);
        // log compensation
        int ncol = vars.length / nx;
        double logP = logCompensation(nx, ncol, dx, result, ans);
        mosseFinalize(ptr); // destroy obj pointer
        return logP;
    }

    public double logCompensation(int nrow, int ncol, double dx, double[] result, double[][] ans) {
        double logP = 0.0;
        int count = 0;
        for (int j = 0; j < ncol; j++) {
            for (int i = 0; i < nrow; i++) {
                ans[i][j] = result[count];
                count++;
            }
        }
        if (ncol > 1 ) {
            double sum = 0.0;
            // sum except first col
            for (int i = 0; i < nrow; i++) {
                for (int j = 1; j < ncol; j++) {
                    sum += ans[i][j] * dx;
                }
            }
            double q = sum;
            // update ans except first col
            for (int i = 0; i < nrow; i++) {
                for (int j = 1; j < ncol; j++) {
                    ans[i][j] = ans[i][j] / q;
                }
            }
            logP = Math.log(q);
        } else {
            logP = 0.0;
        }
        return logP;
    }

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

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new NotImplementedException();
    }

}
