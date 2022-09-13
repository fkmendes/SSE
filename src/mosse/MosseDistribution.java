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

    static {
        System.loadLibrary("test");
    }

    private native long makeMosseFFT(int nx, double dx, int[] array_nd, int flags);

    private native void mosseFinalize(long obj_ptr);

    private native double[] doIntegrateMosse(long obj_ptr, double[] vars, double[] lambda, double[] mu,
            double drift, double diffusion, double[] Q, int nt, double dt, int pad_left, int pad_right);

    private double[][] eProbs;
    private double[][] dProbs;

    public void initAndValidate() {

    }

    public double getLogP(double[] array, int nx, int ncol, double dx, double[][] ans) {
//        double[][] ans = new double[nx][ncol];
        double logP = logCompensation(nx, ncol, dx, array, ans);
        return logP; // return log compensated result
    }

    public double[] doIntegration(int nx, double dx, int[] nd, int flags,
                              double[] vars, double[] lambda, double[] mu,
                              double drift, double diffusion,
                              double[] Q, int nt, double dt_max,
                              int pad_left, int pad_right) {
        // make mosse fft object pointer
        long ptr = makeMosseFFT(nx, dx, nd, flags);
        // integrate using C propagate x and propagate t
        double[] result = doIntegrateMosse(ptr, vars, lambda, mu, drift, diffusion, Q, nt, dt_max, pad_left, pad_right);
        mosseFinalize(ptr); // destroy obj pointer
        return result; // return non logged results
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
