package mosse;

import beast.core.Description;
import beast.core.State;
import beast.evolution.tree.TreeDistribution;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.List;
import java.util.Random;

@Description("Mosse tree prior model")
public class MosseDistribution extends TreeDistribution {

    static {
        System.loadLibrary("test");
    }

    private native long makeMosseFFT(int nx, double dx, int[] array_nd, int flags);

    private native void mosseFinalize(long obj_ptr);

    private native double[] getX(long obj_ptr);

    private native double[] doIntegrateMosse(long obj_ptr, double[] vars, double[] lambda, double[] mu,
            double drift, double diffusion, double[] Q, int nt, double dt, int pad_left, int pad_right);

    public static void main(String[] args) {
        MosseDistribution d = new MosseDistribution();
        int[] arr = {0,1,2,3};
        long ptr = d.makeMosseFFT(1024, 0.1, arr, 0);
        double[] result = d.getX(ptr);

//        double[] vars = {1.0, 1.0, 1.0};
//        double[] lambda = {1.0, 1.0, 1.0};
//        double[] mu = {1.0, 1.0, 1.0};
//        double drift = 0.0;
//        double diffusion = 0.02;
//        double[] Q = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//        int nt = 4;
//        double dt = 0.001;
//        int pad_left = 5;
//        int pad_right = 5;
//        double[] result = d.doIntegrateMosse(ptr, vars, lambda, mu, drift, diffusion, Q, nt, dt, pad_left, pad_right);

        for (double x: result) {
            System.out.print(x + " ");
        }

        d.mosseFinalize(ptr);
    }

    private double[][] eProbs;
    private double[][] dProbs;

    public void initAndValidate() {

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
