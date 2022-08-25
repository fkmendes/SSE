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

    public static void main(String[] args) {

        MosseDistribution d = new MosseDistribution();
        int arr[] = {1,2,3};
        long ptr = d.makeMosseFFT(1024, 0.1, arr, 0);
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

    public void propagateX() {

    }

    public void propagateT() {

    }

    public void fft() {

    }

}
