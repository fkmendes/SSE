package mosse;

import beast.core.Description;
import beast.core.State;
import beast.evolution.tree.TreeDistribution;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.List;
import java.util.Random;

@Description("Mosse tree prior model")
public class MosseDistribution extends TreeDistribution {

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
