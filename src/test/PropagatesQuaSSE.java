package test;

import org.apache.commons.lang3.ArrayUtils;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

import static SSE.SSEUtils.propagateEandDinTQuaLike;

public class PropagatesQuaSSE {

    final static double EPSILON = 1e-6;

    int nUsefulTraitBins, nDimensions;
    double dt;
    double[] birthRate, deathRate;
    double[] esDs, scratch;

    @Test
    public void testPropagateTimeOneChQuaSSETest1() {

        esDs = new double[] { 0.0, 1.0 };
        scratch = new double[] { 0.0, 0.0 };
        birthRate = new double[] { 1.0 };
        deathRate = new double[] { 0.5 };
        dt = 0.01;
        nUsefulTraitBins = 1;
        nDimensions = 2; // 1 for E, 1 for D

        // propagating in place, result left in esDs
        propagateEandDinTQuaLike(esDs, scratch, birthRate, deathRate, dt, nUsefulTraitBins, nDimensions);

        // System.out.println(Arrays.toString(esDs));

        Double[] esDsArrayList = ArrayUtils.toObject(esDs);
        Double[] expected = new Double[] { 0.004962, 0.985161};
        // List<Double> list = Arrays.asList(esDsArrayList);
        assertArrayEquals(esDsArrayList, expected);
        // Assert.assertEquals(2.190298, lnLk1, EPSILON);
    }
}
