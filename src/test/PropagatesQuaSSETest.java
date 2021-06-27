package test;

import SSE.SSEUtils;
import org.jtransforms.fft.DoubleFFT_1D;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.*;
import static SSE.SSEUtils.propagateEandDinTQuaLike;

public class PropagatesQuaSSETest {

    final static double EPSILON = 1e-6;

    int nUsefulTraitBins, nXbins, nLeftFlankBins, nRightFlankBins, nDimensions;
    double dt = 0.01;
    double dx = 0.001;
    double drift, diffusion;
    double[] birthRate, deathRate;
    double[] esDs, scratch;

    @Test
    public void testPropagateTimeOneChQuaSSETest() {

        esDs = new double[] { 0.0, 1.0 };
        scratch = new double[] { 0.0, 0.0 };
        birthRate = new double[] { 1.0 };
        deathRate = new double[] { 0.5 };
        nUsefulTraitBins = 1;
        nDimensions = 2; // 1 for E, 1 for D

        // propagating in place, result left in esDs
        propagateEandDinTQuaLike(esDs, scratch, birthRate, deathRate, dt, nUsefulTraitBins, nDimensions);

        // System.out.println("esDs=" + Arrays.toString(esDs));

        double[] expected = new double[] { 0.004962, 0.985161 };
        assertArrayEquals(expected, esDs, EPSILON);
    }

    /*
     * Makes Normal kernel
     */
    @Test
    public void testMakeNormalKernInPlaceAndFFT() {

        drift = 0.0;
        diffusion = 0.001;
        nLeftFlankBins = nRightFlankBins = 5;
        nXbins = 100;

        double[] fY = new double[200];
        double[] fftY = new double[100];
        double[] expectedStartFy = new double[] { 0.13723362, 0.13054066, 0.11235738, 0.08750402, 0.06166304, 0.03931809, 0.0, 0.0 };
        double[] expectedEndFy = new double[] { 0.0, 0.0, 0.03931809, 0.06166304, 0.08750402, 0.11235738, 0.13054066 };

        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        double[] startFy = Arrays.copyOfRange(fY, 0, 8);
        double[] endFy = Arrays.copyOfRange(fY, 93, 100);

        // System.out.println("fY=" + Arrays.toString(fY));

        DoubleFFT_1D fft = new DoubleFFT_1D(100);
        fft.realForwardFull(fY);

        int j=0;
        for (int i=0; i<fY.length; i+=2) {
            fftY[j] = fY[i];
            j++;
        }

        System.out.println(Arrays.toString(fftY));

        assertArrayEquals(expectedStartFy, startFy, EPSILON);
        assertArrayEquals(expectedEndFy, endFy, EPSILON);
    }

    @Test
    public void testPropagateXOneChQuaSSETest() {

        nXbins = 100;
        nLeftFlankBins = 10;
        nRightFlankBins = 10;
        drift = 0.0;
        diffusion = 1.0;
    }
}
