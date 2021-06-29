package test;

import SSE.SSEUtils;
import org.jtransforms.fft.DoubleFFT_1D;
import org.junit.Test;

import java.util.Arrays;

import static SSE.SSEUtils.*;
import static org.junit.Assert.*;

public class PropagatesQuaSSETest {

    final static double EPSILON = 1e-5;

    DoubleFFT_1D fftForKern, fftForEandD;
    int nUsefulTraitBins, nXbins, nLeftFlankBins, nRightFlankBins, nDimensionsE, nDimensionsD;
    double dt = 0.01;
    double dx = 0.001;
    double drift, diffusion;
    double[] birthRate, deathRate;
    double[][] esDs;
    double[] scratch;

    @Test
    public void testPropagateTimeOneChQuaSSETest() {

        esDs = new double[][] { {0.0}, {1.0} };
        scratch = new double[] { 0.0, 0.0 };
        birthRate = new double[] { 1.0 };
        deathRate = new double[] { 0.5 };
        nUsefulTraitBins = nDimensionsD = 1;

        // propagating in place, result left in esDs
        propagateEandDinTQuaLike(esDs, scratch, birthRate, deathRate, dt, nUsefulTraitBins, nDimensionsE, nDimensionsD);

        // System.out.println("esDs=" + Arrays.toString(esDs));

        double[][] expected = new double[][] { {0.004962}, {0.985161} };
        assertArrayEquals(expected[0], esDs[0], EPSILON);
    }

    /*
     * Makes Normal kernel, then does FFT on kernel, and inverse FFT on the FFT.
     */
    @Test
    public void testMakeNormalKernInPlaceAndFftAndIfft() {

        drift = 0.0;
        diffusion = 0.001;
        nLeftFlankBins = nRightFlankBins = 4;
        nXbins = 48;

        double[] fY = new double[48];

        double[] fftFY = new double[96];
        double[] realFftFy = new double[48];

        double[] ifftFY = new double[96];
        double[] realIfftFy = new double[48];

        // prepare fY
        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        // System.out.println("fY=" + Arrays.toString(fY));

        // fft
        for (int i=0; i<fY.length; i++) {
            fftFY[i] = fY[i];
        }
        fftForKern = new DoubleFFT_1D(48);
        fftForKern.realForwardFull(fftFY);
        everyOtherInPlace(fftFY, realFftFy, true, 1.0);
        // System.out.println(Arrays.toString(realFftFy));

        // ifft
        for (int i=0; i<fftFY.length; i++) {
            ifftFY[i] = fftFY[i];
        }
        fftForKern.complexInverse(ifftFY, false);
        everyOtherInPlace(ifftFY, realIfftFy, true, 1.0);
        // System.out.println(Arrays.toString(realIfftFy));

        double[] expectedStartFy = new double[] { 0.14894618, 0.14168199, 0.12194682, 0.09497228, 0.06692583, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double[] expectedEndFy = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06692583, 0.09497228, 0.12194682, 0.14168199 };
        double[] expectedFftFY = new double[] { 1.000000000, 0.956873917, 0.835109753, 0.655887828, 0.449367562, 0.248270676, 0.081152021, -0.033081907, -0.088189035, -0.090559850 };
        double[] expectedIfftFY = new double[] { 7.149417e+00, 6.800735e+00, 5.853447e+00, 4.558669e+00, 3.212440e+00, 5.957485e-16, 4.549779e-16, 1.366926e-15, 2.766926e-16, 3.584917e-17 };

        assertArrayEquals(expectedStartFy, Arrays.copyOfRange(fY, 0, 10), EPSILON);
        assertArrayEquals(expectedEndFy, Arrays.copyOfRange(fY, 38, 48), EPSILON);
        assertArrayEquals(expectedFftFY, Arrays.copyOfRange(realFftFy, 0, 10), EPSILON);
        assertArrayEquals(expectedIfftFY, Arrays.copyOfRange(realIfftFy, 0, 10), EPSILON);
    }

    @Test
    public void testConvolve() {

        // number of bins must be multiple of 4
        fftForKern = new DoubleFFT_1D(48);
        fftForEandD = new DoubleFFT_1D(48);
        nDimensionsE = nDimensionsD = 1;

        drift = 0.0;
        diffusion = 0.001;
        nLeftFlankBins = nRightFlankBins = 4;
        nXbins = 48;

        double[] fY = new double[96];
        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        fftForKern.realForwardFull(fY);

        // System.out.println(Arrays.toString(fY));

        esDs = new double[2][96];
        double[] randNumbers = new double[] { 8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = 0.0001;
            esDs[1][i] = randNumbers[i];
        }

        // System.out.println(Arrays.toString(esDs[0]));
        // System.out.println(Arrays.toString(esDs[1]));
        double[] realEs = new double[48];
        double[] realDs = new double[48];
        // FFT-ing E's and then D's
        for (int i=0; i<2; i++) {
            fftForEandD.realForwardFull(esDs[i]);
        }

        // just for testing, not used in ifft
        everyOtherInPlace(esDs[0], realEs, true, 1.0);
        everyOtherInPlace(esDs[1], realDs, true, 1.0);

        System.out.println(Arrays.toString(realEs));
        System.out.println(Arrays.toString(realDs));

        for (int i=0, j=0; i<fY.length; i+=2, j++) {
            esDs[0][i] *= fY[i];
            esDs[0][i+1] *= fY[i];
            esDs[1][i] *= fY[i];
            esDs[1][i+1] *= fY[i];
            realEs[j] *= fY[i]; // just for testing, not used in ifft
            realDs[j] *= fY[i]; // just for testing, not used in ifft
        }

        System.out.println(Arrays.toString(realEs));
        System.out.println(Arrays.toString(realDs));

        // inverse-FFT-ing E's and then D's
        for (int i=0; i<2; i++) {
            fftForEandD.complexInverse(esDs[i], false);
        }

        everyOtherInPlace(esDs[0], realEs, true, 1.0/48.0);
        everyOtherInPlace(esDs[1], realDs, true, 1.0/48.0);

        System.out.println(Arrays.toString(realEs));
        System.out.println(Arrays.toString(realDs));
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
