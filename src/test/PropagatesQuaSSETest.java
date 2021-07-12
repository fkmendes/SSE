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
    double[] fY;
    double[][] esDs;
    double[][] scratch;

    @Test
    public void testPropagateTimeOneChQuaSSETest() {

        esDs = new double[2][96];
        scratch = new double[2][96];
        birthRate = new double[48];
        deathRate = new double[48];
        nDimensionsE = nDimensionsD = 1;

        double[] randNumbers = new double[] { 8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = 0.0001;
            esDs[1][i] = randNumbers[i];
            birthRate[i] = 1.0;
            deathRate[i] = 0.5;
        }
        nUsefulTraitBins = 48 - 4 - 4;

        // propagating in place, result left in esDs
        propagateEandDinTQuaSSE(esDs, scratch, birthRate, deathRate, dt, nUsefulTraitBins, nDimensionsD);

        // System.out.println(Arrays.toString(esDs[0]));
        // System.out.println(Arrays.toString(esDs[1]));

        double[] expectedEs = new double[] { 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285, 0.005061285 };
        double[] expectedDs = new double[] { 8.383508, 9.953882, 8.954895, 11.671936, 8.386246, 9.769374, 8.817404, 9.168019, 10.913191, 10.604198 };
        assertArrayEquals(expectedEs, Arrays.copyOfRange(esDs[0], 0, 10), EPSILON);
        assertArrayEquals(expectedDs, Arrays.copyOfRange(esDs[1], 0, 10), EPSILON);
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

        double[] fY = new double[nXbins];

        double[] fftFY = new double[nXbins * 2];
        double[] realFftFy = new double[nXbins * 2];

        double[] ifftFY = new double[nXbins * 2];
        double[] realIfftFy = new double[nXbins * 2];

        // prepare fY
        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        // System.out.println("fY=" + Arrays.toString(fY));

        // fft
        for (int i=0; i<fY.length; i++) {
            fftFY[i] = fY[i];
        }
        fftForKern = new DoubleFFT_1D(nXbins);
        fftForKern.realForwardFull(fftFY);
        everyOtherInPlace(fftFY, realFftFy, nXbins, 0, 0, 1.0);
        // System.out.println(Arrays.toString(realFftFy));

        // ifft
        for (int i=0; i<fftFY.length; i++) {
            ifftFY[i] = fftFY[i];
        }
        fftForKern.complexInverse(ifftFY, false);
        everyOtherInPlace(ifftFY, realIfftFy, nXbins, 0, 0, 1.0);
        // System.out.println(Arrays.toString(realIfftFy));

        double[] expectedStartFy = new double[] { 0.14894618, 0.14168199, 0.12194682, 0.09497228, 0.06692583, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double[] expectedEndFy = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06692583, 0.09497228, 0.12194682, 0.14168199 };
        double[] expectedFftFY = new double[] { 1.000000000, 0.956873917, 0.835109753, 0.655887828, 0.449367562, 0.248270676, 0.081152021, -0.033081907, -0.088189035, -0.090559850 };
        double[] expectedIfftFY = new double[] { 7.149417e+00, 6.800735e+00, 5.853447e+00, 4.558669e+00, 3.212440e+00, 5.957485e-16, 4.549779e-16, 1.366926e-15, 2.766926e-16, 3.584917e-17 };

        assertArrayEquals(expectedStartFy, Arrays.copyOfRange(fY, 0, 10), EPSILON);
        assertArrayEquals(expectedEndFy, Arrays.copyOfRange(fY, 38, nXbins), EPSILON);
        assertArrayEquals(expectedFftFY, Arrays.copyOfRange(realFftFy, 0, 10), EPSILON);
        assertArrayEquals(expectedIfftFY, Arrays.copyOfRange(realIfftFy, 0, 10), EPSILON);
    }

    /*
     * Makes Normal kernel, FFTs it, then convolves with E's and D's.
     */
    @Test
    public void testConvolve() {

        nDimensionsE = nDimensionsD = 1;
        nLeftFlankBins = nRightFlankBins = 4;
        nXbins = 48;
        drift = 0.0;
        diffusion = 0.001;

        // number of bins must be multiple of 4
        fftForKern = new DoubleFFT_1D(nXbins);
        fftForEandD = new DoubleFFT_1D(nXbins);

        fY = new double[nXbins * 2];
        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        // fftForKern.realForwardFull(fY);

        // System.out.println(Arrays.toString(fY));

        // preparing input
        esDs = new double[2][nXbins * 2];
        scratch = new double[2][nXbins * 2];
        double[] randNumbers = new double[] { 8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = scratch[0][i] = 0.0001;
            esDs[1][i] = scratch[1][i] = randNumbers[i];
        }

        SSEUtils.convolveInPlace(scratch, fY, 1, 1, fftForEandD);

        everyOtherInPlace(scratch[0], esDs[0], nXbins, 0, 0, 1.0/nXbins); // grabbing real part and scaling by 1/nXbins
        everyOtherInPlace(scratch[1], esDs[1], nXbins, 0, 0, 1.0/nXbins);

        // System.out.println(Arrays.toString(esDs[0]));
        // System.out.println(Arrays.toString(esDs[1]));

        double[] expectedEs = new double[] { 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 };
        double[] expectedDs = new double[] { 9.73417583936596, 9.63964952722219, 9.58033651944056, 9.61334237373252, 9.70572479120779, 9.84274608587003, 9.93195677948007, 10.0416517201315, 10.14467157742, 10.4371836973515 };
        assertArrayEquals(expectedEs, Arrays.copyOfRange(esDs[0], 0, 10), EPSILON);
        assertArrayEquals(expectedDs, Arrays.copyOfRange(esDs[1], 0, 10), EPSILON);
    }

    /*
     * Calls convolve function, which makes Normal kernel, FFTs it, convolves,
     * and iFFTs it. Then reorganize all elements of esDs by bookkeeping all
     * flanking bins, etc.
     */
    @Test
    public void testPropagateChOneChQuaSSETest() {

        nDimensionsE = nDimensionsD = 1;
        nLeftFlankBins = nRightFlankBins = 4;
        nXbins = 48;
        drift = 0.0;
        diffusion = 0.001;

        esDs = new double[2][nXbins * 2];
        scratch = new double[2][nXbins * 2];

        // number of bins must be multiple of 4
        fftForEandD = new DoubleFFT_1D(nXbins);

        fY = new double[nXbins * 2];
        SSEUtils.makeNormalKernelInPlace(fY, drift, diffusion, nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already

        esDs = new double[2][nXbins * 2];
        scratch = new double[2][nXbins * 2];
        double[] randNumbers = new double[] { 8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = scratch[0][i] = 0.0001;
            esDs[1][i] = scratch[1][i] = randNumbers[i];
        }

        propagateEandDinXQuaLike(esDs, scratch, fY, nXbins, nLeftFlankBins, nRightFlankBins, nDimensionsE, nDimensionsD, fftForEandD);

        double[] expectedEsFirstHalf = new double[] { 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedEsSecondHalf = new double[nXbins];
        double[] expectedDsFirstHalf = new double[] { 8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 9.70572479120779, 9.84274608587003, 9.93195677948007, 10.0416517201315, 10.14467157742, 10.4371836973515, 10.4801274082435, 10.3773781820321, 10.3654783428447, 10.1701178809958, 10.0016160413878, 9.81012635281443, 9.67902839214041, 9.55810533268215, 9.52689648206184, 9.55977301973422, 9.76757793553243, 9.79454433623556, 9.92251398623884, 10.0126014859394, 10.0172029657907, 10.1066796832284, 9.90685464014484, 9.67801788862755, 9.66356740459434, 9.48033275476296, 9.40446378096838, 9.51689441708817, 9.67920291038787, 9.91972188887997, 10.0248943989667, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedDsSecondHalf = new double[nXbins];

        for (int i=0; i<nXbins; i++) {
            expectedEsSecondHalf[i] = expectedDsSecondHalf[i]  = 0.0;
        }

        assertArrayEquals(expectedEsFirstHalf, Arrays.copyOfRange(esDs[0], 0, nXbins), EPSILON);
        assertArrayEquals(expectedEsSecondHalf, Arrays.copyOfRange(esDs[0], nXbins, 2 * nXbins), 0.0);
        assertArrayEquals(expectedDsFirstHalf, Arrays.copyOfRange(esDs[1], 0, nXbins), EPSILON);
        assertArrayEquals(expectedDsSecondHalf, Arrays.copyOfRange(esDs[1], nXbins, 2 * nXbins), 0.0);
    }
}
