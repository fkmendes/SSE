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
    double dx = 0.01; // before I used 0.001
    double drift, diffusion;
    double[] birthRate, deathRate;
    double[] fY;
    double[][] esDs;
    double[][] scratch;

    /*
     * Propagates esDs in time. These are the analytical solutions
     * for the birth-death model.
     */
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
        propagateEandDinTQuaSSEInPlace(esDs, scratch, birthRate, deathRate, dt, nUsefulTraitBins, nDimensionsD);

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
        SSEUtils.makeNormalKernelInPlace(fY, dt * drift, Math.sqrt(dt * diffusion), nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        // System.out.println("fY = " + Arrays.toString(fY));

        // fft
        for (int i=0; i<fY.length; i++) {
            fftFY[i] = fY[i];
        }
        fftForKern = new DoubleFFT_1D(nXbins);
        fftForKern.realForwardFull(fftFY);
        everyOtherInPlace(fftFY, realFftFy, nXbins, 0, 0, 1.0);
        // System.out.println("fft fY = " + Arrays.toString(realFftFy));

        // ifft
        for (int i=0; i<fftFY.length; i++) {
            ifftFY[i] = fftFY[i];
        }
        fftForKern.complexInverse(ifftFY, false);
        everyOtherInPlace(ifftFY, realIfftFy, nXbins, 0, 0, 1.0);
        // System.out.println("ifft fY =" + Arrays.toString(realIfftFy));

        double[] expectedFy = new double[] { 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386 };
        double[] expectedFftFY = new double[] { 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113 };
        double[] expectedIfftFY = new double[] { 47.3617577773852, 0.319121013687545, 9.76198587343688e-08, 7.45931094670027e-17, -4.44089209850063e-16, 1.51614216742138e-16, 4.05992386804242e-16, -1.43599463769065e-16, -4.44089209850063e-16, -5.17654502592123e-16 };

        assertArrayEquals(expectedFy, fY, EPSILON);
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
        SSEUtils.makeNormalKernelInPlace(fY, (dt * drift), Math.sqrt(dt * diffusion), nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        fftForKern.realForwardFull(fY); // first FFT the Normal kernel

        // preparing input
        esDs = new double[2][nXbins * 2];
        scratch = new double[2][nXbins * 2];
        double[] initialValues = new double[] { 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.047681764029297, 0.0886369682387602, 0.1583090316596, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.99454931271489, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924628, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0.0122380386022754, 0.0058389385158292, 0.0026766045152977, 0.0011788613551308, 0, 0, 0, 0, 0 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = scratch[0][i] = 0.0;
            esDs[1][i] = scratch[1][i] = initialValues[i];
        }

        SSEUtils.convolveInPlace(scratch, fY, 1, 1, fftForEandD);
        double[] unnormalizedDs = new double[scratch[1].length];
        for (int i=0, j=0; i<scratch[1].length; i+=2, j++) {
            unnormalizedDs[j] = scratch[1][i];
        }

        everyOtherInPlace(scratch[0], esDs[0], nXbins, 0, 0, 1.0/nXbins); // E's: grabbing real part and scaling by 1/nXbins
        everyOtherInPlace(scratch[1], esDs[1], nXbins, 0, 0, 1.0/nXbins); // D's: grabbing real part and scaling by 1/nXbins

        // System.out.println(Arrays.toString(esDs[0]));
        // System.out.println(Arrays.toString(esDs[1]));

        double[] expectedDsUnnormalized = new double[] { 0.0566871072710171, 0.1290082233227, 0.281301970215054, 0.589342893447835, 1.1863229930959, 2.29444263475987, 4.26373864030162, 7.61277219718144, 13.0597170956171, 21.5259924869188 };
        double[] expectedEs = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double[] expectedDs = new double[] { 0.00118098140147952, 0.00268767131922291, 0.00586045771281363, 0.0122779769468299, 0.0247150623561646, 0.0478008882241639, 0.088827888339617, 0.158599420774613, 0.272077439492024, 0.448458176810809 };

        assertArrayEquals(expectedDsUnnormalized, Arrays.copyOfRange(unnormalizedDs, 0, 10), EPSILON); // comes straight out of convolve
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
        SSEUtils.makeNormalKernelInPlace(fY, (dt * drift), Math.sqrt(dt * diffusion), nXbins, nLeftFlankBins, nRightFlankBins, dx, dt); // normalizes inside already
        fftForEandD.realForwardFull(fY); // first FFT the Normal kernel

        esDs = new double[2][nXbins * 2];
        scratch = new double[2][nXbins * 2];
        double[] initialValues = new double[] { 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.047681764029297, 0.0886369682387602, 0.1583090316596, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.99454931271489, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924628, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0.0122380386022754, 0.0058389385158292, 0.0026766045152977, 0.0011788613551308, 0, 0, 0, 0, 0 };
        for (int i=0; i<esDs[0].length/2; i++) {
            esDs[0][i] = scratch[0][i] = 0.0;
            esDs[1][i] = scratch[1][i] = initialValues[i];
        }

        propagateEandDinXQuaLike(esDs, scratch, fY, nXbins, nLeftFlankBins, nRightFlankBins, nDimensionsE, nDimensionsD, fftForEandD);

        double[] expectedEs = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedDs = new double[] { 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0247150623561646, 0.0478008882241639, 0.088827888339617, 0.158599420774613, 0.272077439492024, 0.448458176810809, 0.710214708266688, 1.0806760140649, 1.57993546382425, 2.21932565164128, 2.99530083804266, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661661, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661661, 5.79327421787606, 4.83940600155166, 3.8841633593627, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.0806760140649, 0.710214708266687, 0.44845817681081, 0.272077439492023, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        assertArrayEquals(expectedEs, Arrays.copyOfRange(esDs[0], 0, nXbins), EPSILON);
        assertArrayEquals(expectedDs, Arrays.copyOfRange(esDs[1], 0, nXbins), EPSILON);
    }
}
