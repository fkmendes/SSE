package SSE;

/*
 * diversitree's (src/quasse-eqs-fftC.c) variables explained:
 *
 * nx = total number of bins for quantitative ch (default = 1024)
 *
 * nkl and nkr:
 *
 * in the normal kernel (giving the pdf of trait value change),
 * specify how many bins on the right (and left, respectively -- and yes, it's reversed!)
 * are non-zero,
 * and
 * in the D and E arrays, they determine how many bins on the left and on the right
 * are not updated by propagate.x()
 *
 * ndat = nx - padding (useful number of bins for quantitative ch, my 'nUsefulTraitBins')
 * nd = number of dimensions for each plan (i.e., number of equations per quantitative ch)
 * x = contains all E's followed by all D's (size 2 * nx; my 'esDs')
 * wrk = scratch space, used to store and restore x values in left and right padding
 *
 * vars = deep copy of E's, used locally in functions, used locally in propagate_t
 * d = deep copy of x skipping E's (contains D's), used locally in propagate_t
 * dd = deep copy of wrk, only used for D's, used locally in propagate_t to save values and restore them later, also for storing terms used in the math
 */

import org.jtransforms.fft.DoubleFFT_1D;

public class SSEUtils {

    private static final double SQRT2PI = Math.sqrt(2 * Math.PI);

    /*
     * Propagates E's and D's, leaving the result in esDs.
     * Equivalent to Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t (src/quasse-eqs-fftC.c)
     *
     * @param   esDs    flat array containing E's followed by D's, for (i) all quantitative ch values, (ii) all equations per ch
     * @param   scratch array for storing/restoring flanking values and other math terms
     * @param   birthRate   birth rates (lambdas), one per quantitative ch
     * @param   birthRate   death rates (mus), one per quantitative ch
     * @param   dt  length of time slice over which to propagate
     * @param   nUsefulTraitBins number of useful (i.e., excluding the # of flanking bins = npad) bins in discretized pdf over quantitative ch values (e.g., morphological ch's, or subst rate)
     * @param   nDimensions number of equations (dimensions in plan) to solve for each quantitative ch
     * @return  no return, leaves result in 'esDs' and 'scratch'
     */
    public static void propagateEandDinTQuaLike(double[][] esDs, double[] scratch, double[] birthRate, double[] deathRate, double dt, int nUsefulTraitBins, int nDimensionsE, int nDimensionsD) {

        // iterating over all continuous character bins (total # = nx), to update E and D for each of those bins
        for (int i = 0; i < nUsefulTraitBins; ++i) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double netDivRate = ithLambda - ithMu;
            double ithZ = Math.exp(dt * netDivRate);
            double ithE = esDs[0][i]; // Ask Xia: note that because i starts at 0, we're grabbing E values in left-padding

            double tmp1 = ithMu - ithLambda * ithE;
            double tmp2 = ithZ * (ithE - 1);

            /* Updating E's */
            esDs[0][i] = (tmp1 + tmp2 * ithMu) / (tmp1 + tmp2 * ithLambda);

            // Saving values for updating D below
            tmp1 = (ithLambda - ithMu) / (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE);
            scratch[i] = ithZ * tmp1 * tmp1;
        }

        /* Updating D's */
        // iterating over dimensions in plan to transform (total # = nd = number of equations for D, say, 4 if A C G and T)
        // one quantitative ch: nd = 2 for QuaSSE (one eq for E, one for D), nd=5 for MoSSE (4 eqs for D, one for E)
        // ithDim starts at 1 because we skip the E's
        for (int ithDim = 1; ithDim < nDimensionsD; ithDim++) {
            // int ithDimStartIdx = nUsefulTraitBins * ithDim; // skipping E's

            // iterating over bins of this dimension
            for (int j = 0; j < nUsefulTraitBins; j++) {
                // int dIdx = ithDimStartIdx + j;
                // int dIdx = ithDim + j;

                if (esDs[ithDim][j] < 0) esDs[ithDim][j] = 0;
                else esDs[ithDim][j] *= scratch[j]; // Ask Xia: why no dimensions here for scratch?
            }
        }
    }

    /*
     * Propagates E's and D's, leaving the result in esDs.
     * Equivalent to Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t (src/quasse-eqs-fftC.c)
     *
     * @param   esDs    flat array containing E's followed by D's, for (i) all quantitative ch values, (ii) all equations per ch
     */
    public static void propagateEandDinXQuaSSE(double[] esDs, double[] fY) {

    }

    public static void makeNormalKernelInPlace(double[] yValues, double drift, double diffusion, int nXbins, int nLeftFlankBins, int nRightFlankBins, double dx, double dt) {
        double total = 0.0;
        double mean = -dt * drift;
        double sd = Math.sqrt(dt * diffusion);

        double x = 0.0;
        for (int i = 0; i <= nRightFlankBins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd); // in diversitree's C code, this is further multiplied by dx (think this is unnecessary, b/c it doesn't change the result!)
            total += yValues[i];
            x += dx;
        }

        for (int i = nRightFlankBins + 1; i < nXbins - nLeftFlankBins; i++) yValues[i] = 0;

        x = -nLeftFlankBins * dx;
        for (int i = nXbins - nLeftFlankBins; i < nXbins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd);
            total += yValues[i];
            x += dx;
        }

        // System.out.println("Unormalized fY = " + Arrays.toString(yValues));
        // System.out.println("Total=" + total);

        for (int i = 0; i <= nRightFlankBins; i++) yValues[i] /= total;
        for (int i = (nXbins - nLeftFlankBins); i < nXbins; i++) yValues[i] /= total;
    }

    public static void convolve(double[][] esDs, double[][] scratch, double[] fY, int nXbins, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {
        double oneOverNXbins = 1.0 / nXbins;
        fft.realForwardFull(fY); // first FFT the Normal kernel

        // doing E's and D's
        for (int ithDim = 1; ithDim < (nDimensionsE+nDimensionsD); ithDim++) {
            fft.realForwardFull(scratch[ithDim]); // FFT for each E and D dimension

            for (int i=0; i < fY.length; i += 2) {
                scratch[ithDim][i] *= fY[i]; // real part
                scratch[ithDim][i+1] *= fY[i]; // complex part
            }

            fft.complexInverse(scratch[ithDim], false); // inverse FFT for each E and D dimension

            everyOtherInPlace(scratch[ithDim], esDs[ithDim], true, oneOverNXbins); // grabbing real part and scaling by 1/nXbins
        }
    }

    /*
     * Not necessary (will be done in makeNormalKernelInPlace)
     */
    public static void normalizeArray(double[] doubleArray) {
        double theSum = 0.0;
        for (int i = 0; i < doubleArray.length; i++) {
            theSum += doubleArray[i];
        }

        for (int i = 0; i < doubleArray.length; i++) {
            doubleArray[i] = doubleArray[i] / theSum;
        }
    }

    public static double getNormalDensity(double x, double mean, double sd) {
        double x0 = x - mean;
        return Math.exp(-x0 * x0 / (2 * sd * sd)) / (sd * SQRT2PI);
    }

    public static void everyOtherInPlace(double[] inArray, double[] outArray, boolean odd, double scaleBy) {
        int k;
        for (int i=0, j=0; i < inArray.length; i+=2, j++) {
            k = i;
            if (!odd) k++;
            if (j < outArray.length) {
                outArray[j] = inArray[k] * scaleBy;
            }
        }
    }
}
