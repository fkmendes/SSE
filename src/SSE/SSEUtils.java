package SSE;

/*
 * diversitree's (src/quasse-eqs-fftC.c) variables explained:
 *
 * nx = total number of bins for quantitative ch (default = 1024)
 * ndat = nx - padding (useful number of bins for quantitative ch, my 'nUsefulTraitBins')
 * nd = number of dimensions for each plan (i.e., number of equations per quantitative ch)
 * x = contains all E's followed by all D's (size 2 * nx; my 'esDs')
 * wrk = scratch space, used to store and restore x values in left and right padding
 *
 * vars = deep copy of E's, used locally in functions, used locally in propagate_t
 * d = deep copy of x skipping E's (contains D's), used locally in propagate_t
 * dd = deep copy of wrk, only used for D's, used locally in propagate_t to save values and restore them later, also for storing terms used in the math
 */

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
    public static void propagateEandDinTQuaLike(double[] esDs, double[] scratch, double[] birthRate, double[] deathRate, double dt, int nUsefulTraitBins, int nDimensions) {

        // iterating over all continuous character bins (total # = nx), to update E and D for each of those bins
        for (int i = 0; i < nUsefulTraitBins; ++i) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double netDivRate = ithLambda - ithMu;
            double ithZ = Math.exp(dt * netDivRate);
            double ithE = esDs[i]; // Ask Xia: note that because i starts at 0, we're grabbing E values in left-padding

            double tmp1 = ithMu - ithLambda * ithE;
            double tmp2 = ithZ * (ithE - 1);

            /* Updating E's */
            esDs[i] = (tmp1 + tmp2 * ithMu) / (tmp1 + tmp2 * ithLambda);

            // Saving values for updating D below
            tmp1 = (ithLambda - ithMu) / (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE);
            scratch[i] = ithZ * tmp1 * tmp1;
        }

        /* Updating D's */
        // iterating over dimensions in plan to transform (total # = nd = number of equations for D, say, 4 if A C G and T)
        // one quantitative ch: nd = 2 for QuaSSE (one eq for E, one for D), nd=5 for MoSSE (4 eqs for D, one for E)
        for (int ithDim = 1; ithDim < nDimensions; ithDim++) {
            int ithDimStartIdx = nUsefulTraitBins * ithDim; // skipping E's

            // iterating over bins of this dimension
            for (int j = 0; j < nUsefulTraitBins; j++) {
                int dIdx = ithDimStartIdx + j;

                if (esDs[dIdx] < 0) esDs[dIdx] = 0;
                else esDs[dIdx] *= scratch[j]; // Ask Xia: why no dimensions here for scratch?
            }
        }
    }

    public static void makeNormalKernelInPlace(double[] yValues, double drift, double diffusion, int nXbins, int nLeftFlankBins, int nRightFlankBins, double dx, double dt) {
        double total = 0.0;
        double mean = -dt * drift;
        double sd = Math.sqrt(dt * diffusion);

        double x = 0.0;
        for (int i=0; i <= nRightFlankBins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd); // in diversitree's C code, this is further multiplied by dx (might need to do this when fft-ing?)
            total += yValues[i];
            x += dx;
        }

        for (int i=nRightFlankBins+1; i < nXbins - nLeftFlankBins; i++) yValues[i] = 0;

        x = -nLeftFlankBins * dx;
        for (int i=nXbins - nLeftFlankBins; i < nXbins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd);
            total += yValues[i];
            x += dx;
        }

        // System.out.println("Unormalized fY = " + Arrays.toString(yValues));
        // System.out.println("Total=" + total);

        for (int i=0; i <= nRightFlankBins; i++) yValues[i] /= total;
        for (int i = (nXbins - nLeftFlankBins); i < nXbins; i++) yValues[i] /= total;
    }

    /*
     * Not necessary (will be done in makeNormalKernelInPlace)
     */
    public static void normalizeArray(double[] doubleArray) {
        double theSum = 0.0;
        for (int i=0; i<doubleArray.length; i++) {
            theSum += doubleArray[i];
        }

        for (int i=0; i<doubleArray.length; i++) {
            doubleArray[i] = doubleArray[i] / theSum;
        }
    }


    public static double getNormalDensity(double x, double mean, double sd) {
        double x0 = x - mean;
        return Math.exp(-x0 * x0 / (2 * sd * sd)) / (sd * SQRT2PI);
    }
}
