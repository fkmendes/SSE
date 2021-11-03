package SSE;

/*
 * diversitree's (src/quasse-eqs-fftC.c) variables explained:
 *
 * nx = total number of bins for quantitative ch (default = 1024)
 *
 * nkl and nkr:
 * in the normal kernel (giving the pdf of trait value change),
 * specify how many bins on the right (and left, respectively -- and yes, it's reversed!)
 * are non-zero,
 *
 * and
 *
 * in the D and E arrays, they determine how many bins on the left and on the right
 * are not updated by propagate.x()
 *
 * npad = nkl + nkr
 *
 * ndat (my 'nUsefulTraitBins') = nx - npad
 *
 * nd = number of dimensions for each plan (i.e., number of equations per quantitative ch)
 * x = flat, contains all E's followed by all D's (size 2 * nx; my 2D-array 'esDs')
 * wrk = scratch space, used to store and restore x values in left and right padding
 *
 * vars = deep copy of E's, used locally in functions, used locally in propagate_t
 * d = deep copy of x skipping E's (contains D's), used locally in propagate_t
 * dd = deep copy of wrk, only used for D's, used locally in propagate_t to save values and restore them later, also for storing terms used in the math
 */

import org.jtransforms.fft.DoubleFFT_1D;
import java.util.Arrays;

public class SSEUtils {

    private static final double SQRT2PI = Math.sqrt(2 * Math.PI);

    /*
     * Propagates E's and D's, leaving the result in esDs.
     * Equivalent to Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t (src/quasse-eqs-fftC.c)
     * Note: only the first nUsefulTraitBins elements in each row of esDs are used (and in birthRate and deathRate);
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratchAtNode 2D-array for storing/restoring flanking values and other math terms
     * @param   birthRate   birth rates (lambdas), one per quantitative ch useful bin
     * @param   deathRate   death rates (mus), one per quantitative ch useful bin
     * @param   dt  length of time slice over which to propagate
     * @param   nUsefulTraitBins number of useful bins (i.e., excluding the # of flanking bins = npad = nLeftFlankBins + nRightFlankbins) in discretized pdf over change in quantitative ch values
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @return  no return, leaves result in 'esDs' and 'scratch'
     */
    public static void propagateEandDinTQuaSSEInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double[] birthRate, double[] deathRate, double dt, int nUsefulTraitBins, int nDimensionsD) {
        // iterating over all continuous character bins (total # = nx), to update E and D for each of those bins
        for (int i = 0; i < nUsefulTraitBins; ++i) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double netDivRate = ithLambda - ithMu;
            double ithZ = Math.exp(dt * netDivRate);
            double ithE = esDsAtNode[0][i]; // Ask Xia: note that because i starts at 0, we're grabbing E values in left-padding

            double tmp1 = ithMu - ithLambda * ithE;
            double tmp2 = ithZ * (ithE - 1);

            /* Updating E's */
            esDsAtNode[0][i] = (tmp1 + tmp2 * ithMu) / (tmp1 + tmp2 * ithLambda);

            // checking against R
            // System.out.println("ithMu = " + ithMu + " ithLambda = " + ithLambda + " ithE = " + ithE + " ithZ = " + ithZ);
            // System.out.println("esDsAtNode[0][i] = " + esDsAtNode[0][i]);

            // Saving values for updating D below
            tmp1 = (ithLambda - ithMu) / (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE);

            // checking against R
            // System.out.println("i = " + i + " numerator (z * r * r) = " + (ithZ * Math.pow((ithLambda - ithMu), 2)));
            // System.out.println("i = " + i + " denominator (z * lambda - mu + (1-z) * lambda * e0) = " + (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE));

            scratchAtNode[1][i] = ithZ * tmp1 * tmp1;

            // checking against R
            // System.out.println("i = " + i + " (numerator/denominator^2) = " + scratchAtNode[1][i]);
        }

        /* Updating D's */
        // iterating over dimensions in plan to transform (total # = nd = number of equations for D, say, 4 if A C G and T)
        // one quantitative ch: nd = 2 for QuaSSE (one eq for E, one for D), nd=5 for MoSSE (4 eqs for D, one for E)
        // ithDim starts at 1 because esDs[0] are the E's

        for (int ithDim=1; ithDim <= nDimensionsD; ithDim++) {
            // int ithDimStartIdx = nUsefulTraitBins * ithDim; // skipping E's

            // iterating over bins of this dimension
            for (int j=0; j < nUsefulTraitBins; j++) {
                if (esDsAtNode[ithDim][j] < 0) esDsAtNode[ithDim][j] = 0;
                else {
                    // checking against R
                    // System.out.println("Scr j=" + j + " scratchAtNode[0][j] = " + scratchAtNode[0][j] + " scratchAtNode[1][j] = " + scratchAtNode[1][j]);
                    // System.out.println("Bef j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                    esDsAtNode[ithDim][j] *= scratchAtNode[1][j];

                    // checking against R
                    // System.out.println("Aft j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                }
            }
        }
    }

    /*
     * Propagate E's and D's along the x-axis (quantitative trait).
     *
     * The math work (FFTs/iFFts and convolution) is done all in convolveInPlace.
     * This function here basically reorganizes the arrays, keeping track
     * of the flanking bins (not letting those be affected by propagate in x),
     * squashing negative numbers back to 0.0, and setting the last
     * (nLeftFlankBins + nRightFlankBins) to 0.0.
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratch 2D-array for storing/restoring flanking values and other math terms
     * @param   fY  Normal kernel giving the density for changes in value of the quantitative trait
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   nLeftFlankBins  how many bins on the right side of kernel are non-zero
     * @param   nRightFlankBins  how many bins on the left side of kernel are non-zero
     * @param   nDimensionsE number of E equations (dimensions in plan) to solve for each quantitative ch
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @param   fft instance of DoubleFFT_1D that will carry out FFT and inverse-FFT
     */
    public static void propagateEandDinXQuaLike(double[][] esDsAtNode, double[][] scratchAtNode, double[] fY, int nXbins, int nLeftFlankBins, int nRightFlankBins, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {

        // fY is already FFTed (in QuaSSEDistribution -> populatefY()), then FFT scratch, inverse-FFT scratch, result is left in scratch
        convolveInPlace(scratchAtNode, fY, nDimensionsE, nDimensionsD, fft);

        // System.out.println("scratchAtNode[1] after convolve: " + Arrays.toString(scratchAtNode[1]));

        int nItems2Copy = nXbins - nLeftFlankBins - nRightFlankBins;
        // System.out.println("nItems2Copy=" + nItems2Copy);
        double scaleBy = 1.0 / nXbins;

        // move stuff from scratch to esDs, making sure left and right flanks keep the original esDs values (central elements come from scratch)
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            // System.out.println("ithDim = " + ithDim);
            everyOtherInPlace(scratchAtNode[ithDim], esDsAtNode[ithDim], nXbins, nLeftFlankBins, nRightFlankBins, 2, scaleBy); // grabbing real part and scaling by 1/nXbins

            for (int i=0; i<nXbins; ++i) {
                // if negative value, set to 0.0
                // if at last (nLeftFlankBins + nRightFlankBins) items, set to 0.0
                if (esDsAtNode[ithDim][i] < 0.0 || i >= (nItems2Copy-1)) esDsAtNode[ithDim][i] = 0.0;
            }
        }

        // System.out.println("esDsAtNode[1] after convolve and scaling: " + Arrays.toString(esDsAtNode[1]));
    }

    /*
     * Carry out two FFTs, one on the Normal kernel, another on the 2D-array containing the E's and D's.
     * Then carry out inverse FFT on the resulting 2D-array, re-ordering elements so that all real elements occupy the first half.
     * Note: the input of FFT uses just the first half of all elements; the output interdigitates real and complex numbers;
     * the input of inverse-FFT is then the interdigitated array, and the output of inverse-FFT remains interdigitated.
     *
     * Function leaves results in scratch (after moving all real elements to the first half of each row of esDs)
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratchAtNode 2D-array for storing/restoring flanking values and other math terms (from a node)
     * @param   fY  Normal kernel giving the density for changes in value of the quantitative trait
     * @param   nDimensionsE number of E equations (dimensions in plan) to solve for each quantitative ch
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @param   fft instance of DoubleFFT_1D that will carry out FFT and inverse-FFT
     */
    public static void convolveInPlace(double[][] scratchAtNode, double[] fY, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {

        //  fft.realForwardFull(fY); // first FFT the Normal kernel

        // doing E's and D's
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("Before FFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            fft.realForwardFull(scratchAtNode[ithDim]); // FFT for each E and D dimension

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After FFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            for (int i = 0; i < fY.length; i += 2) {
                scratchAtNode[ithDim][i] *= fY[i]; // real part
                scratchAtNode[ithDim][i + 1] *= fY[i]; // complex part
            }

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After FFT and * fY scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            fft.complexInverse(scratchAtNode[ithDim], false); // inverse FFT for each E and D dimension

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After iFFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));
        }

        // looking at things
        // System.out.println(Arrays.toString(scratch[0]));
        // System.out.println(Arrays.toString(scratch[1]));
    }

    /*
     * Builds Normal distribution (yValues) where x is the CHANGE in quantitative trait values
     * (not the quantitative trait values themselves!). Under BM, for example, this
     * kernel is centered at 0.0
     *
     * Note: this kernel (yValues) is later FFT-ed in the likelihood class, and then used in the
     * convolution function, and for reasons I do not fully understand, this bell-shaped kernel\
     * needs to be cuts in two, with the left half being placed at the (right-)tail end of the kernel,
     * and the right half at the (left-)head of the kernel
     *
     * 'nLeftFlankBins' and 'nRightFlankBins' are also the number of bins on the left-
     * side (and right-side, after skipping (nLeftFlankBins + nRightFlankBins))
     * of E and D that are not updated by 'propagateEandDinXQuaLike'
     *
     * @param   yValues (= fY) where the result is left; gives the probability density of a given change in quantitative trait value
     * @param   mean    (= changeInXNormalMean = diversitree's drift * -dt) is the mean expected change in quantitative trait value
     * @param   sd  (= changeInXNormalSd = squared root(diversitree's diffusion * dt)) is the standard deviation of the expected change in quantitative trait value
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   nLeftFlankBins  how many bins on the right side of kernel are non-zero
     * @param   nRightFlankBins  how many bins on the left side of kernel are non-zero
     * @param   dx  size of each bin
     */
    public static void makeNormalKernelInPlace(double[] yValues, double mean, double sd, int nXbins, int nLeftFlankBins, int nRightFlankBins, double dx) {

        double total = 0.0;

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

        // System.out.println("fY = " + Arrays.toString(yValues));
    }

    /*
     * (DEPRECATED: will be done in makeNormalKernelInPlace)
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

    /*
     * Return equivalent to R's dnorm(x, mean, sd)
     *
     * @param   x   random variable value
     * @param   mean    Mean of normal distribution
     * @param   sd  Standard deviation of normal distribution
     * @return  density of x under Normal distribution
     */
    public static double getNormalDensity(double x, double mean, double sd) {
        double x0 = x - mean;
        return Math.exp(-x0 * x0 / (2 * sd * sd)) / (sd * SQRT2PI);
    }

    /*
     * Populate toArray in place, getting every other element from
     * fromArray.
     *
     * @param   inArray source array
     * @param   outArray result array to receive every other element from source array
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   skipFirstN  number of discrete quantitative trait bins on the left-side of E and D that are not affected by 'propagateEandDinXQuaLike'
     * @param   skipLastN  number of discrete quantitative trait bins on the right-side of E and D that are not affected by 'propagateEandDinXQuaLike'
     * @param   scaleBy will scale every other element by this
     */
    public static void everyOtherInPlace(double[] fromArray, double[] toArray, int nXbins, int skipFirstN, int skipLastN, int everyOther, double scaleBy) {
        // move these checks to initAndValidate later
        if (skipFirstN == 0 && skipLastN == 0) skipLastN = -1; // this should only happen in debugging, where there are no elements to skip at the start or end
        if (fromArray.length != toArray.length) throw new RuntimeException("Arrays of different size. Exiting...");
        if (toArray.length/nXbins != 2.0) throw new RuntimeException("Size of arrays must be double the number of quantitative ch bins. Exiting...");
        // TODO: instead of multiple of 4, it's actually a power of 2
        if ((2 * nXbins) % 4 != 0.0) throw new RuntimeException("Number of quantitative character bins must be a multiple of 4. Exiting...");
        if ((nXbins - skipLastN - (skipFirstN + skipLastN)) <= skipFirstN) throw new RuntimeException("There are no useful bins left after pushing left and right flanks to the end, on top of right flank. Exiting...");

        for (int i=skipFirstN * 2, j=skipFirstN; i <= (fromArray.length-2); i+=everyOther, j++) {
            // if (j < (nXbins - skipLastN - (skipFirstN + skipLastN) - 1)) { // working
            if (j < (nXbins - skipLastN - (skipFirstN + skipLastN) - 1)) {
                toArray[j] = fromArray[i] * scaleBy;
            }
        }
    }

    public static void hiToLoTransferInPlace(double[] fromArray, double[] toArray, int[] idxs4Transfer) {
        int i = 0;
        for (int j: idxs4Transfer) {
            toArray[i] = fromArray[j];
            i++;
        }
    }

    public static double calculateNormalizationFactor(double[] dS, double binSize) {
        double normalizationFactorFromDs = 0.0;
        for (int i=0; i<dS.length; ++i) normalizationFactorFromDs += dS[i];

        // debugging
        // System.out.println("normalizationFactorFromDs inside calculateNormalizationFactor = " + normalizationFactorFromDs);
        // System.out.println("binSize inside calculateNormalizationFactor = " + binSize);

        normalizationFactorFromDs *= binSize;

        // debugging
        // System.out.println("normalizationFactorFromDs * binSize = " + normalizationFactorFromDs);

        return normalizationFactorFromDs;
    }
}
