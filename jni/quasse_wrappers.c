#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <math.h>
#include <jni.h>
#include <string.h>
#include <assert.h>

#include "rfftw.h"
#include "quasse-eqs-fftC.h"
#include "mosse-eqs-fftC.h"

#define R_D__0 (give_log ? -INFINITY : 0.)
#define M_LN_SQRT_2PI 0.918938533204672741780329736406
#define M_1_SQRT_2PI 0.398942280401432677939946059934

/*
 * Allocate memory for quasse_stash, populate it, return pointer to quasse_stash
 */
quasse_constants_stash* make_quasse_stash(
    int nXbinsLo,
    int nXbinsHi,
    int nUsefulXbinsLo,
    int nUsefulXbinsHi,
    int nLeftFlankBinsLo,
    int nLeftFlankBinsHi,
    int nRightFlankBinsHi,
    int nRightFlankBinsLo,
    double dxBin,
    double *xLo_ptr,
    double *xHi_ptr,
    int nDimensions,
    int *nDnEs_ptr,
    int flags) {

    // nDimensions = total number of required (multi-dimensional) FFT plans, one per each D and each E
    quasse_constants_stash *quasse_cts_stash_ptr = calloc(1, sizeof(quasse_constants_stash));

    int nYLo = (((int)floor(nXbinsLo / 2)) + 1);
    int nYHi = (((int)floor(nXbinsHi / 2)) + 1);
    int i, max_nd = 1;
    for (i=0; i < nDimensions; i++) {
    if (nDnEs_ptr[i] > max_nd) max_nd = nDnEs_ptr[i];
    }

    quasse_cts_stash_ptr->nYLo = nYLo;
    quasse_cts_stash_ptr->nYHi = nYHi;
    quasse_cts_stash_ptr->max_nd = max_nd;

    // basic constant members
    quasse_cts_stash_ptr->nXbinsLo = nXbinsLo;
    quasse_cts_stash_ptr->nXbinsHi = nXbinsHi;
    quasse_cts_stash_ptr->nUsefulXbinsLo = nUsefulXbinsLo;
    quasse_cts_stash_ptr->nUsefulXbinsHi = nUsefulXbinsHi;
    quasse_cts_stash_ptr->dxBin = dxBin;

    // fft-related constant members
    quasse_cts_stash_ptr->nDnEs_ptr = nDnEs_ptr; // pointer to array of doubles
    quasse_cts_stash_ptr->nDimensions = nDimensions;

    // are x and y just "working areas" for FFTing?
    // they are required when setting up fft plans
    // quasse_cts_stash_ptr->fft_x_input_lo = fftw_malloc(max_nd * nXbinsLo * sizeof(double));
    // quasse_cts_stash_ptr->fft_x_input_hi = fftw_malloc(max_nd * nXbinsHi * sizeof(double));
    // quasse_cts_stash_ptr->fft_y_input_lo = fftw_malloc(max_nd * (nYLo + 1) * sizeof(fftw_complex));
    // quasse_cts_stash_ptr->fft_y_input_hi = fftw_malloc(max_nd * (nYHi + 1) * sizeof(fftw_complex));

    quasse_cts_stash_ptr->xLo_ptr = xLo_ptr; // x-ruler in low resolution
    quasse_cts_stash_ptr->xHi_ptr = xHi_ptr; // x-ruler in high resolution

    quasse_cts_stash_ptr->fftLo = (rfftw_plan_real **)calloc(nDimensions, sizeof(rfftw_plan_real *));
    quasse_cts_stash_ptr->fftHi = (rfftw_plan_real **)calloc(nDimensions, sizeof(rfftw_plan_real *));

    quasse_cts_stash_ptr->zLo = (double *)calloc(nXbinsLo, sizeof(double));
    quasse_cts_stash_ptr->wrkLo = (double *)calloc(nXbinsLo, sizeof(double));
    quasse_cts_stash_ptr->wrkdLo = (double *)calloc(max_nd * nXbinsLo, sizeof(double));
    quasse_cts_stash_ptr->zHi = (double *)calloc(nXbinsHi, sizeof(double));
    quasse_cts_stash_ptr->wrkHi = (double *)calloc(nXbinsHi, sizeof(double));
    quasse_cts_stash_ptr->wrkdHi = (double *)calloc(max_nd * nXbinsHi, sizeof(double));

    // making 'nDimensions' fft plans, for both low and high resolution
    for (i=0; i < nDimensions; i++) {
        quasse_cts_stash_ptr->fftLo[i] =
            make_rfftw_plan_real(nDnEs_ptr[i], nXbinsLo, DIR_COLS, quasse_cts_stash_ptr->fft_x_input_lo, quasse_cts_stash_ptr->fft_y_input_lo, flags);

        quasse_cts_stash_ptr->fftHi[i] =
            make_rfftw_plan_real(nDnEs_ptr[i], nXbinsHi, DIR_COLS, quasse_cts_stash_ptr->fft_x_input_hi, quasse_cts_stash_ptr->fft_y_input_hi, flags);
    }

  /* Brownian kernel */
  // quasse_stash_ptr->kern_x = fftw_malloc(nxBin * sizeof(double));
  // quasse_stash_ptr->kern_y = fftw_malloc((ny + 1) * sizeof(fftw_complex));
//  quasse_stash_ptr->kernel =
//      make_rfftw_plan_real(1, nxBin, DIR_COLS, obj->kern_x, obj->kern_y, flags);

  return quasse_cts_stash_ptr;
}



JNIEXPORT jlong JNICALL makeQuasseConstantsStash_JNIWrapper(
    JNIEnv *env,
    jobject thisObject,
    jint nXbinsLo,
    jint nXbinsHi,
    jint nUsefulXbinsLo,
    jint nUsefulXbinsHi,
    jint nLeftFlankBinsLo,
    jint nLeftFlankBinsHi,
    jint nRightFlankBinsLo,
    jint nRightFlankBinsHi,
    jdouble dxBin,
    jdoubleArray xLo_ptr,
    jdoubleArray xHi_ptr,
    jintArray nDnEs_ptr,
    jint flags) {

    jsize nDimensions = (int)((*env)->GetArrayLength(env, nDnEs_ptr));
    assert(nDimensions);

    // allocating memory for C values we will copy content to
    int *nDnEs_ptr_C = malloc(sizeof(int) * nDimensions);
    double *xLo_ptr_C = malloc(sizeof(double) * nXbinsLo);
    double *xHi_ptr_C = malloc(sizeof(double) * nXbinsHi);

    assert(nDnEs_ptr_C);
    assert(xLo_ptr_C);
    assert(xHi_ptr_C);

    // grab Java array content
    jint *const_int_array_content = (*env)->GetIntArrayElements(env, nDnEs_ptr, 0);
    jdouble *const_double_array_content1 = (*env)->GetDoubleArrayElements(env, xLo_ptr, 0);
    jdouble *const_double_array_content2 = (*env)->GetDoubleArrayElements(env, xHi_ptr, 0);

    // copy java arrays to c arrays
    memcpy(nDnEs_ptr_C, const_int_array_content, sizeof(int) * nDimensions);
    (*env)->ReleaseIntArrayElements(env, nDnEs_ptr, const_int_array_content, 0);
    memcpy(xLo_ptr_C, const_double_array_content1, sizeof(double) * nXbinsLo);
    (*env)->ReleaseDoubleArrayElements(env, xLo_ptr_C, const_double_array_content1, 0);
    memcpy(xHi_ptr_C, const_double_array_content2, sizeof(double) * nXbinsHi);
    (*env)->ReleaseDoubleArrayElements(env, xHi_ptr_C, const_double_array_content2, 0);

    // populate quasse_constants_stash
    quasse_constants_stash *quasse_cts_stash_ptr =
        make_quasse_stash(nXbinsLo, nXbinsHi, nUsefulXbinsLo, nUsefulXbinsHi, nLeftFlankBinsLo, nLeftFlankBinsHi, nRightFlankBinsLo, nRightFlankBinsHi, dxBin, double *xLo_ptr, double *xHi_ptr, nDimensions, nDnEs_ptr_C, flags);

  return (jlong)(quasse_cts_stash_ptr);
};




JNIEXPORT jdoubleArray JNICALL populateFYandIntegrateInPlace_JNIWrapper(
    JNIEnv *env,
    jobject thisObject,
    jlong quasse_stash_ptr,
    jdoubleArray esDs, jdoubleArray lambda, jdoubleArray mu,
    jdouble drift,
    jdouble diffusion,
    jint nIntervals,
    jdouble aDt,
    jboolean lowRes) {

    quasse_constants_stash *quasse_cts_stash_ptr = (quasse_constants_stash *)(quasse_stash_ptr);
    assert(quasse_cts_stash_ptr);

    // nkl, nkr, ndat
    int nLeftFlankBins, nRightFlankBins, nUsefulXbins;

    if (lowRes) {
        nLeftFlankBins = quasse_cts_stash_ptr->nLeftFlankBinsLo;
        nRightFlankBins = quasse_cts_stash_ptr->nRightFlankBinsLo;
        nUsefulXbins = quasse_cts_stash_ptr->nUsefulXbinsLo;
    } else {
        nLeftFlankBins = quasse_cts_stash_ptr->nLeftFlankBinsHi;
        nRightFlankBins = quasse_cts_stash_ptr->nRightFlankBinsHi;
        nUsefulXbins = quasse_cts_stash_ptr->nUsefulXbinsHi;
    }

    int nIntervals_C = nIntervals;
    double aDt_C = aDt;
    double drift_C = drift;
    double diffusion_C = diffusion;
};