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

JNIEXPORT jlong JNICALL Java_mosse_MosseDistribution_makeMosseFFT(
    JNIEnv *env, jobject thisObject, jint nx, jdouble dx, jintArray array_nd,
    jint flags) {

  jsize n_fft = (int)((*env)->GetArrayLength(env, array_nd));
  assert(n_fft);

  jint *const_array_body = (*env)->GetIntArrayElements(env, array_nd, 0);
  int *local_nd = malloc(sizeof(int) * n_fft);
  assert(local_nd);
  // copy java array to c array
  memcpy(local_nd, const_array_body, sizeof(int) * n_fft);
  (*env)->ReleaseIntArrayElements(env, array_nd, const_array_body, 0);

  mosse_fft *obj = make_mosse_fft(n_fft, nx, dx, local_nd, flags);

  return (jlong)(obj);
}

// JNIEXPORT void JNICALL Java_mosse_MosseDistribution_propagateT(JNIEnv *env, jobject thisObject, jlong obj_ptr) {
//    mosse_fft *obj = (mosse_fft *)(obj_ptr);
//    int idx = 0;
//    propagate_t_mosse(obj, idx);
// }

JNIEXPORT jdoubleArray JNICALL Java_mosse_MosseDistribution_getX(JNIEnv *env, jobject thisObject, jlong obj_ptr) {
   mosse_fft *obj = (mosse_fft*)(obj_ptr);
   // convert c array to java array and return it
   int size = (obj->max_nd) * (obj->nx);
   jdoubleArray result = (*env)->NewDoubleArray(env, size);
   // move c array structure to java structure 
   (*env)->SetDoubleArrayRegion(env, result, 0, size, obj->x); 
   return result;

}

JNIEXPORT jdoubleArray JNICALL Java_mosse_MosseDistribution_doIntegrateMosse(
    JNIEnv *env, jobject thisObject, jlong mosse_ptr, jdoubleArray vars, jdoubleArray lambda, jdoubleArray mu,
    jdouble drift, jdouble diffusion, jdoubleArray Q, jint nt, jdouble dt, jint pad_left, jint pad_right) {
  
  mosse_fft *obj = (mosse_fft *)(mosse_ptr);
  assert(obj);

  int nkl = pad_left;
  int nkr = pad_right;
  // int ndat = LENGTH(lambda);
  int ndat = (int)((*env)->GetArrayLength(env, lambda));
  double c_dt = dt;
  int c_nt = nt;
  // double *c_lambda = REAL(lambda);
  // double *c_mu = REAL(mu);
  double c_drift = drift;
  double c_diffusion = diffusion;
  int i, idx, nd;

  int len_vars = (int)((*env)->GetArrayLength(env, vars));
  nd = len_vars / obj->nx;

  // setup c lambda array
  jsize n_lambda = (int)((*env)->GetArrayLength(env, lambda));
  jdouble *const_array_body_lambda = (*env)->GetDoubleArrayElements(env, lambda, 0);
  double *c_lambda = malloc(sizeof(double) * n_lambda);
  assert(c_lambda);
  memcpy(c_lambda, const_array_body_lambda, sizeof(double) * n_lambda);
  (*env)->ReleaseDoubleArrayElements(env, lambda, const_array_body_lambda, 0);

  // setup c mu array
  jsize n_mu = (int)((*env)->GetArrayLength(env, mu));
  jdouble *const_array_body_mu = (*env)->GetDoubleArrayElements(env, mu, 0);
  double *c_mu = malloc(sizeof(double) * n_mu);
  assert(c_mu);
  memcpy(c_mu, const_array_body_mu, sizeof(double) * n_mu);
  (*env)->ReleaseDoubleArrayElements(env, mu, const_array_body_mu, 0);

  // setup c vars array
  jsize n_vars = (int)((*env)->GetArrayLength(env, vars));
  jdouble *const_array_body_vars = (*env)->GetDoubleArrayElements(env, vars, 0);
  double *c_vars = malloc(sizeof(double) * n_vars);
  assert(c_vars);
  memcpy(c_vars, const_array_body_vars, sizeof(double) * n_vars);
  (*env)->ReleaseDoubleArrayElements(env, vars, const_array_body_vars, 0);

  // setup Q c array
  jsize n_Q = (int)((*env)->GetArrayLength(env, Q));
  jdouble *const_array_body_Q = (*env)->GetDoubleArrayElements(env, Q, 0);
  double *c_Q = malloc(sizeof(double) * n_Q);
  assert(c_Q);
  memcpy(c_Q, const_array_body_Q, sizeof(double) * n_Q);
  (*env)->ReleaseDoubleArrayElements(env, Q, const_array_body_Q, 0);

  idx = lookup(nd, obj->nd, obj->n_fft);
  if (idx < 0) {
    printf("Failed to find nd = %d\n", nd);
    abort();
  }

  qf_copy_x_mosse(obj, c_vars, nd, 1);

  obj->lambda = c_lambda;
  obj->mu = c_mu;
  obj->Q = c_Q;

  for (i = 0; i < ndat; i++)
    obj->z[i] = exp(c_dt * (c_lambda[i] - c_mu[i]));

  qf_setup_kern_mosse(obj, c_drift, c_diffusion, c_dt, nkl, nkr);

  do_integrate_mosse(obj, c_nt, idx);

  obj->lambda = NULL;
  obj->mu = NULL;

  // PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  int size = obj->nx * nd;
  double *result = (double*)malloc(sizeof(double) * size); 
  qf_copy_x_mosse(obj, result, nd, 0); // copy obj to result

  // copy to java array
  jdoubleArray j_result = (*env)->NewDoubleArray(env, size);
  (*env)->SetDoubleArrayRegion(env, j_result, 0, size, result);

  // free memory
  free(c_lambda);
  free(c_mu);
  free(c_vars);
  free(c_Q);

  return j_result;
}

int lookup(int x, int *v, int len) {
  int i, idx=-1;
  for ( i = 0; i < len; i++ )
    if ( v[i] == x ) {
      idx = i;
      break;
    }

  return idx;
}

double dnorm(double x, double mu, double sigma, int give_log) {

  if (isnan(x) || isnan(mu) || isnan(sigma)) {
    return x + mu + sigma;
  }
  if (!isfinite(sigma)) {
    return R_D__0;
  }
  if (!isfinite(x) && mu == x) {
    return NAN; /* x-mu is NaN */
  }
  if (sigma <= 0) {
    if (sigma < 0) {
      printf("Error!\n");
      return NAN;
    }
    /* sigma == 0 */
    return (x == mu) ? INFINITY : R_D__0;
  }
  x = (x - mu) / sigma;

  if (!isfinite(x)) {
    return R_D__0;
  }
  return (give_log ? -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma))
                   : M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma);
  /* M_1_SQRT_2PI = 1 / sqrt(2 * pi) */
}

rfftw_plan_real *make_rfftw_plan_real(int nd, int nx, int dir, double *x,
                                      fftw_complex *y, int flags) {
  rfftw_plan_real *obj = (rfftw_plan_real *)calloc(1, sizeof(rfftw_plan_real));
  int ny = ((int)floor(nx / 2)) + 1;
  int xstride, xdist, ystride, ydist;

  if (dir == DIR_COLS) {
    xstride = 1;
    ystride = 1;
    xdist = nx;
    ydist = ny;
  } else { /* if ( dir == DIR_ROWS ) { // assume ROWS silently */
    xstride = nd;
    ystride = nd;
    xdist = 1;
    ydist = 1;
  }

  obj->nd = nd;
  obj->nx = nx;
  obj->ny = ny;
  obj->x = x;
  obj->y = y;
  obj->dir = dir;

  /* Consider FFTW_PATIENT, or even FFTW_EXHAUSTIVE */
  obj->plan_f = fftw_plan_many_dft_r2c(1, &nx, nd, obj->x, NULL, xstride, xdist,
                                       obj->y, NULL, ystride, ydist, flags);
  obj->plan_b = fftw_plan_many_dft_c2r(1, &nx, nd, obj->y, NULL, ystride, ydist,
                                       obj->x, NULL, xstride, xdist, flags);
  return obj;
}

void convolve(rfftw_plan_real *obj, fftw_complex *fy) {
  int nx = obj->nx, ny = obj->ny, nd = obj->nd, i, j;
  int nxd = nx * nd;
  double *x = obj->x;
  fftw_complex *y = obj->y;

  fftw_execute(obj->plan_f);

  for (i = 0; i < nd; i++)
    for (j = 0; j < ny; j++, y++)
      (*y) *= fy[j];

  fftw_execute(obj->plan_b);

  for (i = 0; i < nxd; i++)
    x[i] /= nx;
}

/* This does the memory allocation and plans the FFT transforms */
mosse_fft *make_mosse_fft(int n_fft, int nx, double dx, int *nd, int flags) {
  mosse_fft *obj = calloc(1, sizeof(mosse_fft));
  int ny = (((int)floor(nx / 2)) + 1);
  int i, max_nd = 1;
  for (i = 0; i < n_fft; i++)
    if (nd[i] > max_nd)
      max_nd = nd[i];

  obj->n_fft = n_fft;
  obj->nx = nx;
  obj->ny = ny;
  obj->dx = dx;
  obj->nd = nd;

  obj->x = fftw_malloc(max_nd * nx * sizeof(double));
  obj->y = fftw_malloc(max_nd * (ny + 1) * sizeof(fftw_complex));
  obj->max_nd = max_nd;

  obj->z = (double *)calloc(nx, sizeof(double));
  obj->wrk = (double *)calloc(nx, sizeof(double));
  obj->wrkd = (double *)calloc(max_nd * nx, sizeof(double));

  obj->fft = (rfftw_plan_real **)calloc(n_fft, sizeof(rfftw_plan_real *));

  for (i = 0; i < n_fft; i++) {
    obj->fft[i] =
        make_rfftw_plan_real(nd[i], nx, DIR_COLS, obj->x, obj->y, flags);
  }

  /* Brownian kernel */
  obj->kern_x = fftw_malloc(nx * sizeof(double));
  obj->kern_y = fftw_malloc((ny + 1) * sizeof(fftw_complex));
  obj->kernel =
      make_rfftw_plan_real(1, nx, DIR_COLS, obj->kern_x, obj->kern_y, flags);

  return obj;
}

JNIEXPORT void JNICALL Java_mosse_MosseDistribution_mosseFinalize(JNIEnv *env, jobject thisObject, jlong obj_ptr) {
  mosse_fft *obj = (mosse_fft *)(obj_ptr);

  int i;
  /* Rprintf("Cleaning up\n"); */

  for (i = 0; i < obj->n_fft; i++) {
    fftw_destroy_plan(obj->fft[i]->plan_f);
    fftw_destroy_plan(obj->fft[i]->plan_b);
  }
  free(obj->fft);
  free(obj->nd);

  fftw_free(obj->x);
  fftw_free(obj->y);

  free(obj->z);
  free(obj->wrk);
  free(obj->wrkd);

  fftw_destroy_plan(obj->kernel->plan_f);
  fftw_destroy_plan(obj->kernel->plan_b);

  fftw_free(obj->kern_x);
  fftw_free(obj->kern_y);

  free(obj);
}

void qf_copy_x_mosse(mosse_fft *obj, double *x, int nd, int copy_in) {
  int i, n = obj->nx * nd;
  double *fft_x = obj->x;
  if (copy_in) {
    for (i = 0; i < n; i++) {
      fft_x[i] = x[i];
    }
  } else {
    for (i = 0; i < n; i++) {
      x[i] = fft_x[i];
    }
  }
}

void qf_setup_kern_mosse(mosse_fft *obj, double drift, double diffusion,
                         double dt, int nkl, int nkr) {
  const int nx = obj->nx;
  int i;
  double x, *kern_x = obj->kern_x, tot = 0, dx = obj->dx;
  double mean, sd;

  obj->nkl = nkl;
  obj->nkr = nkr;
  obj->npad = nkl + 1 + nkr;
  obj->ndat = nx - obj->npad;
  obj->drift = drift;
  obj->diffusion = diffusion;

  tot = 0;
  mean = -dt * drift;
  sd = sqrt(dt * diffusion);

  for (i = 0, x = 0; i <= nkr; i++, x += dx)
    tot += kern_x[i] = dnorm(x, mean, sd, 0) * dx;
  for (i = nkr + 1; i < nx - nkl; i++)
    kern_x[i] = 0;
  for (i = nx - nkl, x = -nkl * dx; i < nx; i++, x += dx)
    tot += kern_x[i] = dnorm(x, mean, sd, 0) * dx;

  for (i = 0; i <= nkr; i++)
    kern_x[i] /= tot;
  for (i = nx - nkl; i < nx; i++)
    kern_x[i] /= tot;

  fftw_execute(obj->kernel->plan_f);
}

void do_integrate_mosse(mosse_fft *obj, int nt, int idx) {
  int i, nkl = obj->nkl;
  for (i = 0; i < nt; i++) {
    propagate_t_mosse(obj, idx);
    propagate_x_mosse(obj, idx);
    if (isnan(obj->x[nkl])) {
      printf("Integration failure at step %d\n", i);
      abort();
    }
  }
}

/* Lower level functions */
void propagate_t_mosse(mosse_fft *obj, int idx) {
  int ix, id, ik, nx = obj->nx, ndat = obj->ndat, nd = obj->nd[idx],
                  nk = nd - 1, nk2 = nk * nk;
  double *vars = obj->x, *d, *dd = obj->wrk;
  double e, tmp1, tmp2, lambda_x, mu_x, z_x, Q_x, d_x;

  for (ix = 0; ix < ndat; ix++) {
    lambda_x = obj->lambda[ix];
    mu_x = obj->mu[ix];
    z_x = obj->z[ix];
    e = vars[ix];

    /* Update the E values */
    tmp1 = mu_x - lambda_x * e;
    tmp2 = z_x * (e - 1);
    vars[ix] = (tmp1 + tmp2 * mu_x) / (tmp1 + tmp2 * lambda_x);

    tmp1 =
        (lambda_x - mu_x) / (z_x * lambda_x - mu_x + (1 - z_x) * lambda_x * e);
    /* Here is the D scaling factor */
    dd[ix] = z_x * tmp1 * tmp1;
  }

  /* Update the D values */
  for (id = 1; id < nd; id++) {
    d = obj->x + nx * id;

    for (ix = 0; ix < ndat; ix++) {
      if (d[ix] < 0)
        d[ix] = 0;
      else
        d[ix] *= dd[ix];
    }
  }

  for (id = 1; id < nd; id++) {
    d = obj->wrkd + nx * id;

    for (ix = 0; ix < ndat; ix++) {
      d[ix] = 0;

      for (ik = 0; ik < nk; ik++) {
        Q_x = obj->Q[ix * nk2 + nk * (id - 1) + ik];
        d_x = obj->x[nx * (ik + 1) + ix];
        d[ix] += d_x * Q_x;
      }
    }
  }

  for (id = 1; id < nd; id++) {
    for (ix = 0; ix < ndat; ix++)
      obj->x[nx * id + ix] = obj->wrkd[nx * id + ix];
  }
}

void propagate_x_mosse(mosse_fft *obj, int idx) {
  double *x, *dd, *wrk = obj->wrk;
  int i, id, nx = obj->nx;
  int nkl = obj->nkl, nkr = obj->nkr, npad = obj->npad;
  int nd = obj->nd[idx];

  x = obj->x + 0;
  for (i = 0; i < nkl; i++)
    wrk[i] = x[i];
  for (i = nx - npad - nkr; i < nx - npad; i++)
    wrk[i] = x[i];

  convolve(obj->fft[idx], obj->kern_y);

  x = obj->x + 0;
  for (i = 0; i < nkl; i++)
    x[i] = wrk[i];
  for (i = nx - npad - nkr; i < nx - npad; i++)
    x[i] = wrk[i];
  for (i = nx - npad; i < nx; i++)
    x[i] = 0;

  for (id = 1; id < nd; id++) {
    x = obj->x + (obj->nx) * id;
    dd = obj->wrkd + (obj->nx) * id;
    for (i = 0; i < nkl; i++)
      x[i] = dd[i];
    for (i = nx - npad - nkr; i < nx - npad; i++)
      x[i] = dd[i];
    for (i = nx - npad; i < nx; i++)
      x[i] = 0;
  }
}
