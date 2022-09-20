typedef struct {

    int nYLo, nYHi;
    int nXbinsLo, nXbinsHi;
    int nUsefulXbinsLo, nUsefulXbinsHi;
    int nLeftFlankBinsLo, nLeftFlankBinsHi;
    int nRightFlankBinsLo, nRightFlankBinsHi;

    double dxBin;

    int max_nd;

    // fft-related constant members
    int* nDnEs_ptr;
    int nDimensions;

    double* fft_x_input_lo;
    double* fft_x_input_hi;
    fftw_complex* fft_y_input_lo;
    fftw_complex* fft_y_input_hi;

    double* xLo_ptr; // x-ruler in low resolution
    double* xHi_ptr; // x-ruler in high resolution

    rfftw_plan_real **fftLo;
    rfftw_plan_real **fftHi;

    double* zLo;
    double* wrkLo;
    double* wrkdLo;
    double* zHi;
    double* wrkHi;
    double* wrkdHi;
} quasse_constants_stash;

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
    int flags
);
