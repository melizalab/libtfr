/**
 * @file   mtm_impl.h
 * @author C Daniel Meliza
 * @date   Mon Mar  7 2016
 */

#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/**
 * Multi-taper FFT transformation structure. Contains basic FFT
 * parameters as well as pointers to tapers (i.e. windowing
 * functions), output buffers, and the FFTW plan.
 *
 */
struct mfft_s {
        int nfft;        /**< number of points in the transform */
        int npoints;     /**< number of points in the taper(s) */
        int ntapers;     /**< number of tapers */
        double *tapers;  /**< array holding tapers, dim ntapers x npoints */
        double *weights; /**< array holding taper weights, dim ntapers */
        double *buf;     /**< workspace for FFTW, dim ntapers x npoints */
        fftw_plan plan;  /**< FFTW plan */
};

/**
 * The number of frequencies in the spectrum of a real signal
 */
#define SPEC_NFREQ(mtm) (mtm->nfft/2 + 1)

/**
 * The number of frames in a spectrogram that only includes frames
 * with full support in the signal.
 *
 * @param mtm          transform object
 * @param signal_size  number of points in the signal
 * @param step_size    number of points shifted between frames
 */
#define SPEC_NFRAMES(mtm,signal_size,step_size) ((signal_size - mtm->npoints)/step_size + 1)
