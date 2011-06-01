/*
 * @file   tfr.h
 * @author C Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:35:27 2010
 *
 * The libtfrspec library computes multi-tapered time-frequency
 * reassignment spectrograms from real-valued time
 * series. Time-frequency reassignment is a technique for sharpening
 * spectrographic representations of complex time-varying
 * processes. See README for more details on the algorithms used.
 *
 * Transformations are initialized by calling one of the mtm_init
 * functions, which returns an mfft pointer.  A single multi-taper FFT
 * can be calculated by passing the mfft object to mtfft.
 * Spectrograms are calculated by stepping through a long signal in
 * discrete steps and using the transforms of each step to caculate
 * the power and time-frequency displacements.  All of these steps can
 * be combined by using the tfr_spec functions.
 *
 * Copyright C Daniel Meliza 2010.  Licensed for use under GNU
 * General Public License, Version 2.  See COPYING for details.
 */
#ifndef _LIBTFR_H
#define _LIBTFR_H

#define LIBTFR_VERSION "1.0.1"

#ifdef __cplusplus
extern "C" {
#include <complex>
#else
#include <complex.h>
#endif

#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/**
 * Multi-taper FFT transformation structure. Contains basic FFT
 * parameters as well as pointers to tapers (i.e. windowing
 * functions), output buffers, and the FFTW plan.
 *
 * nfft - number of points in the transformm
 * npoints - size of the window/tapers
 * ntapers - number of tapers
 * tapers - array holding tapers, dimension ntapers x npoints
 * lamdas - weights for tapers dimension ntapers
 * buf    - workspace for FFTW, dimension ntapers x npoints
 * plan   - FFTW plan
 */
typedef struct {
        int nfft;
        int npoints;
        int ntapers;
        double *tapers;
        double *lambdas;
        double *buf;
        fftw_plan plan;
} mfft;

/* initialization and destruction functions */

/**
 * Initialize a multitaper mtm transform using preallocated tapers (i.e. with dpss()))
 *
 * Inputs:
 *   nfft - number of points in the transform
 *   npoints - number of points in the tapers (windows)
 *   ntapers - number of tapers
 *   *tapers - pointer to ntapers*npoints array of windowing functions
 *   *lambdas - eigenvalues for tapers; if NULL, assign weight of 1.0 to each taper
 *
 * Returns:
 *   pointer to mfft_params structure (owned by caller)
 *
 * Note:
 *   pointers to tapers and lambdas are now owned by the returned mfft structure
 */
mfft* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);

/**
 * Initialize a mtfft transform using DPSS tapers (i.e. for a standard
 * multitaper transform)
 *
 * Inputs:
 *   nfft - number of points in the transform/dpss tapers
 *   nw   - time-frequency parameter
 *   ntapers - number of tapers to keep
 *
 * Returns:
 *   pointer to mfft structure (owned by caller)
 */
mfft* mtm_init_dpss(int nfft, double nw, int ntapers);

/**
 *  Initialze the mtm engine to compute FFT transforms using the hermitian
 *  function tapers.
 *
 * Inputs:
 *  nfft - the number of points in the fourier transform
 *  npoints - the number of points in the window; controls the time-frequency resolution
 *  order - the maximum order of hermite functions to use. actual # of tapers is 3 times this
 *  tm    - time support for the tapers. If 0 or less, use the default of 6
 *
 * Returns:
 *   pointer to mfft structure
 */
mfft* mtm_init_herm(int nfft, int npoints, int order, double tm);


/**
 * Frees up the mtftt_params structure and dependent data. Note that
 * references to the tapers are considered to be owned by the
 * structure, so if they were calculated elsewhere do not attempt to
 * access them after calling this function.
 */
void mtm_destroy(mfft *mtm);

/* transformation functions */

/**
 * Compute multitaper FFT of a real-valued signal. Note that this can
 * be used for single taper FFTs, if the mfft structure has been
 * initialized with a single window.  The result is stored in the mfft
 * buffer in half-complex format with dimension ntapers x nfft.  Use
 * mtpower or mtcomplex to extract the transformed signal.
 *
 * Inputs:
 *    mtm - parameters for the transform
 *    data - input data (double-precision floating points)
 *    nbins - the number of time points in the signal
 *
 * Returns:
 *    total power in signal (used in computing adaptive power spectra)
 */
double mtfft(mfft *mtm, const double *data, int nbins);

/* spectrogram functions */

/**
 * Compute power spectrum from multiple taper spectrogram.  The
 * 'high-res' method is simply the average of the estimates for each
 * taper weighted by the eigenvalue of the taper.  The 'adaptive'
 * method attempts to fit the contribution from each taper to match
 * the total power in the signal.
 *
 * Inputs:
 *   mtm - mfft structure after running mtfft
 *   sigpow - total power in the signal. If zero or less, uses high-res method
 *
 * Outputs:
 *   pow - power spectral density (linear scale) of the signal. Needs to be
 *         preallocated, with dimensions at least nfft/2 + 1;
 */
void mtpower(const mfft *mtm, double *pow, double sigpower);


/**
 * Export complex multitaper transform of signal.
 *
 * Inputs:
 *   mtm - mfft structure after running mtfft
 *
 * Outputs:
 *   out - complex transform of signal.  Needs to be preallocated with
 *         dimensions at least ntapers by nfft
 */
void mtcomplex(const mfft *mtm, _Complex double *out);

/**
 *  Compute a multitaper spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the PSD for each
 *  frame in the signal.
 *
 * Inputs:
 *  mtm - mfft structure; needs to be initialized with hermite tapers
 *  samples - input signal
 *  nsamples - number of points in input buffer
 *  shift    - number of samples to shift in each frame
 *  adapt    - if true, use adaptive averaging between tapers (otherwise 'high-res')
 *
 * Outputs:
 *  spec     - output spectrogram, with dimension  (nsamples/shift) by (nfft/2+1)
 *             needs to be allocated and zero-filled before calling
 *
 *
 */
void mtm_spec(mfft *mtm, double *spec, const double *samples, int nsamples, int shift, int adapt);

/**
 *  Compute a multitaper complex spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the complex FFT for each
 *  frame in the signal.  If the mfft object is configured for multiple tapers,
 *  these are not averaged.
 *
 * Inputs:
 *  mtm - mfft structure; needs to be initialized with hermite tapers
 *  samples - input signal
 *  nsamples - number of points in input buffer
 *  shift    - number of samples to shift in each frame
 *
 * Outputs:
 *  spec     - output spectrogram, with dimension  (nsamples/shift) by (ntapers) by (nfft)
 *
 *
 */
void mtm_zspec(mfft *mtm, _Complex double *spec, const double *samples, int nsamples, int shift);

/**
 *  Compute a time-frequency reassignment spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the displaced PSD for each
 *  frame in the signal.
 *
 * Inputs:
 *  mtm - mfft structure; needs to be initialized with hermite tapers
 *  samples - input signal
 *  nsamples - number of points in input buffer
 *  k        - which taper to use; -1 for all tapers
 *  shift    - number of samples to shift in each frame
 *  flock    - frequency locking parameter (normalized frequency units)
 *  tlock    - time locking parameter (in frames)
 *  nfreq    - output frequency resolution; if <= 0, defaults to nfft/2+1
 *  fgrid    - output frequency grid; if NULL, defaults to linear scale from 0 to 0.5 (normalized freq)
 *
 * Outputs:
 *  spec     - output spectrogram, with dimension  (nsamples/shift) by (nfft/2+1)
 *             needs to be allocated and zero-filled before calling
 *
 */
void tfr_spec(mfft *mtm, double *spec, const double *samples, int nsamples, int k, int shift,
              double flock, int tlock, int nfreq, const double *fgrid);

/* taper generating functions */

/**
 * Computes the discrete prolate spherical sequences used in the
 * multitaper method power spectrum calculations.
 *
 * Inputs:
 *   npoints   the number of points in the window
 *   nw        the time-bandwidth product. Must be an integer or half-integer
 *             (typical choices are 2, 5/2, 3, 7/2, or 4)
 *   k         how many DPSS vectors to return (up to npoints but k>nw*2-1 are not stable)
 *
 * Outputs:
 *   tapers  - k DPSS sequences in order of decreasing eigenvalue (size k*npoints)
 *   lambdas - k eigenvalues associated with each taper
 *   [outputs need to be allocated]
 *
 * Returns:
 *    0 for success
 *   -1 for invalid parameters (NW >= npoints/2, npoints < 0, k < 0, etc)
 *   -2 for failed eigenvalue solver (shouldn't ever happen)
 */
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);

/**
 *  Computes a set of orthogonal Hermite functions for use in
 *  computing multi-taper reassigned spectrograms
 *
 * Inputs:
 * N - the number of points in the window (must be odd; rounded down)
 * M - the maximum order of the set of functions
 * tm - half-time support
 *
 * Outputs:
 * h - hermite functions (MxN)
 * Dh - first derivative of h (MxN)
 * Th - time multiple of h (MxN)
 *
 * Returns:
 *  The actual number of points in the tapers
 *
 * From the Time-Frequency Toolkit, P. Flandrin & J. Xiao, 2005
 */
int hermf(int N, int M, double tm, double *h, double *Dh, double *Th);

/* reassignment functions */

/**
 * Compute the power spectrum and the time/frequency displacement.
 *
 * Inputs:
 *   mtm - mfft object with computed FFT transforms; assumes that there
 *         are 3x tapers as the order of the multitaper transform (K)
 *
 * Outputs:
 *   q   - power spectrum (NFFT/2+1 x K)
 *   tdispl - time displacements (NFFT/2+1 x K)
 *   fdispl - frequency displacements (NFFT/2+1 x K)
 */
void tfr_displacements(const mfft *mtm, double *q, double *tdispl, double *fdispl);

/**
 *  Assign power from a spectrum to a spectrogram based on time-frequency displacements
 *
 * Inputs:
 *  q - power spectrum (N points)
 *  tdispl  - time displacements (N points)
 *  fdispl  - frequency displacements (N points)
 *  N       - number of points in input spectrums
 *  nfreq   - number of frequency bins in output spectrum
 *  fgrid   - array of output frequency bins (optional; see below)
 *  dt      - spacing between columns of output spectrogram (samples)
 *  qthresh - frequency bins with q<=qthresh are not assigned (unstable)
 *  flock   - maximum frequency displacement (radians; 0.01-0.02 is a good value; 0 to disable)
 *  tminlock  - maximum negative time displacement (number of FRAMES)
 *  tmaxlock  - maximum positive time displacement (number of FRAMES)
 *
 * Outputs:
 *  spec    - output spectrogram (nfreq by >(tmaxlock+tminlock))
 *
 * Note:
 *  The time-frequency reassignment spectrogram is built up through
 *  calls to this function for each time frame.  The spectrum in q
 *  contributes to a range of time bins in spec which is limited by
 *  the tminlock and tmaxlock parameters.  This in turn controls how
 *  the memory pointed to by *spec is accessed.  At the edges of the
 *  spectrogram tminlock and tmaxlock need to be adjusted to avoid
 *  accessing invalid memory locations.  Note that the units are frames,
 *  to make allocating the memory a bit easier.
 *
 *  The bin resolution of the output spectrogram is controlled by the
 *  nfreq and dt parameters.  A spectrogram with arbitrary frequency
 *  bins (e.g. logarithmic) can be generated by specifying an array
 *  fgrid[nfreq], which must contain positive, monotonically
 *  increasing frequency values; energy is assigned to the nearest
 *  value in the grid (i.e. the grid specifies center frequencies)
 */
void tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
                  int N, int nfreq, const double *fgrid,
                  double dt, double qthresh, double flock, int tminlock, int tmaxlock);

#ifdef __cplusplus
}
#endif

#endif
