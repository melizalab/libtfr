/**
 * @file   tfr.h
 * @author C Daniel Meliza
 * @date   Mon Mar  7 2016
 */

/**
 * @mainpage
 * @section intro Introduction
 *
 * The libtfrspec library computes multi-tapered conventional and
 * time-frequency reassignment spectrograms from real-valued time
 * series. Time-frequency reassignment is a technique for sharpening
 * spectrographic representations of complex time-varying
 * processes. See README for more details on the algorithms used.
 *
 * Transforms are started by initializing an @ref mfft structure with
 * the size of the transform and the tapering functions to be
 * used. Transforms using DPSS and Hermite tapers can be initialized
 * with specialized functions:
 * - mtm_init_dpss()
 * - mtm_init_herm()
 *
 * After initializing the transform object, it can be used to perform
 * any number of operations.  In general, the transformed signal is
 * stored in the @ref mfft object and must be extracted.  See the
 * following functions for calculating real or complex spectra:
 * - mtfft()
 * - mtpower()
 * - mtcomplex()
 *
 * Spectrograms are generally calculated by sliding an analysis window
 * along the signal and using it to calculate the spectral density in
 * that window. This is the short-time fourier transform.  The
 * time-frequency reassigned spectrogram improves on this by allowing
 * energy to be reassigned across frames, but the general principle is
 * the same.  The spectrogram functions are fairly straightforward but
 * they do require preallocation of the output array.
 * - mtm_spec()  for standard STFT with one or more tapers
 * - mtm_zspec() for complex STFT with one or more tapers
 * - tfr_spec()  for reassigned spectrogram
 *
 * Other functions in the library are used internally to generate
 * Hermite and DPSS tapers and to do the time-frequency reassignment.
 * They are included in the public interface, though ordinary use
 * should have no need of them.
 *
 * @section Example
 *
 * For an example of how to use the library, see test_tfr.c
 *
 * @section lic License
 * Copyright C Daniel Meliza 2010-2016.  Licensed for use under GNU
 * General Public License, Version 2.  See COPYING for details.
 */
#ifndef _LIBTFR_H
#define _LIBTFR_H

#define LIBTFR_VERSION "2.0.0"

#ifdef __cplusplus
extern "C" {
#include <complex>
#else
#include <complex.h>
#endif

/**
 * Opaque pointer type for multitaper fft transforms
 */
typedef struct mfft_s mfft;

/* initialization and destruction functions */

/**
 * Initialize a multitaper mtm transform and allocate memory for tapers/window
 * functions.
 *
 * @param nfft  number of points in the transform
 * @param npoints  number of points in the tapers (windows)
 * @param ntapers  number of tapers
 * @returns  pointer to mfft_params structure (owned by caller)
 *
 */
mfft * mtm_init(int nfft, int npoints, int ntapers);

/**
 * Initialize a mtfft transform using DPSS tapers
 * (i.e. for a standard multitaper transform)
 *
 * @param nfft   number of points in the transform
 * @param npoints  number of points in the tapers
 * @param nw     time-frequency parameter
 * @param ntapers  number of tapers to keep
 * @returns pointer to mfft structure (owned by caller)
 */
mfft * mtm_init_dpss(int nfft, int npoints, double nw, int ntapers);

/**
 * Initialize mtfft transform for reassigned spectrogram
 * (i.e. using hermitian function tapers)
 *
 * @param nfft  the number of points in the fourier transform
 * @param npoints  the number of points in the window; controls the time-frequency resolution
 * @param order  the maximum order of hermite functions to use. actual # of tapers is 3 times this
 * @param tm     time support for the tapers. If 0 or less, use the default of 6
 * @returns pointer to mfft structure
 */
mfft * mtm_init_herm(int nfft, int npoints, int order, double tm);

/**
 * Copy pre-calculated tapers/window functions (e.g. hanning) into a mtfft
 * transform. Size of arrays must match memory allocated by the transform.
 *
 * @param tapers   pointer to ntapers*npoints array of windowing functions
 * @param weights  weights for tapers; if NULL, assign weight of 1.0 to each taper
 *
 */
void mtm_copy(mfft * mtmh, const double * tapers, const double * weights);

/**
 * Frees up the mfft structure and dependent data.
 *
 * Note that references to the tapers are considered to be owned by
 * the structure, so if they were calculated elsewhere do not attempt
 * to access them after calling this function.
 *
 * @param mtm   the structure to release
 */
void mtm_destroy(mfft * mtm);

/* utility functions */

int mtm_nfft(mfft const * mtm);
int mtm_npoints(mfft const * mtm);
int mtm_ntapers(mfft const * mtm);
int mtm_nreal(mfft const * mtm);
int mtm_nframes(mfft const * mtm, int signal_size, int step_size);
double const * mtm_buffer(mfft const * mtm);
double const * mtm_tapers(mfft const * mtm);


/* transformation functions */

/**
 * Compute multitaper FFT of a real-valued signal.
 *
 * Note that this can be used for single taper FFTs, if the mfft
 * structure has been initialized with a single window.  The result is
 * stored in the mfft buffer in half-complex format with dimension
 * ntapers x nfft.  Use mtpower or mtcomplex to extract the
 * transformed signal.
 *
 * @param mtm  parameters for the transform
 * @param data  input data (double-precision floating points)
 * @param nbins  the number of time points in the signal
 * @returns total power in signal (used in computing adaptive power spectra)
 */
double mtfft(mfft * mtm, double const * data, int nbins);

/* spectrogram functions */

/**
 * Extract power spectrum from multiple taper FFT.
 *
 * The 'high-res' method is simply a weighted average of the estimates for each
 * taper. The 'adaptive' method attempts to fit the contribution from each taper
 * to match the total power in the signal.
 *
 * @param  mtm     mfft structure after running mtfft
 * @param  pow     (output) power spectral density (linear scale) of the signal. Needs to be
 *                 preallocated, with dimensions at least nfft/2 + 1;
 * @param  sigpow  total power in the signal. If zero or less, uses high-res method
 */
void mtpower(mfft const * mtm, double * pow, double sigpow);


/**
 * Extract complex multitaper transform of signal from transform object
 *
 * @param mtm  mfft structure after running mtfft
 * @param out  (output) complex transform of signal. Needs to be preallocated with
 *             dimensions at least ntapers by nfft
 */
void mtcomplex(mfft const * mtm, _Complex double * out);

/**
 *  Compute a multitaper spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the PSD for each
 *  frame in the signal.
 *
 * @param mtm       mfft structure; needs to be initialized with tapers
 * @param samples   input signal
 * @param nsamples  number of points in input buffer
 * @param shift     number of samples to shift in each frame
 * @param adapt     if true, use adaptive averaging between tapers (otherwise 'high-res')
 *
 * @param spec      (output) spectrogram, dimension (nsamples-npoints+1)/shift by nfft/2+1
 *                  needs to be allocated and zero-filled before calling
 */
void mtm_spec(mfft * mtm, double *spec, const double *samples, int nsamples, int shift,
              int adapt);

/**
 *  Compute a multitaper complex spectrogram by stepping through a signal. This
 *  function 'fills' a spectrogram by calculating the complex FFT for each taper
 *  and for each frame in the signal.
 *
 * @param mtm       mfft structure; needs to be initialized with tapers
 * @param samples   input signal
 * @param nsamples  number of points in input buffer
 * @param shift     number of samples to shift in each frame
 *
 * @param spec      (output) spectrogram, dimension (nsamples-npoints+1)/shift
 *                  by (ntapers) by (nfft). Must be allocated and zero-filled.
 */
void mtm_zspec(mfft * mtm, _Complex double *spec, const double *samples, int nsamples,
               int shift);

/**
 *  Compute a time-frequency reassignment spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the displaced PSD for each
 *  frame in the signal.
 *
 * @param mtm       mfft structure; needs to be initialized with hermite tapers
 * @param samples   input signal
 * @param nsamples  number of points in input buffer
 * @param k         which taper to use; -1 for all tapers
 * @param shift     number of samples to shift in each frame
 * @param flock     frequency locking parameter (normalized frequency units)
 * @param tlock     time locking parameter (in frames)
 * @param nfreq     output frequency resolution; if <= 0, defaults to nfft/2+1
 * @param fgrid     output frequency grid; if NULL, defaults to linear scale
 *                  from 0 to 0.5 (normalized freq)
 *
 * @param spec      (output) spectrogram, dimension  (nsamples-npoints+1)/shift by nfft/2+1
 *                  needs to be allocated and zero-filled before calling
 *
 */
void tfr_spec(mfft * mtm, double *spec, const double *samples, int nsamples, int k,
              int shift, double flock, int tlock, int nfreq, const double *fgrid);

/* taper generating functions */

/**
 * Computes discrete prolate spherical sequences.
 * These are used in the multitaper method power spectrum calculations.
 *
 * @param npoints   the number of points in the window
 * @param nw        the time-bandwidth product. Must be an integer or half-integer
 *                  (typical choices are 2, 5/2, 3, 7/2, or 4)
 * @param k         how many DPSS vectors to return (up to npoints but k>nw*2-1 are not stable)

 * @param tapers   (output) k DPSS sequences in order of decreasing eigenvalue (size k*npoints)
 * @param lambda   (output) k eigenvalues associated with each taper
 *
 * @returns 0 for success,
 *          -1 for invalid parameters (NW >= npoints/2, npoints < 0, k < 0, etc),
 *          -2 for failed eigenvalue solver (shouldn't ever happen)
 */
int dpss(double *tapers, double *lambda, int npoints, double nw, int k);

/**
 *  Computes a set of orthogonal Hermite functions.
 *  Used in computing multi-taper reassigned spectrograms
 *
 * @param N  the number of points in the window (must be odd; rounded down)
 * @param M  the maximum order of the set of functions
 * @param tm  half-time support
 *
 * @param h  (output) hermite functions (MxN)
 * @param Dh (output) first derivative of h (MxN)
 * @param Th (output) time multiple of h (MxN)
 *
 * @returns The actual number of points in the tapers
 *
 * From the Time-Frequency Toolkit, P. Flandrin & J. Xiao, 2005
 */
int hermf(int N, int M, double tm, double *h, double *Dh, double *Th);

/* reassignment functions */

/**
 * Compute the power spectrum and the time/frequency displacement.
 *
 * @param mtm  mfft object with computed FFT transforms; assumes that there
 *             are 3x tapers as the order of the multitaper transform (K)
 *
 * @param q       (output) power spectrum (NFFT/2+1 x K)
 * @param tdispl  (output) time displacements (NFFT/2+1 x K)
 * @param fdispl  (output) frequency displacements (NFFT/2+1 x K)
 */
void tfr_displacements(mfft const * mtm, double *q, double *tdispl, double *fdispl);

/**
 *  Assign power from a spectrum to a spectrogram based on time-frequency displacements
 *
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
 *
 * Inputs:
 * @param q        power spectrum (N points)
 * @param tdispl   time displacements (N points)
 * @param fdispl   frequency displacements (N points)
 * @param N        number of points in input spectrums
 * @param nfreq    number of frequency bins in output spectrum
 * @param fgrid    array of output frequency bins (optional; see below)
 * @param dt       spacing between columns of output spectrogram (samples)
 * @param qthresh  frequency bins with q<=qthresh are not assigned (unstable)
 * @param flock    maximum frequency displacement (radians; 0.01-0.02 is a good value; 0 to disable)
 * @param tminlock   maximum negative time displacement (number of FRAMES)
 * @param tmaxlock   maximum positive time displacement (number of FRAMES)
 *
 * @param spec     (output) spectrogram (nfreq by >(tmaxlock+tminlock))
 *                 pre-allocate with zeros
 */
void tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
                  int N, int nfreq, const double *fgrid,
                  double dt, double qthresh, double flock, int tminlock, int tmaxlock);

#ifdef __cplusplus
}
#endif

#endif
