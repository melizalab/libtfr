/*
 * tfr.h
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
 */
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/**
 * Multi-taper FFT transformation structure. Contains basic FFT
 * parameters as well as pointers to tapers (i.e. windowing
 * functions), output buffers, and the FFTW plan
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
 *   *tapers - pointer to npoints*ntapers array of windowing functions
 *   *lambdas - eigenvalues for tapers; if NULL, assign weight of 1.0 to each taper
 *
 * Returns:
 *   pointer to mfft_params structure (owned by caller)
 *
 * Note:
 *   pointers to tapers and lambdas are now owned by the return mtfft structure
 */
mfft* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);

#ifndef NO_LAPACK
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
#endif

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
 * Compute multitaper FFT of a signal. Note that this can be used for
 * single taper FFTs, if the mfft structure has been
 * initialized with a single window
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
 *  spec     - reassigned spectrogram. needs to be allocated and zero-filled before calling
 *
 */  
void mtm_spec(mfft *mtm, double *spec, const double *samples, int nsamples, int shift, int adapt);

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
 *  flock    - frequency locking parameter
 *  tlock    - time locking parameter
 *
 * Outputs:
 *  spec     - reassigned spectrogram. needs to be allocated and zero-filled before calling
 *
 */  
void tfr_spec(mfft *mtm, double *spec, const double *samples, int nsamples, int k, int shift,
	 double flock, int tlock);

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
 *   tapers  - k DPSS sequences in order of decreasing eigenvalue (size npoints*k)
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
 *  nfreq and dt parameters.  
 */
void tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
		  int N, int nfreqout, double dt, double qthresh, double flock, int tminlock, int tmaxlock);

