/*
  Functions to compute time-frequency reassignment spectrograms.
  Asqsembled from various MATLAB sources, including the time-frequency
  toolkit [1], Xiao and Flandrin's work on multitaper reassignment [2]
  and code from Gardner and Magnasco [3].

  The basic principle is to use reassignment to increase the precision
  of the time-frequency localization, essentially by deconvolving the
  spectrogram with the TF representation of the window, recovering the
  center of mass of the spectrotemporal energy.  Reassigned TFRs
  typically show a 'froth' for noise, and strong narrow lines for
  coherent signals like pure tones, chirps, and so forth.  The use of
  multiple tapers reinforces the coherent signals while averaging out
  the froth, giving a very clean spectrogram with optimal precision
  and resolution properties.

  Implementation notes:

  Gardner & Magnasco calculate reassignment based on a different
  algorithm from Xiao and Flandrin.  The latter involves 3 different
  FFT operations on the signal windowed with the hermitian taper
  [h(t)], its derivative [h'(t)], and its time product [t * h(t)].
  The G&M algorithm only uses two FFTs, on the signal windowed with a
  gassian and its time derivative.  If I understand their methods
  correctly, however, this derivation is based on properties of the
  fourier transform of the gaussian, and isn't appropriate for window
  functions based on the Hermitian tapers.

  Therefore, the algorithm is mostly from [2], though I include time
  and frequency locking parameters from [3], which specify how far
  energy is allowed to be reassigned in the TF plane.  Large
  displacements generally arise from numerical errors, so this helps
  to sharpen the lines somewhat. I also included the time/frequency
  interpolation from [3], which can be used to get higher precision
  (at the expense of less averaging) from smaller analysis windows.

  Some fiddling with parameters is necessary to get the best
  spectrograms from a given sort of signal.  Like the window size in
  an STFT, the taper parameters control the time-frequency resolution.
  However, in the reassignment spectrogram the precision
  (i.e. localization) is not affected by the taper size, so the
  effects of taper size will generally only be seen when two coherent
  signals are close to each other in time or frequency.  Nh controls
  the size of the tapers; one can also adjust tm, the time support of
  the tapers, but depending on the number of tapers used, this
  shouldn't get a whole lot smaller than 5.  Increased values of Nh
  result in improved narrowband resolution (i.e. between pure tones)
  but closely spaced clicks can become smeared.  Decreasing Nh
  increases the resolution between broadband components (i.e. clicks)
  but smears closely spaced narrowband components.  The effect of
  smearing can be ameliorated to some extent by adjusting the
  frequency/time locking parameters.

  The frequency zoom parameter can be used to speed up calculation
  quite a bit [3].  Since calculating the multitaper reassigned
  spectrogram takes 3xNtapers FFT operations, smaller FFTs are
  generally better.  The spectrogram can then be binned at a finer
  resolution during reassignment.  These two sets of parameters should
  generate fairly similar results:

  nfft=512, shift=10, tm=6, Nh=257, zoomf=1, zoomt=1  (default)
  nfft=256, shift=10, tm=6, Nh=257, zoomf=2, zoomt=1

  Increasing the order generally reduces the background 'froth', but
  interference between closely spaced components may increase.

  CDM, 8/2008

  [1] http://tftb.nongnu.org/
  [2] http://perso.ens-lyon.fr/patrick.flandrin/multitfr.html
  [3] PNAS 2006, http://web.mit.edu/tgardner/www/Downloads/Entries/2007/10/22_Blue_bird_day_files/ifdv.m 
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mtm_tfr.h"

#ifndef SQR
#define SQR(a) ( (a) * (a) )
#endif

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
int
hermf(int N, int M, double tm, double *h, double *Dh, double *Th)
{
	
	int i, k;
	double dt, *tt, *g, *P, *Htemp;

	// fix even window sizes
	N -= (N % 2) ? 0 : 1;

	tt = (double*)malloc(N*M*sizeof(double));
	g = (double*)malloc(N*sizeof(double));
	P = (double*)malloc(N*(M+1)*sizeof(double));
	Htemp = (double*)malloc(N*(M+1)*sizeof(double));

	dt = 2 * tm / (N-1);
	for (i = 0; i < N; i++) {
		tt[i] = -tm + dt * i;
		g[i] = exp(-tt[i] * tt[i] / 2);
		P[i] = 1.0;
		P[N+i] = 2 * tt[i];
	}

	for (k = 2; k < M+1; k++) {
		for (i = 0; i < N; i++) {
			P[k*N+i] = 2 * tt[i] * P[(k-1)*N+i] - 2 * (k-1) * P[(k-2)*N+i];
		}
	}

	for (k = 0; k < M+1; k++) {
		for (i = 0; i < N; i++) {
			Htemp[k*N+i] = P[k*N+i] * 
				g[i]/sqrt(sqrt(M_PI) * pow(2, k) * gamma(k+1)) *
				sqrt(dt);
		}
	}

	for (k = 0; k < M; k++) {
		for (i = 0; i < N; i++) {
			Dh[k*N+i] = (tt[i] * Htemp[k*N+i] - sqrt(2*(k+1)) * Htemp[(k+1)*N+i])*dt;
			Th[k*N+i] = Htemp[k*N+i] * (-(N-1)/2 + i);
		}
	}


	memcpy(h, Htemp, N*M*sizeof(double));
	free(tt);
	free(g);
	free(P);
	free(Htemp);

	return N;
}

/**
 *  Initialze the mtm engine to compute FFT transforms using the hermitian
 *  function tapers.
 *
 * Inputs:
 *  nfft - the number of points in the fourier transform
 *  npoints - the number of points in the window; controls the time-frequency resolution
 *  order - the maximum order of hermite functions to use. actual # of tapers is 3 times this
 *  tm    - time support for the tapers. If 0 or less, use the default of 6
 */
mfft*
mtm_init_herm(int nfft, int npoints, int order, double tm)
{

	double *tapers = (double*)malloc(npoints*order*3*sizeof(double));
	
	tm = (tm > 0) ? tm : 6;

	npoints = hermf(npoints, order, tm,
			tapers, tapers + order*npoints, tapers + order*npoints*2);

	mfft* mtm = mtm_init(nfft, npoints, order*3, tapers, 0);

	printf("NFFT= %d; npoints= %d; ntapers=%d\n", mtm->nfft, mtm->npoints, mtm->ntapers);
	return mtm;
}

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
void
tfr_displacements(const mfft *mtm, double *q, double *tdispl, double *fdispl)
{

	int i,j;
	int nfft = mtm->nfft;
	int real_count = nfft / 2 + 1;
	int imag_count = (nfft+1) / 2; // not actually the count but the last index
	int K = mtm->ntapers / 3;
	fftw_complex z1,z2,z3;

	for (j = 0; j < K; j++) {
		for (i = 1; i < imag_count; i++) {
			z1 = mtm->buf[j*nfft+i] + mtm->buf[j*nfft+(nfft-i)] * I;
			z2 = mtm->buf[(j*3+1)*nfft+i] + mtm->buf[(j*3+1)*nfft+(nfft-i)] * I;
			z3 = mtm->buf[(j*3+2)*nfft+i] + mtm->buf[(j*3+2)*nfft+(nfft-i)] * I;

			q[j*real_count+i] = cabs(z1) * cabs(z1);
			fdispl[j*real_count+i] =  cimag(z2 / z1 / (2 * M_PI));
			tdispl[j*real_count+i] = creal(z3 / z1);
		}
 		// DC 
 		q[j*real_count] = SQR(mtm->buf[j*nfft]);
		fdispl[j*real_count] = 0.0;
		tdispl[j*real_count] = mtm->buf[(j*3+2)*nfft] / mtm->buf[j*nfft];
		// nyquist
 		if (imag_count < real_count) {
 			i = real_count-1;
 			q[j*real_count+i] = SQR(mtm->buf[j*nfft+i]);
 			fdispl[j*real_count+i] = 0.0; 
 			tdispl[j*real_count+i] = mtm->buf[(j*3+2)*nfft+i] / mtm->buf[j*nfft+i];
 		}
	}
}

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
void
tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
	     int N, int nfreq, double dt, double qthresh, double flock, int tminlock, int tmaxlock)
{

	int f, that, fhat;
	double fref;
	
        for (f = 0; f < N; f++) {
		//spec[f] += q[f];
		fref = (1.0 * f) / N;
		fhat = (int)round((fref - fdispl[f] * 2.0)*nfreq); // note 2xfdisplace for 1-sided psd
		that = (int)round(tdispl[f] / dt);
		//printf("\n%d: %d,%d (%3.3f)", f, fhat, that, q[f]);
		// check that we're in bounds, within locking distance, and above thresh
		if ((fhat < 0) || (fhat >= nfreq))
			continue;
		if (q[f] <= qthresh)
			continue;
		if ((flock > 0) && (fabs(fdispl[f]) > flock))
			continue;
		if ((that > tmaxlock) || (that < -tminlock))
			continue;
                 // make the reassignment
		//printf("- assigned");
		spec[that*nfreq + fhat] += q[f];
	}
}	     

/**
 *  Compute a time-frequency reassignment spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the displaced PSD for each
 *  frame in the signal.
 *
 * Inputs:
 *  mtm - mfft structure; needs to be initialized with hermite tapers
 *  samples - input signal
 *  nsamples - number of points in input buffer
 *  shift    - number of samples to shift in each frame
 *  flock    - frequency locking parameter
 *  tlock    - time locking parameter
 *
 * Outputs:
 *  spec     - reassigned spectrogram. needs to be allocated and zero-filled before calling
 *
 */  
void
tfr_spec(mfft *mtm, double *spec, const short *samples, int nsamples, int shift,
	 double flock, int tlock)
{
	int t;
	int nbins = nsamples / shift;
	int real_count = mtm->nfft / 2 + 1;
	int K = mtm->ntapers / 3;

	double pow = 0.0;
	for (t = 0; t < nsamples; t++)
		pow += (double)samples[t] * samples[t];
	pow /= nsamples;
	printf("Signal: %d samples, %3.4f RMS power\n", nsamples, pow);

	double *q = (double*)malloc(real_count*K*sizeof(double));
	double *td = (double*)malloc(real_count*K*sizeof(double));
	double *fd = (double*)malloc(real_count*K*sizeof(double));

	for (t = 0; t < nbins; t++) {
		mtfft(mtm, samples+(t*shift), nsamples-(t*shift));
		tfr_displacements(mtm, q, td, fd);
		//memcpy(spec+(t*real_count), td, real_count*sizeof(double));
		//memcpy(tdisp+(t*real_count), td, real_count*sizeof(double));
		tfr_reassign(spec+(t*real_count), q, td, fd,
			     real_count, real_count, shift, 1e-6*pow,
			     flock, (t < tlock) ? t : tlock, (t < nbins-tlock) ? tlock : nbins-tlock);
	}
	free(q);
	free(td);
	free(fd);
}
