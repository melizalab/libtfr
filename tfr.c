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

/**
 *  Computes a set of orthogonal Hermite functions for use in
 *  computing multi-taper reassigned spectrograms
 *
 * Inputs:
 * N - the number of points in the window (must be odd)
 * M - the maximum order of the set of functions (default 6)
 * tm - half-time support (default 6)
 *1
 * Outputs:
 * h - hermite functions (MxN)
 * Dh - first derivative of h (MxN)
 * tt - time support of functions (N)
 *
 * From the Time-Frequency Toolkit, P. Flandrin & J. Xiao, 2005
 */
void
hermf(int N, int M, double tm, double *h, double *Dh, double *tt)
{
	
	int i, k;
	double dt, *g, *P, *Htemp;

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
		}
	}

	memcpy(h, Htemp, N*M*sizeof(double));
	free(g);
	free(P);
	free(Htemp);
}
