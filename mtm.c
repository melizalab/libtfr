#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "mtm_tfr.h"

/*
 *  mtm.c - calculate windowed and multitaper FFT transforms
 *
 *  insert some crap here about copyright, etc
 */

/* some LAPACK prototypes */
extern void dsterf_(int *N, double *D, double *E, int *INFO);

extern void dgtsv_(int *N, int *NRHS, 
		   double *DL, double *D, double *DU, double *B, 
		   int *LDB, int *INFO );

#define SINC(A) sin(M_PI * 2.0 * W * (A))/(M_PI * 2.0 * W * (A))

/**
 * Scale a vector by its L2 norm
 *
 * Inputs:
 *   N - number of points
 *   x - vector; altered in place
 */
void
renormalize(int N, double *x)
{
	int i;
	double norm = 0.0;
	for (i = 0; i < N; i++)
		norm += x[i]*x[i];

	norm = sqrt(norm);
	for (i = 0; i < N; i++)
		x[i] /= norm; 
}

/**
 * Compute the self-convolution of a vector using FFT. 
 *
 * Inputs:
 *   N - number of points
 *   x - input vector
 *
 * Outputs:
 *   y - output vector (not allocated; needs to have N points)
 */
void
fftconv(int N, const double *x, double *y)
{
	int i;
	double *X;
	fftw_complex *X1, *X2;
	fftw_plan plan;

	X = (double*)calloc(N * 2, sizeof(double));
	X1 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));
	X2 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));

	memcpy(X, x, sizeof(double)*N);
	plan = fftw_plan_dft_r2c_1d(N*2, X, X1, FFTW_ESTIMATE);
	fftw_execute(plan);

	// flip x and compute again
	for (i = 0; i < N; i++)
		X[i] = x[N-1-i];
	plan = fftw_plan_dft_r2c_1d(N*2, X, X2, FFTW_ESTIMATE);
	fftw_execute(plan);

	for (i = 0; i < N; i++)
		X1[i] *= X2[i];

	// inverse fft
	plan = fftw_plan_dft_c2r_1d(N*2, X1, X, FFTW_ESTIMATE);
	fftw_execute(plan);
	
	for (i = 0; i < N; i++)
		y[i] = X[i] / N / 2;
	fftw_free(X1);
	fftw_free(X2);
	free(X);
	fftw_destroy_plan(plan);
}

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
int
dpss(double *tapers, double *lambda, int npoints, double NW, int k) 
{
	int i, j, m, rv;
	double *d, *sd, *dd1, *dd2, *ee1, *ee2;
	double *taper; 

	double W, ff;

	if ((NW < 0) || (k < 1) || (k >= npoints) || (npoints < 0) || (NW >= npoints/2))
		return -1;

	W = NW/npoints;
	
	d = (double*)malloc(npoints*sizeof(double));
	sd = (double*)malloc(npoints*sizeof(double));
	dd1 = (double*)malloc(npoints*sizeof(double));
	dd2 = (double*)malloc(npoints*sizeof(double));
	ee1 = (double*)malloc((npoints)*sizeof(double)); 
	ee2 = (double*)malloc((npoints)*sizeof(double)); 

	for (i = 0; i < npoints; i++) {
		ff = (npoints - 1 - 2*i);
		d[i] = dd1[i] = 0.25 * cos(2*M_PI*W) * ff * ff;
		sd[i] = ee1[i] = (i+1) * (npoints-(i+1))/2.0;
	}

	// lapack eigenvalue solver; values stored in d in increasing order
	dsterf_(&npoints,d,ee1,&rv);
	if (rv != 0) return -2;
	
	// set up tridiagonal equations:
	for (j = 0; j < k; j++) {
		taper = tapers + j * npoints;  // point into tapers array
		lambda[j] = d[npoints-(j+1)];
		// initialize taper
		for (i = 0; i < npoints; i++)
			taper[i] = sin((j+1) * M_PI * i / (npoints-1));

		for (m = 0; m < 3; m++) {
			// all inputs destroyed by dgtsv
			for (i = 0; i < npoints; i++) {
				dd2[i] = dd1[i] - lambda[j];
				ee1[i] = ee2[i] = sd[i];
			}
			i = 1;
			dgtsv_(&npoints, &i, ee1, dd2, ee2, taper, &npoints, &rv);
			if (rv != 0) return -2;
			renormalize(npoints, taper);
		}

		// fix sign of taper
		if ((j+1) % 2==1) {
			// calculate sum
			ff = 0.0;
			for (i = 0; i < npoints; i++)
				ff += taper[i];
			if (ff < 0.0) {
				for (i = 0; i < npoints; i++)
					taper[i] *= -1;
			}
		}
		else if (taper[2] < 0.0) {
			for (i = 0; i < npoints; i++)
				taper[i] *= -1;
		}

		// calculate lambdas
		fftconv(npoints, taper, dd2);

		ff = 2.0 * W * dd2[npoints-1];  // last point
		for (i = 0; i < npoints-1; i++)
			ff += dd2[i] * 4.0 * W * SINC(npoints-1-i);

		lambda[j] = ff;

	}	
	free(d);
	free(sd);
	free(dd1);
	free(dd2);
	free(ee1);
	free(ee2);
	return 0;
}


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

mtfft_params* 
mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas)
{
	mtfft_params *mtm;
	int *n_array, i;
	fftw_r2r_kind *kind;
	mtm = (mtfft_params*)malloc(sizeof(mtfft_params));

	mtm->nfft = nfft;
	mtm->npoints = npoints;
	mtm->ntapers = ntapers;
	mtm->tapers = tapers;
	if (lambdas)
		mtm->lambdas = lambdas;
	else {
		mtm->lambdas = (double*)malloc(ntapers*sizeof(double));
		for (i = 0; i < ntapers; i++) mtm->lambdas[i] = 1.0;
	}
			

	mtm->buf = (double*)fftw_malloc(nfft*ntapers*sizeof(double));
	//mtm->out_buf = (fftw_complex*)fftw_malloc((nfft/2+1)*ntapers*sizeof(fftw_complex));
	
	// set up fftw plan
	n_array = malloc(sizeof(int)*ntapers);
	kind = malloc(sizeof(int)*ntapers);
	for (i = 0; i < ntapers; i++) {
		n_array[i] = nfft;
		kind[i] = FFTW_R2HC;
	}

	mtm->plan = fftw_plan_many_r2r(1, n_array, ntapers,
				       mtm->buf, NULL, 1, nfft,
				       mtm->buf, NULL, 1, nfft,
				       kind, FFTW_MEASURE);

	free(n_array);
	free(kind);
	return mtm;
}

/**
 * Initialize a mtfft transform and calculate DPSS tapers
 *
 * Inputs:
 *   nfft - number of points in the transform/dpss tapers
 *   nw   - time-frequency parameter
 *   ntapers - number of tapers to keep
 *
 * Returns:
 *   pointer to mfft_params structure (owned by caller)
 */
mtfft_params*
mtm_init_dpss(int nfft, double nw, int ntapers)
{
	double *tapers, *lambdas;
	tapers = (double*)malloc(nfft*ntapers*sizeof(double));
	lambdas = (double*)malloc(nfft*sizeof(double));
	dpss(tapers, lambdas, nfft, nw, ntapers);
	printf("Generated %d tapers of length %d\n", ntapers, nfft);
	int i;
	for (i = 0; i < ntapers; i++)
		printf("%3.6f, ", lambdas[i]);
	printf("\n");
	return mtm_init(nfft, nfft, ntapers, tapers, lambdas);
}

/**
 * Frees up the mtftt_params structure and dependent data. The
 * references to the tapers are considered to be owned by the
 * structure.
 */
void
mtm_destroy(mtfft_params *mtm)
{
	if (mtm->plan) fftw_destroy_plan(mtm->plan);
	if (mtm->tapers) free(mtm->tapers);
	if (mtm->lambdas) free(mtm->lambdas);
	if (mtm->buf) fftw_free(mtm->buf);
	free(mtm);
}

/**
 * Compute multitaper FFT of a signal. Note that this can be used for
 * single taper FFTs, if the mtfft_params structure has been
 * initialized with a single window
 *
 * Inputs:
 *    mtm - parameters for the transform
 *    data - input data (short integers)
 *    nbins - the number of time points in the signal
 *
 * Returns:
 *    total power in signal (used in computing adaptive power spectra)
 */
double
mtfft(mtfft_params *mtm, const short *data, int nbins)
{
	// copy data * tapers to buffer
	int nfft = mtm->nfft;
	int size = mtm->npoints;
	int i,j;
	int nt = (nbins < size) ? nbins : size;
	double pow = 0.0;

	//printf("Windowing data (%d points, %d tapers)\n", nt, mtm->ntapers);
	for (i = 0; i < mtm->ntapers; i++) {
		for (j = 0; j < nt; j++) {
			mtm->buf[j+i*nfft] = mtm->tapers[j+i*size] * (double)data[j];
			pow += data[j] * data[j];
		}
	}
	pow /= mtm->ntapers;
	// zero-pad rest of buffer
	//printf("Zero-pad buffer with %d points\n", mtm->nfft - nt);
	for (i = 0; i < mtm->ntapers; i++) {
		for (j = nt; j < mtm->nfft; j++)
			mtm->buf[j+i*nfft] = 0.0;
	}

	fftw_execute(mtm->plan);

	return pow / nt;
}

/**
 * Compute power spectrum from multiple taper spectrogram.  The
 * 'high-res' method is simply the average of the estimates for each
 * taper weighted by the eigenvalue of the taper.  The 'adaptive'
 * method attempts to fit the contribution from each taper to match
 * the total power in the signal.
 *
 * Inputs:
 *   mtm - mtfft_params structure after running mtfft
 *   sigpow - total power in the signal. If zero or less, uses high-res method
 *
 * Outputs:
 *   pow - power spectral density (linear scale) of the signal. Needs to be
 *         preallocated, with dimensions at least nfft/2 + 1;
 */
void
mtpower(const mtfft_params *mtm, double *pow, double sigpow)
{
	int nfft = mtm->nfft;
	int ntapers = mtm->ntapers;
	int real_count = nfft / 2 + 1;
	int imag_count = (nfft+1) / 2;  // not actually the count but the last index
	int t,n;

	if (sigpow<=0.0 || ntapers==1) {
		memset(pow, 0, real_count*sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				pow[n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t]/ntapers;
			for (n = 1; n < imag_count; n++) {
				pow[n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t]/ntapers;
			}
		}
	}
	else {
		double est, num, den, w;
		double tol, err;
		double *Sk;
		Sk = (double*)calloc(ntapers*real_count, sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				Sk[t*real_count+n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t];
			for (n = 1; n < imag_count; n++)
				Sk[t*real_count+n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t];
			//Sk[t*nfft+n] *= 2;
		}
		// initial guess is average of first two tapers
		err = 0;
		for (n = 0; n < real_count; n++) {
			pow[n] = (Sk[n] + Sk[real_count+n])/2;
			err += abs(pow[n]);
		}

		tol = 0.0005 * sigpow / nfft;
		err /= nfft;
		//printf("err: %3.4g; tol: %3.4g\n", err, tol);
		//for(t = 0; t < ntapers; t++)
		//	printf("%3.4g ", sigpow * (1 - mtm->lambdas[t]));
		//printf("\n");
		while (err > tol) {
			err = 0;
			for (n = 0; n < real_count; n++) {
				est = pow[n];
				num = den = 0;
				//printf("%d: est=%3.4g; ", n, est);
				for (t=0; t < ntapers; t++) {
					w = est / (est * mtm->lambdas[t] + sigpow * (1 - mtm->lambdas[t]));
					w = w * w * mtm->lambdas[t];
					//printf("%3.4g ",Sk[t*real_count+n]);
					num += w * Sk[t*real_count+n];
					den += w;
				}
				pow[n] = num/den;
				err += fabs(num/den-est);
			}
			//printf("err: %3.4g\n", err);
		}
		free(Sk);
	}
	// adjust power for one-sided spectrum
	for (n = 1; n < imag_count; n++)
		pow[n] *= 2;
}

void
getbuffer(const mtfft_params *mtm, double *buf)
{
	memcpy(buf, mtm->buf, mtm->nfft*mtm->ntapers*sizeof(double));
}

void
gettapers(const mtfft_params *mtm, double *buf)
{
	memcpy(buf, mtm->tapers, mtm->npoints*mtm->ntapers*sizeof(double));
}
