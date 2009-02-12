#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "mtm.h"

extern void dstegr_( char *JOBZ, char *RANGE, int *N, double *D, double *E, double *VL, double *VU, 
		     int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, 
		     int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK,
		     int *LIWORK, int *INFO );

extern void dsterf_(int *N, double *D, double *E, int *INFO);

extern void dgttrf_(int *N, double *DL, double *D, double *DU, double *DU2, 
		    int *IPIV, int *INFO );

extern void dgttrs_(char *TRANS, int *N, int *NRHS, 
		    double *DL, double *D, double *DU, double *DU2, 
		    int *IPIV, double *B, int *LDB, int *INFO );

extern void dgtsv_(int *N, int *NRHS, 
		   double *DL, double *D, double *DU, double *B, 
		   int *LDB, int *INFO );

#define SINC(A) sin(M_PI * 2.0 * W * (A))/(M_PI * 2.0 * W * (A))

/*
 * Solve a symmetric tridiagonal system of equations; this doesn't exist 
 * in LAPACK.
 *
 * Inputs:
 *   n - matrix order
 *   e - subdiagonal (length n; first element is ignored)
 *   d - diagonal (length n)  (destroyed in computation)
 *   b - right hand side of equation
 *
 * Outputs:
 *   the solution is returned in b
 */
int
tridisolve(int N, double *e, double *d, double *b)
{
	int j;
	double mu;
        
    	for (j = 0; j < N-1; j++) {
		mu = e[j+1]/d[j];
		d[j+1] = d[j+1] - e[j+1]*mu;
		b[j+1]  = b[j+1] -  b[j]*mu;
	}

	if (fabs(d[N-1]) < DBL_EPSILON)
		return -1;
	else {
	        b[N-1] = b[N-1]/d[N-1];

	        for (j=N-2; j >= 0; j--) {
		       b[j] = (b[j] - e[j+1]*b[j+1])/d[j];
	        }
        }
	return 0;
}

/*
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

/*
 * Compute the self-convolution of a vector using FFT. 
 *
 * Inputs:
 *   N - number of points
 *   x - input vector
 *
 *   y - output vector (not allocated; needs to have N points)
 */
void
fftconv(int N, double *x, double *y)
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

/*
 * Computes the discrete prolate spherical sequences used in the
 * multitaper method power spectrum calculations.
 *
 * Inputs:
 *   npoints   the number of points in the window
 *   nw        the time-bandwidth product. Must be an integer or half-integer
 *             (typical choices are 2, 5/2, 3, 7/2, or 4)
 *   k         how many DPSS vectors to return (up to npoints but k>nw*2-1 are not stable)
 *
 * Outputs (not allocated):
 *   tapers  - k DPSS sequences in order of decreasing eigenvalue (size npoints*k)
 *   lambdas - k eigenvalues associated with each taper
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
	double *d, *sd, *dd1, *dd2, *ee;
	double *taper; 

	double W, ff;

	if ((NW < 0) || (k < 1) || (k >= npoints) || (npoints < 0) || (NW >= npoints/2))
		return -1;

	W = NW/npoints;
	
	d = (double*)malloc(npoints*sizeof(double));
	sd = (double*)malloc(npoints*sizeof(double));
	dd1 = (double*)malloc(npoints*sizeof(double));
	dd2 = (double*)malloc(npoints*sizeof(double));
	// tridisolve ignores first element of subdiagonal
	ee = (double*)malloc((npoints+1)*sizeof(double)); 

	for (i = 0; i < npoints; i++) {
		ff = (npoints - 1 - 2*i);
		d[i] = dd1[i] = 0.25 * cos(2*M_PI*W) * ff * ff;
		sd[i] = ee[i+1] = (i+1) * (npoints-(i+1))/2.0;
	}

	// lapack eigenvalue solver; values stored in d in increasing order
	dsterf_(&npoints,d,sd,&rv);
	if (rv != 0) return -2;
	
	// set up tridiagonal equations:
	for (j = 0; j < k; j++) {
		taper = tapers + j * npoints;  // point into tapers array
		lambda[j] = d[npoints-(j+1)];
		// initialize taper
		for (i = 0; i < npoints; i++)
			taper[i] = sin((j+1) * M_PI * i / (npoints-1));

		for (m = 0; m < 3; m++) {
			// diagonal destroyed by tridisolve
			for (i = 0; i < npoints; i++)
				dd2[i] = dd1[i] - lambda[j];
			if (tridisolve(npoints, ee, dd2, taper) < 0)
				return -2;
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
	free(ee);
	return 0;
}

/*
 * Initialize a multitaper mtm transform using preallocated tapers (i.e. with dpss()))
 * 
 * Inputs:
 *   nfft - number of points in the transform
 *   npoints - number of points in the tapers (windows)
 *   ntapers - number of tapers
 *   *tapers - pointer to npoints*ntapers array of windowing functions
 *   *lambdas - eigenvalues for tapers
 */

mtfft_params* 
init_mtm_prealloc(int nfft, int npoints, int ntapers, double* tapers, double *lambdas)
{
	mtfft_params *mtm;
	mtm = (mtfft_params*)malloc(sizeof(mtfft_params));

	mtm->nfft = nfft;
	mtm->npoints = npoints;
	mtm->ntapers = ntapers;
	mtm->tapers = tapers;
	mtm->lambdas = lambdas;

	mtm->buf = (double*)fftw_malloc(nfft*ntapers*sizeof(double));
	//mtm->out_buf = (fftw_complex*)fftw_malloc((nfft/2+1)*ntapers*sizeof(fftw_complex));
	
	// set up fftw plan
	int *n_array, i;
	fftw_r2r_kind *kind;
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

void
destroy_mtm(mtfft_params *mtm)
{
	free(mtm->tapers);
	free(mtm->lambdas);
	fftw_destroy_plan(mtm->plan);
	fftw_free(mtm->buf);
	free(mtm);
}

/*
 * Compute multitaper FFT of a signal.
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
mtfft(mtfft_params *mtm, short *data, int nbins)
{
	// copy data * tapers to buffer
	int size = mtm->npoints;
	int i,j;
	int nt = (nbins < size) ? nbins : size;
	double pow = 0.0;

	for (i = 0; i < mtm->ntapers; i++) {
		for (j = 0; j < nt; j++) {
			mtm->buf[j+i*size] = mtm->tapers[j+i*size] * data[j];
			pow += data[j] * data[j];
		}
	}
	// zero-pad rest of buffer
	for (i = 0; i < mtm->ntapers; i++) {
		for (j = nt; j < mtm->nfft; j++)
			mtm->buf[j+i*size] = 0.0;
	}

	fftw_execute(mtm->plan);

	return pow;
}

/*
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
mtpower(mtfft_params *mtm, double *pow, double sigpow)
{
	int nfft = mtm->nfft;
	int ntapers = mtm->ntapers;
	int real_count = nfft / 2 + 1;
	int imag_count = (nfft+1) / 2 - 1;
	int t,n;

	if (sigpow<=0.0 || ntapers==1) {
		memset(pow, 0, real_count*sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				pow[n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t]/ntapers;
			for (n = 1; n < imag_count; n++)
				pow[n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t]/ntapers;
		}
	}
	else {
		double est, num, den, w;
		double tol, err;
		double *Sk;
		Sk = (double*)calloc(ntapers*real_count,sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				Sk[t*nfft+n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t];
			for (n = 1; n < imag_count; n++)
				Sk[t*nfft+n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t];
		}
		// initial guess is average of first two tapers
		for (n = 0; n < real_count; n++)
			pow[n] = (Sk[n] + Sk[nfft+n])/2;

		tol = 0.0005 * sigpow;
		err = 0;
		while (err > tol) {
			for (n = 0; n < real_count; n++) {
				est = pow[n];
				num = den = 0;
				for (t=0; t < ntapers; t++) {
					w = est / (est * mtm->lambdas[t] + sigpow * (1 -mtm->lambdas[t]));
					w = w * w * mtm->lambdas[t];
					num += w * Sk[t*nfft+n];
					den += w;
				}
				pow[n] = num/den;
				err += fabs(num/den-est);
			}
		}
		free(Sk);
	}
}
