#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

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
        //const double eps = std::numeric_limits<double>::epsilon();
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

void
fftconv(int N, double *x, double *y)
{
	int i;
	double *X;
	fftw_complex *X1, *X2;
	//double *X1, *X2;
	fftw_plan plan;

	X = (double*)calloc(N * 2, sizeof(double));
	//X1 = (double*)fftw_malloc(N * 2 * sizeof(double));
	//X2 = (double*)fftw_malloc(N * 2 * sizeof(double));
	X1 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));
	X2 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));

	memcpy(X, x, sizeof(double)*N);
	printf("\n");
	plan = fftw_plan_dft_r2c_1d(N*2, X, X1, FFTW_ESTIMATE);
	//plan = fftw_plan_r2r_1d(N*2, X, X1, FFTW_R2HC, FFTW_ESTIMATE);
	fftw_execute(plan);

	// flip x and compute again
	for (i = 0; i < N; i++)
		X[i] = x[N-1-i];
	plan = fftw_plan_dft_r2c_1d(N*2, X, X2, FFTW_ESTIMATE);
	//plan = fftw_plan_r2r_1d(N*2, X, X2, FFTW_R2HC, FFTW_ESTIMATE);
	fftw_execute(plan);

	for (i = 0; i < N; i++)
		X1[i] *= X2[i];

	// inverse fft
	plan = fftw_plan_dft_c2r_1d(N*2, X1, X, FFTW_ESTIMATE);
	//plan = fftw_plan_r2r_1d(N*2, X1, X, FFTW_HC2R, FFTW_ESTIMATE);
	fftw_execute(plan);
	
	for (i = 0; i < N; i++)
		y[i] = X[i] / N / 2;
	//memcpy(y,X,N*sizeof(double));
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
 *   k         how many DPSS vectors to return (up to npoints)
 *
 * Outputs:
 *   tapers  - k DPSS sequences in order of decreasing eigenvalue (size npoints*k)
 *   lambdas - k eigenvalues associated with each taper
 *
 *
 * Returns:
 *    0 for success
 *   -1 for invalid parameters (NW >= npoints/2, npoints < 0, k < 0, etc)
 */
int
dpss(double *tapers, double *lambda, int npoints, double NW, int k) 
{
	int i, j, m, rv;
	double *d, *sd, *dd1, *dd2, *ee, *fftw_buf;
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
	
	// set up tridiagonal equations:
	for (j = 0; j < k; j++) {
		taper = tapers + j * npoints;  // point into tapers array
		lambda[j] = d[npoints-(j+1)];
		printf("taper %d:\n", j);
		// initialize taper
		for (i = 0; i < npoints; i++)
			taper[i] = sin((j+1) * M_PI * i / (npoints-1));

		for (m = 0; m < 3; m++) {
			// diagonal destroyed by tridisolve
			for (i = 0; i < npoints; i++)
				dd2[i] = dd1[i] - lambda[j];
			tridisolve(npoints, ee, dd2, taper);
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

/* 		for (i = 0; i < npoints; i++) { */
/* 			printf("%3.3g ", taper[i]); */
/* 		} */
/* 		printf("\n"); */

		// calculate lambdas
		fftconv(npoints, taper, dd2);
/* 		for (i = 0; i < npoints; i++) { */
/* 			printf("%3.4g, ", dd2[i]); */
/* 		} */
/* 		printf("\n"); */

		//ff = dd2[npoints-1];
		ff = 2.0 * W * dd2[npoints-1];  // last point
		for (i = 0; i < npoints-1; i++)
			//ff += dd2[i];
			ff += dd2[i] * 4.0 * W * SINC(npoints-1-i);
			//printf("%3.5g, ", 4.0 * W * sin(M_PI * 2.0 * W * (npoints-i))/(M_PI * 2.0 * W * (npoints-i)));
			//printf("%3.5g, ", 4.0 * W * SINC(npoints-1-i));

		lambda[j] = ff;
		printf("lambda: %3.6f\n", ff);
	}	
}

int
main(int argc, char** argv)
{
	
	int N = 100;
	double NW = 3.5;
	int k = 5;
	double *tapers, *eigenval;
	
	tapers = (double*)malloc(N*k*sizeof(double));
	eigenval = (double*)malloc(N*sizeof(double));

	dpss(tapers, eigenval, N, NW, k);

}
