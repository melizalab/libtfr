
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "tfr.h"

int N = 256;
int Np = 201;
double NW = 3.5;
int step = 10;
int k = 6;
double tm = 6.0;

void
fmsin(double *val, int N, double fnormin, double fnormax, double period, double t0,
      double fnorm0, double pm1)
{
	double fnormid, delta, phi;
	int t;

	fnormid = 0.5 * (fnormax+fnormin);
	delta = 0.5 * (fnormax-fnormin);
	phi = - copysign(1.0,pm1) * acos((fnorm0 - fnormid)/delta);

	for (t = 0; t < N; t++) {
		val[t] = creal(cexp(I * 2 * M_PI * fnormid * (t - t0) + 
				   delta * period * (sin(2 * M_PI * (t - t0) / period + phi) - sin(phi))));
	}
}

int 
main(int argc, char **argv)
{
	mfft *mtmh;
	double *sig, *psd, *specgram;
	double sigpow;

	printf("* Testing TFR library:\n");
	printf("* N = %d\n", N);
	printf("* shift = %d\n", step);
	printf("* MTM NW = %3.2f\n", NW);
	printf("* TFR Np = %d\n", Np);
	printf("* TFR k = %d\n", k);
	printf("* TFR tm = %3.2f\n", tm);

	sig = (double*)malloc(17590 * sizeof(double));
	psd = (double*)malloc(N * sizeof(double));
	
	fmsin(sig, 17590, 0.15, 0.45, 1024, 256./4, 0.3, -1);

	printf("* Input signal (8300):\n");
	for (int i = 0; i < N; i++)
		printf("%3.2f ", *(sig+8300+i));
	printf("\n");
	
	printf("* Testing MTM PSD:\n");
	mtmh = mtm_init_dpss(N, NW, (int)(NW*2-1));
	sigpow = mtfft(mtmh, sig+8300, N);
	mtpower(mtmh, psd, sigpow);
	mtm_destroy(mtmh);

	for (int i = 0; i < N; i++)
		printf("%3.2f ", psd[i]);
	printf("\n");

	free(psd);
	free(sig);
}

	
	
