#include <fftw3.h>

typedef struct {
	int nfft;
	int npoints;
	int ntapers;
	double *tapers;
	double *lambdas;
	double *buf;
	//fftw_complex *out_buf;
	fftw_plan plan;
} mtfft_params;

/* public functions */
mtfft_params* init_mtm_prealloc(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);
double mtfft(mtfft_params *mtm, short *data, int nbins);
void destroy_mtm(mtfft_params *mtm);
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);


/* internal functions */
int tridisolve(int N, double *e, double *d, double *b);
void renormalize(int N, double *x);
void fftconv(int N, double *x, double *y);
