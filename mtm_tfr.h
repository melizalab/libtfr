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
mtfft_params* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);
mtfft_params* mtm_init_dpss(int nfft, double nw, int ntapers);
double mtfft(mtfft_params *mtm, const short *data, int nbins);
void mtm_destroy(mtfft_params *mtm);
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);
void mtpower(const mtfft_params *mtm, double *out, double sigpower);

int hermf(int N, int M, double tm, double *h, double *Dh, double *Th);
mtfft_params* mtm_init_herm(int nfft, int npoints, int order, double tm);
void tfr_displacements(const mtfft_params *mtm, double *q, double *tdispl, double *fdispl);
void tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
		  int N, int nfreqout, double dt, double qthresh, double flock, int tminlock, int tmaxlock);


/* internal functions */
int tridisolve(int N, const double *e, double *d, double *b);
void renormalize(int N, double *x);
void fftconv(int N, const double *x, double *y);

