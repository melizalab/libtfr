#include <fftw3.h>

typedef struct {
	int nfft;
	int npoints;
	int ntapers;
	double *tapers;
	double *lambdas;
	double *buf;
	fftw_plan plan;
} mfft;

/* public functions */
mfft* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);
mfft* mtm_init_dpss(int nfft, double nw, int ntapers);
double mtfft(mfft *mtm, const short *data, int nbins);
void mtm_destroy(mfft *mtm);
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);
void mtpower(const mfft *mtm, double *out, double sigpower);

int hermf(int N, int M, double tm, double *h, double *Dh, double *Th);
mfft* mtm_init_herm(int nfft, int npoints, int order, double tm);
void tfr_displacements(const mfft *mtm, double *q, double *tdispl, double *fdispl);
void tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
		  int N, int nfreqout, double dt, double qthresh, double flock, int tminlock, int tmaxlock);


/* internal functions */
int tridisolve(int N, const double *e, double *d, double *b);
void renormalize(int N, double *x);
void fftconv(int N, const double *x, double *y);

