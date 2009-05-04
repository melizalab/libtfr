#include <mex.h>
#include "tfr.h"

/*
 * mex (MATLAB C interface) file for computing time-frequency reassignment spectrograms
 */

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m,n;
	int nfft, step, npoints;
	int ntapers = 6;
	int tlock = 5;
	double tm = 6.0;
	double flock=0.01;
	mfft* mtmh;

	if (nlhs < 1)
		return;
	if ((nrhs < 4) || (nrhs > 4) || (nlhs > 1))
		mexErrMsgTxt("SP = TFRSPEC(S, N, step, Np, K, tm, flock, tlock)");
		
	// parse arguments
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	if ((m > 1) && (n > 1))
		mexErrMsgTxt("Input signal must be a 1-D time series");
	// check data type

	nfft = mxGetScalar(prhs[1]);
	step = mxGetScalar(prhs[2]);
	npoints = mxGetScalar(prhs[3]);
	if (nrhs > 4)
		ntapers = mxGetScalar(plhs[4]);
	if (nrhs > 5)
		tm = mxGetScalar(prhs[5]);
	if (nrhs > 6)
		flock = mxGetScalar(prhs[6]);
	if (nrhs > 7)
		tlock = mxGetScalar(prhs[7]);

	mtmh = mtm_init_herm(nfft, npoints, ntapers, tm);
	mtm_destroy(mtmh);
}
