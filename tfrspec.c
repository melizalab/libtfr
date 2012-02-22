/*
 * mex (MATLAB C interface) file for computing time-frequency
 * reassignment spectrograms
 *
 * Copyright C Daniel Meliza 2010.  Licensed for use under GNU
 * General Public License, Version 2.  See COPYING for details.
 */
#include <mex.h>
#include <math.h>
#include "tfr.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        int m,n,nt;
        int nfft, step, npoints;
        int ntapers = 6;
        int tlock = 5;
        double tm = 6.0;
        double flock=0.01;
        double *spec;
        mfft* mtmh;

        mxClassID data_type;
        double *data_double;

        if (nlhs < 1)
                return;
        if ((nrhs < 4) || (nrhs > 4) || (nlhs > 1))
                mexErrMsgTxt("SP = TFRSPEC(S, N, step, Np, K, tm, flock, tlock)");

        /* parse arguments */
        m = mxGetM(prhs[0]);
        n = mxGetN(prhs[0]);
        if ((m > 1) && (n > 1))
                mexErrMsgTxt("Input signal must be a 1-D time series");
        nt = (m > n) ? m : n;
        /* check data type */
        data_type = mxGetClassID(prhs[0]);
        if (!(data_type==mxDOUBLE_CLASS))
                mexErrMsgTxt("Input signal must be double precision");

        /*data_double = mxGetPr(prhs[0]);*/
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

        /* generate transform object */
        mtmh = mtm_init_herm(nfft, npoints, ntapers, tm);

        /* allocate output array */
        plhs[0] = mxCreateDoubleMatrix(SPEC_NFREQ(mtmh), SPEC_NFRAMES(mtmh,nt,step),mxREAL);
        spec = mxGetPr(plhs[0]);

        /* calculate spectrogram; use the right precision fxn */
        data_double = mxGetPr(prhs[0]);
        tfr_spec(mtmh, spec, data_double, nt, -1, step, flock, tlock, -1, 0);

        mtm_destroy(mtmh);
}
