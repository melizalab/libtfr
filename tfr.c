/*
 * tfr.c - functions for calculating time-frequency reassignment spectrograms
 *
 * Hermite tapers and reassignment algorithms adapted from MATLAB code by Xiao and Flandrin,
 * http://perso.ens-lyon.fr/patrick.flandrin/multitfr.html
 * Frequency locking algorithm from Gardner and Magnasco
 * (http://web.mit.edu/tgardner/www/Downloads/Entries/2007/10/22_Blue_bird_day_files/ifdv.m)
 *
 * All other code Copyright C Daniel Meliza 2010.  Licensed for use
 * under GNU General Public License, Version 2.  See COPYING for
 * details.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "tfr.h"
#include "mtm_impl.h"

#ifndef SQR
#define SQR(a) ( (a) * (a) )
#endif

int
hermf(int N, int M, double tm, double *h, double *Dh, double *Th)
{

        int i, k;
        double dt, *tt, *g, *P, *Htemp;

        // fix even window sizes
        N -= (N % 2) ? 0 : 1;

        tt = (double*)malloc(N*M*sizeof(double));
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
                                g[i]/sqrt(sqrt(M_PI) * pow(2, k) * tgamma(k+1)) *
                                sqrt(dt);
                }
        }

        for (k = 0; k < M; k++) {
                for (i = 0; i < N; i++) {
                        Dh[k*N+i] = (tt[i] * Htemp[k*N+i] - sqrt(2*(k+1)) * Htemp[(k+1)*N+i])*dt;
                        Th[k*N+i] = Htemp[k*N+i] * (-(N-1)/2 + i);
                }
        }


        memcpy(h, Htemp, N*M*sizeof(double));
        free(tt);
        free(g);
        free(P);
        free(Htemp);

        return N;
}

mfft *
mtm_init_herm(int nfft, int npoints, int order, double tm)
{
        mfft * mtm = mtm_init(nfft, npoints, order*3);
        tm = (tm > 0) ? tm : 6;
        npoints = hermf(npoints, order, tm,
                        mtm->tapers, mtm->tapers + order*npoints, mtm->tapers + order*npoints*2);

        return mtm;
}

/**
 * Determine closest bin in a grid of arbitrarily spaced
 * frequencies. fgrid is assumed to be monotonically
 * increasing. Transform values to do comparisons on a different scale
 * (i.e. log).  Returns -1 if the test value is out of range.
 *
 */
static int
find_bin(double f, const double *fgrid, int nfreq) {
        double diff;
        int i;
        if (f < *fgrid || f > fgrid[nfreq-1]) return -1;
        for (i = 1; i < nfreq; i++) {
                diff = fgrid[i] - f;
                if (diff >= 0)
                        return ((f-fgrid[i-1]) < diff) ? i-1 : i;
        }
        return -1;
}

void
tfr_displacements(mfft const * mtm, double *q, double *tdispl, double *fdispl)
{

        int i,j;
        int nfft = mtm->nfft;
        int real_count = nfft / 2 + 1;
        int imag_count = (nfft+1) / 2; // not actually the count but the last index
        int K = mtm->ntapers / 3;
        fftw_complex z1,z2,z3;

        for (j = 0; j < K; j++) {
                for (i = 1; i < imag_count; i++) {
                        z1 = mtm->buf[j*nfft+i] + mtm->buf[j*nfft+(nfft-i)] * I;
                        z2 = mtm->buf[(K+j)*nfft+i] + mtm->buf[(K+j)*nfft+(nfft-i)] * I;
                        z3 = mtm->buf[(2*K+j)*nfft+i] + mtm->buf[(2*K+j)*nfft+(nfft-i)] * I;

                        q[j*real_count+i] = cabs(z1) * cabs(z1);
                        fdispl[j*real_count+i] =  cimag(z2 / z1 / (2 * M_PI));
                        tdispl[j*real_count+i] = creal(z3 / z1);
                }
                // DC
                q[j*real_count] = SQR(mtm->buf[j*nfft]);
                fdispl[j*real_count] = 0.0;
                tdispl[j*real_count] = mtm->buf[(2*K+j)*nfft] / mtm->buf[j*nfft];
                // nyquist
                if (imag_count < real_count) {
                        i = real_count-1;
                        q[j*real_count+i] = SQR(mtm->buf[j*nfft+i]);
                        fdispl[j*real_count+i] = 0.0;
                        tdispl[j*real_count+i] = mtm->buf[(2*K+j)*nfft+i] / mtm->buf[j*nfft+i];
                }
        }
}

void
tfr_reassign(double *spec, const double *q, const double *tdispl, const double *fdispl,
             int N, int nfreq, const double *fgrid,
             double dt, double qthresh, double flock, int tminlock, int tmaxlock)
{

        int f, that, fhat;
        double fref;

        for (f = 0; f < N; f++) {
                //spec[f] += q[f];
                fref = (0.5 * f) / N;
                if (fgrid==0) {
                        fhat = (int)round((fref - fdispl[f])*2.0*nfreq);
                        if ((fhat < 0) || (fhat >= nfreq))
                                continue;
                }
                else {
                        fhat = find_bin(fref - fdispl[f], fgrid, nfreq);
                        if (fhat <  0)
                                continue;
                }
                that = (int)round(tdispl[f] / dt);
                //printf("\n%d: %d,%d (%3.3f)", f, fhat, that, q[f]);
                // check that we're in bounds, within locking distance, and above thresh
                if (q[f] <= qthresh)
                        continue;
                if ((flock > 0) && (fabs(fdispl[f]) > flock))
                        continue;
                if ((that > tmaxlock) || (that < -tminlock))
                        continue;
                 // make the reassignment
                //printf("- assigned");
                spec[that*nfreq + fhat] += q[f];
        }
}

void
tfr_spec(mfft * mtm, double *spec, const double *samples, int nsamples, int k, int shift,
         double flock, int tlock, int nfreq, const double *fgrid)
{
        int t,mink = 0;
        int nbins = SPEC_NFRAMES(mtm, nsamples, shift);
        int real_count = SPEC_NFREQ(mtm);
        int K = mtm->ntapers / 3;
        if (nfreq <= 0) nfreq = real_count;

        double pow = 0.0;
        for (t = 0; t < nsamples; t++)
                pow += (double)samples[t] * samples[t];
        pow /= nsamples;
        //printf("Signal: %d samples, %3.4f RMS power\n", nsamples, pow);

        double *q = (double*)malloc(real_count*K*sizeof(double));
        double *td = (double*)malloc(real_count*K*sizeof(double));
        double *fd = (double*)malloc(real_count*K*sizeof(double));

        if (k >= 0) {
                mink = k;
                K = k+1;
        }
        //for (k = mink; k < K; k++) printf("Calculating spectrogram for taper %d\n", k);
        for (t = 0; t < nbins; t++) {
                mtfft(mtm, samples+(t*shift), nsamples-(t*shift));
                tfr_displacements(mtm, q, td, fd);
                for (k = mink; k < K; k++) {
                        tfr_reassign(spec+(t*nfreq),
                                     q+(k*real_count), td+(k*real_count), fd+(k*real_count),
                                     real_count, nfreq, fgrid, shift, 1e-6*pow,
                                     flock*(k+1), (t < tlock) ? t : tlock, (t+tlock >= nbins) ? nbins-t-1 : tlock);
                }
        }
        free(q);
        free(td);
        free(fd);
}
