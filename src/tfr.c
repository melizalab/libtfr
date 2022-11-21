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

// N must be odd
int
hermf(int N, int M, double tm, double *h, double *Dh, double *Th)
{

        int i, k;
        double dt, *tt, *P;

        tt = (double*)malloc(N*sizeof(double));
        P = (double*)malloc(N*2*sizeof(double));

        // tt[] are the x values that the functions P, h, Dh, and Th are evaluated at
        // P is the Hermite polynomial
        // h is the Hermite function, Ψₙ = (2ⁿ n! √π)^-½ * e^(-x²/2) * Pₙ

        // The first Hermite polynomial, P₀, is calculated by explicit expression, the
        // orders after that via recurrence.  P₀ is pre-multiplied by the factor
        // e^(-x²/2), as it more efficient than multiplying when h and Dh are computed
        // and means e^(-x²/2) doesn't need to be saved.

        // Dh, the derivative of h, is calculated as:
        // Ψₙ´ = ((2ⁿ n! √π)^-½ * e^(-x²/2) * Pₙ)´            With P = Hermite polynomial
        //     = f * (e^(-x²/2) * Pₙ)´                        Where f = (2ⁿ n! √π)^-½
        //     = f * (e^(-x²/2) * Pₙ´ - x * e^(-x²/2) * Pₙ)
        //                                                    Using Pₙ´ = Pₙ₋₁ * 2n
        //     = f * (e^(-x²/2) * Pₙ₋₁ * 2n - x * e^(-x²/2) * Pₙ)
        //                                                    Where Pk = e^(-x²/2) * Pₙ
        //     = f * (2n * Pkm1 - x * Pk)

        dt = 2 * tm / (N-1);
        // f is (2^k k! √π)^-½ * √dt, for the current k (now k = 0)
        // It's faster to calculate fₖ through recurrence as fₖ = fₖ₋₁ / √(2k)
        double f = sqrt(dt / sqrt(M_PI));
        for (i = 0; i < N; i++) {
                tt[i] = -tm + dt * i;
                // Include e^(-x²/2) factor here, i.e.. P = Hermite(0) * e^(-x²/2) = 1 * e^(-x²/2)
                P[i] = exp(-tt[i] * tt[i] / 2);

                h[i] = P[i] * f;
                Th[i] = h[i] * (-(N-1)/2 + i);
                // Dh = dt * f * (2k * Pkm1 - tt * Pk), with Pkm1 = 0 and h = f * Pk
                Dh[i] = -dt * tt[i] * h[i];
        }

        // Pkm2 points to Pₖ₋₂ and Pkm1 points to Pₖ₋₁, i.e. the two previous orders of
        // P.  Pk will be the current order, it points to Pkm2, so we overwrite Pkm2 as
        // we calculate the values for Pk.  The pointers are shifted each iteration.
        // This way we only need two rows of P values at once.
        double *Pkm2 = P+N, *Pkm1 = P, *Pk = P+N;

        memset(Pkm2, 0, sizeof(*Pkm2) * N);

        // Note that on the first iteration (k=1), Pkm2 is uninitialized, but it is
        // multiplied by (k-1) = 0, so it doesn't matter.
        for (k = 1; k < M; k++) {
                f /= sqrt(2 * k);

                for (i = 0; i < N; i++) {
                        // Pkm2 is P[k-2,], and Pkm1 is P[k-1,], we will overwrite Pkm2 with P[k,]
                        Pk[i] = 2.0 * (tt[i] * Pkm1[i] - (k-1) * Pkm2[i]);

                        h[k*N+i] = Pk[i] * f;
                        Th[k*N+i] = h[k*N+i] * (-(N-1)/2 + i);
                        Dh[k*N+i] = f * dt * (2 * k * Pkm1[i] - tt[i] * Pk[i]);
                }

                Pkm2 = Pkm1;
                Pkm1 = Pk;
                Pk = Pkm2;
        }

        free(tt);
        free(P);

        return N;
}

mfft *
mtm_init_herm(int nfft, int npoints, int order, double tm)
{
        if ((npoints & 1) == 0)
                return NULL;
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
find_bin(double f, const double *fgrid, int nfreq)
{
        int i = 0, j = nfreq - 1;

        if (f < fgrid[i] || f > fgrid[j]) return -1;

        while (j - i > 1) {
                int k = (i + j) / 2;
                if (f > fgrid[k])
                        i = k;
                else
                        j = k;
        }

        return (f - fgrid[i] < fgrid[j] - f) ? i : j;
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
