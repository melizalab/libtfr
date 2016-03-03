/*
 * @file   test_tfr.c
 * @author C Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:35:27 2010
 *
 * Test program for tfr library
 *
 * Copyright C Daniel Meliza 2010.  Licensed for use under GNU
 * General Public License, Version 2.  See COPYING for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tfr.h"

int npoints = 17590;
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
                complex double phase = 2 * M_PI * fnormid * (t - t0) +
                        delta * period * (sin(2 * M_PI * (t - t0) / period + phi) - sin(phi));
                val[t] = creal(cexp(I * phase));
        }
}

/* output a tab-delimited file */
void
write_file(char const * fn, double *buf, int nrow, int ncol)
{
        FILE *fp = fopen(fn, "wt");
        for (int i = 0; i < nrow; ++i) {
                for (int j = 0; j+1 < ncol; ++j) {
                        fprintf(fp, "%3.4f\t", *buf);
                        ++buf;
                }
                fprintf(fp, "%.6g\n", *buf);
                ++buf;
        }
        fclose(fp);
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

        fmsin(sig, npoints, 0.15, 0.45, 1024, 256./4, 0.3, -1);

        printf("* Input signal to tfr_in.dat\n");
        write_file("tfr_in.dat", sig, npoints, 1);

        mtmh = mtm_init_dpss(N, NW, (int)(NW*2-1));

        printf("* MTM PSD to tfr_out_psd\n");
        psd = (double*)malloc(N * sizeof(double));
        sigpow = mtfft(mtmh, sig+8300, N);
        mtpower(mtmh, psd, sigpow);
        write_file("tfr_out_psd.dat", psd, N/2 + 1, 1);
        free(psd);

        const int l = (npoints - Np + 1) / step;
        printf("* MTM spectrogram to tfr_out_mtm\n");
        specgram = (double*)calloc(l * (N/2+1), sizeof(double));
        mtm_spec(mtmh, specgram, sig, npoints, step, 1);
        write_file("tfr_out_mtm.dat", specgram, l, (N/2+1));
        free(specgram);

        mtm_destroy(mtmh);

        printf("* TFR spectrogram to tfr_out_tfr\n");
        mtmh = mtm_init_herm(N, Np, k, tm);
        specgram = (double*)calloc(l * (N/2+1), sizeof(double));
        tfr_spec(mtmh, specgram, sig, npoints, -1, step, 0.01, 5, 0, NULL);
        write_file("tfr_out_tfr.dat", specgram, l, (N/2+1));
        free(specgram);

        mtm_destroy(mtmh);

        free(sig);
}
