/*
 * @file   test_tfr.c
 * @author C Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:35:27 2010
 *
 * Test program for tfr library. Returns nonzero if any check fails.
 *
 * Pass any argument to also dump the demo spectrograms to .dat files in the
 * working directory; by default nothing is written, so CI leaves no clutter.
 *
 * Copyright C Daniel Meliza 2010.  Licensed for use under GNU
 * General Public License, Version 2.  See COPYING for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// tfr.h deliberately omits <complex.h>, so include it here: this is a plain
// C99 program that uses complex arithmetic and is never built with MSVC.
#include <complex.h>
#include "tfr.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

int npoints = 17590;
int N = 256;
int Np = 201;
double NW = 3.5;
int step = 10;
int k = 6;
double tm = 6.0;

static int failures = 0;
static int checks = 0;

#define CHECK(cond)                                                     \
        do {                                                            \
                ++checks;                                               \
                if (!(cond)) {                                          \
                        fprintf(stderr, "%s:%d: FAIL %s\n",             \
                                __FILE__, __LINE__, #cond);             \
                        ++failures;                                     \
                }                                                       \
        } while (0)

#define CHECK_CLOSE(got, want, tol)                                     \
        do {                                                            \
                ++checks;                                               \
                if (!(fabs((got) - (want)) <= (tol))) {                 \
                        fprintf(stderr,                                 \
                                "%s:%d: FAIL %s == %s (%.17g vs %.17g)\n", \
                                __FILE__, __LINE__, #got, #want,        \
                                (double)(got), (double)(want));         \
                        ++failures;                                     \
                }                                                       \
        } while (0)

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

/* a pure tone, for tests that need a known peak frequency */
static void
tone(double *val, int n, double f0, double fs)
{
        for (int t = 0; t < n; t++)
                val[t] = sin(2 * M_PI * f0 * t / fs);
}

/* output a tab-delimited file */
void
write_file(char const * fn, double *buf, int nrow, int ncol)
{
        FILE *fp = fopen(fn, "wt");
        if (fp == NULL) {
                fprintf(stderr, "unable to open %s for writing\n", fn);
                return;
        }
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

static void
test_mtm_accessors(void)
{
        mfft *mtm = mtm_init_dpss(N, Np, NW, (int)(NW*2-1));
        CHECK(mtm != NULL);
        if (mtm == NULL) return;

        CHECK(mtm_nfft(mtm) == N);
        CHECK(mtm_npoints(mtm) == Np);
        CHECK(mtm_ntapers(mtm) == (int)(NW*2-1));
        CHECK(mtm_nreal(mtm) == N/2 + 1);
        // frames with full support in the signal
        CHECK(mtm_nframes(mtm, npoints, step) == (npoints - Np)/step + 1);
        CHECK(mtm_buffer(mtm) != NULL);
        CHECK(mtm_tapers(mtm) != NULL);

        mtm_destroy(mtm);
        // documented to be a no-op on NULL
        mtm_destroy(NULL);
}

static void
test_dpss(void)
{
        const int np = 256, ntapers = 5;
        double *tapers = (double*)malloc(np * ntapers * sizeof(double));
        double *lambda = (double*)malloc(ntapers * sizeof(double));

        CHECK(dpss(tapers, lambda, np, 3.0, ntapers) == 0);
        for (int t = 0; t < ntapers; t++) {
                // eigenvalues are concentrations, so they sit in (0, 1] and
                // decrease with taper order
                CHECK(lambda[t] > 0.0 && lambda[t] <= 1.0);
                if (t > 0) CHECK(lambda[t] <= lambda[t-1]);
                // each taper is unit norm
                double ss = 0.0;
                for (int j = 0; j < np; j++) ss += tapers[t*np+j] * tapers[t*np+j];
                CHECK_CLOSE(ss, 1.0, 1e-9);
        }
        // nw out of range is rejected
        CHECK(dpss(tapers, lambda, np, -5.0, ntapers) != 0);

        free(tapers);
        free(lambda);
}

/* mtcomplex unpacks FFTW's half-complex buffer into interleaved pairs. Nothing
   in the python extension calls it -- it builds its own equivalent -- so this
   is the only coverage it gets. */
static void
test_mtcomplex(void)
{
        const int nfft = 256, np = 201, ntapers = 5;
        const int imag_count = (nfft + 1) / 2;
        mfft *mtm = mtm_init_dpss(nfft, np, 3.5, ntapers);
        CHECK(mtm != NULL);
        if (mtm == NULL) return;

        double *sig = (double*)malloc(np * sizeof(double));
        tone(sig, np, 1000.0, 8000.0);
        mtfft(mtm, sig, np);

        const int nreal = mtm_nreal(mtm);
        const double *buf = mtm_buffer(mtm);
        cmplx_t *z = (cmplx_t*)malloc(nreal * ntapers * sizeof(cmplx_t));
        mtcomplex(mtm, z);

        for (int t = 0; t < ntapers; t++) {
                for (int n = 0; n < nreal; n++) {
                        // real parts ascend from index 0
                        CHECK_CLOSE(z[t*nreal+n][0], buf[t*nfft+n], 0.0);
                        // imaginary parts descend from index nfft-1; DC and
                        // (for even nfft) Nyquist have none
                        double want = (n > 0 && n < imag_count) ? buf[t*nfft+(nfft-n)] : 0.0;
                        CHECK_CLOSE(z[t*nreal+n][1], want, 0.0);
                }
                CHECK(z[t*nreal][1] == 0.0);
                CHECK(z[t*nreal + nreal-1][1] == 0.0);
        }

        free(z);
        free(sig);
        mtm_destroy(mtm);
}

/* mtm_zspec should equal mtfft+mtcomplex applied frame by frame */
static void
test_mtm_zspec(void)
{
        const int nfft = 256, np = 201, ntapers = 5, shift = 64, nsamples = 2048;
        mfft *mtm = mtm_init_dpss(nfft, np, 3.5, ntapers);
        CHECK(mtm != NULL);
        if (mtm == NULL) return;

        double *sig = (double*)malloc(nsamples * sizeof(double));
        tone(sig, nsamples, 1000.0, 8000.0);

        const int nreal = mtm_nreal(mtm);
        const int nframes = mtm_nframes(mtm, nsamples, shift);
        CHECK(nframes > 1);

        cmplx_t *spec = (cmplx_t*)calloc(nframes * ntapers * nreal, sizeof(cmplx_t));
        mtm_zspec(mtm, spec, sig, nsamples, shift);

        cmplx_t *frame = (cmplx_t*)malloc(ntapers * nreal * sizeof(cmplx_t));
        for (int f = 0; f < nframes; f++) {
                mtfft(mtm, sig + f*shift, nsamples - f*shift);
                mtcomplex(mtm, frame);
                for (int i = 0; i < ntapers * nreal; i++) {
                        CHECK_CLOSE(spec[f*ntapers*nreal + i][0], frame[i][0], 0.0);
                        CHECK_CLOSE(spec[f*ntapers*nreal + i][1], frame[i][1], 0.0);
                }
        }

        free(frame);
        free(spec);
        free(sig);
        mtm_destroy(mtm);
}

/* a pure tone must reassign onto its own frequency bin */
static void
test_tfr_concentration(void)
{
        const double fs = 8000.0, f0 = 1000.0;
        const int nfft = 256, np = 201, order = 6, shift = 64, nsamples = 8192;
        const int nreal = nfft/2 + 1;
        const int expected = (int)(f0 / fs * nfft);

        mfft *mtm = mtm_init_herm(nfft, np, order, tm);
        CHECK(mtm != NULL);
        if (mtm == NULL) return;

        double *sig = (double*)malloc(nsamples * sizeof(double));
        tone(sig, nsamples, f0, fs);
        const int nframes = mtm_nframes(mtm, nsamples, shift);
        // tfr_spec writes [frame][freq]
        double *spec = (double*)calloc(nframes * nreal, sizeof(double));
        tfr_spec(mtm, spec, sig, nsamples, -1, shift, 0.01, 5, 0, NULL);

        double *power = (double*)calloc(nreal, sizeof(double));
        for (int f = 0; f < nframes; f++)
                for (int n = 0; n < nreal; n++)
                        power[n] += spec[f*nreal + n] / nframes;

        int peak = 0;
        double total = 0.0;
        for (int n = 0; n < nreal; n++) {
                if (power[n] > power[peak]) peak = n;
                total += power[n];
                CHECK(isfinite(power[n]));
        }
        CHECK(total > 0.0);
        CHECK(peak == expected);
        // essentially all the energy lands within a bin of the peak
        double near = 0.0;
        for (int n = (peak > 0 ? peak-1 : 0); n <= peak+1 && n < nreal; n++)
                near += power[n];
        CHECK(near / total > 0.99);

        free(power);
        free(spec);
        free(sig);
        mtm_destroy(mtm);
}

/* a silent signal has no power to reassign, and must not produce NaN */
static void
test_tfr_silent(void)
{
        const int nfft = 256, np = 201, order = 6, shift = 64, nsamples = 4096;
        const int nreal = nfft/2 + 1;
        mfft *mtm = mtm_init_herm(nfft, np, order, tm);
        CHECK(mtm != NULL);
        if (mtm == NULL) return;

        double *sig = (double*)calloc(nsamples, sizeof(double));
        const int nframes = mtm_nframes(mtm, nsamples, shift);
        double *spec = (double*)calloc(nframes * nreal, sizeof(double));
        tfr_spec(mtm, spec, sig, nsamples, -1, shift, 0.01, 5, 0, NULL);

        for (int i = 0; i < nframes * nreal; i++) {
                CHECK(isfinite(spec[i]));
                CHECK(spec[i] == 0.0);
        }

        free(spec);
        free(sig);
        mtm_destroy(mtm);
}

/* end-to-end pass over a chirp, covering the stages in combination rather than
   in isolation; optionally dumps the result for eyeballing */
static void
demo(int write_output)
{
        mfft *mtmh;
        double *sig, *psd, *specgram;
        double sigpow;

        sig = (double*)malloc(npoints * sizeof(double));
        // every sample has to be initialized: the transforms below read all
        // npoints of it, not just the first Np
        fmsin(sig, npoints, 0.15, 0.45, 1024, 256./4, 0.3, -1);
        if (write_output) write_file("tfr_in.dat", sig, npoints, 1);

        mtmh = mtm_init_dpss(N, Np, NW, (int)(NW*2-1));
        CHECK(mtmh != NULL);
        if (mtmh == NULL) { free(sig); return; }

        psd = (double*)malloc((N/2 + 1) * sizeof(double));
        sigpow = mtfft(mtmh, sig+8300, N);
        CHECK(sigpow > 0.0);
        mtpower(mtmh, psd, sigpow);
        for (int i = 0; i < N/2 + 1; i++) {
                CHECK(isfinite(psd[i]));
                CHECK(psd[i] >= 0.0);
        }
        if (write_output) write_file("tfr_out_psd.dat", psd, N/2 + 1, 1);
        free(psd);

        const int l = mtm_nframes(mtmh, npoints, step);
        specgram = (double*)calloc(l * (N/2+1), sizeof(double));
        mtm_spec(mtmh, specgram, sig, npoints, step, 1);
        for (int i = 0; i < l * (N/2+1); i++) CHECK(isfinite(specgram[i]));
        if (write_output) write_file("tfr_out_mtm.dat", specgram, l, (N/2+1));
        free(specgram);
        mtm_destroy(mtmh);

        mtmh = mtm_init_herm(N, Np, k, tm);
        CHECK(mtmh != NULL);
        if (mtmh == NULL) { free(sig); return; }
        specgram = (double*)calloc(l * (N/2+1), sizeof(double));
        tfr_spec(mtmh, specgram, sig, npoints, -1, step, 0.01, 5, 0, NULL);
        for (int i = 0; i < l * (N/2+1); i++) CHECK(isfinite(specgram[i]));
        if (write_output) write_file("tfr_out_tfr.dat", specgram, l, (N/2+1));
        free(specgram);
        mtm_destroy(mtmh);

        free(sig);
}

int
main(int argc, char **argv)
{
        printf("* Testing TFR library:\n");
        printf("* N = %d\n", N);
        printf("* shift = %d\n", step);
        printf("* MTM NW = %3.2f\n", NW);
        printf("* TFR Np = %d\n", Np);
        printf("* TFR k = %d\n", k);
        printf("* TFR tm = %3.2f\n", tm);

        test_mtm_accessors();
        test_dpss();
        test_mtcomplex();
        test_mtm_zspec();
        test_tfr_concentration();
        test_tfr_silent();
        demo(argc > 1);

        printf("* %d checks, %d failures\n", checks, failures);
        return failures != 0;
}
