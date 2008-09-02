/*
** sonogram.h
**
** header file for sonogram.c, an interface for obtaining power-spectral
** density information for a 1-dimensional series of sampled data, using
** FFT, Maximum Entropy, or Multitapering methods.
**
** Copyright (C) by Amish S. Dave 2002
**
** This code is distributed under the GNU GENERAL PUBLIC LICENSE (GPL)
** Version 2 (June 1991). See the "COPYING" file distributed with this software
** for more info.
*/

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#define LIBSONO_MAXFFTSIZE 131072

typedef struct
{
    double *wk1, *wk2, *wkm;
    int n, m;
    double *cof;
    double pm;
    double *wprtheta, *wpitheta;
    int theta_array_len;
    float fmin;
    float fmax;
} memspect;

typedef struct
{
    int p_nwin;
    int p_npi;
    float p_f1;
    float p_f2;
    int p_inoise;
    int p_ilog;
    int p_fsmooth;
    int p_inorm;
    int p_ispec;
    int p_ithresh;
    int p_isignal;
} multitaper;

struct sonogram
{
    int fftsize;
    int overlap;
    int ntimebins;
    int nfreqbins;
    int specalloctype;
    int fd;
    float **cachetimebins;
    float *window;
	int windowtype;
    float *spec;
    float *complexbuf;
    char *tempnam;
    int tempfd;
    int linear;
    memspect *ms;
    multitaper *mtm;
    int npoles;
    int method;
    float *psd;
    float *findex;
    int findexlen;
    float fmin;
    float fmax;
    float *fftw_temp;
    float *ffttemppsd;
    int *interp_index;
    int interp_index_size;
#ifdef HAVE_FFTW
	void *fftw_plan;
	int fftw_plan_nsamp;
	double *fftw_data_in;
	double *fftw_data_out;
#endif
	int buftype;
	float dcoff;
};

#define SONO_OPT_NPOLES         1
#define SONO_OPT_SCALE          2
#define   SCALE_LINEAR      1
#define   SCALE_DB          2
#define SONO_OPT_WINDOW         3
#define   WINDOW_HAMMING    1
#define   WINDOW_HANNING    2
#define   WINDOW_BOXCAR     3
#define   WINDOW_TRIANG     4
#define   WINDOW_BLACKMAN   5
#define SONO_OPT_METHOD         4
#define   METHOD_FFT        1
#define   METHOD_MAXENTROPY 2
#define   METHOD_MULTITAPER 3
#define   METHOD_FFTW       4
#define SONO_OPT_MTM_NWIN       5
#define SONO_OPT_MTM_NPI        6
#define SONO_OPT_MTM_INOISE     7
#define SONO_OPT_MTM_ILOG       8
#define SONO_OPT_MTM_FSMOOTH    9
#define SONO_OPT_MTM_INORM     10
#define SONO_OPT_MTM_ISPEC     11
#define SONO_OPT_MTM_ITHRESH   12
#define SONO_OPT_MTM_ISIGNAL   13
#define SONO_OPT_NFREQBINS     14
#define SONO_OPT_BUFTYPE       15
#define   BUFTYPE_LINEAR    1
#define   BUFTYPE_RING      2
#define SONO_OPT_DCOFF         16
#define SONO_OPT_FFTSIZE       17
#define SONO_OPT_OVERLAP       18
#define SONO_OPT_NSAMPLES      19

#define ALLOCTYPE_MALLOC 1
#define ALLOCTYPE_MMAP   2

struct sonogram *create_sonogram(void);
float *calculate_psd_cached_column(struct sonogram *sono, short *samples, int nsamples, int timebin, float fmin, float fmax);
float *calculate_psd_uncached_column(struct sonogram *sono, short *samples, int nsamples, int timebin, float fmin, float fmax);
void free_sonogram(struct sonogram *sono);
int sonogram_setopts(struct sonogram *sono, int option, long value);
void calculate_psd(struct sonogram *sono, short *samples, int nsamples, int xoffset, int psdlen, float *psd, float fmin, float fmax);

memspect *init_memspect(int n, int m);
void free_memspect(memspect *ms);
void mem_calc_coefficients(memspect *ms, float *data);
double mem_eval_spectrum(memspect *ms, double fdt);
void mem_spectrum(memspect *ms, float *spectrum);
void mem_spectrum2(memspect *ms, float *spectrum, int spectrum_len, float fmin, float fmax);


int mtmpsd(float *buf, int nsamples, int p_nwin, int p_npi, float p_f1, float p_f2, int p_inoise, int p_ilog, float p_fsmooth, int p_inorm, int p_ispec, int p_ithresh, int p_isignal);
void my_multitaper_spectrum(float *data, int npoints, int nwin, int p_nwdt, float *psd, int kind);

