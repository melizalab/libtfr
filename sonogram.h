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

#include <fftw3.h>
#include "mtm_tfr.h"

#define LIBSONO_MAXFFTSIZE 131072

struct sonogram
{
	int fftsize;
	int overlap;
	int ntimebins;
	int nfreqbins;
	int specalloctype;
	int fd;
	float **cachetimebins;
	int windowtype;
	float *spec;
	char *tempnam;
	int tempfd;
	int linear;
	mtfft_params *mtm;
	int mtm_ntapers;
	double mtm_nw;
	int mtm_adapt;
	int method;
	double *ffttemppsd;   // workspace for the PSD algorithms
	float *psd;
	float *findex;
	int findexlen;
	float fmin;
	float fmax;
	int *interp_index;
	int interp_index_size;
	int buftype;
	float dcoff;
};

#define SONO_OPT_SCALE          2
#define   SCALE_LINEAR      1
#define   SCALE_DB          2
#define SONO_OPT_WINDOW        3
#define   WINDOW_MULTITAPER  1
#define   WINDOW_HAMMING     2
#define   WINDOW_HANNING     3
#define   WINDOW_BOXCAR      4 
#define   WINDOW_TRIANG      5
#define   WINDOW_BLACKMAN    6
#define SONO_OPT_METHOD        4
#define   METHOD_STFT           1
#define   METHOD_TFR            2
#define SONO_OPT_MTM_NTAPERS   5
#define SONO_OPT_MTM_NW        6
#define SONO_OPT_MTM_ADAPT     7
/* #define SONO_OPT_MTM_INOISE     7 */
/* #define SONO_OPT_MTM_ILOG       8 */
/* #define SONO_OPT_MTM_FSMOOTH    9 */
/* #define SONO_OPT_MTM_INORM     10 */
/* #define SONO_OPT_MTM_ISPEC     11 */
/* #define SONO_OPT_MTM_ITHRESH   12 */
/* #define SONO_OPT_MTM_ISIGNAL   13 */
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

