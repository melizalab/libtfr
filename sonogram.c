/*
** sonogram.c
**
** Interface for obtaining power-spectral density information for a
** 1-dimensional series of sampled data, using FFT, Maximum Entropy,
** or Multitapering methods.
**
** Copyright (C) by Amish S. Dave 2002
**
** This code is distributed under the GNU GENERAL PUBLIC LICENSE (GPL)
** Version 2 (June 1991). See the "COPYING" file distributed with this software
** for more info.
*/


/****************************************************************************************
** Includes
****************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <errno.h>
#include <float.h>

#include "sonogram.h"

#define MAXFFTSIZE LIBSONO_MAXFFTSIZE

#ifdef linux
#ifndef MADV_SEQUENTIAL
#define MADV_SEQUENTIAL
#define madvise(addr,len,flags)
#endif
#endif

#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

#ifndef MINFLOAT
#define MINFLOAT FLT_MIN
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef PI
#define PI      3.1415926537
#endif


/****************************************************************************************
** Prototypes
****************************************************************************************/
static int init_window(int windowtype, double *window, int fftsize);
static void init_hamming_window(double *hammingwindow, int fftsize);
static void init_hanning_window(double *window, int fftsize);
static void init_boxcar_window(double *window, int fftsize);
static void init_triang_window(double *window, int fftsize);
static void init_blackman_window(double *window, int fftsize);

static void calculate_stftpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, double *psd);

static void create_cache(struct sonogram *sono);
static void free_cache(struct sonogram *sono);
static void clear_cache(struct sonogram *sono);


/****************************************************************************************
** Configuration
****************************************************************************************/

/*
** This routine allocates and initialized the sonogram structure which will access
** the FFT data used to draw the sonogram.
** We want to support the following features:
**	1) use of either malloc or mmap to allocate memory
**		for example, a 20kHz, 100s, 256pt, 90% overlap float sonogram needs ~60Mb!!!
**	2) only calculate the FFTs as they are needed.
**		for example, when plotting the above FFT on a 1000 pixel wide window, why
**		calculate the umpteen thousand ffts when you can make do with only 1000 !!
**	3) Eventually, be able to spawn a background thread (or maybe use X11 work functions)
**		to calculate the uncalculated FFTs in the background, before they are needed, thus
**		speeding up zooms and panning.
**	4) NOTE: the time bins are centered on the sample: so, the first time bin is the 
**		result of computing the FFT of a window of data containing the first 1/2 window
**		size samples, with the other half of the window zero padded.  Because this is a
**		special case, this timebin is always computed right here.
**	5) There may be some thrown away data at the end, if nsamples % inc != 0.
**		In this case, the last window is zero padded also.
*/
struct sonogram*
create_sonogram(void)
{
	struct sonogram *sono;

	sono = (struct sonogram *)calloc(1, sizeof(struct sonogram));

	/*
	** THIS MUST BE SET VIA sonogram_setopts()
	*/
	sono->fftsize = 0;

	sono->overlap = 0;
	sono->ntimebins = 0;
	sono->nfreqbins = 0;
	sono->tempfd = -1;
	sono->linear = 0;
	sono->windowtype = WINDOW_HAMMING;
	sono->cachetimebins = NULL;
	sono->spec = NULL;
	sono->fd = -1;
	sono->method = METHOD_STFT;
	sono->mtm = NULL;
	sono->mtm_ntapers = 3;
	sono->mtm_nw = 2.5;
	sono->mtm_adapt = 1;
	sono->psd = NULL;
	sono->fmin = 0.0;
	sono->fmax = 0.5;
	sono->ffttemppsd = NULL;
	sono->interp_index = NULL;
	sono->interp_index_size = 0;
	sono->buftype = BUFTYPE_LINEAR;
	sono->dcoff = 0.0;
	return sono;
}

void 
free_sonogram(struct sonogram *sono)
{
	if (sono->interp_index) free(sono->interp_index);
	if (sono->ffttemppsd) free(sono->ffttemppsd);
	if (sono->spec) free_cache(sono);
	if (sono->mtm) mtm_destroy(sono->mtm);
	if (sono->psd) free(sono->psd);
	/* Zero these out, to trigger segfault on access after free() */
	memset(sono, 0, sizeof(struct sonogram));
	free(sono);
}

int 
sonogram_setopts(struct sonogram *sono, int option, long value)
{
	int inc;
	int clearplan = 0;
	int clearcache = 0;

	switch(option)
	{
		case SONO_OPT_METHOD:
			if (sono->method == value) break;
			sono->method = value;
			clearplan = 1;
			break;
		case SONO_OPT_NSAMPLES:
			inc = ((100.0 - sono->overlap) / 100.0) * sono->fftsize;
			sono->ntimebins = (value + inc - 1) / inc;
			break;
		case SONO_OPT_OVERLAP:
			sono->overlap = value;
			clearcache = 1;
			break;
		case SONO_OPT_FFTSIZE:
			if (sono->fftsize == value) break;
			sono->fftsize = value;
			sono->nfreqbins = sono->fftsize / 2;
			clearcache = 1;
			clearplan = 1;
			break;
		case SONO_OPT_BUFTYPE:
			sono->buftype = value;
			break;
		case SONO_OPT_DCOFF:
			sono->dcoff = (float)value;
			break;
		case SONO_OPT_NFREQBINS:
			if (sono->nfreqbins == value)
				break;
			sono->nfreqbins = value;
			clearcache = 1;
			break;
		case SONO_OPT_SCALE:
			sono->linear = (value == SCALE_LINEAR) ? (1) : (0);
			break;
		case SONO_OPT_WINDOW:
			sono->windowtype = value;
			clearplan = 1;
			break;
		case SONO_OPT_MTM_NTAPERS:
			if (sono->mtm_ntapers == value) break;
			sono->mtm_ntapers = value;
			clearplan = 1;
			break;
		case SONO_OPT_MTM_NW:
			if (sono->mtm_nw == value) break;
			sono->mtm_nw = value;
			clearplan = 1;
			break;
		case SONO_OPT_MTM_ADAPT:
			sono->mtm_adapt = value;
			break;
	}
	if (clearcache) {
		if (sono->spec)
			free_cache(sono);
		sono->spec = NULL;
		sono->cachetimebins = NULL;
	}
	if (clearplan) {
		if (sono->ffttemppsd) free(sono->ffttemppsd);
		sono->ffttemppsd = NULL;
		if (sono->interp_index) free(sono->interp_index);
		sono->interp_index = NULL;
		sono->interp_index_size = 0;

		if (sono->mtm) mtm_destroy(sono->mtm);
		sono->mtm = NULL;

	}
	return 0;
}


float *calculate_psd_cached_column(struct sonogram *sono, short *samples, int nsamples, int timebin, float fmin, float fmax)
{
	int xoffset;
	float *spec, **cachetimebins, *outptr, inc;

	if (sono->spec == NULL)
		create_cache(sono);
	if (sono->spec == NULL)
		return NULL;

	spec = sono->spec;
	cachetimebins = sono->cachetimebins;
	inc = ((100.0 - sono->overlap) / 100.0) * sono->fftsize;
	xoffset = (timebin * inc);

	if ((ABS(fmin - sono->fmin) > 0.001) || (ABS(fmax - sono->fmax) > 0.001))
	{
		sono->fmin = fmin;
		sono->fmax = fmax;
		/*
		** Clear cache references (works even if the cache (sono->spec)
		** is not allocated or being used
		*/
		memset(sono->cachetimebins, 0, sono->ntimebins * sizeof(float *));
		if (sono->interp_index)
			free(sono->interp_index);
		sono->interp_index = NULL;
	}

	if (cachetimebins[timebin] != NULL)
		return cachetimebins[timebin];

	outptr = spec + timebin * sono->nfreqbins;
	calculate_psd(sono, samples, nsamples, xoffset, sono->nfreqbins, outptr, fmin, fmax);
	cachetimebins[timebin] = outptr;
	return outptr;
}

float *calculate_psd_uncached_column(struct sonogram *sono, short *samples, int nsamples, int timebin, float fmin, float fmax)
{
	int sample;
	float *outptr, inc;

	inc = ((100.0 - sono->overlap) / 100.0) * sono->fftsize;
	sample = (timebin * inc);

	//printf("NFFT: %d; Overlap %d\n", sono->fftsize, sono->overlap);
	//printf("Step size: %3.2f; Sample %d\n", inc, sample);
	if (sono->psd == NULL)
	{
		//printf("Allocating %d frequency points\n", sono->nfreqbins);
		sono->psd = (float *)calloc(sono->nfreqbins, sizeof(float));
	}
	outptr = sono->psd;
	calculate_psd(sono, samples, nsamples, sample, sono->nfreqbins, outptr, fmin, fmax);
	return outptr;
}

void calculate_psd(struct sonogram *sono, short *samples, int nsamples, int xoffset, int psdlen, float *psd, float fmin, float fmax)
{
	int y;
	float fres;

	/*
	** Did our y-axis change?
	*/
	if ((ABS(fmin - sono->fmin) > 0.001) || (ABS(fmax - sono->fmax) > 0.001))
	{
		sono->fmin = fmin;
		sono->fmax = fmax;

		/*
		** Clear cache references (works even if the cache
		** is not allocated or being used
		*/
		clear_cache(sono);
		if (sono->interp_index)
			free(sono->interp_index);
		sono->interp_index = NULL;
	}

	/*
	** The following won't happen if we're called from calculate_psd_cached_column...
	** If it did, we'd have problems, since psd would be a pointer into the cache.
	*/
	if (psdlen != sono->nfreqbins)
	{
		if (sono->spec)
			free_cache(sono);
		sono->nfreqbins = psdlen;
	}

	/*
	** Now, we call the FFT routine (which writes over the input)
	*/
	if (sono->method == METHOD_STFT)
	{
		if (sono->ffttemppsd == NULL) {
			//printf("Allocating psd cache, %d samples.\n", sono->fftsize/2+1);
			sono->ffttemppsd = (double *)calloc((sono->fftsize/2)+1, sizeof(double));
		}
		if ((sono->interp_index != NULL) && (sono->interp_index_size != psdlen))
		{
			free(sono->interp_index);
			sono->interp_index = NULL;
		}
		if (sono->interp_index == NULL)
		{
			sono->interp_index = (int *)malloc(psdlen * sizeof(int));
			sono->interp_index_size = psdlen;
			fres = ((sono->fftsize * (fmax - fmin)) - 1.0) / ((float)psdlen - 1);
			for (y=0; y < psdlen; y++)
				sono->interp_index[y] = (fmin * sono->fftsize) + fres * y;
		}
		//printf("Calculating PSD...\n");
		calculate_stftpsd(sono, samples, nsamples, xoffset, sono->ffttemppsd);
		for (y=0; y < psdlen; y++)
			psd[y] = (float)sono->ffttemppsd[sono->interp_index[y]];
	}
	return;
}

static void calculate_stftpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, double *psd)
{
	int i, fftsize;
	double sigpow;

	fftsize = sono->fftsize;

	// initialize tapers/windows and plan if necessary
	if ((sono->mtm == NULL) || fftsize != sono->mtm->nfft) {
		printf("Initializing new plan\n");
		if (sono->mtm) mtm_destroy(sono->mtm);
		if ((sono->windowtype) == WINDOW_MULTITAPER)
			sono->mtm = mtm_init_dpss(fftsize, sono->mtm_nw, sono->mtm_ntapers);
		else {
			double *window = (double*)malloc(fftsize*sizeof(double));
			//printf("Generating window of size %d\n", fftsize);
			init_window(sono->windowtype, window, fftsize);
			//printf("Initializing plan...\n");
			sono->mtm = mtm_init(fftsize, fftsize, 1, window, NULL);
		}
	}

	//printf("Computing FFT\n");
	sigpow = mtfft(sono->mtm, samples+xoffset, nsamples);
	//printf("Frame %d; Signal power: %3.2f\n", xoffset, sigpow);
	if (!sono->mtm_adapt) sigpow = 0.0;
	mtpower(sono->mtm, psd, sigpow);

	if (sono->linear) return;
	for (i = 0; i < (fftsize + 1) / 2; i++)
		psd[i] = 10.0 * log10(psd[i] + DBL_EPSILON);
}


/****************************************************************************************
********   WINDOW FUNCTIONS
****************************************************************************************/
static int init_window(int windowtype, double *window, int fftsize)
{
	switch(windowtype)
	{
		case WINDOW_HAMMING:
			init_hamming_window(window, fftsize);
			break;
		case WINDOW_HANNING:
			init_hanning_window(window, fftsize);
			break;
		case WINDOW_BOXCAR:
			init_boxcar_window(window, fftsize);
			break;
		case WINDOW_TRIANG:
			init_triang_window(window, fftsize);
			break;
		case WINDOW_BLACKMAN:
			init_blackman_window(window, fftsize);
			break;
		default:
			return -1;
			break;
	};
	// normalize the window
	renormalize(fftsize, window);
	return 0;
}

static void init_hamming_window(double *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * M_PI / fftsize;
	for (i=0; i < fftsize; i++)
		window[i] = 0.54 - 0.46*(cos(K*i));
}

static void init_hanning_window(double *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * M_PI / fftsize;
	for (i=0; i < fftsize; i++)
		window[i] = 0.5*(1.0-cos(K*i));
}

static void init_boxcar_window(double *window, int fftsize)
{
	int i;

	for (i=0; i < fftsize; i++)
		window[i] = 1.0;
}

static void init_triang_window(double *window, int fftsize)
{
	int i, half_fftsize;

	half_fftsize = fftsize / 2;
	if ((fftsize % 2) == 0)
	{
		for (i=0; i < half_fftsize; i++)
			window[i] = 1.0 * (i / ((double)(half_fftsize - 1)));
	}
	else
	{
		for (i=0; i <= half_fftsize; i++)
			window[i] = 1.0 * (i / ((double)half_fftsize));
	}
	for (; i < fftsize; i++)
		window[i] = window[fftsize - i - 1];
}

static void init_blackman_window(double *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * M_PI / (fftsize - 1);
	for (i=0; i < fftsize; i++)
		window[i] = (0.42 - 0.5 * cos(K*i) + .08 * cos(2*K*i));
}


/****************************************************************************************
********   UTILITY FUNCTIONS
****************************************************************************************/
static void create_cache(struct sonogram *sono)
{
	size_t len;

	if (sono->spec != NULL)
		free_cache(sono);

	/*
	** For large allocations, create a temporary file and mmap it
	** If this fails, or for small files, malloc the needed space.
	*/
	sono->spec = NULL;
	sono->specalloctype = 0;
	if (sono->ntimebins * sono->nfreqbins > (10000000 / sizeof(float)))
	{
		sono->tempnam = strdup("/tmp/fftXXXXXX");
		if ((sono->tempfd = mkstemp(sono->tempnam)) != -1)
		{
			unlink(sono->tempnam);
			len = sono->ntimebins * sono->nfreqbins * sizeof(float);
			ftruncate(sono->tempfd, len);
			lseek(sono->tempfd, 0, SEEK_SET);
			sono->spec = (float *)mmap(NULL, len, PROT_READ|PROT_WRITE, MAP_PRIVATE, sono->tempfd, 0);
			if ((void *)(sono->spec) != (void *)-1)
			{
				madvise((void *)(sono->spec), len, MADV_SEQUENTIAL);
				sono->specalloctype = ALLOCTYPE_MMAP;
			}
			else
			{
				sono->spec = NULL;
				close(sono->tempfd);
				sono->tempfd = -1;
				sono->specalloctype = 0;
			}
		}
		free(sono->tempnam);
		sono->tempnam = NULL;
	}
	if (sono->spec == NULL)
	{
		if ((sono->spec = (float *)calloc(sono->ntimebins * sono->nfreqbins, sizeof(float))) != NULL)
			sono->specalloctype = ALLOCTYPE_MALLOC;
	}
	if (sono->cachetimebins == NULL)
		sono->cachetimebins = (float **)calloc(sono->ntimebins, sizeof(float *));
}

/*
** This will cause all timebins to be re-calculated when
** needed, instead of returned from the cache
*/
static void clear_cache(struct sonogram *sono)
{
	if (sono->cachetimebins)
		memset(sono->cachetimebins, 0, sono->ntimebins * sizeof(float *));
}

static void free_cache(struct sonogram *sono)
{
	if (sono->spec == NULL)
		return;
	if (sono->specalloctype == ALLOCTYPE_MALLOC)
	{
		free(sono->spec);
	}
	else if (sono->specalloctype == ALLOCTYPE_MMAP)
	{
		munmap((void *)(sono->spec), sono->ntimebins * sono->nfreqbins * sizeof(float));
		if (sono->tempnam != NULL) free(sono->tempnam);
		if (sono->tempfd != -1) close(sono->tempfd);
	}
	sono->spec = NULL;
	if (sono->cachetimebins) free(sono->cachetimebins);
	sono->cachetimebins = NULL;
}

