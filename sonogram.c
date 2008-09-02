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
static int init_window(int windowtype, float *window, int fftsize);
static void init_hamming_window(float *hammingwindow, int fftsize);
static void init_hanning_window(float *window, int fftsize);
static void init_boxcar_window(float *window, int fftsize);
static void init_triang_window(float *window, int fftsize);
static void init_blackman_window(float *window, int fftsize);

#ifdef HAVE_FFTW
static void calculate_fftwpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd);
#endif
static void calculate_fftpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd);
static void calculate_mempsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, int psdlen, float *psd, float fmin, float fmax);
static void calculate_mtmpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd);
static void create_cache(struct sonogram *sono);
static void free_cache(struct sonogram *sono);
static void clear_cache(struct sonogram *sono);
static void window_samples(float *dest, float *window, short *samples, int nsamples, int offset, int fftsize, float dcoff);
static void window_samples_from_ringbuf(float *dest, float *window, short *samples, int nsamples, int offset, int fftsize, float dcoff);
static multitaper *alloc_mtmstruct(void);


/****************************************************************************************
** Configuration
****************************************************************************************/

#define FFTINITROUTINE(fftsize)
extern void rfft(float x[], int length);
#define FFTROUTINE(complexbuf, fftsize)		rfft(complexbuf, fftsize);

#if defined(sun) || defined(linux)
#define log10f(num) ((float)log10((float)(num)))
#endif

#ifndef SQR
#define SQR(a) ( (a) * (a) )
#endif


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
struct sonogram *create_sonogram(void)
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
	sono->window = (float *)calloc(MAXFFTSIZE, sizeof(float));
	sono->cachetimebins = NULL;
	sono->complexbuf = (float *)calloc(2 * MAXFFTSIZE + 1, sizeof(float));
	sono->spec = NULL;
	sono->fd = -1;
	sono->method = METHOD_FFT;
	sono->npoles = 40;
	sono->ms = NULL;
	sono->mtm = alloc_mtmstruct();
	sono->psd = NULL;
	sono->fmin = 0.0;
	sono->fmax = 0.5;
	sono->ffttemppsd = NULL;
	sono->interp_index = NULL;
	sono->interp_index_size = 0;
#ifdef HAVE_FFTW
	sono->fftw_plan = NULL;
	sono->fftw_plan_nsamp = 0;
	sono->fftw_data_in = (double *)fftw_malloc((2 * MAXFFTSIZE + 2) * sizeof(double));
	sono->fftw_data_out = (double *)fftw_malloc((2 * MAXFFTSIZE + 2) * sizeof(double));
#endif
	sono->buftype = BUFTYPE_LINEAR;
	sono->dcoff = 0.0;
	return sono;
}

void free_sonogram(struct sonogram *sono)
{
	if (sono->complexbuf) free(sono->complexbuf);
	if (sono->window) free(sono->window);
	if (sono->interp_index) free(sono->interp_index);
	if (sono->ffttemppsd) free(sono->ffttemppsd);
	if (sono->spec) free_cache(sono);
	if (sono->ms) free_memspect(sono->ms);
	if (sono->mtm) free(sono->mtm);
	if (sono->psd) free(sono->psd);
#ifdef HAVE_FFTW
	if (sono->fftw_data_in) fftw_free(sono->fftw_data_in);
	if (sono->fftw_data_out) fftw_free(sono->fftw_data_out);
	if (sono->fftw_plan) fftw_destroy_plan(sono->fftw_plan);
#endif
	/* Zero these out, to trigger segfault on access after free() */
	memset(sono, 0, sizeof(struct sonogram));
	free(sono);
}

int sonogram_setopts(struct sonogram *sono, int option, long value)
{
	int inc;

	switch(option)
	{
		case SONO_OPT_METHOD:
			sono->method = value;
			break;
		case SONO_OPT_NSAMPLES:
			inc = ((100.0 - sono->overlap) / 100.0) * sono->fftsize;
			sono->ntimebins = (value + inc - 1) / inc;
			break;
		case SONO_OPT_OVERLAP:
			sono->overlap = value;

			/* Free the cache */
			if (sono->spec)
				free_cache(sono);
			sono->spec = NULL;
			sono->cachetimebins = NULL;

			break;
		case SONO_OPT_FFTSIZE:
			sono->fftsize = value;
			sono->nfreqbins = sono->fftsize / 2;

			/* Free the cache */
			if (sono->spec)
				free_cache(sono);
			sono->spec = NULL;
			sono->cachetimebins = NULL;

			if (sono->ffttemppsd) free(sono->ffttemppsd);
			sono->ffttemppsd = NULL;
			if (sono->interp_index) free(sono->interp_index);
			sono->interp_index = NULL;
			sono->interp_index_size = 0;

			init_window(sono->windowtype, sono->window, sono->fftsize);
			FFTINITROUTINE(sono->fftsize);

#ifdef HAVE_FFTW
			if (sono->fftw_plan) fftw_destroy_plan(sono->fftw_plan);
			sono->fftw_plan = NULL;
			sono->fftw_plan_nsamp = 0;
#endif
			if (sono->ms) free_memspect(sono->ms);
			sono->ms = NULL;
			break;
		case SONO_OPT_BUFTYPE:
			sono->buftype = value;
			break;
		case SONO_OPT_DCOFF:
			sono->dcoff = (float)value;
			break;
		case SONO_OPT_NPOLES:
			sono->npoles = (int)value;
			if (sono->ms) free_memspect(sono->ms);
			sono->ms = NULL;
			break;
		case SONO_OPT_NFREQBINS:
			if (sono->nfreqbins == value)
				break;
			if (sono->spec)
				free_cache(sono);
			sono->nfreqbins = value;
			break;
		case SONO_OPT_SCALE:
			sono->linear = (value == SCALE_LINEAR) ? (1) : (0);
			break;
		case SONO_OPT_WINDOW:
			sono->windowtype = value;
			return init_window(sono->windowtype, sono->window, sono->fftsize);
			break;
		case SONO_OPT_MTM_NWIN:
			if (sono->mtm) sono->mtm->p_nwin = value;
			break;
		case SONO_OPT_MTM_NPI:
			if (sono->mtm) sono->mtm->p_npi = value;
			break;
		case SONO_OPT_MTM_INOISE:
			if (sono->mtm) sono->mtm->p_inoise = value;
			break;
		case SONO_OPT_MTM_ILOG:
			if (sono->mtm) sono->mtm->p_ilog = value;
			break;
		case SONO_OPT_MTM_FSMOOTH:
			if (sono->mtm) sono->mtm->p_fsmooth = value;
			break;
		case SONO_OPT_MTM_INORM:
			if (sono->mtm) sono->mtm->p_inorm = value;
			break;
		case SONO_OPT_MTM_ISPEC:
			if (sono->mtm) sono->mtm->p_ispec = value;
			break;
		case SONO_OPT_MTM_ITHRESH:
			if (sono->mtm) sono->mtm->p_ithresh = value;
			break;
		case SONO_OPT_MTM_ISIGNAL:
			if (sono->mtm) sono->mtm->p_isignal = value;
			break;
	}
	return 0;
}


/*
*****
** EXAMPLE:
** fftsize = 4, overlap = 50%, nsamples = 9
** 
** bin		windowed samples
** ---		----------------
** 0		#, #, 0, 1
** 1		0, 1, 2, 3
** 2		2, 3, 4, 5
** 3		4, 5, 6, 7
** 4		6, 7, 8, #
**
** fftsize = 4, overlap = 0%, nsamples = 9
**
** bin		windowed samples
** ---		----------------
** 0		#, #, 0, 1
** 1		2, 3, 4, 5
** 2		6, 7, 8, #
*****/

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

	if (sono->psd == NULL)
	{
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
	if ((sono->method == METHOD_FFT) || (sono->method == METHOD_FFTW))
	{
		if (sono->ffttemppsd == NULL)
			sono->ffttemppsd = (float *)calloc(sono->fftsize / 2, sizeof(float));
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
		if (sono->method == METHOD_FFT)
		{
			calculate_fftpsd(sono, samples, nsamples, xoffset, sono->ffttemppsd);
		}
#ifdef HAVE_FFTW
		else if (sono->method == METHOD_FFTW)
		{
			calculate_fftwpsd(sono, samples, nsamples, xoffset, sono->ffttemppsd);
		}
#endif
		for (y=0; y < psdlen; y++)
			psd[y] = sono->ffttemppsd[sono->interp_index[y]];
	}
	else if (sono->method == METHOD_MULTITAPER)
	{
		if (sono->ffttemppsd == NULL)
			sono->ffttemppsd = (float *)calloc(sono->fftsize / 2, sizeof(float));
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
		calculate_mtmpsd(sono, samples, nsamples, xoffset, sono->ffttemppsd);
		for (y=0; y < psdlen; y++)
			psd[y] = sono->ffttemppsd[sono->interp_index[y]];
	}
	else if (sono->method == METHOD_MAXENTROPY)
	{
		calculate_mempsd(sono, samples, nsamples, xoffset, psdlen, psd, fmin, fmax);
	}
	return;
}

static void calculate_fftpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd)
{
	int i, fftsize;
	float *complexbuf, *ptr1, *ptr2, *outptr;

	complexbuf = sono->complexbuf;
	fftsize = sono->fftsize;

	if (sono->method != METHOD_FFT)
		return;

	if (sono->buftype == BUFTYPE_LINEAR)
		window_samples(complexbuf, sono->window, samples, nsamples, xoffset, fftsize, sono->dcoff);
	else
		window_samples_from_ringbuf(complexbuf, sono->window, samples, nsamples, xoffset, fftsize, sono->dcoff);

	/*
	** Finally, we call the FFT routine (which writes over the input)
	*/
	FFTROUTINE(complexbuf, fftsize);

	/*
	** Now compute the power spectrum
	*/
	outptr = psd;
	if (sono->linear == 0)
	{
		ptr1 = complexbuf;
		ptr2 = ptr1 + 1;
		*(outptr++) = 10.0 * log10f(MINFLOAT + (float)(SQR(*ptr1) / fftsize));
		ptr1 += 2;
		ptr2 += 2;
		for (i=1; i < sono->fftsize / 2; i++)
		{
			*(outptr++) = 10.0 * log10f((float)(MINFLOAT + (SQR(*ptr1) + SQR(*ptr2)) / fftsize));
			ptr1 += 2;
			ptr2 += 2;
		}
	}
	else
	{
		ptr1 = complexbuf;
		ptr2 = ptr1 + 1;
		*(outptr++) = ((float)(SQR(*ptr1) / fftsize));
		ptr1 += 2;
		ptr2 += 2;
		for (i=1; i < sono->fftsize / 2; i++)
		{
			*(outptr++) = ((float)((SQR(*ptr1) + SQR(*ptr2)) / fftsize));
			ptr1 += 2;
			ptr2 += 2;
		}
	}
}

#ifdef HAVE_FFTW
static void calculate_fftwpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd)
{
	int i, fftsize, N;
	float *complexbuf, *outptr;
	double *in_arr;

	if (sono->method != METHOD_FFTW)
		return;

	complexbuf = sono->complexbuf;
	fftsize = sono->fftsize;

	if (sono->buftype == BUFTYPE_LINEAR)
		window_samples(complexbuf, sono->window, samples, nsamples, xoffset, fftsize, sono->dcoff);
	else
		window_samples_from_ringbuf(complexbuf, sono->window, samples, nsamples, xoffset, fftsize, sono->dcoff);

	/*
	** Note that the plan must change if the fftsize changes, or if the arrays change
	** (which they don't)
	*/
	if ((sono->fftw_plan == NULL) || (fftsize != sono->fftw_plan_nsamp))
	{
		if (sono->fftw_plan) fftw_destroy_plan(sono->fftw_plan);
		sono->fftw_plan = fftw_plan_r2r_1d(fftsize, sono->fftw_data_in, sono->fftw_data_out, FFTW_R2HC, FFTW_ESTIMATE);
		sono->fftw_plan_nsamp = fftsize;
	}

	/*
	** Load the data into the input array
	*/
	for (i=0; i < fftsize; i++)
		sono->fftw_data_in[i] = complexbuf[i];

	/*
	** Now, call the FFT routine by running the plan.
	*/
	fftw_execute(sono->fftw_plan);

	/*
	** Now compute the power spectrum
	** Note: the PSD is normally sqrt(re^2 + im^2).
	** 
	** and, decibels are defined as 10 times the log10 of the power(i.e. watts) ratio, i.e.:
	** dB = 10 * log10( power1 / power2).
	** BUT, since our signals are proportional to volts rather than power,
	** (power proportional to volts squared),
	** dB = 20 * log10 (volt1 / volt2)
	**
	** Now, to get the Power spectral density in dB (take volt2 to be 1), we use:
	** dB = 20 * log10 ( sqrt(re^2 + im^2) )
	** which is equivalent to:
	** dB = 20 * log10 ( (re^2 + im^2) ^ (0.5) )
	** which is equivalent to:
	** dB = 20 * 0.5 * log10( (re^2 + im^2) )
	** i.e.
	** dB = 10 * log10( re^2 + im^2 )
	** which is what we do below.
	**
	** The observant will note that we first divide by fftsize.
	** This is just because the FFTW library returns an unnormalized result.
	*/
	N = fftsize;
	outptr = psd;
	in_arr = sono->fftw_data_out;
	*(outptr++) = (MINFLOAT + SQR(in_arr[0])) / (float)N;
	for (i=1; i < (N + 1) / 2; i++)
		*(outptr++) = (MINFLOAT + SQR(in_arr[i]) + SQR(in_arr[N - i])) / (float)N;
	if ((N & 1) == 0) // (if N is even)
		psd[N / 2] = (MINFLOAT + SQR(in_arr[N / 2])) / (float)N;
	if (sono->linear == 0)
	{
		for (i=0; i < (N + 1) / 2; i++)
			psd[i] = 10.0 * log10f(psd[i]);
	}
}
#endif

static void calculate_mempsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, int psdlen, float *psd, float fmin, float fmax)
{
	int i, fftsize;
	float *complexbuf, *outptr, *ptr1;

	complexbuf = sono->complexbuf;
	fftsize = sono->fftsize;

	if (sono->method != METHOD_MAXENTROPY)
		return;

	if (sono->ms == NULL)
		sono->ms = init_memspect(sono->fftsize, sono->npoles);

	/*
	** Fill in the FFT buffer with the samples to transform.
	** The first and last time bin are special-cased because
	** they may need to be zero-padded.
	**
	** NOTE: one of the whole points of the maximum entropy method is that
	** windowing isn't really needed.  But, there are references suggesting
	** that it can be helpful even with this method...  Anyway, we don't
	** really apply a window here - by passing NULL for the window, we
	** effectively use a rectangular/boxcar window regardless of the
	** setting.
	*/
	if (sono->buftype == BUFTYPE_LINEAR)
		window_samples(complexbuf, NULL, samples, nsamples, xoffset, fftsize, sono->dcoff);
	else
		window_samples_from_ringbuf(complexbuf, NULL, samples, nsamples, xoffset, fftsize, sono->dcoff);

	/*
	** The Maximum Entropy code pretty much throws up when it is
	** fed all-zeros.  It works, sort of, but it takes MUCH longer.
	** Thus, I bypass this case, since I already know what the
	** spectrum for all-zeros should look like.
	*/
	for (i=0; i < fftsize; i++)
		if (complexbuf[i] > 0.0000001)
			break;
	if (i == fftsize)
	{
		outptr = psd;
		for (i=0; i < psdlen; i++)
			*(outptr++) = 0.0;
		return;
	}

	mem_calc_coefficients(sono->ms, complexbuf);
	outptr = psd;

#ifdef notdef
	if ((ABS(fmin) > 0.001) || (ABS(fmax - 0.5) > 0.001))
	{
		float f, step;

		/* mem_eval_spectrum can't return 0 */
		step = (fmax - fmin) / psdlen;
		for (f=fmin, i=0; i < psdlen; i++, f += step)
		{
			if (sono->linear)
				*(outptr++) = mem_eval_spectrum(sono->ms, f);
			else
				*(outptr++) = 10.0 * log10f(mem_eval_spectrum(sono->ms, f));
		}
	}
	else
#endif
	{
		/* mem_eval_spectrum can't return 0 */
		mem_spectrum2(sono->ms, complexbuf, psdlen, fmin, fmax);
		ptr1 = complexbuf;
		if (sono->linear)
		{
			for (i=0; i < psdlen; i++)
				*(outptr++) = *(ptr1++);
		}
		else
		{
			for (i=0; i < psdlen; i++)
				*(outptr++) = 10.0 * log10f(*(ptr1++));
		}
	}
}

static void calculate_mtmpsd(struct sonogram *sono, short *samples, int nsamples, int xoffset, float *psd)
{
	int i, fftsize;
	float *complexbuf, *outptr;
	float norm;

	complexbuf = sono->complexbuf;
	fftsize = sono->fftsize;

	if (sono->method != METHOD_MULTITAPER)
		return;

	/*
	** Fill in the input buffer with the samples to transform.
	** The first and last time bin are special-cased because
	** they may need to be zero-padded.
	** Note: we're not really applying a window function here:
	** by passing NULL for the window, we effectively use a boxcar
	** or rectangular window.
	*/
	if (sono->buftype == BUFTYPE_LINEAR)
		window_samples(complexbuf, NULL, samples, nsamples, xoffset, fftsize, sono->dcoff);
	else
		window_samples_from_ringbuf(complexbuf, NULL, samples, nsamples, xoffset, fftsize, sono->dcoff);

#ifdef FORTRAN_MTM
	mtmpsd(complexbuf, fftsize, sono->mtm->p_nwin, sono->mtm->p_npi,
		0.0, 0.0, sono->mtm->p_inoise, sono->mtm->p_ilog,
		sono->mtm->p_fsmooth, sono->mtm->p_inorm, sono->mtm->p_ispec, sono->mtm->p_ithresh,
		sono->mtm->p_isignal);
	norm = fftsize;
#else
	my_multitaper_spectrum(complexbuf, fftsize, sono->mtm->p_nwin, sono->mtm->p_npi, complexbuf, 1 + sono->mtm->p_ispec);
	norm = 1.0;
#endif

	outptr = psd;
	if (sono->linear)
	{
		for (i=0; i < sono->fftsize / 2; i++)
			*(outptr++) = complexbuf[i] / (norm);
	}
	else
	{
		for (i=0; i < sono->fftsize / 2; i++)
			*(outptr++) = 10.0 * log10f(complexbuf[i] / (norm));
	}
}



/****************************************************************************************
********   WINDOW FUNCTIONS
****************************************************************************************/
static int init_window(int windowtype, float *window, int fftsize)
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
	return 0;
}

static void init_hamming_window(float *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * PI / fftsize;
	for (i=0; i < fftsize; i++)
		window[i] = 0.5*(1.0-(float)cos(K*i));
}

static void init_hanning_window(float *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * PI / fftsize;
	for (i=0; i < fftsize; i++)
		window[i] = 0.5*(1.0-(float)cos(K*i));
}

static void init_boxcar_window(float *window, int fftsize)
{
	int i;

	for (i=0; i < fftsize; i++)
		window[i] = 1.0;
}

static void init_triang_window(float *window, int fftsize)
{
	int i, half_fftsize;

	half_fftsize = fftsize / 2;
	if ((fftsize % 2) == 0)
	{
		for (i=0; i < half_fftsize; i++)
			window[i] = 1.0 * (i / ((float)(half_fftsize - 1)));
	}
	else
	{
		for (i=0; i <= half_fftsize; i++)
			window[i] = 1.0 * (i / ((float)half_fftsize));
	}
	for (; i < fftsize; i++)
		window[i] = window[fftsize - i - 1];
}

static void init_blackman_window(float *window, int fftsize)
{
	int i;
	double K;

	K = 2.0 * PI / (fftsize - 1);
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


/*
** Applies a window function to the samples in the linear buffer
** starting at 'samples'.  The number of samples to be windowed
** is 'fftsize'.  These samples will be centered at offset 'offset'
** from the location pointed to by 'samples'.
**
** Here, we adjust so the window is CENTERED on the sample of interest.
** Thus, we subtract half the window size from the sample offset.
** 
** If the new offset is out of bounds, treat as sample of zero value.
*/
static void window_samples(float *dest, float *window, short *samples, int nsamples, int offset, int fftsize, float dcoff)
{
	int i, indx;

	indx = offset - (fftsize / 2);
	for (i=0; i < fftsize; i++)
	{
		if (indx + i < 0)
			*(dest++) = 0.0;
		else if (indx + i >= nsamples)
			*(dest++) = 0.0;
		else if (window != NULL)
			*(dest++) = (float)window[i] * ((float)samples[indx + i] - dcoff);
		else if (window == NULL)
			*(dest++) = ((float)samples[indx + i] - dcoff);
	}
}

/*
** Applies a window function to the samples in the ring buffer
** starting at 'samples'.  The number of samples to be windowed
** is 'fftsize'.  These samples will be centered at offset 'offset'
** from the location pointed to by 'samples'.
**
** NOTE: we adjust so the window is CENTERED on the sample of interest.
** Thus, we subtract half the window size from the sample offset.
*/
static void window_samples_from_ringbuf(float *dest, float *window, short *samples, int nsamples, int offset, int fftsize, float dcoff)
{
	int i, indx;

	indx = (offset - (fftsize / 2) + nsamples) % nsamples;
	if (indx + fftsize < nsamples)
	{
		for (i=0; i < fftsize; i++)
		{
			if (window != NULL)
				*(dest++) = (float)window[i] * ((float)samples[indx + i] - dcoff);
			else
				*(dest++) = ((float)samples[indx + i] - dcoff);
		}
	}
	else
	{
		for (i=0; i < fftsize; i++)
		{
			if (window != NULL)
				*(dest++) = (float)window[i] * ((float)samples[(indx + i) % nsamples] - dcoff);
			else
				*(dest++) = ((float)samples[(indx + i) % nsamples] - dcoff);
		}
	}
}


static multitaper *alloc_mtmstruct(void)
{
	multitaper *mt;

	mt = (multitaper *)malloc(sizeof(multitaper));
	if (mt != NULL)
	{
		mt->p_nwin = -1;
		mt->p_npi = -1;
		mt->p_f1 = -1.0;
		mt->p_f2 = -1.0;
		mt->p_inoise = -1;
		mt->p_ilog = -1;
		mt->p_fsmooth = -1;
		mt->p_inorm = -1;
		mt->p_ispec = -1;
		mt->p_ithresh = -1;
		mt->p_isignal = -1;
	}
	return mt;
}

