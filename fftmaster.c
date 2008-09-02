
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define ERROR_OUTOFMEMORY 2

#define MAXWINDOWSIZE 1024

#ifdef SGICOMPLIB

#ifdef AWROO
#define sfft1dui sfft1di
#endif

#define REALFFT
typedef float FFTREALSIZE;
static FFTREALSIZE workspace[MAXWINDOWSIZE + 15];
static int inittedsize = 0;
#define FFTINITROUTINE(fftsize)				if (inittedsize != fftsize) { sfft1dui(fftsize, workspace); inittedsize = fftsize; }
#define FFTROUTINE(complexbuf, fftsize)		sfft1du(-1, fftsize, complexbuf, 1, workspace);

#else

#define REALFFT
typedef float FFTREALSIZE;
#define FFTINITROUTINE(fftsize)
extern void rfft(FFTREALSIZE x[], int length);
#define FFTROUTINE(complexbuf, fftsize)		rfft(complexbuf, fftsize);

#endif

FFTREALSIZE hammingwindow[MAXWINDOWSIZE];
static int hammingwindowsize = 0;

void init_hamming_window(int fftsize);

#if defined(sun) || defined(linux)
#define log10f(num) ((float)log10((float)(num)))
#endif

long outfft(short *databuf, unsigned long databufsize, long fftsize, long overlap, float **spec_p, long *numtimebins_p)
{
	long status;
	int i, lower, frame;
	short *dataptr;
	long numtimebins, numrows;
	float *spec = NULL, inc, *outptr;
	FFTREALSIZE *ptr1, *ptr2, complexbuf[2 * MAXWINDOWSIZE];

	*spec_p = NULL;
	*numtimebins_p = 0;
	inc = (100.0 - overlap) / 100.0;
	numrows = fftsize / 2;
	numtimebins = databufsize / (fftsize * inc);

	FFTINITROUTINE(fftsize);

	while (numtimebins * (fftsize * inc) + fftsize > databufsize) numtimebins--;
	if (numtimebins == 0) numtimebins = 1;

	if ((spec = (float *)calloc(numtimebins * numrows, sizeof(float))) != NULL)
	{
		init_hamming_window(fftsize);
		for (frame=0; frame < numtimebins; frame++)
		{
			/* Get the appropriate frame, window it, and put it into the fft buffer */
			lower = (int)floor(fftsize * inc * frame);
			dataptr = &(databuf[lower]);
			ptr1 = complexbuf;
			for (i=0; i < fftsize; i++)
			{
#ifdef REALFFT
				*(ptr1++) = (FFTREALSIZE)hammingwindow[i] * *(dataptr++);
#else
				*(ptr1++) = (FFTREALSIZE)hammingwindow[i] * *(dataptr++);
				*(ptr1++) = (FFTREALSIZE)0.0;
#endif
			}

			/* Finally, the FFT (which writes over the input) */
			FFTROUTINE(complexbuf, fftsize);

			/* Power spectrum (writes into spec) */
			outptr = spec + frame*numrows;
			ptr1 = complexbuf;
			ptr2 = ptr1 + 1;
			for (i=0; i < numrows; i++)
			{
				*outptr = 10.0 * log10f((float)(*ptr1 * *ptr1 + *ptr2 * *ptr2));
				outptr++;
				ptr1 += 2;
				ptr2 += 2;
			}
		}
		*numtimebins_p = numtimebins;
		*spec_p = spec;
		status = SUCCESS;
	}
	else status = ERROR_OUTOFMEMORY;

	if ((status != SUCCESS) && (spec != NULL))
		free(spec);
	return status;
}


/*
 * Setup Hamming window
 */
#ifndef PI
#define PI      3.1415926537
#endif

void init_hamming_window(int fftsize)
{
	int i, half_fftsize;
	extern int hammingwindowsize;

	if (hammingwindowsize == fftsize)
		return;
	half_fftsize = fftsize/2;
	for (i=0; i < fftsize; i++)
		hammingwindow[i] = 0.5*(1.0-(FFTREALSIZE)cos(PI*i/(FFTREALSIZE)half_fftsize));
	hammingwindowsize = fftsize;
}


