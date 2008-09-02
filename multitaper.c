/* GRASP: Copyright 1997,1998  Bruce Allen */
/* (C) Copr. 1986-92 Numerical Recipes Software #.3. */

/* 
** INCLUDES
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


/*
** DEFINES
*/
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define DIAG1 0


/*
** GLOBALS
*/
/* range to search for fvalues in */
int imin,imax;


/*
** PROTOTYPES
*/
void realft(float*, unsigned long, int);
int get_pow_2(int inum);
int adwait(double *sqr_spec, double *dcf, double *el, int nwin, int num_freq, double *ares, double *degf, double avar);
void get_F_values(double *sri, int nf, int nwin, float *Fvalue, double *b, float *cest, int imin, int imax);
int hires(double *sqr_spec, double *el, int nwin, int num_freq, double *ares);
int jtridib(int *n, double *eps1, double *d, double *e, double *e2, double *lb, double *ub, int *m11, int *m, double *w, int *ind, int *ierr, double *rv4, double *rv5);
int jtinvit(int *nm, int *n, double *d, double *e, double *e2, int *m, double *w, int *ind, double *z, int *ierr, double *rv1, double *rv2, double *rv3, double *rv4, double *rv6);
void mt_get_spec(float *series, int inum, int klength, float *amp);
void zero_pad(float output[], int start, int olength);
double remove_mean(float x[], int lx);
double fcrit(int npoints, int nwin);
int is_power_of_two(int test);
int empty_window(float *workspace, int nbins, int index, int wwidth);
void get_F_values2(double *sr, double *si, int nf, int nwin, float *Fvalue, double *b, float *cest);
void slepian_tapers(int num_points, int nwin, double *lam, float nwdt, double *tapers, double *tapsum);

#ifdef LINES
void sort_lines(struct removed_lines *line_list, int *num_removed);
#endif


/***************************************************************/
/* This is the critical value of f for npoints number of points in a
   data-set and nwin tapers, according to Thompson's rule of thumb
   (Percival & Walden paragraph at top of pg 513) that the confidence
   should exceed 1-1/n, combined with the explicit formula for
   the appropriate percentage point given on pg 501.
*/
/***************************************************************/
double fcrit(int npoints, int nwin)
{
	double crit;
	double nm1;

	nm1 = nwin - 1.0;
	crit = nm1 * (pow((double)npoints, 1.0 / nm1) - 1.0);
	return crit;
}


/***************************************************************/
/* Checks if the argument is a power of 2 less than 2^32   */
/***************************************************************/
int is_power_of_two(int test)
{
	int i;
	unsigned int k=1;

	for (i=0; i<32; i++)
	{
		if (k==test)
			return 1;
		k+=k;
	}
	return 0;
}


#ifdef LINES
/***************************************************************/
/* Comparison function for ordering structures of type 
   removed_lines on the basis of their F-test values       */
/***************************************************************/
int fvalue_cmp(const void *l1, const void *l2)
{
	float fl1,fl2;

	fl1 = ((struct removed_lines *)l1)->fvalue;
	fl2 = ((struct removed_lines *)l2)->fvalue;
	if (fl1 > fl2)
		return -1;
	else if (fl1 < fl2)
		return 1;
	else
		return 0;
}


/***************************************************************/
/* Comparison function for ordering structures of type 
   removed_lines on the basis of their frequency bin position      */
/***************************************************************/
int index_cmp(const void *l1, const void *l2)
{
	int   il1,il2;

	il1 = ((struct removed_lines *)l1)->index;
	il2 = ((struct removed_lines *)l2)->index;
	if (il1 < il2)
		return -1;
	else if (il1 > il2)
		return 1;
	else
		return 0;
}
#endif


/***************************************************************/
/* A routine to check if there are any other nearby points arising 
   from the same spectral line - if so we only want to remove the
   most significant   */
/***************************************************************/
int empty_window(float *workspace, int nbins, int index, int wwidth)
{
	int i, bin, test = 1;

	/* step through the +- wwidth */
	for (i=-wwidth; i<=wwidth; i++)
	{
		/* if the frequency bin is below 0 or above the end, set correctly */
		bin = index+i;
		if (bin<0)
			bin = 0;
		if (bin>=nbins)
			bin = nbins-1;
		/* test if there are any non-zero entries in workspace */
		if ((workspace[2*bin] != 0.0) || (workspace[2*bin+1] != 0.0))
		{
			test = 0;
			break;
		}
	}
	return test;
}


#ifdef LINES
/**************************************************************/
/* this routine sorts into increasing frequency order and combines any duplicates */
/**************************************************************/
void sort_lines(struct removed_lines *line_list, int *num_removed)
{
	int top, bot, did_remove;

	/* Now sort into increasing frequency index order */
	qsort((void *)line_list, (size_t)*num_removed, (size_t)sizeof(struct removed_lines), index_cmp); 

	/* search/combine duplicates */
	top = 0;
	did_remove = 0;
	while (top<*num_removed)
	{
		bot = top+1;
		while (bot < *num_removed && line_list[bot].index == line_list[top].index)
		{
			line_list[top].re += line_list[bot].re;
			line_list[top].im += line_list[bot].im;
			line_list[bot++].fvalue = 0.0;
			did_remove++;
		}
		top = bot;
	}

	/* sort into value by fvalues -- puts any empty structures at the end of the list */
	qsort((void *)line_list, (size_t)*num_removed, (size_t)sizeof(struct removed_lines), fvalue_cmp); 
	*num_removed -= did_remove;

	/* sort into increasing frequency index order */
	qsort((void *)line_list, (size_t)*num_removed, (size_t)sizeof(struct removed_lines), index_cmp);
	return;
}


/***************************************************************
***************************************************************/
void remove_spectral_lines(float *data, int npoints, int padded_length, float nwdt,
	int nwin, int max_lines, int maxpass, int *num_removed, struct removed_lines *line_list,
	float *mtap_spec_init, float *mtap_spec_final, int outspec, int fimin, int fimax)
{
	int num_freq, i, nlines = 0, bin, wwidth, passes = 0, near_nyquist, did_remove = 0;
	float *dof, *fvalues, *workspace, *cestimate, *spec, fcritical;
	struct removed_lines *work_list;
	void clear(float *array, int n, int spacing);

	/* upper and lower frequency bin limits for F-test (expensive!) */
	imin = fimin;
	imax = fimax;

	/* check that the padded array has a length = integer power of 2 */
	if (!is_power_of_two(padded_length))
	{
		fprintf(stderr, "remove_spectral_lines() %s:%s",__FILE__,__LINE__);
		fprintf(stderr, "Error: the padded_length must");
		fprintf(stderr, "be an integer power of 2 <= 2^31, but is %d \n",padded_length);
		exit(-1);
	}

	/* The total number of frequencies, including DC and Nyquist */
	num_freq = 1+padded_length/2;

	/* The number of degrees of freedom (only useful for adaptive, not high_res tapered spectra */
	dof = (float *)malloc(sizeof(float)*num_freq);

	/* Array of values of the F-test */
	fvalues = (float *)malloc(sizeof(float)*num_freq);

	/* Array to hold the real and imaginary parts of the associated Fourier coefficient */
	cestimate = (float *)malloc(sizeof(float)*2*num_freq); 

	/* Array that will be used to do IFFT to remove the lines */
	workspace = (float *)malloc(sizeof(float)*padded_length);

	/* Working structure to hold frequencies and f-values  */
	work_list = (struct removed_lines *)malloc(sizeof(struct removed_lines)*num_freq);

	if ((dof==NULL) || (fvalues==NULL) || (cestimate==NULL) || (workspace==NULL) || (work_list==NULL))
	{
		fprintf(stderr, "remove_spectral_lines() %s:%d",__FILE__,__LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	/* The critical value of F, below which we do not consider the line component significant */
	/* See Percival and Walden Sec 10.11 and pg 513 (top) for discussion of critical value */
	fcritical = fcrit(npoints, nwin);
	
	/* compute number of bins over which F-test is correlated */
	wwidth = (int)(nwdt*padded_length/npoints);

	*num_removed = 0;
        
	/* run through loop removing lines, ending when no additional lines removed */
	while ((passes == 0) || ((passes < wwidth) && (nlines != 0) && (passes < maxpass)))
	{
		/* subtract the mean value from the data set */
		remove_mean(data,npoints);

		/* if first pass through store spectrum in mtap_spec_init otherwise  */
		/* store it in mtap_spec_final  */
		if (passes++ == 0)
			spec = mtap_spec_init;
		else
			spec = mtap_spec_final;

		/* calculate the multitaper spectrum */
		multitaper_spectrum(data,
			npoints,
			1 /* highres */,
			nwin,
			nwdt,
			4 /* normalization choice */,
			1.0 /* srate */,
			spec,
			dof,
			fvalues,
			padded_length,
			cestimate,
			outspec);

		/* (skip DC, Nyquist,Nyquist-1) looking for frequency bins which pass the F-test */
		nlines = 0;
		near_nyquist = (int)((float)padded_length/(float)npoints);
		for (i=1; i < num_freq-near_nyquist-1; i++)
		{
			if (fvalues[i] > fcritical)
			{
				work_list[nlines].index = i;
				work_list[nlines].fvalue = fvalues[i];
				work_list[nlines].re = cestimate[2*i];
				work_list[nlines].im = cestimate[2*i+1];

				/* pass back estimates of the derivative (f+ - f-)/2 */
				work_list[nlines].dcdbr = 0.5*(cestimate[2*i+2]-cestimate[2*i-2]);
				work_list[nlines].dcdbi = 0.5*(cestimate[2*i+3]-cestimate[2*i-1]);

				/* pass back estimates of the second derivative (f+ - 2f + f-) */
				work_list[nlines].d2cdb2r = cestimate[2*i+2]-2.0*cestimate[2*i]+cestimate[2*i-2];
				work_list[nlines].d2cdb2i = cestimate[2*i+3]-2.0*cestimate[2*i+1]+cestimate[2*i-1];
				nlines++;	
			}		  
		}

		/* Now sort into order of decreasing significance */
		qsort((void *)work_list, (size_t)nlines, (size_t)sizeof(struct removed_lines), fvalue_cmp); 

		/* clear the workspace to hold nonzero value to pass to realft  */
		clear(workspace, padded_length, 1);
		did_remove = 0;

		for (i=0; (i < nlines) && (*num_removed < max_lines); i++)
		{
			/* check to see if we are removing any lines within wwidth bins */
			if (empty_window(workspace, padded_length/2, bin = work_list[i].index,wwidth))
			{
				/* if there are no nearby lines already removed, add lines to remove list */
				line_list[*num_removed] = work_list[i];
				workspace[2*bin] = -2.0*line_list[*num_removed].re;
				workspace[2*bin+1] = -2.0*line_list[*num_removed].im;
				(*num_removed)++;
				did_remove = 1;
			}			
		}

		/* efficiently calculate the sum of the spectral lines to remove */
		/* the conventions are almost exactly Percival and Walden example eqns 20 lines down */
		/* on page 513 but our imaginary part has opposite sign */
		/* because P&W FFT conventions eqn (65ab) are opposite to Numerical Recipes */
		if (did_remove)
		{
			/* subtract the spectral lines */
			realft(workspace-1,padded_length,-1);
			for (i=0;i<npoints;i++)
				data[i]+=workspace[i];
		} 
	}

	/*
	** If we removed no points, final spectrum = initial spectrum
	*/
	if (*num_removed==0) for (i=0;i<num_freq;i++) mtap_spec_final[i] = mtap_spec_init[i];

	/*
	** If last pass through loop changed data, need to find spectrum
	** This will only happen if the number of lines removed is maxlines.  However it
	** may be that the we remove maxlines but didremove==0
	*/
	if (did_remove!=0 && outspec) 
	{
		multitaper_spectrum(data,
			npoints,
			1 /* highres */,
			nwin,
			nwdt,
			4 /* normalization choice */,
			1.0 /* srate */,
			mtap_spec_final,
			dof,
			fvalues,
			padded_length,
			cestimate,
			outspec);
	}

	/* combine any duplicates and sort into increasing frequency bin order */
	sort_lines(line_list,num_removed);

	/* de-allocate the storage space */
	free(dof);
	free(fvalues);
	free(workspace);
	free(work_list);
	free(cestimate);  
	return;
}
#endif

/*--------------------------------------------------------*/
/*----------------mt_get_spec---------------------------*/
	/*
	series=input time series
	inum=length of time series
	klength= number of elements in power spectrum (a power of 2)
	amp=returned power spectrum
	*/
/*--------------------------------------------------------*/
void mt_get_spec(float *series, int inum, int klength, float *amp)
{
	int i, isign=1;

	/* copy amp onto series and apply zero padding to  klength */
	for (i=0; i < inum; i++)
		amp[i] = series[i];

	zero_pad(amp, inum, klength);

	/*
	** Fast Fourier Transform Routine:  here we are using the Numerical
	** Recipes routine realft which returns the fft in the 1-D input
	** array packed as pairs of real numbers. The realft routine
	** requires the input array to start at index=1 so we must decrement
	** the index of amp
	*/
	realft(amp-1, (unsigned long)klength, isign);
	return;
}


/*------------------------------------------------------------*/
/*----------------multitaper_spectrum---------------------------*/
	/*
	data=floating point input time series
	npoints=number of points in data
	kind=flag for choosing hires (1) or adaptive (2) weighting coefficients 
	nwin=number of taper windows to calculate
	nwdt=order of the slepian functions
	inorm=flag for choice of normalization
	dt=sampling interval (time)
	ospec=output spectrum
	dof=degrees of freedom at each frequency
	Fvalues=Ftest value at each frequency estimate
	klen=number of frequecies calculated (power of 2)
	cest is the estimated complex coefficient of a complex exponential (array of length klen+2)
	*/

	/*
	lambda=vector of eigenvalues
	tapsum=sum of each taper, saved for use in adaptive weighting
	tapers= matrix of slepian tapers packed in a 1D double array
	 */
/*------------------------------------------------------------*/
void multitaper_spectrum(float *data, int npoints, int kind, int nwin, float nwdt, int inorm,
	float dt, float *ospec, float *dof, float *Fvalues, int klen, float *cest, int outspec)
{
	int             i, j;
	float          *b;
	int             iwin, kk;
	double          anrm, norm, sqramp, sum;
	double         *ReImSpec,*sqr_spec, *amu;
	float          *amp;
	int             num_freqs, len_taps, num_freq_tap,n1, n2, kf;
	double         *dcf, *degf, avar;
	/*************/
	static int save_npoints = 0,save_nwin = 0;
	static float save_nwdt = 0.0;
	static double   *lambda=NULL,*tapers=NULL,*tapsum=NULL;

	if ((npoints != save_npoints) || (nwin != save_nwin) || (nwdt != save_nwdt))
	{
		/* This means that we must compute tapers, etc */
		len_taps = npoints*nwin;
		lambda = (double *)realloc((void *)lambda,(size_t)nwin*sizeof(double));
		tapsum = (double *)realloc((void *)tapsum,(size_t)nwin*sizeof(double));
		tapers = (double *)realloc((void *)tapers,(size_t)len_taps*sizeof(double));
		if ((lambda == NULL) || (tapsum == NULL) || (tapers == NULL))
		{
			fprintf(stderr, "multitaper_spectrum() %s:%d",__FILE__,__LINE__);
			fprintf(stderr, "Insufficient memory - realloc returned NULL.\n");
			exit(-1);
		}
		save_npoints = npoints;
		save_nwin = nwin;
		save_nwdt = nwdt;

		/* get the set of Slepian tapers  */
		slepian_tapers(npoints, nwin, lambda, nwdt, tapers, tapsum);
	}

	/* choose normalization based on inorm flag  */
	anrm = 1.0;
	switch (inorm)
	{
		case 1:
			anrm = npoints;
			break;
		case 2:
			anrm = 1.0/dt;
			break;
		case 3:
			anrm = sqrt((double) npoints);
			break;
		default:
			anrm = 1.0;
			break;
	}

	num_freqs = 1 + klen/2;
	num_freq_tap = num_freqs * nwin;

	b = (float *)malloc((size_t)npoints * sizeof(float));
	amu = (double *)malloc((size_t)num_freqs * sizeof(double));
	sqr_spec = (double *)malloc((size_t)num_freq_tap * sizeof(double));
	ReImSpec = (double *)malloc(2 * (size_t)num_freq_tap * sizeof(double));
	amp = (float *)malloc((size_t)klen * sizeof(float));
	if ((b==NULL) || (amu==NULL) || (sqr_spec==NULL) || (ReImSpec==NULL) || (amp==NULL))
	{
		fprintf(stderr, "multitaper_spectrum() %s:%d",__FILE__,__LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	/* apply the taper in the loop.  do this nwin times  */
	for (iwin=0; iwin < nwin; iwin++)
	{
		kk = iwin*npoints;
		kf = iwin*num_freqs;

		/* application of iwin-th taper   */
		for (j=0; j < npoints; j++)
			b[j] = data[j] * tapers[kk+j];

		/* calculate the eigenspectrum */
		mt_get_spec(b, npoints, klen, amp);					

		/* check that indices are within bounds */
		if ((2*num_freqs-3 > klen) || (num_freqs-1+kf > num_freq_tap))
		{
			fprintf(stderr, "multitaper_spectrum() %s:%d ",__FILE__,__LINE__);
			fprintf(stderr, "error in index\n");
		}

		if (outspec)
		{
			/* get spectrum from real fourier transform    */
			norm = 1.0/(anrm*anrm);
			sum = 0.0;
			for (i=1; i<num_freqs-1; i++)
			{
				sqramp = SQR(amp[2*i+1])+SQR(amp[2*i]);
				ReImSpec[2*(i+kf)] = amp[2*i];
				ReImSpec[2*(i+kf)+1] = amp[2*i+1];
				sqr_spec[i+kf] = norm*sqramp;
				sum += sqramp;
			}

			/* do DC and Nyquist */
			sqr_spec[kf] = norm*SQR(amp[0]);
			sqr_spec[num_freqs-1+kf] = norm*SQR(amp[1]);
			ReImSpec[2*kf] = amp[0];
			ReImSpec[2*kf+1] = 0.0;
			ReImSpec[2*(num_freqs-1+kf)] = amp[1];
			ReImSpec[2*(num_freqs-1+kf)+1] = 0.0;
			sum += sqr_spec[kf]+sqr_spec[num_freqs-1+kf];
		}
		else
		{
			/* get spectrum from real fourier transform    */
			for (i=2; i<2*num_freqs-2; i++)
			{
				ReImSpec[i+2*kf] = amp[i];
			}

			/* do DC and Nyquist */
			ReImSpec[2*kf] = amp[0];
			ReImSpec[2*kf+1] = 0.0;
			ReImSpec[2*(num_freqs-1+kf)] = amp[1];
			ReImSpec[2*(num_freqs-1+kf)+1] = 0.0;
		}

	}
	free(amp);
	free(b);

	/* choice of hi-res or adaptive weighting for spectra    */
	switch (kind)
	{
		case 1:

			if (outspec)
				hires(sqr_spec, lambda, nwin, num_freqs, amu);
			get_F_values(ReImSpec, num_freqs, nwin, Fvalues, tapsum, cest, imin, imax);

			if (outspec)
			{
				for (i=0; i < num_freqs; i++)
				{
					ospec[i] = amu[i];
					dof[i] = nwin-1;
				}
			}
			else
			{
				for (i=0; i < num_freqs; i++)
					dof[i] = nwin-1;
			}
			break;

		case 2:
			/* get avar=variance */
			n1=0;
			n2=npoints;
			avar=0.0;

			for (i=n1; i < n2; i++)
				avar += (data[i])*(data[i]);
			switch (inorm)
			{
				case 1:
					avar = (avar/npoints)/npoints;
					break;
				case 2:
					avar = avar*dt*dt;
					break;
				case 3:
					avar = avar/npoints;
					break;
				default:
					break;
			}

			dcf = (double *) malloc((size_t) num_freq_tap*sizeof(double));
			degf = (double *) malloc((size_t) num_freqs*sizeof(double));
			if ((dcf==NULL) || (degf==NULL))
			{
				fprintf(stderr, "multitaper_spectrum() %s:%d ",__FILE__,__LINE__);
				fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
				exit(-1);
			}

			adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);
			get_F_values(ReImSpec, num_freqs, nwin, Fvalues, tapsum,cest,imin,imax);

			for (i=0; i < num_freqs; i++)
			{
				ospec[i] = amu[i];
				dof[i] = degf[i];
			}

			free(dcf);
			free(degf);
			break;

		default:
			fprintf(stderr, "multitaper_spectrum() %s:%d ",__FILE__,__LINE__);
			fprintf(stderr, "Do not recognize kind = %d\n",kind);
			exit(-1);
			break;
	}

	/* free up memory and return  */

	free(amu);
	free(sqr_spec);
	free(ReImSpec);
}


/*------------------------------------------------------------*/
/*----------------multitaper_cross_spectrum---------------------------*/
  /*
    data1=floating point input time series #1
    data2=floating point input time series #2
    npoints=number of points in data1 (and data2)
    npadded=padded length (integer power of 2)
    dt=sampling interval (time)
    nwin=number of taper windows to calculate
    nwdt=order of the slepian functions
    ReImSpec12=complex-valued cross-correlation spectrum
  */

  /*
    lambda=vector of eigenvalues
    tapsum=sum of each taper, for adaptive weighting (not used here!!).
    tapers=matrix of slepian tapers packed in a 1D double array
  */
/*------------------------------------------------------------*/
void multitaper_cross_spectrum(float *data1, float *data2, int npoints, int npadded, float dt, int nwin, float nwdt, double *ReImSpec12)
{
	int             i, j;
	int             iwin, kk, kf;
	int             num_freqs, len_taps, num_freq_tap;
	float           a;
	float           *tap_data1, *tap_data2;
	float           *amp1, *amp2;
	double          *ReImSpec1, *ReImSpec2;
	static int      save_npoints=0, save_nwin=0;
	static float    save_nwdt=0.0;
	static double   *lambda=NULL, *tapers=NULL, *tapsum=NULL;

	if ((npoints != save_npoints) || (nwin != save_nwin) || (nwdt != save_nwdt))
	{
		/* This means that we must compute tapers, etc */
		len_taps = npoints*nwin;
		lambda = (double *)realloc((void *)lambda,(size_t)nwin*sizeof(double));
		tapsum = (double *)realloc((void *)tapsum,(size_t)nwin*sizeof(double));
		tapers = (double *)realloc((void *)tapers,(size_t)len_taps*sizeof(double));
		if (lambda==NULL || tapsum==NULL || tapers==NULL)
		{
			fprintf(stderr, "multitaper_cross_spectrum() %s:%d ", __FILE__, __LINE__);
			fprintf(stderr, "Insufficient memory - realloc returned NULL.\n");
			exit(-1);
		}
		save_npoints = npoints;
		save_nwin = nwin;
		save_nwdt = nwdt;
		/* get the set of Slepian tapers  */
		slepian_tapers(npoints, nwin, lambda, nwdt, tapers, tapsum);
	}

	num_freqs = 1+npadded/2;
	num_freq_tap = num_freqs*nwin;

	/* allocate memory */
	tap_data1 = (float *)malloc((size_t)npoints*sizeof(float));
	tap_data2 = (float *)malloc((size_t)npoints*sizeof(float));
	amp1 = (float *)malloc((size_t)npadded*sizeof(float));
	amp2 = (float *)malloc((size_t)npadded*sizeof(float));
	ReImSpec1 = (double *)malloc(2*(size_t)num_freq_tap*sizeof(double));
	ReImSpec2 = (double *)malloc(2*(size_t)num_freq_tap*sizeof(double));
	if (tap_data1==NULL || tap_data2==NULL || amp1==NULL || amp2==NULL || ReImSpec1==NULL || ReImSpec2==NULL)
	{
		fprintf(stderr, "multitaper_cross_spectrum() %s:%d ", __FILE__, __LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	/* apply the taper in the loop.  do this nwin times  */
	for (iwin=0; iwin<nwin; iwin++)
	{
		kk = iwin*npoints;
		kf = iwin*num_freqs;

		/* application of iwin-th taper   */
		for (j=0; j<npoints; j++)
		{
			tap_data1[j] = data1[j]*tapers[kk+j];
			tap_data2[j] = data2[j]*tapers[kk+j];
		}

		/* calculate the eigenspectrum */
		mt_get_spec(tap_data1, npoints, npadded, amp1);
		mt_get_spec(tap_data2, npoints, npadded, amp2);

		/* get spectrum from real fourier transform */
		for (i=1; i<num_freqs-1; i++)
		{
			ReImSpec1[2*(i+kf)] = amp1[2*i];
			ReImSpec1[2*(i+kf)+1] = amp1[2*i+1];

			ReImSpec2[2*(i+kf)] = amp2[2*i];
			ReImSpec2[2*(i+kf)+1] = amp2[2*i+1];
		}

		/* do DC and Nyquist */
		ReImSpec1[2*kf] = amp1[0];
		ReImSpec1[2*kf+1] = 0.0;
		ReImSpec1[2*(num_freqs-1+kf)] = amp1[1];
		ReImSpec1[2*(num_freqs-1+kf)+1] = 0.0;

		ReImSpec2[2*kf] = amp2[0];
		ReImSpec2[2*kf+1] = 0.0;
		ReImSpec2[2*(num_freqs-1+kf)] = amp2[1];
		ReImSpec2[2*(num_freqs-1+kf)+1] = 0.0;
	}

	/* calculate high resolution cross-correlation spectrum */
	for (i=0; i<num_freqs; i++)
	{
		ReImSpec12[2*i] = 0.0;
		ReImSpec12[2*i+1] = 0.0;
	}

	for (iwin=0; iwin<nwin; iwin++)
	{
		a = dt*dt/(lambda[iwin]*nwin);
		kf = iwin*num_freqs;
		for (i=0; i<num_freqs; i++)
		{
			/* real part of cross-correlation spectrum */
			ReImSpec12[2*i] = ReImSpec12[2*i] + a*(ReImSpec1[2*(i+kf)]*ReImSpec2[2*(i+kf)]+ ReImSpec1[2*(i+kf)+1]*ReImSpec2[2*(i+kf)+1]);
			/* imag part of cross-correlation spectrum */
			ReImSpec12[2*i+1] = ReImSpec12[2*i+1] + a*(ReImSpec1[2*(i+kf)]*ReImSpec2[2*(i+kf)+1]- ReImSpec1[2*(i+kf)+1]*ReImSpec2[2*(i+kf)]);
		}
	}

	/* free up memory and return  */
	free(tap_data1);
	free(tap_data2);
	free(amp1);
	free(amp2);
	free(ReImSpec1);
	free(ReImSpec2);
	return;
}


/*--------------------------------------------------------*/
/*--------------------slepian_tapers--------------------------*/
	/*
	get the multitaper slepian functions.
	num_points=number of points in data stream: N in Percival and Walden
	nwin=number of window tapers: K in Percival and Walden (taper sums go k=0 to K-1).
		Should have K<2 nwdt.
	lam= vector of eigenvalues
	nwdt=quantity called NWdt in Percival and Walden: time x freq resolution bandwidth
	tapsum=sum of each taper, saved for use in adaptive weighting
	tapers= matrix of slepian tapers, packed in a 1D double array
	 */
/*--------------------------------------------------------*/
void slepian_tapers(int num_points, int nwin, double *lam, float nwdt, double *tapers, double *tapsum)
{
	int             i, k, kk;
	double          ww, cs, ai, an, eps, rlu, rlb, aa;
	double          dfac, drat, gamma, bh, tapsq, TWOPI, DPI;
	double         *diag, *offdiag, *offsq;
	double         *scratch1, *scratch2, *scratch3, *scratch4, *scratch5;

	/* need to initialize iwflag=0 */
	double          anwdt;
	double         *ell;
	int            *ip;
	double         *evecs;
	long            len;
	int             ierr;
	int             m11;
	DPI = (double) M_PI;
	TWOPI = (double) 2 *DPI;

	/* check that the number of windows is not too large */
	if (nwin>=2.0*nwdt)
	{
		fprintf(stderr, "slepian_tapers() %s:%d ", __FILE__, __LINE__);
		fprintf(stderr, "Warning: number of taper ");
		fprintf(stderr, "windows %d should be LESS than 2 x NWdt = %f\n",nwin,2.0*nwdt);
		fprintf(stderr, "Please see Percival and Walden pg 335.\n");
	}

	anwdt = nwdt;
	an = (double)(num_points);
	ww = (double)(anwdt)/an;	/* this corresponds to P&W's W value  */
	cs = cos(TWOPI*ww);
	len = num_points*nwin;

	ell = (double *)malloc((size_t)nwin*sizeof(double));
	diag = (double *)malloc((size_t)num_points*sizeof(double));
	offdiag = (double *)malloc((size_t)num_points*sizeof(double));
	offsq = (double *)malloc((size_t)num_points*sizeof(double));
	evecs = (double *)malloc((size_t)len*sizeof(double));

	ip = (int *)malloc((size_t)nwin*sizeof(int));
	scratch1 = (double *)malloc((size_t)num_points*sizeof(double));
	scratch2 = (double *)malloc((size_t)num_points*sizeof(double));
	scratch3 = (double *)malloc((size_t)num_points*sizeof(double));
	scratch4 = (double *)malloc((size_t)num_points*sizeof(double));
	scratch5 = (double *)malloc((size_t)num_points*sizeof(double));
	if ((ell==NULL) || (diag==NULL) || (offdiag==NULL) || (offsq==NULL) ||
		(evecs==NULL) || (ip==NULL) || (scratch1==NULL) || (scratch2 == NULL) ||
		(scratch3==NULL) || (scratch4==NULL) || (scratch5==NULL))
	{
		fprintf(stderr, "slepian_tapers() %s:%d ", __FILE__, __LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	/* make the diagonal elements of the tridiag matrix  */
	for (i=0; i < num_points; i++)
	{
		ai = (double) (i);
		diag[i] = -cs*(((an-1.)/2.-ai))*(((an-1.)/2.-ai));
		offdiag[i] = -ai*(an-ai)/2.;
		offsq[i] = offdiag[i]*offdiag[i];
	}
	eps = 1.0e-13;
	m11 = 1;

	/* call the eispac routines to invert the tridiagonal system */
	jtridib(&num_points, &eps, diag, offdiag, offsq, &rlb, &rlu, &m11, &nwin, lam, ip, &ierr, scratch1, scratch2);
#if DIAG1
	fprintf(stderr, "slepian_tapers() %s:%d ",__FILE__,__LINE__);
	fprintf(stderr,  "ierr=%d rlb=%.8f rlu=%.8f\n", ierr, rlb, rlu);
	fprintf(stderr,  "eigenvalues for the eigentapers\n");
	for (k=0; k < nwin; k++)
		fprintf(stderr,  "%.20f ", lam[k]);
	fprintf(stderr,  "\n");
#endif

	jtinvit(&num_points, &num_points, diag, offdiag, offsq, &nwin, lam, ip, evecs, &ierr, scratch1, scratch2, scratch3, scratch4, scratch5);

	free(ip);
	free(scratch1);
	free(scratch2);
	free(scratch3);
	free(scratch4);
	free(scratch5);

	/*
	** we calculate the eigenvalues of the dirichlet-kernel problem i.e.
	** the bandwidth retention factors from slepian 1978 asymptotic
	** formula, gotten from thomson 1982 eq 2.5 supplemented by the
	** asymptotic formula for k near 2n from slepian 1978 eq 61 more
	** precise values of these parameters, perhaps useful in adaptive
	** spectral estimation, can be calculated explicitly using the
	** rayleigh-quotient formulas in thomson (1982) and park et al (1987)
	*/
	dfac = (double) an * DPI * ww;
	drat = (double) 8.0 * dfac;
	dfac = (double) 4.0 * sqrt(DPI * dfac) * exp((double)(-2.0) * dfac);

	for (k=0; k < nwin; k++)
	{
		lam[k] = (double) 1.0 - (double) dfac;
		dfac = dfac * drat/(double) (k+1);
		/* fails as k -> 2n */
	}

	gamma = log((double)8.0 * an * sin((double) 2.0 * DPI * ww)) + (double)0.5772156649;

	for (k=0; k < nwin; k++)
	{
		bh = -2.0 * DPI * (an * ww - (double)(k) / (double)2.0 - (double)0.25)/gamma;
		ell[k] = (double) 1./((double) 1.+exp(DPI*(double) bh));
	}

	for (i=0; i < nwin; i++)
		lam[i] = MAX(ell[i], lam[i]);

	/************************************************************
	c   normalize the eigentapers to preserve power for a white process
	c   i.e. they have rms value unity
	c  tapsum is the average of the eigentaper, should be near zero for
	c  antisymmetric tapers
	************************************************************/

	for (k=0; k < nwin; k++)
	{
		kk = (k)*num_points;
		tapsum[k] = 0.;
		tapsq = 0.;
		for (i=0; i < num_points; i++)
		{
			aa = evecs[i+kk];
			tapers[i+kk] = aa;
			tapsum[k] = tapsum[k] + aa;
			tapsq = tapsq + aa * aa;
		}
		aa = sqrt(tapsq/(double) num_points);
		tapsum[k] = tapsum[k]/aa;

		for (i=0; i < num_points; i++)
		{
			tapers[i+kk] = tapers[i+kk]/aa;
		}
	}

	/* Free Memory */
	free(ell);
	free(diag);
	free(offdiag);
	free(offsq);
	free(evecs);
	return;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*----------------adwait--------------------------------*/
int adwait(double *sqr_spec, double *dcf, double *el, int nwin, int num_freq, double *ares, double *degf, double avar)
{
	/*
	** This version uses thomson's algorithm for calculating the
	** adaptive spectrum estimate
	*/
	double          as, das, tol, a1, scale, ax, fn, fx;
	double         *spw, *bias;
	double          test_tol, dif;
	int             jitter, i,  k, kpoint, jloop;
	float           df;
	/*
	** set tolerance for iterative scheme exit
	*/

#if 0
	fprintf(stderr, "adwait() %s:%d ", __FILE__, __LINE__);
	fprintf(stderr, "test input\n adwait: %d %d %f\n", nwin, num_freq, avar);
	fprintf(stderr, "\n Data=\n");
	for (i=0; i < num_freq; i++)
	{
		fprintf(stderr, "%d %f \n", i, sqr_spec[i]);
	}
#endif


	tol=3.0e-3;
	jitter=0;
	scale=avar;
	/***********************************
	c  we scale the bias by the total variance of the frequency transform
	c  from zero freq to the nyquist
	c  in this application we scale the eigenspectra by the bias in order to avoid
	c  possible floating point overflow
	************************************/
	spw=(double *) malloc((size_t) nwin*sizeof(double));
	bias=(double *) malloc((size_t) nwin*sizeof(double));
	if (spw==NULL || bias==NULL)
	{
		fprintf(stderr, "adwait() %s:%d ", __FILE__, __LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	for (i=0; i < nwin; i++)
	{
		bias[i]=(1.00-el[i]);
	}

	/* START do 100 */
	for (jloop=0; jloop < num_freq; jloop++)
	{
		for (i=0; i < nwin; i++)
		{
			kpoint=jloop+i*num_freq;
			spw[i]=(sqr_spec[kpoint])/scale;
		}
		/********************************************
		** first guess is the average of the two
		** lowest-order eigenspectral estimates
		********************************************/
		as=(spw[0]+spw[1])/2.00;

		/* START do 300 */
		/* c  find coefficients */

		for (k=0; k < 50; k++)
		{
			fn=0.00;
			fx=0.00;

			for (i=0; i < nwin; i++)
			{
				a1=sqrt(el[i])*as/(el[i]*as+bias[i]);
				a1=a1*a1;
				fn=fn+a1*spw[i];
				fx=fx+a1;
			}

			ax=fn/fx;
			dif=ax-as;
			das=fabs(dif);
			// fprintf(stderr, "adwait: jloop=%d k=%d %g %g %g %g\n",jloop,k, fn, fx, ax, das);
			test_tol=das/as;
			if (test_tol < tol)
				break;
			as=ax;
		}

		// fprintf(stderr, "adwait: k=%d test_tol=%f\n",k, test_tol);

		/* flag if iteration does not converge */
		if (k >= 50)
			jitter++;

		ares[jloop] = as*scale;

		/* output sqrt of power spectral density, printing error message if necessary */
		if (ares[jloop] >= 0.0)
		{
			ares[jloop] = sqrt(ares[jloop]);
		}
		else
		{
			fprintf(stderr, "adwait() %s:%d ", __FILE__, __LINE__);
			fprintf(stderr, "Error: trying to take square root of ares[jloop=%d] = %e\n", jloop, ares[jloop]);
			exit(-1);
		}
	
		/* calculate degrees of freedom */
		df=0.0;
		for (i=0; i < nwin; i++)
		{
			kpoint=jloop+i*num_freq;
			dcf[kpoint]=sqrt(el[i])*as/(el[i]*as+bias[i]);
			df=df+dcf[kpoint]*dcf[kpoint];
		}
		/*
		** we normalize degrees of freedom by the weight of the first
		** eigenspectrum this way we never have fewer than two
		** degrees of freedom
		*/
		degf[jloop]=df*2./(dcf[jloop]*dcf[jloop]);
	}

	//fprintf(stderr, "adwait() %s:%d ", __FILE__, __LINE__);
	//fprintf(stderr, "%d failed iterations\n", jitter);
	free(spw);
	free(bias);

	return jitter;
}

/*--------------------------------------------------------*/
/*------------get_F_values----------------------------*/
	/*
	sr contains the real part of the spectrum
	si contains the imaginary part of the spectrum
	nf is the number of frequencies (2^n+1), typically
	nwin is the number of tapers
	Fvalue is the array of f-statistic values
	b is fft of slepian eigentapers at zero freq
	cest is the estimated complex coefficient of a complex exponential (array of length 2 nf)
	*/
/*--------------------------------------------------------*/
void get_F_values2(double *sr, double *si, int nf, int nwin, float *Fvalue, double *b, float *cest)
{
	double sum, sumr, sumi, sum2;
	int i, j, k;
	sum = 0.0;

	for (i=0; i < nwin; i++)
	{
		sum += b[i] * b[i];
	}
	for (i=0; i < nf; i++)
	{
		cest[2*i] = cest[2*i+1] = 0.0;
		for (j=0; j < nwin; j++)
		{
			k = i+j*nf;
			cest[2*i] += sr[k]*b[j];
			cest[2*i+1] += si[k]*b[j];
		}
		cest[2*i] /= sum;
		cest[2*i+1] /= sum;
		sum2 = 0.0;
		for (j=0; j < nwin; j++)
		{
			k = i + j * nf;
			sumr = sr[k] - cest[2*i]*b[j];
			sumi = si[k] - cest[2*i+1]*b[j];
			sum2 = sum2 + sumr * sumr + sumi * sumi;
		}
		Fvalue[i] = (float)(nwin-1) * (SQR(cest[2*i]) + SQR(cest[2*i+1])) * sum / sum2;
	}
	return;
}


/*--------------------------------------------------------*/
/*------------get_F_values----------------------------*/
	/*
	sri contains the real/imag part of the spectrum
	nf is the number of frequencies (2^n+1), typically
	nwin is the number of tapers
	Fvalue is the array of f-statistic values
	b is fft of slepian eigentapers at zero freq
	cest is the estimated complex coefficient of a complex exponential (array of length 2 nf)
	*/
/*--------------------------------------------------------*/
void get_F_values(double *sri, int nf, int nwin, float *Fvalue, double *b, float *cest,int imin,int imax)
{
	double sum,sumr,alpha,cr,ci;
	double *ssri;
	int i, j, twonf;
	unsigned int i2;
	sum = 0.0;

	if (imin<0) imin = 0;
	if (imax>nf) imax = nf;

	/* twice the number of frequencies */
	twonf = 2 * imax;

	/* normalization of Slepian FFT's at DC */
	for (i=0; i<nwin; i++)
	{
		sum += b[i] * b[i];
	}

	/* estimate coefficients */
	alpha = b[0]/sum;
	for (i=2 * imin; i<twonf; i++)
	{
		cest[i] = alpha * sri[i];
	}

	for (j=1; j<nwin; j++)
	{
		alpha = b[j]/sum;
		ssri = sri + 2 * j * nf;
		for (i=2 * imin; i<twonf; i++)
		{
			cest[i] += alpha * ssri[i];
		}
	}

	/* clear Fvalues */
	for (i=0; i<nf; i++)
		Fvalue[i] = 0.0;

	/* store difference between estimated and actual coefficient squares */
	for (j=0; j<nwin; j++)
	{
		alpha = b[j];
		ssri = sri + 2 * j * nf;
		for (i2=2 * imin; i2<twonf; i2++)
		{
			sumr = ssri[i2]-alpha * cest[i2];
			Fvalue[i2/2] += sumr * sumr;
		}
	}

	/* now finish computing the Fvalues */
	sum *= (nwin-1);
	for (i=imin; i<imax; i++)
	{
		cr = cest[2 * i];
		ci = cest[2 * i + 1];
		Fvalue[i] = sum * (cr * cr + ci * ci)/Fvalue[i];
	}

	return;
}

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/*-------------  HIRES  ----------------------------------*/
/*-------------------------------------------------------*/
int hires(double *sqr_spec, double *el, int nwin, int num_freq, double *ares)
{
	int i, j, k, kpoint;
	float a;

	for (j=0; j < num_freq; j++)
		ares[j]=0.0;

	for (i=0; i < nwin; i++)
	{
		k = i * num_freq;
		a = 1.0 / (el[i] * nwin);
		for (j=0; j < num_freq; j++)
		{
			kpoint = j + k;
			ares[j] = ares[j] + a * (sqr_spec[kpoint]);
		}
	}

	for (j=0; j < num_freq; j++)
	{
		if (ares[j] > 0.0)
			ares[j]=sqrt(ares[j]);
		else
		{
			fprintf(stderr, "hires() %s:%d ", __FILE__, __LINE__);
			fprintf(stderr, "Error: trying to take square root of ares[j=%d] = %e\n", j,ares[j]);
			exit(-1);
		}	
	}

	return 1;
}
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*---------------jtinvit----------------------------------*/
/* #define dfabs(x) ffabs(x) */
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)
/* ./ add name=tinvit */
/* ------------------------------------------------------------------ */
int jtinvit(int *nm, int *n, double *d, double *e, double *e2, int *m, double *w, int *ind, double *z, int *ierr, double *rv1, double *rv2, double *rv3, double *rv4, double *rv6)
{
	/* Initialized data */
	static double   machep=1.25e-15;

	/* System generated locals */
	int             z_dim1, z_offset, i1, i2, i3;
	double          d1, d2;

	/* Local variables */
	static double   norm;
	static int      i, j, p, q, r, s;
	static double   u, v, order;
	static int      group;
	static double   x0, x1;
	static int      ii, jj, ip;
	static double   uk, xu;
	static int      tag, its;
	static double   eps2, eps3, eps4;
	static double   rtem;

	/* this subroutine is a translation of the inverse iteration tech- */
	/* nique in the algol procedure tristurm by peters and wilkinson. */
	/* handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). */

	/* this subroutine finds those eigenvectors of a tridiagonal */
	/* symmetric matrix corresponding to specified eigenvalues, */
	/* using inverse iteration. */

	/* on input: */

	/* nm must be set to the row dimension of two-dimensional */
	/* array parameters as declared in the calling program */
	/* dimension statement; */

	/* n is the order of the matrix; */

	/* d contains the diagonal elements of the input matrix; */

	/* e contains the subdiagonal elements of the input matrix */
	/* in its last n-1 positions.  e(1) is arbitrary; */

	/* e2 contains the squares of the corresponding elements of e, */
	/* with zeros corresponding to negligible elements of e. */
	/* e(i) is considered negligible if it is not larger than */
	/* the product of the relative machine precision and the sum */
	/* of the magnitudes of d(i) and d(i-1).  e2(1) must contain */
	/* 0.0d0 if the eigenvalues are in ascending order, or 2.0d0 */
	/* if the eigenvalues are in descending order.  if  bisect, */
	/* tridib, or  imtqlv  has been used to find the eigenvalues, */
	/* their output e2 array is exactly what is expected here; */

	/* m is the number of specified eigenvalues; */

	/*
	*w contains the m eigenvalues in ascending or descending order;
	 */

	/* ind contains in its first m positions the submatrix indices */
	/* associated with the corresponding eigenvalues in w -- */
	/* 1 for eigenvalues belonging to the first submatrix from */
	/*
	*the top, 2 for those belonging to the second submatrix, etc.
	 */

	/* on output: */

	/* all input arrays are unaltered; */

	/* z contains the associated set of orthonormal eigenvectors. */
	/* any vector which fails to converge is set to zero; */

	/* ierr is set to */
	/* zero       for normal return, */
	/* -r         if the eigenvector corresponding to the r-th */
	/* eigenvalue fails to converge in 5 iterations; */

	/* rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays. */

	/* questions and comments should be directed to b. s. garbow, */
	/* applied mathematics division, argonne national laboratory */

	/*
	*------------------------------------------------------------------
	 */

	/* :::::::::: machep is a machine dependent parameter specifying */
	/* the relative precision of floating point arithmetic. */
	/* machep=16.0d0**(-13) for long form arithmetic */
	/* on s360 :::::::::: */
	/* for f_floating dec fortran */
	/* data machep/1.1d-16/ */
	/* for g_floating dec fortran */
	/* Parameter adjustments */
	--rv6;
	--rv4;
	--rv3;
	--rv2;
	--rv1;
	--e2;
	--e;
	--d;
	z_dim1=*nm;
	z_offset=z_dim1+1;
	z -= z_offset;
	--ind;
	--w;

	/* Function Body */

	*ierr=0;
	if (*m == 0) {
		goto L1001;
	}
	tag=0;
	order=1.-e2[1];
	q=0;
	/* :::::::::: establish and process next submatrix :::::::::: */
L100:
	p=q+1;

	i1=*n;
	for (q=p; q <= i1; ++q) {
		if (q == *n) {
			goto L140;
		}
		if (e2[q+1] == 0.) {
			goto L140;
		}
		/* L120: */
	}
	/* :::::::::: find vectors by inverse iteration :::::::::: */
L140:
	++tag;
	s=0;

	i1=*m;
	for (r=1; r <= i1; ++r) {
		if (ind[r] != tag) {
			goto L920;
		}
		its=1;
		x1=w[r];
		if (s != 0) {
			goto L510;
		}
		/* :::::::::: check for isolated root :::::::::: */
		xu=1.;
		if (p != q) {
			goto L490;
		}
		rv6[p]=1.;
		goto L870;
L490:
		norm=(d1=d[p], fabs(d1));
		ip=p+1;

		i2=q;
		for (i=ip; i <= i2; ++i) {
			/* L500: */
			norm=norm+(d1=d[i], fabs(d1))+(d2=e[i], fabs(d2));
		}
		/* :::::::::: eps2 is the criterion for grouping, */
		/* eps3 replaces zero pivots and equal */
		/* roots are modified by eps3, */
		/*
		*eps4 is taken very small to avoid overflow ::::::::: :
		 */
		eps2=norm*.001;
		eps3=machep*norm;
		uk=(double) (q-p+1);
		eps4=uk*eps3;
		uk=eps4/sqrt(uk);
		s=p;
L505:
		group=0;
		goto L520;
		/* :::::::::: look for close or coincident roots :::::::::: */
L510:
		if ((d1=x1-x0, fabs(d1)) >= eps2) {
			goto L505;
		}
		++group;
		if (order*(x1-x0) <= 0.) {
			x1=x0+order*eps3;
		}
		/* :::::::::: elimination with interchanges and */
		/* initialization of vector :::::::::: */
L520:
		v=0.;

		i2=q;
		for (i=p; i <= i2; ++i) {
			rv6[i]=uk;
			if (i == p) {
				goto L560;
			}
			if ((d1=e[i], fabs(d1)) < fabs(u)) {
				goto L540;
			}
			/*
			*:::::::::: warning -- a divide check may occur
			*here if
			 */
			/*
			*e2 array has not been specified correctly ::::::
			*::::
			 */
			xu=u/e[i];
			rv4[i]=xu;
			rv1[i-1]=e[i];
			rv2[i-1]=d[i]-x1;
			rv3[i-1]=0.;
			if (i != q) {
				rv3[i-1]=e[i+1];
			}
			u=v-xu*rv2[i-1];
			v=-xu*rv3[i-1];
			goto L580;
	L540:
			xu=e[i]/u;
			rv4[i]=xu;
			rv1[i-1]=u;
			rv2[i-1]=v;
			rv3[i-1]=0.;
	L560:
			u=d[i]-x1-xu*v;
			if (i != q) {
				v=e[i+1];
			}
	L580:
			;
		}

		if (u == 0.) {
			u=eps3;
		}
		rv1[q]=u;
		rv2[q]=0.;
		rv3[q]=0.;
		/* :::::::::: back substitution */
		/* for i=q step -1 until p do -- :::::::::: */
L600:
		i2=q;
		for (ii=p; ii <= i2; ++ii) {
			i=p+q-ii;
			rtem=rv6[i]-u*rv2[i]-v*rv3[i];
			rv6[i]=(rtem)/rv1[i];
			v=u;
			u=rv6[i];
			/* L620: */
		}
		/* :::::::::: orthogonalize with respect to previous */
		/* members of group :::::::::: */
		if (group == 0) {
			goto L700;
		}
		j=r;

		i2=group;
		for (jj=1; jj <= i2; ++jj) {
	L630:
			--j;
			if (ind[j] != tag) {
				goto L630;
			}
			xu=0.;

			i3=q;
			for (i=p; i <= i3; ++i) {
				/* L640: */
				xu += rv6[i]*z[i+j*z_dim1];
			}

			i3=q;
			for (i=p; i <= i3; ++i) {
				/* L660: */
				rv6[i] -= xu*z[i+j*z_dim1];
			}

			/* L680: */
		}

L700:
		norm=0.;

		i2=q;
		for (i=p; i <= i2; ++i) {
			/* L720: */
			norm += (d1=rv6[i], fabs(d1));
		}

		if (norm >= 1.) {
			goto L840;
		}
		/* :::::::::: forward substitution :::::::::: */
		if (its == 5) {
			goto L830;
		}
		if (norm != 0.) {
			goto L740;
		}
		rv6[s]=eps4;
		++s;
		if (s > q) {
			s=p;
		}
		goto L780;
L740:
		xu=eps4/norm;

		i2=q;
		for (i=p; i <= i2; ++i) {
			/* L760: */
			rv6[i] *= xu;
		}
		/* :::::::::: elimination operations on next vector */
		/* iterate :::::::::: */
L780:
		i2=q;
		for (i=ip; i <= i2; ++i) {
			u=rv6[i];
			/*
			*:::::::::: if rv1(i-1) .eq. e(i), a row
			*interchange
			 */
			/* was performed earlier in the */
			/* triangularization process :::::::::: */
			if (rv1[i-1] != e[i]) {
				goto L800;
			}
			u=rv6[i-1];
			rv6[i-1]=rv6[i];
	L800:
			rv6[i]=u-rv4[i]*rv6[i-1];
			/* L820: */
		}

		++its;
		goto L600;
		/*
		*:::::::::: set error -- non-converged eigenvector
		*::::::::::
		 */
L830:
		*ierr=-r;
		xu=0.;
		goto L870;
		/* :::::::::: normalize so that sum of squares is */
		/* 1 and expand to full order :::::::::: */
L840:
		u=0.;

		i2=q;
		for (i=p; i <= i2; ++i) {
			/* L860: */
			/* Computing 2nd power */
			d1=rv6[i];
			u += d1*d1;
		}

		xu=1./sqrt(u);

L870:
		i2=*n;
		for (i=1; i <= i2; ++i) {
			/* L880: */
			z[i+r*z_dim1]=0.;
		}

		i2=q;
		for (i=p; i <= i2; ++i) {
			/* L900: */
			z[i+r*z_dim1]=rv6[i]*xu;
		}

		x0=x1;
L920:
		;
	}

	if (q < *n) {
		goto L100;
	}
L1001:
	return 0;
	/* :::::::::: last card of tinvit :::::::::: */
}				/* tinvit_ */


/*--------------------------------------------------------*/
/*----------------jtridib---------------------------------*/
/*--------------------------------------------------------*/
int jtridib(int *n, double *eps1, double *d, double *e, double *e2, double *lb, double *ub, int *m11, int *m, double *w, int *ind, int *ierr, double *rv4, double *rv5)
{
	/* Initialized data */

	static double   machep=1.25e-15;

	/* System generated locals */
	int             i1, i2;
	double          d1, d2, d3;

	/* Local variables */
	static int      i, j, k, l, p, q, r, s;
	static double   u, v;
	static int      m1, m2;
	static double   t1, t2, x0, x1;
	static int      m22, ii;
	static double   xu;
	static int      isturm, tag;



	/* this subroutine is a translation of the algol procedure bisect, */
	/* num. math. 9, 386-393(1967) by barth, martin, and wilkinson. */
	/* handbook for auto. comp., vol.ii-linear algebra, 249-256(1971). */

	/* this subroutine finds those eigenvalues of a tridiagonal */
	/* symmetric matrix between specified boundary indices, */
	/* using bisection. */

	/* on input: */

	/* n is the order of the matrix; */

	/* eps1 is an absolute error tolerance for the computed */
	/* eigenvalues.  if the input eps1 is non-positive, */
	/* it is reset for each submatrix to a default value, */
	/* namely, minus the product of the relative machine */
	/* precision and the 1-norm of the submatrix; */

	/* d contains the diagonal elements of the input matrix; */

	/* e contains the subdiagonal elements of the input matrix */
	/* in its last n-1 positions.  e(1) is arbitrary; */

	/* e2 contains the squares of the corresponding elements of e. */
	/* e2(1) is arbitrary; */

	/* m11 specifies the lower boundary index for the desired */
	/* eigenvalues; */

	/* m specifies the number of eigenvalues desired.  the upper */
	/* boundary index m22 is then obtained as m22=m11+m-1. */

	/* on output: */

	/* eps1 is unaltered unless it has been reset to its */
	/* (last) default value; */

	/* d and e are unaltered; */

	/* elements of e2, corresponding to elements of e regarded */
	/* as negligible, have been replaced by zero causing the */
	/* matrix to split into a direct sum of submatrices. */
	/* e2(1) is also set to zero; */

	/* lb and ub define an interval containing exactly the desired */
	/* eigenvalues; */

	/* w contains, in its first m positions, the eigenvalues */
	/* between indices m11 and m22 in ascending order; */

	/* ind contains in its first m positions the submatrix indices */
	/* associated with the corresponding eigenvalues in w -- */
	/* 1 for eigenvalues belonging to the first submatrix from */
	/*
	*the top, 2 for those belonging to the second submatrix, etc.;
	 */

	/* ierr is set to */
	/* zero       for normal return, */
	/* 3*n+1      if multiple eigenvalues at index m11 make */
	/* unique selection impossible, */
	/* 3*n+2      if multiple eigenvalues at index m22 make */
	/* unique selection impossible; */

	/* rv4 and rv5 are temporary storage arrays. */

	/* note that subroutine tql1, imtql1, or tqlrat is generally faster */
	/* than tridib, if more than n/4 eigenvalues are to be found. */

	/* questions and comments should be directed to b. s. garbow, */
	/* applied mathematics division, argonne national laboratory */

	/*
	*------------------------------------------------------------------
	 */

	/* :::::::::: machep is a machine dependent parameter specifying */
	/* the relative precision of floating point arithmetic. */
	/* machep=16.0d0**(-13) for long form arithmetic */
	/* on s360 :::::::::: */
	/* for f_floating dec fortran */
	/* data machep/1.1d-16/ */
	/* for g_floating dec fortran */
	/* Parameter adjustments */
	--rv5;
	--rv4;
	--e2;
	--e;
	--d;
	--ind;
	--w;

	/* Function Body */

	*ierr=0;
	tag=0;
	xu=d[1];
	x0=d[1];
	u=0.;
	/* :::::::::: look for small sub-diagonal entries and determine an */
	/* interval containing all the eigenvalues :::::::::: */
	i1=*n;
	for (i=1; i <= i1; ++i) {
		x1=u;
		u=0.;
		if (i != *n) {
			u=(d1=e[i+1], fabs(d1));
		}
		/* Computing MIN */
		d1=d[i]-(x1+u);
		xu=min(d1, xu);
		/* Computing MAX */
		d1=d[i]+(x1+u);
		x0=max(d1, x0);
		if (i == 1) {
			goto L20;
		}
		if ((d1=e[i], fabs(d1)) > machep*((d2=d[i], fabs(d2))+(
						 d3=d[i-1], fabs(d3)))) {
			goto L40;
		}
L20:
		e2[i]=0.;
L40:
		;
	}

	/* Computing MAX */
	d1=fabs(xu), d2=fabs(x0);
	x1=max(d1, d2)*machep*(double) (*n);
	xu -= x1;
	t1=xu;
	x0 += x1;
	t2=x0;
	/* :::::::::: determine an interval containing exactly */
	/* the desired eigenvalues :::::::::: */
	p=1;
	q=*n;
	m1=*m11-1;
	if (m1 == 0) {
		goto L75;
	}
	isturm=1;
L50:
	v=x1;
	x1=xu+(x0-xu)*.5;
	if (x1 == v) {
		goto L980;
	}
	goto L320;
L60:
	if ((i1=s-m1) < 0) {
		goto L65;
	} else if (i1 == 0) {
		goto L73;
	} else {
		goto L70;
	}
L65:
	xu=x1;
	goto L50;
L70:
	x0=x1;
	goto L50;
L73:
	xu=x1;
	t1=x1;
L75:
	m22=m1+*m;
	if (m22 == *n) {
		goto L90;
	}
	x0=t2;
	isturm=2;
	goto L50;
L80:
	if ((i1=s-m22) < 0) {
		goto L65;
	} else if (i1 == 0) {
		goto L85;
	} else {
		goto L70;
	}
L85:
	t2=x1;
L90:
	q=0;
	r=0;
	/* :::::::::: establish and process next submatrix, refining */
	/* interval by the gerschgorin bounds :::::::::: */
L100:
	if (r == *m) {
		goto L1001;
	}
	++tag;
	p=q+1;
	xu=d[p];
	x0=d[p];
	u=0.;

	i1=*n;
	for (q=p; q <= i1; ++q) {
		x1=u;
		u=0.;
		v=0.;
		if (q == *n) {
			goto L110;
		}
		u=(d1=e[q+1], fabs(d1));
		v=e2[q+1];
L110:
		/* Computing MIN */
		d1=d[q]-(x1+u);
		xu=min(d1, xu);
		/* Computing MAX */
		d1=d[q]+(x1+u);
		x0=max(d1, x0);
		if (v == 0.) {
			goto L140;
		}
		/* L120: */
	}

L140:
	/* Computing MAX */
	d1=fabs(xu), d2=fabs(x0);
	x1=max(d1, d2)*machep;
	if (*eps1 <= 0.) {
		*eps1=-x1;
	}
	if (p != q) {
		goto L180;
	}
	/* :::::::::: check for isolated root within interval :::::::::: */
	if (t1 > d[p] || d[p] >= t2) {
		goto L940;
	}
	m1=p;
	m2=p;
	rv5[p]=d[p];
	goto L900;
L180:
	x1 *= (double) (q-p+1);
	/* Computing MAX */
	d1=t1, d2=xu-x1;
	*lb=max(d1, d2);
	/* Computing MIN */
	d1=t2, d2=x0+x1;
	*ub=min(d1, d2);
	x1=*lb;
	isturm=3;
	goto L320;
L200:
	m1=s+1;
	x1=*ub;
	isturm=4;
	goto L320;
L220:
	m2=s;
	if (m1 > m2) {
		goto L940;
	}
	/* :::::::::: find roots by bisection :::::::::: */
	x0=*ub;
	isturm=5;

	i1=m2;
	for (i=m1; i <= i1; ++i) {
		rv5[i]=*ub;
		rv4[i]=*lb;
		/* L240: */
	}
	/* :::::::::: loop for k-th eigenvalue */
	/* for k=m2 step -1 until m1 do -- */
	/*
	*(-do- not used to legalize -computed go to-) ::::::::::
	 */
	k=m2;
L250:
	xu=*lb;
	/* :::::::::: for i=k step -1 until m1 do -- :::::::::: */
	i1=k;
	for (ii=m1; ii <= i1; ++ii) {
		i=m1+k-ii;
		if (xu >= rv4[i]) {
			goto L260;
		}
		xu=rv4[i];
		goto L280;
L260:
		;
	}

L280:
	if (x0 > rv5[k]) {
		x0=rv5[k];
	}
	/* :::::::::: next bisection step :::::::::: */
L300:
	x1=(xu+x0)*.5;
	if (x0-xu <= machep*2.*(fabs(xu)+fabs(x0))+fabs(*eps1)) {
		goto L420;
	}
	/* :::::::::: in-line procedure for sturm sequence :::::::::: */
L320:
	s=p-1;
	u=1.;

	i1=q;
	for (i=p; i <= i1; ++i) {
		if (u != 0.) {
			goto L325;
		}
		v=(d1=e[i], fabs(d1))/machep;
		if (e2[i] == 0.) {
			v=0.;
		}
		goto L330;
L325:
		v=e2[i]/u;
L330:
		u=d[i]-x1-v;
		if (u < 0.) {
			++s;
		}
		/* L340: */
	}

	switch ((int) isturm) {
	case 1:
		goto L60;
	case 2:
		goto L80;
	case 3:
		goto L200;
	case 4:
		goto L220;
	case 5:
		goto L360;
	}
	/* :::::::::: refine intervals :::::::::: */
L360:
	if (s >= k) {
		goto L400;
	}
	xu=x1;
	if (s >= m1) {
		goto L380;
	}
	rv4[m1]=x1;
	goto L300;
L380:
	rv4[s+1]=x1;
	if (rv5[s] > x1) {
		rv5[s]=x1;
	}
	goto L300;
L400:
	x0=x1;
	goto L300;
	/* :::::::::: k-th eigenvalue found :::::::::: */
L420:
	rv5[k]=x1;
	--k;
	if (k >= m1) {
		goto L250;
	}
	/* :::::::::: order eigenvalues tagged with their */
	/* submatrix associations :::::::::: */
L900:
	s=r;
	r=r+m2-m1+1;
	j=1;
	k=m1;

	i1=r;
	for (l=1; l <= i1; ++l) {
		if (j > s) {
			goto L910;
		}
		if (k > m2) {
			goto L940;
		}
		if (rv5[k] >= w[l]) {
			goto L915;
		}
		i2=s;
		for (ii=j; ii <= i2; ++ii) {
			i=l+s-ii;
			w[i+1]=w[i];
			ind[i+1]=ind[i];
			/* L905: */
		}

L910:
		w[l]=rv5[k];
		ind[l]=tag;
		++k;
		goto L920;
L915:
		++j;
L920:
		;
	}

L940:
	if (q < *n) {
		goto L100;
	}
	goto L1001;
	/* :::::::::: set error -- interval cannot be found containing */
	/* exactly the desired eigenvalues :::::::::: */
L980:
	*ierr=*n*3+isturm;
L1001:
	*lb=t1;
	*ub=t2;
	return 0;
	/* :::::::::: last card of tridib :::::::::: */
}				/* tridib_ */


/*--------------------------------------------------------*/
/*---------------get_pow_2----------------------------*/
/*--------------------------------------------------------*/
int get_pow_2(int inum)
{
	int j, klength;
	/* find smallest power of 2 that encompasses the data */

	for (j=1; pow((double) 2, (double) j) < inum; j++);
	klength=pow((double) 2, (double) j);
	return klength;
}


/*--------------------------------------------------------*/
/*-----------remove_mean----------------------------*/
/*--------------------------------------------------------*/
double remove_mean(float x[], int lx)
{
	int k;
	double mean;

	mean = 0.0;
	if (lx < 2)
	{
		mean=x[0];
		x[0]=0.;
		return mean;
	}

	for (k=0; k < lx; k++)
	{
		mean=x[k]+mean;
	}
	mean = mean / (float)lx;

	for (k=0; k < lx; k++)
	{
		x[k] = x[k] - mean;
	}
	return mean;
}


/*********************************************************/
/*--------------------------------------------------------*/
/*-----------zero_pad----------------------------------*/
/*--------------------------------------------------------*/
void zero_pad(float output[], int start, int olength)
{
	int i;

	for (i=start; i < olength; i++)
	{
		output[i] = 0.0;
	}
}



/*********************************************************
** Functions from numerical recipes
** 
**********************************************************/
void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i, i1, i2, i3, i4, np3;
	float c1=0.5, c2, h1r, h1i, h2r, h2i;
	double wr, wi, wpr, wpi, wtemp, theta;

	theta = 3.141592653589793/(double) (n>>1);
	if (isign == 1)
	{
		c2 = -0.5;
		four1(data,n>>1,1);
	}
	else
	{
		c2 = 0.5;
		theta = -theta;
	}
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0+wpr;
	wi = wpi;
	np3 = n+3;
	for (i=2;i<=(n>>2);i++)
	{
		i4 = 1+(i3 = np3-(i2 = 1+(i1 = i+i-1)));
		h1r = c1*(data[i1]+data[i3]);
		h1i = c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i = c2*(data[i1]-data[i3]);
		data[i1] = h1r+wr*h2r-wi*h2i;
		data[i2] = h1i+wr*h2i+wi*h2r;
		data[i3] = h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr = (wtemp = wr)*wpr-wi*wpi+wr;
		wi = wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1)
	{
		data[1] = (h1r = data[1])+data[2];
		data[2] = h1r-data[2];
	}
	else
	{
		data[1] = c1*((h1r = data[1])+data[2]);
		data[2] = c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i=1;i<n;i+=2)
	{
		if (j > i)
		{
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax)
	{
		istep = mmax << 1;
		theta = isign*(6.28318530717959/mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1;m<mmax;m+=2)
		{
			for (i=m;i<=n;i+=istep)
			{
				j = i+mmax;
				tempr = wr*data[j]-wi*data[j+1];
				tempi = wr*data[j+1]+wi*data[j];
				data[j] = data[i]-tempr;
				data[j+1] = data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
	}
}
#undef SWAP


/*********************************************************************************
** my_multitaper_spectrum
**
** data=floating point input time series
** npoints=number of points in data
** kind=flag for choosing hires (1) or adaptive (2) weighting coefficients 
** nwin=number of taper windows to calculate
** nwdt=order of the slepian functions
** inorm=flag for choice of normalization
** dt=sampling interval (time)
** psd=output spectrum
** dof=degrees of freedom at each frequency
** Fvalues=Ftest value at each frequency estimate
** klen=number of frequecies calculated (power of 2)
** cest is the estimated complex coefficient of a complex exponential (array of length klen+2)
**
** lambda=vector of eigenvalues
** tapsum=sum of each taper, saved for use in adaptive weighting
** tapers= matrix of slepian tapers packed in a 1D double array
**
*********************************************************************************/
void my_multitaper_spectrum(float *data, int npoints, int nwin, int p_nwdt, float *psd, int kind)
{
	int             i, j;
	float          *b;
	int             iwin, kk;
	double          anrm, norm, sqramp, sum;
	double         *ReImSpec,*sqr_spec, *amu;
	float          *amp;
	int             num_freqs, len_taps, num_freq_tap,n1, n2, kf;
	double         *dcf, *degf, avar;
	/*************/
	static int save_npoints = 0,save_nwin = 0;
	static float save_nwdt = 0.0;
	static double   *lambda=NULL,*tapers=NULL,*tapsum=NULL;
	float nwdt = p_nwdt;
	int inorm = 0;
	float dt = 1.0;
	int klen = npoints;

	if ((npoints != save_npoints) || (nwin != save_nwin) || (nwdt != save_nwdt))
	{
		/* This means that we must compute tapers, etc */
		len_taps = npoints*nwin;
		lambda = (double *)realloc((void *)lambda,(size_t)nwin*sizeof(double));
		tapsum = (double *)realloc((void *)tapsum,(size_t)nwin*sizeof(double));
		tapers = (double *)realloc((void *)tapers,(size_t)len_taps*sizeof(double));
		if ((lambda == NULL) || (tapsum == NULL) || (tapers == NULL))
		{
			fprintf(stderr, "multitaper_spectrum() %s:%d",__FILE__,__LINE__);
			fprintf(stderr, "Insufficient memory - realloc returned NULL.\n");
			exit(-1);
		}
		save_npoints = npoints;
		save_nwin = nwin;
		save_nwdt = nwdt;

		/* get the set of Slepian tapers  */
		slepian_tapers(npoints, nwin, lambda, nwdt, tapers, tapsum);
	}

	/* choose normalization based on inorm flag  */
	anrm = 1.0;
	switch (inorm)
	{
		case 1:
			anrm = npoints;
			break;
		case 2:
			anrm = 1.0/dt;
			break;
		case 3:
			anrm = sqrt((double) npoints);
			break;
		default:
			anrm = 1.0;
			break;
	}

	num_freqs = 1 + klen/2;
	num_freq_tap = num_freqs * nwin;

	b = (float *)malloc((size_t)npoints * sizeof(float));
	amu = (double *)malloc((size_t)num_freqs * sizeof(double));
	sqr_spec = (double *)malloc((size_t)num_freq_tap * sizeof(double));
	ReImSpec = (double *)malloc(2 * (size_t)num_freq_tap * sizeof(double));
	amp = (float *)malloc((size_t)klen * sizeof(float));
	if ((b==NULL) || (amu==NULL) || (sqr_spec==NULL) || (ReImSpec==NULL) || (amp==NULL))
	{
		fprintf(stderr, "multitaper_spectrum() %s:%d",__FILE__,__LINE__);
		fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
		exit(-1);
	}

	/* apply the taper in the loop.  do this nwin times  */
	for (iwin=0; iwin < nwin; iwin++)
	{
		kk = iwin*npoints;
		kf = iwin*num_freqs;

		/* application of iwin-th taper   */
		for (j=0; j < npoints; j++)
			b[j] = data[j] * tapers[kk+j];

		/* calculate the eigenspectrum */
		mt_get_spec(b, npoints, klen, amp);					

		/* check that indices are within bounds */
		if ((2*num_freqs-3 > klen) || (num_freqs-1+kf > num_freq_tap))
		{
			fprintf(stderr, "multitaper_spectrum() %s:%d ",__FILE__,__LINE__);
			fprintf(stderr, "error in index\n");
		}

		/* get spectrum from real fourier transform    */
		norm = 1.0/(anrm*anrm);
		sum = 0.0;
		for (i=1; i<num_freqs-1; i++)
		{
			sqramp = SQR(amp[2*i+1])+SQR(amp[2*i]);
			ReImSpec[2*(i+kf)] = amp[2*i];
			ReImSpec[2*(i+kf)+1] = amp[2*i+1];
			sqr_spec[i+kf] = norm*sqramp;
			sum += sqramp;
		}

		/* do DC and Nyquist */
		sqr_spec[kf] = norm*SQR(amp[0]);
		sqr_spec[num_freqs-1+kf] = norm*SQR(amp[1]);
		ReImSpec[2*kf] = amp[0];
		ReImSpec[2*kf+1] = 0.0;
		ReImSpec[2*(num_freqs-1+kf)] = amp[1];
		ReImSpec[2*(num_freqs-1+kf)+1] = 0.0;
		sum += sqr_spec[kf]+sqr_spec[num_freqs-1+kf];
	}
	free(amp);
	free(b);

	switch (kind)
	{
		case 1:
			hires(sqr_spec, lambda, nwin, num_freqs, amu);
			for (i=0; i < num_freqs; i++)
				psd[i] = amu[i];
		case 2:
			/* get avar=variance */
			n1=0;
			n2=npoints;
			avar=0.0;

			for (i=n1; i < n2; i++)
				avar += (data[i])*(data[i]);
			switch (inorm)
			{
				case 1:
					avar = (avar/npoints)/npoints;
					break;
				case 2:
					avar = avar*dt*dt;
					break;
				case 3:
					avar = avar/npoints;
					break;
				default:
					break;
			}

			dcf = (double *) malloc((size_t) num_freq_tap*sizeof(double));
			degf = (double *) malloc((size_t) num_freqs*sizeof(double));
			if ((dcf==NULL) || (degf==NULL))
			{
				fprintf(stderr, "multitaper_spectrum() %s:%d ",__FILE__,__LINE__);
				fprintf(stderr, "Insufficient memory - malloc returned NULL.\n");
				exit(-1);
			}

			adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);

			for (i=0; i < num_freqs; i++)
				psd[i] = amu[i];

			free(dcf);
			free(degf);
			break;
	}

	/* free up memory and return  */
	free(amu);
	free(sqr_spec);
	free(ReImSpec);
}

