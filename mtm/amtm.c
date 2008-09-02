#include <stdio.h>
#include <math.h>
#include <g2c.h>

#define maxlen 50000
#define zero 0.0
#define maxsignal 40
#define nlim 32768

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )


#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifdef notdef
int main(int argc, char *argv[])
{
	char name[256];
	int i;

	real fsignal[maxsignal + 1];
	real confsignal[maxsignal + 1];
	integer irec[maxsignal + 1];
	real a[maxlen + 1];
	real dummy[maxlen + 1];
	real dt,anpi,f1,f2,f0;
	integer nscan,nsignals;
	integer nf0,nbnd,iwhich;
	real fconf[7];
	real conf[5];
	real specraw0[nlim + 1];
	real specmed0[nlim + 1];
	real specresh0[nlim + 1];
	real harmonic0[nlim + 1];
	real specback0[nlim + 1];
	real ftest[nlim + 1];
	real whiteraw0;
	real whiterob0;
	real rhoraw0;
	real rhorob0;
	real tauraw0;
	real taurob0;
	real fsmooth;
	real df;
	real demn;
	real fray, fny, bndwdth, frange;
	real flow, fhigh;

	/* Set defaults */
	integer npi = 2;
	integer nwin = 3;
	integer inorm = 0;
	integer ispec = 1;
	integer iresh = 1;
	integer ithresh = 3;
	real fthresh = 95.0;
	integer inoise = 0;
	integer ismooth = 0;
	integer ilog = 0;
	integer isignal = 0;
	logical fset = FALSE;
	logical smooset = FALSE;

	for (i=1; i < argc; i++)
	{
		if (strcmp(argv[i], "-w") == 0)
		{
			nwin = strtol(argv[++i], NULL, 10);
		}
		else if (strcmp(argv[i], "-b") == 0)
		{
			npi = strtol(argv[++i], NULL, 10);
		}
		else if (strcmp(argv[i], "-f") == 0)
		{
			fset = TRUE;
			f1 = strtod(argv[++i], NULL);
			f2 = strtod(argv[++i], NULL);
		}
		else if (strcmp(argv[i], "-null") == 0)
		{
			if (strstr(argv[i+1], "red")) inoise = 0;
			if (strstr(argv[i+1], "white")) inoise = 1;
			if (strstr(argv[i+1], "lw")) inoise = 2;
			if (strstr(argv[i+1], ":m")) ilog = 0;
			if (strstr(argv[i+1], ":l")) ilog = 1;
			i++;
		}
		else if (strcmp(argv[i], "-m") == 0)
		{
			smooset = TRUE;
			ismooth = 1;
			fsmooth = strtod(argv[++i], NULL);
		}
		else if (strcmp(argv[i], "-n") == 0)
		{
			inorm = strtol(argv[++i], NULL, 10);
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			char *ch = argv[++i];
			if (*ch == 'h') ispec = 0;
			if (*ch == 'a') ispec = 1;
			if (strlen(ch) > 1)
			{
				ithresh = strtol(ch + 1, NULL, 10);
				iresh = 1;
			}
		}
		else if (strcmp(argv[i], "-s") == 0)
		{
			isignal = strtol(argv[++i], NULL, 10);
		}
		else if (strcmp(argv[i], "-") != 0)
		{
			strcpy(name, argv[i]);
		}
		else
		{
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
		}
	}

	/*
	** Read in data
	*****  Call Get1dField (name, ncols, nscan, a, dummy)
	*/
	{
		FILE *fp;
		char buf[256];
		int nlines = 0;

		fp = fopen(name, "r");
		while (fgets(buf, 255, fp) != NULL)
		{
			sscanf(buf, "%f", &(a[nlines]));
			dummy[nlines] = nlines;
			nlines++;
		}
		nscan = nlines;
	}

	/*
	** Calculate Rayleigh and Nyquist frequencies
	*/
	dt = (dummy[nscan - 1] - dummy[0]) / (nscan - 1);
	for (i=0; i < nscan; i++)
		demn += a[i];
	demn = demn / nscan;
	fray = 1.0 / (nscan * dt);
	fny = 0.5 / dt;
	bndwdth = 2.0 * npi * fray;

	if (!fset)
	{
		f1 = 0.0;
		f2 = fny;
	}
	frange = f2 - f1;

	if (!smooset) fsmooth = frange / 6.0;
	if ((fsmooth > (frange / 2.0)) || (fsmooth < bndwdth))
	{
		fsmooth = MIN(MAX(fsmooth, bndwdth), frange / 2.0);
		fprintf(stderr, "Your frequency-smoothing was out of bounds!\n");
		fprintf(stderr, "It has been reset to %f\n", fsmooth);
	}

	for (i=0; i < maxsignal; i++)
	{
		irec[i] = 0;
		fsignal[i] = 0.0;
		confsignal[i] = 0.0;
	}
	nsignals = 0;

	if (isignal == 2)
	{
		inoise = 2;
	}

	if (iresh == 1)
	{
		switch (++ithresh)
		{
			case 1: fthresh = 50.0; break;
			case 2: fthresh = 90.0; break;
			case 3: fthresh = 95.0; break;
			case 4: fthresh = 99.0; break;
			case 5: fthresh = 99.5; break;
			case 6: fthresh = 99.9; break;
		}
	}

	fprintf(stderr, "1) variance/resolution tradeoff:\n");
	fprintf(stderr, "   * resolution = %d f_R\n", npi);
	fprintf(stderr, "   * number of tapers = %d\n", nwin);
	fprintf(stderr, "2) null hypothesis: ");
	switch(inoise)
	{
		case 0: fprintf(stderr, "red noise\n"); break;
		case 1: fprintf(stderr, "white noise\n"); break;
		default: fprintf(stderr, "locally-white noise\n"); break;
	}
	if (isignal != 2)
	{
		if (ismooth == 0)
			fprintf(stderr, "   * raw noise background estimation\n");
		else
		{
			fprintf(stderr, "   * robust noise background estimation\n");
			fprintf(stderr, "   * & fsmooth = %f\n", fsmooth);
			if (inoise == 0)
			{
				if (ilog == 0)
					fprintf(stderr, "   (min misfit noise background)\n");
				else
					fprintf(stderr, "   (min log-misfit noise background)\n");
			}
		}
	}
	if (isignal == 0)
		fprintf(stderr, "3) signal assumption: harmonic or quasiperiodic\n");
	else if (isignal == 1)
		fprintf(stderr, "3) signal assumption: quasiperiodic\n");
	else
		fprintf(stderr, "3) signal assumption: harmonic\n");

	fprintf(stderr, "4) spectrum:\n");
	if (ispec == 1)
		fprintf(stderr, "    * adaptive estimate\n");
	else
		fprintf(stderr, "    * high-resolution estimate\n");

	if (inorm == 0)
		fprintf(stderr, "    * no normalization\n");
	else if (inorm == 1)
		fprintf(stderr, "    * normalize by N\n");
	else
		fprintf(stderr, "    * normalize by 1/dt\n");

	if (iresh == 0)
		fprintf(stderr, "    * unreshaped\n");
	else
		fprintf(stderr, "    * reshaped at threshold = %f\n", fthresh);

	fprintf(stderr, "    * frequency range: %f:%f\n", f1, f2);

	fprintf(stderr, "5) signals to reconstruct: %d\n", nsignals);
	fprintf(stderr, "6) display:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");

	spec_(a, &nwin, &npi, &dt, &nscan, &f1, &f2, &inoise, &ismooth, &ilog,
		&fsmooth, &isignal, &ispec, &inorm, &iresh, &ithresh, &nf0, &df,
		specraw0, specmed0, specresh0, harmonic0, specback0,
		ftest, conf, fconf, &nbnd, &whiteraw0, &whiterob0, &rhoraw0,
		&rhorob0, &tauraw0, &taurob0);

	findpeaks_(&flow, &fhigh, &isignal,
		&nf0, &df, &nbnd, specraw0, specmed0, specresh0, harmonic0,
		specback0, ftest,
		conf, fconf,
		fsignal, confsignal, &nsignals);

	if (nsignals > 0)
	{
		fprintf(stderr, "Significant Peaks\n");
		fprintf(stderr, "signal frequency signif%\n");
		for (i=0; i < nsignals; i++)
		{
			if (fsignal[i] < bndwdth)
				fprintf(stderr, "(%2d)   TREND       %2d%\n", i, fsignal[i], (int)((confsignal[i] + 0.5)));
			else
				fprintf(stderr, "(%2d)   %7.4f       %2d%\n", i, fsignal[i], (int)((confsignal[i] + 0.5)));
		}
	}

	for (i=0; i < nf0 - 1; i++)
	{
		printf("%9f\n", specraw0[i]);
	}
}
#endif


int mtmpsd(float *buf, int nsamples,
	int p_nwin,
	int p_npi,
	float p_f1,
	float p_f2,
	int p_inoise,
	int p_ilog,
	float p_fsmooth,
	int p_inorm,
	int p_ispec,
	int p_ithresh,
	int p_isignal)
{
	int i;
	real fsignal[maxsignal + 1], confsignal[maxsignal + 1];
	real a[maxlen + 1], dummy[maxlen + 1];
	real dt,anpi,f1,f2,f0;
	real fconf[6 + 1], conf[4 + 1];
	real specraw0[nlim + 1], specmed0[nlim + 1], specresh0[nlim + 1];
	real harmonic0[nlim + 1], specback0[nlim + 1], ftest[nlim + 1];
	real whiteraw0, whiterob0, rhoraw0, rhorob0, tauraw0, taurob0, fsmooth;
	real df, demn, fray, fny, bndwdth, frange, flow, fhigh;
	integer irec[maxsignal + 1];
	integer nscan,nsignals, nf0,nbnd,iwhich;

	/* Set defaults */
	integer npi = 2;
	integer nwin = 3;
	integer inorm = 0;
	integer ispec = 1;
	integer iresh = 1;
	integer ithresh = 3;
	real fthresh = 95.0;
	integer inoise = 0;
	integer ismooth = 0;
	integer ilog = 0;
	integer isignal = 0;
	logical fset = FALSE;
	logical smooset = FALSE;

	if (p_nwin > 0) nwin = p_nwin;
	if (p_npi > 0) npi = p_npi;
	if (p_f1 < p_f2)
	{
		fset = TRUE;
		f1 = p_f1;
		f2 = p_f2;
	}
	if (p_inoise >= 0) inoise = p_inoise;
	if (p_ilog >= 0) ilog = p_ilog;
	if (p_fsmooth >= 0) { fsmooth = p_fsmooth; ismooth = 1; }
	if (p_inorm >= 0) inorm = p_inorm;
	if (p_ispec >= 0) ispec = p_ispec;
	if (p_ithresh >= 0) { ithresh = p_ithresh; iresh = 1; }
	if (p_isignal >= 0) isignal = p_isignal;


	/*
	** Read in data
	*/
	for (i=0; i < nsamples; i++)
	{
		dummy[i] = i;
		a[i] = buf[i];
	}
	nscan = nsamples;

	/*
	** Calculate Rayleigh and Nyquist frequencies
	*/
	dt = (dummy[nscan - 1] - dummy[0]) / (nscan - 1);
	for (i=0; i < nscan; i++)
		demn += a[i];
	demn = demn / nscan;
	fray = 1.0 / (nscan * dt);
	fny = 0.5 / dt;
	bndwdth = 2.0 * npi * fray;

	if (!fset)
	{
		f1 = 0.0;
		f2 = fny;
	}
	frange = f2 - f1;

	if (!smooset) fsmooth = frange / 6.0;
	if ((fsmooth > (frange / 2.0)) || (fsmooth < bndwdth))
	{
		fsmooth = MIN(MAX(fsmooth, bndwdth), frange / 2.0);
		fprintf(stderr, "Your frequency-smoothing was out of bounds!\n");
		fprintf(stderr, "It has been reset to %f\n", fsmooth);
	}

	for (i=0; i < maxsignal; i++)
	{
		irec[i] = 0;
		fsignal[i] = 0.0;
		confsignal[i] = 0.0;
	}
	nsignals = 0;

	if (isignal == 2)
	{
		inoise = 2;
	}

	if (iresh == 1)
	{
		switch (++ithresh)
		{
			case 1: fthresh = 50.0; break;
			case 2: fthresh = 90.0; break;
			case 3: fthresh = 95.0; break;
			case 4: fthresh = 99.0; break;
			case 5: fthresh = 99.5; break;
			case 6: fthresh = 99.9; break;
		}
	}

#ifdef notdef
	fprintf(stderr, "1) variance/resolution tradeoff:\n");
	fprintf(stderr, "   * resolution = %d f_R\n", npi);
	fprintf(stderr, "   * number of tapers = %d\n", nwin);
	fprintf(stderr, "2) null hypothesis: ");
	switch(inoise)
	{
		case 0: fprintf(stderr, "red noise\n"); break;
		case 1: fprintf(stderr, "white noise\n"); break;
		default: fprintf(stderr, "locally-white noise\n"); break;
	}
	if (isignal != 2)
	{
		if (ismooth == 0)
			fprintf(stderr, "   * raw noise background estimation\n");
		else
		{
			fprintf(stderr, "   * robust noise background estimation\n");
			fprintf(stderr, "   * & fsmooth = %f\n", fsmooth);
			if (inoise == 0)
			{
				if (ilog == 0)
					fprintf(stderr, "   (min misfit noise background)\n");
				else
					fprintf(stderr, "   (min log-misfit noise background)\n");
			}
		}
	}
	if (isignal == 0)
		fprintf(stderr, "3) signal assumption: harmonic or quasiperiodic\n");
	else if (isignal == 1)
		fprintf(stderr, "3) signal assumption: quasiperiodic\n");
	else
		fprintf(stderr, "3) signal assumption: harmonic\n");

	fprintf(stderr, "4) spectrum:\n");
	if (ispec == 1)
		fprintf(stderr, "    * adaptive estimate\n");
	else
		fprintf(stderr, "    * high-resolution estimate\n");

	if (inorm == 0)
		fprintf(stderr, "    * no normalization\n");
	else if (inorm == 1)
		fprintf(stderr, "    * normalize by N\n");
	else
		fprintf(stderr, "    * normalize by 1/dt\n");

	if (iresh == 0)
		fprintf(stderr, "    * unreshaped\n");
	else
		fprintf(stderr, "    * reshaped at threshold = %f\n", fthresh);

	fprintf(stderr, "    * frequency range: %f:%f\n", f1, f2);

	fprintf(stderr, "5) signals to reconstruct: %d\n", nsignals);
	fprintf(stderr, "6) display:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
#endif

	spec_(a, &nwin, &npi, &dt, &nscan, &f1, &f2, &inoise, &ismooth, &ilog,
		&fsmooth, &isignal, &ispec, &inorm, &iresh, &ithresh, &nf0, &df,
		specraw0, specmed0, specresh0, harmonic0, specback0,
		ftest, conf, fconf, &nbnd, &whiteraw0, &whiterob0, &rhoraw0,
		&rhorob0, &tauraw0, &taurob0);

#ifdef notdef
	findpeaks_(&flow, &fhigh, &isignal,
		&nf0, &df, &nbnd, specraw0, specmed0, specresh0, harmonic0,
		specback0, ftest,
		conf, fconf,
		fsignal, confsignal, &nsignals);

	if (nsignals > 0)
	{
		fprintf(stderr, "Significant Peaks\n");
		fprintf(stderr, "signal frequency signif%\n");
		for (i=0; i < nsignals; i++)
		{
			if (fsignal[i] < bndwdth)
				fprintf(stderr, "(%2d)   TREND       %2d%\n", i, fsignal[i], (int)((confsignal[i] + 0.5)));
			else
				fprintf(stderr, "(%2d)   %7.4f       %2d%\n", i, fsignal[i], (int)((confsignal[i] + 0.5)));
		}
	}
#endif

	if (nf0 - 1 > nsamples)
	{
		fprintf(stderr, "WARNING: mtm unexpected return size!\n");
		nf0 = nsamples + 1;
	}
	buf[0] = 0.0;
	for (i=0; i < nf0 - 1; i++)
	{
		buf[i+ 1] = specraw0[i];
	}
}

