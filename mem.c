#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sonogram.h"

#ifndef SQR
#define SQR(a) ( (a) * (a) )
#endif

memspect *init_memspect(int n, int m)
{
	memspect *ms;
	int j;

	if ((ms = (memspect *)malloc(sizeof(memspect))) != NULL)
	{
		ms->n = n;
		ms->m = m;
		ms->wk1 = (double *)malloc((n+1) * sizeof(double));
		ms->wk2 = (double *)malloc((n+1) * sizeof(double));
		ms->wkm = (double *)malloc((m+1) * sizeof(double));
		ms->cof = (double *)malloc((m+1) * sizeof(double));
		ms->wprtheta = (double *)malloc(n * sizeof(double));
		ms->wpitheta = (double *)malloc(n * sizeof(double));
		for (j=0; j < n / 2; j++)
		{
			ms->wprtheta[j] = cos(M_PI * 2.0 * ((double)j / (double)n));
			ms->wpitheta[j] = sin(M_PI * 2.0 * ((double)j / (double)n));
		}
		ms->fmin = 0.0;
		ms->fmax = 0.5;
		ms->theta_array_len = n / 2;
	}
	return ms;
}

void free_memspect(memspect *ms)
{
	free(ms->wk1);
	free(ms->wk2);
	free(ms->wkm);
	free(ms->cof);
	free(ms->wprtheta);
	free(ms->wpitheta);
	free(ms);
}

void mem_calc_coefficients(memspect *ms, float *data)
{
	int i, j, k;
	double p = 0.0, *wk1, *wk2, *wkm;
	int n, m;
	double *pm, *cof;

	wk1 = ms->wk1;
	wk2 = ms->wk2;
	wkm = ms->wkm;
	cof = ms->cof;
	pm = &(ms->pm);
	n = ms->n;
	m = ms->m;

	memset((void *)wk1, 0, (n + 1) * sizeof(double));
	memset((void *)wk2, 0, (n + 1) * sizeof(double));
	memset((void *)wkm, 0, (m + 1) * sizeof(double));
	memset((void *)cof, 0, (m + 1) * sizeof(double));
	ms->pm = 0.0;

	for (j=1; j <= n; j++) p += SQR(data[j-1]);
	*pm = p / n;

	wk1[1] = data[0];
	wk2[n-1] = data[n-1];
	for (j=2; j <= n-1; j++)
	{
		wk1[j] = data[j-1];
		wk2[j-1] = data[j-1];
	}

	for (k=1; k <= m; k++)
	{
		double num = 0.0, denom = 0.0;

		for (j=1; j <= (n-k); j++)
		{
			num += wk1[j] * wk2[j];
			denom += SQR(wk1[j]) + SQR(wk2[j]);
		}
		cof[k] = 2.0 * num / denom;
		*pm *= (1.0 - SQR(cof[k]));
		for (i=1; i <= (k-1); i++)
			cof[i] = wkm[i] - cof[k] * wkm[k-i];

		if (k==m)
			break;

		for (i=1; i <= k; i++) wkm[i] = cof[i];
		for (j=1; j <= (n-k-1); j++)
		{
			wk1[j] -= wkm[k] * wk2[j];
			wk2[j] = wk2[j+1] - wkm[k] * wk1[j+1];
		}
	}

	return;
}

double mem_eval_spectrum(memspect *ms, double fdt)
{
	int m;
	double *cof, pm;
	int i;
	double sumr = 1.0;
	double sumi = 0.0;
	double wr = 1.0;
	double wi = 0.0;
	double wpr, wpi, wtemp, theta;

	m = ms->m;
	cof = ms->cof;
	pm = ms->pm;

	theta = M_PI * 2 * fdt;
	wpr = cos(theta);
	wpi = sin(theta);

	for (i=1; i <= m; i++)
	{
		wr = (wtemp = wr) * wpr - wi * wpi;
		wi = wi * wpr + wtemp * wpi;
		sumr -= cof[i] * wr;
		sumi -= cof[i] * wi;
	}
	return pm / (SQR(sumr) + SQR(sumi));
}


void mem_spectrum(memspect *ms, float *spectrum)
{
	int m, n;
	double *cof, pm;
	int i, j;
	double sumr = 1.0;
	double sumi = 0.0;
	double wr = 1.0;
	double wi = 0.0;
	double wpr, wpi, wtemp;

	m = ms->m;
	n = ms->n;
	cof = ms->cof;
	pm = ms->pm;

	for (j=0; j < n / 2; j++)
	{
		sumr = wr = 1.0;
		sumi = wi = 0.0;
		wpr = ms->wprtheta[j];
		wpi = ms->wpitheta[j];
		for (i=1; i <= m; i++)
		{
			wr = (wtemp = wr) * wpr - wi * wpi;
			wi = wi * wpr + wtemp * wpi;
			sumr -= cof[i] * wr;
			sumi -= cof[i] * wi;
		}
		spectrum[j] = pm / (SQR(sumr) + SQR(sumi));
	}
}


void mem_spectrum2(memspect *ms, float *spectrum, int spectrum_len, float fmin, float fmax)
{
	int m, n;
	double *cof, pm;
	int i, j;
	double sumr = 1.0;
	double sumi = 0.0;
	double wr = 1.0;
	double wi = 0.0;
	double wpr, wpi, wtemp;
	double fdt, step;

	m = ms->m;
	n = ms->n;
	cof = ms->cof;
	pm = ms->pm;

	if ((ms->theta_array_len != spectrum_len) || (ms->fmin != fmin) || (ms->fmax != fmax))
	{
		ms->fmin = fmin;
		ms->fmax = fmax;
		ms->wprtheta = (double *)realloc(ms->wprtheta, spectrum_len * sizeof(double));
		ms->wpitheta = (double *)realloc(ms->wpitheta, spectrum_len * sizeof(double));
		ms->theta_array_len = spectrum_len;
		step = (fmax - fmin) / (spectrum_len);
		for (fdt = fmin, j=0; j < spectrum_len; j++, fdt += step)
		{
			ms->wprtheta[j] = cos(M_PI * 2.0 * fdt);
			ms->wpitheta[j] = sin(M_PI * 2.0 * fdt);
		}
	}

	for (j=0; j < spectrum_len; j++)
	{
		sumr = wr = 1.0;
		sumi = wi = 0.0;
		wpr = ms->wprtheta[j];
		wpi = ms->wpitheta[j];
		for (i=1; i <= m; i++)
		{
			wr = (wtemp = wr) * wpr - wi * wpi;
			wi = wi * wpr + wtemp * wpi;
			sumr -= cof[i] * wr;
			sumi -= cof[i] * wi;
		}
		spectrum[j] = pm / (SQR(sumr) + SQR(sumi));
	}
}

