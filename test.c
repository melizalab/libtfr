
#include <stdlib.h>
#include <stdio.h>
#include <pcmio.h>
#include "sonogram.h"

#define MIN(a, b) ( ((a) < (b)) ? (a) : (b) )

int
main(int argv, char **argc)
{
	short *data;
	int nsamples;
	PCMFILE *fp = pcm_open("A7.pcm","r");
	int i =0;
	float *psd;

	pcm_read(fp, &data, &nsamples);
	printf("Number of samples: %d\n", nsamples);

	struct sonogram *sono = create_sonogram();

	printf("Setting options\n");
	sonogram_setopts(sono, SONO_OPT_FFTSIZE, 256);
	sonogram_setopts(sono, SONO_OPT_OVERLAP, 50);
	sonogram_setopts(sono, SONO_OPT_NSAMPLES, nsamples);
	sonogram_setopts(sono, SONO_OPT_METHOD, METHOD_STFT);
	sonogram_setopts(sono, SONO_OPT_WINDOW, WINDOW_HANNING);
	
	sonogram_setopts(sono, SONO_OPT_WINDOW, WINDOW_MULTITAPER);
	sonogram_setopts(sono, SONO_OPT_MTM_NW, 3.5);
	sonogram_setopts(sono, SONO_OPT_MTM_NTAPERS, 5);
	sonogram_setopts(sono, SONO_OPT_MTM_ADAPT, 0);
	printf("Set options\n");
	
	FILE *out = fopen("testsono.bin","w");

	double sigpow = 0.0;
	for (i = 0; i < nsamples; i++)
		sigpow += (double)data[i] * data[i];
	sigpow /= nsamples;

	for (i = 0; i < 100; i++) {
		psd = calculate_psd_cached_column(sono, data, nsamples, i, 0, 0.5);
		fwrite(psd, sizeof(float), 129, out);
	}

	fclose(out);
	pcm_close(fp);
	return (0);
}
	
