
cdef extern from "tfr.h":
    ctypedef struct mfft:
        pass

    int dpss(double *, double *, int, double, int)
    int hermf(int N, int M, double tm, double *h, double *Dh, double *Th)

    mfft * mtm_init(int nfft, int npoints, int ntapers)
    mfft * mtm_init_dpss(int nfft, double nw, int ntapers)
    mfft * mtm_init_herm(int nfft, int npoints, int order, double tm)
    void mtm_copy(mfft * mtm, const double *, const double *)
    void mtm_destroy(mfft * mtm)

    int mtm_nfft(const mfft * mtm)
    int mtm_npoints(const mfft * mtm)
    int mtm_ntapers(const mfft * mtm)
    int mtm_nframes(const mfft * mtm, int signal_size, int step_size)
    const double * mtm_buffer(const mfft * mtm)

    double mtfft(mfft * mtm, const double *data, int nbins)

    void tfr_spec(mfft * mtm, double * spec, const double *samples, int nsamples,
                  int k, int shift, double flock, int tlock, int nfreq, const double *fgrid)
