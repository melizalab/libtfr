
cdef extern from "tfr.h":
    ctypedef struct mfft:
        pass

    int dpss(double *, double *, int, double, int) nogil
    int hermf(int N, int M, double tm, double *h, double *Dh, double *Th) nogil

    mfft * mtm_init(int nfft, int npoints, int ntapers)
    mfft * mtm_init_dpss(int nfft, int npoints, double nw, int ntapers)
    mfft * mtm_init_herm(int nfft, int npoints, int order, double tm)
    void mtm_copy(mfft * mtm, const double *, const double *)
    void mtm_destroy(mfft * mtm)

    int mtm_nfft(const mfft * mtm) nogil
    int mtm_npoints(const mfft * mtm) nogil
    int mtm_ntapers(const mfft * mtm) nogil
    int mtm_nreal(const mfft * mtm) nogil
    int mtm_nframes(const mfft * mtm, int signal_size, int step_size) nogil
    const double * mtm_buffer(const mfft * mtm) nogil
    const double * mtm_tapers(const mfft * mtm) nogil

    double mtfft(mfft * mtm, const double * data, int nbins) nogil
    void mtpower(const mfft * mtm, double * pow, double sigpow) nogil

    void mtm_spec(mfft * mtm, double *spec, const double *samples, int nsamples,
                  int shift, int adapt) nogil


    void tfr_spec(mfft * mtm, double * spec, const double *samples, int nsamples,
                  int k, int shift, double flock, int tlock, int nfreq, const double *fgrid) nogil
