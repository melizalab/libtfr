# -*- coding: utf-8 -*-
# -*- mode: cython -*-
"""Interface to libtfr spectrogram library using numpy.

Spectrograms are returned as 2D arrays with frequency indexed by row and time by
column. Signals are assumed to be real; therefore real power spectrograms with a
transform size of N have N/2+1 rows, corresponding to frequencies from 0 to
Nyquist. The number of time points in the spectrogram is (M - W + 1)/S, where M
is the length of the signal, S is the shift (S), and the analysis window size is
W (this may be less than or equal to N). Only time points corresponding to
window positions that completely overlap with the signal are returned.

Copyright C Daniel Meliza 2010-2016.  Licensed for use under GNU
General Public License, Version 2.  See COPYING for details.

"""
from cython cimport view, boundscheck
import numpy as nx
cimport numpy as nx
nx.import_array()
cimport tfr

DTYPE = nx.double
ctypedef nx.double_t DTYPE_t
CTYPE = nx.complex128
ctypedef nx.complex128_t CTYPE_t


cdef class mfft:
    """
    Computes multi-tapered transforms of real signals. Instantiate with factory
    functions.

    """
    cdef tfr.mfft * _mfft
    def __init__(self):
        raise TypeError("This class cannot be directly instantiated")
    def __dealloc__(self):
        if self._mfft is not NULL:
            tfr.mtm_destroy(self._mfft)

    @property
    def ntapers(self):
        return tfr.mtm_ntapers(self._mfft)
    @property
    def nfft(self):
        return tfr.mtm_nfft(self._mfft)
    @property
    def npoints(self):
        return tfr.mtm_npoints(self._mfft)
    @property
    def tapers(self):
        """A copy of the transform object's tapers, dimension (ntapers, npoints)"""
        cdef nx.npy_intp dims[2]
        cdef nx.ndarray out
        dims[0] = tfr.mtm_ntapers(self._mfft)
        dims[1] = tfr.mtm_npoints(self._mfft)
        out = nx.PyArray_SimpleNewFromData(2, dims,
                                           nx.NPY_DOUBLE, tfr.mtm_tapers(self._mfft))
        return out.copy()


    def mtfft(self, s not None):
        """
        Computes complex multitaper FFT of a real-valued signal.

        s - input data (1D time series)
        returns array of complex numbers, dimension (nfft, ntapers)
        """
        cdef double[:] data = nx.asarray(s).astype(DTYPE)
        cdef size_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef size_t real_count = tfr.mtm_nreal(self._mfft)
        cdef nx.ndarray[CTYPE_t, ndim=2] out = nx.zeros((ntapers, real_count), dtype=CTYPE)
        tfr.mtfft(self._mfft, &data[0], data.size);
        hc2cmplx(self._mfft, out)
        return out.T

    def mtpsd(self, s not None, adapt=True):
        """Compute PSD of a signal using multitaper methods

        s -  input data (1D time series)
        adapt - compute adaptive spectrum (default True)

        @returns  N/2+1 1D real power spectrum density
        """
        cdef double[:] data = nx.asarray(s).astype(DTYPE)
        cdef size_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef nx.ndarray[DTYPE_t, ndim=1] spec = nx.zeros(nfreq, dtype=DTYPE)
        cdef double sigpow = tfr.mtfft(self._mfft, &data[0], data.size)
        if not adapt:
            sigpow = 0.0
        tfr.mtpower(self._mfft, &spec[0], sigpow)
        return spec

    def mtspec(self, s not None, int step, adapt=True):
        """Compute spectrogram of a signal using multitaper methods

        s -     input data (1D time series)
        step -  number of samples to step between frames
        adapt - compute adaptive spectrum (default True)

        @returns real power spectrogram, dim (N/2+1, L)
        """
        cdef double[:] data = nx.asarray(s).astype(DTYPE)
        cdef size_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef size_t nt = tfr.mtm_nframes(self._mfft, data.size, step)
        cdef nx.ndarray[DTYPE_t, ndim=2] spec = nx.zeros((nt, nfreq), dtype=DTYPE)
        tfr.mtm_spec(self._mfft, &spec[0,0], &data[0], data.size, step, adapt)
        return spec.T

    def mtstft(self, s not None, int step):
        """Compute STFT of a signal using multiple tapers

        s -     input data (1D time series)
        step -  number of samples to step between frames

        @returns complex STFT, dim (N/2+1, L, ntapers)
        """
        cdef double[:] data = nx.asarray(s).astype(DTYPE)
        cdef size_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef size_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef size_t nt = tfr.mtm_nframes(self._mfft, data.size, step)
        cdef nx.ndarray[CTYPE_t, ndim=3] out = nx.zeros((nt, ntapers, nfreq), dtype=CTYPE)
        cdef size_t t
        for t in range(nt):
            tfr.mtfft(self._mfft, &data[t*step], data.size - t*step)
            hc2cmplx(self._mfft, out[t,:,:])
        return out.transpose(2, 0, 1)


def mfft_dpss(int nfft, double nw, int ntapers, int npoints):
    """
    Initializes a mfft transform using DPSS tapers (i.e. for a standard
    multitaper transform)

    nfft -     number of points in the transform/dpss tapers
    nw -       time-frequency parameter
    ntapers -  number of tapers to generate
    """
    cdef mfft instance = mfft.__new__(mfft)
    instance._mfft = tfr.mtm_init_dpss(nfft, npoints, nw, ntapers)
    return instance


def mfft_precalc(int nfft, tapers not None, weights=None):
    """
    Copy pre-calculated tapers/window functions (e.g. hanning) into a mtfft

    nfft -    number of points in the transform/dpss tapers
    tapers -  array with tapers, either (npoints,) or (ntapers,npoints)
    weights - array with weights for the tapers, or None to give equal weight
    """
    cdef mfft instance = mfft.__new__(mfft)
    cdef int npoints
    cdef int ntapers
    tapers = nx.asarray(tapers)
    if tapers.ndim == 1:
        tapers.shape = (1, tapers.size)
    ntapers,npoints = tapers.shape
    cdef nx.ndarray[DTYPE_t, ndim=2] tapers_c = tapers.astype(DTYPE)

    cdef nx.ndarray[DTYPE_t, ndim=1] weights_c
    if weights is None:
        weights_c = nx.ones(ntapers, dtype=DTYPE)
    elif weights.size != ntapers:
        raise ValueError("Number of weights does not match number of tapers")
    else:
        weights_c = nx.asarray(weights).astype(DTYPE)

    instance._mfft = tfr.mtm_init(nfft, npoints, ntapers)
    tfr.mtm_copy(instance._mfft, &tapers_c[0,0], &weights_c[0])
    return instance


def tfr_spec(s not None, int N, int step, int Np, int K=6,
             double tm=6.0, double flock=0.01, int tlock=5, fgrid=None):
    """
    Compute time-frequency reassignment spectrogram of input signal s

    s - input signal (real)
    N - number of frequency points
    step - step size (in time points)
    Np - window size (should be <= N)
    K - number of tapers to use (default 6)
    tm - time support of tapers (default 6.0)
    flock - frequency locking parameter; power is not reassigned
            more than this value (normalized frequency; default 0.01)
    tlock - time locking parameter (in frames; default 5)
    fgrid - output frequency bins: monotonically increasing
            (default linear scale with N points; Nyquist is 1.0)

    returns an N/2+1 by L power spectrogram, or if fgrid is specified,
    fgrid.size by L
    """
    # coerce data to proper type
    cdef double[:] samples = nx.asarray(s).astype(DTYPE)

    # generate/convert frequency grid
    cdef int nfreq = N/2 + 1
    cdef nx.ndarray[DTYPE_t, ndim=1] fgrid_cast
    cdef double * fgridp = NULL
    if fgrid is not None:
        fgrid_cast = nx.asarray(fgrid).astype(DTYPE)
        fgridp = &fgrid_cast[0]
        nfreq = fgrid_cast.size

    # initialize transform
    cdef tfr.mfft * mtmh = tfr.mtm_init_herm(N, Np, K, tm)

    # allocate output array
    cdef size_t nt = tfr.mtm_nframes(mtmh, samples.size, step)
    cdef nx.ndarray[DTYPE_t, ndim=2] out = nx.zeros((nt, nfreq), dtype=DTYPE)

    tfr.tfr_spec(mtmh, &out[0,0], &samples[0], samples.size, -1,
                 step, flock, tlock, nfreq, fgridp);
    tfr.mtm_destroy(mtmh)

    return out.T


def hermf(int N, int M=6, double tm=6.0):
    """
    Computes a set of orthogonal Hermite functions for use in computing
    multi-taper reassigned spectrograms

    @param N      the number of points in the window (must be odd)
    @param M      the maximum order of the set of functions (default 6)
    @param tm     half-time support (default 6)

    @returns  hermite functions (MxN), first derivative of h (MxN), time-multiple of h (MxN)
    """
    cdef nx.ndarray[DTYPE_t, ndim=2] h = nx.empty((M, N), dtype=DTYPE)
    cdef nx.ndarray[DTYPE_t, ndim=2] Dh = nx.empty((M, N), dtype=DTYPE)
    cdef nx.ndarray[DTYPE_t, ndim=2] Th = nx.empty((M, N), dtype=DTYPE)
    tfr.hermf(N, M, tm, &h[0,0], &Dh[0,0], &Th[0,0])
    return (h, Dh, Th)


def dpss(int N, double NW, int k):
    """
    Computes the discrete prolate spherical sequences used in the
    multitaper method power spectrum calculations.

    @param npoints   the number of points in the window
    @param mtm_p     the time-bandwidth product. Must be an integer or half-integer
                     (typical choices are 2, 5/2, 3, 7/2, or 4)
    @param k         the number of DPSS vectors to generate. Must be less than
                     npoints, but k > NW*2 - 1 are not stable

    @returns (2D array of tapers, shape (k,npoints),
              1D array of concentration values, length k)
    """
    cdef nx.ndarray[DTYPE_t, ndim=2] tapers = nx.empty((k, N), dtype=DTYPE)
    cdef nx.ndarray[DTYPE_t, ndim=1] lambdas = nx.empty(k, dtype=DTYPE)

    rv = tfr.dpss(&tapers[0,0], &lambdas[0], N, NW, k)
    if rv == 0:
        return tapers, lambdas
    elif rv == -1:
        raise ValueError("Invalid DPSS parameters")
    elif rv == -2:
        raise RuntimeError("Eigenvalue solver failed")
    else:
        raise RuntimeError("Unknown error")


@boundscheck(False)
cdef void hc2cmplx(tfr.mfft * mtm, CTYPE_t[:,:] out) nogil:
    """Copy data from workspace of mtm object into a complex array"""
    cdef size_t nfft = tfr.mtm_nfft(mtm)
    cdef size_t ntapers = tfr.mtm_ntapers(mtm)
    cdef size_t real_count = nfft / 2 + 1
    cdef size_t imag_count = (nfft + 1) / 2
    cdef size_t t, n
    cdef double * buf = tfr.mtm_buffer(mtm)

    for t in range(ntapers):
        for n in range(0, real_count):
            out[t, n] = buf[t*nfft+n]
        for n in range(1, imag_count):
            out[t, n] += buf[t*nfft+(nfft-n)] * 1j


### Utility functions
def log_fgrid(double fmin, double fmax, int N, Fs=None):
    """
    Generates a logarithmic frequency grid between fmin and fmax

    @param fmin  first frequency
    @param fmax  last frequency
    @param N     number of points
    @param base  log base
    @param Fs    set to a positive value to convert values to relative frequencies
    """
    from numpy import log, logspace, e
    lfmin, lfmax = log((fmin, fmax))
    out = logspace(lfmin, lfmax, N, base=e)
    if Fs is not None:
        assert Fs > fmax, "Fs must be greater than fmax"
        return out / Fs
    else:
        return out


# utility functions
def fgrid(double Fs, int nfft, fpass):
    """
    Calculate the frequency grid associated with an fft computation

    @param Fs        sampling frequency associated with the data
    @param nfft      number of points in fft
    @param fpass     upper and lower frequencies of interest,
                     [fmin fmax), in same units as Fs

    @returns (1D array of frequencies,
              1D array indexing frequencies in the full frequency grid)

    Example:

    If Fs=1000, and nfft=1048, an fft calculation of a real signal
    generates 512 frequencies between 0 and 500 (i.e. Fs/2) Hz. Now if
    fpass=(0,100), findx will contain the indices in the frequency grid
    corresponding to frequencies < 100 Hz.

    From Chronux 1_50
    """
    from numpy import arange, abs

    cdef double df = Fs / nfft
    f = arange(0, Fs, df)  # all possible frequencies

    if len(fpass) != 1:
        findx = ((f >= fpass[0]) & (f < fpass[-1])).nonzero()[0]
    else:
        findx = abs(f - fpass).argmin()

    return f[findx], findx


def tgrid(S, double Fs, int shift):
    """
    Calculate the time grid associated with an STFT. Note that
    spectrograms generated by libtfr do not include frames that extend
    beyond the edge of the signal.  These can be excluded if one knows
    the size of the analysis taper(s), but this size may be adjusted
    silently by the functions that generate Hermitian and DPSS tapers
    in order to ensure that the number of points is even or odd.  This
    function will accept either a scalar or an ndarray as its first
    argument: in the former case the full time grid (including
    unsupported frames) will be returned; in the latter, the time grid
    will be truncated to match the length of the spectrogram.

    @param S         length of signal (in samples) OR the 2D spectrogram array (freq x time)
    @param Fs        sampling frequency associated with the data
    @param shift     number of samples shifted between data frames

    @returns a 1D array of frame start times
    """
    from numpy import arange, ndarray
    if isinstance(S, ndarray):
        if S.ndim == 1: raise ValueError, "Input must be a scalar or a 2D numpy array"
        return arange(0, 1. / Fs * S.shape[1] * shift, 1. / Fs * shift)
    else:
        return arange(0, 1. / Fs * S, 1. / Fs * shift)


def dynamic_range(nx.ndarray S, double dB):
    """
    Compress a spectrogram's dynamic range by thresholding all values
    dB less than the peak of S (linear scale).

    @param S    the input spectrogram or spectrum
    @param dB   the dynamic range to rescale to

    @returns copy of S after thresholding
    """
    from numpy import log10, where
    cdef double smax = S.max()
    cdef double thresh = 10 ** (log10(smax) - dB / 10.)
    return where(S >= thresh, S, thresh)
