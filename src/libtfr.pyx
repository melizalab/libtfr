# -*- coding: utf-8 -*-
# -*- mode: cython -*-
#cython: infer_types=True, language_level=3
"""Interface to libtfr spectrogram library using numpy.

Spectrograms are returned as 2D arrays with frequency indexed by row and time by
column. Signals are assumed to be real; therefore real power spectrograms with a
transform size of N have N/2+1 rows, corresponding to frequencies from 0 to
Nyquist. The number of time points in the spectrogram is (M - W + 1)/S, where M
is the length of the signal, S is the shift (S), and the analysis window size is
W (this may be less than or equal to N). Only time points corresponding to
window positions that completely overlap with the signal are returned. Power
spectra and spectrograms are not normalized; they are just the absolute values
of the complex FFT results. To get the power spectrum, divide by the square of
the sum of the window. To get the power spectral density, divide by the sampling
rate times the sum of the square of the window function.

Copyright C Daniel Meliza 2010-2016.  Licensed for use under GNU
General Public License, Version 2.  See COPYING for details.

"""
cimport cython
from cython.view cimport array as cvarray
cimport tfr

import numpy as np

ctypedef double complex cmplx_t

ITYPE = np.int32
DTYPE = np.double
CTYPE = np.complex128

__version__ = "2.1.9"

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

    @staticmethod
    cdef mfft from_ptr(tfr.mfft *ptr):
         if ptr is NULL:
             raise MemoryError
         cdef mfft instance = mfft.__new__(mfft)
         instance._mfft = ptr
         return instance

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
    def nreal(self):
        return tfr.mtm_nreal(self._mfft)

    @property
    def tapers(self):
        """A copy of the transform object's tapers, dimension (ntapers, npoints)"""
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        cdef double [:, :] arr_view = <double[:ntapers, :npoints]>tfr.mtm_tapers(self._mfft)
        # allocate empty array and copy
        out = np.empty((ntapers, npoints), dtype=DTYPE)
        out[...] = arr_view
        return out

    def tapers_fft(self, double scale):
        """The FFT of the transform object's tapers, dimension (ntapers, nreal) """
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t real_count = tfr.mtm_nreal(self._mfft)
        tfr.mtm_tapers_fft(self._mfft, scale)
        out = np.zeros((ntapers, real_count), dtype=CTYPE)
        hc2cmplx(self._mfft, out)
        return out


    def tapers_interpolate(self, double[:] t, double t0, double dt):
        """
        Interpolate the transform object's tapers at specified values. The time
        support for the tapers is specified by a start time and a sampling interval.

        t - input data (array of times)
        t0 - the start time of the tapers
        dt - the sampling interval of the tapers
        returns a 2D array, dimension (ntapers, t.size)
        """
        cdef Py_ssize_t ntimes = t.size
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        out = np.zeros((ntapers, ntimes), dtype=DTYPE)
        cdef double[:, :] out_view = out
        tfr.mtm_tapers_interp(self._mfft, &out_view[0, 0], &t[0], ntimes, t0, dt)
        return out

    def mtfft(self, s not None):
        """
        Computes complex multitaper FFT of a real-valued signal.

        s - input data (1D time series)
        returns array of complex numbers, dimension (nreal, ntapers)
        """
        # this allows the caller to use any kind of array as input
        cdef const double[:] data = np.asarray(s).astype(DTYPE)
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t real_count = tfr.mtm_nreal(self._mfft)
        out = np.empty((ntapers, real_count), dtype=CTYPE)
        tfr.mtfft(self._mfft, &data[0], data.size);
        hc2cmplx(self._mfft, out)
        return out.T

    def mtfft_pt(self, t not None, double dt, double t0):
        """
        Computes complex multitaper FFT of a point process

        times - input data (1D time series)
        t0, dt - define the support of the window
        returns array of complex numbers, dimension (nreal, ntapers)

        """
        # this algorithm could be further cythonized
        times = np.asarray(t).astype(DTYPE)
        cdef Py_ssize_t nevents = times.size
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        cdef Py_ssize_t nfft = tfr.mtm_nfft(self._mfft)
        cdef double Fs = 1 / dt
        # finite size correction
        cdef double Msp = nevents / npoints
        H = self.tapers_fft(1.0).T
        if nevents == 0:
            return np.zeros_like(H)
        # apply tapers to input times
        ht = self.tapers_interpolate(times, t0, dt)
        # exp(2 pi i omega t)
        f, idx = fgrid(Fs, nfft)
        w = -2j * np.pi * f
        Y = np.exp(np.outer(w, times - t0))
        # integrate  exp (2 pi i omega t) * ht with matrix multiplication
        return np.dot(Y, ht.T) - H * Msp


    def mtpsd(self, s not None, adapt=True):
        """Compute PSD of a signal using multitaper methods

        s -  input data (1D time series)
        adapt - if more than one taper, compute adaptive spectrum (default True)

        @returns  N/2+1 1D real power spectrum, not normalized
        """
        cdef const double[:] data = np.asarray(s).astype(DTYPE)
        cdef Py_ssize_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef double sigpow = tfr.mtfft(self._mfft, &data[0], data.size)
        if not adapt:
            sigpow = 0.0
        spec = np.empty(nfreq, dtype=DTYPE)
        cdef double[:] spec_view = spec
        tfr.mtpower(self._mfft, &spec_view[0], sigpow)
        return spec

    def mtspec(self, s not None, int step, adapt=True):
        """Compute spectrogram of a signal using multitaper methods

        s -     input data (1D time series)
        step -  number of samples to step between frames
        adapt - if more than one taper, compute adaptive spectrum (default True)

        @returns real power spectrogram, dim (N/2+1, L), not normalized
        """
        cdef const double[:] data = np.asarray(s).astype(DTYPE)
        cdef Py_ssize_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef Py_ssize_t nt = tfr.mtm_nframes(self._mfft, data.size, step)
        spec = np.zeros((nt, nfreq), dtype=DTYPE)
        cdef double[:, :] spec_view = spec
        tfr.mtm_spec(self._mfft, &spec_view[0,0], &data[0], data.size, step, adapt)
        return spec.T

    def mtstft(self, s not None, int step):
        """Compute STFT of a signal using multiple tapers

        s -     input data (1D time series)
        step -  number of samples to step between frames

        @returns complex STFT, dim (N/2+1, L, ntapers)
        """
        cdef const double[:] data = np.asarray(s).astype(DTYPE)
        cdef Py_ssize_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t nt = tfr.mtm_nframes(self._mfft, data.size, step)
        out = np.empty((nt, ntapers, nfreq), dtype=CTYPE)
        cdef cmplx_t[:, :, :] out_view = out
        cdef Py_ssize_t t
        for t in range(nt):
            tfr.mtfft(self._mfft, &data[t*step], data.size - t*step)
            hc2cmplx(self._mfft, out[t,:,:])
        return out.transpose(2, 0, 1)

    def mtstft_pt(self, t not None, double dt, double step, double t0, double tN):
        """
        Computes complex multitaper STFT of a point process

        times - input data (1D time series)
        dt    - implied sampling rate of the signal (determines frequency resolution)
        step  - time interval between frames
        t0, tN, - start and stop of the signal
        returns array of complex numbers, dimension (nreal, ntapers)
        """
        cdef unsigned int i
        cdef double tw0, Msp
        times = np.asarray(t).astype(DTYPE)
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        cdef Py_ssize_t nfft = tfr.mtm_nfft(self._mfft)
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef double window = npoints * dt
        cdef double Fs = 1 / dt
        cdef Py_ssize_t nframes = int((tN - t0 - window) / step) + 1
        # exp(2 pi i omega t)
        f = fgrid(Fs, nfft)[0]
        w = -2j * np.pi * f
        H = self.tapers_fft(1.0).T
        J = np.zeros((f.size, nframes, ntapers), dtype=CTYPE)
        Nsp = np.zeros(nframes, dtype=ITYPE)
        for i in range(nframes):
            tw0 = t0 + i * step
            # this is not very efficient
            idx = (times >= tw0) & (times < (tw0 + window))
            events = times[idx]
            Nsp[i] = events.size
            if Nsp[i] == 0:
                J[:, i, :] = 0.0
            else:
                Msp = 1. * Nsp[i] / npoints
                ht = self.tapers_interpolate(events, tw0, dt)
                Y = np.exp(np.outer(w, events - tw0))
                J[:, i, :] = np.dot(Y, ht.T) - H * Msp
        return J, Nsp


def mfft_dpss(int nfft, double nw, int ntapers, int npoints=0):
    """
    Initializes a mfft transform using DPSS tapers (i.e. for a standard
    multitaper transform)

    nfft -     number of points in the transform/dpss tapers
    nw -       time-frequency parameter
    ntapers -  number of tapers to generate
    npoints -  number of points in the taper
    """
    if npoints <=0:
        npoints = nfft
    cdef mfft instance = mfft.from_ptr(tfr.mtm_init_dpss(nfft, npoints, nw, ntapers))
    return instance


def mfft_precalc(int nfft, tapers not None, weights=None):
    """
    Copy pre-calculated tapers/window functions (e.g. hanning) into a mtfft

    nfft -    number of points in the transform/dpss tapers
    tapers -  array with tapers, either (npoints,) or (ntapers,npoints)
    weights - array with weights for the tapers, or None to give equal weight
    """
    cdef int npoints
    cdef int ntapers
    tapers = np.asarray(tapers).astype(DTYPE)
    if tapers.ndim == 1:
        tapers.shape = (1, tapers.size)
    ntapers, npoints = tapers.shape
    cdef double[:, :] tapers_view = tapers

    if weights is None:
        weights = np.ones(ntapers, dtype=DTYPE)
    elif weights.size != ntapers:
        raise ValueError("Number of weights does not match number of tapers")
    else:
        weights = np.asarray(weights).astype(DTYPE)
    cdef double[:] weights_view = weights

    cdef mfft instance = mfft.from_ptr(tfr.mtm_init(nfft, npoints, ntapers))
    tfr.mtm_copy(instance._mfft, &tapers_view[0,0], &weights_view[0])
    return instance


def tfr_spec(s not None, int N, int step, int Np, int K=6,
             double tm=6.0, double flock=0.01, int tlock=5, fgrid=None):
    """
    Compute time-frequency reassignment spectrogram of input signal s

    s - input signal (real)
    N - number of frequency points
    step - step size (in time points)
    Np - window size (should be <= N, must be odd)
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

    if Np > N:
        raise ValueError("Np must be less or equal to N")
    if Np % 2 == 0:
        raise ValueError("Np must be odd")

    # coerce data to proper type
    cdef double[:] samples = np.asarray(s).astype(DTYPE)

    # generate/convert frequency grid
    cdef int nfreq = N//2 + 1
    cdef double[:] fgrid_view
    cdef double * fgridp = NULL
    if fgrid is not None:
        fgrid = np.asarray(fgrid).astype(DTYPE)
        fgrid_view = fgrid
        fgridp = &fgrid_view[0]
        nfreq = fgrid.size

    # initialize transform
    cdef tfr.mfft * mtmh = tfr.mtm_init_herm(N, Np, K, tm)
    if mtmh is NULL:
        raise RuntimeError(f"Arguments ({N}, {Np}, {K}, {tm}) rejected by mtm_init_herm or other error")

    # allocate output array
    cdef Py_ssize_t nt = tfr.mtm_nframes(mtmh, samples.size, step)
    out = np.zeros((nt, nfreq), dtype=DTYPE)
    cdef double[:, :] out_view = out
    tfr.tfr_spec(mtmh, &out_view[0,0], &samples[0], samples.size, -1,
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
    if N % 2 == 0:
        raise ValueError("N must be odd")
    h = np.empty((M, N), dtype=DTYPE)
    cdef double[:, :] h_view = h
    Dh = np.empty((M, N), dtype=DTYPE)
    cdef double[:, :] Dh_view = Dh
    Th = np.empty((M, N), dtype=DTYPE)
    cdef double[:, :] Th_view = Th
    tfr.hermf(N, M, tm, &h_view[0,0], &Dh_view[0,0], &Th_view[0,0])
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
    tapers = np.empty((k, N), dtype=DTYPE)
    cdef double[:, :] tapers_view = tapers
    lambdas = np.empty(k, dtype=DTYPE)
    cdef double[:] lambdas_view = lambdas

    rv = tfr.dpss(&tapers_view[0,0], &lambdas_view[0], N, NW, k)
    if rv == 0:
        return tapers, lambdas
    elif rv == -1:
        raise ValueError("Invalid DPSS parameters")
    elif rv == -2:
        raise RuntimeError("Eigenvalue solver failed")
    else:
        raise RuntimeError("Unknown error")


cdef void hc2cmplx(tfr.mfft * mtm, cmplx_t[:,:] out) noexcept nogil:
    """Copy data from workspace of mtm object into a complex array"""
    cdef Py_ssize_t nfft = tfr.mtm_nfft(mtm)
    cdef Py_ssize_t ntapers = tfr.mtm_ntapers(mtm)
    cdef Py_ssize_t real_count = nfft // 2 + 1
    cdef Py_ssize_t imag_count = (nfft + 1) // 2
    cdef Py_ssize_t t, n
    cdef const double * buf = tfr.mtm_buffer(mtm)

    with cython.boundscheck(False), cython.wraparound(False):
        for t in range(ntapers):
            for n in range(0, real_count):
                out[t, n] = buf[t*nfft+n]
            for n in range(1, imag_count):
                out[t, n] = out[t, n] + buf[t*nfft+(nfft-n)] * 1j


### Utility functions: not much benefit to cython as written but nice to have in
### the same module
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


def fgrid(double Fs, int nfft, fpass=None):
    """
    Calculate the frequency grid associated with an fft computation

    @param Fs        sampling frequency associated with the data
    @param nfft      number of points in fft
    @param fpass     upper and lower frequencies of interest,
                     [fmin fmax), in same units as Fs. If not supplied,
                     returns all the frequencies up to Nyquist.

    @returns (1D array of frequencies,
              1D array indexing frequencies in the full frequency grid)

    Example:

    If Fs=1000, and nfft=1048, an fft calculation of a real signal
    generates 512 frequencies between 0 and 500 (i.e. Fs/2) Hz. Now if
    fpass=(0,100), findx will contain the indices in the frequency grid
    corresponding to frequencies < 100 Hz.

    From Chronux 1_50
    """
    cdef double df = Fs / nfft
    f = np.arange(0, Fs, df)  # all possible frequencies

    if fpass is not None:
        f1, f2 = fpass
        findx = ((f >= f1) & (f < f2)).nonzero()[0]
    else:
        findx = (f <= Fs / 2).nonzero()[0]

    return f[findx], findx


def tgrid(S not None, double Fs, int shift):
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
    if isinstance(S, np.ndarray):
        if S.ndim == 1: raise ValueError, "Input must be a scalar or a 2D numpy array"
        return np.arange(0, 1. / Fs * S.shape[1] * shift, 1. / Fs * shift)
    else:
        return np.arange(0, 1. / Fs * S, 1. / Fs * shift)


def dynamic_range(S not None, double dB):
    """
    Compress a spectrogram's dynamic range by thresholding all values
    dB less than the peak of S (linear scale).

    @param S    the input spectrogram or spectrum
    @param dB   the dynamic range to rescale to

    @returns copy of S after thresholding
    """
    cdef double smax = S.max()
    cdef double thresh = 10 ** (np.log10(smax) - dB / 10.)
    return np.where(S >= thresh, S, thresh)
