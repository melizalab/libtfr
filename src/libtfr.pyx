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
from importlib.metadata import version as _version
from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
cimport tfr

import numpy as np

# pdoc renders docstrings as markdown; this selects the Args/Returns dialect
__docformat__ = "google"

__all__ = [
    "mfft",
    "mfft_dpss",
    "mfft_precalc",
    "tfr_spec",
    "hermf",
    "dpss",
    "log_fgrid",
    "fgrid",
    "tgrid",
    "dynamic_range",
    "ITYPE",
    "DTYPE",
    "CTYPE",
]

ctypedef double complex cmplx_t

ITYPE = np.int32
DTYPE = np.double
CTYPE = np.complex128

try:
    __version__ = _version("libtfr")
except _PackageNotFoundError:
    __version__ = "unknown"

cdef class mfft:
    """
    Computes multi-tapered transforms of real signals. Instantiate with factory
    functions.

    """
    cdef tfr.mfft * _mfft

    def __init__(self):
        """Not callable. Use `mfft_dpss` or `mfft_precalc` to construct a transform."""
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
        """The number of tapers used by the transform."""
        return tfr.mtm_ntapers(self._mfft)

    @property
    def nfft(self):
        """The number of points in the transform."""
        return tfr.mtm_nfft(self._mfft)

    @property
    def npoints(self):
        """The number of points in each taper, which may be less than nfft."""
        return tfr.mtm_npoints(self._mfft)

    @property
    def nreal(self):
        """The number of frequency bins in the spectrum of a real signal, nfft/2 + 1."""
        return tfr.mtm_nreal(self._mfft)

    @property
    def tapers(self):
        """A copy of the transform's tapers, dimension (ntapers, npoints)."""
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        cdef double [:, :] arr_view = <double[:ntapers, :npoints]>tfr.mtm_tapers(self._mfft)
        # allocate empty array and copy
        out = np.empty((ntapers, npoints), dtype=DTYPE)
        out[...] = arr_view
        return out

    def tapers_fft(self, double scale):
        """The FFT of the transform's tapers.

        Args:
            scale: positive factor to rescale the tapers by before transforming

        Returns:
            Complex array, dimension (ntapers, nreal).
        """
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t real_count = tfr.mtm_nreal(self._mfft)
        tfr.mtm_tapers_fft(self._mfft, scale)
        out = np.zeros((ntapers, real_count), dtype=CTYPE)
        hc2cmplx(self._mfft, out)
        return out


    def tapers_interpolate(self, double[:] t, double t0, double dt):
        """Interpolate the transform's tapers at specified times.

        The time support of the tapers is given by a start time and a sampling
        interval.

        Args:
            t: times at which to interpolate
            t0: start time of the tapers
            dt: sampling interval of the tapers

        Returns:
            2D array, dimension (ntapers, t.size).
        """
        cdef Py_ssize_t ntimes = t.size
        cdef Py_ssize_t ntapers = tfr.mtm_ntapers(self._mfft)
        cdef Py_ssize_t npoints = tfr.mtm_npoints(self._mfft)
        out = np.zeros((ntapers, ntimes), dtype=DTYPE)
        cdef double[:, :] out_view = out
        tfr.mtm_tapers_interp(self._mfft, &out_view[0, 0], &t[0], ntimes, t0, dt)
        return out

    def mtfft(self, s not None):
        """Compute the complex multitaper FFT of a real-valued signal.

        Args:
            s: input data (1D time series)

        Returns:
            Complex array, dimension (nreal, ntapers).
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
        """Compute the complex multitaper FFT of a point process.

        Args:
            t: event times (1D array)
            dt: sampling interval of the window
            t0: start time of the window

        Returns:
            Complex array, dimension (nreal, ntapers).
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
        """Compute the power spectral density of a signal using multitaper methods.

        Args:
            s: input data (1D time series)
            adapt: with more than one taper, compute the adaptive spectrum

        Returns:
            1D real power spectrum of length nreal. Not normalized; see the
            module docstring for how to scale it.
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
        """Compute the spectrogram of a signal using multitaper methods.

        Args:
            s: input data (1D time series)
            step: number of samples to advance between frames
            adapt: with more than one taper, compute the adaptive spectrum

        Returns:
            Real power spectrogram, dimension (nreal, nframes). Not normalized;
            see the module docstring for how to scale it.
        """
        cdef const double[:] data = np.asarray(s).astype(DTYPE)
        cdef Py_ssize_t nfreq = tfr.mtm_nreal(self._mfft)
        cdef Py_ssize_t nt = tfr.mtm_nframes(self._mfft, data.size, step)
        spec = np.zeros((nt, nfreq), dtype=DTYPE)
        cdef double[:, :] spec_view = spec
        tfr.mtm_spec(self._mfft, &spec_view[0,0], &data[0], data.size, step, adapt)
        return spec.T

    def mtstft(self, s not None, int step):
        """Compute the short-time Fourier transform of a signal using multiple tapers.

        Args:
            s: input data (1D time series)
            step: number of samples to advance between frames

        Returns:
            Complex STFT, dimension (nreal, nframes, ntapers).
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
        """Compute the complex multitaper STFT of a point process.

        Args:
            t: event times (1D array)
            dt: implied sampling interval of the signal, which sets the
                frequency resolution
            step: time interval between frames
            t0: start time of the signal
            tN: stop time of the signal

        Returns:
            A tuple of the complex STFT, dimension (nreal, nframes, ntapers),
            and the number of events in each frame, length nframes.
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
    """Initialize an mfft transform using DPSS tapers.

    This is the standard multitaper transform.

    Args:
        nfft: number of points in the transform
        nw: time-bandwidth parameter
        ntapers: number of tapers to generate
        npoints: number of points in each taper; defaults to nfft

    Returns:
        An `mfft` transform object.
    """
    if npoints <=0:
        npoints = nfft
    cdef mfft instance = mfft.from_ptr(tfr.mtm_init_dpss(nfft, npoints, nw, ntapers))
    return instance


def mfft_precalc(int nfft, tapers not None, weights=None):
    """Copy pre-calculated tapers or window functions into an mfft transform.

    Use this for window functions libtfr does not generate itself, such as a
    Hanning window.

    Args:
        nfft: number of points in the transform
        tapers: taper array, either (npoints,) or (ntapers, npoints)
        weights: weight for each taper; defaults to equal weights

    Returns:
        An `mfft` transform object.

    Raises:
        ValueError: if the number of weights does not match the number of tapers
    """
    cdef int npoints
    cdef int ntapers
    tapers = np.asarray(tapers).astype(DTYPE)
    if tapers.ndim == 1:
        # reshape rather than assigning .shape: in-place shape assignment is
        # deprecated as of numpy 2.5
        tapers = tapers.reshape(1, tapers.size)
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
    """Compute the time-frequency reassignment spectrogram of a signal.

    Args:
        s: input signal (real)
        N: number of frequency points
        step: number of samples to advance between frames
        Np: window size; must be odd and no larger than N
        K: number of tapers to use
        tm: time support of the tapers
        flock: frequency locking parameter. Power is not reassigned further
            than this, in normalized frequency
        tlock: time locking parameter, in frames
        fgrid: output frequency bins, monotonically increasing. Defaults to a
            linear scale with N points, with Nyquist at 1.0

    Returns:
        Power spectrogram, dimension (N/2+1, nframes), or (fgrid.size, nframes)
        if fgrid is given.

    Raises:
        ValueError: if Np is larger than N, or is even
        RuntimeError: if the transform cannot be initialized
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
    """Compute a set of orthogonal Hermite functions.

    These are the tapers used for multi-taper reassigned spectrograms.

    Args:
        N: number of points in the window; must be odd
        M: maximum order of the set of functions
        tm: half-time support

    Returns:
        A tuple of three (M, N) arrays: the Hermite functions, their first
        derivatives, and their time multiples.

    Raises:
        ValueError: if N is even
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
    """Compute discrete prolate spheroidal sequences.

    These are the tapers used for multitaper power spectrum estimation.

    Args:
        N: number of points in the window
        NW: time-bandwidth product. Must be an integer or half-integer;
            typical choices are 2, 5/2, 3, 7/2 or 4
        k: number of DPSS vectors to generate. Must be less than N, and
            vectors beyond NW*2 - 1 are not numerically stable

    Returns:
        A tuple of the tapers, shape (k, N), and their concentration values,
        length k.

    Raises:
        ValueError: if the parameters are invalid
        RuntimeError: if the eigenvalue solver fails
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
    """Generate a logarithmic frequency grid between fmin and fmax.

    Args:
        fmin: first frequency
        fmax: last frequency
        N: number of points
        Fs: sampling frequency. If given, the grid is returned as relative
            frequencies; must be greater than fmax

    Returns:
        1D array of N frequencies.
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
    """Calculate the frequency grid associated with an FFT computation.

    Args:
        Fs: sampling frequency of the data
        nfft: number of points in the FFT
        fpass: lower and upper frequencies of interest, [fmin, fmax), in the
            same units as Fs. Defaults to all frequencies up to Nyquist

    Returns:
        A tuple of the frequencies and the indices of those frequencies within
        the full frequency grid.

    Example:
        With `Fs=1000` and `nfft=1048`, an FFT of a real signal generates 512
        frequencies between 0 and 500 Hz. With `fpass=(0, 100)`, the returned
        indices are those of the frequencies below 100 Hz.

    Adapted from Chronux 1_50.
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
    """Calculate the time grid associated with an STFT.

    Spectrograms generated by libtfr omit frames that extend past the edge of
    the signal. Those frames can only be excluded here if the taper size is
    known, and that size may be adjusted silently by the Hermite and DPSS taper
    functions to keep the number of points even or odd. So the two forms of
    this function differ: given a signal length, the full grid is returned,
    including unsupported frames; given a spectrogram, the grid is truncated to
    match it.

    Args:
        S: length of the signal in samples, or the 2D spectrogram array
            (frequency x time)
        Fs: sampling frequency of the data
        shift: number of samples shifted between frames

    Returns:
        1D array of frame start times.

    Raises:
        ValueError: if S is a 1D array, which is ambiguous
    """
    if isinstance(S, np.ndarray):
        if S.ndim == 1:
            raise ValueError("Input must be a scalar or a 2D numpy array")
        return np.arange(0, 1. / Fs * S.shape[1] * shift, 1. / Fs * shift)
    else:
        return np.arange(0, 1. / Fs * S, 1. / Fs * shift)


def dynamic_range(S not None, double dB):
    """Compress a spectrogram's dynamic range.

    Values more than `dB` below the peak of `S` are clamped to that threshold.
    `S` is on a linear scale.

    Args:
        S: input spectrogram or spectrum
        dB: dynamic range to rescale to

    Returns:
        A copy of S after thresholding.
    """
    cdef double smax = S.max()
    cdef double thresh = 10 ** (np.log10(smax) - dB / 10.)
    return np.where(S >= thresh, S, thresh)
