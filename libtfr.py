 # -*- coding: utf-8 -*-
"""
Interface to libtfr spectrogram library using numpy.

Spectrograms are returned as 2D arrays with frequency indexed by row
and time by column.  Signals are assumed to be real; therefore real
power spectrograms with a transform size of N have N/2+1 rows,
corresponding to frequencies from 0 to Nyquist.  Complex spectrograms
have N rows.  The number of time points in the spectrogram is (M - W +
1)/S, where M is the length of the signal, S is the shift (S), and the
analysis window size is W (this may be less than or equal to N).  Only
time points corresponding to window positions that completely overlap
with the signal are returned.

Copyright C Daniel Meliza 2010-2012.  Licensed for use under GNU
General Public License, Version 2.  See COPYING for details.
"""

import _libtfr
__version__ = _libtfr.__version__


def tfr_spec(s, N, step, Np, K=6, tm=6.0, flock=0.01, tlock=5, fgrid=None):
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
    return _libtfr.tfr_spec(s, N, step, Np, K, tm, flock, tlock, fgrid)


def stft(s, w, step, N=0, complex=False):
    """
    Compute short time Fourier transforms with arbitrary window functions

    s - input signal (real)
    w - window function, N-vector or KxN array (for multiple tapers)
    step - step size (in samples)
    N - size of FFT transform (defaults to size of window)
    complex - if True, return complex results (default False)

    returns N/2+1 by L spectrogram (or N by L by K for complex)
    """
    S = _libtfr.stft(s, w, step, N, complex)
    # reorder axes here
    return S.swapaxes(1, 2) if complex else S


def hermf(N, order=6, tm=6.0):
    """
    Computes a set of orthogonal Hermite functions for use in computing
    multi-taper reassigned spectrograms

    @param N      the number of points in the window (must be odd)
    @param order  the maximum order of the set of functions (default 6)
    @param tm     half-time support (default 6)

    @returns  hermite functions (MxN), first derivative of h (MxN), time-multiple of h (MxN)
    """
    return _libtfr.hermf(N, order, tm)


if hasattr(_libtfr, 'mtm_spec'):
    def mtm_spec(s, N, step, NW, k=0, adapt=True):
        """Compute multitaper spectrogram using DPSS tapers

        @param s  input signal
        @param N  number of frequency points (i.e. window size)
        @param step  step size
        @param NW time-frequency product
        @param k  number of tapers (default NW*2-1)
        @param adapt - compute adaptive spectrogram (default True)

        @returns N/2+1 by L power spectrogram
        """
        return _libtfr.mtm_spec(s, N, step, NW, k, adapt)


if hasattr(_libtfr, 'mtm_psd'):
    def mtm_psd(s, NW, k=0, adapt=True):
        """Compute PSD of a signal using multitaper methods

        @param s  input signal
        @param NW  time-frequency product
        @param k  number of tapers (default NW*2-1)
        @param adapt - compute adaptive spectrum (default True)

        @returns  N/2+1 1D real power spectrum density
        """
        return _libtfr.mtm_psd(s, NW, k, adapt)


if hasattr(_libtfr, 'mtfft'):
    def mtfft(s, NW, k=0, N=0):
        """Compute multitaper transform of a signal

        @param s input signal
        @param NW time-frequency product
        @param k number of tapers (default NW*2-1)
        @param N number of points in FFT (default s.size)

        @returns 2D complex array, dimension N x k
        """
        return _libtfr.mtfft(s, NW, k, N)


if hasattr(_libtfr, 'dpss'):
    def dpss(N, NW, k):
        """
        Computes the discrete prolate spherical sequences used in the
        multitaper method power spectrum calculations.

        @param npoints   the number of points in the window
        @param mtm_p     the time-bandwidth product. Must be an integer or half-integer
                         (typical choices are 2, 5/2, 3, 7/2, or 4)
        @param k         If a scalar, returns the 1:k DPSS vectors
                         If a 2-ple, returns the k[0]:k[1] DPSS vectors
                         Default is to return all vectors

        @returns (2D array of eigenvectors, shape (n,npoints),
                  2D array of eigenvalues, length n = (mtm_p * 2 - 1))
        """
        return _libtfr.dpss(N, NW, k)


def log_fgrid(fmin, fmax, N, Fs=None):
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
    if Fs:
        assert Fs > fmax, "Fs must be greater than fmax"
        return out / Fs
    else:
        return out


# utility functions
def fgrid(Fs, nfft, fpass):
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

    df = float(Fs) / nfft
    f = arange(0, Fs, df)  # all possible frequencies

    if len(fpass) != 1:
        findx = ((f >= fpass[0]) & (f < fpass[-1])).nonzero()[0]
    else:
        findx = abs(f - fpass).argmin()

    return f[findx], findx


def tgrid(S, Fs, shift):
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


def dynamic_range(S, dB):
    """
    Compress a spectrogram's dynamic range by thresholding all values
    dB less than the peak of S (linear scale).

    @param S    the input spectrogram or spectrum
    @param dB   the dynamic range to rescale to

    @returns copy of S after thresholding
    """
    from numpy import log10, where
    smax = S.max()
    thresh = 10 ** (log10(smax) - dB / 10.)
    return where(S >= thresh, S, thresh)
