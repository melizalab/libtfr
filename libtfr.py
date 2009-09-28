#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""
Interface to libtfrspec library using numpy. 

Copyright C.D. Meliza, 2009
dmeliza@uchicago.edu
"""

import _libtfr

#__version__ = _libtfr.__version__

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

    returns an N/2+1 by L power spectrogram (L = length(s) / step), or if
    fgrid is specified, fgrid.size by L
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
    return S.swapaxes(1,2) if complex else S

def hermf(N, order=6, tm=6.0):
    """
    Computes a set of orthogonal Hermite functions for use in computing
    multi-taper reassigned spectrograms

    N - the number of points in the window (must be odd)
    order - the maximum order of the set of functions (default 6)
    tm - half-time support (default 6)

    Returns: h, Dh, Th
    h - hermite functions (MxN)
    Dh - first derivative of h (MxN)
    Th - time-multiple of h (MxN)

    """
    return _libtfr.hermf(N, order, tm)

if hasattr(_libtfr,'mtm_spec'):
    def mtm_spec(s, N, step, NW, k=0, adapt=True):
	"""
	Compute multitaper spectrogram using DPSS tapers

	s - input signal
	N - number of frequency points (i.e. window size)
	step - step size
	NW - time-frequency product
	k - number of tapers (default NW*2-1)
	adapt - compute adaptive spectrogram (default True)

	returns an N/2+1 by L power spectrogram
	"""
	return _libtfr.mtm_spec(s, N, step, NW, k, adapt)

if hasattr(_libtfr,'mtm_psd'):
    def mtm_psd(s, NW, k=0, adapt=True):
	"""
	Compute PSD of a signal using multitaper methods

	s - input signal
	NW - time-frequency product
	k - number of tapers (default NW*2-1)
	adapt - compute adaptive spectrum (default True)

	returns an N/2+1 real power spectrum density
	"""
	return _libtfr.mtm_psd(s, NW, k, adapt)


if hasattr(_libtfr,'mtm_psd'):
    def mtfft(s, NW, k=0, N=0):
	"""
	Compute multitaper transform of a signal

	s - input signal
	NW - time-frequency product
	k - number of tapers (default NW*2-1)
        N - number of points in FFT (default s.size)

	returns an NxK array of complex numbers
	"""
	return _libtfr.mtfft(s, NW, k, N)


if hasattr(_libtfr,'dpss'):
    def dpss(N, NW, k):
	"""
	Computes the discrete prolate spherical sequences used in the
	multitaper method power spectrum calculations.

	npoints   the number of points in the window
	mtm_p     the time-bandwidth product. Must be an integer or half-integer
		  (typical choices are 2, 5/2, 3, 7/2, or 4)
	k         If a scalar, returns the 1:k DPSS vectors
		  If a 2-ple, returns the k[0]:k[1] DPSS vectors
		  Default is to return all vectors

	returns:
	e - 2D array of eigenvectors, shape (n,npoints)
	v - 2D array of eigenvalues, length n = (mtm_p * 2 - 1)
	"""
	return _libtfr.dpss(N,NW,k)


def log_fgrid(fmin, fmax, N, Fs=None):
    """
    Generates a logarithmic frequency grid between fmin and fmax

    fmin - first frequency
    fmax  - last frequency
    N      - number of points
    base   - log base
    Fs     - set to a positive value to convert values
             to relative frequencies
    """
    from numpy import log, logspace, e
    lfmin, lfmax = log((fmin, fmax)) 
    out = logspace(lfmin, lfmax, N, base=e)
    if Fs:
        assert Fs > fmax, "Fs must be greater than fmax"
        return out / Fs
    else:
        return out
