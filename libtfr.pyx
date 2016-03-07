# -*- coding: utf-8 -*-
# -*- mode: cython -*-
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

Copyright C Daniel Meliza 2010-2016.  Licensed for use under GNU
General Public License, Version 2.  See COPYING for details.
"""
from cython cimport view
import numpy as nx
cimport numpy as nx
nx.import_array()
cimport tfr

DTYPE = nx.double
ctypedef nx.double_t DTYPE_t
CTYPE = nx.complex128
ctypedef nx.complex128_t CTYPE_t


# def class mfft:
#     cdef tfr.mfft * _mfft
#     def __cinit__(self):
#         self._mfft = mtm_init


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
    cdef nx.ndarray[DTYPE_t, ndim=1] samples = nx.asarray(s).astype(DTYPE)

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
    nt = tfr.mtm_nframes(mtmh, samples.size, step)
    cdef nx.ndarray[DTYPE_t, ndim=2] out = nx.zeros((nfreq, nt), dtype=DTYPE, order='F')

    tfr.tfr_spec(mtmh, &out[0,0], &samples[0], samples.size, -1,
                 step, flock, tlock, nfreq, fgridp);
    tfr.mtm_destroy(mtmh)

    return out


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


def mtfft(s not None, double NW, int k=0, N=None):
    """Compute multitaper transform of a signal

    @param s 1D input signal
    @param NW time-frequency product
    @param k number of tapers (default NW*2-1)
    @param N number of points in FFT (default s.size)

    @returns 2D complex array, dimension N x k
    """
    cdef tfr.mfft * mtm
    cdef int npoints

    # coerce data to proper type
    cdef nx.ndarray[DTYPE_t, ndim=1] samples = nx.asarray(s).astype(DTYPE)

    if k < 1:
        k = <int>(NW*2)-1;
    if N is None:
        npoints = samples.size
    else:
        npoints = N

    mtmh = tfr.mtm_init_dpss(npoints, NW, k)
    tfr.mtfft(mtmh, &samples[0], npoints);
    cdef nx.ndarray spec = hc2cmplx(mtmh)
    tfr.mtm_destroy(mtmh)
    return spec


cdef hc2cmplx(tfr.mfft * mtm):
    """Copy data from workspace of mtm object into a complex array"""
    cdef int nfft = tfr.mtm_nfft(mtm)
    cdef int ntapers = tfr.mtm_ntapers(mtm)
    cdef int real_count = nfft / 2 + 1
    cdef int imag_count = (nfft + 1) / 2
    cdef size_t t, n
    cdef double * buf = tfr.mtm_buffer(mtm)

    # allocate output array
    cdef nx.ndarray[CTYPE_t, ndim=2] out = nx.zeros((ntapers, real_count), dtype=CTYPE)
    for t in range(ntapers):
        for n in range(0, real_count):
            out[t, n].real = buf[t*nfft+n]
        for n in range(1, imag_count):
            out[t, n].imag = buf[t*nfft+(nfft-n)]
    return out
