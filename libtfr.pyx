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

import numpy as nx
cimport numpy as nx
nx.import_array()

ctypedef nx.double_t DTYPE_t
cimport tfr

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
    pass


def dpss(int N, float NW, int k):
    """
    Computes the discrete prolate spherical sequences used in the
    multitaper method power spectrum calculations.

    @param npoints   the number of points in the window
    @param mtm_p     the time-bandwidth product. Must be an integer or half-integer
                     (typical choices are 2, 5/2, 3, 7/2, or 4)
    @param k         Returns the 1:k DPSS vectors
                     If a 2-ple, returns the k[0]:k[1] DPSS vectors
                     Default is to return all vectors

    @returns (2D array of tapers, shape (k,npoints),
              2D array of concentration values, length k
    """
    cdef nx.ndarray[DTYPE_t, ndim=2] tapers = nx.empty((k, N), dtype='d')
    cdef nx.ndarray[DTYPE_t, ndim=1] lambdas = nx.empty(k, dtype='d')

    rv = tfr.dpss(&tapers[0,0], &lambdas[0], N, NW, k)
    if rv == 0:
        return tapers, lambdas
    elif rv == 1:
        raise ValueError("Invalid DPSS parameters")
    else:
        raise RuntimeError("Eigenvalue solver failed.")
