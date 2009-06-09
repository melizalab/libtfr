#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""
Interface to libtfrspec library using numpy

Copyright C.D. Meliza, 2009 (except for fmsin)
dmeliza@uchicago.edu
"""

import numpy as nx
import _libtfr

#__version__ = _libtfr.__version__

def tfr_spec(s, N, step, Np, K=6, tm=6.0, flock=0.01, tlock=5):
    """
    Compute time-frequency reassignment spectrogram of input signal s

    s - input signal (short integer, float, and double supported)
    N - number of frequency points
    step - step size (in time points)
    Np - window size (should be <= N)
    K - number of tapers to use (default 6)
    tm - time support of tapers (default 6.0)
    flock - frequency locking parameter; power is not reassigned
            more than this value (in Hz; default 0.01)
    tlock - time locking parameter (in frames; default 5)

    returns an N/2+1 by L power spectrogram (L = length(s) / step)
    """
    return _libtfr.tfr_spec(s, N, step, Np, K, tm, flock, tlock)

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
    def mtfft(s, NW, k=0):
        """
        Compute multitaper transform of a signal

        s - input signal
        NW - time-frequency product
        k - number of tapers (default NW*2-1)

        returns an NxK array of complex numbers
        """
        return _libtfr.mtfft(s, NW, k)

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
        e - 2D array of eigenvectors, shape (npoints, n)
        v - 2D array of eigenvalues, length n = (mtm_p * 2 - 1)
        """
        return _libtfr.dpss(N,NW,k)


def fmsin(N, fnormin=0.05, fnormax=0.45, period=None, t0=None, fnorm0=0.25, pm1=1):
    """
    Signal with sinusoidal frequency modulation.

    generates a frequency modulation with a sinusoidal frequency.
    This sinusoidal modulation is designed such that the instantaneous
    frequency at time T0 is equal to FNORM0, and the ambiguity between
    increasing or decreasing frequency is solved by PM1.
	
    N       : number of points.
    FNORMIN : smallest normalized frequency          (default: 0.05) 
    FNORMAX : highest normalized frequency           (default: 0.45)
    PERIOD  : period of the sinusoidal fm            (default: N   )
    T0      : time reference for the phase           (default: N/2 )
    FNORM0  : normalized frequency at time T0        (default: 0.25)
    PM1     : frequency direction at T0 (-1 or +1)	 (default: +1  )

    Returns:
    Y       : signal
    IFLAW   : its instantaneous frequency law
	
    Example: 
    z,i=fmsin(140,0.05,0.45,100,20,0.3,-1.0)
	
    Original MATLAB code F. Auger, July 1995.
    (note: Licensed under GPL; see main LICENSE file)
    """

    if period==None:
        period = N
    if t0==None:
        t0 = N/2
    pm1 = nx.sign(pm1)

    fnormid=0.5*(fnormax+fnormin);
    delta  =0.5*(fnormax-fnormin);
    phi    =-pm1*nx.arccos((fnorm0-fnormid)/delta);
    time   =nx.arange(1,N)-t0;
    phase  =2*nx.pi*fnormid*time+delta*period*(nx.sin(2*nx.pi*time/period+phi)-nx.sin(phi));
    y      =nx.exp(1j*phase)
    iflaw  =fnormid+delta*nx.cos(2*nx.pi*time/period+phi);

    return y,iflaw

if __name__=="__main__":
    N = 256
    NW = 3.5
    step = 10
    k = 6
    tm = 6.0
    Np = 201

    # generate a nice dynamic signal
    siglen = 17590
    ss,iff = fmsin(siglen, .15, 0.45, 1024, 256/4,0.3,-1)
    s = ss.real + nx.random.randn(ss.size)/2

    for i in range(100):
        print i
        E,V   = dpss(N, NW, k)
        mpsd  = mtm_psd(s[8300:8600], NW)
        J     = mtfft(s[8300:8600], NW)
        mspec = mtm_spec(s, N, step, NW)
        tspec = tfr_spec(s, N, step, Np, k, tm)
