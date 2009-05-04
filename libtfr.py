#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""
Interface to libtfrspec library using numpy and ctypes
"""

import numpy as nx
import ctypes as cx
from numpy.ctypeslib import ndpointer

ltfr = cx.cdll.LoadLibrary('libtfrspec.dylib')

ltfr.mtm_init.restype = cx.c_void_p
ltfr.mtm_init.argtypes = [cx.c_int, cx.c_int, cx.c_int,
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]
ltfr.mtm_destroy.argtypes = [cx.c_void_p]

ltfr.mtfft.restype = cx.c_double
ltfr.mtfft.argtypes = [cx.c_void_p, ndpointer(dtype='h'), cx.c_int]
ltfr.mtfft_float.restype = cx.c_double
ltfr.mtfft_float.argtypes = [cx.c_void_p, ndpointer(dtype='f'), cx.c_int]
ltfr.mtfft_double.restype = cx.c_double
ltfr.mtfft_double.argtypes = [cx.c_void_p, ndpointer(dtype='d'), cx.c_int]

ltfr.mtpower.argtypes = [cx.c_void_p, ndpointer(dtype='d'), cx.c_double]

ltfr.mtm_spec.argtypes = [cx.c_void_p, ndpointer(dtype='d'), ndpointer(dtype='h'),
                           cx.c_int, cx.c_int, cx.c_int]
ltfr.mtm_spec_float.argtypes = [cx.c_void_p, ndpointer(dtype='d'), ndpointer(dtype='f'),
                           cx.c_int, cx.c_int, cx.c_int]
ltfr.mtm_spec_double.argtypes = [cx.c_void_p, ndpointer(dtype='d'), ndpointer(dtype='d'),
                           cx.c_int, cx.c_int, cx.c_int]

ltfr.hermf.argtypes = [cx.c_int, cx.c_int, cx.c_double,
                        ndpointer(dtype='d'), ndpointer(dtype='d'), ndpointer(dtype='d')]
ltfr.mtm_init_herm.restype = cx.c_void_p
ltfr.mtm_init_herm.argtypes = [cx.c_int, cx.c_int, cx.c_int, cx.c_double]

ltfr.tfr_displacements.argtypes = [cx.c_void_p, ndpointer(dtype='d'),
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]

ltfr.tfr_reassign.argtypes = [ndpointer(dtype='d'),ndpointer(dtype='d'),
                               ndpointer(dtype='d'),ndpointer(dtype='d'),
                               cx.c_int, cx.c_int, cx.c_double, cx.c_double, cx.c_double,
                               cx.c_int, cx.c_int]

ltfr.tfr_spec.argtypes = [cx.c_void_p, ndpointer(dtype='d'),ndpointer(dtype='h'),
                           cx.c_int, cx.c_int, cx.c_int, cx.c_double, cx.c_int]
ltfr.tfr_spec_float.argtypes = [cx.c_void_p, ndpointer(dtype='d'),ndpointer(dtype='f'),
                           cx.c_int, cx.c_int, cx.c_int, cx.c_double, cx.c_int]
ltfr.tfr_spec_double.argtypes = [cx.c_void_p, ndpointer(dtype='d'),ndpointer(dtype='d'),
                           cx.c_int, cx.c_int, cx.c_int, cx.c_double, cx.c_int]

# Note: these functions only get built if there's a functioning LAPACK
ltfr.dpss.restypes = None
ltfr.dpss.argtypes = [ndpointer(dtype='d'), ndpointer(dtype='d'),
                       cx.c_int, cx.c_double, cx.c_int]
ltfr.mtm_init_dpss.restype = cx.c_void_p
ltfr.mtm_init_dpss.argtypes = [cx.c_int, cx.c_double, cx.c_int]

def tfr_spec(s, N, step, Np, K=6, tm=6.0, flock=0.01, tlock=5):
    """
    Compute time-frequency reassignment spectrogram of input signal s

    s - input signal (short integer, float, and double supported)
    N - number of frequency points
    step - step size (in time points)
    Np - window size (should be <= N)
    K - number of tapers to use
    tm - time support of tapers
    flock - frequency locking parameter; power is not reassigned
            more than this value (in Hz)
    tlock - time locking parameter (in frames)

    returns an N/2+1 by L power spectrogram (L = length(s) / step)
    """
    if s.dtype=='h':
        fun = ltfr.tfr_spec
    elif s.dtype=='f':
        fun = ltfr.tfr_spec_float
    elif s.dtype=='d':
        fun = ltfr.tfr_spec_double
    else:
        raise TypeError, "libtfrspec doesn't support dtype %s" % s.dtype

    mtmh = ltfr.mtm_init_herm(N,Np,K,tm)
    spec = nx.zeros((N/2+1,s.size / step), order='F')
    fun(mtmh, spec, s, s.size, -1, step, flock, tlock)
    ltfr.mtm_destroy(mtmh)
    return spec

def mtm_spec(s, N, step, NW, k=None, adapt=True):
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
    if k==None:
        k = int(NW*2-1)
    if s.dtype=='h':
        fun = ltfr.mtm_spec
    elif s.dtype=='f':
        fun = ltfr.mtm_spec_float
    elif s.dtype=='d':
        fun = ltfr.mtm_spec_double
    else:
        raise TypeError, "libtfrspec doesn't support dtype %s" % s.dtype

    mtm = ltfr.mtm_init_dpss(N, NW, k)
    spec = nx.zeros((N/2+1,s.size / step), order='F')
    fun(mtm, spec, s, s.size, step, adapt)
    ltfr.mtm_destroy(mtm)
    return spec

def mtm_psd(s, NW, k=None, adapt=True):
    """
    Compute PSD of a signal using multitaper methods

    s - input signal
    NW - time-frequency product
    k - number of tapers (default NW*2-1)
    adapt - compute adaptive spectrogram (default True)

    returns an N/2+1 real power spectrum density
    """
    N = s.size
    if k==None:
        k = int(NW*2-1)    
    if s.dtype=='h':
        fun = ltfr.mtfft
    elif s.dtype=='f':
        fun = ltfr.mtfft_float
    elif s.dtype=='d':
        fun = ltfr.mtfft_double
    else:
        raise TypeError, "libtfrspec doesn't support dtype %s" % s.dtype

    mtm = ltfr.mtm_init_dpss(N, NW, k)
    sigpow = fun(mtm, s, N)
    out = nx.zeros(N/2+1)
    if adapt:
        ltfr.mtpower(mtm, out, sigpow)
    else:
        ltfr.mtpower(mtm, out, 0.0)

    ltfr.mtm_destroy(mtm)
    return out

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

    mpsd  = mtm_psd(s[8300:8600], NW)
    mspec = mtm_spec(s, N, step, NW)
    tspec = tfr_spec(s, N, step, Np, k, tm)
    
