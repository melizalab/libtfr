#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- mode: python -*-
"""
Just some quick tests of the libtfr python interface.

Copyright (C) 2009 Daniel Meliza <dmeliza@meliza-laptop-1.uchicago.edu>
Created 2009-06-09
"""

import numpy as nx
from scipy.signal import get_window
from libtfr import *

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
    PM1     : frequency direction at T0 (-1 or +1)       (default: +1  )

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

    w = get_window('hamming',N)

    for i in range(1):
	print i
	h,Dh,Th = hermf(N, k, tm)
	E,V   = dpss(N, NW, k)
	mpsd  = mtm_psd(s[8300:8600], NW)
	J     = mtfft(s[8300:8600], NW)
        spec  = stft(s, w, step)
	mspec = mtm_spec(s, N, step, NW)
	tspec = tfr_spec(s, N, step, Np, k, tm)


# Variables:
# End: