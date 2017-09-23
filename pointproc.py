# -*- mode: python -*-
# -*- coding: utf-8 -*-
""" Pure python implementation of point process fft """
from __future__ import division
from __future__ import unicode_literals

import numpy as np
import libtfr

def interpolate(y, t, t0, dt):
    """Linearly interpolate y(t) given y defined on evenly spaced grid"""
    nt = t.size
    ny, nd = y.shape
    out = np.zeros((nt, nd), dtype='d')
    for i, v in enumerate(t):
        ti = (v - t0) / dt
        i0 = int(ti)
        # assume y = 0 outside support
        if i0 + 1 == ny:
            out[i] = y[i0]
        elif i0 >= 0 and i0 < ny:
            tid = (ti - np.floor(ti))
            out[i] = y[i0] * (1 - tid) + y[i0 + 1] * tid
    return out


def mtfft(transform, events, start, stop):
    data = np.asarray(events, dtype='d')
    nevents = data.size
    dt = (stop - start) / transform.npoints
    Fs = 1. / dt
    # a finite size correction
    H = transform.tapers_fft(np.sqrt(Fs)).T
    Msp = nevents / transform.npoints
    # this is effectively events * tapers
    h = transform.tapers
    X = interpolate(h.T, data, start, dt)
    # exp(2 pi i omega t)
    f, idx = libtfr.fgrid(Fs, transform.nfft, (0, Fs/2))
    w = -2j * np.pi * f
    Y = np.exp(np.outer(w, data - start))
    # integrate over t as a matrix multiplication
    return np.dot(Y, X) - H * Msp
