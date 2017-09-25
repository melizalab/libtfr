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
    # h = transform.tapers
    # X = interpolate(h.T, data, start, dt)
    X = transform.tapers_interpolate(data, start, dt)
    # exp(2 pi i omega t)
    f, idx = libtfr.fgrid(Fs, transform.nfft, (0, Fs/2))
    w = -2j * np.pi * f
    Y = np.exp(np.outer(w, data - start))
    # integrate over t as a matrix multiplication
    return np.dot(Y, X.T) - H * Msp


def mtstft_pt(transform, t, window, step, t0, tN):
    """
    Computes complex multitaper STFT of a point process

    times - input data (1D time series)
    window, step - determine window size and t0, dt - define the support of the window
    returns array of complex numbers, dimension (nreal, ntapers)
    """
    times = np.asarray(t).astype('d')
    dt = window / transform.npoints
    Fs = 1 / dt
    tgrid = np.arange(t0, tN - window, step)
    nframes = tgrid.size
    # exp(2 pi i omega t)
    f = libtfr.fgrid(Fs, transform.nfft)[0]
    w = -2j * np.pi * f
    H = transform.tapers_fft(np.sqrt(Fs)).T
    J = np.zeros((f.size, nframes, transform.ntapers), dtype=np.complex128)
    Nsp = np.zeros(nframes, dtype='i')
    for i in range(nframes):
        tw0 = t0 + i * step
        # this is not very efficient
        idx = (times >= tw0) & (times < (tw0 + window))
        events = times[idx]
        Nsp[i] = events.size
        Msp = 1. * Nsp[i] / transform.npoints
        # apply tapers to input times
        ht = transform.tapers_interpolate(events, tw0, dt)
        Y = np.exp(np.outer(w, events - tw0))
        J[:, i, :] = np.dot(Y, ht.T) - H * Msp
    return J, Nsp
