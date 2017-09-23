# -*- mode: python -*-
# -*- coding: utf-8 -*-
""" Pure python implementation of point process fft """
from __future__ import division
from __future__ import unicode_literals

import numpy as np
from scipy.interpolate import interp1d
import libtfr

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
    t = np.arange(start, stop, dt)
    interpolator = interp1d(t, h.T, axis=0)
    X = interpolator(data)
    # exp(2 pi i omega t)
    f, idx = libtfr.fgrid(Fs, transform.nfft, (0, Fs/2))
    w = -2j * np.pi * f
    Y = np.exp(np.outer(w, data - start))
    # integrate over t as a matrix multiplication
    return np.dot(Y, X) - H * Msp
