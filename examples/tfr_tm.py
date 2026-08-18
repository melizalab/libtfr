# -*- mode: python -*-
"""
A worked example of the libtfr python interface.

Builds a signal whose instantaneous frequency is known, runs it through each of
the transforms, and checks that the reassigned spectrogram recovers the true
frequency law. Run it directly; pass a repeat count to use it as a crude timing
loop.

Copyright (C) 2009-2026 C Daniel Meliza
Created 2009-06-09

SPDX-License-Identifier: GPL-2.0-or-later
"""

import libtfr
import numpy as nx


def fmsin(N, fnormin=0.05, fnormax=0.45, period=None, t0=None, fnorm0=0.25, pm1=1):
    """
    Signal with sinusoidal frequency modulation.

    The instantaneous frequency at time t0 is fnorm0; pm1 (-1 or +1) resolves
    whether the frequency is rising or falling there.

    Returns the signal and its instantaneous frequency law, both length N-1.

    Original MATLAB code F. Auger, July 1995, licensed under the GPL.
    """
    if period is None:
        period = N
    if t0 is None:
        t0 = N / 2
    pm1 = nx.sign(pm1)

    fnormid = 0.5 * (fnormax + fnormin)
    delta = 0.5 * (fnormax - fnormin)
    phi = -pm1 * nx.arccos((fnorm0 - fnormid) / delta)
    time = nx.arange(1, N) - t0
    phase = 2 * nx.pi * fnormid * time + delta * period * (
        nx.sin(2 * nx.pi * time / period + phi) - nx.sin(phi)
    )
    y = nx.exp(1j * phase)
    iflaw = fnormid + delta * nx.cos(2 * nx.pi * time / period + phi)

    return y, iflaw


if __name__ == "__main__":
    import sys

    N = 256  # transform size
    NW = 3.5  # time-bandwidth product
    step = 10  # samples between frames
    k = 6  # number of tapers
    tm = 6.0  # time support of the hermite tapers
    Np = 201  # reassignment window size (odd)

    nloop = 1 if len(sys.argv) == 1 else int(sys.argv[1])

    # a sinusoidal FM sweep, plus noise, with the true frequency law in hand
    signal, iflaw = fmsin(17590, 0.15, 0.45, 1024, 256 / 4, 0.3, -1)
    s = signal.real + nx.random.randn(signal.size) / 2

    for i in range(nloop):
        if nloop > 1:
            print(f"loop {i}")

        # tapers can be generated on their own. hermf needs an odd window,
        # so it takes Np, the reassignment window, rather than the fft size
        h, Dh, Th = libtfr.hermf(Np, k, tm)
        E, V = libtfr.dpss(N, NW, k)

        # most transforms hang off a transform object, which owns the tapers
        D = libtfr.mfft_dpss(N, NW, k)
        psd = D.mtpsd(s[8300:8600])
        J = D.mtfft(s[8300:8600])
        mspec = D.mtspec(s, step)

        # a single arbitrary window, e.g. a plain hamming STFT
        H = libtfr.mfft_precalc(N, nx.hamming(N))
        spec = H.mtstft(s, step)

        # reassigned spectrograms, on linear, zoomed and logarithmic grids
        tspec = libtfr.tfr_spec(s, N, step, Np, k, tm)
        tspec_zoom = libtfr.tfr_spec(
            s, N, step, Np, k, tm, fgrid=nx.linspace(0.1, 0.475, 512)
        )
        tspec_log = libtfr.tfr_spec(
            s, N, step, Np, k, tm, fgrid=libtfr.log_fgrid(0.1, 0.45, 256)
        )

    print(f"hermite tapers   {h.shape}")
    print(f"dpss tapers      {E.shape}, concentrations {V.min():.4f}-{V.max():.4f}")
    print(f"multitaper psd   {psd.shape}")
    print(f"multitaper fft   {J.shape}")
    print(f"multitaper spec  {mspec.shape}")
    print(f"hamming stft     {spec.shape}")
    print(f"reassigned spec  {tspec.shape}")
    print(f"  zoomed         {tspec_zoom.shape}")
    print(f"  log frequency  {tspec_log.shape}")

    # the reassigned spectrogram should track the frequency law it was built
    # from: compare the peak frequency of each frame against the true value at
    # the centre of that frame
    frames = nx.arange(tspec.shape[1]) * step + Np // 2
    frames = frames[frames < iflaw.size]
    peak = tspec[:, : frames.size].argmax(0) / N
    error = nx.abs(peak - iflaw[frames])
    print(f"\nreassigned peak vs true frequency: median error {nx.median(error):.5f}")
