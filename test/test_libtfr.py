# -*- coding: utf-8 -*-
# -*- mode: python -*-
from __future__ import division
import sys

if sys.hexversion < 0x03000000:
    import unittest2 as unittest
else:
    import unittest

import numpy as np

import libtfr

def fmsin(N, fnormin, fnormax, period, t0, fnorm0, pm1):
    from numpy import sign, sin, arccos, exp, arange, pi, real
    fnormid = 0.5 * (fnormax+fnormin)
    delta = 0.5 * (fnormax-fnormin)
    phi = - sign(pm1) * arccos((fnorm0 - fnormid)/delta)
    t = arange(0.0, N)
    phase = 2 * pi * fnormid * (t - t0) + delta * period * (sin(2 * pi * (t - t0) / period + phi) - sin(phi))
    return real(exp(1j * phase))

def ppt(sig):
    from numpy import exp, random
    p = exp(sig - 1)
    events = (p > random.uniform(size=p.size)).nonzero()[0].astype('d')
    # jitter
    events += random.uniform(low=-0.25, high=0.25, size=events.size)
    return events

sig = fmsin(17664, 0.15, 0.45, 1024, 256./4, 0.3, -1)
events = ppt(sig)

# these values are from matlab's implementation of DPSS. Note that Prieto's library
# http://wwwprof.uniandes.edu.co/~gprieto/software/mwlib.html gives slightly
# different values. Could probably pick something more useful than the means
dpss_vals = [((128, 3, 5),
              (0.999999867187541, 0.999990855025110, 0.999717094982790, 0.994937064876468, 0.946249296585016),
              (0.066431706597466,   0.000000000000000, 0.044471301509568, 0.000000000000000, 0.033807993521897)),
             ((65536, 3.5, 6),
              (0.999999993658772, 0.999999484460332, 0.999980764501928, 0.999568509118410, 0.993676443756899, 0.941057523072299),
              (0.002829381466133, 0.000000000000020, 0.001915014524847, 0.000000000000030, 0.001539836330905, -0.000000000000026))]

# these values are from Xiao and Flandrin's matlab implementation
hermf_vals = [((201, 6, 6.0),
               (0.038241135791969, 0.000000000000000, 0.027040563146761, 0.000000000000000, 0.023417748319781, -0.000000000000000),
               (0.001800000000000, 0.005400000000000, 0.008999999999995, 0.012599999999893, 0.016199999998273, 0.019799999978400))]


class TestTapers(unittest.TestCase):

    def test_dpss(self):
        from numpy import ones_like
        for A, C, M in dpss_vals:
            with self.subTest(args=A, concentrations=C, means=M):
                E, V = libtfr.dpss(*A)
                self.assertTrue(np.allclose(V, C))
                self.assertTrue(np.allclose(E.mean(1), M))
                self.assertTrue(np.allclose((E**2).sum(1), ones_like(V)))

    def test_dpss_bad_args(self):
        with self.assertRaises(ValueError):
            E,V = libtfr.dpss(128, -5, 3)

    def test_hermf(self):
        for A, H, D in hermf_vals:
            with self.subTest(args=A, hmeans=H, dnorms=D):
                h, Dh, Th = libtfr.hermf(*A)
                self.assertTrue(np.allclose(h.mean(1), H))
                self.assertTrue(np.allclose((Dh**2).sum(1), D))


class TestTransforms(unittest.TestCase):
    def test_dpss_fft(self):
        from numpy import sqrt, fft
        for args, C, M in dpss_vals:
            with self.subTest(args=args, concentrations=C, means=M):
                D = libtfr.mfft_dpss(args[0], args[1], args[2], args[0])
                E = D.tapers
                Z = D.tapers_fft(1.0)
                self.assertTupleEqual(Z.shape, (args[2], args[0]//2 + 1))
                self.assertTrue(np.allclose(Z, fft.fft(E, axis=1)[:, :Z.shape[1]]))

    # these tests simply assert that the returned arrays have the correct shape and type
    def test_tfr(self):
        nfft = 256
        Np = 201
        shift = 10
        K = 6
        tm = 6.0
        flock = 0.01
        tlock = 5
        Z = libtfr.tfr_spec(sig, nfft, shift, Np, K, tm, flock, tlock)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, (sig.size - Np)// shift + 1))
        self.assertEqual(Z.dtype, libtfr.DTYPE)


    def test_dpss_mtfft(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtfft(sig)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, ntapers))
        self.assertEqual(Z.dtype, libtfr.CTYPE)


    def test_dpss_mtfft_pt_noevents(self):
        from numpy import zeros_like
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        J = D.mtfft_pt([], 1, 0)
        self.assertTupleEqual(J.shape, (nfft//2 + 1, ntapers))
        self.assertEqual(J.dtype, libtfr.CTYPE)
        self.assertTrue(np.allclose(J, zeros_like(J)))


    def test_dpss_mtfft_pt(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        J = D.mtfft_pt(events, 1, 0)
        self.assertTupleEqual(J.shape, (nfft//2 + 1, ntapers))
        self.assertEqual(J.dtype, libtfr.CTYPE)


    def test_dpss_mtpsd(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtpsd(sig)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1,))
        self.assertEqual(Z.dtype, libtfr.DTYPE)


    def test_dpss_mtspec(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtspec(sig, shift)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, nframes))
        self.assertEqual(Z.dtype, libtfr.DTYPE)


    def test_dpss_mtstft(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtstft(sig, shift)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, nframes, ntapers))
        self.assertEqual(Z.dtype, libtfr.CTYPE)


    def test_dpss_mtstft_pt_noevents(self):
        from numpy import zeros_like
        events = []
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z, Nsp = D.mtstft_pt(events, 1, shift, 0, sig.size)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, nframes, ntapers))
        self.assertEqual(Nsp.size, nframes)
        self.assertEqual(Z.dtype, libtfr.CTYPE)
        self.assertTrue(np.allclose(Z, zeros_like(Z)))


    def test_dpss_mtstft_pt(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z, Nsp = D.mtstft_pt(events, 1, shift, 0, sig.size)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, nframes, ntapers))
        self.assertEqual(Nsp.size, nframes)
        self.assertEqual(Z.dtype, libtfr.CTYPE)


    def test_hanning_mtstft(self):
        from numpy import hanning
        nfft = 256
        shift = 10
        window = hanning(nfft - 50)
        nframes = (sig.size - window.size) // shift + 1
        D = libtfr.mfft_precalc(nfft, window)
        Z = D.mtstft(sig, shift)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1, nframes, 1))
        self.assertEqual(Z.dtype, libtfr.CTYPE)


    def test_precalc_psd(self):
        nfft = 256
        E,V = libtfr.dpss(200, 3, 5)
        D = libtfr.mfft_precalc(nfft, E, V)
        self.assertTrue(np.allclose(E, D.tapers))
        Z = D.mtpsd(sig)
        self.assertTupleEqual(Z.shape, (nfft//2 + 1,))
        self.assertEqual(Z.dtype, libtfr.DTYPE)


class TestUtility(unittest.TestCase):

    def test_fgrid(self):
        Fs = 100
        nfft = 256
        f, idx = libtfr.fgrid(Fs, nfft)
        self.assertEqual(f.size, idx.size)
        self.assertEqual(f[-1], Fs / 2)
        f, idx = libtfr.fgrid(Fs, nfft, (10, 40))
        self.assertTrue(f[0] >= 10)
        self.assertTrue(f[-1] <= 40)


    def test_tgrid(self):
        nfft = 256
        shift = 10
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtstft(sig, shift)
        tgrid1 = libtfr.tgrid(sig.size, 1, shift)
        tgrid2 = libtfr.tgrid(Z, 1, shift)
        #assert_array_equal(tgrid1, tgrid2)


    def test_interpolation(self):
        from numpy import interp, arange
        nfft1 = 256
        nfft2 = nfft1 * 2
        ntapers = 5
        D1 = libtfr.mfft_dpss(nfft1, 3, ntapers, nfft1)
        t1 = arange(0, nfft1, 1)
        t2 = arange(0, nfft1, nfft1 / nfft2)
        h1_interp = D1.tapers_interpolate(t2, 0, 1)
        self.assertTupleEqual(h1_interp.shape, (ntapers, nfft2))
        for i in range(ntapers):
            self.assertTrue(np.allclose(h1_interp[i], interp(t2, t1, D1.tapers[i])))
