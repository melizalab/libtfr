# -*- mode: python -*-

import libtfr
import numpy as np
import pytest


def fmsin(N, fnormin, fnormax, period, t0, fnorm0, pm1):
    from numpy import arange, arccos, exp, pi, real, sign, sin

    fnormid = 0.5 * (fnormax + fnormin)
    delta = 0.5 * (fnormax - fnormin)
    phi = -sign(pm1) * arccos((fnorm0 - fnormid) / delta)
    t = arange(0.0, N)
    phase = 2 * pi * fnormid * (t - t0) + delta * period * (
        sin(2 * pi * (t - t0) / period + phi) - sin(phi)
    )
    return real(exp(1j * phase))


def ppt(sig, seed=20240818):
    from numpy import exp
    from numpy.random import default_rng

    # seeded so a CI failure is reproducible; the point-process tests only
    # check shapes and dtypes, so the particular draw does not matter
    rng = default_rng(seed)
    p = exp(sig - 1)
    events = (p > rng.uniform(size=p.size)).nonzero()[0].astype("d")
    # jitter
    events += rng.uniform(low=-0.25, high=0.25, size=events.size)
    return events


sig = fmsin(17664, 0.15, 0.45, 1024, 256.0 / 4, 0.3, -1)
events = ppt(sig)

# these values are from matlab's implementation of DPSS. Note that Prieto's library
# http://wwwprof.uniandes.edu.co/~gprieto/software/mwlib.html gives slightly
# different values. Could probably pick something more useful than the means
dpss_vals = [
    (
        (128, 3, 5),
        (
            0.999999867187541,
            0.999990855025110,
            0.999717094982790,
            0.994937064876468,
            0.946249296585016,
        ),
        (
            0.066431706597466,
            0.000000000000000,
            0.044471301509568,
            0.000000000000000,
            0.033807993521897,
        ),
    ),
    (
        (65536, 3.5, 6),
        (
            0.999999993658772,
            0.999999484460332,
            0.999980764501928,
            0.999568509118410,
            0.993676443756899,
            0.941057523072299,
        ),
        (
            0.002829381466133,
            0.000000000000020,
            0.001915014524847,
            0.000000000000030,
            0.001539836330905,
            -0.000000000000026,
        ),
    ),
]

# these values are from Xiao and Flandrin's matlab implementation
hermf_vals = [
    (
        (201, 6, 6.0),
        (
            0.038241135791969,
            0.000000000000000,
            0.027040563146761,
            0.000000000000000,
            0.023417748319781,
            -0.000000000000000,
        ),
        (
            0.001800000000000,
            0.005400000000000,
            0.008999999999995,
            0.012599999999893,
            0.016199999998273,
            0.019799999978400,
        ),
    )
]

# (nfft, step, Np, K) combinations for the reassignment tests. Np must be odd
# and no larger than nfft.
tfr_geometries = [
    (256, 64, 201, 6),
    (512, 64, 255, 5),
    (128, 32, 101, 4),
]

# signals carrying no power at all, in the various ways that can happen
no_power_signals = [
    ("zeros", np.zeros(4096)),
    ("denormal", np.full(4096, 1e-300)),
    ("negative_zero", np.full(4096, -0.0)),
]


class TestTapers:
    @pytest.mark.parametrize("args, concentrations, means", dpss_vals)
    def test_dpss(self, args, concentrations, means):
        from numpy import ones_like

        E, V = libtfr.dpss(*args)
        assert np.allclose(V, concentrations)
        assert np.allclose(E.mean(1), means)
        assert np.allclose((E**2).sum(1), ones_like(V))

    def test_dpss_bad_args(self):
        with pytest.raises(ValueError):
            _ = libtfr.dpss(128, -5, 3)

    @pytest.mark.parametrize("args, hmeans, dnorms", hermf_vals)
    def test_hermf(self, args, hmeans, dnorms):
        h, Dh, _Th = libtfr.hermf(*args)
        assert np.allclose(h.mean(1), hmeans)
        assert np.allclose((Dh**2).sum(1), dnorms)


class TestTransforms:
    @pytest.mark.parametrize("args, concentrations, means", dpss_vals)
    def test_dpss_fft(self, args, concentrations, means):
        from numpy import fft

        D = libtfr.mfft_dpss(args[0], args[1], args[2], args[0])
        E = D.tapers
        Z = D.tapers_fft(1.0)
        assert Z.shape == (args[2], args[0] // 2 + 1)
        assert np.allclose(Z, fft.fft(E, axis=1)[:, : Z.shape[1]])

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
        assert Z.shape == (nfft // 2 + 1, (sig.size - Np) // shift + 1)
        assert Z.dtype == libtfr.DTYPE

    @pytest.mark.parametrize("nfft, step, Np, K", tfr_geometries)
    @pytest.mark.parametrize(
        "signal",
        [s for _, s in no_power_signals],
        ids=[name for name, _ in no_power_signals],
    )
    def test_tfr_degenerate_input(self, signal, nfft, step, Np, K):
        """Signals with no power give an all-zero spectrogram, never NaN/inf.

        tfr_displacements divides by z1, the transform of the first taper,
        which is exactly zero for a digitally silent frame and underflows to
        zero for a denormal one. tfr_reassign currently insulates the output
        from that: a non-finite displacement produces an out-of-range target
        bin, and q is below threshold as well, so the frame is dropped on two
        independent counts. This test pins the resulting API contract -- it
        does not constrain how the division itself is implemented, since the
        intermediate value never reaches the caller.
        """
        Z = libtfr.tfr_spec(signal, nfft, step, Np, K, 6.0, 0.01, 5)
        assert np.isfinite(Z).all(), "non-finite values in spectrogram"
        # no power in, no power out
        assert (Z == 0).all()

    @pytest.mark.parametrize("nfft, step, Np, K", tfr_geometries)
    def test_tfr_constant_input(self, nfft, step, Np, K):
        """A constant signal does carry power, and must survive intact."""
        Z = libtfr.tfr_spec(np.ones(4096), nfft, step, Np, K, 6.0, 0.01, 5)
        assert np.isfinite(Z).all()
        assert Z.max() > 0.0

    @pytest.mark.parametrize("nfft, step, Np, K", tfr_geometries)
    @pytest.mark.parametrize("f0", [500.0, 1000.0, 1500.0, 2000.0])
    def test_tfr_reassignment_concentration(self, f0, nfft, step, Np, K):
        """A pure tone must reassign onto its own frequency bin.

        This is the point of reassignment, and it is the only numerical check
        on tfr_displacements -- test_tfr above covers shape and dtype only.
        Asserting the structure (which bin, how concentrated) rather than
        stored reference values keeps it stable across platforms, where the
        floating-point results differ in the last ulp.
        """
        Fs = 8000.0
        t = np.arange(8192) / Fs
        Z = libtfr.tfr_spec(np.sin(2 * np.pi * f0 * t), nfft, step, Np, K, 6.0, 0.01, 5)
        power = Z.mean(1)
        freqs = np.arange(Z.shape[0]) / nfft * Fs
        peak = power.argmax()
        assert freqs[peak] == f0
        # essentially all the energy lands within a bin of the peak
        concentration = power[max(0, peak - 1) : peak + 2].sum() / power.sum()
        assert concentration > 0.99

    def test_dpss_mtfft(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtfft(sig)
        assert Z.shape == (nfft // 2 + 1, ntapers)
        assert Z.dtype == libtfr.CTYPE

    def test_dpss_mtfft_pt_noevents(self):
        from numpy import zeros_like

        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        J = D.mtfft_pt([], 1, 0)
        assert J.shape == (nfft // 2 + 1, ntapers)
        assert J.dtype == libtfr.CTYPE
        assert np.allclose(J, zeros_like(J))

    def test_dpss_mtfft_pt(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        J = D.mtfft_pt(events, 1, 0)
        assert J.shape == (nfft // 2 + 1, ntapers)
        assert J.dtype == libtfr.CTYPE

    def test_dpss_mtpsd(self):
        nfft = sig.size
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtpsd(sig)
        assert Z.shape == (nfft // 2 + 1,)
        assert Z.dtype == libtfr.DTYPE

    def test_dpss_mtspec(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtspec(sig, shift)
        assert Z.shape == (nfft // 2 + 1, nframes)
        assert Z.dtype == libtfr.DTYPE

    def test_dpss_mtstft(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtstft(sig, shift)
        assert Z.shape == (nfft // 2 + 1, nframes, ntapers)
        assert Z.dtype == libtfr.CTYPE

    def test_dpss_mtstft_pt_noevents(self):
        from numpy import zeros_like

        events = []
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z, Nsp = D.mtstft_pt(events, 1, shift, 0, sig.size)
        assert Z.shape == (nfft // 2 + 1, nframes, ntapers)
        assert Nsp.size == nframes
        assert Z.dtype == libtfr.CTYPE
        assert np.allclose(Z, zeros_like(Z))

    def test_dpss_mtstft_pt(self):
        nfft = 256
        shift = 10
        ntapers = 5
        nframes = (sig.size - nfft) // shift + 1
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z, Nsp = D.mtstft_pt(events, 1, shift, 0, sig.size)
        assert Z.shape == (nfft // 2 + 1, nframes, ntapers)
        assert Nsp.size == nframes
        assert Z.dtype == libtfr.CTYPE

    def test_hanning_mtstft(self):
        from numpy import hanning

        nfft = 256
        shift = 10
        window = hanning(nfft - 50)
        nframes = (sig.size - window.size) // shift + 1
        D = libtfr.mfft_precalc(nfft, window)
        Z = D.mtstft(sig, shift)
        assert Z.shape == (nfft // 2 + 1, nframes, 1)
        assert Z.dtype == libtfr.CTYPE

    def test_precalc_psd(self):
        nfft = 256
        E, V = libtfr.dpss(200, 3, 5)
        D = libtfr.mfft_precalc(nfft, E, V)
        assert np.allclose(E, D.tapers)
        Z = D.mtpsd(sig)
        assert Z.shape == (nfft // 2 + 1,)
        assert Z.dtype == libtfr.DTYPE


class TestUtility:
    def test_fgrid(self):
        Fs = 100
        nfft = 256
        f, idx = libtfr.fgrid(Fs, nfft)
        assert f.size == idx.size
        assert f[-1] == Fs / 2
        f, idx = libtfr.fgrid(Fs, nfft, (10, 40))
        assert f[0] >= 10
        assert f[-1] <= 40

    def test_tgrid(self):
        nfft = 256
        shift = 10
        ntapers = 5
        D = libtfr.mfft_dpss(nfft, 3, ntapers, nfft)
        Z = D.mtstft(sig, shift)
        _tgrid1 = libtfr.tgrid(sig.size, 1, shift)
        _tgrid2 = libtfr.tgrid(Z, 1, shift)
        # assert_array_equal(tgrid1, tgrid2)

    @pytest.mark.parametrize("i", range(5))
    def test_interpolation(self, i):
        from numpy import arange, interp

        nfft1 = 256
        nfft2 = nfft1 * 2
        ntapers = 5
        D1 = libtfr.mfft_dpss(nfft1, 3, ntapers, nfft1)
        t1 = arange(0, nfft1, 1)
        t2 = arange(0, nfft1, nfft1 / nfft2)
        h1_interp = D1.tapers_interpolate(t2, 0, 1)
        assert h1_interp.shape == (ntapers, nfft2)
        assert np.allclose(h1_interp[i], interp(t2, t1, D1.tapers[i]))
