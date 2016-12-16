
Libtfr is a library for calculating multi-taper time-frequency reassignment (TFR) spectrograms. Time-frequency reassignment is a method that makes use of the instantaneous frequency and phase values in a spectrogram to 'deconvolve' the image, and can yield substantially sharper spectrograms with better signal-noise resolution than conventional windowed spectrograms (i.e. short-time Fourier transforms) or multi-taper spectrograms (using the discrete prolate spherical sequences).

The library will also calculate conventional windowed spectrograms and multitaper spectrograms.

Libtfr has C and Python APIs. The Python package is compatible with versions 2.7, 3.2, 3.3, and 3.4. It is known to work on Linux and OS X, but has not been tested on Windows.

[![Build Status](https://travis-ci.org/melizalab/libtfr.png?branch=master)](https://travis-ci.org/melizalab/libtfr)

## Quick start

To use the python module, you'll need to install some system dependencies first. On Debian:

```bash
sudo apt-get install libfftw3-dev liblapack-dev
```

On OS X with Macports:

```bash
sudo port install fftw-3
```

To install the python module from source:

```bash
pip install -r requirements.txt
python setup.py install
```

To install from PyPI:

```bash
pip install pkgconfig libtfr
```

Windows wheels with statically linked FFTW and LAPACK libraries have kindly been developed by [carlkl](https://github.com/carlkl). Install with `pip install -i https://pypi.anaconda.org/carlkl/simple libtfr`

To compute a time-frequency reassignment spectrogram in Python:

```python
import libtfr
nfft = 256
Np = 201
shift = 10
K = 6
tm = 6.0
flock = 0.01
tlock = 5

# load signal of dimension (npoints,)
signal = ...
S = libtfr.tfr_spec(signal, nfft, shift, Np, K, tm, flock, tlock)
```

See below for more information on the parameters.

### Mulitaper Spectral Analysis

Libtfr can also calculate multitaper and standard windowed Fourier transforms. For example, discrete prolate spherical sequences can be used to obtain multiple independent estimates of spectral statistics while maintaining optimal time-frequency tradeoffs (Prieto et al, 2007). The Python interface for MT calculations is somewhat different:

```python
import libtfr

# load signal of dimension (npoints,)
signal = ...

# generate a transform object with size equal to signal length and 5 tapers
D = libtfr.mfft_dpss(npoints, 3, 5)
# complex windowed FFTs, one per taper
Z = D.mtfft(signal)
# power spectrum density estimate using adaptive method to average tapers
P = D.mtpsd(signal)

# generate a transform object with size 512 samples and 5 tapers for short-time analysis
D = libtfr.mfft_dpss(512, 3, 5)
# complex STFT with frame shift of 10 samples
Z = D.mtstft(signal, 10)
# spectrogram with frame shift of 10 samples
P = D.mtspec(signal, 10)

# generate a transform object with 200-sample hanning window padded to 256 samples
from numpy import hanning
D = libtfr.mfft_precalc(256, hanning(200))
# spectrogram with frame shift of 10 samples
P = D.mtspec(signal, 10)
```

### C Library

To build the C library you will also need to have [scons](http://www.scons.org) installed. You may need to edit the SConstruct file to make sure it points to the correct LAPACK libraries. To build the shared library:

    scons lib

To install the libraries and header (default to `/usr/local/lib` and `/usr/local/include`):

    scons install

A small test program, *test_tfr*, can be built with `scons test`. The program generates a signal with sinusoidally modulated frequency and then calculates a multitaper PSD, a multitaper spectrogram, and a time-frequency reassigned spectrogram. The results are output in ASCII format to `tfr_in.dat`, `tfr_out_psd.dat`, `tfr_out_mtm.dat`, and `tfr_out_tfr.dat`.

See `test_tfr.c` for an example of how to use the C API.

### Documentation

The C header `tfr.h` and python module `libtfr.pyx` are both extensively documented.

### Algorithm and usage notes

The software was assembled from various MATLAB sources, including the time-frequency toolkit, Xiao and Flandrin's work on multitaper reassignment, and code from Gardner and Magnasco.

The basic principle is to use reassignment to increase the precision of the time-frequency localization, essentially by deconvolving the spectrogram with the TF representation of the window, recovering the center of mass of the spectrotemporal energy. Reassigned TFRs typically show a 'froth' for noise, and strong narrow lines for coherent signals like pure tones, chirps, and so forth. The use of multiple tapers reinforces the coherent signals while averaging out the froth, giving a very clean spectrogram with optimal precision and resolution properties.

Gardner & Magnasco calculate reassignment based on a different algorithm from Xiao and Flandrin. The latter involves 3 different FFT operations on the signal windowed with the hermitian taper *h(t)*, its derivative *h'(t)*, and its time product *t × h(t)*. The G&M algorithm only uses two FFTs, on the signal windowed with a Gaussian and its time derivative. If I understand their methods correctly, however, this derivation is based on properties of the fourier transform of the gaussian, and isn't appropriate for window functions based on the Hermitian tapers, which have more optimal distribution of energy over the TF plane (i.e., it takes fewer Hermitian tapers than Gaussian tapers to achieve the same quality spectrogram)

Therefore, the algorithm is mostly from X&F, though I include time and frequency locking parameters from G&M, which specify how far energy is allowed to be reassigned in the TF plane. Large displacements generally arise from numerical errors, so this helps to sharpen the lines somewhat. I also included the time/frequency interpolation from , which can be used to get higher precision (at the expense of less averaging) from smaller analysis windows.

Some fiddling with parameters is necessary to get the best spectrograms from a given sort of signal. Like the window size in an STFT, the taper parameters control the time-frequency resolution. However, in the reassignment spectrogram the precision (i.e. localization) is not affected by the taper size, so the effects of taper size will generally only be seen when two coherent signals are close to each other in time or frequency. `Nh` controls the size of the tapers; one can also adjust `tm`, the time support of the tapers, but depending on the number of tapers used, this shouldn't get a whole lot smaller than 5. Increased values of `Nh` result in improved narrowband resolution (i.e. between pure tones) but closely spaced clicks can become smeared. Decreasing `Nh` increases the resolution between broadband components (i.e. clicks) but smears closely spaced narrowband components. The effect of smearing can be ameliorated to some extent by adjusting the frequency/time locking parameters.

The frequency zoom parameter can be used to speed up calculation quite a bit. Since calculating the multitaper reassigned spectrogram takes 3xNtapers FFT operations, smaller FFTs are generally better. The spectrogram can then be binned at a finer resolution during reassignment. These two sets of parameters should generate fairly similar results:

    nfft=512, shift=10, tm=6, Nh=257, zoomf=1, zoomt=1 (default)
    nfft=256, shift=10, tm=6, Nh=257, zoomf=2, zoomt=1

Increasing the order generally reduces the background 'froth', but interference between closely spaced components may increase.

Additional improvements in resolution may be achievable averaging across different window sizes, or by using other averaging methods (i.e. as in Xiao and Flandrin)

### License

libtfr was written by C Daniel Meliza (dmeliza@uchicago.edu) and is licensed under the Gnu Public License (GPL) version 2; see COPYING for details.

some code is adapted from chronux (<http://www.chronux.org>), by Partha Mitra and Hemant Bokil, also licensed under GPL version 2

THE PROGRAMS ARE PROVIDED "AS IS" WITHOUT WARRANTY OF MERCANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. IN NO EVENT SHALL THE UNIVERSITY OF CHICAGO OR DR. MELIZA BE LIABLE FOR ANY DIRECT OR CONSEQUENTIAL DAMAGES RESULTING FROM USE OF THE PROGRAMS. THE USER BEARS THE ENTIRE RISK FOR USE OF THE PROGRAMS.

### References

* Time-frequency toolkit: <http://tftb.nongnu.org/>
* Xiao, J. & Flandrin, P. Multitaper Time-Frequency Reassignment for Nonstationary Spectrum Estimation and Chirp Enhancement Signal Processing, IEEE Transactions on, Signal Processing, IEEE Transactions on, 2007, 55, 2851-2860 code: <http://perso.ens-lyon.fr/patrick.flandrin/multitfr.html>
* Gardner, T. J. & Magnasco, M. O. Sparse time-frequency representations. Proc. Natl. Acad. Sci. U S A, 2006, 103, 6094-6099 code: <http://web.mit.edu/tgardner/www/Downloads/Entries/2007/10/22_Blue_bird_day_files/ifdv.m>
* Prieto, G.A., Parker, R. L., Thomson D. J., Vernon, F. L., & Graham, R. L. Reducing the bias of multitaper spectrum estimates. Geophys. J. Int. 2007, doi: 10.1111/j.1365-246X.2007.03592.x.
