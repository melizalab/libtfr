#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- mode: python -*-
from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension
import os, sys

ext_libs = ['fftw3']
ext_incl = []

if hasattr(os, 'uname'):
    system = os.uname()[0]
else:
    system = 'Windows'

if system == 'Darwin':
    ext_libs.append('lapack')
    ext_incl.append('/opt/local/include')
elif system == 'Linux':
    ext_libs.append('lapack')

cls_txt = """
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License (GPL)
Programming Language :: Python
Topic :: Scientific/Engineering
Operating System :: Unix
Operating System :: POSIX :: Linux
Operating System :: MacOS :: MacOS X
Natural Language :: English
"""

short_desc = "Calculates multi-taper windowed and time-frequency reassignment spectrograms"

long_desc = """
libtfr provides high-performance C and Python libraries for
calculating multitaper time-frequency reassignment (TFR) spectrograms
as well as conventional STFTs.  TFR spectrograms have enhanced
signal-to-noise characteristics and can provide very precise spectral
estimates under many conditions. The library requires FFTW for the underlying
FFT transformations.
"""

setup(
    name = 'libtfr',
    version = "1.0.3",
    py_modules = ['libtfr'],
    ext_modules = [Extension('_libtfr',
                             sources=['libtfr.c','tfr.c','mtm.c'],
                             libraries=ext_libs,
                             extra_compile_args=['-std=c99'],
                             include_dirs=ext_incl)],

    description = short_desc,
    long_description = long_desc,
    author = 'C Daniel Meliza',
    author_email = '"dan" at the domain "meliza.org"',
    maintainer = 'C Daniel Meliza',
    maintainer_email = '"dan" at the domain "meliza.org"',
    url = 'http://dmeliza.github.com/libtfr',
    download_url = 'https://github.com/downloads/dmeliza/libtfr',
)


# Variables:
# End:
