#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -*- mode: python -*-
import sys
if sys.hexversion < 0x02070000:
    raise RuntimeError("Python 2.7 or higher required")

# setuptools 0.7+ doesn't play nice with distribute, so try to use existing
# package if possible
try:
    from setuptools import setup, Extension
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, Extension

try:
    from Cython.Distutils import build_ext
    SUFFIX = '.pyx'
except ImportError:
    from distutils.command.build_ext import build_ext
    SUFFIX = '.c'

import numpy

# --- Distutils setup and metadata --------------------------------------------
VERSION = '2.0.1'

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

import pkgconfig
compiler_settings = pkgconfig.parse("fftw3")
compiler_settings['include_dirs'].append(numpy.get_include())
compiler_settings['libraries'].append('lapack')
if sys.platform == 'darwin':
    compiler_settings['include_dirs'].append('/opt/local/include')
compiler_settings = dict((k,list(v)) for k,v in compiler_settings.items())

sources = ['tfr.c','mtm.c', 'libtfr' + SUFFIX,]

setup(
    name= 'libtfr',
    version= VERSION,
    ext_modules= [Extension('libtfr',
                             sources=sources,
                            **compiler_settings)],
    cmdclass={'build_ext': build_ext},
    description= short_desc,
    long_description= long_desc,
    author= 'C Daniel Meliza',
    author_email= 'dan@meliza.org',
    maintainer= 'C Daniel Meliza',
    maintainer_email= 'dan@meliza.org',
    url= 'http://melizalab.github.com/libtfr',
    download_url= 'https://github.com/downloads/melizalab/libtfr',
    setup_requires=["pkgconfig>=1.2"],
    zip_safe= False,
    test_suite='nose.collector'
)


# Variables:
# End:
