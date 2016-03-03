#!/usr/bin/env python
 # -*- coding: utf-8 -*-
# -*- mode: python -*-
from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension
import os, sys

sources = ['libtfr.c','tfr.c','mtm.c']
ext_libs = []
ext_libdirs = []
ext_incl = []
compiler_defines = []
compiler_args = []
package_data = []

if hasattr(os, 'uname'):
    system = os.uname()[0]
else:
    system = 'Windows'

if system == 'Darwin':
    ext_libs.extend(('lapack','fftw3'))
    ext_incl.append('/opt/local/include')
elif system == 'Linux':
    ext_libs.extend(('lapack','fftw3'))
elif system == 'Windows':
    ext_incl.append('fftw3')
    ext_libs.append('fftw3-3')
    ext_libdirs.append('fftw3')
    compiler_defines.append(('NO_LAPACK', None))
    package_data.append(('', ['fftw3/libfftw3-3.dll']))

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
    name= 'libtfr',
    version= "1.0.6",
    py_modules= ['libtfr'],
    ext_modules= [Extension('_libtfr',
                             sources=sources,
                             define_macros=compiler_defines,
                             libraries=ext_libs,
                             extra_compile_args=compiler_args,
                             include_dirs=ext_incl,
                             library_dirs=ext_libdirs)],

    description= short_desc,
    long_description= long_desc,
    author= 'C Daniel Meliza',
    author_email= '"dan" at the domain "meliza.org"',
    maintainer= 'C Daniel Meliza',
    maintainer_email= '"dan" at the domain "meliza.org"',
    url= 'http://melizalab.github.com/libtfr',
    download_url= 'https://github.com/downloads/melizalab/libtfr',
    data_files= package_data,
    zip_safe= False
)


# Variables:
# End:
