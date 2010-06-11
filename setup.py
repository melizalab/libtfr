#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- mode: python -*-
from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension
import os,sys

ext_libs = ['fftw3']
ext_incl = []

if hasattr(os,'uname'):
    system = os.uname()[0]
else:
    system = 'Windows'

if system=='Darwin':
    ext_libs.append('lapack')
    ext_incl.append('/opt/local/include')
elif system=='Linux':
    ext_libs.extend(('atlas','cblas','f77blas','lapack'))
    
setup(
    name = 'libtfr',
    version = "1.0.0",  
    py_modules = ['libtfr'],
    ext_modules = [Extension('_libtfr',
                             sources=['libtfr.c','tfr.c','mtm.c'],
                             libraries=ext_libs,
                             include_dirs=ext_incl)],

    install_requires = ["numpy>=1.3"],

    description = "Python package for calculating tfr and mtm spectrograms",
    
    author = "CD Meliza",
    maintainer = "CD Meliza",
    maintainer_email = "dmeliza@uchicago.edu",
)


# Variables:
# End:
