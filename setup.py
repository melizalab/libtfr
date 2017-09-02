#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -*- mode: python -*-
import sys
if sys.hexversion < 0x02070000:
    raise RuntimeError("Python 2.7 or higher required")
from setuptools import setup, Extension

try:
    from Cython.Distutils import build_ext
    SUFFIX = '.pyx'
except ImportError:
    from distutils.command.build_ext import build_ext
    SUFFIX = '.c'


# --- Distutils setup and metadata --------------------------------------------
VERSION = '2.0.2'

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


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    from setuptools import distutils
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except distutils.errors.CompileError:
            return False
    return True


# these helper functions are used to postpone imports until dependencies are installed
class BuildExt(build_ext):
    def build_extensions(self):
        import numpy
        import pkgconfig
        compiler_settings = pkgconfig.parse("fftw3")
        compiler_settings['include_dirs'].append(numpy.get_include())
        c_opts = []
        if has_flag(self.compiler, '-ffast-math'):
            c_opts.append('-ffast-math')
        for ext in self.extensions:
            for k, v in compiler_settings.items():
                getattr(ext, k).extend(v)
            ext.extra_compile_args.extend(c_opts)
        build_ext.build_extensions(self)


# import pkgconfig
# compiler_settings = pkgconfig.parse("fftw3")
# compiler_settings['include_dirs'].append(numpy.get_include())
# compiler_settings['libraries'].append('lapack')
# if sys.platform == 'darwin':
#     compiler_settings['include_dirs'].append('/opt/local/include')
# compiler_settings = dict((k, list(v)) for k, v in compiler_settings.items())

sources = ['tfr.c', 'mtm.c', 'libtfr' + SUFFIX]

setup(
    name='libtfr',
    version=VERSION,
    ext_modules=[Extension('libtfr', sources=sources, libraries=['lapack'])],
    cmdclass={'build_ext': BuildExt},
    description=short_desc,
    long_description=long_desc,
    author='C Daniel Meliza',
    author_email='dan@meliza.org',
    maintainer='C Daniel Meliza',
    maintainer_email='dan@meliza.org',
    url='https://melizalab.github.io/libtfr/',
    download_url='https://github.com/melizalab/libtfr/archive/2.0.1.tar.gz',
    setup_requires=["pkgconfig>=1.2", "numpy>=1.10"],
    zip_safe=False,
    test_suite='nose.collector'
)


# Variables:
# End:
