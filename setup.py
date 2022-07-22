#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -*- mode: python -*-
from setuptools import setup, Extension
import os
import sys
if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >= 3.6 required.")

from Cython.Distutils import build_ext

# --- Distutils setup and metadata --------------------------------------------
VERSION = '2.1.6'

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


class BuildExt(build_ext):
    def build_extensions(self):
        import numpy
        import pkgconfig
        compiler_settings = pkgconfig.parse("fftw3")
        compiler_settings['include_dirs'].insert(0, "include")
        compiler_settings['include_dirs'].append(numpy.get_include())
        if os.environ.get('STATIC_LAPACK', False):
            # For the manylinux wheels: we have to build and statically
            # link to lapack and blas.
            compiler_settings['extra_link_args'].extend(("/usr/src/lapack/liblapack.a",
                                                         "/usr/src/lapack/librefblas.a",
                                                         "-lgfortran"))
        else:
            compiler_settings['libraries'].extend(("lapack", "m"))
        c_opts = []
        if has_flag(self.compiler, '-ffast-math'):
            c_opts.append('-ffast-math')
        for ext in self.extensions:
            for k, v in compiler_settings.items():
                setattr(ext, k, v)
            ext.extra_compile_args.extend(c_opts)
        build_ext.build_extensions(self)


sources = ['src/tfr.c', 'src/mtm.c', 'src/libtfr.pyx']

setup(
    name='libtfr',
    version=VERSION,
    ext_modules=[Extension('libtfr', sources=sources)],
    cmdclass={'build_ext': BuildExt},
    description=short_desc,
    long_description=long_desc,
    classifiers=[x for x in cls_txt.split("\n") if x],
    author='C Daniel Meliza',
    maintainer='C Daniel Meliza',
    url='https://melizalab.github.io/libtfr/',
    install_requires=["numpy>=1.19.5"],
    zip_safe=False,
    test_suite='tests'
)


# Variables:
# End:
