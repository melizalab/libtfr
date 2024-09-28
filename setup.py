#!/usr/bin/env python
# -*- mode: python -*-
import os

from Cython.Build import cythonize
from pkg_resources import get_build_platform
from setuptools import Extension, setup


def get_include_dirs():
    include_dirs = [os.path.join(os.getcwd(), "include")]
    build_platform = get_build_platform()
    if build_platform.startswith("linux"):
        include_dirs.append("/usr/include")
    elif build_platform in ("win32", "win-amd64"):
        include_dirs.append(os.path.join(os.getcwd(), "fftw"))
    elif build_platform.startswith("freebsd"):
        include_dirs.append("/usr/local/include")
    elif build_platform.startswith("macosx"):
        # macports
        include_dirs.append("/opt/local/include")
        # homebrew
        include_dirs.append("/opt/homebrew/include")
        include_dirs.append("/usr/local/include")
    return include_dirs


def get_lib_dirs():
    build_platform = get_build_platform()
    lib_dirs = []
    if build_platform.startswith("macosx"):
        # macports
        lib_dirs.append("/opt/local/lib")
        # homebrew
        lib_dirs.append("/opt/homebrew/lib")
        lib_dirs.append("/usr/local/lib")
    return lib_dirs


compiler_settings = {
    "include_dirs": get_include_dirs(),
    "libraries": ["fftw3", "lapack", "m"],
    "library_dirs": get_lib_dirs(),
    "extra_link_args": [],
}

sources = ["src/tfr.c", "src/mtm.c", "src/libtfr.pyx"]
extensions = [Extension("libtfr", sources=sources, **compiler_settings)]

setup(
    name="libtfr",
    ext_modules=cythonize(extensions),
    # ext_modules=[Extension("libtfr", sources=sources)],
    # cmdclass={"build_ext": BuildExt},
)


# Variables:
# End:
