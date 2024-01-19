#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -*- mode: python -*-
from setuptools import setup, Extension
import os
import platform
from pkg_resources import get_build_platform

from Cython.Distutils import build_ext
from Cython.Build import cythonize


def get_include_dirs():
    include_dirs = [os.path.join(os.getcwd(), "include")]
    build_platform = get_build_platform()
    if build_platform.startswith("linux"):
        include_dirs.append("/usr/include")
    elif build_platform == in ("win32", "win-amd64"):
        include_dirs.append(os.path.join(os.getcwd(), "fftw"))
    elif build_platform.startswith("freebsd"):
        include_dirs.append("/usr/local/include")
    elif build_platform.startswith("macosx"):
        # macports
        include_dirs.append("/opt/local/include")
        arch = platform.machine()
        if arch == "arm64":
            include_dirs.append("/opt/homebrew/include")
        elif arch == "x86_64":
            include_dirs.append("/usr/local/include")
    return include_dirs


compiler_settings = {
    "include_dirs": get_include_dirs(),
    "libraries": ["fftw3", "lapack", "m"],
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
