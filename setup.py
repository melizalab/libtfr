#!/usr/bin/env python
# -*- mode: python -*-
import os
import platform

from Cython.Build import cythonize
from sysconfig import get_platform
from setuptools import Extension, setup


def get_include_dirs():
    include_dirs = [os.path.join(os.getcwd(), "include")]
    build_platform = get_platform()
    if build_platform.startswith("linux"):
        include_dirs.append("/usr/include")
    elif build_platform in ("win32", "win-amd64"):
        fftw_path = os.environ.get("FFTW_PATH")
        if fftw_path:
            include_dirs.append(os.path.join(fftw_path, "include"))
        else:
            include_dirs.append(os.path.join(os.getcwd(), "fftw"))
        lapack_path = os.environ.get("LAPACK_PATH")
        if lapack_path:
            include_dirs.append(os.path.join(lapack_path, "include"))
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
    build_platform = get_platform()
    lib_dirs = []
    if build_platform.startswith("macosx"):
        # macports
        lib_dirs.append("/opt/local/lib")
        # homebrew
        lib_dirs.append("/opt/homebrew/lib")
        lib_dirs.append("/usr/local/lib")
    elif build_platform in ("win32", "win-amd64"):
        fftw_path = os.environ.get("FFTW_PATH")
        if fftw_path:
            lib_dirs.append(os.path.join(fftw_path, "lib"))
        lapack_path = os.environ.get("LAPACK_PATH")
        if lapack_path:
            lib_dirs.append(os.path.join(lapack_path, "lib"))
    return lib_dirs


def get_link_libraries():
    build_platform = get_platform()
    if build_platform in ("win32", "win-amd64"):
        link_libraries = ["fftw3", "openblas"]
    else:
        link_libraries = ["fftw3", "lapack", "m"]
    return link_libraries


def get_define_macros():
    # Defined on every platform, not just Windows: MSVC's C compiler has no C99
    # _Complex, so fftw_complex must be a double[2]. Keeping that unconditional
    # means there is one code path for the complex math rather than a Windows
    # path that only gets exercised in CI.
    return [("FFTW_NO_Complex", "1")]

compiler_settings = {
    "include_dirs": get_include_dirs(),
    "library_dirs": get_lib_dirs(),
    "libraries": get_link_libraries(),
    "define_macros": get_define_macros(),
    "extra_link_args": [],
}

cython_directives = {}

if platform.python_implementation() == 'PyPy':
    cython_directives.update({
        "legacy_implicit_noexcept": True,
        "c_api_binop_methods": False,
    })

sources = ["src/tfr.c", "src/mtm.c", "src/libtfr.pyx"]
extensions = [Extension("libtfr", sources=sources, **compiler_settings)]
ext_modules = cythonize(extensions, compiler_directives=cython_directives)
setup(name="libtfr", ext_modules=ext_modules)


# Variables:
# End:
