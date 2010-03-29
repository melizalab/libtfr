#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- mode: python -*-
"""
pyext.py - tool chain for building python extension modules

Copyright (C) 2010 Daniel Meliza <dmeliza@dylan.uchicago.edu>
Created 2010-03-29
"""
import os
import SCons
from SCons.Builder import Builder
from SCons.Action import Action

import distutils.sysconfig


def generate(env,**kw):
    pybase = distutils.sysconfig.get_python_lib(plat_specific=1)
    pyvars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS', 'CCSHARED', 'LDSHARED', 'SO')
    (cc, cxx, opt, basecflags, ccshared, ldshared, so_ext) = [x if x!=None else '' for x in pyvars]
    nxdir = os.path.join(pybase, 'numpy/core/include')

    def build_module(env, target, csources, sharedobjects):
        pycppath = [env.get('CPPPATH'), distutils.sysconfig.get_python_inc(), nxdir],
        pycppflags = [env.get('CPPFLAGS'), env.Split(basecflags), env.Split(opt)]
        obj = env.SharedObject(csources, CC=cc, CPPPATH=pycppath, CPPFLAGS=pycppflags)
        return env.SharedLibrary(target, [obj + sharedobjects], SHLINK=ldshared, SHLINKFLAGS=[],
                                 SHLIBPREFIX="", SHLIBSUFFIX=so_ext)

    def install_module(env, sources):
        return env.Install(pybase, sources)

    env.AddMethod(build_module, "PyExt")
    env.AddMethod(install_module, "PyInstall")


def exists(env):
    return env['BUILDERS'].has_key('PyExt')
    
    

# Variables:
# End:
