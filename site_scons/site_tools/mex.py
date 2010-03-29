# mex.py:  Matlab extension builder
# Joe VanAndel, vanandel@ucar.edu, 2010/1/15


import os
import re
import SCons
from SCons.Builder import Builder
from SCons.Action import Action
from subprocess import Popen,PIPE


def generate(env):
##     _options = env.Variables()
##     _options.Add('MEXC',
##                  'Path to mex command, or else "mex" if unset.')
##     _options.Update(env)
    # find mex if MEXC not set
    if not env.get('MEXC'):
        extra_paths = [ '/usr/bin' ]
        if env.has_key('OPT_PREFIX'):
            extra_paths.append("%s/bin" % env['OPT_PREFIX'])
        opts = ['el4','el3','ws3','fc4','fc3','fc2']
        extra_paths.extend([ "/net/opt_lnx/local_%s/bin" % o for o in opts])
        env['MEXC'] = env.WhereIs('mex', extra_paths)
    if not env.get('MEXC'):
        SCons.Warnings.warn(SCons.Warnings.DependencyWarning, "Could not find mex program.")
        return
    # matlab *should* be in the same directory; use it to check extension
    try:
        cmd = [os.path.join(os.path.split(env['MEXC'])[0],'matlab'), '-nodisplay', '-nojvm' ]
        p1 = Popen(cmd,stdin=PIPE,stdout=PIPE)
        os.write(p1.stdin.fileno(), "mexext\n")
        os.write(p1.stdin.fileno(), "quit\n")
        # invoke grep to retrieve the extension
        p2 = Popen(['grep', 'mex'], stdin=p1.stdout, stdout=PIPE)
        env['MEX_EXT']  = p2.communicate()[0][:-1]

        bld = Builder(action = '$MEXC $SOURCES -o $TARGET',
                      suffix=env['MEX_EXT'])
        env['BUILDERS']['MEX'] = bld
    except OSError, e:
        # couldn't find matlab most likely
        SCons.Warnings.warn(SCons.Warnings.DependencyWarning, "Could not find matlab program.")

def exists(env):
    print "exists?"
    return env['BUILDERS'].has_key('MEX')
