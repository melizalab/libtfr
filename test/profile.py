# -*- mode: python -*-

#import pstats, cProfile

import libtfr
from numpy.random import randn

sig = randn(17590)
D = libtfr.dpss(256, 3, 5)

#cProfile.runctx("D.mtspec(sig, 10)", globals(), locals(), "Profile.prof")
#s = pstats.Stats("Profile.prof")
#s.strip_dirs().sort_stats("time").print_stats()
