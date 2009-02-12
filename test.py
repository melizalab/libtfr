
import numpy as nx
import ctypes as cx
from numpy.ctypeslib import ndpointer

lsono = cx.cdll.LoadLibrary('libsono2.dylib')

lsono.dpss.restypes = None
lsono.dpss.argtypes = [ndpointer(dtype='d'), ndpointer(dtype='d'),
                       cx.c_int, cx.c_double, cx.c_int]

lsono.init_mtm_prealloc.restype = cx.c_void_p
lsono.init_mtm_prealloc.argtypes = [cx.c_int, cx.c_int, cx.c_int,
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]

lsono.mtfft.restype = cx.c_double
lsono.mtfft.argtypes = [cx.c_void_p, ndpointer(dtype='h'), cx.c_int]

lsono.mtpower.argtypes = [cx.c_void_p, ndpointer(dtype='d'), cx.c_double]

N = 100
NW = 3.5
k = 5

# allocate storage
tapers = nx.zeros(N*k)
lambdas = nx.zeros(k)
test = (nx.random.randn(N) * 100).astype('h')
out = nx.zeros(N)

# compute fft
lsono.dpss(tapers, lambdas, N, NW, k)
mtm = lsono.init_mtm_prealloc(N, N, k, tapers, lambdas)
sigpow = lsono.mtfft(mtm, test, N)
lsono.mtpower(mtm, out, sigpow)

# test against python code
## from dlab import signalproc
## lambda2,tapers2 = signalproc.dpss(N,NW,k)
## J,f = signalproc.mtfft(test,nfft=N,mtm_p=NW,tapers=k)
## out2 = nx.absolute(J**2).mean(1)
