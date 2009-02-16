
import numpy as nx
import ctypes as cx
from numpy.ctypeslib import ndpointer

lsono = cx.cdll.LoadLibrary('libsono2.dylib')
#lsono1 = cx.cdll.LoadLibrary('libsono.dylib')

lsono.dpss.restypes = None
lsono.dpss.argtypes = [ndpointer(dtype='d'), ndpointer(dtype='d'),
                       cx.c_int, cx.c_double, cx.c_int]

lsono.mtm_init.restype = cx.c_void_p
lsono.mtm_init.argtypes = [cx.c_int, cx.c_int, cx.c_int,
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]
lsono.mtm_init_dpss.restype = cx.c_void_p
lsono.mtm_init_dpss.argtypes = [cx.c_int, cx.c_double, cx.c_int]

lsono.mtfft.restype = cx.c_double
lsono.mtfft.argtypes = [cx.c_void_p, ndpointer(dtype='h'), cx.c_int]

lsono.mtpower.argtypes = [cx.c_void_p, ndpointer(dtype='d'), cx.c_double]

lsono.create_sonogram.restype = cx.c_void_p
lsono.sonogram_setopts.argtypes = [cx.c_void_p, cx.c_int, cx.c_long]

lsono.getbuffer.argtypes = [cx.c_void_p, ndpointer(dtype='d')]

lsono.hermf.argtypes = [cx.c_int, cx.c_int, cx.c_double,
                        ndpointer(dtype='d'), ndpointer(dtype='d'), ndpointer(dtype='d')]
lsono.mtm_init_herm.restype = cx.c_void_p
lsono.mtm_init_herm.argtypes = [cx.c_int, cx.c_int, cx.c_int, cx.c_double]

lsono.tfr_displacements.argtypes = [cx.c_void_p, ndpointer(dtype='d'),
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]

def hc2cmplx(X):

    N = X.shape[0]
    real_count = N / 2 + 1;
    imag_count = (N+1) / 2;
    out = X[:real_count].astype('D')
    out[1:imag_count].imag = X[:(N-imag_count):-1]
    return out


if __name__=="__main__":
    N = 256
    NW = 3.5
    k = 5

    # allocate storage
    #tapers = nx.zeros(N*k)
    #lambdas = nx.zeros(k)
    from dlab import pcmio
    s = pcmio.sndfile('A7.pcm').read()
    #test = s[:N]  #(nx.random.randn(N) * 100).astype('h')
    test = s[6800:(6800+N)]
    out = nx.zeros(N/2+1)
    out2 = nx.zeros_like(out)
    out3 = nx.zeros_like(out)

    # compute fft
    #lsono.dpss(tapers, lambdas, N, NW, k)
    #mtm = lsono.mtm_init(N, N, k, tapers, lambdas)
    mtm = lsono.mtm_init_dpss(N, NW, k)
    #for i in range(100):
    sigpow = lsono.mtfft(mtm, test, N)
    lsono.mtpower(mtm, out, sigpow)
    lsono.mtpower(mtm, out2, 0.0)
    #lsono.destroy_mtm(mtm)

    # test with single window
    window = nx.hanning(N)
    ll = nx.ones(1)
    mtm2 = lsono.mtm_init(N, N, 1, window, ll)
    #for i in range(100):
    sigpow = lsono.mtfft(mtm2, test, N)
    lsono.mtpower(mtm2, out3, 0.0)

    # test against python code
##     from dlab import signalproc
##     lambda2,tapers2 = signalproc.dpss(N,NW,k)
##     J,f = signalproc.mtfft(test,tapers=tapers2)
##     out4 = nx.absolute(J**2).mean(1)

    # test time-freq reassignemnt

    k = 6
    tm = 6.0
    Np = 201
    h = nx.zeros((Np,k),order='F')
    Dh = nx.zeros_like(h)
    Th = nx.zeros_like(h)

    lsono.hermf(Np,k,tm,h,Dh,Th)

##     from dlab import tfr
##     h1,Dh1,tt1 = tfr.hermf(N,k,tm)
    mtmh = lsono.mtm_init_herm(N,Np,k,tm)
    sigpow = lsono.mtfft(mtmh,test,Np)

    buf = nx.zeros((N,k*3),order='F')
    lsono.getbuffer(mtmh, buf)

    q = nx.zeros((N/2+1,k),order='F')
    td = nx.zeros_like(q)
    fd = nx.zeros_like(q)
    lsono.tfr_displacements(mtmh, q, td, fd)

    #Z = hc2cmplx(buf[:,:6])
    #lsono.tfr_cmplx(mtmh,Y)


