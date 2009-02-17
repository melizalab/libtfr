
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
lsono.gettapers.argtypes = [cx.c_void_p, ndpointer(dtype='d')]

lsono.hermf.argtypes = [cx.c_int, cx.c_int, cx.c_double,
                        ndpointer(dtype='d'), ndpointer(dtype='d'), ndpointer(dtype='d')]
lsono.mtm_init_herm.restype = cx.c_void_p
lsono.mtm_init_herm.argtypes = [cx.c_int, cx.c_int, cx.c_int, cx.c_double]

lsono.tfr_displacements.argtypes = [cx.c_void_p, ndpointer(dtype='d'),
                                    ndpointer(dtype='d'), ndpointer(dtype='d')]

lsono.tfr_reassign.argtypes = [ndpointer(dtype='d'),ndpointer(dtype='d'),
                               ndpointer(dtype='d'),ndpointer(dtype='d'),
                               cx.c_int, cx.c_int, cx.c_double, cx.c_double, cx.c_double,
                               cx.c_int, cx.c_int]

lsono.tfr_spec.argtypes = [cx.c_void_p, ndpointer(dtype='d'),ndpointer(dtype='h'),
                           cx.c_int, cx.c_int, cx.c_double, cx.c_int]

lpck = cx.cdll.LoadLibrary('/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib')
lpck.dptsv.argtypes = [cx.POINTER(cx.c_int), cx.POINTER(cx.c_int), ndpointer(dtype='d'),
                        ndpointer(dtype='d'), ndpointer(dtype='d'),
                        cx.POINTER(cx.c_int), cx.POINTER(cx.c_int)]
lpck.dgtsv.argtypes = [cx.POINTER(cx.c_int), cx.POINTER(cx.c_int), 
                       ndpointer(dtype='d'),ndpointer(dtype='d'),
                       ndpointer(dtype='d'), ndpointer(dtype='d'),
                       cx.POINTER(cx.c_int), cx.POINTER(cx.c_int)]
def hc2cmplx(X):
    N = X.shape[0]
    real_count = N / 2 + 1;
    imag_count = (N+1) / 2;
    out = X[:real_count].astype('D')
    out[1:imag_count].imag = X[:(N-imag_count):-1]
    return out

def mwindow(x,w):
    out = nx.zeros(x.size,'d')
    out[:w.size] = x[:w.size] * w
    return out

def alt_tridisolve(ee,d,e):
    n = cx.c_int(d.size)
    nrhs = cx.c_int(1)
    ldb = cx.c_int(d.size)
    info = cx.c_int(0)

    #lpck.dptsv_(cx.byref(n), cx.byref(nrhs), d, ee, e, cx.byref(ldb), cx.byref(info))
    #lpck.dptsv(n, nrhs, d, ee, e, ldb, info)
    lpck.dgtsv(n,nrhs,ee,d,ee.copy(),e,ldb,info)
    return info

if __name__=="__main__":
    N = 256
    NW = 3.5
    k = 5

##     from dlab import signalproc
##     ee,d,e = signalproc.dpss(128,3.5)
##     ee1,d1,e1 = signalproc.dpss(128,3.5)
##     ee = ee[1:]
##     alt_tridisolve(ee,d,e)

    # allocate storage
    #tapers = nx.zeros(N*k)
    #lambdas = nx.zeros(k)
    from dlab import pcmio
    s = pcmio.sndfile('A7.pcm').read()
    #test = s[:N]  #(nx.random.randn(N) * 100).astype('h')
    test = s[8300:(8300+N)]
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

    k = 1
    tm = 6.0
    Np = 201
    tapers = nx.zeros((Np,k*3),order='F')
    #Dh = nx.zeros_like(h)
    #Th = nx.zeros_like(h)

    #lsono.hermf(Np,k,tm,h,Dh,Th)
    mtmh = lsono.mtm_init_herm(N,Np,k,tm)
    lsono.gettapers(mtmh, tapers)
##     h = tapers[:,0]
##     Dh = tapers[:,1]
##     Th = tapers[:,2]

    from dlab import tfr
    h1,Dh1,Th1 = tfr.hermf(Np,k,tm)
    #q1,te1,fe1 = tfr.tfrrsph(test,256,h1[0],Dh1[0],Th1[0])
    spec2 = tfr.tfrrsph(s,256,h1[0],Dh1[0],Th1[0])
    

##     S = tfr.sfft.fft(mwindow(test,h), N)
##     tf2 = tfr.sfft.fft(mwindow(test,Th), N)
##     tf3 = tfr.sfft.fft(mwindow(test,Dh), N)

##     td1 = nx.real(tf2 / S) 
##     fd1 = nx.imag(tf3 / S / (2 * nx.pi))

    sigpow = lsono.mtfft(mtmh,test,Np)

##     buf = nx.zeros((N,k*3),order='F')
##     lsono.getbuffer(mtmh, buf)
##     Z = hc2cmplx(buf)
##     td2 = nx.real(Z[:,2] / Z[:,0])
##     fd2 = nx.imag(Z[:,1] / Z[:,0] / (2 * nx.pi))

    q = nx.zeros((N/2+1,k),order='F')
    td = nx.zeros_like(q)
    fd = nx.zeros_like(q)
    lsono.tfr_displacements(mtmh, q, td, fd)

    qthresh = 299341e-6
    spec = nx.zeros((N/2+1,s.size / 10), order='F')
    #lsono.tfr_reassign(spec, q, td, fd, 129, 129, 10, qthresh, 0.01, 5, 5)
                      
    lsono.tfr_spec(mtmh, spec, s, s.size, 10, 0.01, 5)
