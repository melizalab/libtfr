      subroutine spec(a0,nwin0,ip,dt0,nscan0,
     $   flow,fhigh,
     $   inoise0,ismooth,ilog0,fsmooth,
     $   isignal,ispec,inorm,iresh,iper,
     $   nf0,df,specraw0,specmed0,specresh0,harmonic0,specback0,
     $   ftest,conf,fconf,nbnd,
     $   whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0)
c
c     (c) Michael Mann
c
c     $Id: spec.f,v 1.5 1997/04/29 00:30:41 weibel Exp $
c
c     MTM procedure of Mann and Lees (1996) adapted fromm original MTM code
c     of J. Park, to perform multiple null-hypothesis testing, both harmonic
c     and quasiperiodic signal detection, and robust background noise
c     estimation.
c
      integer maxwin,maxlen,maxdof,nlim,nmed
      real big,small,tol
      parameter (maxlen=50000)
      parameter (maxwin=8,maxdof=2*maxwin,nlim=32768)
      parameter (nmed=2000)
      parameter (big=1.0e+36,small=1.0e-12,tol=1e-4)
c
      real dt0
      real a0(maxlen)
      integer ip,nscan0,nwin0,ilog0
      integer inoise,inoise0
      integer nbnd,nbnd0
      integer idone(nlim)
      real rho,rho0
      real MSE
      real*8 avar,el
      real val(nmed)
      real demean(nlim),specmed
      real specraw0(nlim)
      real specresh0(nlim),harmonic0(nlim),specback0(nlim),
     $      specmed0(nlim)
      real rednoise(nlim),rednoise0(nlim)
      real ftest(nlim)
      real base(nlim)
      integer iharm(nlim)
      real fconf(6)
      real conf(4)
      real adum,sumpeak
      complex zed
c
c     note that following limits are hard-wired into the code:
c       max # tapers = 8
c       max # frequency paris = 16400
c       max # points in dataseries = 32800
c
c
      common/taperz/ta(16400,8),tai(16400,8),dcf(16400,8),
     $   amu(16400,8)
      common/staple/tad(32800),tad1(32800)
      common/stap2/tas(32800,16)
      common/misc/npad
      common/npiprol/anpi
      common/data/a(maxlen)
      common/work/b(32800)
      common/work2/fspec(16400)
      common/BLK1/white,specmed(nlim),fnynew,ddf,f0
      common/BLK2/ilog,nfreq
c
      dimension dum(35)
      real *8 b1(32800),b2(32800)
      dimension reim(2),el(10)
      equivalence (zed,reim),(iy,dum)
c
      real chisqdata(5,16)
      data chisqdata /  0.455, 2.706,  3.841,  6.635, 10.827,
     2                  0.693, 2.303,  2.996,  4.605, 6.908,
     3                  0.789, 2.084,  2.605,  3.780, 5.423,
     4                  0.839, 1.945,  2.372,  3.319, 4.617,
     5                  0.870, 1.847,  2.214,  3.017, 4.102,
     6                  0.891, 1.774,  2.099,  2.802, 3.743,
     7                  0.907, 1.717,  2.010,  2.639, 3.475,
     8                  0.918, 1.670,  1.938,  2.511, 3.266,
     9                  0.927, 1.632,  1.880,  2.407, 3.097,
     $                  0.934, 1.599,  1.831,  2.321, 2.959,
     1                  0.940, 1.570,  1.789,  2.248, 2.842,
     2                  0.945, 1.546,  1.752,  2.185, 2.742,
     3                  0.949, 1.524,  1.720,  2.130, 2.656,
     4                  0.953, 1.505,  1.692,  2.082, 2.580,
     5                  0.956, 1.487,  1.666,  2.039, 2.513,
     6                  0.959, 1.471,  1.644,  2.000, 2.453
     $                /

      real ftestdata(6,7)
      data ftestdata  /
     2                 1.00,  9.00, 19.0, 99.0,  199.0, 99.0,
     4                 0.828, 4.32, 6.94, 18.0,  26.28, 61.25,
     6                 0.780, 3.46, 5.14, 10.92, 14.54, 27.0,
     8                 0.757, 3.11, 4.46,  8.65, 11.04, 18.49,
     $                 0.743, 2.92, 4.1,   7.56,  9.43, 14.91,
     2                 0.735, 2.81, 3.89,  6.93,  8.51, 12.97,
     4                 0.729, 2.73, 3.74,  6.51,  7.92, 11.78
     $              /

c
c     define some constants
c
      pi=3.14159265358979
      radian=180./pi
      tiny = 1e-6
c
c     asign parameters from main
c
      nscan = nscan0
      do i=1,maxlen
         a(i)=a0(i)
      end do
      npi = ip
      anpi=float(ip)
      nwin=nwin0
      dt = dt0
      nscan = nscan0
      ilog = ilog0
c
c     zero pad to first power of 2 > # data points
c     constrain minimum padding to 128 points...
c
      npts = nscan
      npad=npts-1
      ij=0
c
c     set zero padding to first power of 2 greater than
c     number of data points
c
c     impose minimum value of 128
c
1000  npad=npad/2
      ij=ij+1
      if(npad.ne.0) go to 1000
c
      npad=max(2**ij,128)
c
c     determine frequency limits, associated DFT points, etc.
c
      fmin=flow
      fmax=fhigh
      fnynew = fhigh
      ddf=1.0/(npad*dt)
      df = ddf
      eps=.1*ddf
      n1=(fmin+eps)/ddf
      n2=((fmax+eps)/ddf)
      fmin=n1*ddf
      fmax=n2*ddf
      nf=n2-n1+1
      nf0 = nf
      nfreq=nf
      f0 = fmin
      nmin=n1*2+1
      nmax=n2*2+2
c
c     determine Rayleigh, Nyquist frequencies, bandwidths, etc.
c
      fny = 0.5/dt
      fray=1.0/(npts*dt)
      bndwdth = 2.0*anpi*fray
      halfwdth = bndwdth/2.0
      nbnd = bndwdth/ddf
      nbnd0 = fray/ddf
      timefreq = 2.0*anpi
c
c     demean the series
c
      demn = 0.0
      do i=1,nscan
         demn=demn+a(i)
      end do
      demn = demn/float(nscan)
      do i=1,nscan
         demean(i)=a(i)-demn
      end do
c
c
c     determine raw noise parameters
c
c
      inoise = inoise0
      if (inoise.gt.0) then
         rho0 = 0.0
      else
c
c      determine raw lag 1 autocorrelation coefficient
c
       var = 0.0
       do i=1,nscan
         var = var + demean(i)**2
       end do
       var = var/float(nscan)
       sd = sqrt(var)
       c1 = 0.0
       icount = 0
       do i=2,nscan
         icount = icount + 1
         c1 = c1 + demean(i-1)*demean(i)
       end do
       c1 = c1/float(icount)
       rho0 = c1/var
      endif
c
c     determine chi-square values for confidence level
c     determination of spectrum
c
      idofs = 2*nwin
      if (idofs .le. maxdof) then      
         do i=1,4
            conf(i) = chisqdata(i,idofs-1)
         enddo
      endif
c
      iftest = 1
      jper = 3
c
c     pick f-test data if appropriate
c
      if (nwin.eq.1) iftest = 0
      if (iftest.eq.1) then
         i = nwin-1
         if (i .le. 7) then
            do j = 1,6
               fconf(j) = ftestdata(j,i)
            enddo
         endif
         thresh = fconf(iper)
      endif
c
c     determine amplitude theshold for reshaping
c
      afact = 0.0
      if ((iresh.eq.1).and.(isignal.ne.2)) then
         afact = conf(1)
         if (iper.eq.2) afact = conf(2)
         if (iper.eq.3) afact = conf(3)
         if (iper.gt.3) afact = conf(4)
      endif
c
c     create output files
c
c     DETERMINE SPECTRUM
c
c     calculate the eigentapers
c
      call taper2(npts,nwin,el)
c
c     normalization:
c
c     mult by dt if we wish absolute spectral estimate
c         e.g. analysis of time-limited signal
c     divide by npts if we wish amplitude spectrum per unit time
c
      if (inorm.eq.2) then  !  amp spec (time limited)
        anrm=1.0/dt
      else
        if(inorm.eq.1) then  ! amp spect per unit time
          anrm=float(npts)
        else
          anrm=1.
        endif
      endif
c
c     perform convolution of the series with each datataper
c
      do iwin=1,nwin
        do i=1,nscan
          b1(i)=demean(i)*tas(i,iwin)
          b2(i)=0.0
        end do
        j=nscan+1
        do i=j,npad
          b1(i)=0.0
          b2(i)=0.0
        end do
        call fft2(b1,b2,npad)
        sum=0.
        j = 0
        do i=1,npad/2
          j = j+1
          ta(j,iwin)=b1(i)/anrm
          tai(j,iwin)=b2(i)/anrm
          b(i) = b1(i)
          b(i+1) = b2(i)
          sumi=b1(i)*b1(i)+b2(i)*b2(i)
          sum=sum+sumi
        end do
      end do
c
      if (ispec.eq.0) then
c
c       calculate "high-resolution" spectrum
c
        call hires(ta,tai,el,nwin,nf-1,amu)
c
c       determine psd as squared amplitude spectrum
c       dissallow zero values
c
        do i=1,nf-1
          amu(i,1)=abs(amu(i,1))**2
        end do
      endif
c
      if (ispec.eq.1) then
c
c       calculate adaptively weighted spectrum
c
        avar=0.d0
        do i=1,nf-1
          avar=avar+demean(i)*demean(i)
        end do
c
c       avar is a factor that scales the bias factor in the adaptive
c       weighting scheme.
c
        if (inorm.eq.1) then ! amp spect per unit time
          avar=avar/(npts*npts)
        elseif(inorm.eq.2) then  ! absolute amp spect
          avar=avar*dt*dt
        endif
c
        call adwait(ta,tai,dcf,el,nwin,nf-1,amu,amu(1,2),avar)
c
c       normalized psd, disallow negative values
c
        do i=1,nf-1
          amu(i,1)=abs(amu(i,1))
        end do
      endif
c
c     initialize raw, reshaped, and harmonic spectra as the psd estimated
c     above, and the f-test as null
c
      do i=1,nf-1
         iharm(i)=0
         specraw0(i)=amu(i,1)
         specresh0(i)=amu(i,1)
         harmonic0(i)=amu(i,1)
         idone(i)=0
      end do
c
      if (iftest.eq.1) then
c
c
c       perform f-test for phase coherence if indicated
c
        call regre(ta,tai,nf-1,nwin,amu)
c
c
        do i=1,nf-1
           ftest(i)=amu(i,3)
        end do
      endif
c
c     determine the average (white) power level of the raw spectrum
c
      white0 = 0.0
      do i=1,nf-1
         white0 = white0+specraw0(i)
      end do
      white0 = white0/float(nf-1)
c
c     set default noise parameters as raw noise parameters
c
      rho = rho0
      white = white0
c
      if (ismooth.eq.1) then
c
c      determine median smoothed spectrum and associated
c     "robust" average (white) power level
c
       white = 0.0
       nsmooth = fsmooth/ddf
       do j=1,nf-1
         if1 = j-nsmooth/2
         if2 = j+nsmooth/2
         if (if1.lt.1) then
           if1=1
         endif
         if (if2.gt.nf-1) then
           if2=nf-1
         endif
         nblk = if2-if1+1
         do i=1,nblk
            val(i)=specraw0(if1+i-1)
         end do
c
c        sort spectrum in this block
c
         kmid = (nblk+1)/2
         do kk1=1,nblk
            do kk2=kk1+1,nblk
               if (val(kk2).lt.val(kk1)) then
                  adum = val(kk2)
                  val(kk2)=val(kk1)
                  val(kk1)=adum
               endif
            end do
         end do
         specmed(j)=val(kmid)
         specmed0(j)=specmed(j)
         white = white + specmed(j)
       end do
       white = white/float(nf-1)
c
       if (inoise.ne.2) then
c
c         Unless user has selected the "local white" assumption,
c         we attempt to fit a parametric (white or red) noise
c         background robustly
c
c         determine best fit spectrum of the form
c
c         rednoise = white*(1.0-rho**2)/
c                     (1.0-2.0*rho*cos(arg)+rho**2)
c
c         to the median-smoother.
c
c         note that the "white noise" is  a trivial case
c
c        "white" is the robust white noise power level
c         estimated above from the median smoothed
c         spectrum
c
c         do a global search of the interval [0,1)
c         to find the optimal rho as determined by minimum MSE
c
          amin = big
          do rho1=0.0,0.999,0.001
           amiss = MSE(rho1)
           if (amiss.lt.amin) then
            rho=rho1
            amin = amiss
           endif
          end do
       endif
      endif
c
c     determine raw and robust estimates of the decorrelation
c     timescale of the noise
c
      tau = 1e+16
      tau0 = 1e+16
      if (rho.gt.0.0) tau = -dt/log(rho)
      if (rho0.gt.0.0) tau0 = -dt/log(rho0)
c
c     determine the noise background
c
      do i=1,nf-1
         ff = fmin+(i-1)*ddf
         freq_norm = ff/fnynew
         arg = freq_norm*pi
         rednoise0(i) = white0*(1.0-rho0**2)/
     $      (1.0-2.0*rho0*cos(arg)+rho0**2)
         rednoise(i) = white*(1.0-rho**2)/
     $      (1.0-2.0*rho*cos(arg)+rho**2)
c
c        if the "harmonic signal" option is selected,
c        the estimated null "base" spectrum is zero,
c        otherwise it is the noise background determined above
c
         if (isignal.eq.2) then
            base(i)=0.0
         else
           if (inoise.lt.2) then
             if (ismooth.eq.1) then
                base(i) = rednoise(i)
             else
                base(i) = rednoise0(i)
             endif
           else
             base(i) = specmed(i)*conf(1)
           endif
           specback0(i)=base(i)
         endif
      end do
c
c     now perform reshaping procedure if indicated
c
      if (iresh.eq.1) then
c
c       do a poor mans reshaping -detect harmonic peaks and then
c       interpolate the continuous spectrum across the effected
c       bandwidth  -only reshape if harmonic peak is greater than
c       the significance level in terms of overall power that was
c       indicated for the F-test harmonic detection procedure
c
c       note that for harmonic signal assumption, base=0 and all
c       signficant f-test peaks will be reshaped
c
c       note: reshaping is only done at frequencies outside
c       the secular band
c
        do i=nbnd,nf-1
           if ((ftest(i).gt.thresh).
     $         and.(specraw0(i).gt.afact*base(i))) iharm(i)=1
        end do
c
        do i=nbnd,nf-1
c
c           determine frequency points at boarder of the Rayleigh bandwidth
c           and the full spectral bandwidth
c
c           if frequency falls in a band within reshaping was already
c           performed, we skip
c
            if (idone(i).ne.1) then
c
            ipre0 = i-nbnd0/2-1
            iaft0 = i+nbnd0/2+1
            ipre = i-nbnd/2
            iaft = i+nbnd/2
            if (ipre.lt.nbnd+1) ipre=nbnd+1
            if (iaft.gt.nf-2) iaft=nf-2
            if (ipre0.lt.nbnd) ipre0=nbnd
            if (iaft0.gt.nf-1) iaft0=nf-1
c
c           reshaped spectrum is estimating by assuming that the
c           spectrum is continuous across the Rayleigh spectral
c           bandwidth within which a periodic signal is detected.
c           this gap is linearly interpolated across that bandwidth
c
c           "harmonic0" is the sum of the reshaped "continuous background"
c           and the detected spectral line. The line component is estimated
c           by assuming that the total power within the reshaped region
c           is represented by the sum of the continuous background and
c           a narrow peak with the Rayleigh resolution
c
            sumpeak = 0.0
            if (iharm(i).eq.1) then
             do j=ipre,iaft
              specresh0(j)=0.5*(specraw0(ipre-1)+specraw0(iaft+1))
              harmonic0(j)=specresh0(j)
             end do
             do j=ipre,iaft
               idone(j)=1
               sumpeak = sumpeak+specraw0(j)
             end do
             do j=ipre0,iaft0
                harmonic0(j)=sumpeak/(iaft0-ipre0+1)
             end do
            endif

            endif

        end do
      endif
c
c     for periodic signal assumption, the estimated noise background
c     is simply the reshaped (ie, estimated continuous) spectrum
c
      if (isignal.eq.2) then
         do i=1,nf-1
            specback0(i)=specresh0(i)
         end do
      endif
c
      whiteraw0=white0
      whiterob0=white
      rhoraw0=rho0
      rhorob0=rho
      tauraw0 = tau0
      taurob0 = tau
c
8888  continue
c
      return
      end
c
      real function MSE(rho)
c
      parameter (nlim=32768)
      COMMON /BLK1/white,specmed(nlim),fnynew,ddf,f0
      COMMON /BLK2/ilog,nfreq
      real rho
      real pie,dff,freq_norm,ff,arg,small,rednoise,val1,val2
      pie = 3.1415926
      small = 1e-12
      dff = 0.0
      do j=1,nfreq
         ff = f0+(j-1)*ddf
         freq_norm = ff/fnynew
         arg = freq_norm*pie
         rednoise = white*(1.0-rho**2)/
     $      (1.0-2.0*rho*cos(arg)+rho**2)
         if (ilog.eq.0) then
           dff = dff + (specmed(j)-rednoise)**2
         else
           val1 = abs(specmed(j))
           val2 = abs(rednoise)
           if (val1.lt.small) val1=small
           if (val2.lt.small) val2=small
           dff = dff +
     $     (log(val1)-log(val2))**2
         endif
      end do
      MSE = dff
      return
      end
c
      subroutine hires(ta,tai,el,nwin,nf,ares)
      real*8 el
      dimension ta(16400,1),tai(16400,1),el(1),ares(1)
      do j=1,nf
        ares(j)=0.
      end do
      do i=1,nwin
        a=1./(el(i)*nwin)
        do j=1,nf
          ares(j)=ares(j)+a*(ta(j,i)*ta(j,i)+tai(j,i)*tai(j,i))
        end do
      end do
      do j=1,nf
        ares(j)=sqrt(ares(j))
      end do
      return
      end
c
c
      subroutine adwait(ta,tai,dcf,el,nwin,nf,ares,degf,avar)
c
c  this version uses Thomson's algorithm for calculating
c  the adaptive spectrum estimate
c
      real*8 avar,spw,as,das,tol,el,a1,bias,scale,ax,fn,fx
      dimension ta(16400,1),tai(16400,1),el(1),ares(1),degf(1)
      dimension spw(10),bias(10),dcf(16400,1)
c
c  set tolerance for iterative scheme exit
c
      tol=3.d-4
      jitter=0
      scale=avar
c
c  we scale the bias by the total variance of the frequency transform
c  from zero freq to the nyquist
c  in this application we scale the eigenspectra by the bias in order to avoid
c  possible floating point overflow
c
      do 200 i=1,nwin
  200 bias(i)=(1.d0-el(i))
      do 100 j=1,nf
        do 150 i=1,nwin
  150   spw(i)=(ta(j,i)*ta(j,i)+tai(j,i)*tai(j,i))/scale
c
c  first guess is the average of the two lowest-order eigenspectral estimates
c
        as=(spw(1)+spw(2))/2.d0
        do 300 k=1,20
c
c  find coefficients
c
          fn=0.d0
          fx=0.d0
          do 350 i=1,nwin
            a1=dsqrt(el(i))*as/(el(i)*as+bias(i))
            a1=a1*a1
            fn=fn+a1*spw(i)
            fx=fx+a1
  350     continue
          ax=fn/fx
          das=dabs(ax-as)
          if(das/as.lt.tol) go to 400
  300   as=ax
c
c  flag if iteration does not converge
c
      jitter=jitter+1
  400 continue
      ares(j)=as*scale
c
c  calculate degrees of freedom
c
      df=0.
      do 450 i=1,nwin
      dcf(j,i)=dsqrt(el(i))*as/(el(i)*as+bias(i))
  450 df=df+dcf(j,i)*dcf(j,i)
c
c  we normalize degrees of freedom by the weight of the first eigenspectrum
c  this way we never have fewer than two degrees of freedom
c
  100 degf(j)=df*2./(dcf(j,1)*dcf(j,1))
      return
      end
c
      subroutine regre(sr,si,nf,nwin,amu)
c
      real b
      real *8 junk
      dimension sr(16400,1),si(16400,1),amu(16400,8)
      common/tapsum/b(10),junk(10)
c
c     "b" is the DFT of Slepian eigentapers at zero frequency
c     "sr" and "si" are the eigenspectra
c     "amu" contains line frequency estimates and f-test parameter
c
      sum=0.
      do i=1,nwin
        sum=sum+b(i)*b(i)
      end do
      do i=1,nf
        amu(i,1)=0.
        amu(i,2)=0.
        do j=1,nwin
          amu(i,1)=amu(i,1)+sr(i,j)*b(j)
          amu(i,2)=amu(i,2)+si(i,j)*b(j)
        end do
        amu(i,1)=amu(i,1)/sum
        amu(i,2)=amu(i,2)/sum
        sum2=0.
        do j=1,nwin
          sumr=sr(i,j)-amu(i,1)*b(j)
          sumi=si(i,j)-amu(i,2)*b(j)
          sum2=sum2+sumr*sumr+sumi*sumi
        end do
        amu(i,3)=(nwin-1)*(amu(i,2)**2+amu(i,1)**2)*sum/sum2
      end do
      return
      end
c
c
      subroutine taper2(n,nwin,el)
c
c     generate slepian tapers
c     "ta" is a real*4 array
c
c     written by J. Park
c
      real*8 el,a,z,pi,ww,cs,ai,an,eps,rlu,rlb
      real*8 dfac,drat,gamma,bh,ell
      real dump
      common/npiprol/anpi
      common/tapsum/tapsum(10),ell(10)
      common/misc/npad
      common/work/ip(32800)
      common/taperz/z(65536),dump(393728)
      common/stap2/ta(32800,16)
      dimension a(32800,8),el(10)
      data pi/3.14159265358979d0/
      equivalence (a(1,1),ta(1,1))
      an=dfloat(n)
      ww=dble(anpi)/an
      cs=dcos(2.d0*pi*ww)
c
c     note:
c
c initialize matrix for eispack subroutine
c this matrix is not the bandwidth retention factor matrix A (for minimizing
c spectral leakage) described in various multitaper papers
c --e.g. Thomson (1982) [proc ieee], Park et al (1987) [jgr].
c the bandwidth matrix x returns a cluster of eigenvectors
c (which are the slepian tapers) with eigenvalues very close to unity.
c since the spacing of eigenvalues can be comparable to machine precision,
c numerical instability can occur for time-bandwidth product n .ge. 6
c (that is, 6pi-prolate tapers)
c also, the bandwidth retention matrix A is toeplitz, but full, and it is not
c feasible to calculate tapers explicitly for long time series (in an earlier
c code i interpolated an m-point taper from the eigenvectors of a 128x128
c bandwidth retention matrix). the following is cribbed from a 1978 paper
c by Slepian (Prolate spheroidal wave functions, Fourier analysis
c and uncertainty, Bell System Tech Journal, v57, pp1371-1430, 1978)
c which gives a three-term recursion that is satisfied by the
c slepian tapers. the three-term recursion can be manipulated to show that the
c slepian tapers are the eigenvectors of a tridiagonal matrix a, with
c well-spaced eigenvalues that are (practically speaking) unrelated to the
c eigenvalues of x. using this matrix a we obtain both numerical stability
c and speed (tridiagonal matrices can be decomposed rapidly enough to
c calculate m-point tapers on the fly).
c
      do i=0,n-1
        ai=dfloat(i)
        a(i+1,1)=-cs*((an-1.d0)/2.d0-ai)**2
        a(i+1,2)=-ai*(an-ai)/2.d0
c
c       next statement is eispack routine tridib, see its documentation
c
        a(i+1,3)=a(i+1,2)**2
      end do
      eps=1.e-13
      m11=1
      call tridib(n,eps,a(1,1),a(1,2),a(1,3),rlb,rlu,m11,nwin,el,ip,
     x       ierr,a(1,4),a(1,5))
      call tinvit(n,n,a(1,1),a(1,2),a(1,3),nwin,el,ip,z,ierr,
     x            a(1,4),a(1,5),a(1,6),a(1,7),a(1,8))
c
c  note:
c
c  we calculate the eigenvalues of the dirichlet-kernel problem
c  i.e. the bandwidth retention factors
c  from slepian 1978 asymptotic formula, gotten from thomson 1982 eq 2.5
c  supplemented by the asymptotic formula for k near 2n from Slepian (1978)
c  more precise values of these parameters, perhaps useful in adaptive
c  spectral estimation, can be calculated explicitly using the
c  rayleigh-quotient formulas in Thomson (1982) and Park et al (1987)
c
      dfac=an*pi*ww
      drat=8.d0*dfac
      dfac=4.d0*dsqrt(pi*dfac)*dexp(-2.d0*dfac)
      do k=1,nwin
        el(k)=1.d0-dfac
        dfac=dfac*drat/k  ! is this correct formula? yes,but fails as k -> 2n
      end do
      gamma=dlog(8.d0*an*dsin(2.d0*pi*ww))+0.5772156649d0
      do k=1,nwin
        bh=-2.d0*pi*(an*ww-dfloat(k-1)/2.d0-.25d0)/gamma
        ell(k)=1.d0/(1.d0+dexp(pi*bh))
      end do
      do i=1,nwin
        el(i)=dmax1(ell(i),el(i))
      end do
c
c  normalize the eigentapers to preserve power for a white process
c  i.e. they have rms value unity
c  "tapsum" is the average of the eigentaper, should be near zero for
c  antisymmetric tapers
c
      do k=1,nwin
        kk=(k-1)*n
        tapsum(k)=0.
        tapsq=0.
        do i=1,n
          aa=z(kk+i)
          ta(i,k)=aa
          tapsum(k)=tapsum(k)+aa
          tapsq=tapsq+aa*aa
        end do
        aa=sqrt(tapsq/n)
        tapsum(k)=tapsum(k)/aa
        do i=1,n
          ta(i,k)=ta(i,k)/aa
        end do
      end do
c
c  the real FFT will preserve amplitudes with zeropadding
c  for example, a(i)=1.,i=1,100 will transform at zero freq
c  to b(f=0)=100 no matter how much zero padding is done
c  therefore we need not doctor the taper normalization,
c  but wait until the fft to force the desired units
c
      return
      end

