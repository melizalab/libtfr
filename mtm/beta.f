c
c  Beta and related functions
c


      Real Function BetaInc(a,b,x)
c
c  Incomplete beta function Ix(a,b).
c  Based on the algorithms of Press et. al, Numerical Recipes, 1st ed.,
c  Sec. 6.3.
c
      if (x.lt.0. .or. x.gt.1.) write(0,*) 'BetaInc: bad argument x'
      if (x.eq.0. .or. x.eq.1.) then
        bt = 0.
      else
        bt = exp(Gammalog(a+b)-Gammalog(a)-Gammalog(b)
     *      +a*log(x)+b*log(1.-x))
      end if
      if (x .lt. (a+1.)/(a+b+2.)) then
        BetaInc = bt*Betafrac(a,b,x)/a
      else
        BetaInc = 1.-bt*Betafrac(b,a,1.-x)/b
      end if
      return
      end


      Real Function Betafrac(a,b,x)
c
c  algorithm:
c    compute the continued fraction expansion of Ix(a,b)
c      Do the recursion relation in two steps
c      Renormalize to prevent over- or underflows
c
      parameter (epsilon=3.e-7)

      am = 1.
      bm = 1.
      az = 1.
      bz = 1.-(a+b)*x/(a+1.)
      do m = 1, 100
        em = m
        twom = em+em
        d = em*(b-em)*x/((a-1.d0+twom)*(a+twom))
        Aodd = az+d*am
        Bodd = bz+d*bm
        d = -(a+em)*(a+b+em)*x/((a+twom)*(a+1.+twom))
        Aeven = Aodd+d*az
        Beven = Bodd+d*bz
        aold = az
        am = Aodd/Beven
        bm = Bodd/Beven
        az = Aeven/Beven
        bz = 1.
        if (abs(az-aold) .lt. epsilon*abs(az)) go to 9
      end do
      write(0,*) 'Too many iterations.  Perhaps a or b too big'
 9    Betafrac = az
      return
      end
