c
c  Gamma and related functions
c


      real function Gamma(z)
      real z
      Gamma = exp(Gammalog(z))
      return
      end
      

c Log de la fonction gamma, estimee a partir de la formule (6.1.5)
c du Numerical Recipes de Press et al. (1988), derivee de l'approximation
c de C. Lanczos, J. SIAM, Num. Anal ser B, v 1, p. 86.
c  (transcription TRIVIALE)

      real function Gammalog(z)
      real z
      double precision x, prod, sum, cof(6), pi
      data cof /76.180091729406, -86.505320327112, 24.014098222230,
     $     -1.231739516140, 0.120858003e-2, -0.536382e-5/
      data pi /3.14159265358/
      logical small
c
c  for z < 1, use the reflection formula (Num. Rec., form. 6.1.4)
      if (z .lt. 1.d0) then
         small = .true.
         zz = 2.d0 - z
      else
         small = .false.
         zz = z
      end if

      x = zz - 1.d0
      prod = x+5.5d0
      prod = (x+0.5d0)*log(prod) - prod
      sum = 1.000000000178
      do i=1,6
         x = x+1.d0
         sum = sum + cof(i)/x
      enddo
      Gl = log(2.50662827465d0*sum) + prod

      if (small) then
         zz = 1. - z
         Gammalog = log(pi*zz) - Gl - log(sin(pi*zz))
      else
         Gammalog = Gl
      end if

      return
      end


c
c  Qgamma - compute the incomplete gamma function Q(a,x) 
c             Appropriate for x > 0.  Preferred for x > a+1
c
c  Based on eq. 6.2.6 of Press et al.
c
c     Gam(a,x) =  e**-x * x**a (1/x+ (1-a)/1+ 1/x+ (2-a)/1+ 2/x+ ...)
c
      Real Function Qgamma(a,x)
      real a, x
      double precision OddA,OddB,EvenA,EvenB,even,odd,fodd,feven
      external Pgamma

      if (x .eq. 0.d0) then
         Qgamma = 1.0
         return
      end if
      if (x .lt. a + 1.d0) then
         Qgamma = 1. - Pgamma(a,x)
         return
      endif

c  Solve by continued fractions
c
      epsilon = 1.d-6

      OddA = 1
      EvenA = 1
      OddB = x
      EvenB = x + 1 - a
      do odd = 1.d0, 49.d0
         even = odd+1.d0 - a
         OddA  = x*EvenA + odd*OddA
         OddB  = x*EvenB + odd*OddB
         EvenA =    OddA + even*EvenA 
         EvenB =    OddB + even*EvenB
         fodd  = OddA/OddB
         feven = EvenA/EvenB
         if (abs(1.d0 - fodd/feven) .lt. epsilon) go to 9
      end do
 9    continue
      Qgamma = exp(a*log(x) -x + log(feven) - gammalog(a))
      return

 10                         format (4g18.8)

      end


      real function Pgamma(a,x)
c
c  Incomplete gamma function P(a,x)
c  This algorithm is based on eqs. 6.2.5  and 6.1.3 of Press et al.
c
c  Bugs:
c      Might cause infinite recursion with Qgamma.
c
      real a,x
      double precision sum, term, factor, epsilon
      external Qgamma

      if (x .eq. 0.d0) then
         Pgamma = 0.d0
         return
      end if
      if (x .ge. a+1.d0) then 
         Pgamma = 1. - Qgamma(a,x)
         return
      end if

c     Use the infinite series representation
      denom = a
      epsilon = 1.d-7
      term = 1.d0/denom
      sum = term
      do i = 1,99
         denom = denom + 1.d0
         term = term * x / denom
         sum = sum + term
         if (abs(term/sum) .lt. epsilon) go to 9
      enddo
 9    continue
      Pgamma = exp(a*log(x) - x + log(sum) - gammalog(a))
      return
      end
