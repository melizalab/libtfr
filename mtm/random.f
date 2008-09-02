      Function PUNI (idu)
c     __________________________________________________________________
c     initiate by calling rstart with i,j,k,l = integers but not all = one
c                                               1 <=  i,j,k <= 178
c                                               1 <=  l <= 168
c    
c     from marsaglia, zaman, and tsang (1990), "Toward a universal 
c     random number generator", Stats & Prob Letters, v. 8, p. 35-39.
c
c     Uniform variate generator with return period of about 2**144 !!
c     ___________________________________________________________________

      real u(97)
      common /set1/ u,c,cd,cm,ip,jp
 
      uni=u(ip)-u(jp)
      if (uni.lt.0.) uni=uni+1.
      u(ip)=uni
      ip=ip-1
       if (ip.eq.0) ip=97
      jp=jp-1
      if (jp.eq.0) jp=97
      c=c-cd
      if (c.lt.0.) c=c+cm
      uni=uni-c
      if(uni.lt.0) uni=uni+1.
      PUNI=uni

      return
      end


       Subroutine RSTART (i,j,k,l)
c      _________________________________________________________________
c
c      initialization for marsaglia, zaman, and tsang uniform
c      variate generator in function PUNI
c
c  WARNING:  On DEC computers (at least), calling rstart with fixed arguments
c  WARNING:  will lead to generation of nonrandom variates by subsequent calls
c  WARNING:  to PUNI()
c  WARNING:
c  WARNING:  Therefore, NEVER EVER
c  WARNING:
c  WARNING:     call RSTART (12,34,56,78)
c  WARNING:
c  WARNING:  instead make sure to call rstart in the form that follows
c  WARNING:
c  WARNING:     i1=12 (or whatever seed you want)
c  WARNING:     i2=34
c  WARNING:     i3=56
c  WARNING:     i4=78
c  WARNING:     call RSTART(i1,i2,i3,i4)
c      _________________________________________________________________

       real u(97)
       common /set1/ u,c,cd,cm,ip,jp

       do 2 ii=1,97
       s=0.
       t=0.5
       do 3 jj=1,24
        m=mod(mod(i*j,179)*k,179)
        i=j
        j=k
        k=m
        l=mod(53*l+1,169)
        if (mod(l*m,64).ge.32) s=s+t
  3     t=0.5*t
  2     u(ii)=s
        c=362436./16777216.
        cd=7654321./16777216.
        cm=16777213./16777216.
        ip=97
        jp=33

        return
        end


       Function GNNORM (iran)
c      _________________________________________________________________
c      gaussian random number fcn -  randub/box-muller
c      revision and extension of jim slack's gaussv.    wk 8/75. 4/78. 6/83.
c      uses single precision and calls only randuv, not randub.
c
c      from usgs xlib.f77:  kirby, wh, 1983, computer routines for
c      probability distributions, random numbers, and related functions:
c      us geological survey water-resources investigations report
c      83-4257, page 6.
c      __________________________________________________________________
 
      dimension w(2)
c      data twopi/6.2831853/, fneg2/-2./ !  ,  angl/1e31/, none/.true./
      data twopi/6.2831853/, fneg2/-2./, angl/1e31/
 
 10   w(1)=PUNI(iran)
      w(2)=PUNI(iran)
      dist=sqrt(fneg2*alog(w(2)))
      angl=twopi*w(1)
      GNNORM=dist*cos(angl)

      return
      end
