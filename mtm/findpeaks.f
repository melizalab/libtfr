      subroutine findpeaks(flow,fhigh,isignal,
     $   nf0,df,nbnd,specraw0,specmed0,specresh0,harmonic0,
     $   specback0,ftest,
     $   conf,fconf,
     $   fsignal,confsignal,nsignals)
c
      parameter (maxsignal=40,nlarge=300)
      parameter (nlim=32768)
      real fsignal(maxsignal),confsignal(maxsignal)
      real fsig(nlarge),confsig(nlarge),sigrat(nlarge)
      integer null(nlim)
      integer nbnd,isignal
      real fconf(6),conf(4),sig(3)
      real ratio,ratmax,ratio0,dummy
      real specraw0(nlim),specmed0(nlim),specresh0(nlim),
     $   harmonic0(nlim),specback0(nlim),ftest(nlim)
c      real whiteraw0,whiterob0,rhoraw0,rhorob0,tauraw0,taurob0
c
c     determine central frequency of all signals above the 90%
c     level (or the maxsignal most significant signals)
c
      nsignals = 0
      do i=1,nf0
         null(i)=0
      end do
      sig(1)=99.0
      sig(2)=95.0
      sig(3)=90.0
c
c     loop search in order of decreasing signficance:99,95,90% conf
c
      do k=1,3
       if (isignal.lt.2) then
        thresh = conf(5-k)
       else
        thresh = fconf(5-k)
       endif
       do i=1,nf0-1
         ff = flow+float(i-1)*df
c
         if (isignal.lt.2) then
            ratio0 = specraw0(i)/specback0(i)
         else
            ratio0 = ftest(i)
         endif
         if ((ratio0.gt.thresh).and.(null(i).eq.0)) then
c
           nsignals=nsignals+1
c
c          if max signal capacity exceeded, quit
c
           if (nsignals.gt.nlarge) then
              nsignals=nsignals-1
              goto 888
           endif
c
c          peak identified -- search for "center" of peak
c          within the spectral bandwidth
c
           ratmax = ratio0
           imax = i
           fmax = ff
           do j=max(1,i-nbnd),min(i+nbnd,nf0-1)
              f0 = flow+float(j-1)*df
              if (isignal.lt.2) then
                ratio = specraw0(j)/specback0(j)
              else
                ratio = ftest(j)
              endif
              if (ratio.gt.ratmax) then
                imax = j
                ratmax=ratio
                fmax = f0
              endif
           end do
c
c          center (significance-wise) of peak identified
c
c          if the center of the peak now falls within one
c          bandwidth of a previously recognized peak, we
c          must throw it out
c
           if (null(imax).eq.1) then
              nsignals = nsignals-1
              goto 444
           endif
c
c          otherwise, we happily add to the list
c
           fsig(nsignals)=fmax
           confsig(nsignals)=sig(k)
           sigrat(nsignals)=ratmax
c
c          now block out the bandwidth from future peak identification
c
           do j=max(1,imax-nbnd),min(imax+nbnd,nf0-1)
              null(j)=1
           end do
c
         endif
444      continue
       end do
      end do 
888   continue
c
c     sort signals in terms of significance
c
      do i=1,nsignals
         do j=i+1,nsignals
            if (sigrat(j).gt.sigrat(i)) then
              dummy = sigrat(j)
              sigrat(j)=sigrat(i)
              sigrat(i) = dummy
              dummy = fsig(j)
              fsig(j)=fsig(i)
              fsig(i) = dummy
              dummy =  confsig(j)
              confsig(j)=confsig(i)
              confsig(i) = dummy
            endif
         end do
      end do
c
c     now store only the first "maxsignal" of these
c
      nsignals=min(maxsignal,nsignals)
c
c     sort the remaining peaks by increasing frequency
c
      do i=1,nsignals
         do j=i+1,nsignals
            if (fsig(j).lt.fsig(i)) then
              dummy = fsig(i)
              fsig(i)=fsig(j)
              fsig(j)=dummy
              dummy = confsig(i)
              confsig(i)=confsig(j)
              confsig(j)=dummy
            endif
         end do
      end do
      do i=1,nsignals         
        fsignal(i)=fsig(i)
        confsignal(i)=confsig(i)
      end do
      return
      end
