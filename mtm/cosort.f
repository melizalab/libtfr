      Subroutine COSORT (a,NM,N,b)
c
c     Sort the array A in descending order by straight insertion.
c     Sort the corresponding vectors in array B.
c
      real a(NM),b(NM,1)

      do i = 2,N
         do j = i,2,-1
            if (a(j) .gt. a(j-1)) then
               call SWAP(a(j),a(j-1))
               call BIGSWAP(b(1,j),b(1,j-1),n)
            else
               go to 2
            endif
         end do
 2       continue
      end do

      return
      end


      Subroutine INDEXSORT (a,index,n)
c
c     Create an array INDEX of integers from 1 to n, 
c     sorted such that a(index(j)) <= a(index(j+1))
c     Shell sort is used here.
c
      integer n, index(n), d, flag
      real a(n)

      do i = 1,N
         index(i) = i
      enddo
      d = n
 2    d = (d+1)/2
 3    flag = 0
      do i = 1, n-d
         if (a(index(i)) .gt. a(index(i+d))) then
            call INTSWAP(index(i),index(i+d))
            flag = 1
         end if
      end do
      if (flag .eq. 1 .or. d .gt. 1) goto 2
      return
      end



      Subroutine BIGSWAP(a1,a2,N)
      real a1(N),a2(N)
      do i = 1,N
         call SWAP(a1(i),a2(i))
      end do
      return
      end


      Subroutine SWAP(x,y)
      real x,y
      z = x
      x = y
      y = z
      return
      end

      Subroutine INTSWAP(x,y)
      integer x,y
      z = x
      x = y
      y = z
      return
      end

