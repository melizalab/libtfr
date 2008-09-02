c
c $Id: fieldio.F,v 1.8 1997/02/14 17:49:57 weibel Exp $
c


      Subroutine Get1dField ( file, ncols, npoints, y, x )
c
c get1dfield - get a column of formatted data from a file.
c        return the number of columns, number of points, and the data.
c        if the file has only 1 column, return the data in y, the index in x.
c        otherwise, return column 1 in x and column 2 in y
c
c  file  -  name of the input file to open.  
c           If file = 'stdin' or blank, read from the standard input
c           If file does not exist, return ncols, npoints = 0.
c
      character file*(*), record*1024, c
      integer ncols, npoints, luin
      real x(1), y(1)

c initialize the data counter, and open the input file
      npoints = 1
      ncols = 0
      if ( file .eq. 'stdin' .or. file .eq. ' ') then
         luin = 5
      else
         luin = lu_handy()
         open (unit=luin, file=file, status='OLD', err=9)
      end if

c skip comments     
 1    read (luin, '(a)', end=9) record
      c = record(1:1)
      if ((c.ge.'A'.and.c.le.'z') .or. c.eq.'%' .or. c.eq.'#') go to 1
     
c The current record is the first data point.  Count the number of words in it.
c   Then read the rest of the file
      ncols = nwords(record,256)

      if (ncols .eq. 1) then
         read (record, *) y(npoints)
 24      npoints = npoints+1
         read (luin, *, end=9) y(npoints)
         x(npoints) = npoints-1
         go to 24
      endif
      if (ncols .gt. 1) then
         read (record, *) x(npoints), y(npoints)
 26      npoints = npoints+1
         read (luin, *, end=9) x(npoints), y(npoints)
         go to 26
      endif

 9    npoints = npoints-1
      close(luin)
      return
      end



c
c nwords  -  count the number of words in a string.
c            a word is defined here as any string that contains no spaces, tabs
c            or commas
c
      integer function nwords ( string, len)
      logical inword, comma
      character*(*) string
      character c, tab
      tab = '	'

      inword = .false.
      comma = .false.
      nwords = 0
      do i=1,len
         c=(string(i:i+1))
         if (c .eq. ' ' .or. c .eq. tab .or. c .eq. ',') then
           inword = .false.
           if ( c .eq. ',' ) then
             if (comma) nwords = nwords + 1
             comma = .true.
           endif
         else if (.not. inword) then
           inword = .true.
           comma = .false.
           nwords = nwords + 1
         end if
      end do
      return
      end


c
c find a logical unit that is not being used
c
      integer function lu_handy()
      logical in_use
      do lu=1,99
         lu_handy = lu
         inquire (lu, opened=in_use)
         if (.not.in_use) return
      end do
      lu_handy = -1
      end

#ifdef AIX
c
c A substitute version of exit for IBM RS6000 systems running AIX
c
      subroutine exit(status)
      integer status
      call exit_(status)
      end

      integer function lnblnk(string)
      character*(*) string
      lnblnk = len(trim(string))
      return
      end
#endif
