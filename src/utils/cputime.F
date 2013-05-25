      subroutine cputime(time)
      double precision time
c     *********
c
c     Subroutine cputime
c
c     This subroutine is used to determine cpu time.
c
c     The subroutine statement is
c
c       subroutine cputime(time)
c
c     where
c
c       time is a double precision variable.
c         On input time need not be specified.
c         On output time specifies the cputime time.
c
c     MINPACK-2 Project. March 2000.
c     Argonne National Laboratory.
c
c     **********
#if defined rs6000
      integer mclock

      time = dble(mclock())/100.0d0

#else
      real tarray(2)
      real etime

      time = etime(tarray)
      time = dble(tarray(1))
#endif

      end
