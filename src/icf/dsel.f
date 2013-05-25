      subroutine dsel(n,x,k)
      integer n, k
      double precision x(n)
c     **********
c
c     Subroutine dsel
c
c     Given an array x of length n, this subroutine permutes
c     the elements of the array x so that 
c    
c       x(i) <= x(k),  1 <= i <= k,
c       x(k) <= x(i),  k <= i <= n.
c
c     In other words, the smallest k elements of x are x(i), i = 1,...,k, 
c     and x(k) is the kth smallest element.
c     
c     The subroutine statement is
c
c       subroutine dsel(n,x,k)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of keys.
c         On exit n is unchanged.
c
c       x is a double precision array of length n.
c         On entry x is the array to be permuted.
c         On exit x is permuted so that the smallest k elements of x are
c            x(i), i = 1,...,k, and x(k) is the kth smallest element.
c
c       k is an integer.
c         On entry k specifes the kth largest element.
c         On exit k is unchanged.
c
c     MINPACK-2 Project. March 1998.
c     Argonne National Laboratory.
c     William D. Kastak, Chih-Jen Lin, and Jorge J. More'.
c
c     **********
      integer i, l, lc, lp, m, p, p1, p2, p3, u
      double precision swap

      if (n .le. 1 .or. k .le. 0 .or. k .gt. n) return

      u = n
      l = 1
      lc = n
      lp = 2*n

c     Start of iteration loop.

      do while (l .lt. u)

c        Choose the partition as the median of the elements in
c        positions l+s*(u-l) for s = 0, 0.25, 0.5, 0.75, 1.
c        Move the partition element into position l.

         p1 = (u+3*l)/4
         p2 = (u+l)/2
         p3 = (3*u+l)/4

c        Order the elements in positions l and p1.

         if (x(l) .gt. x(p1)) then
            swap = x(l)
            x(l) = x(p1)
            x(p1) = swap
            end if

c        Order the elements in positions p2 and p3.

         if (x(p2) .gt. x(p3)) then
            swap = x(p2)
            x(p2) = x(p3)
            x(p3) = swap
            end if

c        Swap the larger of the elements in positions p1
c        and p3, with the element in position u, and reorder
c        the first two pairs of elements as necessary.

         if (x(p3) .gt. x(p1)) then
            swap = x(p3)
            x(p3) = x(u)
            x(u) = swap
            if (x(p2) .gt. x(p3)) then
               swap = x(p2)
               x(p2) = x(p3)
               x(p3) = swap
               end if
         else
            swap = x(p1)
            x(p1) = x(u)
            x(u) = swap
            if (x(l) .gt. x(p1)) then
               swap = x(l)
               x(l) = x(p1)
               x(p1) = swap
               end if
            end if

c        We have now permuted the array x so that
c
c          x(l) <= x(p1), x(p2) <= x(p3), mxx(x(p1),x(p3)) <= x(u).
c
c        Find the third largest element of the four remaining
c        elements (the median), and place in position l.

         if (x(p1) .gt. x(p3)) then
            if (x(l) .le. x(p3)) then
               swap = x(l)
               x(l) = x(p3)
               x(p3) = swap
               end if
         else
            if (x(p2) .le. x(p1)) then
               swap = x(l)
               x(l) = x(p1)
               x(p1) = swap
            else
               swap = x(l)
               x(l) = x(p2)
               x(p2) = swap
               end if
            end if

c        Partition the array about the element in position l.

         m = l
         do i = l+1, u
            if (x(i) .lt. x(l)) then
               m = m + 1
               swap = x(m)
               x(m) = x(i)
               x(i) = swap
               end if
         end do

c        Move the partition element into position m.

         swap = x(l)
         x(l) = x(m)
         x(m) = swap

c        Adjust the values of l and u.

         if (k .ge. m) l = m + 1
         if (k .le. m) u = m - 1

c        Check for multiple medians if the length of the subarray
c        has not decreased by 1/3 after two consecutive iterations.

         if (3*(u-l) .gt. 2*lp .and. k .gt. m) then

c           Partition the remaining elements into those elements
c           equal to x(m), and those greater than x(m). Adjust
c           the values of l and u.

            p = m
            do i = m+1, u
               if (x(i) .eq. x(m)) then
                  p = p + 1
                  swap = x(p)
                  x(p) = x(i)
                  x(i) = swap
                  end if
            end do
            l = p + 1
            if (k .le. p) u = p - 1
            end if

c        Update the length indicators for the subarray.

         lp = lc
         lc = u-l

      end do

      return

      end
