      subroutine dicfs(n,nnz,a,adiag,acol_ptr,arow_ind,
     +                 l,d,lcol_ptr,lrow_ind,
     +                 p,alpha,droptol,minpiv,iwa,wa1,wa2)
      integer n, nnz, p
      integer acol_ptr(n+1), arow_ind(nnz)
      integer lcol_ptr(n+1), lrow_ind(nnz+n*p)
      integer iwa(3*n)
      double precision alpha, minpiv, droptol
      double precision wa1(n), wa2(n)
      double precision a(nnz), adiag(n), l(nnz+n*p), d(n)
c     *********
c
c     Subroutine dicfs
c
c     Given a symmetric matrix A in compressed column storage, this
c     subroutine computes an incomplete LDL factorization of A + alpha*D
c     where alpha > 0 is a shift and D is a diagonal matrix whose
c     entries are the 2-norm of the columns of A multiplied by the sign
c     of the diagonal element of A.
c     dominique.orban@gerad.ca based on Lin and More's original ICFS.
c
c     The subroutine statement is
c
c       subroutine dicfs(n,nnz,a,adiag,acol_ptr,arow_ind,
c                        l,d,lcol_ptr,lrow_ind,
c                        p,alpha,iwa,wa1,wa2)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of A.
c         On exit n is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz is the number of nonzeros in the strict lower
c            triangular part of A.
c         On exit nnz is unchanged.
c
c       a is a double precision array of dimension nnz.
c         On entry a must contain the strict lower triangular part
c            of A in compressed column storage.
c         On exit a is unchanged.
c
c       adiag is a double precision array of dimension n.
c         On entry adiag must contain the diagonal elements of A.
c         On exit adiag is unchanged.
c
c       acol_ptr is an integer array of dimension n + 1.
c         On entry acol_ptr must contain pointers to the columns of A.
c            The nonzeros in column j of A must be in positions
c            acol_ptr(j), ... , acol_ptr(j+1) - 1.
c         On exit acol_ptr is unchanged.
c
c       arow_ind is an integer array of dimension nnz.
c         On entry arow_ind must contain row indices for the strict
c            lower triangular part of A in compressed column storage.
c         On exit arow_ind is unchanged.
c
c       l is a double precision array of dimension nnz+n*p.
c         On entry l need not be specified.
c         On exit l contains the strict lower triangular part
c            of L in compressed column storage.
c
c       d is a double precision array of dimension n.
c         On entry d need not be specified.
c         On exit d contains the diagonal elements of D.
c
c       lcol_ptr is an integer array of dimension n + 1.
c         On entry lcol_ptr need not be specified.
c         On exit lcol_ptr contains pointers to the columns of L.
c            The nonzeros in column j of L are in the
c            lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.
c
c       lrow_ind is an integer array of dimension nnz+n*p.
c         On entry lrow_ind need not be specified.
c         On exit lrow_ind contains row indices for the strict lower
c            triangular part of L in compressed column storage.
c
c       p is an integer variable.
c         On entry p specifes the amount of memory available for the
c            incomplete Cholesky factorization.
c         On exit p is unchanged.
c
c       alpha is a double precision variable.
c         On entry alpha is the initial guess of the shift.
c         On exit alpha is final shift
c
c       minpiv is a double precision variable.
c         On entry minpiv is the minimum allowed pivot magnitude.
c           If a pivot of magnitude smaller than minpiv is generated,
c           the shift will be increased and a new factorization will
c           be attempted.
c         On exit minpiv is unchanged.
c
c       droptol is a double precision variable.
c         On entry droptol is the drop tolerance. Any entry in L
c           with magnitude inferior to droptol will be discarded.
c         On exit droptol is unchanged.
c
c       iwa is an integer work array of dimesnion 3*n.
c
c       wa1 is a double precision work array of dimension n.
c
c       wa2 is a double precision work array of dimension n.
c
c     Subprograms called
c
c       MINPACK-2  ......  dicf
c
c     MINPACK-2 Project. October 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      integer nbmax
      parameter(nbmax=1)
      double precision alpham, nbfactor
      parameter(alpham=1.0d-3,nbfactor=512)
      double precision zero, one, two
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)

      integer i, info, j, k, nb
      double precision alphas

      external dicf

c     Compute the l2 norms of the columns of A
c     and largest element in absolute value.

      do i = 1, n
         wa1(i) = adiag(i)**2
      end do
      do j = 1, n
         do i = acol_ptr(j), acol_ptr(j+1)-1
            k = arow_ind(i)
            wa1(j) = wa1(j) + a(i)**2
            wa1(k) = wa1(k) + a(i)**2
         end do
      end do
      do j = 1, n
         wa1(j) = sqrt(wa1(j))
      end do

c     Compute the scaling matrix S.

      do i = 1, n
         if (wa1(i) .gt. zero) then
            wa2(i) = one/sqrt(wa1(i))
         else
            wa2(i) = one
         endif
      end do

c     Compute the initial shift.

      alphas = alpham
      if (alpha .gt. zero) alpha = max(alpha, alphas)
      do i = 1, n
         if (adiag(i) .eq. zero) then
            alpha = max(alpha, alphas)
         end if
      end do

c     Check values of minpiv and droptol.

      if (minpiv .lt. zero) minpiv = zero
      if (droptol .lt. zero) droptol = zero

c     Search for an acceptable shift. During the search we decrease
c     the lower bound alphas until we determine a lower bound that
c     is not acceptable. We then increase the shift.
c     The lower bound is decreased by nbfactor at most nbmax times.

      nb = 1
      do while (1 .eq. 1)

c     Copy the sparsity structure of A into L.

         do i = 1, n+1
            lcol_ptr(i) = acol_ptr(i)
         end do
         do i = 1, nnz
            lrow_ind(i) = arow_ind(i)
         end do

c     Scale A and store in the lower triangular matrix L.

         do j = 1, n
            d(j) = adiag(j)*(wa2(j)**2)
            if (adiag(j) .gt. zero) then
               d(j) = d(j) + alpha
            else
               d(j) = d(j) - alpha
            endif
         end do
         do j = 1, n
            do i = acol_ptr(j), acol_ptr(j+1)-1
               l(i) = a(i)*wa2(j)*wa2(arow_ind(i))
            end do
         end do

c     Attempt an incomplete factorization.

         call dicf(n,nnz,l,d,lcol_ptr,lrow_ind,p,droptol,minpiv,info,
     +             iwa(1),iwa(n+1),iwa(2*n+1),wa1)

c        If the factorization exists, then test for termination.
c        Otherwise increment the shift.

         if (info .ge. 0) then

c           If the shift is at the lower bound, reduce the shift.
c           Otherwise undo the scaling of L and exit.

            if (alpha .eq. alphas .and. nb .lt. nbmax) then
               alphas = alphas/nbfactor
               alpha = alphas
               nb = nb + 1
            else

c     Undo the scaling.

               do i = 1, n
                  d(i) = d(i) / (wa2(i)**2)
               end do
               do j = 1, n
                  do i = lcol_ptr(j), lcol_ptr(j+1)-1
                     l(i) = l(i) * wa2(j) / wa2(lrow_ind(i))
                  end do
               end do
               return
            endif

         else
            alpha = max(two*alpha,alphas)
         end if
      end do

      end
