/*
The calling syntax :

   [L, Ldiag] = icfmex(n, lA, Adiag, p) ;

*/
#if ARCH == aix-4
#define dicfs_ dicfs
#endif

#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

/* For versions of the Matlab API prior to 7.3 */
#if (MX_API_VER < 0x07030000)
typedef int mwSize;
typedef int mwIndex;
#endif

/* Prototypes */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void quicksort_icfs(mwIndex numbers[], double values[],
                    mwIndex low, mwIndex up);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int n, p ;
   double *a, *l;
   double *adiag, *ldiag;
   mwIndex *acol_ptr, *arow_ind ;
   mwIndex *lcol_ptr, *lrow_ind ;

   /* Working arrays */

   double *w1, *w2 ;
   int *iw ;

   /* Local variables */

   int i, nnz ;
   double alpha ;

   /* Right hand side */

   n = (int) mxGetScalar(prhs[0]);
   a = (double*) mxGetData(prhs[1]);
   arow_ind = mxGetIr(prhs[1]);
   acol_ptr = mxGetJc(prhs[1]);
   adiag = (double*) mxGetData(prhs[2]) ;
   p = (int) mxGetScalar(prhs[3]) ;

   nnz = (int) acol_ptr[n] ;

#ifdef MXDEBUG
   mexPrintf("Calling icfmex with n=%d, nnz=%d, p=%d\n", n, nnz, p);
   mexPrintf("%d input args and %d output args\n", nrhs, nlhs);
   mexPrintf("arow_ind = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) arow_ind[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("acol_ptr = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) acol_ptr[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("adiag = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) adiag[i]);
   mexPrintf(" ... ]\n");
#endif

   /* Left-hand side */

   if (plhs[0] != NULL) mxDestroyArray(plhs[0]) ;
   plhs[0] = mxCreateSparse((mwSize) n, (mwSize) n, (mwSize) (nnz+n*p), mxREAL) ;

   l =  (double*) mxGetData(plhs[0]);
   lrow_ind = mxGetIr(plhs[0]);
   lcol_ptr = mxGetJc(plhs[0]);

   if (plhs[1] != NULL) mxDestroyArray(plhs[1]) ;
   plhs[1] = mxCreateDoubleMatrix((mwSize) n, (mwSize) 1, mxREAL) ;
   ldiag = (double*) mxGetData(plhs[1]) ;

   /* Allocate working arrays */

   if (( (w1 = (double*) mxCalloc(n, sizeof(double))) == (double *)NULL ) ||
       ( (w2 = (double*) mxCalloc(n, sizeof(double))) == (double *)NULL ) ||
       ( (iw = (int*) mxCalloc(3*n, sizeof(int))) == (int *)NULL )) {
     mexErrMsgTxt("Not enough memory\n") ;
   }

   /* Change a's i j into fortran index */

   for (i = 0; i < nnz; i++) arow_ind[i] = arow_ind[i] + 1 ;
   for (i = 0; i <= n; i++)  acol_ptr[i] = acol_ptr[i] + 1 ;

#ifdef MXDEBUG
   mexPrintf("After conversion to Fortran indexing:\n");
   mexPrintf("arow_ind = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) arow_ind[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("acol_ptr = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) acol_ptr[i]);
   mexPrintf(" ... ]\n");
#endif

   /* Call icf fortran subroutine */

   alpha = 0.0 ;

#ifdef MXDEBUG
   mexPrintf("Calling dicfs\n");
#endif

   dicfs_(&n, &nnz, a, adiag, (int*) acol_ptr, (int*) arow_ind,
	  l, ldiag, (int*) lcol_ptr, (int*) lrow_ind, &p, &alpha,
	  iw, w1, w2) ;

#ifdef MXDEBUG
   mexPrintf("Returning from dicfs with alpha=%7.1e\n", alpha);
   mexPrintf("lrow_ind = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) lrow_ind[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("lcol_ptr = [ ");
   for (i=0; i<5; i++) mexPrintf("%d ", (int) lcol_ptr[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("l = [ ");
   for (i=0; i<5; i++) mexPrintf("%8.1e ", l[i]);
   mexPrintf(" ... ]\n");
   mexPrintf("ldiag = [ ");
   for (i=0; i<5; i++) mexPrintf("%8.1e ", ldiag[i]);
   mexPrintf(" ... ]\n");
#endif

   /* Sort each segment of lrow_ind in ascending order
      and keep l synchronized. */

#ifdef MXDEBUG
   mexPrintf("Sorting L\n");
#endif

   for (i = 0; i < n; i++)
      quicksort_icfs(lrow_ind, l, lcol_ptr[i], lcol_ptr[i+1]-1);

#ifdef MXDEBUG
   mexPrintf("Restoring C indexing\n");
#endif

   /* Change a's i j into C index */
   for (i = 0; i < nnz; i++)
      arow_ind[i] = arow_ind[i] - 1 ;
   for (i = 0; i < lcol_ptr[n]; i++)
      lrow_ind[i] = lrow_ind[i] - 1 ;
   for (i = 0; i <= n; i++) {
      acol_ptr[i] = acol_ptr[i] - 1 ;
      lcol_ptr[i] = lcol_ptr[i] - 1 ;
   }

#ifdef MXDEBUG
   mexPrintf("Freeing memory\n");
#endif

   /* Free memory */

   mxFree((void*) w1) ;
   mxFree((void*) w2) ;
   mxFree((void*) iw) ;
} /* mexFunction */


void quicksort_icfs(mwIndex numbers[], double values[],
                    mwIndex low, mwIndex up) {
 int current, low_current, up_current;
 double dcurrent;

 low_current = low;
 up_current = up;
 current = numbers[low];
 dcurrent = values[low];
 while (low < up)
 {
   while ((numbers[up] >= current) && (low < up))
     up--;
   if (low != up)
   {
     numbers[low] = numbers[up];
     values[low] = values[up];
     low++;
   }
   while ((numbers[low] <= current) && (low < up))
     low++;
   if (low != up)
   {
     numbers[up] = numbers[low];
     values[up] = values[low];
     up--;
   }
 }
 numbers[low] = current;
 values[low] = dcurrent;
 current = low;
 low = low_current;
 up = up_current;
 if (low < current)
   quicksort_icfs(numbers, values, low, current-1);
 if (up > current)
   quicksort_icfs(numbers, values, current+1, up);
}
