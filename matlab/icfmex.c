/*
The calling syntax :

   [L, Ldiag] = icfmex(lA, Adiag, p);

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

/* Safeguard against C++ symbol mangling */
#ifdef __cplusplus
extern "C" {
#endif

   /* Prototypes */
   void mexFunction(int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[]);
   void dicfs_(int *n, int *nnz,
               double *a, double *adiag, int *acol_ptr, int *arow_ind,
               double *l, double *ldiag, int *lcol_ptr, int *lrow_ind,
               int *p, double *alpha,
               int *iw, double *w1, double *w2);
   void quicksort_icfs(int numbers[], double values[], int low, int up);

   /* Main entry point */
   void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
   {
      int n, p;
      double *a, *l;
      double *adiag, *ldiag;
      mwIndex *acol_ptr, *arow_ind;
      int *acol_ptr_int, *arow_ind_int;
      mwIndex *lcol_ptr, *lrow_ind;
      int *lcol_ptr_int, *lrow_ind_int;

      /* Working arrays */

      double *w1, *w2;
      int *iw;

      /* Local variables */

      int i, nnz;
      double alpha = 0.0;

      /* Right hand side: input args */

      /* Check data type of input argument  */
      if (!(mxIsDouble(prhs[0])))
         mexErrMsgIdAndTxt("MATLAB:icfmex:inputNotDouble",
               "Input argument 0 must be of type double.");

      if (!(mxIsSparse(prhs[0])))
        mexErrMsgIdAndTxt("MATLAB:icfmex:inputNotSparse",
                      "Input argument 0 must be sparse.");

      if (mxGetNumberOfDimensions(prhs[0]) != 2)
         mexErrMsgIdAndTxt("MATLAB:icfmex:inputNot2D",
          "Input argument 0 must be two dimensional.");

      n = (int)mxGetN(prhs[0]);

      if (n != (int)mxGetM(prhs[0]))
         mexErrMsgIdAndTxt("MATLAB:icfmex:inputNotSquare",
                       "Input argument 0 must be square.");

      a = mxGetPr(prhs[0]);                  /* Length nnz (0 ... nnz-1) */
      arow_ind = mxGetIr(prhs[0]);           /* Length nnz (0 ... nnz-1) */
      acol_ptr = mxGetJc(prhs[0]);           /* Length n+1 (0 ... n)     */
      nnz = (int)acol_ptr[n];

      if (!(mxIsDouble(prhs[1])))
         mexErrMsgIdAndTxt("MATLAB:icfmex:inputNotDouble",
               "Input argument 1 must be of type double.");

      adiag = (double*)mxGetData(prhs[1]);   /* Length n   (0 ... n-1)   */

      p = (int)mxGetScalar(prhs[2]);

      /* Because mwIndex differs from int, we are required to COPY the
         input matrix into int* arrays. Isn't Matlab the tool?
         Anyways. Account for Fortran indexing. */
      if (! (arow_ind_int = (int *)mxCalloc(nnz, sizeof(int))))
         mexErrMsgTxt("Not enough memory for arow_ind");
      if (! (acol_ptr_int = (int *)mxCalloc(n+1, sizeof(int))))
         mexErrMsgTxt("Not enough memory for acol_ptr");

      for (i = 0; i < nnz; i++) {
         if (arow_ind[i] >= INT_MAX)
            mexErrMsgTxt("Indexing exceeds INT limits.");
         arow_ind_int[i] = (int)arow_ind[i] + 1;
      }

      for (i = 0; i <= n; i++) {
         if (acol_ptr[i] >= INT_MAX)
            mexErrMsgTxt("Indexing exceeds INT limits.");
         acol_ptr_int[i] = (int)acol_ptr[i] + 1;
      }

   #ifdef MXDEBUG
      mexPrintf("Calling icfmex with n=%d, nnz=%d, p=%d\n", n, nnz, p);
      mexPrintf("%d input args and %d output args\n", nrhs, nlhs);
      mexPrintf("arow_ind = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", (int)arow_ind[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", (int)arow_ind[nnz-1-i]);
      mexPrintf("]\n");

      mexPrintf("acol_ptr = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", (int)acol_ptr[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", (int)acol_ptr[n-i]);
      mexPrintf("]\n");

      mexPrintf("adiag = [ ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", adiag[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", adiag[n-1-i]);
      mexPrintf("]\n");

      mexPrintf("a = [ ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", a[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", a[nnz-1-i]);
      mexPrintf("]\n");
   #endif

      /* Left-hand side: output args */

      plhs[0] = mxCreateSparse((mwSize)n, (mwSize)n, (mwSize)(nnz+n*p), mxREAL);
      l = (double *)mxGetPr(plhs[0]);
      lrow_ind = mxGetIr(plhs[0]);
      lcol_ptr = mxGetJc(plhs[0]);

      /* Allocate int* arrays to pass to Fortran. These will be copied
         into lrow_ind and lcol_ptr later. */
      if (! (lrow_ind_int = (int *)mxCalloc(nnz+n*p, sizeof(int))))
         mexErrMsgTxt("Not enough memory for lrow_ind");
      if (! (lcol_ptr_int = (int *)mxCalloc(n+1, sizeof(int))))
         mexErrMsgTxt("Not enough memory for lcol_ptr");

      plhs[1] = mxCreateDoubleMatrix((mwSize)n, (mwSize)1, mxREAL);
      ldiag = (double*)mxGetData(plhs[1]);

      /* Allocate working arrays */
      if (( (w1 = (double*)mxCalloc(n, sizeof(double))) == (double *)NULL ) ||
          ( (w2 = (double*)mxCalloc(n, sizeof(double))) == (double *)NULL ) ||
          ( (iw = (int*)mxCalloc(3*n, sizeof(int))) == (int *)NULL )) {
        mexErrMsgTxt("Not enough memory for temporary arrays");
      }

   #ifdef MXDEBUG
      mexPrintf("After conversion to Fortran indexing:\n");
      mexPrintf("arow_ind = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", arow_ind_int[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", arow_ind_int[nnz-1-i]);
      mexPrintf("]\n");

      mexPrintf("acol_ptr = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", acol_ptr_int[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", acol_ptr_int[n-i]);
      mexPrintf("]\n");
   #endif

      /* Call icf fortran subroutine */

   #ifdef MXDEBUG
      mexPrintf("Calling dicfs\n");
   #endif

      dicfs_(&n, &nnz,
             a, adiag, acol_ptr_int, arow_ind_int,
   	       l, ldiag, lcol_ptr_int, lrow_ind_int,
             &p, &alpha,
   	       iw, w1, w2);

      int nnzl = lcol_ptr_int[n] - 1;  /* -1 accounts for 1-based indexing */

   #ifdef MXDEBUG
      mexPrintf("Returning from dicfs with alpha=%7.1e, nnzl=%d\n", alpha, nnzl);
      mexPrintf("lrow_ind = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", lrow_ind_int[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", lrow_ind_int[nnzl-1-i]);
      mexPrintf("]\n");

      mexPrintf("lcol_ptr = [ ");
      for (i=0; i<3; i++) mexPrintf("%d ", lcol_ptr_int[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%d ", lcol_ptr_int[n-i]);
      mexPrintf("]\n");

      mexPrintf("ldiag = [ ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", ldiag[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", ldiag[n-1-i]);
      mexPrintf("]\n");

      mexPrintf("l = [ ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", l[i]);
      mexPrintf(" ... ");
      for (i=0; i<3; i++) mexPrintf("%8.1e ", l[nnzl-1-i]);
      mexPrintf("]\n");
   #endif

      /* Sort each segment of lrow_ind in ascending order and sync l. */
      for (i = 0; i < n; i++) {
   #ifdef MXDEBUG
      mexPrintf("Sorting lrow_ind[%d : %d]\n", lcol_ptr_int[i], lcol_ptr_int[i+1]-1);
   #endif
         quicksort_icfs(lrow_ind_int, l, lcol_ptr_int[i], lcol_ptr_int[i+1]-1);
      }

      /* Copy lcol_ptr_int and lrow_ind_int into mwIndex* arrays.
         Account for C indexing. */
      for (i = 0; i < nnzl; i++)
         lrow_ind[i] = (mwIndex)(lrow_ind_int[i] - 1);

      for (i = 0; i <= n; i++)
         lcol_ptr[i] = (mwIndex)(lcol_ptr_int[i] - 1);

   #ifdef MXDEBUG
      mexPrintf("Freeing memory\n");
   #endif

      /* Free memory */
      mxFree(arow_ind_int);
      mxFree(acol_ptr_int);
      mxFree(lrow_ind_int);
      mxFree(lcol_ptr_int);
      mxFree(w1);
      mxFree(w2);
      mxFree(iw);
   } /* mexFunction */


   void quicksort_icfs(int numbers[], double values[],
                       int low, int up) {
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

/* Safeguard against C++ symbol mangling */
#ifdef __cplusplus
}
#endif
