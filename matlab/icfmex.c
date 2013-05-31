/*
Calling syntax :

   [L, Ldiag] = icfmex(lA, Adiag, p);

*/
#if ARCH == aix-4
#define dicfs_ dicfs
#endif

#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

#define ICFS_MAX(a,b) ((a) > (b) ? (a) : (b))
#define ICFS_MIN(a,b) ((a) < (b) ? (a) : (b))

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
  void print_int_array(const char *name, int *array, int len,
                       const char *fmt);
  void print_double_array(const char *name, double *array, int len,
                          const char *fmt);

   /* Main entry point */
   void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
   {
      int n = 0, p = 0;
      double *a = NULL, *l = NULL;
      double *adiag = NULL, *ldiag = NULL;
      mwIndex *acol_ptr = NULL, *arow_ind = NULL;
      int *acol_ptr_int = NULL, *arow_ind_int = NULL;
      mwIndex *lcol_ptr = NULL, *lrow_ind = NULL;
      int *lcol_ptr_int = NULL, *lrow_ind_int = NULL;

      /* Working arrays */

      double *w1 = NULL, *w2 = NULL;
      int *iw = NULL;

      /* Local variables */

      int i = 0, nnz = 0;
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
      if (n >= INT_MAX)
        mexErrMsgTxt("Indexing exceeds INT limits.");

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
      if (! (arow_ind_int = (int *)mxCalloc(ICFS_MAX(nnz,1), sizeof(int))))
         mexErrMsgTxt("Not enough memory for arow_ind");
      if (! (acol_ptr_int = (int *)mxCalloc(n+1, sizeof(int))))
         mexErrMsgTxt("Not enough memory for acol_ptr");

      for (i = 0; i < nnz; i++) {
         if (arow_ind[i] >= INT_MAX)
            mexErrMsgTxt("Indexing exceeds INT limits.");
         arow_ind_int[i] = (int)arow_ind[i] + 1;
      }

      for (i = 0; i <= n; i++)
         acol_ptr_int[i] = (int)acol_ptr[i] + 1;

  #ifdef MXDEBUG
      mexPrintf("Calling icfmex with n=%d, nnz=%d, p=%d\n", n, nnz, p);
      mexPrintf("%d input args and %d output args\n", nrhs, nlhs);
      mexPrintf("After conversion to Fortran indexing:\n");
      print_int_array("arow_ind", arow_ind_int, nnz, "%d ");
      print_int_array("acol_ptr", acol_ptr_int, n+1, "%d ");
      print_double_array("adiag", adiag, n, "%8.1e ");
      print_double_array("a", a, nnz, "%8.1e ");
   #endif

      /* Left-hand side: output args */

      plhs[0] = mxCreateSparse((mwSize)n, (mwSize)n,
                               (mwSize)ICFS_MAX(nnz+n*p,1), mxREAL);
      l = (double *)mxGetPr(plhs[0]);
      lrow_ind = mxGetIr(plhs[0]);
      lcol_ptr = mxGetJc(plhs[0]);

      /* Allocate int* arrays to pass to Fortran. These will be copied
         into lrow_ind and lcol_ptr later. */
      if (! (lrow_ind_int = (int *)mxCalloc(ICFS_MAX(nnz+n*p,1), sizeof(int))))
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
      print_int_array("lrow_ind", lrow_ind_int, nnzl, "%d ");
      print_int_array("lcol_ptr", lcol_ptr_int, n+1, "%d ");
      print_double_array("ldiag", ldiag, n, "%8.1e ");
      print_double_array("l", l, nnzl, "%8.1e ");
   #endif

      /* Sort each segment of lrow_ind in ascending order and sync l. */
      for (i = 0; i < n; i++) {
   #ifdef MXDEBUG
        mexPrintf("Sorting lrow_ind[%d : %d]\n", lcol_ptr_int[i], lcol_ptr_int[i+1]-1);
   #endif
        quicksort_icfs(lrow_ind_int, l, lcol_ptr_int[i]-1, lcol_ptr_int[i+1]-2);
      }

      /* Copy lcol_ptr_int and lrow_ind_int into mwIndex* arrays.
         Account for C indexing. */
      if (nnzl < nnz + n*p) {
        if (! (mxRealloc((void *)lrow_ind, nnzl * sizeof(mwIndex))))
          mexErrMsgTxt("Insufficient memory to reallocate lrow_ind");
        if (! (mxRealloc((void *)l, nnzl * sizeof(double))))
          mexErrMsgTxt("Insufficient memory to reallocate l");
      }
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
      quicksort_icfs(numbers, values, low, current-1);    if (up > current)
      quicksort_icfs(numbers, values, current+1, up);
   }

   /* Debug helper */

   #define NELMAX 3    /* Number of bordering array elements to print */

   void print_int_array(const char *name, int *array, int len,
                        const char *fmt) {

      int i;
      int nel = ICFS_MIN(len, NELMAX);

      mexPrintf(name);
      mexPrintf(" = [ ");
      for (i = 0; i < nel; i++) mexPrintf(fmt, array[i]);

      if (len > 2*nel) mexPrintf(" ... ");

      if (len > nel)
        for (i = 0; i < ICFS_MIN(len-nel, NELMAX); i++)
          mexPrintf(fmt, array[len-1-i]);

      mexPrintf("]\n");

   }

   void print_double_array(const char *name, double *array, int len,
                           const char *fmt) {

      int i;
      int nel = ICFS_MIN(len, NELMAX);

      mexPrintf(name);
      mexPrintf(" = [ ");
      for (i = 0; i < nel; i++) mexPrintf(fmt, array[i]);

      if (len > 2*nel) mexPrintf(" ... ");

      if (len > nel)
        for (i = 0; i < ICFS_MIN(len-nel, NELMAX); i++)
          mexPrintf(fmt, array[len-1-i]);

      mexPrintf("]\n");

   }

/* Safeguard against C++ symbol mangling */
#ifdef __cplusplus
}
#endif
