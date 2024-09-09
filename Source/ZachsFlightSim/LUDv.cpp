/*==========================================================================
   LUD.c

   This source file contains functions for decomposing and solving a linear
   system of equations, i.e., Ax = b, using the LU decomposition method.

   Functions used:
      LUDecomp - this functions performs LU decomposition of a matrix
      LUSolve  - this function solves the LU decomposed system
===========================================================================*/
/*--------------------------- Include Files -------------------------------*/
#include "LUDv.h"
#include <cmath>
#define ZERO 1.0e-12
using namespace std;

/*------------------------- Function Prototypes -----------------------------*/
double LUDecomp( double** a, int n, double* pvt);
void LUSolve( double** a, int n, double* b,double* pvt);


/*===========================================================================
   LUDecomp

   Description:
      Given an n x n matrix A[0..n-1][0..n-1], this function replaces
      it by the LU decomposition of itself.  This method uses
      Gaussian Elimination techniques to convert matrix A into a
      combination of an upper triangular matrix and a permuted version
      of a lower triangular matrix.  This function is used in
      combination with function LUSolve() to solve systems of equations
      or to invert a matrix.  The determinant of A can also be obtained
      on output from LUDecomp() as follows:
         det| A |  = pvt[n-1] * A[0][0] * A[1][1] * ... * A[n-1][n-1]

   Inputs:
      A - the two-dimensional matrix to be decomposed
      n - the order of matrix A (e.g., the size of A is n x n)

   Outputs:
      pvt - the pivot vector--indicates which rows were swapped
             during decomposition
             pvt[k] = the index of the k-th pivot row
             pvt[n-1] = +1 or -1 depending on whether the number of
                         row interchanges was even or odd, respectively

   Return value:
      The determinant of A is returned.

   Programmer: Jeffrey A. Talbert (based on the FORTRAN version) 01 Oct 1989
   Revised by: D. Snyder (c++ vector class version) 18 Feb 2005
---------------------------------------------------------------------------*/
//double LUDecomp( vector< vector<double> >& a, int n, vector<int>& pvt)

double LUDecomp( double** a, int n, double* pvt)
{
   double anorm, t, detA;
   int nm1, i, j, k, kp1, m;

   pvt[n-1] = 1;

   /* begin processing */
   nm1 = n-1;
   anorm = 0.0;
   for (j=0; j<n; j++) {
      t = 0.0;
      for(i=0; i<n; i++) {
         t += fabs(a[i][j]);
      }
      if(t > anorm) anorm = t;
   }

   /* perform Gaussian Elimination with partial pivoting */
   for(k = 0; k < nm1; k++) {
      kp1 = k+1;

      /* find the pivot element */
      m = k;
      for(i = kp1; i < n; i++) {
         if(fabs(a[i][k]) > fabs(a[m][k])) m = i;
      }
      pvt[k] = m;
      if(m != k) pvt[ n-1 ] = -pvt[ n-1 ];
      t = a[m][k];
      a[m][k] = a[k][k];
      a[k][k] = t;

      /* skip the rest if the pivot is zero */
      if(t != 0.0) {
         /* compute the multipliers */
         for(i = kp1; i < n; i++)
            a[i][k] = -a[i][k]/t;

         /* interchange and eliminate by columns */
         for(j = kp1; j < n; j++) {
            t = a[m][j];
            a[m][j] = a[k][j];
            a[k][j] = t;
            if(t != 0.0) {
               for(i = kp1; i < n; i++)
                  a[i][j] += a[i][k]*t;
            }
         }
      }
   }

   /* Compute and return the determinant of A */
   detA = pvt[n-1];
   for(i = 0; i < n; i++)
      detA *= a[i][i];
   return detA;
}

/*--------------------------------------------------------------------------
   LUSolve

   Description:
      This function solves a set of "n" linear equations in the form
      Ax = b.  The matrix A is input as the LU decomposition of itself,
      determined by the function LUDecomp().  Also, the vector pvt(),
      which is returned from LUDecomp(), as the permutation vector.  The
      vector b is the right-hand side vector on input, and returns with
      the solution vector x on output.  A, n, and pvt are not modified
      by this routine and can be left in place for successive calls with
      different right-hand sides b.

   Inputs:
      A - the LU decomposition of the original matrix A
      n - the order of matrix A (e.g., the size of A is n x n)
      b - the right-hand side vector for the eq'n Ax = b
      pvt - the pivot vector--indicates which rows were swapped
             during decomposition
             pvt[k] = the index of the k-th pivot row
             pvt[n-1] = +1 or -1 depending on whether the number of
                         row interchanges was even or odd, respectively

   Modifies:
      b - the solution vector for the equation Ax = b

   Programmer: Jeffrey A. Talbert (based on the FORTRAN version) 01 Oct 1989
   Revised by:
---------------------------------------------------------------------------*/
void LUSolve( double** a, int n, double* b,double* pvt)

//void LUSolve( vector< vector<double> >& a, int n, vector<double>& b,
 //             vector<int>& pvt)
{
   double t;
   int kb, nm1, kp1, i, k, m;

   /* check for a 1-by-1 matrix */
   if(n == 1) {
      b[0] = b[0]/a[0][0];
      return;
   }

   /* forward elimination */
   nm1 = n-1;
   for(k = 0; k < nm1; k++) {
      kp1 = k+1;
      m = pvt[k];
      t = b[m];
      b[m] = b[k];
      b[k] = t;
      for(i = kp1; i < n; i++)
         b[i] = b[i]+a[i][k]*t;
   }

   /* backward substitution */
   for(kb = 1; kb < n; kb++) {
      k = n-kb;
      b[k] = b[k]/a[k][k];
      t = -b[k];
      for(i = 0; i < k; i++)
         b[i] = b[i]+a[i][k]*t;
   }
   b[0] = b[0]/a[0][0];
   return;
}