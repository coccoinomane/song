#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#define ALTERNATING_SIGN(m) ((m)%2 == 0 ? 1 : -1)

/*
   Recursive definition of determinate using expansion by minors.
   Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
*/
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         // det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         det += ALTERNATING_SIGN(j1+2) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}


/*
   Find the cofactor matrix of a square matrix
   Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         // b[i][j] = pow(-1.0,i+j+2.0) * det;
         b[i][j] = ALTERNATING_SIGN(i+j+2) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}


/*
   Transpose of a square matrix, do it in place
   Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}


/* 
   Inverse of a square matrix
*/
void InverseMatrix(double **in,int n,double **out)
{

  int i,j;
  
  /* Compute determinant */
  double det = Determinant(in, n);

  /* Compute cofactor matrix */
  double ** cofactor = malloc(n*sizeof(double *));
  for (i=0;i<n;i++)
    cofactor[i] = malloc(n*sizeof(double));
  CoFactor(in, n, cofactor);

  /* Compute inverse matrix (note that we take the transpose of cofactor) */
  for (i=0;i<n;i++)
     for (j=0;j<n;j++)
        out[i][j] = cofactor[j][i]/det;
  
  for (i=0;i<n;i++)
     free(cofactor[i]);
  free(cofactor);

}



int main (int argc, char const *argv[])
{

  int n = 3;
  int i,j;

  double ** M = malloc(n*sizeof(double *));
  for (i=0;i<n;i++)
    M[i] = malloc(n*sizeof(double));

  double ** M_INV = malloc(n*sizeof(double *));
  for (i=0;i<n;i++)
    M_INV[i] = malloc(n*sizeof(double));

  M[0][0] = +1.; M[0][1] = +2.; M[0][2] = +3.;
  M[1][0] = +2.; M[1][1] = +0.; M[1][2] = +2.;
  M[2][0] = +3.; M[2][1] = -1.; M[2][2] = -4.;

  printf("~~~ Original matrix ~~~\n");
  for (i=0;i<n;i++) {
     for (j=0;j<n;j++) {
       printf ("%8g ", M[i][j]);
     }
     printf("\n");
  }
  printf("\n");

  double det = Determinant(M, n);
  InverseMatrix(M, n, M_INV);

  printf("~~~ Determinant ~~~\n");  
  printf("det = %g\n", det);
  printf("\n");

  printf("~~~ Inverted matrix ~~~\n");
  for (i=0;i<n;i++) {
     for (j=0;j<n;j++) {
       printf ("%8g ", M_INV[i][j]);
     }
     printf("\n");
  }
  printf("\n");
  
  return 0;
}














