/** @file print_matrices.c 
 * Guido Walter Pettinari 23/11/2011
 *
 * Print to screen the matrices needed to compute the multipole decomposition
 * of Boltzmann equation.  These are the Chi matrices.
 *
 */
 
#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions (1st-order) */
  struct perturbs2 pt2;       /* for source functions (2nd-order) */  
  struct threej tj;           /* for 3j symbols (only needed at 2nd-order) */    
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */


  // Print arguments (debug)
  // int jj=0;
  // printf ("argc = %d\n", argc);
  // for (jj=0; jj<argc; ++jj)
  //   printf("argv[%d] = %s\n", jj, argv[jj]);


  // Parse arguments.  First two arguments must be .ini and .pre files
  if (argc != 5) {
    printf("usage:    %s <params.ini> <params.pre> <l> <m>\n", argv[0]);
    return _FAILURE_;
  }
  int l = atoi(argv[3]);
  int m = atoi(argv[4]);
  if (abs(m)>l) {
    printf("Condition abs(m)>l not respected.\n");
    return _FAILURE_;
  }

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&pt2,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (threej_init(&pr,&tj) == _FAILURE_) {
    printf("\n\nError in threej_init \n =>%s\n",tj.error_message);
    return _FAILURE_;
  }
  
  // The array index corresponding to 'm' is just m-MIN(m), that is m+l
  int index_m = m+l;
  
  // Print some info
  printf("# l = %d,  m = %d, index_m = %d\n", l, m, index_m);  
  
  int i=0, j=0, k=0, s=0;  

  if (l==1) {
    for(i=0; i<tj.xi_size; ++i)
      printf("%g %+gi\n", creal(tj.xi[index_m][i]), cimag(tj.xi[index_m][i]));
  }

  if (l==2) {
    for(i=0; i<3; ++i) {

      // Print a line
      for(j=0; j<3; ++j) {
        complex double element = tj.chi_ij[index_m][3*i + j];
        printf("%g %+gi  |", creal(element), cimag(element));
      }
      printf("\n");    

    }
  }

  if (l==3) {
    for(i=0; i<3; ++i) {

      // Separation between matrices
      printf("i = %d\n", i);
      
      for(j=0; j<3; ++j) {
      
        // Print a line
        printf("      ");

        for(k=0; k<3; ++k) {
          complex double element = tj.chi_ijk[index_m][3*3*i+3*j+k];
          printf(" %g %+gi |", creal(element), cimag(element));
          // printf("%+10g + %+10gi  |", creal(tj.chi_ijk[index_m][3*3*i+3*j+k]), cimag(tj.chi_ijk[index_m][3*3*i+3*j+k]));          
        }
        printf("\n");    

      } // j
    } // i
  } // if

  if (l==4) {
    for(i=0; i<3; ++i) {

      // Separation between matrices of matrices
      printf("i = %d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", i);
      
      for(j=0; j<3; ++j) {

        // Separation between matrices
        printf("       j = %d\n", j);
      
        for(k=0; k<3; ++k) {

        // Print a line
        printf("              ");

          for (s=0; s<3; ++s) {
            complex double element = tj.chi_ijks[index_m][3*3*3*i+3*3*j+3*k+s];
            printf(" %g %+gi |", creal(element), cimag(element));
          }
          printf("\n");    
          
        } // s
      } // k
    } // j
  } // i


  return _SUCCESS_;

}
