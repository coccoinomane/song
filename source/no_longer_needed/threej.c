/** @file bessel.c Documented Bessel module.
 *
 * Guido W. Pettinari, 2.11.2011
 *
 * This module loads the Wigner 3j symbols needed to compute the Boltzmann
 * hierarchy at second order (either read from file or computed from scratch).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# threej_init() at the beginning (anytime after input_init() and before transfer_init())
 * -# bessel_at_x() at any time for computing a value j_l(x) at any x by interpolation
 * -# threej_free() at the end

 */

#include "threej.h"

int threej_init(
 		    struct precision * ppr,
		    struct bessels * pbs,
 		    struct threej * p3j
 		    )
{
 		      
  /*  By default, we use the same l and x samplings as in the Bessel module
    (used at 1st order).  */
  p3j->l_size = pbs->l_size;
  p3j->l = pbs->l;
  p3j->l_max = pbs->l_max;
  p3j->x_max = pbs->x_max;
  p3j->x_size = pbs->x_size;
  p3j->x_step = pbs->x_step;


  

  /* This is not really needed, is it? */
  // class_call(
  //   threej_fill_matrices(p3j),
  //   p3j->error_message,
  //   p3j->error_message);
 		      
  return _SUCCESS_;

}




int LOS_function_for_lm (
 		    struct precision * ppr,
		    struct bessels * pbs,
 		    struct threej * p3j
 		    )
{
  
  
  return _SUCCESS_;
  
}



























int threej_fill_matrices(
        struct threej * p3j
        ) {


  // Temporary variables
  int i=0, j=0, k=0, s=0;
  int m=0, m1=0, m2=0;
  int l=0, l1=0, l2=0;
  int m1_min, m1_max;
  double * result;
  int result_size;
  complex double increment=0.;
  complex double sum_so_far=0.;

  // Shortcuts to the matrices
  double complex ** xi; 
  double complex ** chi_ij; 
  double complex ** chi_ijk; 
  double complex ** chi_ijks;

  // Write down size of the matrices
  p3j->xi_size = 3;
  p3j->m1_size = (2*1+1);
  p3j->chi_ij_size = 3*3;
  p3j->m2_size = (2*2+1);
  p3j->chi_ijk_size = 3*3*3;
  p3j->m3_size = (2*3+1);  
  p3j->chi_ijks_size = 3*3*3*3;  
  p3j->m4_size = (2*4+1);  

  // Allocate memory for the matrices.  Note that they are allocated with calloc
  // since their elements are defined by sums (and therefore aregoing to be
  // accumulated).
  class_alloc(xi, p3j->m1_size*sizeof(double complex*), p3j->error_message);
  for(i=0; i<p3j->m1_size; ++i)
    class_calloc(xi[i], p3j->xi_size, sizeof(double complex), p3j->error_message);

  class_alloc(chi_ij, p3j->m2_size*sizeof(double complex*), p3j->error_message);
  for(i=0; i<p3j->m2_size; ++i)
    class_calloc(chi_ij[i], p3j->chi_ij_size, sizeof(double complex), p3j->error_message);

  class_alloc(chi_ijk, p3j->m3_size*sizeof(double complex*), p3j->error_message);
  for(i=0; i<p3j->m3_size; ++i)
    class_calloc(chi_ijk[i], p3j->chi_ijk_size, sizeof(double complex), p3j->error_message);

  class_alloc(chi_ijks, p3j->m4_size*sizeof(double complex*), p3j->error_message);
  for(i=0; i<p3j->m4_size; ++i)
    class_calloc(chi_ijks[i], p3j->chi_ijks_size, sizeof(double complex), p3j->error_message);  


  // ==============
  // = Compute xi =
  // ==============
  // There is nothing to compute here, since the 3-vectors xi[m] are the
  // the directions of the spherical harmonic Y1m (n), i.e. the so-called
  // spherical basis.  The expression of the three xi's is taken from Beneke
  // and Fidler (2010), eq. A.12
  xi[0][0] = 1./sqrt(2.) * (+1);
  xi[0][1] = 1./sqrt(2.) * (+IM);
  xi[0][2] = 0.;    

  xi[1][0] = 0.;
  xi[1][1] = 0.;
  xi[1][2] = 1.;    

  xi[2][0] = 1./sqrt(2.) * (-1);
  xi[2][1] = 1./sqrt(2.) * (+IM);
  xi[2][2] = 0.;    


  // ====================
  // =   Compute chi_ij =
  // ====================
  // The l=2 matrix is obtained starting from the l=1 one (that is xi).  This is a
  // common pattern, as we shall see.
  l = 2;
  l1 = 1;
  l2 = 1;

  for(m=-l; m<=l; ++m) {  
    
    // These will contain the 3j for all possible values of 'm1'.
    // The calculated 'm1' values are between the output parameters
    // m1_min and m1_max, which correspond to m1_min = fmax(-l1, m-l2),
    // m1_max = fmin(l1, m+l2);
    class_call(
      threej_m2 ( p3j,                    // In
                  l,l1,l2,-m,             // In
                  &m1_min, &m1_max,       // Out
                  &result,                // Out
                  &result_size),          // Out
              p3j->error_message,
              p3j->error_message);


    for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
        
        // Sum on m1
        for(m1=m1_min; m1<=m1_max; ++m1) {
        
          m2 = m-m1;

          // When addressing an array with an 'm' index, always keep in mind that you should
          // always use index_m = m - MIN(m), which in most cases is equal to m - (-l) = m+l
          chi_ij[m+l][3*i + j] += sqrt(10./3.) * pow(-1,-m) * xi[m1+l1][i] * xi[m2+l2][j] * result[m1-m1_min];


        } // End of sum on m1

      } // End of loop on j
    } // End of loop on i

    free(result);

  } // End of loop on m



  // ==================
  // = Compute xi_ijk =
  // ==================
  // The chi_ijk matrix is obtained starting from xi and chi_ij
  l = 3;
  l1 = 2;
  l2 = 1;
  
  for(m=-l; m<=l; ++m) {  

    // Compute 3j
    class_call(
      threej_m2 ( p3j,                    // In
                  l,l1,l2,-m,             // In
                  &m1_min, &m1_max,       // Out
                  &result,                // Out
                  &result_size),          // Out
              p3j->error_message,
              p3j->error_message);

    for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
        for(k=0; k<3; ++k) {

          // Sum on m1
          for(m1=m1_min; m1<=m1_max; ++m1) {
          
            m2 = m-m1;

            // Increment to the ijk,m element of the matrix.
            increment = - sqrt(21./5.) * pow(-1,-m) * chi_ij[m1+l1][3*i+j] * xi[m2+l2][k] * result[m1-m1_min];
            
            // Numerical noise makes it impossible to have exactly zero when subtracting equal
            // numbers.  This is exactly what happens here.  We get some values to have a very tiny
            // imaginary part (tipically around 1e-17) when they should be zero.  So, we manually set
            // to zero the imaginary part of the element being accumulated if the two elements
            // of the difference are different less than 1 part in 1e16.
            sum_so_far = chi_ijk[m+l][3*3*i+3*j+k];            
            if( (cimag(sum_so_far) != 0.) && ((fabs(1. + cimag(increment) / cimag(sum_so_far))) < DBL_SMALL_ILON*100) ) {
              if (p3j->threej_verbose>1) printf("l=%d, m=%d, i,j,k = (%d,%d,%d): set imaginary part to zero\n", l, m, i, j, k);
              chi_ijk[m+l][3*3*i+3*j+k] = creal(sum_so_far);
            } 
            else {
              chi_ijk[m+l][3*3*i+3*j+k] += increment;
            }
                      
          } // End of sum on m1

        } // End of loop on k
      } // End of loop on j
    } // End of loop on i

    free(result);

  } // End of loop on m
    
    
    
  // ===================
  // = Compute xi_ijks =
  // ===================
  // The chi_ijks matrix is obtained starting from chi_ik only
  l = 4;
  l1 = 2;
  l2 = 2;
  
  for(m=-l; m<=l; ++m) {  

    // Compute 3j
    class_call(
      threej_m2 ( p3j,                    // In
                  l,l1,l2,-m,             // In
                  &m1_min, &m1_max,       // Out
                  &result,                // Out
                  &result_size),          // Out
              p3j->error_message,
              p3j->error_message);

    for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
        for(k=0; k<3; ++k) {
          for(s=0; s<3; ++s) {
          
            // Sum on m1
            for(m1=m1_min; m1<=m1_max; ++m1) {
          
              m2 = m-m1;
          
              // Increment to the ijks,m element of the matrix
              increment = 18./sqrt(70.) * pow(-1,-m) * chi_ij[m1+l1][3*i+j] * chi_ij[m2+l2][3*k+s] * result[m1-m1_min];

              // Account for numerical noise (see comment for the chi_ijk case)
              sum_so_far = chi_ijks[m+l][3*3*3*i+3*3*j+3*k+s];            
              if( (cimag(sum_so_far) != 0.) && ((fabs(1. + cimag(increment) / cimag(sum_so_far))) < DBL_SMALL_ILON*100) ) {
                if (p3j->threej_verbose>1) printf("l=%d, m=%d, i,j,k,s = (%d,%d,%d,%d): set imaginary part to zero\n", l, m, i, j, k, s);
                chi_ijks[m+l][3*3*3*i+3*3*j+3*k+s] = creal(sum_so_far);
              } 
              else {
                chi_ijks[m+l][3*3*3*i+3*3*j+3*k+s] += increment;
              }
                        
            } // End of sum on m1

          } // End of loop on s
        } // End of loop on k
      } // End of loop on j
    } // End of loop on i

    free(result);

  } // End of loop on m
    
    
  // Make the structure matrices point to our aliases
  p3j->xi = xi;
  p3j->chi_ij = chi_ij;
  p3j->chi_ijk = chi_ijk;
  p3j->chi_ijks = chi_ijks; 

  return _SUCCESS_;

}

