/** 
 * Definitions for some support function to the perturbation2 module
 */

#ifndef __SONG_TOOLS__
#define __SONG_TOOLS__

#include "common.h"

#ifndef _SPLINE_NATURAL_
#define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
#endif

#ifndef _SPLINE_EST_DERIV_
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */
#endif

/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif



  /** 
   * Check whether a triad of numbers (l1,l2,l3) satisfies the triangular condition
   * |l1-l2| <= l3 <= l1+l2
   */
  int is_triangular_int (int l1, int l2, int l3);
  int is_triangular_double (double l1, double l2, double l3);


  /** 
   * Compute the 3j symbol
   * (    l1     l2     l3   )
   * ( -m2-m3    m2     m3   )
   * 
   * The result shall be stored into 'threej'.
   *
   */
  int threej_single(
         int l1, int l2, int l3, int m2, int m3, // In
         double *threej,                         // Out
         ErrorMsg errmsg       
         );



  /** 
   *  Compute the 3j symbol
   *  (  l1    l2     l3   )
   *  (  m1    m2  -m1-m2  )
   *  for all allowed values of 'm2'.
   *  
   *  This function has double precision.
   *
   */
  int threej_m2(
         int l1, int l2, int l3, int m1,    // In
         int *m2_min, int *m2_max,          // Out, 'm2' limits for which result was computed
         double **result,                   // Out
         int *result_size,                  // Out
         ErrorMsg errmsg         
         );

  /** 
   *  Compute the 3j symbol
   *  (    l1     l2     l3   )
   *  ( -m2-m3    m2     m3   )
   *  for all allowed values of 'L1'.
   *  
   *  The result shall be stored into the 'result' array, that needs to be preallocated
   *  to contain at least l2+l3 - abs(l2-l3) + 1 elements.
   *  
   *  This function has double precision.
   *
   */
  int threej_l1(
         int l2, int l3, int m2, int m3,    // In
         int *l1_min, int *l1_max,          // Out, 'm2' limits for which result was computed
         double **result,                   // Out, has size equal to 'result_size' (should be preallocated)
         int *result_size,                  // Out
         ErrorMsg errmsg
         );



  /** 
   * Compute the 6j symbol
   * {   l1     l2     l3   }
   * {   l4     l5     l6   }
   * for all allowed values of 'l1'.
   * 
   * The result shall be stored into the 'result' array, that needs to be preallocated
   * to contain at least MIN(L2+L3,L5+L6) - MAX(ABS(L2-L3),ABS(L5-L6)) + 1 elements.
   *
   */

  int sixj_l1(
         int l2, int l3, int l4, int l5, int l6,    // In
         int *l1_min, int *l1_max,                  // Out, 'l1' limits for which the result is computed
         double **result,                           // Out, has size equal to 'result_size' (should be preallocated)
         int *result_size,                          // Out, equal to l1_max + l1_min + 1.
         ErrorMsg errmsg       
         );






  /**      
   *  Given the arguments l >= 0.0D0, x >= 0.0D0, compute
   *  the value of the spherical Bessel function
   *
   *        j_l( x )
   *
   *  for values of l going from l to l+N-1.
   *
   *  The result shall be stored into the 'result' array, that needs to be preallocated
   *  to contain at least N elements.
   *
   *  This function has single precision.  In order to use double
   *  precision, change float -> double and us dbesj_ instead of
   *  besj_.
   *
   */
   int besselj_l1(
          float l,                          // In
          float x,                          // In
          int N,                             // In
          float **result,                   // Out (array of size N)
          ErrorMsg errmsg
          );

  /**      
   *  Given the arguments l >= 0.0D0, x >= 0.0D0, compute
   *  the value of the Bessel function
   *
   *        J_l( x )
   *
   *  for values of l going from l to l+N-1.
   *
   *  This function has single precision.  In order to use double
   *  precision, change float -> double and us dbesj_ instead of
   *  besj_.
   *
   */
   int besselJ_l1(
          float l,                          // In
          float x,                          // In
          int N,                            // In
          float * result,                   // Out (array of size N)
          ErrorMsg errmsg
          );
 
  /**
   * Compute spherical Bessel function j_l(x) for a given l and x.
   *
   * Inspired from Numerical Recipies. This is the same as the function
   * bessel_j in the Bessel structure, but without requiring pbs as an
   * argument, so that it can be called from anywhere.
   */
  double spherical_bessel_j(
         int l,
         double x
         );





  /*   Function from the Slatec library to compute the 3j symbol for
     all the allowed values of 'm2'.  IMPORTANT: this is a Fortran
     function contained in the file tools/slatec_3j_f90.f90     */
  void drc3jm_ (double *l1, double *l2, double *l3,
                double *m1, double *m2_min, double *m2_max,
                double *thrcof, int *ndim, int *ier);
	       

  /*   Function from the Slatec library to compute the 3j symbol for
     all the allowed values of 'l1'.  IMPORTANT: this is a Fortran
     function contained in the file tools/slatec_3j_f90.f90     */
  void drc3jj_ (double *l2, double *l3,
                double *m2, double *m3, double *l1_min, double *l1_max,
                double *thrcof, int *ndim, int *ier);

  /*    Function from the SLATEC library to compute the J Bessel function J_l ( x )
    for values of l going from l to l+N-1.  IMPORTANT: this is a Fortran
    function contained in the file tools/slatec_3j_f90.f90   */
  void dbesj_ (double *x, double *l, int *N, double *result, int *NZ); // double precision
  void besj_ (float *x, float *l, int *N, float *result, int *NZ);     // single precision




  /* Associated Legendre polynomials P_LM, from Numerical Recipes, Third edition, pag.294, Press et al. 2002. */
  double plegendre_lm(int l, int m, double x);

  /* Associated Legendre polynomials, rescaled so that they do not include the (1-x)^(m/2) factor */
  double plegendre_lm_rescaled(int l, int m, double x);
  double plegendre_lm_rescaled_analytically (int l, int m, double x);

  /* Legendre polynomials P_L, from alglib-3.6.0 */
  double plegendre (int n, double x);

  /* Function to calculate second derivatives of ppt->sources at the nodes.
    This function reflects the very specific indexing pattern of ppt->sources, and therefore
    is not easily recycleable outside CLASS. */
  int spline_sources_derivs(
			     double * x, /* vector of size tau_size */
			     int tau_size,
			     double *** y_array, /* array of size tau_size*tp_size with elements 
						  										y_array[index_tau*tp_size+index_tp] */
			     int tp_size,   
			     double *** ddy_array, /* array of size tau_size*tp_size */
			     short spline_mode,
           int index_mode,
           int index_ic,           
           int index_k,
           int k_size,
			     ErrorMsg errmsg
           );
           
  int spline_sources_derivs_two_levels(
  			     double * x, /* vector of size tau_size */
  			     int tau_size,
  			     double ** y_array,
  			     int tp_size,   
  			     double ** ddy_array,
  			     short spline_mode,
  			     ErrorMsg errmsg
             );

  int spline_sources_interpolate_two_levels(
  			     double * x_array,
  			     int tau_size,
  			     double ** y_array,
  			     double ** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
  			     ErrorMsg errmsg
             );

  int spline_sources_interpolate_growing_closeby(
  			     double * x_array,
  			     int tau_size,
  			     double *** y_array,
  			     double *** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
             int index_mode,
             int index_ic,           
             int index_k,
             int k_size,
  			     ErrorMsg errmsg
  			     );
           
           
  int spline_sources_interpolate_two_levels_growing_closeby(
  			     double * x_array,
  			     int tau_size,
  			     double ** y_array,
  			     double ** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
  			     ErrorMsg errmsg
             );

  /* Function to interpolate ppt->sources at a specific time.  The result will be stored
  in the vector 'result' which should already be allocated.  The function will write
  ppt->tp_size doubles inside result.
  This function reflects the very specific indexing pattern of ppt->sources, and therefore
  is not easily recycleable outside CLASS. */
  int spline_sources_interpolate(
  			     double * x_array,
  			     int x_size,
  			     double *** y_array,
  			     double *** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
             int index_mode,
             int index_ic,           
             int index_k,
             int k_size,
  			     ErrorMsg errmsg
             );

  int array_interpolate_spline_fake(
			       double * x_array,
			       int n_lines,
			       double * array,
			       double * array_splined,
			       int n_columns,
			       double x,
			       int * last_index,
			       double * result,
			       int result_size, /** from 1 to n_columns */
			       ErrorMsg errmsg);
	


  /* Fill an array with logarithmically spaced points. */
  int log_space (double * xx, double x_min, double x_max, int n_points);
  
  /* Fill an array with linearly spaced points */
  int lin_space (double * xx, double x_min, double x_max, int n_points);

  /* Matrix inversion, credits to Christopher M. Brown
  (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html) */
  double Determinant(double **a,int n);
  void CoFactor(double **a,int n,double **b);
  void Transpose(double **a,int n);
  void InverseMatrix(double **in,int n,double **out);




  // ===============================================================
  // =                    Coupling coefficients                    =
  // ===============================================================



  /* C coupling coefficients.  They appear in the Boltzmann hierarchy for the photon
  temperature (see eqs. A.11 and 2.18 of Beneke & Fidler 2011). */
  double coupling_c_plus (int l, int m1, int m);
  double coupling_c_minus (int l, int m1, int m);

  /* D coupling coefficients.  They appear in the Boltzmann hierarchy for the photon
  polarization (see eqs. 2.17 and 2.19 of B&F 2011 and eq. 141 of B&F 2010). */
  double coupling_d_zero (int l, int m1, int m);
  double coupling_d_minus (int l, int m1, int m);
  double coupling_d_plus (int l, int m1, int m);    

  /* Generic coupling coefficient, obtained as the product of two 3j-symbols */
  int coupling_general (
    int l2, int l3, int m1, int F,
    double * three_j_000, /* should be preallocated with at least l2_max doubles */
    int three_j_000_size,
    double * three_j_mmm, /* should be preallocated with at least m1_max doubles */
    int three_j_mmm_size,
    int * l1_min, int * l1_max,
    int * m2_min, int * m2_max,
    double ** result,     /* should be preallocated with at least l2_max*m1_max doubles */
    ErrorMsg errmsg 
    );


  // ===============================================================
  // =                  Multipole related functions                =
  // ===============================================================
  
  /* Return the index corresponding to an l,m pair, considering m in [0,min(l,m_max)] */ 
  int multipole2offset_l_m(int l, int m, int m_max);

  /* Return the number of elements in the l,m hierarchy with l in [0,l_max] and m in [0,min(l,m_max)] */
  int size_l_m(int l_max, int m_max);

  /* Return the index corresponding to an l,m pair, considering that the M-index is constrained to the
 	  values contained in m_vec.  Note that if m_vec = (0,1,2, ..., m_max), this function should give the
 		same result as multipole2offset_l_m(L, M, m_max). */ 
	int multipole2offset_l_indexm(int L, int M, int * m_vec, int m_size);

  /* Inverse function of multipole2offset_l_indexm. Given an offset, l_max and a list of m's,
  determine the unique (l,m) pair associated with that offset. The output is written into
  L, index_M. */
  int offset2multipole_l_indexm (int offset, int l_max, int * m_vec, int m_size,
                                 int * L, int * index_M);


  /* Return the number of elements in the l,m hierarchy considering that the m-index is constrained to the
 	  values contained in m_vec */
	int size_l_indexm(int l_max, int * m_vec, int m_size);

  /* Return the index corresponding to an L,M pair, considering that the L-index and M-index are
    constrained to the values contained in l_vec and m_vec, respectively. Note that if
    l_vec = (0,1,2,...,l_max) and m_vec = (0,1,2, ..., m_max), this function should give the
	  same result as multipole2offset_l_m(L, M, m_max). */
  int multipole2offset_indexl_indexm(int L, int M, int * l_vec, int l_size, int * m_vec, int m_size);

  /* Return the number of elements in the l,m hierarchy considering that the m-index is constrained to the
 	  values contained in m_vec */
	int size_indexl_indexm(int * l_vec, int l_size, int * m_vec, int m_size);

  /* Inverse function of multipole2offset_indexl_indexm. Given an offset and two lists of possible l's and m's,
    determine the unique (l,m) pair associated with that offset. The output is written into
    L, M, index_L, index_M */
  int offset2multipole_indexl_indexm(int offset, int * l_vec, int l_size, int * m_vec, int m_size,
                       int * index_L, int * index_M);

  /* Index the massive hierarchy with the tree indices n,l,m.  'n' can be any positive integer,
    l should be a positive integer smaller than 'n', and 'm' should be a positive integer smaller
    than 'l'. */
  int multipole2offset_unconstrained_n_l_m(int n, int l, int m, int l_max, int m_max);

  /* This function is the same as 'multipole2offset_unconstrained_n_l_m', but it only allows l to be equal to 0 or n. */
  int multipole2offset_n_l_m(int n, int l, int m, int l_max, int m_max);
  int size_n_l_m(int n_max, int l_max, int m_max);

  /* This function is the same as 'multipole2offset_n_l_m', but it only allows m to be in the list m_vec */
  int multipole2offset_n_l_indexm(int N, int L, int M, int l_max, int * m_vec, int m_size);
  int size_n_l_indexm (int n_max, int l_max, int * m_vec, int m_size);






#ifdef __cplusplus
}
#endif

#endif
