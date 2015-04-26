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
  

  // ==============================================================================
  // =                                3J and 6J symbols                           =
  // ==============================================================================

  int threej_single(
         int l1, int l2, int l3, int m2, int m3, // In
         double *threej,                         // Out
         ErrorMsg errmsg       
         );

  int threej_ratio_L_recursive (
        int l1, int l2, int l3, int N,     // In
        double *result,                    // Out, should be allocated with M+1 elements
        ErrorMsg errmsg
        );
        
  int threej_ratio_L (
        int l1, int l2, int l3,            // In
        int N1, int N2, int N3,            // In
        double *result,                    // Out
        ErrorMsg errmsg
        );

  int threej_ratio_L1 (
        int l1, int l2, int l3, int N,     // In
        double *result,                    // Out
        ErrorMsg errmsg
        );

  int threej_ratio_M_recursive (
        int l1, int l2, int l3, int M,     // In
        double *result,                    // Out, should be allocated with M+1 elements
        ErrorMsg errmsg
        );
        
  int threej_ratio_M (
        int l1, int l2, int l3, int M,     // In
        double *result,                    // Out
        ErrorMsg errmsg
        );
        
  int threej_A (
        int l1, int l2, int l3,
        int m1,
        double *result,
        ErrorMsg errmsg
        );

  int threej_B (
        int l1, int l2, int l3,
        int m1, int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_C (
        int l1, int l2, int l3,
        int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_D (
        int l1, int l2, int l3,
        int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_A_factor (
        int l1, int l2, int l3,
        int n,
        double *result,
        ErrorMsg errmsg
        );

    

  // ======================================================================
  // =                           Bessel functions                         =
  // ======================================================================

  double spherical_bessel_j(
         int l,
         double x
         );
         

  // ====================================================================================
  // =                                Coupling factors                                  =
  // ====================================================================================

  double coupling_c_plus (int l, int m1, int m);
  double coupling_c_minus (int l, int m1, int m);
  double coupling_d_zero (int l, int m1, int m);
  double coupling_d_minus (int l, int m1, int m);
  double coupling_d_plus (int l, int m1, int m);    

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
    

  // ============================================================================================
  // =                                 Legendre polynomials                                     =
  // ============================================================================================

  double plegendre_lm(int l, int m, double x);
  double plegendre_lm_rescaled(int l, int m, double x);
  double plegendre_lm_rescaled_analytically (int l, int m, double x);
  double plegendre (int n, double x);


  // ============================================================================================
  // =                             Multipole related functions                                  =
  // ============================================================================================
  
  int multipole2offset_l_m(int l, int m, int m_max);
  int size_l_m(int l_max, int m_max);
	int multipole2offset_l_indexm(int L, int M, int * m_vec, int m_size);
  int offset2multipole_l_indexm (int offset, int l_max, int * m_vec, int m_size,
                                 int * L, int * index_M);
	int size_l_indexm(int l_max, int * m_vec, int m_size);
  int multipole2offset_indexl_indexm(int L, int M, int * l_vec, int l_size, int * m_vec, int m_size);
	int size_indexl_indexm(int * l_vec, int l_size, int * m_vec, int m_size);
  int offset2multipole_indexl_indexm(int offset, int * l_vec, int l_size, int * m_vec, int m_size,
                       int * index_L, int * index_M);
  int multipole2offset_unconstrained_n_l_m(int n, int l, int m, int l_max, int m_max);
  int multipole2offset_n_l_m(int n, int l, int m, int l_max, int m_max);
  int size_n_l_m(int n_max, int l_max, int m_max);
  int multipole2offset_n_l_indexm(int N, int L, int M, int l_max, int * m_vec, int m_size);
  int size_n_l_indexm (int n_max, int l_max, int * m_vec, int m_size);
  

  // =====================================================================================
  // =                                     Interpolation                                 =
  // =====================================================================================

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
	


  // ============================================================================
  // =                                 Sampling related                         =
  // ============================================================================ 

  int log_space (double * xx, double x_min, double x_max, int n_points);
  int lin_space (double * xx, double x_min, double x_max, int n_points);



  // ===================================================================================
  // =                                Matrix operations                                =
  // ===================================================================================

  double Determinant(double **a,int n);
  void CoFactor(double **a,int n,double **b);
  void Transpose(double **a,int n);
  void InverseMatrix(double **in,int n,double **out);
  void PrintMatrix(double **in,int n);
  

  // ================================================================================
  // =                              Assert functions                                =
  // ================================================================================

  int is_triangular_int (int l1, int l2, int l3);
  int is_triangular_double (double l1, double l2, double l3);


  // ================================================================================
  // =                               Fortran functions                              =
  // ================================================================================
  
  /* The following functions are Fortran procedures from the Slatec library
  in the file tools/slatec_3j_f90.f90 */
  
  /**
   * Function from the Slatec library to compute the 3j symbol for
   * all the allowed values of 'm2'. 
   */
  void drc3jm_ (double *l1, double *l2, double *l3,
                double *m1, double *m2_min, double *m2_max,
                double *thrcof, int *ndim, int *ier);
	       

  /** 
   * Function from the Slatec library to compute the 3j symbol for
   * all the allowed values of 'l1'.  IMPORTANT: this is a Fortran
   * function contained in the file tools/slatec_3j_f90.f90
   */
  void drc3jj_ (double *l2, double *l3,
                double *m2, double *m3, double *l1_min, double *l1_max,
                double *thrcof, int *ndim, int *ier);

  /**
   * Function from the SLATEC library to compute the J Bessel function J_l ( x )
   * for values of l going from l to l+N-1.  IMPORTANT: this is a Fortran
   * function contained in the file tools/slatec_3j_f90.f90
   */
  void dbesj_ (double *x, double *l, int *N, double *result, int *NZ); // double precision
  void besj_ (float *x, float *l, int *N, float *result, int *NZ);     // single precision


  // ========================================================================
  // =                                 Misc                                 =
  // ========================================================================

  double identity_double (double x);

  int sign_int (int x);

  int ordering_int (
        int * n,            /* In */
        int * ordering,     /* Out */
        ErrorMsg errmsg
        );

  int reorder_int (
        int * n,             /* In/Out */
        int * ordering,      /* In */
        ErrorMsg errmsg
        );


#ifdef __cplusplus
}
#endif


// ===============================================================================
// =                              Preprocessor macros                            =
// ===============================================================================

/* Needed by spherical_bessel_j */
#define _GAMMA1_ 2.6789385347077476336556
#define _GAMMA2_ 1.3541179394264004169452


#endif
