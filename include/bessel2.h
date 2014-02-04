/** @file bessel2.h Documented includes for Bessel module */

#ifndef __BESSEL2__
#define __BESSEL2__

#include "bessel.h"
#include "common2.h"
#include "perturbations2.h"


/* What projection function to compute? */
enum projection_function_types {
  J_TT,
  J_EE,
  J_EB
};


/**
 * Structure containing everything about spherical bessel functions 
 * that other modules need to know.
 *
 * Once initialized by bessel2_init(), contains table of
 * spherical Bessel functions \f$ j_l(x) \f$).
 */

struct bessels2 {


  /* The following definitions are used to compute and store the geometrical tables needed for
  the line of sight integration at second order.  At second order we have mode coupling, hence
  we need a combination of Bessel functions and 3j-symbols rather than just Bessel functions.
  We store such coefficients for each value of (l,m) and argument x inside the table double ***** J.
  There is an extra level which indexes different kinds of J's: TT, EE, EB. The formulas for these
  projection functions can be found in eq. B.12 of Beneke and Fidler (2010). */
  
  
  // *** Flags
  short has_J_TT;              // Whether to compute and store the intensity J's or not
  short has_J_EE;              // Whether to compute and store the polarisation J's or not
  short has_J_EB;              // Whether to compute and store the mixing J's or not
  short extend_l1_using_m;     // Should we compute extend the pbs2->l1 array to take into account high values of m?


  // *** Domains
  long int count_allocated_Js;   // Number of doubles stored in the projection function arrays.
                                 
  double j_l1_cut;               // Value of j_l1(x) below which it is approximated by zero (in the region x << l)
  double J_Llm_cut;              // Value of J_Llm(x) below which it is approximated by zero (in the region x << l)
                                 
  int L_max;                     // Maximum value of L for which to compute the J_Llm(x)                                 
  int L_size;                    // Number of L's to compute J_Llm(x) for.  This is just ppr->l_max_los+1
  int * L;                       // Array containing the L's to be computed J_Llm(x) for
                                 
  int m_size;                    // Number of m's to compute J_Llm(x) for.  This is just ppr->m_size.
  int * m;                       // Array containing the m's to be computed J_Llm(x) for.  This is just ppr->m.
                                 
  double * xx;                   // xx is the grid where J_Llm(x) is sampled.  To each (L,l,m) it corresponds a different starting index in x.
  double xx_step;                // Linear step dx for sampling the J_Llm Bessel functions
  int xx_size;                   // Size of xx. This is determined by pbs->x_max and pbs->x_step.
  double xx_max;                 // Maximum value of xx (always multiple of xx-step).  Determined in input.c



  // *** Arrays and variables related to the projection functions
  
  int index_J_TT;                // do we need the temperature projection functions?
  int index_J_EE;                // do we need the E-mode polarization projection functions?
  int index_J_EB;                // do we need the projection functions for the mixing of polarisations?
  int J_size;                    // number of projection functions to compute (so far, TT, EE, BB, EB, BE)
  
  int **** index_xmin_J;         // index_x_min_J[index_J][index_L][index_l][index_m] is the index of the xx array where J_Llm(x) starts to be non-negligible
  int **** x_size_J;             // x_size_J[index_J][index_L][index_l][index_m] is the number of x values we sample J_Llm(x) in
  double **** x_min_J;           // x_min_J[index_J][index_L][index_l][index_m] is the first x where you have a non-negliglible value for J_Llm(x)
  int x_size_max_J;              // maximum value of x_size_J[index_J][index_L][index_l][index_m] over L,l,m
                                 
  double ***** J_Llm_x;          // J_Llm_x[index_J][index_L][index_l][index_m][index_x] with x=x_min_J[index_J][index_L][index_l][index_m]+x_step*index_x
  double ***** ddJ_Llm_x;        // Same indexing as J_Llm_x, used for spline interpolation

  short * has_allocated_J;       // was the memory for the index_J projection functions allocated?


  // *** Sampling of j_l1                                 
  int * l1;                      // A multipole list that includes all points in pbs->l, plus more needed in the computation of J_Llm(x)
  int * index_l1;                // pbs2->index_l1[l1] is the index of 'l' inside pbs2->l1. If 'l1' is not contained in pbs2->l1, then pbs2->index_l1[l1]=-1
  int l1_size;                   
  double ** j_l1;                /*  Spherical Bessel function j_l(x), indexed as pbs->j_l1[index_l1][index_x].  
                                    The l1 level of pbs->j_l1 should be addressed the same way as pbs->l1[index_l1], while
                                    the x level should be addressed as 'index_x - pbs->index_xmin_l1[index_l1]', where index_x
                                    is the index of pbs->xx. */
  double ** ddj_l1;              // Same indexing as j_l1, used for spline interpolation
  int * index_xmin_l1;           // index_xmin_l1[index_l1] is the index of pbs->xx where j_l1(x) starts to be non-negligible
  int * x_size_l1;               // x_size_l1[index_l1] is the number of x values we sample j_l1(x) in
  double * x_min_l1;             // x_min_l1[index_l1] is the first x where you have a non-negliglible value for j_l1(x)


  // *** Technical parameters
  short bessels2_verbose;        // verbosity flag (none if set to zero)
  ErrorMsg error_message;        // zone for writing error messages


};



/**
 * Structure containing the result of the computation of the spherical Bessel functions
 *
 *          j_l1(x)
 *
 * and of the 3j-symbols
 *
 *      (    l     l1      L   )
 *      (    0      0      0   )   ,
 *
 *      (    l     l1      L   )
 *      (   -m      0      m   )
 *
 *
 * for all the needed values of 'l1', and for a given set of (L,l,m) indices.
 *
 * This structure is required by bessel2_J_Llm.  We fill it inside bessel2_J_for_Llm.
 *
 */
struct J_Llm_data {

  float * bessels;            // Temporary l1-array to hold the spherical Bessel funcions j_l1(x) neeeded to compute J_Llm(x)
  double * first_3j;          // Temporary l1-array to hold the 3j-symbol (l l1 L)(0  0  0) neeeded to compute J_Llm(x)
  double * second_3j;         // Temporary l1-array to hold the 3j-symbol (l l1 L)(m  0 -m) neeeded to compute J_Llm(x)
  
  int l1_size;                // Allowed values of l1 (top left index in the 3j's)
  int l1_min;                 // Minimum allowed values of l1 (top left index in the 3j's)
  int l1_max;                 // Maximum allowed values of l1 (top left index in the 3j's)
  
  int index_l1_min;           // Position of l1_min in the array pbs->l1
  
};




/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int bessel2_init(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      );

  int bessel2_free(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      );

  int bessel2_get_l1_list(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      );

  int bessel2_J_Llm_at_x(
      struct bessels * pbs,
      struct bessels2 * pbs2,
      double x,
      int index_L,
      int index_l,
      int index_m,
      double * J_Llm_x
      );

  int bessel2_J_for_Llm(
       struct precision * ppr,
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       int index_J,
       int index_L,
       int index_l,
       int index_m
       );

  int bessel2_J_Llm(
         struct precision2 * ppr2,
         struct bessels * pbs,
         struct bessels2 * pbs2,
         enum projection_function_types projection_function,
         int L,
         int l,
         int m,
         int index_x,
         struct J_Llm_data * bessel_3j_data,
         double * J_Llm_x
         );

  int bessel2_l1_at_x(
      struct bessels2 * pbs2,
      double x,
      int index_l1,
      double * j_l1
      );

  int bessel2_l1_at_x_linear(
      struct bessels2 * pbs2,
      double x,
      int index_l1,
      double * j_l1
      );
    
  int bessel2_j_for_l1(
         struct precision * ppr,
         struct precision2 * ppr2,
         struct bessels * pbs,
         struct bessels2 * pbs2,
         int index_l1
         );

  int bessel2_get_xx_list(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      );

    
  int bessel2_convolution(
      struct precision * ppr,
      struct bessels2 * pbs2,
      double * kk,
      double * delta_kk,
      int k_size,
      double * f,
      double * g,
      int index_l,
      double r,
      double * integral,
      ErrorMsg error_message
      );
    
#ifdef __cplusplus
}
#endif



#endif
