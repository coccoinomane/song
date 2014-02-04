/** @file threej.h Documented includes for ThreeJ module */

#ifndef __THREEJ__
#define __THREEJ__

#include "common.h"
#include "arrays.h"
#include "bessel.h"

#include "complex.h"
#undef I                  /* 'I' may already be taken... */
#define IM _Complex_I     /* ... hence, use 'IM' to denote sqrt(-1) */


/**
 * Structure containing everything about 3j symbols and multipole matrices
 * that other modules need to know.
 *
 * Once initialized by threej_init(), contains table of
 * spherical Bessel functions \f$ j_l(x) \f$).
 */

struct threej {

  
  // ======================
  // = Spherical matrices =
  // ======================
	// Matrices used to compute, respectively, the linear, quadratic, cubic, and quartic terms in the
	// photon direction in Boltzmann equation.  The matrices have, respectively, 2, 3 and 4 
	// indices that range from 1 to 3, plus a 'm' index that ranges from -l to l, where
	// l is the number of indices

  // **** XI *****
  double complex ** xi;
  int xi_size;
  int m1_size;

  // **** chi_ij *****
	double complex ** chi_ij;
  int chi_ij_size;
  int m2_size;
	// chi_ij[m][ 3*i + j]

  // **** chi_ijk *****
	double complex ** chi_ijk;
  int chi_ijk_size;
  int m3_size;  
	// chi_ijk[m][ 3*3*i + 3*j + k]

  // **** chi_ijks *****
	double complex ** chi_ijks;
  int chi_ijks_size;
  int m4_size;  
	// chi_ijks[m][ 3*3*3*i + 3*3*j + 3*k + s]







  /** @name - input parameters initialized by user in input module 
      (all other quantitites are computed in this module, given these 
      parameters and the content of the 'precision' structure) */
  
  //@{

  int l_max; /**< last l value */

  double x_max; /**< maximum value of x (always multiple of x-step) */

  double x_step; /**< step dx for sampling Bessel functions */

  short bessel_always_recompute; /**< if set to true, 3j are never read from / written in files */

 //@}

  /** @name - parameters defining uniquely the exact content of the Bessel table
      (hence when reading a file, will compare these values with needed value
      in order to take the decision to recompute or not) */

  //@{

  int l_size; /**< number of multipole values */
  int * l; /**< list of multipole values, l[index_l] */

  double j_cut; /**< value of \f$ j_l \f$ below which it is approximated by zero (in the region x << l) */
  int has_dj;   /**< set to true means j_l'(x) also need to be stored */

 //@}

  /** @name - Bessel table, and arrays necessary for reading it */

  //@{

  int * x_size; /**< x_size[index_l] is the number of x values for l[index_l]; hence *x_min[index_l]+x_step*(x_size[index_l]-1) = x_max */

  int x_size_max;

  double ** buffer; /**< buffer[index_l] is a pointer towards a memory zone containing x_min, all j_l(x) and all j_l''(x) for each value of l */

  double ** x_min; /**< x_min[index_l] is a pointer towards the minimum value of x for l[index_l], given j_cut; always a multiple of x_step */

  double ** j; /* (j[index_l])[index_x] is \f$ j_l(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */ 

  double ** ddj; /* (ddj[index_l])[index_x] \f$ j_l''(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x (in view of spline interpolation) */ 

  double ** dj; /* (dj[index_l])[index_x] is \f$ j_l'(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x */ 

  double ** dddj; /* (dddj[index_l])[index_x] \f$ j_l'''(x) \f$ for l[index_l] and x=x_min[index_l]+x_step*index_x (in view of spline interpolation) */

  //@}

  /** @name - technical parameters */

  //@{

  short threej_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  // int bessel_at_x(
  //      struct threej * p3j,
  //      double x,
  //      int l,
  //      double * j
  //      );

  int threej_init(
		  struct precision * ppr,
		  struct bessels * pbs,
		  struct threej * p3j
		  );

  int threej_free(
		  struct threej * p3j
		  );

  int threej_get_l_list(
			struct precision * ppr,
			struct threej * p3j
			);

  // int bessel_j_for_l(
  //         struct precision * ppr,
  //         struct threej * p3j,
  //         int index_l
  //         );

  int threej_fill_matrices(
	       struct threej * p3j
         );


    
#ifdef __cplusplus
}
#endif


#endif
