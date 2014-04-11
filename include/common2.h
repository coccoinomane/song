/** @file common2.h Generic libraries, parameters and functions used in SONG. */

#include "common.h"
#include "song_tools.h"      /* A bunch of functions needed by SONG */
#include "slatec_3j_C.h"


#ifndef __COMMON2__
#define __COMMON2__

#define _SONG_VERSION_ "v1.0"


/* The differential system at second-order has to be solved for a set of three wavemodes
(k1,k2,k3) whereby k3 has to be in the range |k1-k2|<=k3<=k1+k1. When k3 is too close
to the boundaries, numerical instabilities might arise such as nan's or larger than one
sines and cosines. In order to avoid that, we define here a safety distance between
k3 and the bounds. This safety distance is going to correspond to the largest scale
probed by SONG. Using k_scalar_min_tau0_2nd_order=1e-3, that correspond to k_min=1e-8,
it seems that setting _MIN_K3_DISTANCE_=1e-10 is ok. */
#define _MIN_K3_DISTANCE_ 1e-10
#define _MIN_K3_RATIO_ 100

/* The following macros are used to index many arrays in the code.  The idea is that
the (L,M) multipole is found at y[monopole_g + lm(L,M)].

We define a similar function to index the massive hierarchy, with l_max set to 2
since we retain only the n=0,1,2 beta-moments for baryons and cold dark matter. */

/* Used to index ppt2->sources, radiation species */
#define lm(l,m) ppt2->lm_array[l][ppr2->index_m[m]]
#define lm_bis(l,m) ppt2->lm_array[l][ppr2->index_m[m]]

/* Used to index ppt2->sources, massive species */
#define nlm(n,l,m) ppt2->nlm_array[n][l][ppr2->index_m[m]]
#define nlm_bis(n,l,m) ppt2->nlm_array[n][l][ppr2->index_m[m]]

/* Used to index the ppw2->rotation_1 and ppw2->rotation_2 arrays */
#define lm_quad(l,m) ppt2->lm_array_quad[l][m]

/* Used to index ptr2->transfer */
#define lm_cls(index_l,index_m) ptr2->lm_array[index_l][index_m]


/* Maximum number of azimuthal numbers that can be computed. We set it so that the
factorial of 'm' never overflows, assuming a limit of 10^30. The factorial of m is
needed in bispectra2.c */
#define _MAX_NUM_AZIMUTHAL_ 14


/**
 * All precision parameters for the second-order part of SONG.
 *  
 */
struct precision2
{
  

  /* Tolerance for the integration of the 2nd-order system.  This parameter goes
    directly into the evolver as the parameter 'rtol' */
  double tol_perturb_integration_2nd_order;
  

  /* How many multipoles should we evolve at second-order? */
  int l_max_g_2nd_order; 
  int l_max_pol_g_2nd_order; 
  int l_max_ur_2nd_order;
  int l_max_g_ten_2nd_order;  
  int l_max_pol_g_ten_2nd_order;

  /* How many multipoles should we keep in the quadratic sources? */
  int l_max_g_quadsources; 
  int l_max_pol_g_quadsources; 
  int l_max_ur_quadsources;
  int l_max_g_ten_quadsources;  
  int l_max_pol_g_ten_quadsources;
  int m_max_quadsources;
  
  /* How many multipoles should we keep in the line of sight integration? */
  int l_max_los_t;          
  int l_max_los_p;
  int l_max_los;

  int l_max_los_quadratic_t;
  int l_max_los_quadratic_p;
  int l_max_los_quadratic;

  /* Array containing the 'm' values for which we need to solve the system.
  m=0 -> scalar modes
  m=1 -> vector modes
  m=2 -> tensor modes  */
  int m[_MAX_NUM_AZIMUTHAL_];
  int m_size;
  
  /* Maximum 'm' contained in ppr2->m */
  int m_max_2nd_order;

  /* index_m_max[l] is the maximum allowed m (in ppr2->m) for a given l */
  int * index_m_max;

  /* Logical array of size ppr2->m_max_2nd_order. If m is contained in ppr2->m,
  then ppr2->compute[m] == _TRUE_ */
  int compute_m[_MAX_NUM_AZIMUTHAL_];

  /* ppr2->index_m[M] contains the index of 'M' inside ppr2->m. If M is not contained
  in ppr2->m, then ppr2->index_m[M]=-1. */
  int * index_m;
    

  // ==============================
  // =       Time samplings       =
  // ==============================
  
  /* Frequency of the time-sampling for the second-order line-of-sight sources */
  double perturb_sampling_stepsize_2nd_order;
  
  /* When should we start to evolve the system? */
  double start_small_k_at_tau_c_over_tau_h_2nd_order;
  double start_large_k_at_tau_h_over_tau_k_2nd_order;

  /* Time sampling for the line of sight integration */
  double tau_step_trans_2nd_order;  

  

  // ==================================
  // =           k1-k2 sampling       =
  // ==================================

  /* N.B. The k-sampling for the second-order sources is in the precision structure instead of the perturbs2 structure
    because we need the perturbs1 structure to access it. */


  // *** Parameters needed to reproduce CLASS sampling at first order
  double k_scalar_min_tau0_2nd_order; /**< number defining k_min for the computation of scalar Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */
  double k_scalar_max_tau0_over_l_max_2nd_order; /**< number defining k_max for the computation of scalar Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */
  double k_scalar_step_sub_2nd_order; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_scalar_linstep_super_2nd_order; /**< step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  
  double k_scalar_logstep_super_2nd_order; /**< logarithmic step in k space, pure number, used to determine the very small k-values */  
  double k_scalar_step_transition_2nd_order; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */



  // *** Parameters for the custom lin/log sampling

  /* Scalars */
  double k_min_scalars, cosk1k2_min_scalars;
  double k_max_scalars, cosk1k2_max_scalars;
  int k_size_scalars, cosk1k2_size_scalars;



  // ====================================
  // =           k3 sampling            =
  // ====================================
  int k3_size_min; /* Minimum number of grid points for any (k1,k2) pair, used when 'k3_sampling' is set to smart */
  int k3_size;     /* Fixed number of grid points for any (k1,k2) pair, used when 'k3_sampling' is set to lin or log */

  double k_step_trans_scalars_2nd_order; /* Upper bound on linear sampling step in k space for the transfer functions, in units of one period of acoustic oscillation at decoupling (usually chosen to be between k_scalar_step_sub and k_scalar_step_super) */


  // ========================================
  // =            Interpolations            =
  // ========================================

  /* How to interpolate the sources in the line of sight integral? */
  enum interpolation_methods sources_time_interpolation;
  enum interpolation_methods sources_k3_interpolation;



	// =============================
	// =       Bessels       =
	// =============================

  double bessel_j_cut_2nd_order;  /* Value of j_l1(x) below which it is approximated by zero (in the region x << l) */
	double bessel_J_cut_2nd_order;	/* Value of J_Llm(x) below which it is approximated by zero (in the region x << l) */
  double bessel_x_step_2nd_order; /* Linear step dx for sampling spherical Bessel functions j_l1(x) and functions J_Llm(x) */



	// ==========================
	// =       Bispectrum       =
	// ==========================






	// ==========================
	// =           Misc         =
	// ==========================

  ErrorMsg error_message;         /* Zone for writing error messages */
  short store_transfers_to_disk;  /* Should we store the transfer functions to disk? */
  short load_transfers_from_disk; /* Should we load the transfer functions from disk? */
  short store_sources_to_disk;    /* Should we store the source functions to disk? */
  short load_sources_from_disk;   /* Should we load the source functions from disk? */


};  /* end of struct precision2 declaration */

#endif
