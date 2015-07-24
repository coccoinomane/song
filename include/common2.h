/** @file common2.h
 * Include file that defines the precision structure and generic macros
 * used in SONG.
 */

#include "common.h"
#include "song_tools.h"
#include "slatec_3j_C.h"

#ifndef __COMMON2__
#define __COMMON2__

#define _SONG_VERSION_ "v1.0-beta1"

/* Maximum number of azimuthal numbers that can be computed. We set it so that the
factorial of 'm' never overflows, assuming a limit of 10^30. The factorial of m is
needed in bispectra2.c */
#define _MAX_NUM_AZIMUTHAL_ 14


/**
 * All precision parameters for the second-order part of SONG. 
 */
struct precision2
{

  /* Tolerance for the integration of the 2nd-order system. This parameter goes
  directly into the evolver as the parameter 'rtol' */
  double tol_perturb_integration_song;
  
  /* How many multipoles should we evolve at second-order? */
  int l_max_g; 
  int l_max_pol_g;
  int l_max_ur;
  int l_max_g_ten;  
  int l_max_pol_g_ten;

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

  /* Logical array of size ppr2->m_max_2nd_order. If m is contained in ppr2->m,
  then ppr2->compute[m] == _TRUE_ */
  int compute_m[_MAX_NUM_AZIMUTHAL_];

  /* ppr2->index_m[M] contains the index of 'M' inside ppr2->m. If M is not contained
  in ppr2->m, then ppr2->index_m[M]=-1. */
  int index_m[_MAX_NUM_AZIMUTHAL_];
    
  /* index_m_max[l] is the maximum allowed m (in ppr2->m) for a given l */
  int * index_m_max;



  // ==============================
  // =       Time samplings       =
  // ==============================
  
  /** Time at which the second-order system will start being evolved. By default it
  is zero, which means that the start time will be determined automatically. */
  double custom_tau_start_evolution; 
  
  /** Parameter used to determine the start time of the differential system at second
  order; set to a value smaller than one to start evolving the k-modes in the tight
  coupling regime. */
  double start_small_k_at_tau_c_over_tau_h_song; 

  /** Parameter used to determine the start time of the differential system at second
  order; set to a value smaller than one to start evolving the k-modes when they are
  outside the horizon. */
  double start_large_k_at_tau_h_over_tau_k_song; 

  /** Time step for the sampling of the second-order line-of-sight sources, in
  units of the timescale involved in the sources. For example, for the CMB, the
  timescale is set by the variations in the visibility function, while for the power
  spectrum it is simply the Hubble time. */
  double perturb_sampling_stepsize_song;
  
  /** Time sampling for the line of sight integration, to be used when the step of
  the sources, perturb_sampling_stepsize_song, is not enough to capture the frequent
  oscillations of the Bessel functions. */
  double tau_step_trans_song;  



  // ==================================
  // =           k1-k2 sampling       =
  // ==================================

  /* Parameters for the k-sampling at second order, using the same algorithm as CLASS */

  double k_min_tau0; /**< number defining k_min for the computation of scalar Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */
  double k_max_tau0_over_l_max; /**< number defining k_max for the computation of scalar Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */
  double k_step_sub; /**< linear step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */
  double k_step_super; /**< linear step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  
  double k_logstep_super; /**< logarithmic step in k space, pure number, used to determine the very small k-values */  
  double k_step_transition; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */

  /* Parameters for the k-sampling at second order, using a linear or logarithmic sampling */

  double k_min_custom;
  double k_max_custom;
  int k_size_custom;



  // ====================================
  // =           k3 sampling            =
  // ====================================

  int k3_size_min; /**< Minimum number of grid points for any (k1,k2) pair,
                      used when 'k3_sampling' is set to smart */
  int k3_size;     /**< Fixed number of grid points for any (k1,k2) pair,
                      used when 'k3_sampling' is set to lin or log */
  double q_linstep_song; /**< Upper bound on linear sampling step in k space for the transfer functions,
                            in units of one period of acoustic oscillation (2*pi/(tau0-tau_rec)) */



  // ========================================
  // =            Interpolations            =
  // ========================================

  /* How to interpolate the sources in the line of sight integral? */
  enum interpolation_methods sources_time_interpolation;
  enum interpolation_methods sources_k3_interpolation;



	// =============================
	// =       Bessels       =
	// =============================

  double bessel_j_cut_song;  /* Value of j_l1(x) below which it is approximated by zero (in the region x << l) */
	double bessel_J_cut_song;	/* Value of J_Llm(x) below which it is approximated by zero (in the region x << l) */
  double bessel_x_step_song; /* Linear step dx for sampling spherical Bessel functions j_l1(x) and functions J_Llm(x) */



	// ==========================
	// =       Bispectrum       =
	// ==========================






	// ==========================
	// =           Misc         =
	// ==========================

  ErrorMsg error_message;         /**< Zone for writing error messages */
  short store_transfers_to_disk;  /**< Should we store the transfer functions to disk? */
  short load_transfers_from_disk; /**< Should we load the transfer functions from disk? */
  short store_sources_to_disk;    /**< Should we store the source functions to disk? */
  short load_sources_from_disk;   /**< Should we load the source functions from disk? */
  short old_run; /**< set to _TRUE_ if the run was stored with a version of SONG smaller than 1.0 */

};  /* end of struct precision2 declaration */

#endif
