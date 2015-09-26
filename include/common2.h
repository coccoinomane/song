/** @file common2.h
 * Include file that defines the precision structure and generic macros
 * used in SONG.
 */

#include "common.h"
#include "song_tools.h"
#include "slatec_3j_C.h"

#ifndef __COMMON2__
#define __COMMON2__

#define _SONG_VERSION_ "v1.0-beta2"
#define _SONG_URL_ "https://github.com/coccoinomane/song"

/** Maximum number of azimuthal numbers that can be computed. We set it so that the
factorial of m never overflows, assuming a limit of 10^30. The factorial of m is
needed in bispectra2.c */
#define _MAX_NUM_AZIMUTHAL_ 14


/**
 * All precision parameters for the second-order part of SONG. 
 */
struct precision2
{

  /** Tolerance for the integration of the second-order system. This parameter goes
  directly into the evolver as the parameter rtol */
  double tol_perturb_integration_song;


  // ====================================================================================
  // =                                   Multipoles                                     =
  // ====================================================================================

  int l_max_g;         /**< Number of multipoles to evolve for the photon intensity */
  int l_max_pol_g;     /**< Number of multipoles to evolve for the photon polarisation */
  int l_max_ur;        /**< Number of multipoles to evolve for the neutrinos */
  int l_max_boltzmann; /**< maximum number of multipoles to keep in any of the second-order Boltzmann hierarchies */

  int l_max_g_quadsources;      /**< Number of multipoles to include in the quadratic part of the Boltzmann hierarchy for the photon intensity. Set to -1 to include all quadratic sources up to ppr2->l_max_g. */
  int l_max_pol_g_quadsources;  /**< Number of multipoles to include in the quadratic part of the Boltzmann hierarchy for the photon polarisation. Set to -1 to include all quadratic sources up to ppr2->l_max_pol_g. */
  int l_max_ur_quadsources;     /**< Number of multipoles to include in the quadratic part of the Boltzmann hierarchy for the neutrinos. Set to -1 to include all quadratic sources up to ppr2->l_max_ur. */
  
  int l_max_los_t;  /**< Number of multipoles to keep in the line of sight sources for the photon intensity */
  int l_max_los_p;  /**< Number of multipoles to keep in the line of sight sources for the photon polarisation */
  int l_max_los;    /**< Maximum number of multipoles to keep in the line of sight sources for any species */

  int l_max_los_quadratic_t;  /**< Number of multipoles to keep in the quadratic part of line of sight sources for the photon intensity */
  int l_max_los_quadratic_p;  /**< Number of multipoles to keep in the quadratic part of line of sight sources for the photon polarisation */
  int l_max_los_quadratic;    /**< Maximum number of multipoles to keep in the quadratic part of the line of sight sources for any species */

  /** Array containing the 'm' values for which we need to solve the system.
  m=0  ->  scalar modes
  m=1  ->  vector modes
  m=2  ->  tensor modes
  m=3+ ->  no name because the metric does not have these modes */
  int m[_MAX_NUM_AZIMUTHAL_];
  int m_size; /**< Number of azimuthal modes to evolve */
  int m_max_2nd_order;  /**< Maximum m contained in ppr2->m */

  /** Logical array of size ppr2->m_max_2nd_order; if compute_m[M]==_TRUE_, then SONG will solve
  the Boltzmann equation for that azimuthal mode M. In other words, if M is contained in ppr2->m,
  then ppr2->compute[m] == _TRUE_. */
  int compute_m[_MAX_NUM_AZIMUTHAL_];  

  /** If the azimuthal mode M is in the list of modes computed by SONG, ppr2->m, then
  ppr2->index_m[M] is the index of M inside ppr2->m. Otherwise, ppr2->index_m[M]=-1. */
  int index_m[_MAX_NUM_AZIMUTHAL_];
    
  /** For a given multipole L, ppr2->index_m_max[L] is index in ppr2->m associated to the
  maximum allowed M in ppr2->m. For example, if SONG is going to compute scalar (m=0) and
  tensor (m=2) modes, then:
    - ppr2->index_m_max[0]=0
    - ppr2->index_m_max[1]=0
    - ppr2->index_m_max[2]=1
    - ppr2->index_m_max[L]=1 for L>2.
  This array is used to loop over l and m. */
  int * index_m_max;



  // ====================================================================================
  // =                                 Time samplings                                   =
  // ====================================================================================
  
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
  
  /** If larger than one, this parameter adds points to the time sampling of the
  CMB source function at late times. It appears as a multiplicative factor for the 
  Hubble rate in perturbations2_timesampling_for_sources(). The idea is to make
  the integration grid in the transfer module more dense in order to capture the
  oscillations of the projection functions at late times. This parameter is ignored
  if the CMB is not asked. */
  double perturb_sampling_late_time_boost;
  
  /** Frequency of the time sampling for the line of sight integration. This parameter
  is used when the time-sampling of the sources, ppt2->tau_sampling, is not dense enough
  to capture the frequent oscillations of the Bessel functions. */
  double tau_linstep_song;  



  // ====================================================================================
  // =                                   k1-k2 sampling                                 =
  // ====================================================================================

  /**
   * Parameters for the k-sampling of the CMB observables at second-order, such as
   * the angular power spectrum C_l and the CMB bispectrum B_l1_l2_l3.
   */
  //@{
  double k_min_tau0; /**< number defining k_min for the computation of scalar Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one */

  double k_max_tau0_over_l_max; /**< number defining k_max for the computation of scalar Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two */

  double k_step_sub; /**< linear step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling */

  double k_step_super; /**< linear step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling */  

  double k_logstep_super; /**< logarithmic step in k space, pure number, used to determine the very small k-values */  

  double k_step_transition; /**< dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision. */
  //@}

  /**
   * Parameters for the k-sampling of the LSS (large scale structure) observables
   * at second-order, such as the power spectrum P(k) and the matter bispectrum
   * B_k1_k2_k3.
   *
   * All these parameters are overridden by the CMB ones for all k values smaller
   * than k_max_cmb = l_scalar_max * k_max_tau0_over_l_max / tau0.
   */
  //@{
  double k_per_decade_for_pk; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region*/

  double k_per_decade_for_bao; /**< if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)*/

  double k_bao_center; /**< in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak) */

  double k_bao_width; /**< in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak) */
  //@}

  /* Parameters for the k-sampling at second order, using a linear or logarithmic sampling */

  double k_min_custom; /**< User-provided minimum value for the wavemode k (in ppt2->k); ignored unless ppt2->k_sampling==lin_k_sampling or log_k_sampling */
  double k_max_custom; /**< User-provided maximum value for the wavemode k (in ppt2->k); ignored unless ppt2->k_sampling==lin_k_sampling or log_k_sampling */
  int k_size_custom; /**< User-provided size of the k vector ppt2->k; ignored unless ppt2->k_sampling==lin_k_sampling or log_k_sampling */



  // ====================================================================================
  // =                                    k3 sampling                                   =
  // ====================================================================================

  int k3_size_min;        /**< Minimum number of grid points for any (k1,k2) pair,
                          used when 'k3_sampling' is set to smart */
  int k3_size;            /**< Fixed number of grid points for any (k1,k2) pair,
                          used when 'k3_sampling' is set to lin or log */
  double q_linstep_song;  /**< Upper bound on linear sampling step in k space for the
                          transfer functions, in units of one period of acoustic oscillation,
                          2*pi/(tau0-tau_rec) */



  // ====================================================================================
  // =                                   Interpolation                                  =
  // ====================================================================================

  enum interpolation_methods sources_time_interpolation; /**< Method for the time-interpolation of the line of sight sources */
  enum interpolation_methods sources_k3_interpolation;   /**< Method for the k3-interpolation of the line of sight sources */



  // ====================================================================================
	// =                                      Bessels                                     =
  // ====================================================================================

  double bessel_j_cut_song;  /* Value of j_l1(x) below which it is approximated by zero (in the region x << l) */
	double bessel_J_cut_song;	/* Value of J_Llm(x) below which it is approximated by zero (in the region x << l) */
  double bessel_x_step_song; /* Linear step dx for sampling spherical Bessel functions j_l1(x) and functions J_Llm(x) */



  // ====================================================================================
	// =                                        Misc                                      =
  // ====================================================================================

  ErrorMsg error_message;         /**< Zone for writing error messages */
  short store_transfers_to_disk;  /**< Should we store the transfer functions to disk? */
  short load_transfers_from_disk; /**< Should we load the transfer functions from disk? */
  short store_sources_to_disk;    /**< Should we store the source functions to disk? */
  short load_sources_from_disk;   /**< Should we load the source functions from disk? */
  short old_run; /**< set to _TRUE_ if the run was stored with a version of SONG smaller than 1.0 */

};  /* end of struct precision2 declaration */

#endif
