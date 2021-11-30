/** @file perturbations2.h Documented header file for the perturbations2.c module */

#ifndef __PERTURBATIONS2__
#define __PERTURBATIONS2__

#include "common2.h"
#include "input.h"
#include "perturbations.h"
#include "perturbations2_macros.h"


// ======================================================================================
// =                                      Options                                       =
// ======================================================================================

/**
 * Implementations of the tight coupling approximation.
 *
 * The TCA allows to truncate the photon hierarchy during the tight coupling regime,
 * using the fact that kappa_dot is larger than any other scale in the system.
 *
 * In SONG we use a conservative approach and turn off the TCA approximation
 * as soon as any of the three wavemodes satisfies k*tau_c >> 1.
 */
enum tca2_method {
  tca2_none,               /**< Turn off the TCA approximation */
  tca2_first_order_pitrou  /**< Use the second-order generalisation of the TCA1 scheme
                           from Pitrou 2011, http://arxiv.org/abs/1012.0546 */
};


/**
 * Implementations of the radiation-streaming approximation.
 *
 * The RSA provides an efficient way to track the density contrast and the velocity
 * of the photons and neutrinos, while neglecting the shear and the higher moments.
 *
 * In SONG, we turn the RSA on when k3 is subhorizon, k3*tau >> 1, regardless of
 * k1 and k2.
 */
enum rsa2_method {
  rsa2_none,          /**< Turn off the RSA approximation */
  rsa2_MD,            /**< Implement the RSA approximation assuming no collisions (i.e. no reionisation effects) */
  rsa2_MD_with_reio   /**< Implement the RSA approximation including collisions up to first-order in tau_c=1/kappa_dot */
};


/**
 * Implementations of the ultra relativistic fluid approximation.
 * 
 * NOT IMPLEMENTED YET!
 */
enum ufa2_method {
  ufa2_none            /**< Turn off the UFA approximation */
};


/**
 * Implementations of the no-radiation approximation.
 *
 * This is a blunt approximation where we treat photons, neutrinos and all massive
 * species as perfect fluids well after matter radiation equality.
 */
enum nra2_method {
  nra2_none,            /**< No NRA approximation */
  nra2_all,             /**< Switch off all multipoles for all relativistic species */
  nra2_fluid            /**< Switch off l>1 multipoles for all relativistic species */
};


/**
 * On and off flags for the NRA (no radiation approximation)
 */
enum nra_flags {
  nra_off,
  nra_on
};


/**
 * Possible sampling methods for the ppt2->k array.
 */
enum sources2_k_sampling {
  lin_k_sampling,                  /**< Linear k sampling */
  log_k_sampling,                  /**< Logarithmic k sampling */
  smart_sources_k_sampling         /**< Smart k-sampling, that is, logarithmic + linear */
};


/**
 * Possible sampling methods for the ppt2->k3[index_k1][index_k2] array.
 */
enum sources2_k3_sampling {

  lin_k3_sampling,         /**< Linear k3 sampling */

  log_k3_sampling,         /**< Logarithmic k3 sampling */

  smart_k3_sampling,       /**< Smart k3-sampling k, see perturb2_get_k3_list() */

  sym_k3_sampling,         /**< Adopt a symmetric sampling in (k1,k2,k3) by means of a transformation whereby
                           the three wavemodes end up sharing the same grid. This sampling allows to describe
                           the physical Fourier space as a cube despite the triangular condition. The symmetric
                           sampling is crucial to compute the power spectrum P(k) of second-order perturbations,
                           as it allows to sample the edges of the triangular condition with a logarithmic power.
                           SONG's default sampling, the 'smart' sampling, will always fail to describe P(k) at
                           some small enough value of k. In some cases (eg. the magnetic field at early times), the
                           failure happens quicker because of a specific structure of the kernel, but also for 
                           the simple matter power spectrum you need the symmetric sampling at small k. Note that
                           the symmetric sampling can be considered the Fourier-space version of the (l1,l2,l3)
                           sampling adopted in Fergusson et al 2009 (http://arxiv.org/abs/0812.3413) for the cubic
                           interpolation of the bispectrum. */
         
  theta12_k3_sampling      /**< Linear sampling in the angle between k1 and k2 */
};


/**
 * Which quadratic sources should we interpolate? Options for the 
 * perturb2_quadratic_sources_at_tau() function.
 */
enum quadratic_source_interpolation {
  interpolate_total,         /**< Interpolate the quadratic part of the full Boltzmann equation
                             (collision term + Liouville term) */
  interpolate_collision,     /**< Interpolate the quadratic part of the collision term */
  interpolate_d_total,       /**< Interpolate the conformal time derivative of the quadratic part of
                             the full Boltzmann equation (collision term + Liouville term). Usable only
                             if ppt2->compute_quadsources_derivatives==_TRUE_. */
  interpolate_d_collision,   /**< Interpolate the conformal time derivative of the quadratic part of
                             the collision term. Usable only if ppt2->compute_quadsources_derivatives==_TRUE_. */
  interpolate_dd_total,      /**< Interpolate the conformal-time second derivative of the quadratic part
                             of the full Boltzmann equation (collision term + Liouville term). Usable only
                             if ppt2->compute_quadsources_derivatives==_TRUE_. */
  interpolate_dd_collision   /**< Interpolate the conformal-time second derivative of the quadratic part
                             of the collision term. Usable only if ppt2->compute_quadsources_derivatives==_TRUE_. */
};

/**
 * Which quadratic sources should we compute? Options for the
 * perturb2_quadratic_sources() function.
 */
enum quadratic_source_computation {
  compute_total_and_collision,  /**< compute the quadratic part of the full Boltzmann equation
                                (collision term + Liouville term) in pvec_quadsources and that
                                of the collision term alone in pvec_quadcollision. */
  compute_only_liouville,       /**< compute only the quadratic part of the Liouville term in
                                pvec_quadsources */
  compute_only_collision,       /**< compute only the quadratic part of the collision term in
                                pvec_quadcollision */
  compute_only_loss_term,       /**< compute only the collisional loss term in pvec_quadcollision */ 
  compute_only_gain_term        /**< compute only the collisional gain term in pvec_quadcollision */ 
};


// ======================================================================================
// =                                     Definitions                                    =
// ======================================================================================

/**
 * Maximum number of equations evolved in the differential system.
 *
 * Feel free to increase it, it is just a memory parameter.
 */
#define _MAX_NUM_EQUATIONS_ 4096

/**
 * Exclude edges of the triangular condition on (k1,k2,k3).
 *
 * The differential system at second-order has to be solved for a set of three wavemodes
 * (k1,k2,k3) whereby k3 has to be in the range |k1-k2|<=k3<=k1+k1. When k3 is too close
 * to the boundaries, numerical instabilities might arise such as nan's or larger-than-one
 * sines and cosines. In order to avoid that, we define here the safety distance between
 * k3 and the bounds. This safety distance is going to correspond to the largest scale
 * probed by SONG. The largest angular scale probed by SONG is l_min=2, which corresponds
 * to k~l_min/tau_0^-1~1e-4. Therefore, setting _MIN_K3_DISTANCE_ to anything smaller
 * than 1e-6 is quite safe.
 */
#define _MIN_K3_DISTANCE_ 1e-10

/**
 * Exclude squeezed configurations (TODO: remove)
 *
 * The differential system dies when k3 is much smaller than k1+k2. These
 * configurations are irrelevant, so we set a minimum ratio between k1+k2
 * and k3.
 */
#define _MIN_K3_RATIO_ 100



// ======================================================================================
// =                              Perturbations structure                               =
// ======================================================================================

struct perturbs2
{

  // ------------------------------------------------------------------------------------
  // -                                    Output flags                                  -
  // ------------------------------------------------------------------------------------

  /**
   * Flags that the perturbations2.c module inherits form the input2.c module, and will be
   * read by other modules as well.
   */
  //@{

  short has_perturbations2;      /**< Do we need second-order perturbations at all? */

  short has_cmb_bispectra;       /**< Do we need to compute the second-order bispectrum of the CMB? */

  short has_cmb_spectra;         /**< Do we need to compute the second-order spectrum (C_l and Pk) of the CMB? */

  short has_cls;                 /**< Do we need to compute the second-order C_l? */

  short has_cmb_temperature;     /**< Do we need to compute spectra or bispectra for the CMB temperature? */

  short has_cmb_polarization_e;  /**< Do we need to compute spectra or bispectra for the CMB E-modes? */

  short has_cmb_polarization_b;  /**< Do we need to compute spectra or bispectra for the CMB B-modes? */

  short has_bk_delta_cdm;        /**< Do we need to compute the bispectrum delta_cdm(k1,k2,k3,tau) of the density contrast of the cold dark matter component? */

  short has_pk_delta_cdm;        /**< Do we need to compute the power spectrum of the density contrast of the cold dark matter component? */

  short has_bk_delta_b;          /**< Do we need to compute the bispectrum delta_b(k1,k2,k3,tau) of the density contrast of baryons? */

  short has_pk_delta_b;          /**< Do we need to compute the power spectrum of the density contrast of baryons? */

  short has_pk_magnetic;         /**< Do we need to compute the power spectrum of the magnetic field generated at recombination?
                                 For reference, see Fidler, Pettinari & Pitrou 2015. If you are interested in an optimal result
                                 at early times and large scales, choose the symmetric k-sampling for the sources by setting
                                 sources2_k3_sampling=sym_sampling, and make the quadsources time sampling finer by decreasing the
                                 parameter perturb_sampling_stepsize_for_quadsources (0.01 is already good) */
  //@}
  

  // ------------------------------------------------------------------------------------
  // -                                Initial conditions                                -
  // ------------------------------------------------------------------------------------

  short has_ad;                       /**< Adiabatic initial conditions */
  short has_ad_first_order;           /**< Adiabatic initial conditions, with no quadratic sources. This is a useful 
                                      debug tool when paired with quadratic_sources=no, because then SONG results
                                      must match CLASS results. See perturb2_initial_conditions() for further
                                      details. */
  short has_zero_ic;                  /**< Vanishing initial conditions, for debug purposes */  
  short has_unphysical_ic;            /**< Custom initial conditions, for debug purposes */
  double primordial_local_fnl_phi;    /**< Amount of primordial non-Gaussianity of the local type */

  int index_ic_first_order;           /**< Index of the initial conditions used at first order; it is one of the
                                      ppt->index_ic_X indices. */




  // ------------------------------------------------------------------------------------
  // -                                Differential system                               -
  // ------------------------------------------------------------------------------------
  
  short has_polarization2;                  /**< Include the polarised hierarchies in the differential sytem? */  

  short has_quadratic_liouville;            /**< Include the quadratic sources in the Liouville operator? */      
  short has_quadratic_collision;            /**< Include the quadratic sources in the photon-baryon collision term? */      

  short has_perfect_baryons;                /**< Treat baryons as a pressureless perfect fluid with no shear? */
  short has_perfect_cdm;                    /**< Treat cold dark matter as a pressureless perfect fluid with no shear? */

  short has_perturbed_recombination_stz;    /**< Include the perturbation to the fraction of free electrons when computing the collision term? */
  short perturbed_recombination_use_approx; /**< Compute the perturbed fraction of free electrons using the approximation in eq. 3.23 of
                                            Senatore et al. 2009 (http://arxiv.org/abs/0812.3652)? */

  /**
   * Should we include the quadratic sources in the differential system?
   * 
   * If you run SONG without quadratic sources, you are effectively running a first-order
   * code. Tranfer functions will be contant in k1 and k2, and will match those from CLASS.
   * In general, SONG power spectra will match CLASS ones. For this reason, setting
   * quadratic_sources=no is a great debugging tool.
   *
   * The intrinsic bispectrum, on the other hand, will match the local bispectrum, because
   * in the bispectrum module we effectively integrate three first-order transfer functions
   * together with two power spectra P(k1)*P(k2). This integral corresponds exactly to the
   * local bispectrum with f_nl = ppt2->primordial_local_fnl_phi. For more info on this
   * "local limit", see Sec. 6.5.3 of http://arxiv.org/abs/1405.2280.
   * 
   * Note that in absence of quadratic sources, the initial conditions need to have a
   * non-vanishing amplitude, lest all perturbations vanish. This is achieved by 
   * setting a nonzero value for primordial_local_fnl_phi.
  */
  short has_quadratic_sources;


  /**
   * Should we stop sampling the line of sight sources at the end of recombination?
   *
   * The CMB is sourced by recombination effects (Sachs-Wolfe, doppler, early ISW)
   * and late-time effects (reionisation and late ISW). If the user requested only
   * recombination effects, this flag will be true and the evolution of the 
   * differential system will stop shortly after recombination.
   *
   * Computing only the collisional sources is a major speed-up, because it allows
   * SONG to evolve the differential system and to compute the transfer functions
   * only up to the end of recombination.
   *
   * The time marking the end of recombination is
   *
   *   ppt2->tau_sampling[ppt2->index_tau_end_of_recombination],
   *
   * where ppt2->index_tau_end_of_recombination is computed based on the parameter
   * recombination_max_to_end_ratio, in the function perturb2_timesampling_for_sources().
   *
   * This flag is ignored if the user asked for a custom time sampling for the sources
   * (custom_time_sampling_song_sources=yes), while it is always true if the user
   * explicitly set only_recombination=yes in the parameter file.
   */
  int has_only_recombination;


  /**
   * Should we start sampling the line of sight sources at the beginning of
   * reionisation?
   *
   * This flag is useful to investigate and debug reionisation effects.
   *
   * The time marking the beginning of reionisation is
   *
   *   ppt2->tau_sampling[ppt2->index_tau_reio_start],
   *
   * This flag is ignored if the user asked for a custom time sampling for the sources
   * (custom_time_sampling_song_sources=yes).
   */
  int has_only_reionisation;


  /**
   * Debug flag: should we include only the loss part of the collision term in 
   * the second order line of sight sources?
   */
  int has_only_loss_term;


  /**
   * Debug flag: should we include only the gain part of the collision term in 
   * the second-order line of sight sources?
   */
  int has_only_gain_term;


  /**
   * Which equation should we use to evolve the curvature potential phi in Newtonian gauge?
   *
   * The current options are:
   *
   *  - poisson, to use the time-time Einstein equation (eq. 5.2 of http://arxiv.org/abs/1405.2280).
   *
   *  - longitudinal, to use the space-time Einstein equation (eq. 3.98).
   *
   *  - huang, to use a combination of the trace and time-time Einstein equations, as in Huang
   *    2012 (http://arxiv.org/abs/1201.5961).
   * 
   * In principle, we could compute phi directly from psi by using the constraint equation
   * obtained by combining the time-time and space-time Einstein equations (eq. 5.7 of
   * http://arxiv.org/abs/0812.3652). This option however turns up to be be numerically
   * unstable at early times, both at first and second order.
   */
  enum phi_equation phi_eq;


  /**
   * Should we rescale the CMB source function in view of bispectrum integration?
   * 
   * In order to compute the bispectrum integral, it is useful to rescale the line of sight
   * sources by a sin(theta_1)^(-m) factor, where theta_1 is the angle between \vec{k1} and
   * \vec{k}. The benefit from this rescaling is explained in Sec. 6.2.1.2 of
   * http://arxiv.org/abs/1405.2280.
   *
   * Note that, in SONG instead of rescaling the source function directly, we rescale the
   * quadratic sources by the same factor, with the same result.
   *
   * The rescaling is asymmetrical in k1<->k2 because it contains theta_1 but not theta_2.
   * It follows that the rescaling breaks the (k1,k2) symmetry of the transfer functions.
   * The effect of exchanging k1 with k2 in the rescaled source function, however, is very
   * simple (eq. 6.47):
   *
   * S_rescaled(k2,k1,k3) = S_rescaled(k1,k2,k3) *  (-1)^m * (k2/k1)^m
   * 
   * This follows from the identity k2 * sin(theta_2) = k1 * sin(theta_1), which is valid
   * for our choice of wavmode geometry (eq. B.4), and from the relation between the 
   * unrescaled source functions (eq. B.12):
   *
   * S(k2,k1,k3) = S(k1,k2,k3) *  (-1)^m.
   */
  short rescale_cmb_sources;


  // ------------------------------------------------------------------------------------
  // -                               Line of sight sources                              -
  // ------------------------------------------------------------------------------------

  /**
   * Array containing the source function S_lm(k1,k2,k3,tau) for all required types (T,E,B)
   *
   * The ppt2->sources array should be indexed as follows:
   *  
   *    ppt2->sources  [index_tp2]
   *                   [index_k1]
   *                   [index_k2]
   *                   [index_tau*k3_size + index_k3]
   *
   * - index_tp2 is a composite index that includes both the field (X=T,E,B...) and the
   *   multipole (l,m); it is expanded as index_tp2 = ptr2->index_tp2_X + lm(l,m).
   * - index_k1 goes from 0 to ppt2->k_size-1.
   * - index_k2 goes from 0 to ppt2->k_size-index_k1-1 due to symmetry reasons.
   * - index_k3 goes from 0 to ppt2->k3_size[index_k1][index_k2]-1.
   * - k3_size is short for ppt2->k3_size[index_k1][index_k2].
   */
  double **** sources;
  
  short * has_allocated_sources;  /**< If has_allocated_sources[index_k1]==_TRUE_,
                                  then ppt2->sources[index_k1] is fully allocated */

  /**
   * Flags for internal use in the perturbations2.c module.
   */
  //@{

  short has_source_T;          /**< Should we store in ppt2->sources the source function for the CMB temperature? */
  short has_source_E;          /**< Should we store in ppt2->sources the source function for the CMB E-polarization? */
  short has_source_B;          /**< Should we store in ppt2->sources the source function for the CMB B-polarization? */
  short has_source_delta_cdm;  /**< Should we store in ppt2->sources the density contrast of cold dark matter? */
  short has_source_M;          /**< Should we store in ppt2->sources the source function for the magnetic field generated
                               at recombination? Includes only the m=1 dipole. For reference, see Fidler, Pettinari &
                               Pitrou 2015. */

  short has_cmb;               /**< Do we need CMB-related sources at all? (e.g. photon temperature or polarisation) */
  short has_lss;               /**< Do we need sources related to the large scale structure? (e.g. lensing potential) ? */  


  /* - Collisional sources */

  short has_pure_scattering_in_los;   /**< Include the purely second-order scattering terms in the line-of-sight sources? */
  short has_quad_scattering_in_los;   /**< Include the quadratic scattering terms in the line-of-sight sources? */
  

  /* - Metric sources */

  short has_pure_metric_in_los;     /**< Shall we include the purely second-order metric terms in the line-of-sight sources,
                                    as they appear in the Boltzmann equation (eq. 5.112 of http://arxiv.org/abs/1405.2280)?
                                    If this flag is turned on, no integration by parts will be performed, and the SW and ISW
                                    effects will be both included. If either of the SW and ISW flags are turned on, this flag is 
                                    overridden. */

  short has_quad_metric_in_los;     /**< Shall we include the terms quadratic in the metric in the line-of-sight sources? The
                                    exact terms that are included vary depending on whether we are using exponential potentials;
                                    compare eq. 4.97 and 4.100 of http://arxiv.org/abs/1405.2280. */

  short has_sw;                     /**< Include the Sachs-Wolfe term in the line-of-sight sources for the photon monopole?
                                    The SW term has the form g*psi, where g is the visibility function and psi is the time
                                    potential in Newtonian gauge. The SW term arises together with the ISW term after
                                    integrating by parts the 4*k*psi contribution to the photon dipole. For this reason,
                                    the flag ppt2->has_sw overrides the flag ppt2->has_pure_metric_in_los, lest SONG double
                                    counts the 4*k*psi contribution. */

  short has_isw;                    /**< Include the integrated Sachs-Wolfe term in the line-of-sight sources for the photon
                                    monopole? The ISW term has the form e^-kappa*(phi'+psi'), where kappa is the opacity and
                                    phi and psi are the curvature and time potentials in Newtonian gauge, respectively. The ISW
                                    term arises together with the SW term after integrating by parts the 4*k*psi contribution
                                    to the photon dipole. For this reason, the flag ppt2->has_isw overrides the flag
                                    ppt2->has_pure_metric_in_los, lest SONG double counts the 4*k*psi contribution. */

  short only_early_isw;             /**< Should we neglect the contribution to the integrated Sachs-Wolfe effect from late
                                    times? */

  short use_exponential_potentials; /**< Should we use exponential potentials in constructing the line of sight sources?
                                    The exponential potentials are defined as psi_e=psi(1-psi) and phi_e=phi(1+phi);
                                    they differ from the regular potentials by a quadratic contribution (see eq. 3.21
                                    of http://arxiv.org/abs/1405.2280). If you include all metric sources in the source
                                    function (that is, if you set has_pure_metric_in_los=yes and has_quad_metric_in_los=yes),
                                    using exponential potentials will not affect the result. Turn this flag on, together
                                    with the SW and ISW effects, to match the analytical approximation for the squeezed
                                    bispectrum for small l. */


  /* - Propagation sources */

  short has_time_delay_in_los;      /**< Include the time-delay term for photons in the line-of-sight sources? */

  short has_redshift_in_los;        /**< Include the redshift term for photons in the line-of-sight sources? This is an obsolete
                                    flag, as we now use the delta_tilde transformation to deal with the redshift term; refer 
                                    to Sec. 5.5.3 of http://arxiv.org/abs/1405.2280 for more details on the transformation. */

  short has_lensing_in_los;         /**< Include the lensing term for photons in the line-of-sight sources? */

  short use_delta_tilde_in_los;     /**< Shall we implement the delta_tilde transformation to deal with the redshift term?
                                    Refer to Sec. 5.5.3 of http://arxiv.org/abs/1405.2280 for details on the transformation. */
  //@}


	short use_test_source;   	/**< If true, use an hard-coded test source for debug purposes */

  int index_tp2_T;             /**< Beginning of the photon intensity hierarchy in the ppt2->sources array */
  int index_tp2_E;             /**< Beginning of the photon E-mode hierarchy in the ppt2->sources array */
  int index_tp2_B;             /**< Beginning of the photon B-mode hierarchy in the ppt2->sources array */
  int index_tp2_delta_cdm;     /**< Index for the second-order density contrast of cold dark matter in the ppt2->sources array */
  int index_tp2_M;             /**< Index for the magnetic field source in the ppt2->sources array */

  int n_sources_T;           /**< Number of sources to be computed for photon temperature */
  int n_sources_E;           /**< Number of sources to be computed for photon E-polarization */
  int n_sources_B;           /**< Number of sources to be computed for photon B-polarization */

  int n_nonzero_sources_E;   /**< Number of nonzero sources to be computed for photon E-polarization; this
                             is basically ppt2->n_sources_E minus the l=0 and l=1 modes */

  int n_nonzero_sources_B;   /**< Number of nonzero sources to be computed for photon B-polarization; this
                             is basically ppt2->n_sources_B minus the l=0, l=1 and m=0 modes */

  int tp2_size;              /**< Number of source types that we need to compute */

  char (*tp2_labels)[_MAX_LENGTH_LABEL_];  /**< Labels of the various source types. For example, ppt2->tp2_labels[index_tp2_phi]
                                           is equal to the string "phi". Useful for printing out results. */


  
  // ------------------------------------------------------------------------------------
  // -                                     CMB fields                                   -
  // ------------------------------------------------------------------------------------

  /**
   * Field and parity of the CMB perturbations.
   *
   * The temperature field has different spin (S=0) than the polarisation fields (S=2).
   * Similarly, the temperature and E fields have different parity (even) than the B
   * field (odd).
   *
   * These disparities translate into slightly different geometrical properties. An
   * important example is the harmonic decomposition of the product of two perturbations,
   * which is crucial for the delta_tilde transformation.
   *
   * Here we define indices to keep track of field type (T, E, B, Rayleigh...) and of its
   * parity (even or odd). 
   */
  //@{  
  int index_pf_t; /**< Index denoting the temperature CMB field */
  int index_pf_e; /**< Index denoting the E-mode polarisation CMB field */
  int index_pf_b; /**< Index denoting the B-mode polarisation CMB field */
  int pf_size;    /**< Number of fields (T,E,B,...) considered in the current run */

  int field_parity[_MAX_NUM_FIELDS_];  /**< ppt2->field_parity[ppt2->index_pf_x] with x=t,e,b is the parity of the
                                       field x; it is either EVEN or ODD. The size of the array is ppt2->pf_size. */

  char (*pf_labels)[_MAX_LENGTH_LABEL_]; /**< Labels of the fields; for example, if temperature is requested,
                                         ppt2->pf_labels[ppt2->index_pf_t] is "t" */
  //@}

  
  // ------------------------------------------------------------------------------------
  // -                                     Multipoles                                   -
  // ------------------------------------------------------------------------------------

  /** 
   * Array containing the values of the azimuthal index m for which we will solve the
   * system and compute the source function.
   *
   * SONG can compute any m-mode, not only scalar scalar (m=0), vector (m=1) and tensor
   * (m=2) modes. Note however that CMB perturbations will be progressively suppressed
   * as m increases because of tight coupling.
   *
   * Certain observables are sourced only by non-scalar modes, for example, the B-modes
   * of photon polarisation, the velocity vorticity, and the magnetic fields.
   *
   * The ppt2->m array cannot contain negative values. SONG will include the negative 
   * contributions automatically using symmetry relations, e.g. T_lm = (-1)^m T_l-m.
   * 
   * The ppt2->m array is just a copy of ppr2->m, which in turn is a copy of the parameter
   * modes_song in the ini file.
   *
   */
  int * m;
  int m_size; /**< Size of the ppt2->m array */

  /** 
   * Number of additional multipoles l required to compute the quadratic sources in
   * the Boltzmann equation.
   * 
   * The Boltzmann equation couples neighbouring moments, so that the time derivative
   * of the Boltzmann moment Delta_l depends on Delta_(l-1) and Delta_(l+1).
   *
   * The coupling holds also for the quadratic sources, with the difference that it
   * can extend farther than the neighbouring moments.
   *
   * The parameter ppt2->lm_extra encodes how far the coupling extends: the time
   * derivative of the second-order Delta_l depends on the first-order perturbations
   * with Delta_(l-lm_extra) and Delta_(l+lm_extra).
   * 
   * The value of lm_extra is a measure of the angular complexity of the quadratic
   * sources. It is gauge dependent: for the Newtonian gauge lm_extra=1 while for
   * the synchronous gauge lm_extra=3. This can be seen from eq. 3.29 of Naruko
   * et al. 2013 (http://arxiv.org/abs/1304.6929): the Newtonian gauge has a simple
   * spatial metric (h_ij) while the synchronous gauge, with its h_ij~k_i*k_j
   * dependence, gives rise to terms quartic in n (the photon direction). These
   * quartic terms, once expanded in multipoles space, are what gives the larger
   * ppt2->lm_extra for the synchronous gauge.
   * 
   * We use ppt2->lm_extra to determine the size of the rotation vectors
   * ppw2->rotation_1 and ppw2->rotation_2.
   */
  int lm_extra;
  
  int largest_l;   /**< The largest multipole l for a second-order perturbation in this module */
  
  int largest_l_quad; /**< The largest multipole l for a first-order perturbation in this module, given
                      by ppt2->largest_l+ppt2->lm_extra; see documentation for ppt2->lm_extra.  */

  /**
   * Array used to access the second-order massless hierarchies in the perturbations2.c module.
   *
   * In SONG, we decompose the angular dependence of the perturbations into spherical harmonics.
   * This allows to treat the Boltzmann equation numerically by expressing it as a hierarchy of
   * equations that depend on the harmonic indices (l,m). The perturbations themselves form a
   * hierarchy of moments, each with its own (l,m) indices.
   * 
   * To write this information in terms of C arrays, we need to create a correspondence between the
   * (l,m) indices and a sequence of integers starting from zero. This is the same problem of
   * compressing a square matrix in a 1D array, with the difference that here the matrix is ragged,
   * because the m direction depends on the considered l.
   * 
   * The multipole index l can assumes all the values between 0 and l_max, the maximum multipole
   * required in this module. The azimuthal index m, on the other hand, has two constraints:
   * it needs to be smaller or equal than l, and it must be contained in ppt2->m.
   * 
   * We encode all this information in the array ppt2->lm_array[l][index_m]. To ease the notation,
   * we define a shorthand, lm(l,m). This will be used ubiquitously in perturbations2.c to access
   * the perturbations in many arrays; for example, the (l,m) multipole of the photon hierarchy is
   * acessed as:
   * 
   *   ppw2->pv->y[ppw2->pv->index_pt2_monopole_g + lm(l,m)]
   *   ppw2->pv->dy[ppw2->pv->index_pt2_monopole_g + lm(l,m)]
   * 
   * while the line of sight source for the (l,m) multipole of the E-modes is in
   * 
   *   ppt2->sources[ppt2>index_tp2_E + lm(l,m)].
   * 
   * Be careful accessing lm(l,m). If either l or m do not satisfy the constraints (0<=l<=l_max,
   * m<=l, m in ppt2->m), you will encounter segmentation faults or unexpected behaviour. The
   * lm(l,m) macro should be used only inside constrained (l,m) loops or inside if blocks.
   */
  int ** lm_array;


  /**
   * Array used to access the first-order massless hierarchies in the perturbations2.c
   * module.
   *
   * We define a new array rather than recycling ppt2->lm_array because the first-order
   * hierarchies need to extend further than the second-order ones.
   * 
   * Used to address the ppw2->rotation_1 and ppw2->rotation_2 arrays.
   *
   * See documentation for ppt2->lm_array for further details.
   */
  int ** lm_array_quad;  


  /**
   * Array used to access the massive hierarchies in the perturbations2.c module.
   *
   * The massive species (baryons and cold dark matter) are expanded both in velocity moments
   * and in spherical harmonics. Therefore, their hierarchies have a velocity index n in addition
   * to the multipole indices l and m.
   * 
   * To index the massive hierarchies, we use the array ppt2->nlm_array[n][l][index_m], with
   * nlm(n,l,m) as a shorthand. For example, the (n,l,m) multipole of the baryon hierarchy is
   * acessed as:
   * 
   *   ppw2->pv->y[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)]
   *   ppw2->pv->dy[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)].
   * 
   * As for now, n is capped to n=2, or n=1 when using the perfect fluid approximation.
   * 
   * See documentation for ppt2->lm_array for further details.
   *
   * For more details on the velocity expansion of the massive species, please refer to
   * Sec. 5.3.1.3 of http://arxiv.org/abs/1405.2280.
   */
  int *** nlm_array;


  /**
   * Coupling coefficients C and D in harmonic space.
   * 
   * The coupling coefficients link multipoles with adjacent l-values in the Boltzmann
   * equation.
   * 
   * They are a convenient rearrangement of the 3j symbols that arise from the spherical
   * harmonics expansion of products in the Boltzmann equation; they are present both
   * at first and second order.
   *
   * The coupling coefficients C and D are usually accessed via the shorthands
   * defined in perturbations2_macros.h.
   * 
   * Refer to Sec. A.4.1 of http://arxiv.org/abs/1405.2280 and eq. 141 of Beneke &
   * Fidler 2010 for further details on C and D.
   */
  //@{
  double ** c_minus;
  double ** c_plus;
  double ** d_minus;
  double ** d_plus;
  double ** d_zero;
  //@}
  

  /**
   * General coupling coefficients in harmonic space.
   *
   * In harmonic space, the product of two real-space perturbations X*Y reads as a weighted
   * sum over X(l1,m1) and Y(l2,m2). We refer to the weights in such sum as the coupling
   * coefficients.
   * 
   * If one of the two perturbations is a dipole (X_lm=0 unless l=1) the coupling coefficients
   * are just C and D, which we have defined above. If instead X and Y are two arbitrary photon
   * perturbations, the coupling coefficients take the following form:
   * 
   *       prefactor * (-1)^m * (2*l+1) * ( l1  l2  l ) * (  l1  l2   l )
   *                                      ( 0   F  -F )   (  m1  m2  -m )
   * 
   * where F=2 polarisation and F=0 otherwise, and the prefactor is shown below.
   * 
   * The coupling coefficients will be crucial to compute the spherical decomposition of
   * the delta-squared term in the delta_tilde transformation. They are indexed as:
   *
   *   ppt2->coupling_coefficients[ppt2->index_pf_t]
   *                              [lm(l,m)]
   *                              [l1]
   *                              [m1+ppt2->l1_max]
   *                              [l2]
   * 
   * and they vanish for configurations where the triangular inequality between l, l1 and l2
   * is not met.
   *
   * FORMAL DEFINITION
   * 
   * In general, a function P[ab] in helicity space, can be decomposed in multipole space
   * into its intensity (I), E-mode polarisation (E) and B-mode polarisation (B) parts according
   * to eq. 2.10 of arXiv:1401.3296:
   * 
   *   P[-+](l,m) = P_E(l,m) - i P_B(l,m)
   *   P[+-](l,m) = P_E(l,m) + i P_B(l,m)
   *   1/2 * (P[++] + P[--])(l,m)  =  P_I(l,m)
   * 
   * If P[ab] is a product in helicity space, 
   * 
   *   P[ab] = 1/2 (X[ac]*Y[cb] + Y[ac]*X[cb]),
   * 
   * then its I, E and B components can be inferred, respectively, by eq. 3.9, 3.7 and 3.6 of
   * the same reference:
   * 
   *   P_I(l,m) ->  + 0.5 * i^L * ( l1   l2 | l ) * (  l1   l2  | l )  
   *                              ( 0    0  | 0 )   (  m1   m2  | m ) 
   *                      * [ X_I(l1,m1) Y_I(l2,m2) + Y_I(l1,m1) X_I(l2,m2) ] 
   * 
   *   P_E(l,m) ->  + 1 * i^L * ( l1   l2 | l ) * (  l1   l2  | l )  
   *                            ( 0    2  | 2 )   (  m1   m2  | m ) 
   *                    * [ X_I(l1,m1) Y_E(l2,m2) + Y_I(l1,m1) X_E(l2,m2) ]  for L even, 0 otherwise
   * 
   *   P_B(l,m) ->  + 1 * i^(L-1) * ( l1   l2 | l  ) * (  l1   l2  | l )  
   *                                ( 0    2  | 2  )   (  m1   m2  | m ) 
   *                    * [ X_I(l1,m1) Y_E(l2,m2) + Y_I(l1,m1) X_E(l2,m2) ]  for L odd, 0 otherwise
   * 
   * where: 
   * 
   *  - We are neglecting the first-order B-modes.
   *  - The symbol (    | ) denotes a Clebsch-Gordan symbol.
   *  - L = l-l1-l2.
   *  - The indices (l1,l2,m1,m2) are summed.
   *  - The indices (l,m) are free.
   *  - The sums in P_I and P_E vanish for odd values of L because the intensity and E-mode fields
   *    have even parity. The B-mode field instead has odd parity. This property also ensures
   *    that both the i^L and i^(L-1) factors are real-valued.
   *  - With respect to arXiv:1401.3296 we have a -2 factor; this counters the fact that the expressions
   *    in that reference were for the quadratic term -1/2*delta*delta rather than for a generic X*Y.
   *    
   * In what follows, we compute the general coupling coefficients given by
   * 
   *       prefactor * ( l1   l2 |  l  ) * (  l1   l2 |  l  )
   *                   ( 0    F  |  F  )   (  m1   m2 |  m  ) ,
   * 
   * for F=0 (intensity) or F=2 (E and B-mode polarisation). We store the result in the array
   * ppt2->coupling_coefficients[index_pf][lm][l1][m1-m1_min][l2], where index_pf refers to
   * the considered field (I,E,B...). Note that, in terms of 3j symbols, the coefficients read: 
   *       
   *       prefactor * (-1)^m * (2*l+1) * ( l1   l2   l  ) * (  l1   l2   l  )
   *                                      ( 0    F   -F  )   (  m1   m2  -m  )
   * 
   * The prefactor for the I,E,B fields can be read from the equations shown above:
   * 
   *       prefactor(I) = +1 * i^L
   *       prefactor(E) = +2 * i^L
   *       prefactor(B) = +2 * i^(L-1)
   * 
   * where L=l-l1-l2. For the even-parity fields (I and E) L is always even, while
   * for the odd-parity ones (B) it is odd. Therefore, the prefactors are always
   * real-valued. The polarisation prefactors are 2 insteaed of 1 because they get
   * contributions from both E and B-modes contribute. 
   * 
   * These coefficients are needed for the delta_tilde transformation that absorbs the
   * redshift term of Boltzmann equation in the polarised case (see Sec. 3.10 of
   * arXiv:1401.3296 and Huang and Vernizzi 2013a). The coefficients always multiply two
   * first-order perturbations (Eqs. 3.6, 3.7 and 3.9 of arXiv:1401.3296):
   *
   *   - delta_I(l1)*delta_I(l2) in the case of the temperature;
   *   - delta_I(l1)*delta_E(l2) in the case of the both E and B-modes.
   *
   * Note that there is no delta_B because we neglect the first-order B-modes.
   *
   */
  double ***** coupling_coefficients;

  int l1_max;   /**< Maximum l1 for which we have stored the coupling coefficients in ppt2->coupling_coefficients.
                It is determined by ppr2->l_max_los_quadratic and is used the truncate the summation for the
                delta_tilde transformation. */

  int l2_max;   /**< Maximum l2 for which we have stored the coupling coefficients in ppt2->coupling_coefficients.
                It is determined by ppr2->l_max_los_quadratic and is used the truncate the summation for the
                delta_tilde transformation. */
  

  // ------------------------------------------------------------------------------------
  // -                                    k-sampling                                    -
  // ------------------------------------------------------------------------------------

  /**
   * Variables linked to the sampling in Fourier space (k1,k2,k3) of the line of sight
   * sources.
   *
   * At second order, the transfer functions depend on three Fourier wavemodes rather than one,
   * due to mode coupling. As a result, we solve the Boltzmann-Einstein system on a 3D grid
   * in (k1,k2,k), where k1 and k2 are the magnitudes of the dummy wavevectors and k (also
   * called k3 in the code) is the magnitude of their sum. 
   *
   * The dummy waveveoctors k1 and k2 share a common k-sampling (ppt2->k) while k is sampled
   * differently for each (k1,k2) pair (ppt2->k3[index_k1][index_k3]). Given a k1 value, we
   * compute the perturbations only for k2<=k1, and obtain the remaining values using symmetry
   * relations.
   */
  //@{
  enum sources2_k_sampling k_sampling;  /**< What method should we use to determine the sampling of k1 and k2? That is, to
                                        fill ppt2->k? For the available options, see documentation of sources2_k_sampling. */

  double * k;      /**< Array containing the sampling in the k1 and k2 wavemodes, with size ppt2->k_size. Note that, for a
                   given value of k1, k2 will be computed only for those values of ppt2->k that are larger than k1. */

  int k_size;      /**< Size of ppt2->k */
                                                                  
  double k_min;    /**< Minimum k that will be ever used in SONG */
  double k_max;    /**< Maximum k that will be ever used in SONG */
                                                                  
  enum sources2_k3_sampling k3_sampling;   /**< What method should we use to determine the sampling of k3? That is, to fill
                                           ppt2->k3[index_k1][index_k2]? For the available options, see documentation of
                                           sources2_k3_sampling. */

  double *** k3;   /**< Given a (k1,k2) pair from ppt2->k, ppt2->k3[index_k1][index_k2] is 
                   the corresponding k3-grid, with size ppt2->k3_size[index_k1][index_k2]. */

  int ** k3_size;  /**< Given a (k1,k2) pair from ppt2->k, ppt2->k3_size[index_k1][index_k2] is 
                   the size of the corresponding k3-grid. */
  
  int ** index_k3_min;   /**< Given a (k1,k2) pair from ppt2->k, ppt2->index_k3_min[index_k1][index_k2] is the
                         index in ppt2->k at the left of the first element in ppt2->k3[index_k1][index_k2]. It 
                         is allocated only when the k3 sampling is set to smart. */
    
  double k_max_for_pk;  /**< Maximum value of k in 1/Mpc in P(k) */
  //@}


  // -----------------------------------------------------------------------------------
  // -                                 Time sampling                                   -
  // -----------------------------------------------------------------------------------

  /**
   * Variables linked to the sampling in conformal time (tau) of the line of sight
   * sources and of the quadratic sources.
   *
   * In SONG we deal with two time samplings:
   *
   * - ppt2->tau_sampling, for the line-of-sight sources in ppt2->sources.
   * - ppt->tau_sampling_quadsources, for the first-order perturbations in ppt->quadsources
   *   and for the quadratic sources in ppw2->quadsources_table.
   *
   * The time sampling for the line-of-sight sources needs to capture the time evolution of
   * all the physical effects relevant to the CMB: Compton collisions at recombination and
   * reionisation, Sachs-Wolfe effect, integrated Sachs-Wolfe effect, and propagation effects
   * such as time-delay and lensing. Each has a different timescale of evolution, which we
   * capture in ppt2->tau_sampling using CLASS algorithm, described in Sec. 5.3.2.3 of
   * http://arxiv.org/abs/1405.2280. Note that ppt2->tau_sampling will also be used as the
   * integration grid for the line of sight integral in the transfer2.c module.
   *
   * The quadratic sources act as sources for the second-order differential system. They
   * follow the natural timescale of the differential system, which is given by the
   * Hubble conformal time, tau_h=1/(a*H). Since this timescale is proportional to the 
   * conformal time tau, the timesampling for the quadratic sources is basically a logarithmic
   * sampling in tau.
   * 
   * Both the sources and quadsources time samplings are computed in
   * perturb2_timesampling_for_sources().
   */
  //@{

  double * tau_sampling; /**< Array with the time values where the line-of-sight sources will be computed */
  int tau_size; /**< Size of ppt2->tau_sampling */

  int index_tau_end_of_recombination; /**< Index in ppt2->tau_sampling that marks the end of recombination.
                                      Defined only if has_only_recombination==_TRUE_. */

  double z_end_of_recombination; /**< Redshift that marks the end of recombination. Defined only if has_only_recombination==_TRUE_. */

  int index_tau_rec; /**< Index in ppt2->tau_sampling where the visibility function peaks at recombination */

  int index_tau_reio_start; /**< Index in ppt2->tau_sampling where reionisation starts */

  /** Value of g/g(tau_rec) when to stop sampling the line of sight sources, where g is the
  visibility function. For example, if set to 100, then the last conformal time used
  to sample the sources will satisfy g(tau)=100*g(tau_rec). This parameter is ignored
  when the user asks for ISW or other late-time effects, because in that case the time
  sampling needs to go all the way to today. */
  double recombination_max_to_end_ratio;

  short  has_custom_timesampling;     /**< Should we use use the user-given parameters, rather than SONG algorithm,
                                      to compute the time sampling of the line-of-sight sources? Useful for debugging. */
  double custom_tau_ini;                          /**< Initial time for the time sampling of the line of sight sources
                                                  (used only if ppt2->has_custom_timesampling==_TRUE_) */
  double custom_tau_end;                          /**< Final time for the time sampling of the line of sight sources
                                                  (used only if ppt2->has_custom_timesampling==_TRUE_) */
  int    custom_tau_size;                         /**< Number of points in the time sampling of the line of sight sources
                                                  (used only if ppt2->has_custom_timesampling==_TRUE_) */
  enum   sources_tau_samplings custom_tau_mode;   /**< Sampling method for the time grid of the line of sight sources
                                                  (used only if ppt2->has_custom_timesampling==_TRUE_). See documentation
                                                  for sources_tau_samplings in perturbations.c. */

  //@}


  // ------------------------------------------------------------------------------------
  // -                                 Approximations                                   -
  // ------------------------------------------------------------------------------------
    
  /**
   * Flags and parameters related to the approximations at second order.
   * 
   * In certain limits, the Einstein-Boltzmann differential system can be simplified
   * and the hierarchies truncated, thus reducing the computational time needed to
   * solve it.
   *
   * CLASS adopts a number of such approximations. The idea is to split the integration
   * range of the differential system into as many intervals as the number of active
   * approximations. To each time interval, it corresponds a different number of equations
   * to be evolved, depending on the approximation for that interval.
   *
   * In SONG, we adopt the same mechanism. So far, we have implemented the following
   * approximations:
   *
   * - The tight coupling approximation (TCA), which allows to truncate the photon
   *   hierarchies before recombination. To implement it at second order, we generalise
   *   the formalism of Pitrou 2011 (http://arxiv.org/abs/1012.0546).
   *
   * - The radiation streaming approximation (RSA), which allows to truncate the
   *   neutrino and photon hiearchies for subhorizon modes. To implement it at
   *   second order, we generalise the formalism of Blas, Lesgourgues and Tram 2011
   *   (Sec. 4 of http://arxiv.org/abs/1104.2933).
   *
   * - The no-radiation approximation (NRA), a poor man's version of the RSA,
   *   whereby we treat photons and neutrinos as perfect fluids with no shear after
   *   the matter radiation equality.
   */
  //@{

  int tight_coupling_approximation;    /**< Which scheme to adopt for the tight coupling approximation at second order?
                                       See documentation for tca2_method for details. */
  
  double tight_coupling_trigger_tau_c_over_tau_h;    /**< Switch off the tight-coupling approximation when tau_c/tau_H is larger
                                                     than this parameter. The condition tau_c/tau_H<<1 is exactly the condition
                                                     whereby photons and electrons are tightly coupled. */

  double tight_coupling_trigger_tau_c_over_tau_k;    /**< Switch off tight-coupling approximation when tau_c/tau_k is smaller than
                                                     this parameter. The idea is that tight coupling does not hold for scales smaller
                                                     than the mean interaction distance. Note that we adopt a conservative approach and
                                                     choose k=MAX(k1,k2,k3). In the future we might relax this constraint and choose k=k3,
                                                     as we have done for the RSA (TODO)*/

  int radiation_streaming_approximation;  /**< Which scheme to adopt for the radiation streaming approximation at second order?
                                          See documentation for rsa2_method for details. */
                                        
  double radiation_streaming_trigger_tau_over_tau_k; /**< Switch on the radiation-streaming approximation when k3*tau is larger than
                                                     this parameter. The idea is that for subhorizon modes we can neglect 
                                                     all photon and neutrino moments higher than the dipole. */

  double radiation_streaming_trigger_tau_c_over_tau; /**< Switch on the radiation-streaming approximation when tau_c/tau is larger
                                                     than this parameter. The idea is that we should turn on the approximation well
                                                     after recombination (and well into the matter era...) */

  int no_radiation_approximation;  /**< Which scheme to adopt for the no-radiation approximation at second order?
                                   See documentation for nra2_method for details. */
    
  double no_radiation_approximation_rho_m_over_rho_r; /**< Switch on the no-radiation approximation when rho_matter/rho_radiation is
                                                      larger than this parameter. The idea is that we should wait until radiation
                                                      becomes irrelevant before messing with the radiation species. */

  int ur_fluid_approximation;               /**< NOT IMPLEMENTED YET */
  double ur_fluid_trigger_tau_over_tau_k;   /**< NOT IMPLEMENTED YET */


  //@}


  // ------------------------------------------------------------------------------------
  // -                                 Disk storage                                     -
  // ------------------------------------------------------------------------------------

  char sources_dir[_FILENAMESIZE_]; /**< Directory containing the line-of-sight sources. If it already exists,
                                    and ppr2->load_sources_from_disk==_TRUE_, the sources will be read from this
                                    folder into the array ppt2->sources. If it does not exist, and
                                    ppr2->store_sources_to_disk==_TRUE_, the sources will be first computed and then
                                    written to this folder from the array ppt2->sources. Either way, the directory contains
                                    one binary file for each value of k1, for a total of ppt2->k_size files. The file
                                    corresponding to index_k1 is located at ppt2->sources_paths[index_k1]; its stream
                                    is in ppt2->sources_files[index_k1]. */

  char ** sources_paths; /**< sources_paths[index_k1] is the path to the file with the line-of-sight sources
                         for k1=ppt2->k[index_k1]. Used only if ppr2->store_sources_to_disk==_TRUE_ or
                         ppr2->load_sources_from_disk==_TRUE_. */

  FILE ** sources_files; /**< sources_paths[index_k1] is the pointer to the file with the line-of-sight sources
                         for k1=ppt2->k[index_k1]. Used only if ppr2->store_sources_to_disk==_TRUE_ or
                         ppr2->load_sources_from_disk==_TRUE_. */

  FILE * sources_status_file;                 /**< NOT IMPLEMENTED YET */
  char sources_status_path[_FILENAMESIZE_];   /**< NOT IMPLEMENTED YET */


  // ------------------------------------------------------------------------------------
  // -                               Debug parameters                                   -
  // ------------------------------------------------------------------------------------

  ErrorMsg error_message;             /**< String where to write error messages */
  short perturbations2_verbose;       /**< Flag regulating the amount of information sent to standard output (none if set to zero) */

  long int count_allocated_sources;   /**< Number of allocated entries in ppt2->sources */
  long int count_memorised_sources;   /**< Number of used entries of ppt2->sources */
  
  long int count_k_configurations;    /**< Number of k-modes for which we shall solve the differential system */

  short stop_at_perturbations1;    /**< If _TRUE_, SONG will stop execution after having run the perturbations.c 
                                      module. Useful to debug the first-order transfer functions at recombination. */
  short stop_at_perturbations2;    /**< If _TRUE_, SONG will stop execution after having run the perturbations2.c
                                      module. Useful to debug the second-order transfer functions at recombination. */

  short compute_quadsources_derivatives; /**< Should we compute the first, third and fourth derivatives of the quadratic
                                         sources? Useful for debugging the RSA approximation. */

  /**
   * Parameters related to the creation of output files
   */
  //@{

  int file_verbose; /**< How much information should we include in the perturbations output files? */

  /* - k output files */

  int k_out_size; /**< Number of (k1,k2,k3) triplets where to output the perturbations (default=0) */

  int k_index_out_size;  /**< Number of (index_k1,index_k2,index_k3) triplets where to output the
                         perturbations (default=0) */
  double * k1_out; /**< List of k1 values where perturbation output is requested,
                   with size k_out_size; filled in input2.c */
  double * k2_out; /**< List of k2 values where perturbation output is requested,
                   with size k_out_size; filled in input2.c */
  double * k3_out; /**< List of k3 values where perturbation output is requested,
                   with size k_out_size; filled in input2.c */
  int * k1_index_out; /**< List of k1 indices in ppt2->k where perturbation output is requested,
                      with size k_index_out_size; filled in input2 */
  int * k2_index_out; /**< List of k2 indices in ppt2->k where perturbation output is requested,
                      with size k_index_out_size; filled in input2.c */
  int * k3_index_out; /**< List of k3 indices in ppt2->k3[index_k1][index_k2] where perturbation
                      output is requested, with size k_index_out_size; filled in input2.c */
  int * index_k1_out; /**< index_k1_out[index_k_output] is the index in ppt2->k corresponding to
                      k1=k1_out[index_k_output]; filled in perturbations2.c */
  int * index_k2_out; /**< index_k2_out[index_k_output] is the index in ppt2->k corresponding to
                      k2=k2_out[index_k_output]; filled in perturbations2.c */
  int * index_k3_out; /**< index_k3_out[index_k_output] is the index in ppt2->k[index_k1][index_k2]
                      corresponding to k3=k3_out[index_k_output]. If k3 does not satisfy the triangular
                      condition, its location in the array will contain -1. Filled in perturbations2.c */
  char (*k_out_paths)[_FILENAMESIZE_]; /**< Path of the ASCII files that will contain the perturbations as a function
                                       of tau at the desired (k1,k2,k3) values; filled in the input2.c module */
  FILE * *k_out_files; /**< ASCII file that will contain the perturbations as a function
                       of tau at the desired (k1,k2,k3) values; filled in the input2.c module */
  char (*k_out_paths_sources)[_FILENAMESIZE_]; /**< Path of the binary files that will contain the source function as a function
                                               of (k3,tau) at the desired (k1,k2) values; filled in the input2.c module */
  FILE * * k_out_files_sources; /**< Binary files that will contain the sources as a function
                               of (k3,tau) at the desired (k1,k2) values; filled in the input2.c module */
  int * k_out_data_byte; /**< k_out_data_byte[index_k_tau] is the location in the binary file 
                         k_out_files_sources[index_k_tau] of the first data block with the elements
                         of ppt2->sources. In practice, this tells the functions that want to
                         access the output files where the actual data starts, ignoring the 
                         preceding accessory data (eg. tau grid, k3 grid) */
  short output_class_perturbations; /**< If _TRUE_, output the first-order perturbations for all the values contained in ppt2->k1_out
                                    and ppt2->k2_out */
  short output_quadratic_sources;   /**< If _TRUE_, output the quadratic sources of the second-order differential system for all
                                    the values contained in ppt2->k1_out and ppt2->k2_out */

  char (*k_out_paths_quad)[_FILENAMESIZE_]; /**< Path of the ASCII files that will contain the quadratic sources as a function
                                            of tau at the desired (k1,k2,k3) values; filled in the input2.c module */
  FILE * *k_out_files_quad; /**< ASCII file that will contain the quadratic sources as a function
                            of tau at the desired (k1,k2,k3) values; filled in the input2.c module */
    
  char k_out_swap_message[_MAX_INFO_SIZE_]; /**< Message to be printed to the output files when the user asks for a (k1,k2) pair with k2>k1 */
  
  int * k_out_was_swapped; /**< Logical array to keep track whether a k_out configuration had k1 and k2 swapped */
  

  /* - tau output files */

  int tau_out_size; /**< Number of tau values where to output the perturbations (default=0) */
  int z_out_size;   /**< Number of z values where to output the perturbations (default=0) */
  double * tau_out; /**< List of tau values where perturbation output is requested,
                    with size tau_out_size; filled in input2.c. V */
  double * z_out; /**< List of z values where perturbation output is requested,
                  with size z_out_size; filled in input2.c */
  int * index_tau_out; /**< index_tau_out[index_tau_output] is the index in ppt2->tau_sampling corresponding to
                       tau=tau_out[index_tau_output]; filled in perturbations2.c */
  char (**tau_out_paths)[_FILENAMESIZE_]; /**< Path of the ASCII files that will contain the perturbations as a function
                                          of k3 at the desired (k1,k2,tau) values; filled in the input2.c module and
                                          indexed as [index_k_out][index_tau_out]. */
  FILE * (**tau_out_files); /**< ASCII file that will contain the perturbations as a function
                            of k3 at the desired (k1,k2,tau) values; filled in the input2.c module
                            and indexed as [index_k_out][index_tau_out] */
  char (*tau_out_paths_sources)[_FILENAMESIZE_]; /**< Path of the binary files that will contain the source function as a function
                                                 of (k1,k2,k3) at the desired tau values; filled in the input2.c module */
  FILE * *tau_out_files_sources; /**< Binary files that will contain the sources as a function
                                 of (k1,k2,k3) at the desired tau values; filled in the input2.c module */

  char tau_out_reduction_message[_MAX_INFO_SIZE_]; /**< Message to be printed to the output files when the user asks for a time too large or a redshift too low*/
  
  int * tau_out_was_reduced; /**< Logical array to keep track whether a tau_out configuration had tau reduced (or z increased) */


  /**
   * Should SONG compute only a specific set of wavemodes?
   *
   * If _TRUE_, SONG will only compute the wavemodes in k1_out, k2_out and k3_out,
   * and stop execution after the source function has been computed.
   *
   * This is a quick way to output the perturbations and the source function from
   * the perturbations2.c module without having to do a full run.
   */
  short k_out_mode;

  //@}


};



/**
 * Workspace passed between the functions in the perturbations2.c module, containing
 * useful information about the state of the second-order differential system. 
 *
 * This workspace is shared between the functions called by perturb2_solve(). It is
 * updated at each time step and for each (k1,k2,k3) triplet.
 * 
 * The most important data contained in the workspace are:
 *
 * - Geometrical quantities such as scalar products, rotation coefficients and
 *   coupling coefficients.
 *
 * - Arrays for the interpolation in time of the background and thermodynamical
 *   quantities, the first-order perturbations and the quadratic sources.
 *
 * - The density, velocity and shear of the various species.
 *
 * - The value of approximation-dependent perturbations, such as the photon shear
 *   in the tight coupling approximation.
 *
 * - The vector of evolved perturbations, pv, which contains the state of the
 *   differential system at the current time.
 *
 */
struct perturb2_workspace 
{

  // ====================================================================================
  // =                                Evolution variables                               =
  // ====================================================================================
  
  struct perturb2_vector * pv;  /**< Pointer to vector of evolved perturbations containing the
                                state of the differential system at the current time */

  double tau_start_evolution;   /**< Conformal time when we start evolving the current wavemode */
  
  int l_max_g;                /**< Number of multipoles to evolve in the second-order Boltzmann hierarchy for the
                              photon intensity. The highest this number, the smallest the noise from numerical 
                              reflection (see Sec. 5.3.1.2 of http://arxiv.org/abs/1405.2280). */
  int n_hierarchy_g;          /**< Number of equations to evolve in the Boltzmann hierarchy for the photon intensity */

  int l_max_pol_g;            /**< Number of multipoles to evolve in the second-order Boltzmann hierarchy for the
                              photon polarisation. The highest this number, the smallest the noise from numerical 
                              reflection (see Sec. 5.3.1.2 of http://arxiv.org/abs/1405.2280). */
  int n_hierarchy_pol_g;      /**< Number of equations to evolve in the Boltzmann hierarchy for the photon polarisation */

  int n_hierarchy_b;          /**< Number of equations to evolve in the Boltzmann hierarchy for the baryons */
  
  int n_hierarchy_cdm;        /**< Number of equations to evolve in the Boltzmann hierarchy for the cold dark matter */

  int l_max_ur;               /**< Number of multipoles to evolve in the second-order Boltzmann hierarchy for the
                              neutrino intensity. The highest this number, the smallest the noise from numerical 
                              reflection (see Sec. 5.3.1.2 of http://arxiv.org/abs/1405.2280). */
  int n_hierarchy_ur;         /**< Number of equations to evolve in the Boltzmann hierarchy for the neutrino intensity */
  


  
  // ====================================================================================
  // =                                Geometric variables                               =
  // ====================================================================================

  /**
   * Variables related to the geometry of the Fourier wavemodes.
   *
   * At second order, the transfer functions depend on three Fourier wavemodes rather than
   * one, due to mode coupling. Statistical isotropy allows us to reduce the degress of
   * freedom associated to the three vectors to three degrees of freedom, which we
   * choose to be their magnitudes (k1,k2,k3). In particular, we chose the k3 vector
   * to be aligned with the z-axis.
   *
   * Here we store the value of (k1,k2,k3) along with many derived quantities: angles,
   * sines, cosines, scalar products, tensorial products. 
   */
  //@{

  int index_k1;   /**< index in ppt->k of the k1 wavemode currently being evolved */
  int index_k2;   /**< index in ppt->k of the k2 wavemode currently being evolved */
  int index_k3;   /**< index in ppt->k3[index_k1][index_k2] of the k3 wavemode currently being evolved */
  double k1;      /**< magnitude of the k1 vector currently being evolved */
  double k2;      /**< magnitude of the k2 vector currently being evolved */
  double k;       /**< magnitude of the k3 vector currently being evolved */
  double k_sq;    /**< square of the magnitude of the k3 vector currently being evolved */
  
  double cosk1k;   /**< cosine of the angle between k1 and k, assuming k is aligned with the z axis */
  double cosk2k;   /**< cosine of the angle between k2 and k, assuming k is aligned with the z axis */
  double cosk1k2;  /**< cosine of the angle between k1 and k2, assuming k is aligned with the z axis */
  double sink1k;   /**< sine of the angle between k1 and k, assuming k is aligned with the z axis */
  double sink2k;   /**< sine of the angle between k2 and k, assuming k is aligned with the z axis */
  double theta_1;  /**< angle between k and k1 in radians */
  double theta_2;  /**< angle between k and k2 in radians */

  double k1_m[3];  /**< spherical components of the k1 wavemode, defined as k1_m[m+1] = xi[m]_i k1^i, where
                   xi[m]^i is the spherical basis described in sec. A.3.1 of http://arxiv.org/abs/1405.2280. */
  double k2_m[3];  /**< spherical components of the k2 wavemode, defined as k2_m[m+1] = xi[m]_i k2^i, where
                   xi[m]^i is the spherical basis described in sec. A.3.1 of http://arxiv.org/abs/1405.2280. */

  double k1_dot_k2;  /**< scalar products between the k1 and k2 vectors */

  double k1_ten_k2[5];  /**< tensorial product between k1 and k2, defined as k1_ten_k2[m+2] = X[m]^ij k1_i k2_j, where
                        X[m]^ij are the projection matrix described in Sec. A.3.2 of http://arxiv.org/abs/1405.2280. */
  double k1_ten_k1[5];  /**< tensorial product between k1 and k1, defined as k1_ten_k1[m+2] = X[m]^ij k1_i k1_j, where
                        X[m]^ij are the projection matrix described in Sec. A.3.2 of http://arxiv.org/abs/1405.2280. */
  double k2_ten_k2[5];  /**< tensorial product between k2 and k2, defined as k2_ten_k2[m+2] = X[m]^ij k2_i k2_j, where
                        X[m]^ij are the projection matrix described in Sec. A.3.2 of http://arxiv.org/abs/1405.2280. */

  /* Rotation coefficients needed to implement the scalar linear perturbations into
  our second-order system. These arrays depend on the (l,m) multipole considered
  and on the cosine of the angle between k and either k1 or k2. They are indexed as
  ppw2->rotation_1[lm_quad(l,m)]. */
  double *rotation_1;           /**< Rotation coefficients for k1, indexed as ppw2->rotation_1[lm_quad(l,m)]; for more 
                                details, refer to perturb2_geometrical_corner(). */
  double *rotation_2;           /**< Rotation coefficients for k2, indexed as ppw2->rotation_2[lm_quad(l,m)]; for more 
                                details, refer to perturb2_geometrical_corner(). */
  double *rotation_1_minus;     /**< Rotation coefficients for k1 (negative m), indexed as ppw2->rotation_1_minus[lm_quad(l,m)];
                                for more details, refer to perturb2_geometrical_corner(). */
  double *rotation_2_minus;     /**< Rotation coefficients for k2 (negative m), indexed as ppw2->rotation_2_minus[lm_quad(l,m)];
                                for more details, refer to perturb2_geometrical_corner(). */

  //@}



  // ====================================================================================
  // =                               Coupling coefficients                              =
  // ====================================================================================

  /**
   * Summed coupling coefficients.
   *
   * In the quadratic part of the second-order Boltzmann equation, the coupling
   * coefficients appear in a sum over m1 with two first-order perturbations
   * (see Sec. 4.6.2 of http://arxiv.org/abs/1405.2280).
   *
   * The m-dependent part of a first-order perturbation, assuming no initial
   * vector or tensor modes, can be factored in a geometrical rotation coefficient
   * (see documentation for rot_1 and rot_2).
   *   
   * This means that we can factor out the perturbations from the sum over m1.
   * Then, the sum will contain geometrical factors, which can be computed and
   * summed beforehand. We compute these factors and the whole sum in
   * perturb2_geometrical_corner(), and access the result with these macros.
   * 
   * The summed coupling coefficients are usually accessed via the shorthands
   * defined in perturbations2_macros.h.
   */
  //@{

  /* Intensity couplings */
  double * c_minus_product_12;
  double * c_minus_product_21;
  double * c_plus_product_12;
  double * c_plus_product_21;
  double * c_minus_product_11;
  double * c_minus_product_22;
  double * c_plus_product_11;
  double * c_plus_product_22;

  double * r_minus_product_12;
  double * r_minus_product_21;
  double * r_plus_product_12;
  double * r_plus_product_21;

  /* E-mode polarization couplings */
  double * d_minus_product_12;
  double * d_minus_product_21;
  double * d_plus_product_12;
  double * d_plus_product_21;
  double * d_minus_product_11;
  double * d_minus_product_22;
  double * d_plus_product_11;
  double * d_plus_product_22;

  double * k_minus_product_12;
  double * k_minus_product_21;
  double * k_plus_product_12;
  double * k_plus_product_21;
  double * k_minus_product_11;
  double * k_minus_product_22;
  double * k_plus_product_11;
  double * k_plus_product_22;
  
  /* B-mode polarization couplings */
  double * d_zero_product_12;
  double * d_zero_product_21;
  double * d_zero_product_11;
  double * d_zero_product_22;

  double * k_zero_product_12;
  double * k_zero_product_21;
  double * k_zero_product_11;
  double * k_zero_product_22;
  
  //@}



  // =============================================================================================
  // =                                      Metric indices                                       =
  // =============================================================================================

  /**
   * Indices for accessing the metric variables.
   *
   * We compute the metric variables separately from the other perturbations, in
   * perturb2_einstein(), and store their values in the vector ppw2->pvecmetric(),
   * which is accessed via the index_mt2_XXX indices defined here.
   * 
   * Refer to sec. 3.3 of http://arxiv.org/abs/1405.2280 for details on the metric
   * variables used by SONG, and Sec. 5.3.1.1 for the equations governing them.
   */

  /* Newtonian gauge */
  int index_mt2_psi;                      /**< Time potential psi of Newtonian gauge */
  int index_mt2_phi_prime;                /**< (d phi/d tau), with phi curvature potential of Newtonian gauge  */
  int index_mt2_phi_prime_prime;          /**< (d^2 phi/d tau^2) with phi curvature potential of Newtonian gauge.
                                          Used only if ppt2->phi_eq==huang. */
  int index_mt2_phi_prime_poisson;        /**< (d phi/d tau) with phi curvature potential of Newtonian gauge, using
                                          the time-time Einstein equation (eq. 5.2 of http://arxiv.org/abs/1405.2280). */
  int index_mt2_phi_prime_longitudinal;   /**< (d phi/d tau) with phi curvature potential of Newtonian gauge, using
                                          the space-time Einstein equation (eq. 3.98 of http://arxiv.org/abs/1405.2280). */
  int index_mt2_omega_m1_prime;           /**< vector mode of the metric in Newtonian gauge */
  int index_mt2_gamma_m2_prime_prime;     /**< tensor mode of the metric in Newtonian gauge */           

  /* Synchronous gauge */
  int index_mt2_h_prime;         /**< (d h/d tau) in synchronous gauge */
  int index_mt2_h_prime_prime;   /**< (d^2 h/d tau^2) in synchronous gauge */
  int index_mt2_eta_prime;       /**< (d eta/d tau) in synchronous gauge */
  int index_mt2_alpha_prime;     /**< (d alpha/d tau) in synchronous gauge, where alpha=(h'+6 eta')/(2 k^2) */

  int mt2_size;                  /**< Size of the metric vector ppw2->pvecmetric */

 

  // ====================================================================================
  // =                                 Quadratic sources                                =
  // ====================================================================================

  double ** quadsources_table; /**< Table that will contain all the quadratic sources needed by
                               SONG to solve the differential system for the current wavemode.
                               Indexed as ppw2->quadsources_table[index_qs2_XXX][index_tau] */

  double ** d_quadsources_table;  /**< First-order time derivative of quadsources_table,
                                  needed for some approximations */

  double ** dd_quadsources_table; /**< Second-order time derivative of quadsources_table,
                                  needed for spline interpolation */
  
  double ** ddd_quadsources_table; /**< Third-order time derivative of quadsources_table,
                                  needed for spline interpolation of the first-derivative */
  
  double ** dddd_quadsources_table; /**< Fourth-order time derivative of quadsources_table,
                                    needed for spline interpolation of the second-derivative */
  
  double ** quadcollision_table; /**< Table that will contain the quadratic part of the collision
                                 term. Needed to compute the tight coupling approximation.
                                 Indexed as ppw2->quadsources_table[index_qs2_XXX][index_tau] */

  double ** d_quadcollision_table; /**< First-order time derivative of quadcollision_table,
                                   needed for some approximations */

  double ** dd_quadcollision_table; /**< Second-order time derivative of quadcollision_table,
                                    needed for spline interpolation */

  double ** ddd_quadcollision_table; /**< Third-order time derivative of quadcollision_table,
                                     needed for spline interpolation of the first-derivative */
  
  double ** dddd_quadcollision_table; /**< Fourth-order time derivative of quadcollision_table,
                                      needed for spline interpolation of the second-derivative */
  
  int qs2_size; /**< Size of ppw2->quadsources_table and ppw2->pvec_quadsources */

  /* Quadratic sources for the metric, Newtonian gauge */
  int index_qs2_psi;
  int index_qs2_psi_prime;
  int index_qs2_phi_prime;
  int index_qs2_phi_prime_prime;               
  int index_qs2_phi_prime_poisson;      
  int index_qs2_phi_prime_longitudinal;
  int index_qs2_omega_m1_prime;
  int index_qs2_gamma_m2_prime_prime;

  /* Quadratic sources for the metric, synchronous gauge */
  int index_qs2_h_prime;
  int index_qs2_h_prime_prime;
  int index_qs2_eta_prime;
  int index_qs2_alpha_prime;  

  /* Quadratic sources for each matter species */
  int index_qs2_monopole_g;
  int index_qs2_monopole_E;
  int index_qs2_monopole_B;
  int index_qs2_monopole_b;
  int index_qs2_monopole_cdm;
  int index_qs2_monopole_ur;

  /* Quadratic sources for the magnetic field */
  int index_qs2_M;

  /* Other useful quadratic sources. */
  int index_qs2_vv_g;    /**< (velocity potential)^2 of photons */
  int index_qs2_vv_ur;   /**< (velocity potential)^2 of neutrinos */
  int index_qs2_dd_b;    /**< (density contrast)^2 of baryons */
  int index_qs2_vv_b;    /**< (velocity potential)^2 of baryons */
  int index_qs2_vv_cdm;  /**< (velocity potential)^2 of CDM */



  // ====================================================================================
  // =                                Time interpolation                                =
  // ====================================================================================

  double * pvecback;              /**< Interpolated values of the background quantitites at the current time tau; indexed by the pba->index_ba_XXX indices. */
  double * pvecthermo;            /**< Interpolated values of the thermodynamics quantitites at the current time tau; indexed by the pth->index_th_XXX indices. */
  double * pvecmetric;            /**< Metric quantitites at the current time tau; indexed by the ppw2->index_mt2_XXX indices. */
  double * pvec_quadsources;      /**< Interpolated/computed values of the quadratic sources at the current time tau; indexed
                                  by ppw2->index_qs2_XXX indices. Filled with perturb2_quadratic_sources() (for computation) or
                                  with perturb2_quadratic_sources_at_tau() (for interpolation). */
  double * pvec_quadcollision;    /**< Interpolated/computed values of the quadratic collisional sources at the current time tau; indexed by ppw2->index_qs2_XXX indices. */

  double * pvec_d_quadsources;    /**< Interpolated values of the conformal-time derivatives of the quadratic sources at the current time tau; indexed by ppw2->index_qs2_XXX indices. */
  double * pvec_d_quadcollision;  /**< Interpolated values of the conformal-time derivatives of the quadratic collisional sources at the current time tau; indexed by ppw2->index_qs2_XXX indices. */
  double * pvec_dd_quadsources;   /**< Interpolated values of the conformal-time second derivatives of the quadratic sources at the current time tau; indexed by ppw2->index_qs2_XXX indices. */
  double * pvec_dd_quadcollision; /**< Interpolated values of the conformal-time second derivatives of the quadratic collisional sources at the current time tau; indexed by ppw2->index_qs2_XXX indices. */

  double * pvec_sources1;       /**< Interpolated values of the first-order perturbations in k1 and tau; indexed by
                                ppt->index_qs_XXX indices and filled by perturb_song_sources_at_tau() */
  double * pvec_sources2;       /**< Interpolated values of the first-order perturbations in k2 and tau; indexed by
                                ppt->index_qs_XXX indices and filled by perturb_song_sources_at_tau() */

  int last_index_back;         /**< Keep track of the row where we last accessed the background table, useful to
                               speed up the interpolation performed by background_at_tau() */
  int last_index_thermo;       /**< Keep track of the row where we last accessed the thermodynamics table, useful to
                               speed up the interpolation performed by thermodynamics_at_z() */
  int last_index_sources;      /**< Keep track of the row where we last accessed the ppt->quadsources table, useful to
                               speed up the interpolation performed by perturb_song_sources_at_tau() */



  // ====================================================================================
  // =                                Fluid variables                                   =
  // ====================================================================================

  /** Photon fluid variables; filled & explained in perturb2_fluid_variables() */
  //@{
  double delta_g_1, delta_g_2, u_g_1[2], u_g_2[2];
  double v_dot_v_g, v_ten_v_g[3];
  double delta_g, delta_g_adiab, u_g[2], pressure_g, shear_g[3];
  //@}

  /** Baryon fluid variables; filled & explained in perturb2_fluid_variables() */
  //@{
  double delta_b_1, delta_b_2, delta_b_1_prime, delta_b_2_prime;
  double u_b_1[2], u_b_2[2], u_b_1_prime[2], u_b_2_prime[2];
  double v_dot_v_b, v_ten_v_b[3], v_ten_v_b_prime[3];
  double delta_b, u_b[2], pressure_b, shear_b[3];
  //@}

  /** Cold dark matter fluid variables; filled & explained in perturb2_fluid_variables() */
  //@{
  double delta_cdm_1, delta_cdm_2, u_cdm_1[2], u_cdm_2[2];
  double v_dot_v_cdm, v_ten_v_cdm[3];
  double delta_cdm, u_cdm[2], pressure_cdm, shear_cdm[3];
  //@}
  
  /** Ultra relativistic species fluid variables; filled & explained in perturb2_fluid_variables() */
  //@{
  double delta_ur_1, delta_ur_2, u_ur_1[2], u_ur_2[2];
  double v_dot_v_ur, v_ten_v_ur[3];
  double delta_ur, u_ur[2], pressure_ur, shear_ur[3];
  //@}
  


  // ====================================================================================
  // =                                  Approximations                                  =
  // ====================================================================================

  int * approx;        /**< Logical array of active approximations at a given time: approx[ppw2->index_ap2_XXX] */
  int ap2_size;        /**< Number of possible approximation schemes and size of ppw2->approx */
  
  int index_ap2_tca;   /**< Index for the tight-coupling approximation in the array ppw2->approx */
  int index_ap2_rsa;   /**< Index for the radiation streaming approximation in the array ppw2->approx */
  int index_ap2_ufa;   /**< Index for the ultra-relativistic fluid approximation in the array ppw2->approx */
  int index_ap2_nra;   /**< Index for the no-radiation approximation in the array ppw2->approx */

  int n_active_approximations; /**< Number of approximations active for the current (k1,k2,k3) and time tau */

  double I_1m_tca1[2];     /**< Photon intensity dipole in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double u_g_tca1[2];      /**< Photon velocity u_g[m]=i*v_g[m] in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double I_2m_tca0[3];     /**< Photon quadrupole in the tight coupling approximations, neglecting O(tau_c) terms */
  double I_2m_tca1[3];     /**< Photon quadrupole in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double shear_g_tca1[3];  /**< Photon shear in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double Pi_tca0[3];       /**< Pi factor in the tight coupling approximations, neglecting O(tau_c) terms */
  double Pi_tca1[3];       /**< Pi factor in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double E_2m_tca1[3];     /**< E-mode quadrupole in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double B_2m_tca1[3];     /**< B-mode quadrupole in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double U_slip_tca1[2];   /**< Velocity slip U[m]=u_b[m]-u_g[m] in the tight coupling approximations, neglecting O(tau_c)^2 terms */
  double C_1m_tca1[2];     /**< Purely second-order part of the photon dipole collision term,
                           C_1m = kappa_dot * (4/3*b_11m-I_1m), in the tight coupling approximations,
                           neglecting O(tau_c)^2 terms. */

  double I_00_rsa;       /**< Value of the photon intensity monopole in the radiation streaming approximations */
  double I_1m_rsa[2];    /**< Value of the photon intensity dipole in the radiation streaming approximations */
  double N_00_rsa;       /**< Value of the neutrino monopole in the radiation streaming approximations */
  double N_1m_rsa[2];    /**< Value of the neutrino dipole in the radiation streaming approximations */
  double delta_g_rsa;    /**< Value of the photon density contrast in the radiation streaming approximations */
  double u_g_rsa[2];     /**< Value of the photon velociy in the radiation streaming approximations */
  double delta_ur_rsa;   /**< Value of the neutrino density contrast in the radiation streaming approximations */
  double u_ur_rsa[2];    /**< Value of the neutrino velocity in the radiation streaming approximations */


  double I_00;         /**< The photon intensity monopole fed to the evolver (depends on RSA) */
  double N_00;         /**< The neutrino monopole fed to the evolver (depends on RSA) */
  double N_1m[2];      /**< The neutrino dipole fed to the evolver (depends on RSA) */
  double I_1m[2];      /**< The photon intensity dipole fed to the evolver (depends on TCA and RSA) */
  double I_2m[3];      /**< The photon intensity quadrupole fed to the evolver (depends on TCA) */
  double E_2m[3];      /**< The E-polarisation quadrupole fed to the evolver (depends on TCA) */
  double B_2m[3];      /**< The B-polarisation quadrupole fed to the evolver (depends on TCA) */
  double C_1m[2];      /**< Purely second-order part of the photon dipole collision term,
                       C_1m = kappa_dot * (4/3*b_11m-I_1m); the collision term for baryons is
                       the same times -r=-rho_g/rho_b. */
  double b_200;        /**< The baryon "pressure" fed to the evolver (depends on perfect fluid approximation) */
  double b_22m[3];     /**< The baryon quadrupole fed to the evolver (depends on perfect fluid approximation) */
  double cdm_200;      /**< The CDM "pressure" fed to the evolver (depends on perfect fluid approximation) */
  double cdm_22m[3];   /**< The CDM quadrupole fed to the evolver (depends on perfect fluid approximation) */



  // ====================================================================================
  // =                                    Debug                                         =
  // ====================================================================================
  
  int derivs_calls;   /**< Counter to keep track of how many times the function perturb2_derivs()
                      has been called for the considered set of (k1,k2,k3). */

  int sources_calls;  /**< Counter to keep track of how many times the function perturb2_sources()
                      has been called for the considered set of (k1,k2,k3). */

  long int n_steps;   /**< Number of steps taken by the differential system so far for the active
                      (k1,k2,k3) mode. Computed only for those triplets belonging to k1_out,
                      k2_out and k3_out. */

  char info [_MAX_INFO_SIZE_];   /**< String with information on the wavemode that is being integrated */
  
  
  /**
   * Function used to output intermediate values from the differential system.
   *
   * This function will be given as an argument to generic_evolver and is  used only for
   * debug and output purposes. By default it is set to NULL, which means that it is never
   * called.  It is called only if the evolved (k1,k2,k3) triplet corresponds to one of
   * the triplets requested via the parameters k1_out, k2_out and k3_out.
   */
  int (*print_function)(double x, int index_x, double y[], double dy[],
    void *parameters_and_workspace, ErrorMsg error_message);

  int index_k_out;    /**< Variable passed to print_function() to decide whether to append data to
                      the k_out files. If negative, the perturbations will not be outputted to file
                      for the current (k1,k2) pair. If positive, index_k_out is the index inside
                      k1_out and k2_out of the current (k1,k2) pair. */
  
  int index_k_out_for_tau_out;  /**< Variable passed to print_function() to decide whether to append data to
                                the tau_out files. If negative, the perturbations will not be outputted to the
                                tau_out files for the current (k1,k2) pair. If positive, it is the index inside
                                k1_out and k2_out of the current (k1,k2) pair. The difference with respect to 
                                index_k_out is that index_k_out_for_tau_out is constant for triplets with
                                identical values of (k1_out,k2_out), thus allowing SONG to output data to a
                                single tau_out file for these triplets. */
  
  char file_header[3*_MAX_INFO_SIZE_]; /**< Header to be added on top of the output files */  

  

};





/**
 * Structure with the state of the differential system.
 *
 * This structure is what is filled and read by the differential evolver.
 * It is also used to set the initial conditions and to set which equations
 * to evolve.
 * 
 * The state of the diffential system is given by:
 *
 * - The list of perturbations currenty being evolved, indexed by the 
 *   indices index_pt2_XXX.
 *
 * - The values of the evolved perturbations, stored in y.
 *
 * - The derivatives of the perturbations with respect to conformal time,
 *   stored in dy.
 * 
 * Note that the number of perturbations to be evolved, and hence the size of
 * y and dy, vary depending on the active approximation schemes.
 *
 */

struct perturb2_vector
{

  // ====================================================================================
  // =                                 Number of equations                              =
  // ====================================================================================

  /**
   * Paremeters used to determine the number of equation to evolve.
   *
   * The number of evolved equations varies depending on the active approximation. We
   * split the integration of a given wavemode in time intervals. Each time interval
   * corresponds to a different approximation (or lack thereof) and will be integrated
   * by the evolver separately.
   *
   * The following variables reflect the varying number of evolved equations; they
   * will be updated at the beginning of every time interval by the
   * perturb2_vector_init() function.
   *
   * The massive species (baryons and cold dark matter) need to be described by three
   * indices (n,l,m) rather than just (l,m). The n parameter arises from the velocity
   * expansion of the distribution function. For massless particles, n is irrelevant
   * because they have constant velocity. We call the (n,l,m) moments as beta-moments;
   * they are a natural way to treat massive particles in the Boltzmann formalism.
   * For more details on the beta moments, please refer to sec. 5.3.1.3 of
   * http://arxiv.org/abs/1405.2280; see also Lewis and Challinor 2002.
   */
  //@{

  int l_max_g;               /**< Largest multipole evolved in the photon intensity hierarchy */
  int l_max_pol_g;           /**< Largest multipole evolved in the photon polarisation hierarchies (both E and B-modes) */
  int l_max_ur;              /**< Largest multipole evolved in the neutrino hierarchy */
  int l_max_b;               /**< Largest multipole evolved in the baryon hierarchy */
  int l_max_cdm;             /**< Largest multipole evolved in the cold dark matter hierarchy */

  int n_hierarchy_g;         /**< Number of equations evolved for the photon temperature Boltzmann hierarchy */
  int n_hierarchy_pol_g;     /**< Number of equations evolved for each polarisation hierarchy */
  int n_hierarchy_b;         /**< Number of equations evolved for the baryon Boltzmann hierarchy */
  int n_hierarchy_cdm;       /**< Number of equations evolved for the cold dark matter Boltzmann hierarchy */
  int n_hierarchy_ur;        /**< Number of equations evolved for the neutrino Boltzmann hierarchy */

  int n_max_b;     /**< Largest velocity moment evolved for the baryon hierarchy;
                   if has_perfect_baryons==_TRUE_, it is equal to 1, otherwise 2. */
  int n_max_cdm;   /**< Largest velocity moment evolved for the cold dark matter hierarchy;
                   if has_perfect_baryons==_TRUE_, it is equal to 1, otherwise 2. */

  short use_closure_g;     /**< Should we use the closure relations for the photon intensity? Usually turned off if an approximation is turned on. */
  short use_closure_pol_g; /**< Should we use the closure relations for the polarisation? Usually turned off if an approximation is turned on. */
  short use_closure_ur;    /**< Should we use the closure relations for the neutrinos? Usually turned off if an approximation is turned on. */

  //@}



  // ====================================================================================
  // =                               Evolved perturbations                              =
  // ====================================================================================

  /* Radiation hierarchies */
  int index_pt2_monopole_g;       /**< Index inside y and dy for the photon temperature hierarchy */
  int index_pt2_monopole_E;       /**< Index inside y and dy for the photon E-mode polarization hierarchy */  
  int index_pt2_monopole_B;       /**< Index inside y and dy for the photon B-mode polarization hierarchy */    
  int index_pt2_monopole_ur;      /**< Index inside y and dy for the neutrino hierarchy */
  
  /* Matter hierarchies */
  int index_pt2_monopole_b;       /**< Index inside y and dy for the baryon hierarchy */
  int index_pt2_monopole_cdm;     /**< Index inside y and dy for the cold dark matter hierarchy */

  /* Metric variables */
  int index_pt2_phi;              /**< Index inside y and dy for the Newtonian gauge curvature potential phi */
  int index_pt2_phi_prime;        /**< Index inside y and dy for the the time derivative of the Newtonian gauge curvature potential phi */
  int index_pt2_omega_m1;         /**< Index inside y and dy for the Newtonian gauge vector potential omega_[m=1] */
  int index_pt2_gamma_m2;         /**< Index inside y and dy for the Newtonian gauge tensor potential gamma_[m=2] */
  int index_pt2_gamma_m2_prime;   /**< Index inside y and dy for the time derivative of gamma_[m=2] */
  int index_pt2_eta;              /**< Index inside y and dy for the synchronous gauge metric perturbation eta */

  /* Magnetic fields */
  int index_pt2_M;                /**< Index inside y and dy for the the magnetic field */

  int pt2_size;                   /**< Number of evolved perturbations in the considered time interval, and size of the
                                  perturbation vectors y and dy */

  /** Array of perturbations to be integrated. It is filled by the evolver at the end of
  each time step. */
  double * y;

  /** Vector containing the time-derivative of the evolved perturbations contained in y.
  It is filled by each call of perturb2_derivs(). */
  double * dy;               

  /** Boolean array specifying which perturbations enter in the calculation of the
  source functions. Only the marked perturbations will be interpolated at the
  times requested in ppt2->tau_sampling. */
  int * used_in_sources;

  /** Array of strings that contain the labels of the various evolved perturbations.
  For debug purposes only. */
  char (*pt2_labels)[_MAX_LENGTH_LABEL_];


};




/**
 * Structure with all the data needed by the evolver.
 *
 * The differential evolver is a general purpose function that does not know anything
 * about SONG or cosmology. We pass such information via a void pointer to the following
 * structure, which contains all the information needed to run perturb2_derivs() and
 * perturb2_sources().
 */

struct perturb2_parameters_and_workspace {

  struct precision * ppr;               /**< Pointer to the precision structure */
  struct precision2 * ppr2;             /**< Pointer to the precision2 structure */
  struct background * pba;              /**< Pointer to the background structure */
  struct thermo * pth;                  /**< Pointer to the thermodynamics structure */
  struct perturbs * ppt;                /**< Pointer to the perturbation structure */
  struct perturbs2 * ppt2;              /**< Pointer to the second-order perturbation structure */
  struct perturb2_workspace * ppw2;     /**< Pointer to the computation workspace */
    
};





/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
  extern "C" {
#endif

    int perturb2_init(
         struct precision * ppr,
         struct precision2 * ppr2,
         struct background * pba,
         struct thermo * pth,
         struct perturbs * ppt,
         struct perturbs2 * ppt2
         );

    int perturb2_free(
         struct precision2 * ppr2,
         struct perturbs2 * ppt2
         );
         

    int perturb2_indices_of_perturbs(
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background * pba,
            struct thermo * pth,
            struct perturbs * ppt,
            struct perturbs2 * ppt2
            );

    int perturb2_timesampling_for_sources(
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2
           );

    int perturb2_start_time_evolution (
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background * pba,
            struct thermo * pth,
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double k,
            double * tau_ini
            );

    int perturb2_end_of_recombination (
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2
          );
    
    int perturb2_get_lm_lists(
         struct precision * ppr,
         struct precision2 * ppr2,
         struct background * pba,
         struct thermo * pth,
         struct perturbs * ppt,
         struct perturbs2 * ppt2
         );


    int perturb2_get_k_lists(
         struct precision * ppr,
         struct precision2 * ppr2,
         struct background * pba,
         struct thermo * pth,
         struct perturbs * ppt,
         struct perturbs2 * ppt2
         );

    int perturb2_workspace_init(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2,                
          struct perturb2_workspace * ppw2
          );

    int perturb2_workspace_init_quadratic_sources(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2,
          struct perturb2_workspace * ppw2
          );


    int perturb2_workspace_free(
          struct perturbs2 * ppt2,
          struct background * pba,
          struct perturb2_workspace * ppw2
          );

    int perturb2_vector_init(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2,
          double tau,
          struct perturb2_workspace * ppw2,
          int * old_approx
          );

    int perturb2_vector_free(
          struct perturb2_vector * pv
          );

    int perturb2_geometrical_corner(
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background * pba,
            struct thermo * pth,
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            int index_k1,
            int index_k2,
            int index_k3,                    
            struct perturb2_workspace * ppw2
            );


    int perturb2_solve(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct perturbs2 * ppt2,
          int index_k1,
          int index_k2,
          int index_k3,                    
          struct perturb2_workspace * ppw2
          );

    int perturb2_initial_conditions(
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           struct perturb2_workspace * ppw2
           );

    int perturb2_timescale(
          double tau,
          void * parameters_and_workspace,
          double * timescale,
          ErrorMsg error_message
          );


  int perturb2_find_approximation_switches(
					struct precision * ppr,
					struct precision2 * ppr2,
					struct background * pba,
					struct thermo * pth,
					struct perturbs * ppt,
					struct perturbs2 * ppt2,
					struct perturb2_workspace * ppw2,
					double tau_ini,
					double tau_end,
					double precision,
					int interval_number,
					int * interval_number_of,
					double * interval_limit, /* interval_limit[index_interval] (already allocated) */
					int ** interval_approx   /* interval_approx[index_interval][index_ap] (already allocated) */
          );


    int perturb2_find_approximation_number(
					  struct precision * ppr,
					  struct precision2 * ppr2,
					  struct background * pba,
					  struct thermo * pth,
					  struct perturbs * ppt,
					  struct perturbs2 * ppt2,
					  struct perturb2_workspace * ppw2,
					  double tau_ini,
					  double tau_end,
					  int * interval_number,
					  int * interval_number_of
					  );


    int perturb2_approximations(
         struct precision * ppr,
         struct precision2 * ppr2,
         struct background * pba,
         struct thermo * pth,
         struct perturbs * ppt,
         struct perturbs2 * ppt2,
         double tau,
         struct perturb2_workspace * ppw2
         );

    int perturb2_workspace_at_tau (
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           double * y,
           struct perturb2_workspace * ppw2
           );
       
    int perturb2_tca_variables (
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           double * y,
           struct perturb2_workspace * ppw2
           );

    int perturb2_rsa_variables (
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           double * y,
           struct perturb2_workspace * ppw2
           );

    int perturb2_fluid_variables (
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           double * y,
           struct perturb2_workspace * ppw2
           );

    int perturb2_einstein(
       struct precision * ppr,
       struct precision2 * ppr2,
       struct background * pba,
       struct thermo * pth,
       struct perturbs * ppt,
       struct perturbs2 * ppt2,
       double tau,
       double * y,
       struct perturb2_workspace * ppw2
       );

    int perturb2_compute_psi_prime(
           struct precision * ppr,
           struct precision2 * ppr2,
           struct background * pba,
           struct thermo * pth,
           struct perturbs * ppt,
           struct perturbs2 * ppt2,
           double tau,
           double * y,
           double * dy,
           double * psi_prime,
           struct perturb2_workspace * ppw2);

    int perturb2_sources(
          double tau,
          double * y,
          double * dy,
          int index_tau,
          void * parameters_and_workspace,
          ErrorMsg error_message
          );

    int perturb2_save_perturbations(
          double tau,
          int index_tau,
          double * y,
          double * dy,
          void * parameters_and_workspace,
          ErrorMsg error_message
          );

    int perturb2_print_variables(double tau,
        double * y,
        double * dy,
        int index_tau,        
        void * parameters_and_workspace,
        ErrorMsg error_message
        );
  
    int what_if_ndf15_fails(int (*derivs)(double x, 
              double * y, 
              double * dy, 
              void * parameters_and_workspace,
              ErrorMsg error_message),
            double x_ini,
            double x_end,
            double * y, 
            int * used_in_output,
            int y_size,
            void * parameters_and_workspace_for_derivs,
            double tolerance, 
            double minimum_variation,
            int (*evaluate_timescale)(double x, 
              void * parameters_and_workspace,
              double * timescale,
              ErrorMsg error_message),
            double timestep_over_timescale,
            double * x_sampling,
            int x_size,
            int (*output)(double x,
              double y[],
              double dy[],
              int index_x,
              void * parameters_and_workspace,
              ErrorMsg error_message),
            int (*print_variables)(double x,
                 double y[], 
                 double dy[],
                 void * parameters_and_workspace,
                 ErrorMsg error_message),
            ErrorMsg error_message);

    int perturb2_derivs(
           double tau,
           double * y,
           double * dy,
           void * parameters_and_workspace,
           ErrorMsg error_message
           );

    int perturb2_quadratic_sources(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,            
          struct perturbs * ppt,
          struct perturbs2 * ppt2,
          int index_tau,
          double tau,
          int what_to_compute,
          double * pvec_quadsources,
          double * pvec_quadcollision,
          struct perturb2_workspace * ppw2
          );

    int perturb2_quadratic_sources_for_k1k2k(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct background * pba,
          struct thermo * pth,            
          struct perturbs * ppt,
          struct perturbs2 * ppt2,
          struct perturb2_workspace * ppw2
          );

    int perturb2_quadratic_sources_at_tau(
            struct precision * ppr,
            struct precision2 * ppr2,
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double tau,
            int what_to_interpolate,
            struct perturb2_workspace * ppw2
            );

    int perturb2_quadratic_sources_at_tau_linear(
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double tau,
            int what_to_interpolate,
            struct perturb2_workspace * ppw2
            );

    int perturb2_quadratic_sources_at_tau_spline(
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double tau,
            int what_to_interpolate,
            struct perturb2_workspace * ppw2
            );

    int perturb2_wavemode_info(
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background * pba,
            struct thermo * pth,
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            struct perturb2_workspace * ppw2
            );

    int perturb2_output(
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background * pba,
            struct perturbs * ppt,
            struct perturbs2 * ppt2
            );

    int perturb2_store_sources_to_disk(
            struct perturbs2 * ppt2,
            int index_k1
            );

    int perturb2_store_sources_k3_tau(
            struct perturbs2 * ppt2,
            int index_tp2,
            int index_k1,
            int index_k2,
            char * filepath,
            FILE * output_stream
            );

    int perturb2_load_sources_from_disk(
            struct perturbs2 * ppt2,
            int index_k1
            );

    int perturb2_load_sources_k3_tau(
            struct perturbs2 * ppt2,
            int index_tp2,
            int index_k1,
            int index_k2,
            char * filepath,
            FILE * input_stream
            );

    int perturb2_allocate_k1_level(
         struct perturbs2 * ppt2,
         int index_k1
         );

    int perturb2_free_k1_level(
         struct perturbs2 * ppt2,
         int index_k1
         );


#ifdef __cplusplus
  }
#endif

/**************************************************************/
  
  

#endif