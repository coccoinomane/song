#ifndef __PERTURBATIONS2__
#define __PERTURBATIONS2__

#include "perturbations2_macros.h"
#include "perturbations.h"
#include "common2.h"


// ======================================================================================
// =                                  Enum structures                                   =
// ======================================================================================

/**
 * Implementations of the tight coupling approximation.
 */
enum tca2_method {
  tca2_none,           /**< No TCA approximation */
  tca2_zero            /**< Zero-order TCA approximation (NOT IMPLEMENTED YET) */
};

/**
 * Implementations of the radiation streaming approximation. RSA provides a 
 * way to track the oscillations of the photon monopole and dipole at
 * subhorizon scales, all the way to today.
 */
enum rsa2_method {
  rsa2_none,          /**< No RSA approximation */
};

/**
 * Implementations of the ultra relativistic fluid approximation. RSA allows
 * to track the oscillations of the neutrino monopole and dipole at
 * subhorizon scales.
 */
enum ufa2_method {
  ufa2_none            /**< No UFA approximation */
};

/**
 * No radiation approximation. This is a blunt approximation where, well after
 * matter radiation equality, we either stop evolving massless species or we
 * approximate them with a perfect fluid.
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
 * Which Einstein equation to use for the derivative of the Newtonian potential phi?
 */
enum phi_prime_equation {
  poisson,      /**< Use Poisson equation (eq 5.2 of http://arxiv.org/abs/1405.2280) */
  longitudinal  /**< Use the longitudinal equation (eq 3.98 of http://arxiv.org/abs/1405.2280) */
};

/**
 * Possible sampling methods for the ppt2->k array
 */
enum sources2_k_sampling {
  lin_k_sampling,                  /**< Linear k sampling */
  log_k_sampling,                  /**< Logarithmic k sampling */
  class_sources_k_sampling,        /**< k sampling adopted in perturb_get_k_list */
  smart_sources_k_sampling         /**< smart k-sampling, logarithmic + linear */
};


/**
 * Possible sampling methods for the ppt2->k3[index_k1][index_k2] array
 */
enum sources2_k3_sampling {
  lin_k3_sampling,                  /**< Linear k sampling */
  log_k3_sampling,                  /**< Logarithmic k sampling */
  smart_k3_sampling,                /**< k sampling adopted in perturb_get_k3_list */
  theta12_k3_sampling,              /**< Sampling linear in the angle between k1 and k2 */
  theta13_k3_sampling               /**< Sampling linear in the angle between k1 and k3 */
};


// ---------------------------------------------------------------------------------
// -                             Precision parameters                              -
// ---------------------------------------------------------------------------------

/* The differential system at second-order has to be solved for a set of three wavemodes
(k1,k2,k3) whereby k3 has to be in the range |k1-k2|<=k3<=k1+k1. When k3 is too close
to the boundaries, numerical instabilities might arise such as nan's or larger than one
sines and cosines. In order to avoid that, we define here a safety distance between
k3 and the bounds. This safety distance is going to correspond to the largest scale
probed by SONG. Using k_min_tau0=1e-3, that corresponds to k_min=1e-8,
it seems that setting _MIN_K3_DISTANCE_=1e-10 is ok. */
#define _MIN_K3_DISTANCE_ 1e-10
#define _MIN_K3_RATIO_ 100



struct perturbs2
{

  // ====================================================================================
  // =                                    Output flags                                  =
  // ====================================================================================

  /* Do we need to compute 2nd order perturbations at all? */
  short has_perturbations2;

  /* These flags are initialised in the input2.c module and are used by all other modules to determine
  which probes (T,E,B...) to compute */
  short has_cmb_temperature;          /* Do we need to compute spectra or bispectra for the CMB temperature? */
  short has_cmb_polarization_e;       /* Do we need to compute spectra or bispectra for the CMB E-modes? */
  short has_cmb_polarization_b;       /* Do we need to compute spectra or bispectra for the CMB B-modes? */
  short has_pk_matter;                /* Do we need the second-order matter Fourier spectrum? */

  /* TODO: I am not sure we really need these two */
  short has_cls;                      /* Do we need any harmonic space spectrum C_l?*/
  short has_bispectra;                /* Do we need any harmonic space bispectrum */



  // ========================================================================================
  // =                                Initial conditions                                    =
  // ========================================================================================

  /* Initial condition flags */
  short has_ad;                      /* Do we need adiabatic mode? */
  short has_ad_first_order;          /* Use first-order initial conditions for the adiabatic modes */
  short has_zero_ic;                 /* Do we need adiabatic mode, with vanishing initial conditions? */  
  short has_unphysical_ic;           /* Custom initial conditions to be specified in perturb2_intial_conditions */


  /* Primordial non-Gaussianity */
  double primordial_local_fnl_phi;         /* Amount of primordial non-Gaussianity of the local type */


  /* Index pointing to the first-order initial conditions needed to solve the second order system. For the
  time being it is set to adiabatic initial conditions only. */
  int index_ic_first_order;



  // =======================================================================================
  // =                                Differential system                                  =
  // =======================================================================================
  

  // *** Option flags for the differential system (initialized in the input module)
  short has_polarization2;                  /* Shall we evolve the photon polarization hierarchy at second-order? */  
  short has_quadratic_sources;              /* Shall we include the quadratic sources in the 2nd-order system at all? */  
  short has_quadratic_liouville;            /* Shall we include the quadratic sources in the Liouville operator? */      
  short has_quadratic_collision;            /* Shall we include the quadratic sources in the photon-baryon collision term? */      
  short has_perfect_baryons;                /* Shall we treat baryons as a pressureless perfect fluid? */
  short has_perfect_cdm;                    /* Shall we treat cold dark matter as a pressureless perfect fluid? */

  short has_perturbed_recombination_stz;    /* Shall we use the perturbed fraction of free electrons? */
  int perturbed_recombination_use_approx;   /* Shall we use the approximation in eq. 3.23 of Senatore et al. 2009? */


  /* Is it ok to evolve the system only up to recombination time?  This is _TRUE_ only if both
  'has_lensing_in_los' and 'has_metric_in_los' are _FALSE_, i.e. if the only contribution to the
  line-of-sight integral comes from visibility function terms.  */
  int has_recombination_only;

  /* Variable that controls whether to use the longitudinal or Poisson equation to compute
  the derivative of the curvature potential (phi_prime) in Newtonian gauge. */
  enum phi_prime_equation phi_prime_eq;


  /* In order to compute the bispectrum integral, it is useful to rescale the line of sight sources
  by a 1/sin(theta_1)^m factor, where theta_1 is the angle between \vec{k1} and \vec{k3}. This
  breaks the (k1,k2) symmetry of the transfer functions, T(k1,k2,k3)=T(k2,k1,k3), which in principle
  is an issue because we compute the transfer functions only or k1>k2. However, since
  k1*sin(theta_1)=k2*sin(theta_2), we can still compute only the transfer functions with k1>k2 and
  obtain the k2>k1 ones using the following identity:
  
  T_rescaled(k2,k1,k3) = T(k2,k1,k3) / sin(theta(k1,k2,k3))^m * (k2/k1)^m
  
  which follows from these other identities:

  T_rescaled(k1,k2,k3) = T(k1,k2,k3) / sin(theta(k1,k2,k3))^m
  T_rescaled(k2,k1,k3) = T(k2,k1,k3) / sin(theta(k2,k1,k3))^m
  T(k1,k2,k3) = T(k2,k1,k3)
  k2 * sin(theta(k2,k1,k3)) = k1 * sin(theta(k1,k2,k3)),
  T_rescaled(k2,k1,k3) = T(k2,k1,k3) / sin(theta(k1,k2,k3))^m * (k2/k1)^m
  
  where T is the symmetric transfer function. */
  short rescale_quadsources;




  // ==========================================================================================
  // =                                 Line of sight sources                                  =
  // ==========================================================================================


  /* Multi-array containing the source functions for each value of k1,k2,k3 and tau. The main task
  of the perturbations2 module is to fill such array.
    
      ppt2->sources [index_type]
                    [index_k1]
                    [index_k2]
                    [index_tau*k3_size + index_k3]

  Due to symmetry properties, 'index_k2' runs from 0 to index_k1, while
  'index_k3' runs from 0 to ppt2->k3_size[index_k1][index_k2]  */
  double **** sources;
  
  /* Logical array. If the index_k1 position is true, then ppt2->sources[index_k1] is allocated */
  short * has_allocated_sources;


  /* For which probes should we compute the sources? These are internal flags used only in 
  the perturbations2 module, contrary to the output flags. */
  short has_source_T;                 /* Do we need source for CMB temperature? */
  short has_source_E;                 /* Do we need source for CMB E-polarization? */
  short has_source_B;                 /* Do we need source for CMB B-polarization? */
  short has_source_g;                 /* Do we need source for gravitational potential? */

  short has_cmb;                      /* Do we need CMB-related sources (temperature, polarization) ? */
  short has_lss;                      /* Do we need LSS-related sources (lensing potential, ...) ? */  

  /* Line of sight sources from the collision term */
  short has_pure_scattering_in_los;           /* Shall we include the purely second-order scattering terms the LOS sources? */
  short has_photon_monopole_in_los;           /* Shall we include the g*I_0_0 term in the line-of-sight sources? */
  short has_quad_scattering_in_los;           /* Shall we include the 'baryon velocity times multipole' terms in the LOS sources? */
  
  /* Line of sight sources from the Liouville term */
  short has_metric_in_los;                    /* Shall we include the metric terms in the LOS sources? */
  short has_quad_metric_in_los;               /* Shall we include the metric*metric terms in the LOS sources? */

  short has_time_delay_in_los;                /* Shall we include the time-delay term for photons in the LOS sources? */
  short has_redshift_in_los;                  /* Shall we include the redshift term for photons in the LOS sources? */
  short has_lensing_in_los;                   /* Shall we include the lensing term for photons in the LOS sources? */

  short use_delta_tilde_in_los;               /* Shall we speficy the LOS sources for Delta-Delta*Delta/2 (see Huang & Vernizzi 2012, Fidler, Pettinari et al. 2014)? */

  /* Sachs-Wolfe (SW) and integrated Sachs-Wolfe (ISW) sources. These sources arise from the integration
  by parts of the purely second-order redshift term. */
  short has_sw;                     /* Shall we include the g*psi Sachs-Wolfe term in the LOS sources? */
  short has_isw;                    /* Shall we include the e^-kappa*( phi'+psi') ISW term in the LOS sources? */
  short use_exponential_potentials; /* Use exponential potentials in constructing the line of sight sources? (see eq. 3.21 of my thesis) */
  short only_early_isw;             /* Include only the early ISW, no late ISW */


  /* ~~~~~~~~~~~~~       INTEGRATION BY PARTS        ~~~~~~~~~~~~~~~~
   *
   * The following discussion is based on sec. 4.2 of http://arxiv.org/abs/1302.0832.
   * 
   *   The source terms forming the line of sight integral can be included as they appear, or they can be integrated
   * by parts.  Applying integration by parts to a term in the source S_l,m generates two terms
   * involving the source on the previous multipole, S_(l-1),m. The first of such terms is peaked at recombination, while the
   * second one consists of the time derivative of the original term. The advantage of the technique is evident when
   * the time derivative of the original term is very small. This is the case of the split between the Sachs-Wolfe (g*psi)
   * and integrated Sachs-Wolfe effects (e^-kappa (phi' + psi')), with the second vanishing for most of the
   * evolution of the Universe. This allows to split nicely the ISW effect into an early part (that cannot be separated
   * from the SW effect) and a late ISW, which is important in the presence of, for instance, dark energy.
   * 
   *   It is important to note that integration by parts does not affect the final result, if one integrates the
   * sources up to today. Only in this case the surface term vanishes.
   * 
   *   At second order, we have many more terms in the line of sight sources. Their time derivative does not
   * vanish at later times, not even in the case of the gravitational potentials (see, for example, eqs. 2.4
   * and 2.5 of http://iopscience.iop.org/1475-7516/2009/08/029/). Performing integration by parts is still
   * useful if one adopts the flat-sky approximation to compute the bispectrum (as done in CMBQuick) because
   * in that formalism it is important that all effects are peaked at recombination.
   * 
   *   The has_integration_by_parts_of_los flag is completely dependent on the flags ppt2->has_sw
   * and ppt2->has_integration_by_parts_of_los, and it is set in the input module.
   */


	/* If true, use an hard-coded test source */
	short use_test_source;

  /* Indices running on types (temperature, polarization, lensing, ...) */
  int index_tp2_T;                    /**< Index value for photon temperature */
  int index_tp2_E;                    /**< Index value for photon E-polarization */
  int index_tp2_B;                    /**< Index value for photon B-polarization */
  int index_tp2_g;                    /**< Index value for gravitational potential */
  int n_sources_T;                    /**< Number of sources to be computed for photon temperature */
  int n_sources_E;                    /**< Number of sources to be computed for photon E-polarization */
  int n_sources_B;                    /**< Number of sources to be computed for photon B-polarization */
  int tp2_size;                       /**< Number of source types that we need to compute */

  /**< Array of strings that contain the labels of the various source types
  For example, tp2_labels[index_tp2_phi] is equal to the string "phi" */
  char ** tp2_labels;

  int n_nonzero_sources_E;   /**< Number of nonzero sources to be computed for photon E-polarization; this is basically
                                  ppt2->n_sources_E minus the l=0 and l=1 modes */
  int n_nonzero_sources_B;   /**< Number of nonzero sources to be computed for photon B-polarization; this is basically
                                  ppt2->n_sources_B minus the l=0, l=1 and m=0 modes */
  
  // ===============================================================================
  // =                                  CMB fields                                 =
  // =============================================================================== 

  /* Indices running on the CMB fields (temperature, E-modes, B-modes, Rayleigh...) considered
  in the perturbation module. We will define similar indices also in the bispectrum module. It
  might seem as a redundancy, but we need to do so as in the polarisation module one has consider
  both types of polarisation at all times (E and B mix), while in the bispectrum module one can
  ask for, e.g., only the E-mode polarisation bispectrum. */
  int index_pf_t;
  int index_pf_e;
  int index_pf_b;
  int pf_size;

  /* Parity of the considered field. Even parity (T,E) is represented by zero, odd parity (B)
  by 1. Indexed as field_parity[index_pf] */
  int field_parity[_MAX_NUM_FIELDS_];

  /* Array of strings that contain the text labels of the various fields */
  char pf_labels[_MAX_NUM_FIELDS_][_MAX_LENGTH_LABEL_]; /* T,E,B... */


  // ==============================================================================
  // =                                  Multipoles                                =
  // ==============================================================================


  /*  Array containing the 'm' values for which we need to solve the system.
     m=0 -> scalar modes
     m=1 -> vector modes
     m=2 -> tensor modes
     m>2 -> modes that do not couple with the metric.  */
  int * m;
  int m_size;
  
    
  /* Flags to denote whether ppt2->m contains m=0, m=1 or m=2 */
  short has_scalars;
  short has_vectors;
  short has_tensors;


  /* Correspondance between (type,l,m) and index_tt */
  int * corresponding_l;         /* corresponding_l[index_tp2] is the multipole 'l' corresponding to index_tp2 */
  int * corresponding_index_m;   /* corresponding_index_m[index_tp2] is the index in ppt2->m corresponding to index_tp2 */
  int * index_monopole;          /* index_monopole[index_tp2] is the index in ppt2->sources of the monopole corresponding to index_tp2 */


  /* 'lm_extra' is an index of the complexity of the angular dependence of the qudratic
  terms in Boltzmann equation.  It is equal to 1 for Newtonian gauge and 3 for
  synchronous gauge.  It is defined so that the (L,M) multipole of the quadratic terms 
  depends on multipoles with (L+-lm_extra, M+-lm_extra).  We use 'lm_extra' to fill and
  index the rotation vectors ppw2->rotation_1 and ppw2->rotation_2. */
  int lm_extra;
  
  /* 'largest_l' is the maximum 'l' requested for the various species at second-order, as
  chosen by the user in the parameter files. */
  int largest_l;
  
  /* 'largest_l_quad' is the highest l-index we shall ever need to use.  This is given by
  whatever is bigger between the number of multipoles evolved at first order, and 'largest_l'
  plus the 'lm_extra' contribution (see above). */
  int largest_l_quad;

  /* 'lm_array' contains the index associated with a given (l,m) couple. The idea is that
  the (l,m) multipole should be addressed from its array, say 'y' as
       y[index_monopole + lm_array[l][index_m]]. */
  int ** lm_array;

  /* 'nlm_array' contains the index associated with a given (n,l,m) triplet. The idea is that
  the (n,l,m) multipole should be addressed from its array, say 'y' as
       y[index_monopole + nlm_array[n][l][index_m]]. */
  int *** nlm_array;


  /* 'lm_array_quad' has the same use of 'lm_array', but it runs up to 'l_max+lm_extra'
  and 'm_max + lm_extra'.  It is used to index the rotation coefficients. */
  int ** lm_array_quad;  

  /* Coupling coefficients C and D (see eq. 141 of Beneke and Fidler 2010) */
  double ** c_minus;
  double ** c_plus;
  double ** d_minus;
  double ** d_plus;
  double ** d_zero;
  
  /* Coupling coefficients in multipole space, given approximately by the product of two Clebsch-Gordan
  symbols. Indexed as ppt2->coupling_coefficients[index_pf][lm(l,m)][l1][m1+ppt2->l1_max][l2]. They
  vanish for configurations where the triangular inequality between l3, l1 and l2 is not met.
  We define separate coefficients for temperature and polarization because the latter is spin-dependent
  (compare Eqs. 3.6, 3.7 and 3.9 of arXiv:1401.3296). */
  double ***** coupling_coefficients;

  /* Maximum values of l1 and l2 where we have computed the Gaunt coefficients. They are
  both determined by ppr2->l_max_los_quadratic and are used the truncate the summation
  for the delta_tilde transformation. */
  int l1_max;
  int l2_max;
  

  // ===============================================================================
  // =                                 k-sampling                                  =
  // ===============================================================================


  /* At second order, we solve the Boltzmann-Einstein system of equations for
  a 3D grid of parameters given by k1 (magnitude of k1 wavevector), k2 (magnitude
  of k2 wavevector) and k3 (magnitude the k1+k2 wavevectors). */

  /* The values that k1 and k2 can assume are controlled by the k_min, k_max,
  k_size parameters that are specified by the user in the .ini file. Both k1
  and k2 take their values from the ppt2->k[index_k] vector that
  is logarithimically sampled from k_min, k_max, k_size.  k2, however, is
  defined so that it can only be larger or equal than k1:  k2>=k1, as the
  equations respect the k1<->k2 symmetry.  This is why we have defined two
  separate ppt2->k1 and ppt2->k2 vectors, instead of using ppt2->k for both
  k1 and k2. */
  
  enum sources2_k_sampling k_sampling;        /* lin, log or default CLASS sampling for ppt2->k? */
  double * k;                                 /* Array containing the magnitudes of the k1 and k2 wavemodes */
  int k_size;                                 /* Size of ppt2->k for what concerns the k1 and k2 sampling */
                                                                  
                                                                  
  enum sources2_k3_sampling k3_sampling;      /* lin, log or default CLASS sampling for ppt2->k3? */
  double *** k3;
  int ** k3_size;
  
  /* Index in ppt2->k corresponding to ppt2->k3[index_k1][index_k2]. It is allocated only
  when the k3 sampling is set to smart. */
  int ** index_k3_min;




  // ===================================================================================
  // =                                 Time sampling                                   =
  // ===================================================================================


  /* Vector that contains the time values where the line-of-sight sources will be computed */
  double * tau_sampling;
  int tau_size;

  /* Time at which the second-order system will start being evolved */
  double tau_start_evolution;

  /* Value of g/g(tau_rec) when to stop sampling the line of sight sources, where g is the
  visibility function. For example, if set to 100, then the last conformal time where
  we will sample the sources will satisfy g(tau)/g(tau_rec)=100. This parameter is overridden
  when the user asks for ISW or other late-time effects, because in that case the sampling
  goes all the way to today */
  double recombination_max_to_end_ratio;

  /* Should we compute the 2nd-order sources at the times specified by the user? */
  short  has_custom_timesampling;          
  double custom_tau_ini;                                    /* Initial time for the custom sampling of the source terms */
  double custom_tau_end;                                    /* Final time for the custom sampling of the source terms */
  int    custom_tau_size;                                      /* Number of points where to sample the source terms (custom timesampling) */
  enum   sources_tau_samplings custom_tau_mode;              /* lin, log or default CLASS sampling for ppt2->sources? */


  /* Should we adjsut the time sampling of the 1st-order line-of-sight sources to the one of the
  2nd-order sources? If _TRUE_, the first-order Cl's will include line-of-sight effects up to
  the same time specified for the 2nd-order sources. */
  short match_final_time_los;
    




  // =========================================================================================
  // =                                    Approximations                                     =
  // =========================================================================================

  
  /* For simplicity, we store the 2nd-order approximation flags here (in the ppt2 structure) rather
  than in the precision one.  (This is not the case for the first-order perturbations in CLASS).
  The following parameters are documented in the declaration of the ppr structure inside the
  file common.h */
  int tight_coupling_approximation;
  double tight_coupling_trigger_tau_c_over_tau_h;
  double tight_coupling_trigger_tau_c_over_tau_k;

  int radiation_streaming_approximation;
  double radiation_streaming_trigger_tau_over_tau_k;

  int ur_fluid_approximation;
  double ur_fluid_trigger_tau_over_tau_k;

  int no_radiation_approximation;
  double no_radiation_approximation_rho_m_over_rho_r;






  // =================================================================================
  // =                        Storage of intermediate results                        =
  // =================================================================================

  /* Files where the sources will be stored (one file for each k1) */
  char sources_dir[_FILENAMESIZE_];
  FILE ** sources_files;
  char ** sources_paths;

  /* File that will keep track how how many sources files have been succesfully written */
  FILE * sources_status_file;
  char sources_status_path[_FILENAMESIZE_];






  // ====================================================================================
  // =                                 Debug parameters                                 =
  // ====================================================================================

  ErrorMsg error_message; /**< String where to write error messages */
  short perturbations2_verbose; /**< Flag regulating the amount of information sent to standard output (none if set to zero) */

  long int count_allocated_sources;   /**< Number of allocated entries in ppt2->sources */
  long int count_memorised_sources;   /**< Number of used entries of ppt2->sources */
  
  long int count_k_configurations;   /**< Number of k-modes for which we shall solve the differential system */

  short has_early_transfers1_only; /**< If _TRUE_, SONG will compute only the first-order early transfer functions and
                                        do not care about the other flags. Useful for debugging. */
  short has_early_transfers2_only; /**< If _TRUE_, SONG will compute only the second-order early transfer functions and
                                        do not care about the other flags. Useful for debugging. */


  /* - Parameters related to the creation of debug files */
  short has_debug_files;        /**< Shall we dump to file the intermediate results such as evolved transfer functions and quadratic sources? */

  char transfers_filename[_FILENAMESIZE_];       /**< Path to the file that will contain the early transfer functions for the different species */
  char quadsources_filename[_FILENAMESIZE_];     /**< Path to the file that will contain the quadratic sources for the different species */
  char quadliouville_filename[_FILENAMESIZE_];   /**< Path to the file that will contain the quadratic sources of the Liouville perator for the different species */
  char quadcollision_filename[_FILENAMESIZE_];   /**< Path to the file that will contain the quadratic sources of the collision term for the different species */

  FILE * transfers_file;         /**< File that will contain the early transfer functions for the different species */
  FILE * quadsources_file;       /**< File that will contain the quadratic sources for the different species */ 
  FILE * quadliouville_file;     /**< File that will contain the quadratic sources of the Liouville perator for the different species */ 
  FILE * quadcollision_file;     /**< File that will contain the quadratic sources of the collision term for the different species */ 

  int index_k1_debug;      /**< Which k1 value should we dump to file? */
  int index_k2_debug;      /**< Which k2 value should we dump to file? */
  int index_k3_debug;      /**< Which k3 value should we dump to file? */
  int l_max_debug;         /**< For any hierarchy, how many 'l' values should we dump to file? */

};




/**
 * Workspace containing, among other things, the value at a given time
 * of all background/perturbed quantitites, as well as their indices.
 *
 * The workspace has to be thought as containing everything that is needed
 * for the evolution of a given wavemode-set (k1,k2,k3). There will
 * be one such structure for each thread (in case of parallel computing).
 */
struct perturb2_workspace 
{


  // ***************       Geometrical Variables       *********************
  
  /* Wavemodes for which we are solving the differential system. */
  int index_k1;
  int index_k2;
  int index_k3;
  double k1, k2;
  
  /* Norm of the k3 vector, which is the sum of the k1 and k2 wavemodes. Inside the
    differential system, we shall call it simply 'k' */
  double k, k_sq;
  
  /* Cosine of the angles between k1 and k, between k2 and k, and between k1 and k2. These are
  obtained assuming the vector k is aligned with the zenith (z axis) */
  double cosk1k, cosk2k, cosk1k2, sink1k, sink2k;
  
  /* Angles between k and k1, and k and k2 */
  double theta_1, theta_2;

  /* Scalar products between the three wavemodes */
  double k1_dot_k2;                    /* scalarProduct(k1,k2) */
  double k1_dot_k2_rescaled;

  /* Spherical components of the wavemodes.  These are obtained using k[m] = xi[m]_i k^i,
  where xi[m]^i is the spherical basis (see Beneke & Fidler 2010, eqs. 68 and A.12) */
  double k1_m[3];                           /* The m-th component is given by k1_m[m+1] */
  double k2_m[3];                           /* The m-th component is given by k2_m[m+1] */

  /* Tensorial product between k1 and k2, that is, X[0]^ij k1_i k2_j, where
  X[0]^ij is the projection matrix.  It appears in the (2,2,m) hierachy for
  cold matter, and in general whenever one wants to convert tensor equations
  into m=0,1,2 multipoles. It is also present when converting beta-moments to
  fluid variables. Use k1_ten_k2[m+2] to access the m-th component. */
  double k1_ten_k2[5];
  double k1_ten_k1[5];
  double k2_ten_k2[5];


  /* Array containing the rotation coefficients needed to relate the first-order quantities
   (computed in the module perturbations.c) to our second-order system.  These factors depend
   on the lm multipole considered and on the cosine of the angle between k and either k1 or k2. */
  double *rotation_1;                         /* Rotation factors for k1 */
  double *rotation_2;                         /* Rotation factors for k2 */
  double *rotation_1_minus;                   /* Rotation factors for k1 (negative m's) */
  double *rotation_2_minus;                   /* Rotation factors for k2 (negative m's) */



  // *** Inner products weighted by the coupling coefficients

  /* Sum over 'm' of the coupling coefficients (c_minus, c_plus, c_zero, etc.) with the rotation
    coefficients defined above.  It is not strictly necessary to precompute these product arrays,
    but it saves a lot of computational time as they do not depend on time and can be computed
    once for each wavemode-set we evolve */
  
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



  // ******************            Indices for the metric          *********************

  /* All possible useful indices for those metric perturbations which are not integrated
  over time, but just inferred from Einstein equations (_mt_" stands for "metric") */

  /* Newtonian gauge */
  int index_mt2_psi;                      /* psi in newtonian gauge */
  int index_mt2_phi_prime;                /* (d phi/d tau) in newtonian gauge, will set to be equal to one of the below indices. */
  int index_mt2_phi_prime_poisson;        /* (d phi/d tau) in newtonian gauge, using the Poisson equation */
  int index_mt2_phi_prime_longitudinal;   /* (d phi/d tau) in newtonian gauge, using the longitudinal equation */
  int index_mt2_omega_m1_prime;           /* vector mode of the metric in Newtonian gauge */
  int index_mt2_gamma_m2_prime_prime;     /* tensor mode of the metric in Newtonian gauge */           


  /* Synchronous gauge */
  int index_mt2_h_prime;     /**< h' (wrt conf. time) in synchronous gauge */
  int index_mt2_h_prime_prime; /**< h'' (wrt conf. time) in synchronous gauge */
  int index_mt2_eta_prime;   /**< eta' (wrt conf. time) in synchronous gauge */
  int index_mt2_alpha_prime; /**< (d \f$ \alpha \f$/d conf.time) in synchronous gauge, where \f$ \alpha = (h' + 6 \eta') / (2 k^2) \f$ */

  int mt2_size;              /**< size of metric perturbation vector */

 



  // ******************            Indices for the quadratic sources          *********************

  int qs2_size;

  /* Quadratic sources for the metric, Newtonian gauge */
  int index_qs2_psi;
  int index_qs2_psi_prime;
  int index_qs2_phi_prime;               
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

  /* Other useful quadratic sources. */
  int index_qs2_dd_b;    /* Quadratic density of baryons */
  int index_qs2_vv_b;    /* Quadratic velocity of baryons */
  int index_qs2_vv_cdm;  /* Quadratic velocity of CDM */

  /* Constants needed to assign the indices */
  int l_max_g;
  int n_hierarchy_g;
  int l_max_pol_g;
  int n_hierarchy_pol_g;
  int n_hierarchy_b;
  int n_hierarchy_cdm;    
  int l_max_ur;
  int n_hierarchy_ur;
  
  
  // ******************            Indices for the source terms          *********************
  
  /* Indices for the ppw2->source_term_table, used to perform the integration by parts of the line-of-sight sources.
  We shall temporarily store in ppw2->source_term_table the results obtained from the differential system.
  We cannot use the ppt2->sources array for this task, as we don't want to keep these intermediate results. The
  advantage of using this buffer array ppw2->source_term_table is that we can obtain its derivatives once the evolution
  of the system is over. */
   
  int st2_size;                       /* Number of variables stored in pptw->source_term_table */ 

  double * source_term_table;         /* Indexed as source_term_table[index_tau*ppw2->st2_size+index_st2] */


  // ******************            Tables for the quadratic sources          *********************

  double ** quadsources_table;

  /* These tables will contain the second-derivatives of the aboves, in view of spline interpolation */
  double ** dd_quadsources_table;



  // *****************         Arrays Filled at Each Evolution Step         ************************

  double * pvecback;                           /* Background quantitites */
  double * pvecthermo;                         /* Thermodynamics quantitites */
  double * pvecmetric;                         /* Metric quantitites */
  double * pvec_quadsources;                   /* Quadratic sources */
  double * pvec_quadcollision;                 /* Quadratic collisional sources */

  /* 1st-order perturbations in k1 and k2.  They are filled at each call of
  'perturb2_quadratic_sources' and 'perturb2_initial_conditions' */
  double * pvec_sources1;       
  double * pvec_sources2;

  /* Pointer to vector of integrated perturbations and their time-derivatives. */
  struct perturb2_vector * pv;

  /* The interpolation functions background_at_tau(), thermodynamics_at_z(), perturb_sources_at_tau()
    keep memory of the last point called through the following indices. */
  int inter_mode;
  int last_index_back;
  int last_index_thermo;
  int last_index_sources;



  // ************************  Approximations  *******************************

  /* Approximations used at a given time */

  int index_ap2_tca;         /* Index for tight-coupling approximation */
  int index_ap2_rsa;         /* Index for radiation streaming approximation */
  int index_ap2_ufa;         /* Index for ultra-relativistic fluid approximation */
  int index_ap2_nra;         /* Index for no-radiation approximation */
                              
  int * approx;             /* Array of approximation flags holding at a given time: approx[index_ap] */
  int ap2_size;              /* Number of relevant approximations for a given mode */





  // **********************  Baryon and CDM fluids   **************************
  
  /* Value of the n=2 beta-moments for the baryon and CDM fluids, assuming they are perfect fluids */
  double b_200;
  double b_22m[3];
  double cdm_200;
  double cdm_22m[3];


  // **********************         Debug parameters         *************************

  int derivs_count;   /**< Counter to keep track of how many times the function perturb2_derivs has been called
                           for the considered set of (k1,k2,k3). */

  /** Function used to output intermediate values from the differential system.  This
  function will be given as an argument to 'generic_evolver' and is  used only for debug
  purposes.  By default it is set to NULL, which means that it is never called.  It is
  called only if ppt2->has_debug_files==_TRUE_ and if the evolved wavemode corresponds
  to the one requested through ppt2->index_k1_debug, ppt2->index_k2_debug and
  ppt2->index_k3_debug. This function will be called for each time step in the
  differential system. */
  int (*print_function)(double x, double y[], double dy[], void *parameters_and_workspace, ErrorMsg error_message);
  
  /** String that contains information on the wavemode that is being currently integrated.
  such information is printed to the debug files if ppt2->perturbations2_verbose is high
  enough */
  char info [4096];

  long int n_steps;   /**< Number of steps taken by the differential system so far for the active (k1,k2,k3) mode.
                           Computed only if has_debug_files==_TRUE_ */

};





/**
 * Structure containing the indices and the values of the perturbation
 * variables which are integrated over time (as well as their
 * time-derivatives). For a given wavenumber, the size of these
 * vectors changes when the approximation scheme changes.
 *
 * The vector structure is what is filled and read by the differential
 * evolver. It is also used to set the initial conditions and to set
 * which equations to evolve.
 *
 */

struct perturb2_vector
{

  // ================================================
  // =                    Shortcuts                 =
  // ================================================

  /* Note that the number of evolved equations changes during the evaluation of a single (k1,k2,k3) mode.
  This happens, for example, whenever the approximation scheme changes at some time.  Hence, the following
  'l_max' and 'n_hierarchy' variables are not related to the particular (k1,k2,k3) mode, but to the
  time interval. */
  int l_max_g;                  /* Max 'l' momentum in Boltzmann hierarchy */
  int l_max_pol_g;              /* Max 'l' momentum in Boltzmann hierarchy */
  int l_max_ur;                 /* Max 'l' momentum in Boltzmann hierarchy */
  int n_hierarchy_g;            /* Number of photon temperature hierarchy equations */
  int n_hierarchy_pol_g;        /* Number of photon polarization hierarchy equations (same for E and B-modes) */  
  int n_hierarchy_b;            /* Number of baryon hierarchy equations */
  int n_hierarchy_cdm;          /* Number of cold dark matter hierarchy equations */    
  int n_hierarchy_ur;           /* Number of neutrino hierarchy equations */

  /* Needed to address the baryon and CDM hierarchies with a loop over n,l,m.
  For perfect fluids, n_max=l_max=1, otherwise =2. */
  int n_max_b;                  
  int l_max_b;
  int n_max_cdm;
  int l_max_cdm;



  // ==============================================
  // = Indices of evolved second-order quantities = 
  // ==============================================

  // ***** Photon hierarchies
  int index_pt2_monopole_g;                   /* Photon temperature hierarchy starts here */
  int index_pt2_monopole_E;                   /* Photon E-mode polarization hierarchy starts here */  
  int index_pt2_monopole_B;                   /* Photon B-mode polarization hierarchy starts here */    

  // *** Baryons hierarchy
  int index_pt2_monopole_b;

  // ***** Cold dark matter hierarchy
  int index_pt2_monopole_cdm;

  // ***** Neutrino hierarchy
  int index_pt2_monopole_ur;
      
  // ***** Metric variables
  int index_pt2_eta;                        /* Synchronous gauge metric perturbation eta */
  int index_pt2_phi;                        /* Newtonian gauge potential phi */
  int index_pt2_omega_m1;                   /* Newtonian gauge vector potential omega_[m=1] */
  int index_pt2_gamma_m2;                   /* Newtonian gauge tensor potential gamma_[m=2] */
  int index_pt2_gamma_m2_prime;             /* the equation for gamma_[m=2] is a second order DE */

  int pt2_size;                             /* Size of perturbation vector, equivalent to the number of evolved quantities */

  // =================
  // =    Arrays     =
  // =================

  /* Array of perturbations to be integrated. It is filled at the end of each iteration of the evolver. */
  double * y;

  /* Vector containing the time-derivative of the pv->y. It is filled by each call of perturb2_derivs. */
  double * dy;               

  /* Boolean array specifying which perturbations enter in the calculation of the source functions. Only
    the marked perturbations will be interpolated at the times requested in ppt2->tau_sampling. */
  int * used_in_sources;


  // =================
  // =     Misc      =
  // =================

 
};




/**
 * Structure pointing towards all what the function that perturb2_derivs
 * needs to know: fixed input parameters and indices contained in the
 * various structures, workspace, etc.
*/ 
struct perturb2_parameters_and_workspace {

  struct precision * ppr;                         /* Pointer to the precision structure */
  struct precision2 * ppr2;                       /* Pointer to the precision2 structure */
  struct background * pba;                        /* Pointer to the background structure */
  struct thermo * pth;                            /* Pointer to the thermodynamics structure */
  struct perturbs * ppt;                          /* Pointer to the perturbation structure */
  struct perturbs2 * ppt2;                        /* Pointer to the 2nd-order perturbation structure */
  struct perturb2_workspace * ppw2;               /* Worspace defined above */
    
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
          int * pa_old
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

    int perturb2_save_early_transfers(double tau,
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
          int index_tau, /* if negative, use interpolation with tau, below */
          double tau,
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
            struct perturb2_workspace * ppw2
            );

    int perturb2_quadratic_sources_at_tau_linear(
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double tau,
            struct perturb2_workspace * ppw2
            );

    int perturb2_quadratic_sources_at_tau_spline(
            struct perturbs * ppt,
            struct perturbs2 * ppt2,
            double tau,
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

    int perturb2_store_sources_to_disk(
            struct perturbs2 * ppt2,
            int index_k1
            );

    int perturb2_load_sources_from_disk(
            struct perturbs2 * ppt2,
            int index_k1
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