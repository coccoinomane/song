/** @file thermodynamics.h Documented includes for thermodynamics module */

#ifndef __THERMODYNAMICS__
#define __THERMODYNAMICS__

#include "background.h"
//#include "arrays.h"
//#include "helium.h"
//#include "hydrogen.h"

/**
 * List of possible recombination algorithms.
 */

enum recombination_algorithm {
  recfast, 
  hyrec
};

/**
 * List of possible reionization schemes.
 */

enum reionization_parametrization {
  reio_none, /**< no reionization */
  reio_camb,  /**< reionization parameterized like in CAMB */
  reio_bins_tanh  /**< binned reionization history with tanh inteprolation between bins */ 
};

/**
 * Is the input parameter the reionization redshift or optical depth?
 */

enum reionization_z_or_tau {
  reio_z,  /**< input = redshift */
  reio_tau /**< input = tau */
};

/**
 * Two useful smooth step functions, for smoothing transitions in recfast.
 */

#define f1(x) (-0.75*x*(x*x/3.-1.)+0.5)  /* goes from 0 to 1 when x goes from -1 to 1 */
#define f2(x) (x*x*(0.5-x/3.)*6.)        /* goes from 0 to 1 when x goes from  0 to 1 */

/**
 * All thermodynamics parameters and evolution that other modules need to know.
 *
 * Once initialized by thermodynamics_init(), contains all the
 * necessary information on the thermodynamics, and in particular, a
 * table of thermodynamical quantities as a function of the redshift,
 * used for interpolation in other modules.
 */

struct thermo 
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantitites are computed in this module, given these parameters
   *   and the content of the 'precision' and 'background' structures) */
  
  //@{

  double YHe;  /**< \f$ Y_{He} \f$ : primordial helium fraction */

  enum recombination_algorithm recombination; /**< reionization code */

  enum reionization_parametrization reio_parametrization; /**< reionization scheme */

  enum reionization_z_or_tau reio_z_or_tau; /**< is the input parameter the reionization redshift or optical depth? */

  double tau_reio; /**< if above set to tau, input value of reionization optical depth */

  double z_reio;   /**< if above set to z,   input value of reionization redshift */

  short compute_cb2_derivatives; /**< do we want to include in computation derivatives of baryon sound speed? */

  // *** MY MODIFICATIONS ***
  short has_rayleigh_scattering;  /** Do we want to compute the interaction rate and visibility function of Rayleigh scattering?
                                      Note that the frequency dependence is not included at this stage, so that the only
                                      difference with the usual interaction rate is that the number density of neutral 
                                      Hydrogen is used instead of that of free electrons. */

  double rayleigh_frequency;                /**< frequency where to compute the Rayleigh scattering */
  
  // *** END OF MY MODIFICATIONS ***
    

  /** parameters for reio_camb */

  double reionization_width; /**< width of H reionization */

  double reionization_exponent; /**< shape of H reionization */

  double helium_fullreio_redshift; /**< redshift for of helium reionization */

  double helium_fullreio_width; /**< width of helium reionization */

  /** parameters for reio_bins_tanh */

  int binned_reio_num; /**< with how many bins de we want to describe reionization? */

  double * binned_reio_z; /**< central z value for each bin */

  double * binned_reio_xe; /**< imposed x_e(z) value at center of each bin */

  double binned_reio_step_sharpness; /**< sharpness of tanh() step interpolating between binned values */

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities stored in table */

  //@{

  int index_th_xe;            /**< ionization fraction \f$ x_e \f$ */
  int index_th_Tb;            /**< baryon temperature \f$ T_b \f$ */
  int index_th_cb2;           /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_th_dcb2;          /**< derivative wrt conformal time of squared baryon sound speed \f$ d [c_b^2] / d \tau \f$ (only computed if some non-mininmal tight-coupling schemes is requested) */
  int index_th_ddcb2;         /**< second derivative wrt conformal time of squared baryon sound speed  \f$ d^2 [c_b^2] / d \tau^2 \f$ (only computed if some non0-minimal tight-coupling schemes is requested) */
  int index_th_rate;          /**< maximum variation rate of \f$ exp^{-\kappa}, g and (d g / d \tau), used for computing integration step in perturbation module */

  int index_th_dkappa;        /**< total scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_ddkappa;       /**< total scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa;      /**< total scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa;  /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g;             /**< total visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg;           /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */

  // *** MY MODIFICATIONS ***
  int index_th_thomson_dkappa;        /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_thomson_ddkappa;       /**< Thomson scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_thomson_dddkappa;      /**< Thomson scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_thomson_exp_m_kappa;   /**< \f$ exp^{-\kappa} \f$ for Thomson scattering */
  int index_th_thomson_g;             /**< Thomson visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_thomson_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_thomson_ddg;            /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */

  int index_th_rayleigh_dkappa;        /**< Rayleigh scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_rayleigh_ddkappa;       /**< Rayleigh scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_rayleigh_dddkappa;      /**< Rayleigh scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_rayleigh_exp_m_kappa;   /**< \f$ exp^{-\kappa} \f$ for Rayleigh scattering */
  int index_th_rayleigh_g;             /**< Rayleigh visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_rayleigh_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_rayleigh_ddg;            /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  // *** END OF MY MODIFICATIONS ***

  int th_size;                /**< size of thermodynamics vector */ 


  // *** MY MODIFICATIONS ***

  /* The perturbed fraction of free electrons, delta_Xe, is needed by the second-order module. To compute it,
    (see Senatore et al. 2009 or Pitrou et al. 2010) we need the derivatives of a function, Q, defined
    in Senatore et al. 2009. These derivatives are background quantities that can be obtained in
    the thermodynamics module. */
  short has_perturbed_recombination;

  /* From which value of x should we start integrating the system to obtain delta_Xe? Should be larger than 34 to avoid
    the discontinuity in the Q function. Note that x = eps / T, where eps is the ionization energy of the
    hydrogen atom, and T is the temperature of the Universe (T = T_cmb / a). */
  double perturbed_recombination_turnx;
  
  /* Derivatives of the Q function (N.B. following CMBquick example, our Q function is actually Q/n_e
    in Senatore et al. 2009) */
  int index_th_dQ_dx;          /* Derivative with respect to x = H ionization energy / T */
  int index_th_dQ_dX;          /* Derivative with respect to the fraction of free electrons */
  int index_th_dQ_dn;          /* Derivative with respect to the baryon density */
  int index_th_dQ_dH;          /* Derivative with respect to the Hubble factor */
  int index_th_Q;              /* Background value of the Q function */

  double dQ_dx;
  double dQ_dX;
  double dQ_dn;
  double dQ_dH;
  double Q;   

  /* Quantities needed to compute the approximated perturbed recombination (i.e. eq. 3.23 of Senatore et
    al. 2009) */
  int index_th_dxe;            /**< ionization fraction first derivative */
  int index_th_ddxe;           /**< ionization fraction second derivative */  


  /* Do we want to include in computation derivatives of ionization fraction? */
  short compute_xe_derivatives;

  // *** END OF MY MODIFICATIONS ***

  //@}

  /** @name - thermodynamics interpolation tables */

  //@{

  int tt_size; /**< number of lines (redshift steps) in the tables */
  double * z_table; /**< vector z_table[index_z] with values of redshift (vector of size tt_size) */
  double * thermodynamics_table; /**< table thermodynamics_table[index_z*pth->tt_size+pba->index_th] with all other quantities (array of size th_size*tt_size) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2thermodynamics_dz2_table; /**< table d2thermodynamics_dz2_table[index_z*pth->tt_size+pba->index_th] with values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */

  //@}

  /** @name - redshift, conformal time and sound horizon at recombination */

  //@{

  double z_rec;   /**< z at which the visibility reaches its maximum (= recombination redshift) */
  // *** MY MODIFICATIONS ***
  double z_rec_thomson;     /**< z at which the Thomson visibility reaches its maximum */
  // *** END MY MODIFICATIONS ***
  double tau_rec; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec;  /**< comoving sound horizon at recombination */
  double ds_rec;  /**< physical sound horizon at recombination */
  double da_rec;  /**< angular diameter distance to recombination */

  //@}

  /** @name - redshift, conformal time and sound horizon at recombination */

  //@{

  double tau_free_streaming;   /**< minimum value of tau at which free-streaming approximation can be switched on */

  //@}

  /** @name - initial conformal time at which thermodynamical variables have been be integrated */

  //@{

  double tau_ini;

  //@}

/** @name - total number density of electrons today (free or not) */

  //@{

  double n_e;

  //@}

  /** 
   *@name - some flags needed for thermodynamics functions
   */
  
  //@{

  short inter_normal;  /**< flag for calling thermodynamics_at_z and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short thermodynamics_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */ 

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};

/**
 * Temporary structure where all the recombination history is defined and stored.
 *
 * This structure is used internally by the thermodynamics module, 
 * but never passed to other modules.
 */ 

struct recombination {

  /** @name - indices of vector of thermodynamics variables related to recombination */

  //@{

  int index_re_z;          /**< redshift \f$ z \f$ */
  int index_re_xe;         /**< ionization fraction \f$ x_e \f$ */
  int index_re_Tb;         /**< baryon temperature \f$ T_b \f$ */
  int index_re_cb2;        /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_re_dkappadtau; /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  // *** MY MODIFICATIONS ***
  int index_re_rayleigh_dkappadtau; /**< Rayleigh scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) (with no frequency dependence)*/
  // *** END OF MY MODIFICATIONS ***
  int re_size;             /**< size of this vector */

  //@}

  /** @name - table of the above variables at each redshift, and number of redshits */

  //@{

  int rt_size; /**< number of lines (redshift steps) in the table */
  double * recombination_table; /**< table recombination_table[index_z*preco->re_size+index_re] with all other quantities (array of size preco->rt_size*preco->re_size) */

  //@}

  /** @name - recfast parameters needing to be passed to
      thermodynamics_derivs_with_recfast() routine */

  //@{

  double CDB;
  double CR;
  double CK;
  double CL;
  double CT;
  double fHe;
  double CDB_He;
  double CK_He;
  double CL_He;
  double fu;
  double H_frac;
  double Tnow;
  double Nnow;
  double Bfact;
  double CB1;
  double CB1_He1;
  double CB1_He2;
  double H0; 

  //@}

};

/**
 * Temporary structure where all the reionization history is defined and stored.
 *
 * This structure is used internally by the thermodynamics module, 
 * but never passed to other modules.
 */ 

struct reionization {

  /** @name - indices of vector of thermodynamics variables related to reionization */

  //@{

  int index_re_z;          /**< redshift \f$ z \f$ */
  int index_re_xe;         /**< ionization fraction \f$ x_e \f$ */
  int index_re_Tb;         /**< baryon temperature \f$ T_b \f$ */
  int index_re_cb2;        /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_re_dkappadtau; /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_re_dkappadz;   /**< Thomson scattering rate with respect to redshift \f$ d \kappa / d z\f$ (units 1/Mpc) */
  int index_re_d3kappadz3; /**< second derivative of previous quantity with respect to redshift */
  // *** MY MODIFICATIONS ***
  int index_re_rayleigh_dkappadtau; /**< Rayleigh scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) (with no frequency dependence)*/
  // *** END OF MY MODIFICATIONS ***
  int re_size;             /**< size of this vector */

  //@}

  /** @name - table of the above variables at each redshift, and number of redshits */

  //@{

  int rt_size;                 /**< number of lines (redshift steps) in the table */
  double * reionization_table; /**< table reionization_table[index_z*preio->re_size+index_re] with all other quantities (array of size preio->rt_size*preio->re_size) */

  //@}

  /** @name - reionization optical depth infered from reionization history */

  //@{

  double reionization_optical_depth;

  //@}

  /** @name - indices describing input parameters used in the definition of the various possible functions x_e(z) */

  //@{

  /** parameters used by reio_camb */

  int index_reio_redshift;  /**< hydrogen reionization redshift */
  int index_reio_exponent;  /**< an exponent used in the function x_e(z) in the reio_camb scheme */
  int index_reio_width;     /**< a width defining the duration of hydrogen reionization in the reio_camb scheme */
  int index_reio_xe_before; /**< ionization fraction at redshift 'reio_start' */
  int index_reio_xe_after;  /**< ionization fraction after full reionization */
  int index_helium_fullreio_fraction; /**< helium full reionization fraction infered from primordial helium fraction */
  int index_helium_fullreio_redshift; /**< helium full reionization redshift */
  int index_helium_fullreio_width;    /**< a width defining the duration of helium full reionization in the reio_camb scheme */

  /** parameters used by reio_bins_tanh */

  int reio_num_z;
  int index_reio_first_z;
  int index_reio_first_xe;
  int index_reio_step_sharpness;

  /** parameters used by all schemes */

  int index_reio_start;     /**< redshift above which hydrogen reionization neglected */

  //@}

  /** @name - vector of such parameters, and its size */

  double * reionization_parameters;
  int reio_num_params;

  //@}

  /** @name - index of line in recombination table corresponding to first line of reionization table */

  //@{

  int index_reco_when_reio_start;

  //@}

};

/** 
 * temporary  parameters and workspace passed to the thermodynamics_derivs function 
 */

struct thermodynamics_parameters_and_workspace {

  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;
  struct precision * ppr;
  struct recombination * preco;

  /* workspace */
  double * pvecback;

};

/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int thermodynamics_at_z(
        struct background * pba,
        struct thermo * pth,
        double z,
        short inter_mode,
        int * last_index,
        double * pvecback,
        double * pvecthermo
        );

  int thermodynamics_init(
        struct precision * ppr,
        struct background * pba,
        struct thermo * pth
        );

  int thermodynamics_free(
        struct thermo * pthermo
        );

  int thermodynamics_indices(
           struct thermo * pthermo,
           struct recombination * preco,
           struct reionization * preio
           );

  // *** MY MODIFICATIONS ***
  
  /* Compute functions needed to obtain the first-order fraction of free electrons */
  int thermodynamics_compute_Q_derivatives(
            struct precision * ppr,
            struct background * pba,
            struct thermo * pth
            );

  int thermodynamics_compute_Q_derivatives_at_z(
            // struct precision * ppr,
            struct background * pba,
            struct thermo * pth,
            double X_e,
            double z
            );

  // *** MY MODIFICATIONS ***

  int thermodynamics_helium_from_bbn(
             struct precision * ppr,
             struct background * pba,
             struct thermo * pth
             );

  int thermodynamics_reionization_function(
             double z,
             struct thermo * pth,
             struct reionization * preio,
             double * xe
             );

  int thermodynamics_reionization(
          struct precision * ppr,
          struct background * pba,
          struct thermo * pth,
          struct recombination * preco,
          struct reionization * preio,
          double * pvecback
          );

  int thermodynamics_reionization_sample(
           struct precision * ppr,
           struct background * pba,
           struct thermo * pth,
           struct recombination * preco,
           struct reionization * preio,
           double * pvecback
           );
  
  int thermodynamics_get_xe_before_reionization(struct precision * ppr,
            struct thermo * pth,
            struct recombination * preco,
            double z,
            double * xe);

  int thermodynamics_recombination(
           struct precision * ppr,
           struct background * pba,
           struct thermo * pth,
           struct recombination * prec,
           double * pvecback
           );

  int thermodynamics_recombination_with_hyrec(
            struct precision * ppr,
            struct background * pba,
            struct thermo * pth,
            struct recombination * prec,
            double * pvecback
            );

  int thermodynamics_recombination_with_recfast(
            struct precision * ppr,
            struct background * pba,
            struct thermo * pth,
            struct recombination * prec,
            double * pvecback
            );

  int thermodynamics_derivs_with_recfast(
           double z,
           double * y,
           double * dy,
           void * fixed_parameters,
           ErrorMsg error_message
           );

  int thermodynamics_merge_reco_and_reio(
           struct precision * ppr,
           struct thermo * pth,
           struct recombination * preco,
           struct reionization * preio
           );
    
#ifdef __cplusplus
}
#endif

/**************************************************************/

/**
 * @name some flags
 */

//@{

#define _BBN_ -1

//@}

/**  
 * @name Some basic constants needed by RECFAST:
 */

//@{

#define _m_e_ 9.10938215e-31  /**< electron mass in Kg */
#define _m_p_ 1.672621637e-27 /**< proton mass in Kg */
#define _m_H_ 1.673575e-27    /**< Hydrogen mass in Kg */
#define _not4_ 3.9715         /**< Helium to Hydrogen mass ratio */
#define _sigma_ 6.6524616e-29 /**< Thomson cross-section in m^2 */

//@}

/**  
 * @name Some specific constants needed by RECFAST:
 */

//@{

#define _RECFAST_INTEG_SIZE_ 3

#define _Lambda_ 8.2245809
#define _Lambda_He_ 51.3
#define _L_H_ion_ 1.096787737e7
#define _L_H_alpha_ 8.225916453e6
#define _L_He1_ion_ 1.98310772e7
#define _L_He2_ion_ 4.389088863e7
#define _L_He_2s_ 1.66277434e7
#define _L_He_2p_ 1.71134891e7
#define _A2P_s_   1.798287e9     /*updated like in recfast 1.4*/
#define _A2P_t_   177.58e0       /*updated like in recfast 1.4*/
#define _L_He_2Pt_  1.690871466e7  /*updated like in recfast 1.4*/
#define _L_He_2St_  1.5985597526e7 /*updated like in recfast 1.4*/
#define _L_He2St_ion_ 3.8454693845e6 /*updated like in recfast 1.4*/
#define _sigma_He_2Ps_  1.436289e-22   /*updated like in recfast 1.4*/
#define _sigma_He_2Pt_  1.484872e-22   /*updated like in recfast 1.4*/

//@}

/**  
 * @name Some specific constants needed by recfast_derivs:
 */

//@{

#define _a_PPB_ 4.309
#define _b_PPB_ -0.6166
#define _c_PPB_ 0.6703
#define _d_PPB_ 0.5300
#define _T_0_ pow(10.,0.477121)   /* from recfast 1.4 */
#define _a_VF_ pow(10.,-16.744)
#define _b_VF_ 0.711
#define _T_1_ pow(10.,5.114)
#define _a_trip_ pow(10.,-16.306) /* from recfast 1.4 */
#define _b_trip_ 0.761            /* from recfast 1.4 */

//@}

/**  
 * @name Some limits imposed on cosmological parameter values:
 */

//@{

#define _YHE_BIG_ 0.5      /**< maximal \f$ Y_{He} \f$ */
#define _YHE_SMALL_ 0.01   /**< minimal \f$ Y_{He} \f$ */
#define _Z_REC_MAX_ 2000.
#define _Z_REC_MIN_ 500.

//@}



// *** MY MODIFICATIONS ***
/**  
 * @name Some constants needed to compute perturbed recombination, mostly copied from HyRec.
 */

//@{

// #define _EI_   13.598286071938324              /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */
#define _EI_   13.605698                          /* Hydrogen ionization energy in eV, from CMBquick */

/* Energy differences between excited levels of hydrogen -- used often */
#define _E21_  10.198714553953742
#define _E31_  12.087365397278509
#define _E41_  12.748393192442178
#define _E32_  1.8886508433247664
#define _E42_  2.5496786384884356

#define _hPc_       1.239841874331e-04   /* hc in eV cm */
#define _mH_        0.93878299831e9      /* Hydrogen atom mass in eV/c^2 */ 
#define _kBoltz_    8.617343e-5          /* Boltzmann constant in eV/K */
#define _L2s1s_     8.2206               /* 2s -> 1s two-photon decay rate in s^{-1} (Labzowsky et al 2005) */

#define _FINE_             7.2973525698e-3   /* Fine structure constant */
#define _s_in_inverse_eV_  1.51926751e15     /* A second in eV^-1 */
#define _m_in_eV_          5.06773093e6      /* An eV in meters, this is the above divided the speed of light in m/s */
#define _m_e_eV_           5.10998910e5       /* Electron mass in eV */

//@}

/**  
 * @name Physical constant needed to compute the Rayleigh scattering rate
 */

//@{

#define _NU_EFF_ 3.1016927457e6      /**< Effective frequency for the Rayleigh scattering, in Ghz (see arXiv:1307.8148) */

//@}


// *** END OF MY MODIFICATIONS ***



#endif
