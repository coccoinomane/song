/** @file print_background.c 
 * Created by Guido W. Pettinari on 01.09.2011
 * Last edited by Guido W. Pettinari on 01.09.2011
 *
 * Print to screen the quantities derived in the background module.
 *
 * usage:     print_background <ini file> <pre file>
 *
 * These are the quantities that will be printed to screen.
 *    
 *    1 CONFORMAL TIME TAU
 *    2 int index_bg_a;             // scale factor 
 *    3 int index_bg_H;             // Hubble parameter in Mpc^{-1} 
 *    4 int index_bg_H_prime;       // its derivative w.r.t. conformal time 
 *    5 int index_bg_rho_g;         // photon density 
 *    6 int index_bg_rho_b;         // baryon density 
 *    7 int index_bg_rho_cdm;       // cdm density 
 *    8 int index_bg_rho_lambda;    // cosmological constant density 
 *    * int index_bg_rho_fld;       // fluid with constant w density 
 *    9 int index_bg_rho_ur;        // relativistic neutrinos/relics density 
 *                                
 *    * int index_bg_rho_ncdm1;     // density of first ncdm species (others contiguous) 
 *    * int index_bg_p_ncdm1;       // pressure of first ncdm species (others contiguous) 
 *    * int index_bg_pseudo_p_ncdm1;// another statistical momentum useful in ncdma approximation 
 *                                
 *    10 int index_bg_Omega_r;       // relativistic density fraction (\f$ \Omega_{\gamma} + \Omega_{\nu r} \f$) 
 *    11 int index_bg_rho_crit;      // critical density 
 *    12 int index_bg_Omega_m;       // non-relativistic density fraction (\f$ \Omega_b + \Omega_cdm + \Omega_{\nu nr} \f$) 
 *    13 int index_bg_conf_distance; // conformal distance (from us) 
 *    14 int index_bg_time;          // proper (cosmological) time 
 *    15 int index_bg_rs;            // comoving sound horizon 
 *    16 Hc
 *    17 log10(a/a_eq)
 * 
 *  Copy and paste:
 *  tau a  H  H_prime  rho_g  rho_b  rho_cdm  rho_lambda  rho_ur  Omega_r  rho_crit  Omega_m  conf_distance  ang_distance  lum_distance  time  rs Hc y
 */
 
 
 
 

 
 

#include "class.h"


int main(int argc, char **argv) {

  struct precision pr;        /* precision parameters (1st-order) */
  struct precision2 pr2;      /* precision parameters (2nd-order) */
  struct background ba;       /* cosmological background */
  struct thermo th;           /* thermodynamics */
  struct perturbs pt;         /* source functions (1st-order) */
  struct perturbs2 pt2;       /* source functions (2nd-order) */  
  struct bessels bs;          /* bessel functions (1st-order) */
  struct bessels2 bs2;        /* bessel functions (2nd-order) */
  struct transfers tr;        /* transfer functions (1st-order) */
  struct transfers2 tr2;      /* transfer functions (2nd-order) */
  struct primordial pm;       /* primordial spectra */
  struct spectra sp;          /* output spectra (1st-order) */
  struct bispectra bi;        /* bispectra */
  struct fisher fi;           /* fisher matrix */
  struct nonlinear nl;        /* non-linear spectra */
  struct lensing le;          /* lensed spectra */
  struct output op;           /* output files */
  ErrorMsg errmsg;            /* error messages */

  // ===================
  // = Parse arguments =
  // ===================
  if (argc > 3) {
    printf ("usage:     %s <ini file> [<pre file>]\n", argv[0]);
    return _FAILURE_;
  }


  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }



  // =======================================
  // = Compute background & thermodynamics =
  // =======================================
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }


  /* Some useful information */
  printf ("a_equality = %g\n", ba.a_eq);

  // Time sampling
  int tau_size = ba.bt_size;
  double * tau_sampling = ba.tau_table;
  
  // Running indices
  int index_tau;
  int index_bg;           // Index on the type of background quantity
  double tau, var, a;
  
  // Vector that will contain background quantities (we are only interested in 'a')
  int dump;
  double * pvecback = malloc(ba.bg_size*sizeof(double));
  // double * pvecthermo = malloc(ba.bg_size_short*sizeof(double));  

  
  // ==============================
  // = Print the background table =
  // ==============================

  /* Loop on time */
  for (index_tau = 0; index_tau < tau_size; ++index_tau) {
    
    /* Take only a limited number of points */
    if ((index_tau%4!=0)&&(index_tau!=tau_size-1)) continue;
  
    tau = tau_sampling[index_tau];
  
    // Find H in conformal time, i.e. a_prime_over_a
    double a = ba.background_table[index_tau*ba.bg_size+ba.index_bg_a];
    double a_prime_over_a = a*ba.background_table[index_tau*ba.bg_size+ba.index_bg_H];
        
    // First column is time.  The plus indicates that even plus signs will be printed.
    fprintf (stderr, "%+6e ", tau);
  
    // Columns from 2 to ba.bg_size+1 are the background quantities
    for (index_bg = 0; index_bg < ba.bg_size; ++index_bg) {
     
      var = ba.background_table[index_tau*ba.bg_size+index_bg];
         
      fprintf (stderr, "%+6e ", var);
    }
    
    // The column before the last is a_prime_over_a
    fprintf (stderr, "%+6e ", a_prime_over_a);
    
    // Last column is the scale factor normalized to equality, i.e. y as used by Pitrou et al.
    fprintf (stderr, "%+6e\n", log10(a/ba.a_eq));
  }
  
  
  /****** all calculations done, now free the structures ******/
  
  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }
  
  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
