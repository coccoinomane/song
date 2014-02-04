/** @file print_background.c 
 * Created by Guido W. Pettinari on 01.09.2011
 * Last edited by Guido W. Pettinari on 01.09.2011
 *
 * Print to screen the quantities derived in the thermodynamics module.
 *
 * usage:     print_thermodynamics <ini file> <pre file>
 *
 * These are the quantities that will be printed to screen.
 *    
 *    1 REDSHIFT Z
 *    2 CONFORMAL TIME TAU
 *    3 SCALE FACTOR A
 *    4 log10(a/a_eq)
 *    5 int index_th_xe;            // ionization fraction \f$ x_e \f$ 
 *    6 int index_th_dxe;           // derivative of ionization fraction wrt to time
 *    7 int index_th_ddxe;          // second derivative of ionization fraction wrt to time
 *    8 int index_th_dkappa;        // Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) 
 *    9 int index_th_ddkappa;       // scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ 
 *    10 int index_th_dddkappa;      // scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ 
 *    11 int index_th_exp_m_kappa;   // \f$ exp^{-\kappa} \f$ 
 *    12 int index_th_g;             // visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ 
 *    13 int index_th_dg;            // visibility function derivative \f$ (d g / d \tau) \f$ 
 *    14 int index_th_ddg;           // visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ 
 *    15 int index_th_Tb;            // baryon temperature \f$ T_b \f$ 
 *    16 int index_th_cb2;           // squared baryon sound speed \f$ c_b^2 \f$ 
 *    17 int index_th_dcb2;          // derivative wrt conformal time of squared baryon sound speed
 *    18 int index_th_ddcb2;         // second derivative wrt conformal time of squared baryon sound speed
 *    19 int index_th_rate;          // maximum variation rate of \f$ exp^{-\kappa}, g and (d g / d \tau)
 * 
 *
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


  // Set verbosity to zero (we want to be able to send the output to a file and plot it)
  ba.background_verbose = 0;
  th.thermodynamics_verbose = 0;

  // We want to print the derivatives of the baryons sound of speed
  th.compute_cb2_derivatives = _TRUE_;

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


  // Time sampling
  int z_size = th.tt_size;
  double * z_sampling = th.z_table;
  
  // Running indices
  int index_z;
  int index_th;           // Index on the type of background quantity
  
  // Vector that will contain background quantities (we are only interested in 'a')
  int dump;
  double * pvecback = malloc(ba.bg_size*sizeof(double));
  double * pvecthermo = malloc(th.th_size*sizeof(double));  


  // =========================================
  // = Find era of matter-radiation equality =
  // =========================================
  //// Using bisection, compute the scale factor at equality.
  double a_equality = ba.a_eq;
  double tol_tau_approx = 1e-6;      // Desired precision in tau
  
  // ==================================
  // = Print the thermodynamics table =
  // ==================================
  // Add rows (loop on redshift)
  for (index_z = 0; index_z < z_size; ++index_z) {
  

    /* Consider only one point every tot */
    if (index_z%64!=0)
      continue;
    
  
    double z = th.z_table[index_z];
  
    // =================================
    // = Get scale factor and tau at z =
    // =================================
    double tau;
    class_call(background_tau_of_z(&ba,	z, &tau),
      th.error_message,
      errmsg);
    
    // Extract values of background quantities at tau(z)
    class_call(background_at_tau(
             &ba,
             tau,
             ba.long_info,
             ba.inter_normal,
             &dump,
             pvecback
             ),
  	     ba.error_message,
  	     errmsg);

    // Value of the scale factor at z
    double a = pvecback[ba.index_bg_a];
    double Hc = a*pvecback[ba.index_bg_H];

    // =================
    // = Print columns =
    // =================

    // First column is redshift, second is tau, third is scale factor.
    fprintf (stderr, "%+6e %+6e %+6e %+6e ", z, tau, a, log10(a/a_equality));

    // Columns from 5 to th.th_size+3 are the thermodynamical quantities
    for (index_th = 0; index_th < th.th_size; ++index_th) {

      double var = th.thermodynamics_table[index_z*th.th_size+index_th];
      fprintf (stderr, "%+6e ", var);
    }
    
    // Add the ratio between the scattering rate and the Hubble time
    double kappa_dot = th.thermodynamics_table[index_z*th.th_size+th.index_th_dkappa];
    fprintf (stderr, "%+6e ", Hc/kappa_dot);

    // Scattering rate
    fprintf (stderr, "%+6e ", kappa_dot);

    // Conformal hubble  
    fprintf (stderr, "%+6e ", Hc);
    
    fprintf (stderr, "\n");

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
