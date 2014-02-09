/** @file print_kernel_2nd_order.c 
 * Created by Guido W. Pettinari on 08.09.2011
 * Last edited by Guido W. Pettinari on 08.09.2011
 *
 * Print to screen the time and cosk1k2 evolution of a perturbation
 * along with conformal time (1st column), scale factor (second column), and
 * the scale factor normalized to equality log10(a/a_eq) (third column).
 *
 * The columns from 4 to 4+pt2.cosk1k2_size-1 contain the time evolution of the perturbation
 * for cosk1k2 = pt2.cosk1k2(index_cosk1k2).
 *
 *
 * usage:     print_kernel <ini file> [<pre file>] <index_mode> <index_k1> <index_k2> <index_tp>
 *
 */
 
#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions (1st-order) */
  struct perturbs2 pt2;       /* for source functions (2nd-order) */  
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */             
  ErrorMsg errmsg;            /* for error messages */



  // ===================
  // = Parse arguments =
  // ===================
  // Print arguments (debug)
  // int ii=0;
  // for (ii=0; ii<argc; ++ii)
  //   printf("argv[%d] = %s\n", ii, argv[ii]);

  // Fixed indices
  int index_ic = 0;
  
  // Parse arguments
  int index_mode, index_k1, index_k2, index_tp;
  if (argc == 6) {
    index_mode = atoi(argv[2]);
    index_k1 = atoi(argv[3]);
    index_k2 = atoi(argv[4]);
    index_tp = atoi(argv[5]);
  }
  else if (argc == 7) {
    index_mode = atoi(argv[3]);
    index_k1 = atoi(argv[4]);
    index_k2 = atoi(argv[5]);
    index_tp = atoi(argv[6]);
  }
  else {
    printf ("usage:     %s <ini file> [<pre file>] <index_mode> <index_k1> <index_k2> <index_tp>\n", argv[0]);
    return _FAILURE_;
  }


  // ===========================
  // = Calculate perturbations =
  // ===========================
  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&pt2,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  // Set verbosity to zero (we want to be able to send the output to a file and plot it)
  ba.background_verbose = 0;
  th.thermodynamics_verbose = 0;
  pt.perturbations_verbose = 0;
  pt2.perturbations2_verbose = 0;  

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb2_init(&pr,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  // Account for an index_tp out of bounds
  int tp2_size = pt2.tp2_size[index_mode];  
  if((index_tp > tp2_size-1) || (index_tp < 0) ) {
    printf("Wrong type index, returning _FAILED_\n");
    return _FAILURE_;
  }

  // =========================================
  // = Find era of matter-radiation equality =
  // =========================================
  // Using bisection, compute the scale factor at equality.
  // We need it to define the column associated to y=a/a_eq. 
  double a_equality;
  double tol_tau_approx = 1e-6;      // Desired precision in tau
  
  class_call(background_epoch_of_equality(
          &ba,
          tol_tau_approx,
          &a_equality),
        ba.error_message,
        errmsg);
  
  // Print some debug
  double h = ba.h;
  printf("#### cosk1k2 dependence of the variable %s ####\n", pt2.tp2_labels[index_mode][index_tp]);
  printf("# Cosmological parameters:\n");
  printf("# Omega_b = %g, Tcmb = %g, Omega_cdm = %g, omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.Omega0_b, ba.Tcmb, ba.Omega0_cdm, ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  printf("# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);
  printf("# Found a_equality = %g  by bisection\n", a_equality);

  // =================
  // = Peform checks =
  // =================
  // Check that the index_mode is correct
  if ((index_mode >= pt2.md2_size) || (index_mode < 0)) {
    printf ("ERROR: index_mode=%d is out of range. It should be between 0 and %d.\n", index_mode, pt2.md2_size-1);
    return _FAILURE_;
  }

  // Sizes
  int k_size = pt2.k_size[index_mode];
  int cosk1k2_size = pt2.cosk1k2_size[index_mode];

  // Account for overshooting of the wavemodes
  if(index_k1 > k_size-1) index_k1 = k_size-1;
  if(index_k1 < 0) index_k1 = 0;
  if(index_k2 > k_size-1) index_k2 = k_size-1;
  if(index_k2 < 0) index_k2 = 0;


  // =================
  // = Print sources =
  // =================    
  // Running indices
  int index_tau;
  int index_cosk1k2;
  int tau_size = pt2.tau_size;
  double tau, a, cosk1k2;

  // Print information on 'k'
  double k1 = pt2.k[index_mode][index_k1];
  double k2 = pt2.k[index_mode][index_k2];
  printf("# k1 = %g, k2 = %g\n", k1, k2);

  // Vector that will contain background quantities (we are only interested in 'a')
  double * pvecback = malloc(ba.bg_size_short*sizeof(double));

  // Some debug info  
  printf("# index_mode = %d (%s mode)\n", index_mode, pt2.md2_labels[index_mode]);
  printf("# index_ic = %d\n", index_ic);    

  // Print the used cosk1k2 values
  printf("# Will print the following cosk1k2 values:\n");
  printf("# ");
  for (index_cosk1k2 = 0; index_cosk1k2 < cosk1k2_size; ++index_cosk1k2)
    printf("%g ", pt2.cosk1k2[index_mode][index_cosk1k2]);
  printf("\n");    

  // First row contains the labels of the different types
  printf ("%13s ","tau");
  // Second column is the scale factor
  printf ("%13s ","a");
  // Third column is the scale factor normalized to equality, i.e. y as used by Pitrou et al.
  printf ("%13s ","y");
  // The other columns contains the time evolution of a certain cosk1k2
  for (index_cosk1k2 = 0; index_cosk1k2 < cosk1k2_size; ++index_cosk1k2) {
    printf ("%13g ", pt2.cosk1k2[index_mode][index_cosk1k2]);
  }
  printf("\n");

  // =======================
  // = Memorize the kernel =
  // =======================
  // Vector that will contain the kernel (see eq. 45 of Bernardeau, Colombi, et al. 2002
  // and eq. 70 of Pitrou et al. 2008)
  double * kernel = malloc(cosk1k2_size*sizeof(double));
  
  // Vector that will contain delta today (see eq 71 of Pitrou et al. 2008)
  double * delta_today = malloc(cosk1k2_size*sizeof(double));

  // Find delta at first order in k1 and k2
  int index_tau_today = pt.tau_size - 1;
  double delta_1_today = pt.sources[index_mode][index_ic*pt.tp_size[index_mode] + pt.index_tp_delta_cdm]
    [index_tau_today*pt.k_size[index_mode] + index_k1];
  double delta_2_today = pt.sources[index_mode][index_ic*pt.tp_size[index_mode] + pt.index_tp_delta_cdm]
    [index_tau_today*pt.k_size[index_mode] + index_k2];

  // Useful factors
  double k1_sq = k1*k1;
  double k2_sq = k2*k2;
  
  for (index_cosk1k2 = 0; index_cosk1k2 < cosk1k2_size; ++index_cosk1k2) {
    double k1_dot_k2 = k1*k2*pt2.cosk1k2[index_mode][index_cosk1k2];
    kernel[index_cosk1k2] = 5./7. + 1./2. * k1_dot_k2 * (1./k1_sq + 1./k2_sq) + 2./7. * k1_dot_k2*k1_dot_k2/(k1_sq*k2_sq);
    delta_today[index_cosk1k2] = 2. * kernel[index_cosk1k2] * delta_1_today * delta_2_today;
    // delta_today[index_cosk1k2] = delta_1_today * delta_2_today;
  }

  // ================
  // = LOOP ON TIME =
  // ================
  // I am adding a row each iteration of the cycle on tau
  for (index_tau = 0; index_tau < tau_size; ++index_tau) {
  
    tau = pt2.tau_sampling[index_tau];
    
    // Extract value of the scale factor
    int dump;
    background_at_tau(
         &ba,
         tau,
         short_info,
         normal,
         &dump,
         pvecback
         );
    a = pvecback[ba.index_bg_a];     
    
    // First columns are tau, a, a/aeq.  The plus indicates that even plus signs will be printed.
    printf ("%+6e %+6e %+6e ", tau, a, log10(a/a_equality));
    
    // Columns from 4 to cosk1k2_size+3 are the sources
    for (index_cosk1k2 = 0; index_cosk1k2 < cosk1k2_size; ++index_cosk1k2) {
      double var = pt2.sources[index_mode]
              [index_ic*tp2_size + index_tp]
              [index_tau * k_size * k_size * cosk1k2_size +
              + index_k1 * k_size * cosk1k2_size
              + index_k2 * cosk1k2_size
              + index_cosk1k2
              ];
      printf ("%+6e ", var);
    }

    printf ("\n");
  }
  
  printf("%+6e %+6e %+6e ", 0., 0., 0.);
  for (index_cosk1k2 = 0; index_cosk1k2 < cosk1k2_size; ++index_cosk1k2)
    printf("%+6e ", delta_today[index_cosk1k2]);
    
  printf ("\n");

  /****** all calculations done, now free the structures ******/


  if (perturb2_free(&pr2,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pr,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

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
