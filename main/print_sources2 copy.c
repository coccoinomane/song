/** @file print_sources2.c 
 * Created by Guido W. Pettinari on 10.08.2011
 * Last edited by Guido W. Pettinari on 03.07.2012
 *
 * Print to screen the sources array inside the perturbs2 structure,
 * along with conformal time (1st column), scale factor (second column), and
 * the scale factor normalized to equality a/a_eq (third column).
 *
 * usage:     print_sources2 <ini file> [<pre file>] <index_k1> <index_k2> <index_cosk1k2>
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
  struct transfers tr;        /* for transfer functions (1st-order) */
  struct transfers2 tr2;      /* for transfer functions (2nd-order) */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra (1st-order) */
  struct bispectra bi;        /* for output spectra (2nd-order) */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */



  // ===================
  // = Parse arguments =
  // ===================

  /* We introduce the n_args variable to differentiate between CLASS arguments and the arguments
    for this function */
  int n_args = 3;
  int index_k1, index_k2, index_cosk1k2;

	/* CLASS can accept either one or two arguments */
  if (argc == 2 + n_args) {
    index_k1 = atoi(argv[2]);
    index_k2 = atoi(argv[3]);
    index_cosk1k2 = atoi(argv[4]);
  }
  else if (argc == 3 + n_args) {
    index_k1 = atoi(argv[3]);
    index_k2 = atoi(argv[4]);
    index_cosk1k2 = atoi(argv[5]);
  }
  else {
    printf ("usage:     %s <ini file> <pre file> <index_k1> <index_k2> <index_cosk1k2>\n", argv[0]);
    printf ("           %s <run_directory> <index_k1> <index_k2> <index_cosk1k2>\n", argv[0]);
    return _FAILURE_;
  }

  /* Check that index_k2 is correct */
  if (index_k2 < index_k1) {
    printf ("ERROR: index_k2=%d must be larger than index_k1=%d.\n", index_k2, index_k1);
    return _FAILURE_;
  }


  // ===========================
  // = Calculate perturbations =
  // ===========================
  /* Decrease the argument counter.  The reason is that CLASS should be fed only its default arguments, that is
    the parameter files and, optionally, the run directory */
  argc -= n_args;
  
  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&pt2,&bs,&tr,&tr2,&pm,&sp,&bi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  /* Compute the second-order sources only no matter what is written in params.ini */
  pt2.has_perturbations2 = _TRUE_;
	pt.has_cl_cmb_temperature = pt2.has_cl_cmb_temperature = _TRUE_;


  /* Set verbosity to zero (we want to be able to send the output to a file and plot it) */
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

  

  
  
  /* Information about the cosmological model */
  double h = ba.h;
  double a_equality = ba.a_eq;  
  printf("# Cosmological parameters:\n");
  printf("# tau0 = %g, a_equality = %g, Omega_b = %g, Tcmb = %g, Omega_cdm = %g, omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.conformal_age, ba.a_eq, ba.Omega0_b, ba.T_cmb, ba.Omega0_cdm, ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  printf("# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);

  // =================
  // = Peform checks =
  // =================

  /* Sizes associated to the non-running indices in the ******sources table */
  int k_size = pt2.k_size;  
  int cosk1k2_size = pt2.cosk1k2_size;

  /* Account for overshooting of the wavemodes */
  if(index_k1 > k_size-1) index_k1 = k_size-1;
  if(index_k1 < 0) index_k1 = 0;
  if(index_k2 > k_size-index_k1-1) index_k2 = k_size-index_k1-1;
  if(index_k2 < 0) index_k2 = 0;
  if(index_cosk1k2 > cosk1k2_size-1) index_cosk1k2 = cosk1k2_size-1;
  if(index_cosk1k2 < 0) index_cosk1k2 = 0;


  /* Load sources from disk if they were previously stored.  This can be true either because we are loading
    them from a precomputed run, or because we stored them in this run. */
  if ( (pt2.load_sources_from_disk == _TRUE_) || (pt2.store_sources_to_disk == _TRUE_) ) {
    /* Allocate and fill k1 level of the sources */
    if (perturb2_load_sources_from_disk(&pt2, index_k1) == _FAILURE_) {
      printf("\n\nError in perturb2_load_sources_from_disk \n=>%s\n", pt2.error_message);
      return _FAILURE_;
    }       
  }
  else if (pr.load_run == _TRUE_) {
    printf("Error: you are trying to load the sources from a run that has not them.\n");
    return _FAILURE_;
  }

  // =================
  // = Print sources =
  // =================    
  /* Number of columns to print apart from the tau and 'a' columns */
  int tp2_size = pt2.tp2_size;
  
  /* Running indices */
  int index_tau;
  int index_type;
  int tau_size = pt2.tau_size;
  double tau, a, Hc;
  
  /* Print information on wavemodes */
  double k1 = pt2.k[index_k1];
  double k2 = pt2.k[index_k2];
  double cosk1k2 = pt2.cosk1k2[index_cosk1k2];    
  double k = sqrt(k1*k1 + k2*k2 + 2*k1*k2*cosk1k2);
  double cosk1k = (k1 + k2*cosk1k2)/k;
  double cosk2k = (k2 + k1*cosk1k2)/k;
  printf("# k1 = %g, k2 = %g, cosk1k2 = %g, k = %g, cosk1k = %g, cosk2k = %g\n", k1, k2, cosk1k2, k, cosk1k, cosk2k);
  
  /* Vector that will contain background quantities (we are only interested in 'a') */
  double * pvecback = malloc(ba.bg_size_short*sizeof(double));
  
  /* Some debug info */
  printf("# gauge = ");
  if (pt.gauge == newtonian)
    printf("Newtonian gauge\n");
  if (pt.gauge == synchronous)
    printf("synchronous gauge\n");

  printf("# Time-sampling of quadsources with %d points from tau=%g to %g\n",
    pt.tau_size_quadsources, pt.tau_sampling[0], pt.tau_sampling[pt.tau_size_quadsources-1]);
  printf("# Time-sampling of sources with %d points from tau=%g to %g\n",
    pt2.tau_size, pt2.tau_sampling[0], pt2.tau_sampling[pt2.tau_size-1]);
  printf("# Number of source types: %d\n", tp2_size);    
  
  /* Running index used to number the columns */
  int index_print=1;
  
  /* First row contains the labels of the different types */
  printf ("%11s(%03d) ", "tau", index_print++);          // Time
  printf ("%11s(%03d) ", "a",index_print++);             // Conformal expansion rate 
  printf ("%11s(%03d) ", "y",index_print++);             // Scale factor normalized to equality
  
  for (index_type = 0; index_type < tp2_size; ++index_type)
    printf ("%11s(%03d) ", pt2.tp2_labels[index_type], index_print++);
  
  printf("\n");
  
  // ================
  // = LOOP ON TIME =
  // ================
  for (index_tau = 0; index_tau < pt2.tau_size; ++index_tau) {
    
    tau = pt2.tau_sampling[index_tau];
    
    /* Extract value of the scale factor */
    int dump;
    background_at_tau(
         &ba,
         tau,
         ba.short_info,
         ba.inter_normal,
         &dump,
         pvecback
         );
     a = pvecback[ba.index_bg_a];
     Hc = a*pvecback[ba.index_bg_H];
    
    /* Conformal time.  The plus indicates that even plus signs will be printed. */
    printf ("%+16e ", tau);
  
    /* Scale factor */
    printf ("%+16e ", a);
  
    /* Scale factor normalized to equality, i.e. y as used by Pitrou et al.  */
    printf ("%+16e ", log10(a/a_equality));
      
    /* Columns from 4 to tp2_size+3 are the sources */
    for (index_type = 0; index_type < tp2_size; ++index_type) {
      double var = pt2.sources[index_type][index_k1][index_k2-index_k1][index_cosk1k2][index_tau];
      printf ("%+16e ", var);
    }
  
    printf ("\n");
  }
  
  
  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

  /* Free the memory associated with the line-of-sight sources for the considered k1 */
  if ( (pt2.load_sources_from_disk == _TRUE_) || (pt2.store_sources_to_disk == _TRUE_) )
    if (perturb2_free_k1_level(&pt2, index_k1) == _FAILURE_) {
      printf("\n\nError in perturb2_free_k1_level \n=>%s\n",pt2.error_message);
      return _FAILURE_;    
    }
  
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
