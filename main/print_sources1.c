/** @file print_sources1.c 
 * Created by Guido W. Pettinari on 17.06.2011
 * Last edited by Guido W. Pettinari on 26.02.2015
 *
 * Print to standard error either the first-order LOS sources from
 * CLASS (contained in pt.sources) or the first-order perturbations
 * needed for SONG (contained in pt.quadsources), according
 * to whether the output requested in the .ini file requires a first
 * or a second-order computation, respectively. 
 *
 * First three columns are conformal time, scale factor, and
 * scale factor normalized to equality a/a_eq.
 *
 * usage:     print_sources1 <ini file> [<pre file>] <index k>"
 *
 */
 
#include "song.h"

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


  /* Fixed indices */
  int index_mode = 0;
  int index_ic = 0;


  // ===================
  // = Parse arguments =
  // ===================

  /* We introduce the n_args variable to differentiate between CLASS arguments and the arguments
    for this function */
  int n_args = 1;
  int index_k;


  if (argc == 2 + n_args) {
    struct stat st;
    stat (argv[1], &st);
    int is_dir = (S_ISDIR (st.st_mode) != 0);
    if (!is_dir) {
      printf ("ERROR: when giving two arguments, the first one ('%s') should be a run directory", argv[1]);
      return _FAILURE_;
    }
    index_k = atoi(argv[2]);
  }
  else if (argc == 3 + n_args) {
    index_k = atoi(argv[3]);
  }
  else {
    printf ("usage:     %s <ini file> <pre file> <index k>\n", argv[0]);
    printf ("           %s <run_directory> <index k>\n", argv[0]);
    return _FAILURE_;
  }

  // ====================================================================
  // =                       Compute perturbations                      =
  // ====================================================================
  
  /* Decrease the argument counter.  The reason is that CLASS should be fed only its default arguments, that is
    the parameter files and, optionally, the run directory */
  argc -= n_args;

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (pt.has_perturbations2 == _TRUE_) {
    if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
      printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
      return _FAILURE_;
    }
    
    /* Compute only the first-order transfer functions, no matter what is specified in the input files */
    pt2.has_early_transfers1_only = _TRUE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }
  
  if (pt.has_perturbations2 == _FALSE_) {
    if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
      printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
      return _FAILURE_;
    }
  }
  else {
    if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
      printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
      return _FAILURE_;
    }
  }

  /* Run some checks */
  if (pt.has_scalars == _FALSE_) {
    printf ("ERROR: only scalar modes supported by this function\n");
    return _FAILURE_;
  }
    
  if ((index_k<0) || (index_k>=pt.k_size[index_mode])) {
    printf ("ERROR: index_k should be between index_k=%d (k=%g) and index_k=%d (k=%g)\n",
    0, pt.k[index_mode][0], pt.k_size[index_mode]-1, pt.k[index_mode][pt.k_size[index_mode]-1]);
    return _FAILURE_;
  }
    

  // ======================================================
  // =              LOS of quadratic sources?             =
  // ======================================================

  /* Should we show the line-of-sight sources or the quadratic sources? (CHOOSE HERE) */
  short tabulate_los_sources = _FALSE_;

  /* If we are in standard 1st-order mode, we can only access the line-of-sight sources */
  if (pt2.has_perturbations2 == _FALSE_)
    tabulate_los_sources = _TRUE_;
  
  
  
  /* Information about the cosmological model */
  double h = ba.h;
  double a_equality = ba.a_eq;  
  fprintf (stderr, "# Cosmological parameters:\n");
  fprintf (stderr, "# tau0 = %g, a_equality = %g, Omega_b = %g, Tcmb = %g, Omega_cdm = %g, omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.conformal_age, ba.a_eq, ba.Omega0_b, ba.T_cmb, ba.Omega0_cdm, ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  fprintf (stderr, "# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);


  /* Sizes associated to the non-running indices in the ***sources table */
  int k_size = pt.k_size[index_mode];
  int tp_size = pt.tp_size[index_mode];

  int qs_size;
  if (pt.has_perturbations2 == _TRUE_)
    qs_size = pt.qs_size[index_mode];

  /* Account for overshooting of 'k' */
  if(index_k > k_size-1) index_k = k_size-1;
  if(index_k < 0) index_k = 0;


  // =================
  // = Print sources =
  // =================    
  /* Running indices */
  int index_tau;
  int index_tp;
  double tau, var, a;
  
  /* Print information on 'k' */
  double k = pt.k[index_mode][index_k];
  printf("# k = %g\n", pt.k[index_mode][index_k]);
  fprintf (stderr, "# k = %g\n", pt.k[index_mode][index_k]);

  /* Vector that will contain background quantities (we are only interested in 'a') */
  double * pvecback = malloc(ba.bg_size_short*sizeof(double));

  /* Number of rows that will be printed */
  int tau_size;
  double tau_ini, tau_end;
  if (tabulate_los_sources == _TRUE_) {
    tau_size = pt.tau_size;
    tau_ini = pt.tau_sampling[0];
    tau_end = pt.tau_sampling[pt.tau_size-1];
  }
  else {
    tau_size = pt.tau_size_quadsources;
    tau_ini = pt.tau_sampling_quadsources[0];
    tau_end = pt.tau_sampling_quadsources[pt.tau_size_quadsources-1];
  };

  /* Some debug info */
  if (tabulate_los_sources == _TRUE_)
    fprintf (stderr, "# Number of source types tp_size=%d\n", tp_size);
  else
    fprintf (stderr, "# Number of source types qs_size=%d\n", qs_size);
  fprintf (stderr, "# Will print %d rows from tau=%g to %g\n", tau_size, tau_ini, tau_end);
  fprintf (stderr, "# Considered mode index_mode=%d\n", index_mode);
  fprintf (stderr, "# Considered initial condition index_ic=%d\n", index_ic);

  /* Running index used to number the columns */
  int index_print=1;


  /* First row contains the labels of the different types */
  fprintf (stderr, "%11s(%03d) ", "tau", index_print++);     // Conformal time
  fprintf (stderr, "%11s(%03d) ", "a", index_print++);       // Scale factor   
  fprintf (stderr, "%11s(%03d) ", "y", index_print++);       // Scale factor normalized to equality
  fprintf (stderr, "%11s(%03d) ", "k_tau", index_print++);   // Scale times Conformal time    

  
  char label[32];

  /* Labels of the line-of-sight sources */
  if (tabulate_los_sources == _TRUE_) {
    for (index_tp = 0; index_tp < tp_size; ++index_tp) {
      sprintf(label, "S_%d", index_tp);
      fprintf (stderr, "%11s(%03d) ", label, index_print++);
    }
  }
  /* Labels for the second-order quadsources */
  else {
    for (index_tp = 0; index_tp < qs_size; ++index_tp)
      fprintf (stderr, "%11s(%03d) ", pt.qs_labels[index_mode][index_tp], index_print++);
  }

  fprintf (stderr, "\n");

  // ==========================================================
  // =                     Loop on Time                       =
  // ==========================================================
  for (index_tau = 0; index_tau < tau_size; ++index_tau) {

    tau = (tabulate_los_sources==_TRUE_ ? pt.tau_sampling[index_tau] : pt.tau_sampling_quadsources[index_tau]);
    
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
    
    fprintf (stderr, "%+16g ", tau);
    fprintf (stderr, "%+16e ", a);
    fprintf (stderr, "%+16e ", log10(a/a_equality));
    fprintf (stderr, "%+16e ", tau*k);


    /* Line of sight sources */
    if (tabulate_los_sources == _TRUE_) {
      for (index_tp = 0; index_tp < tp_size; ++index_tp) {
     
        var = pt.sources[index_mode][index_ic*tp_size + index_tp][index_tau*k_size + index_k];
    
        fprintf (stderr, "%+16e ", var);
      }
    }
    /* Sources for the second-order system */
    else {
      for (index_tp = 0; index_tp < qs_size; ++index_tp) {
     
        var = pt.quadsources[index_mode][index_ic*qs_size + index_tp][index_tau*k_size + index_k];

        fprintf (stderr, "%+16e ", var);
      }
    }
    
    fprintf (stderr, "\n");
    
  } // end of for(index_tau)
  

  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

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
