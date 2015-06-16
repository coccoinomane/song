/** @file print_sources1.c 
 *
 * Print to standard error either the first-order LOS sources from
 * CLASS (contained in pt.sources) or the first-order perturbations
 * needed for SONG (contained in pt.quadsources), according
 * to whether the output requested in the .ini file requires a first
 * or a second-order computation, respectively.
 *
 * For example, if you write 'output=tCl', the function will print the
 * the line-of-sight sources at first order (those that will be later
 * used to compute the first-order C_l). If you instead use
 * 'output=tCl2' or 'output=transfers2', the function will print the
 * sources needed by SONG to solve the second-order differential system
 * (those contained in ppt->quadsource).
 * 
 * The first three columns of the output are conformal time, scale factor,
 * and scale factor normalized to equality (y=a/a_eq), respectively.
 *
 * usage:     print_sources1 <ini file> <pre file> <index k>"
 *            print_sources1 <run_directory> <index k>
 *
 * Created by Guido W. Pettinari on 17.06.2011
 * Last edited by Guido W. Pettinari on 20.05.2015
 */
 
#include "song.h"

int main(int argc, char **argv) {

  struct precision pr;        /* precision parameters (1st-order) */
  struct precision2 pr2;      /* precision parameters (2nd-order) */
  struct background ba;       /* cosmological background */
  struct thermo th;           /* thermodynamics */
  struct perturbs pt;         /* source functions (1st-order) */
  struct perturbs2 pt2;       /* source functions (2nd-order) */  
  struct transfers tr;        /* transfer functions (1st-order) */
  struct bessels bs;          /* bessel functions (1st-order) */
  struct bessels2 bs2;        /* bessel functions (2nd-order) */
  struct transfers2 tr2;      /* transfer functions (2nd-order) */
  struct primordial pm;       /* primordial spectra */
  struct spectra sp;          /* output spectra (1st-order) */
  struct nonlinear nl;        /* non-linear spectra */
  struct lensing le;          /* lensed spectra */
  struct bispectra bi;        /* bispectra */
  struct fisher fi;           /* fisher matrix */
  struct output op;           /* output files */
  ErrorMsg errmsg;            /* error messages */


  // ===================================================================================
  // =                                 Parse arguments                                 =
  // ===================================================================================

  /* We introduce the n_args variable to differentiate between CLASS arguments and the arguments
  for this function */
  int n_args = 1;
  int index_k;

	/* CLASS/SONG can accept either one argument... */
  if (argc == 2 + n_args) {
    struct stat st;
    stat (argv[1], &st);
    int is_dir = (S_ISDIR (st.st_mode) != 0);
    if (!is_dir) {
      printf ("ERROR: when giving two arguments, the first one ('%s') should be a run directory\n", argv[1]);
      return _FAILURE_;
    }
    index_k = atoi(argv[2]);
  }
  /* ... or two arguments */
  else if (argc == 3 + n_args) {
    index_k = atoi(argv[3]);
  }
  else {
    printf ("usage:     %s <ini file> <pre file> <index k>\n", argv[0]);
    printf ("           %s <run_directory> <index k>\n", argv[0]);
    return _FAILURE_;
  }


  // ===================================================================================
  // =                               Compute perturbations                             =
  // ===================================================================================
  
  /* Decrease the argument counter. The reason is that CLASS should be fed only its
  default arguments, that is the parameter files and, optionally, the run directory */
  argc -= n_args;

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,
    &sp,&nl,&le,&bs,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (pt.has_perturbations2 == _TRUE_) {

    if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&tr,&bs,&bs2,&tr2,&pm,
      &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
      printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
      return _FAILURE_;
    }
    
    /* Compute only the first-order early transfer functions, no matter what is specified in
    the input files */
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
  int index_md = pt.index_md_scalars;
    
  if (pt.has_ad == _FALSE_) {
    printf ("ERROR: only adiabatic modes supported by this function\n");
    return _FAILURE_;
  }
  int index_ic = pt.index_ic_ad;
  
  if ((index_k<0) || (index_k>=pt.k_size[index_md])) {
    printf ("ERROR: index_k should be between index_k=%d (k=%g) and index_k=%d (k=%g)\n",
    0, pt.k[index_md][0], pt.k_size[index_md]-1, pt.k[index_md][pt.k_size[index_md]-1]);
    return _FAILURE_;
  }
    

  // ===================================================================================
  // =                             First or second order?                              =
  // ===================================================================================

  /* Should we show the line-of-sight sources or the quadratic sources? */
  short tabulate_los_sources = _FALSE_;

  /* If we are in standard 1st-order mode, we can only access the line-of-sight sources */
  if (pt2.has_perturbations2 == _FALSE_)
    tabulate_los_sources = _TRUE_;
  
  /* Information about the cosmological model */
  double h = ba.h;
  double a_equality = ba.a_eq;  
  fprintf (stderr, "# Cosmological parameters:\n");
  fprintf (stderr, "# tau0 = %g, a_equality = %g, Omega_b = %g, Tcmb = %g, Omega_cdm = %g\n",
    ba.conformal_age, ba.a_eq, ba.Omega0_b, ba.T_cmb, ba.Omega0_cdm);
  fprintf (stderr, "# omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  fprintf (stderr, "# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);

  /* Info about the used gauge at first order */
  fprintf (stderr, "# gauge = ");
  if (pt.gauge == newtonian)
    fprintf (stderr, "Newtonian gauge\n");
  if (pt.gauge == synchronous)
    fprintf (stderr, "synchronous gauge\n");


  /* Sizes associated to the non-running indices in the ***sources table */
  int k_size = pt.k_size[index_md];
  int tp_size = pt.tp_size[index_md];

  int qs_size;
  if (pt.has_perturbations2 == _TRUE_)
    qs_size = pt.qs_size[index_md];

  /* Account for overshooting of 'k' */
  if(index_k > k_size-1) index_k = k_size-1;
  if(index_k < 0) index_k = 0;


  // ===================================================================================
  // =                                   Print sources                                 =
  // ===================================================================================

  /* Print information on 'k' */
  double k = pt.k[index_md][index_k];
  printf("# k = %g (index_k=%d)\n", pt.k[index_md][index_k], index_k);
  fprintf (stderr, "# k = %g (index_k=%d)\n", pt.k[index_md][index_k], index_k);

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
  fprintf (stderr, "# Considered mode: index_md=%d\n", index_md);
  fprintf (stderr, "# Considered initial condition: index_ic=%d\n", index_ic);

  /* Running index used to number the columns */
  int index_print=1;

  /* First row contains the labels of the different types */
  fprintf (stderr, "%11s(%03d) ", "tau", index_print++);     // Conformal time
  fprintf (stderr, "%11s(%03d) ", "a", index_print++);       // Scale factor   
  fprintf (stderr, "%11s(%03d) ", "y", index_print++);       // Scale factor normalized to equality
  fprintf (stderr, "%11s(%03d) ", "k_tau", index_print++);   // Scale times Conformal time    
  
  /* Build labels for the line-of-sight sources */
  char label[32];
  if (tabulate_los_sources == _TRUE_) {
    for (int index_tp = 0; index_tp < tp_size; ++index_tp) {
      sprintf(label, "S_%d", index_tp);
      fprintf (stderr, "%11s(%03d) ", label, index_print++);
    }
  }
  /* Labels for the second-order quadsources */
  else {
    for (int index_tp = 0; index_tp < qs_size; ++index_tp)
      fprintf (stderr, "%11s(%03d) ", pt.qs_labels[index_md][index_tp], index_print++);
  }

  fprintf (stderr, "\n");


  // ===================================================================================
  // =                                   Loop on time                                  =
  // ===================================================================================

  for (int index_tau = 0; index_tau < tau_size; ++index_tau) {

    double tau = (tabulate_los_sources==_TRUE_ ? pt.tau_sampling[index_tau] : pt.tau_sampling_quadsources[index_tau]);
    
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
    double a = pvecback[ba.index_bg_a];      
    
    fprintf (stderr, "%+16g ", tau);
    fprintf (stderr, "%+16e ", a);
    fprintf (stderr, "%+16e ", log10(a/a_equality));
    fprintf (stderr, "%+16e ", tau*k);

    /* Line of sight sources */
    if (tabulate_los_sources == _TRUE_) {
      for (int index_tp = 0; index_tp < tp_size; ++index_tp) {
        double var = pt.sources[index_md][index_ic*tp_size + index_tp][index_tau*k_size + index_k];
        fprintf (stderr, "%+16e ", var);
      }
    }
    /* Sources for the second-order system */
    else {
      for (int index_tp = 0; index_tp < qs_size; ++index_tp) {
        double var = pt.quadsources[index_md][index_ic*qs_size + index_tp][index_tau*k_size + index_k];
        fprintf (stderr, "%+16e ", var);
      }
    }
    
    fprintf (stderr, "\n");
    
  } // end of for(index_tau)
  

  // =====================================================================================
  // =                                    Free memory                                    =
  // =====================================================================================

  if (pt.has_perturbations2 == _TRUE_) {
    if (perturb2_free(&pr2, &pt2) == _FAILURE_) {
      printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
      return _FAILURE_;
    }
  }

  if (perturb_free(&pt) == _FAILURE_) {
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
