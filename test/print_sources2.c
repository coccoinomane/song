/** @file print_sources2.c 
 *
 * Created by Guido W. Pettinari on 10.08.2011
 * Last edited by Guido W. Pettinari on 20.05.2015
 *
 * Print to screen the sources array inside the perturbs2 structure,
 * along with conformal time (first column), scale factor (second column),
 * and the scale factor normalized to equality a/a_eq (third column).
 *
 * usage:     print_sources2 <ini file> <pre file> <variable to print> <index_1> <index_2> <index_3>
 *            print_sources2 <run_directory> <variable to print> <index_1> <index_2> <index_3>
 *
 * <variable to print> can be either 'k2', 'k' or 'tau'.
 * Specify a negative time index to show the sources at the peak of recombination.
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

  /* We introduce the n_args variable to differentiate between CLASS arguments (either 1 or 2) and
  the arguments for this function (n_args) */
  int n_args = 4;
  
  /* Variable that will be printted, read from command line, either "k2", "k", "tau" */
  char variable_to_print[32];

  /* Account for the different inputs accepted by CLASS (either 1 or 2 arguments)  */
  int first_arg;

  if (argc == 2 + n_args) {
    first_arg = 2;
  }
  else if (argc == 3 + n_args) {
    first_arg = 3;
  }
  else {
    printf ("\n");
    printf ("usage:  %s <ini file> <pre file> <variable to print> <index_1> <index_2> <index_3>\n", argv[0]);
    printf ("        %s <run_directory> <variable to print> <index_1> <index_2> <index_3>\n", argv[0]);
    printf ("\n");
    printf ("<variable to print> can be either 'k2', 'k3' or 'tau'\n");
    printf ("The index arguments vary according to what you choose for <variable to print>.\n");
    printf ("Here are the possible combinations:\n");
    printf ("        %s <ini file> <pre file> tau index_k1 index_k2 index_k3\n", argv[0]);
    printf ("        %s <ini file> <pre file> k3 index_k1 index_k2 index_tau\n", argv[0]);
    printf ("        %s <ini file> <pre file> k2 index_k1 index_k3 index_tau\n", argv[0]);
    printf ("Specify a negative time index to show the sources at the peak of recombination.\n");
    return _FAILURE_;
  }


  /* Decide which variable to print according to the first argument */
  strcpy(variable_to_print, argv[first_arg]);
  short PRINT_K2 = _FALSE_;
  short PRINT_K3 = _FALSE_;
  short PRINT_TAU = _FALSE_;
  double * pvecback, a, * pvecthermo, kappa_dot, exp_minus_kappa;

  /* Indices to access pt2.sources[index_type][index_k1][index_k2][index_tau*k3_size+index_k3] */
  int index_k1, index_k2, index_k3, index_tau;
  double k1, k2, k3, cosk1k2, tau;

  if (strcmp(variable_to_print, "k2") == 0) {
    
    PRINT_K2 = _TRUE_;
    index_k1 = atoi(argv[first_arg+1]);
    index_k3 = atoi(argv[first_arg+2]);
    index_tau = atoi(argv[first_arg+3]);
  }
  else if ((strcmp(variable_to_print, "k") == 0)
    || (strcmp(variable_to_print, "k3") == 0)
    || (strcmp(variable_to_print, "cos") == 0)) {
    
    PRINT_K3 = _TRUE_;
    index_k1 = atoi(argv[first_arg+1]);
    index_k2 = atoi(argv[first_arg+2]);
    index_tau = atoi(argv[first_arg+3]);
  }
  else if ((strcmp(variable_to_print, "tau") == 0)
    || (strcmp(variable_to_print, "time") == 0)) {
    
    PRINT_TAU = _TRUE_;
    index_k1 = atoi(argv[first_arg+1]);
    index_k2 = atoi(argv[first_arg+2]);
    index_k3 = atoi(argv[first_arg+3]);
  }
  else {
    printf ("ERROR: The argument 'variable_to_print' must be one between 'k2', 'k', 'tau'\n");
    return _FAILURE_;
  }

  /* Check that index_k2 is correct */
  if (PRINT_K2 == _FALSE_) {
    if (index_k2 > index_k1) {
      printf ("ERROR: index_k1=%d must be larger than index_k2=%d.\n", index_k1, index_k2);
      return _FAILURE_;
    }
  }
  

  // ===============================================================================
  // =                           Compute perturbations                             =
  // ===============================================================================

  /* Decrease the argument counter. The reason is that CLASS should be fed only its
  default arguments, that is the parameter files and, optionally, the run directory */
  argc -= n_args;

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,
    &sp,&nl,&le,&bs,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&tr,&bs,&bs2,&tr2,&pm,
    &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (pt2.has_perturbations2 == _FALSE_) {
    printf ("\n\nSpecify an output requiring a second-order computation (e.g. out=tCl2)\n");
    return _FAILURE_;
  }

  /* Print the source functions without the artificial rescaling factor sin(theta)^m */
  pt2.rescale_quadsources = _FALSE_;

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  
  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }
  
  /* Now that we loaded the background module, we can allocate the background vector */
  if (PRINT_TAU) {
    pvecback = malloc(ba.bg_size*sizeof(double));
    pvecthermo = malloc(th.th_size*sizeof(double));
  }

  /* Print a custom column, and a la mierda */
  // int index_type = pt2.index_tp2_T;
  // int k3_size = pt2.k3_size[index_k1][index_k2];
  // for (index_tau=0; index_tau < pt2.tau_size; ++index_tau)
  //   fprintf (stderr, "%g %g\n", pt2.tau_sampling[index_tau],
  //     pt2.sources[pt2.index_tp2_T][index_k1][index_k2][index_tau*k3_size + index_k3]);
  // class_stop("","");
  
  
  // =====================================================================
  // =                            Perform checks                         =
  // =====================================================================
  
  /* Sizes associated to the non-running indices in the ****sources table */
  int tau_size = pt2.tau_size;
  int k_size = pt2.k_size;  
  int k3_size;
  
  /* Account for overshooting of the indices */
  if(index_k1 > k_size-1) index_k1 = k_size-1;
  if(index_k1 < 0) index_k1 = 0;
  
  if (PRINT_K2 == _FALSE_) {
    if(index_k2 > k_size-1) index_k2 = k_size-1;
    if(index_k2 < 0) index_k2 = 0;
  }

  if (PRINT_K3 == _FALSE_) {
    k3_size = pt2.k3_size[index_k1][index_k2];
    if(index_k3 > k3_size-1) index_k3 = k3_size-1;
    if(index_k3 < 0) index_k3 = 0;
  }

  if (PRINT_TAU == _FALSE_) {
    
    /* Fix sources at recombination time if no time is specified */
    if (index_tau < 0) {
      index_tau = 0;
      while (pt2.tau_sampling[index_tau] < th.tau_rec)
        index_tau++;
    }
    /* Otherwise, check bounds of specified time */
    else {
      if(index_tau > tau_size-1) index_tau = tau_size-1;
      if(index_tau < 0) index_tau = 0;
    }
  }
  
  
  // ==========================================================================
  // =                           Load sources from disk                       =
  // ==========================================================================
  
  /* Load sources from disk if they were previously stored.  This can be true either because we are loading
  them from a precomputed run, or because we stored them in this run. */
  if ( (pr2.load_sources_from_disk == _TRUE_) || (pr2.store_sources_to_disk == _TRUE_) ) {
    if (perturb2_load_sources_from_disk(&pt2, index_k1) == _FAILURE_) {
      printf("\n\nError in perturb2_load_sources_from_disk \n=>%s\n", pt2.error_message);
      return _FAILURE_;
    }       
  }
  else if (pr.load_run == _TRUE_) {
    printf("Error: you are trying to load the sources from a run that has not them.\n");
    return _FAILURE_;
  }
  
  
  // ===================================================================
  // =                           Print debug info                      =
  // ===================================================================
  
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
    
    
  /* Extract values for the quantities that are fixed */
  k1 = pt2.k[index_k1];
  fprintf (stderr, "# Fixing k1 = %g\n", k1);
  printf ("# Fixing k1 = %g\n", k1);
  
  if (PRINT_K2 == _FALSE_) {
    k2 = pt2.k[index_k2];
    fprintf (stderr, "# Fixing k2 = %g\n", k2);
    printf ("# Fixing k2 = %g\n", k2);
  }
  
  // int * k3_indices;
  // double * k3_values;
  
  if (PRINT_K3 == _FALSE_) {
  
    /* The grid in k3 depends on the value of k2. Here we choose the value of k3 which is equal to k2,
    because we know it is in the grid for any k2 value. */
    /* TODO: I think the whole print_k2 business is broken! */
    if (PRINT_K2 == _TRUE_) {
  
      // k3_indices = malloc (k_size*sizeof(int));
      // k3_values = malloc (k_size*sizeof(double));
  
      if (pt2.k3_sampling == smart_k3_sampling) {
  
        k3 = pt2.k[index_k3];
        fprintf (stderr, "# Fixing k3 = %g\n", k3);
        printf ("# Fixing k3 = %g\n", k3);

        // double desired_k3 = pt2.k[index_k2];
        // 
        // for (index_k2 = 0; index_k2 < k_size; ++index_k2) {
        //   
        //   index_k3 = 0;
        // 
        //   if (index_k1 >= index_k2) {
        //     k3 = pt2.k3[index_k1][index_k2][0];
        //     k3_size = pt2.k3_size[index_k1][index_k2];
        //   }
        //   else {
        //     k3 = pt2.k3[index_k2][index_k1][0];
        //     k3_size = pt2.k3_size[index_k2][index_k1];
        //   }
        // 
        //   /* Skip the current k2 value if the k3 sampling is not drawn from pt2.k */
        //   if (k3_size == pr.k3_size_min) {
        //     index_k3 = -1;
        //     continue;
        //   }
        // 
        //   while (k3 < desired_k3) {
        // 
        //     if (index_k1 >= index_k2)
        //       k3 = pt2.k3[index_k1][index_k2][index_k3];
        //     else
        //       k3 = pt2.k3[index_k2][index_k1][index_k3];
        // 
        //     index_k3++;
        //   }
        // 
        //   k3_indices[index_k2] = index_k3;
        //   k3_values[index_k2] = k3;
        // 
        //   /* Some debug */
        //   printf ("k2=%g, k3=%g\n", pt2.k[index_k2], k3);
        // 
        // } // end of for(index_k2)
  
      } // end of if smart sampling
      else {
        printf ("Function not-implemented for k3_sampling != smart_k3_sampling.\n");
        return _SUCCESS_;
      }
    } // end of if print_k2==_TRUE_

    /* When k1 and k2 are both fixed, there is no ambiguity in fixing also the index of k3 */
    else if (PRINT_TAU == _TRUE_) {
      k3 = pt2.k3[index_k1][index_k2][index_k3];
      fprintf (stderr, "# Fixing k3 = %g\n", k3);
      printf ("# Fixing k3 = %g\n", k3);
    }
  } // end of if(print_cos==_FALSE_)
  
  if (PRINT_TAU == _FALSE_) {
    tau = pt2.tau_sampling[index_tau];
    fprintf (stderr, "# Fixing tau = %g\n", tau);
    printf ("# Fixing tau = %g\n", tau);
  }
  
  /* Print information on what we are going to print */
  if (PRINT_K2 == _TRUE_) {
    fprintf (stderr, "# k2 spans %d values from k2=%g to %g\n", index_k1+1, pt2.k[0], pt2.k[index_k1]);
    fprintf (stderr, "# k1 = %g, tau = %g\n", k1, tau);
  }
  else if (PRINT_K3 == _TRUE_) {
    k3_size = pt2.k3_size[index_k1][index_k2];
    double k_min = pt2.k3[index_k1][index_k2][0];
    double k_max = pt2.k3[index_k1][index_k2][k3_size-1];
    double cosk1k2_min = (k_min*k_min - k1*k1 - k2*k2)/(2.*k1*k2);
    double cosk1k2_max = (k_max*k_max - k1*k1 - k2*k2)/(2.*k1*k2);
    fprintf (stderr, "# k spans %d values from k=%g to %g, equivalent to cosk1k2 in (%g,%g)\n",
      k3_size, k_min, k_max, cosk1k2_min, cosk1k2_max);
    fprintf (stderr, "# k1 = %g, k2 = %g, tau = %g\n", k1, k2, tau);
  }
  else if (PRINT_TAU == _TRUE_) {
    k3 = pt2.k3[index_k1][index_k2][index_k3];
    cosk1k2 = (k3*k3 - k1*k1 - k2*k2)/(2.*k1*k2);
    double cosk1k = (k1 + k2*cosk1k2)/k3;
    double cosk2k = (k2 + k1*cosk1k2)/k3;
    fprintf (stderr, "# tau spans %d values from k1=%g to %g\n",
      tau_size, pt2.tau_sampling[0], pt2.tau_sampling[tau_size-1]);
    fprintf (stderr, "# k1 = %g, k2 = %g, k3 = %g, cosk1k = %g, cosk2k = %g, cosk1k2 = %g\n",
      k1, k2, cosk1k2, k3, cosk1k, cosk2k);
  }
  
      
  /* Info about the time sampling */
  fprintf (stderr, "# Time-sampling of quadsources with %d points from tau=%g to %g\n",
    pt.tau_size_quadsources, pt.tau_sampling_quadsources[0], pt.tau_sampling_quadsources[pt.tau_size_quadsources-1]);
  
  if (PRINT_TAU == _FALSE_)
    fprintf (stderr, "# Time-sampling of sources with %d points from tau=%g to %g\n",
      pt2.tau_size, pt2.tau_sampling[0], pt2.tau_sampling[pt2.tau_size-1]);
  
  /* Print number of columns */
  int tp2_size = pt2.tp2_size;
  fprintf (stderr, "# Number of source types: %d\n", tp2_size);
  
  
  
  // ==========================================================
  // =                      Print labels                      =
  // ==========================================================
  
  
  
  /* Running index used to number the columns */
  int index_print=1;
  
  /* First row contains the labels of the different types */
  if (PRINT_K2 == _TRUE_) {
    fprintf (stderr, "%11s(%03d) ", "k2", index_print++);
  }
  else if (PRINT_K3 == _TRUE_) {
    fprintf (stderr, "%11s(%03d) ", "k", index_print++);
    fprintf (stderr, "%11s(%03d) ", "cosk1k2", index_print++);
  } 
  else if (PRINT_TAU == _TRUE_) {
    fprintf (stderr, "%11s(%03d) ", "tau", index_print++);           // Time
    fprintf (stderr, "%11s(%03d) ", "a", index_print++);             // Conformal expansion rate 
    fprintf (stderr, "%11s(%03d) ", "y", index_print++);             // Scale factor normalized to equality
  }
  
  for (int index_type = 0; index_type < tp2_size; ++index_type)
    fprintf (stderr, "%11s(%03d) ", pt2.tp2_labels[index_type], index_print++);
  
  fprintf (stderr, "%11s(%03d) ", "kappa_dot", index_print++);
  fprintf (stderr, "%11s(%03d) ", "exp_m_kappa", index_print++);
  
  fprintf (stderr, "\n");
  
  
  // =======================================================
  // =                     Print rows                      =
  // =======================================================
  
  /* Determine the number of rows to print according to the chosen variable */
  int running_index, n_rows;
  
  if (PRINT_K2 == _TRUE_)
    n_rows = index_k1+1;
  else if (PRINT_K3 == _TRUE_)
    n_rows = k3_size;
  else if (PRINT_TAU == _TRUE_)
    n_rows = tau_size;
  
  
  // ***   Loop on rows   ***
  
  double source;
  int INDEX_K3;
  int WRITE_ZERO;
  
  for (running_index = 0; running_index < n_rows; ++running_index) {
    
    WRITE_ZERO = _FALSE_;
  
    /* The size of the k3 grid is always needed because the tau and k3 levels are mixed inside pt2.sources */
    if (index_k1 >= index_k2)
      k3_size = pt2.k3_size[index_k1][index_k2];
    else
      k3_size = pt2.k3_size[index_k2][index_k1];
  
    /* TODO: I think the whole print_k2 business is broken! */
    if (PRINT_K2 == _TRUE_) {

      index_k2 = running_index;
      k2 = pt2.k[index_k2];
      fprintf (stderr, "%+16e ", k2);
      
      if (pt2.k3_sampling == smart_k3_sampling) { 
      
        int index_k3_min;
      
        if (index_k1 >= index_k2)
          index_k3_min = pt2.index_k3_min[index_k1][index_k2];
        else
          index_k3_min = pt2.index_k3_min[index_k2][index_k1];
                
        /* Skip those k2 points where k3 is not drawn from pt2.k */
        if (index_k3_min < 0)
          WRITE_ZERO = _TRUE_;
      
        INDEX_K3 = index_k3 - index_k3_min;
      
        /* Skip those points that do not satisfy the triangular condition */
        if (INDEX_K3 < 0)
          WRITE_ZERO = _TRUE_;
        
        /* Some debug */
        // fprintf (stderr, "%+16d ", k3_size);
        // fprintf (stderr, "%+16d ", index_k3_min);
        // fprintf (stderr, "%+16d ", INDEX_K3);
      }
      
    }
    else if (PRINT_K3 == _TRUE_) {
     
      INDEX_K3 = running_index;
      k3 = pt2.k3[index_k1][index_k2][INDEX_K3];
      cosk1k2 = (k3*k3 - k1*k1 - k2*k2)/(2.*k1*k2);
      
      fprintf (stderr, "%+16e ", k3);
      fprintf (stderr, "%+16e ", cosk1k2);
    }
    else if (PRINT_TAU == _TRUE_) {
      
      INDEX_K3 = index_k3;
      index_tau = running_index;
      tau = pt2.tau_sampling[index_tau];
  
      /* Extract value of the scale factor */
      int dump;
      background_at_tau (&ba, tau, ba.short_info, ba.inter_normal, &dump, pvecback);
      a = pvecback[ba.index_bg_a];

      /* Extract kappa_dot and exp_minus_kappa */
      thermodynamics_at_z(&ba, &th, 1/a-1, th.inter_normal, &dump, pvecback, pvecthermo);
      kappa_dot = pvecthermo[th.index_th_dkappa];
      exp_minus_kappa = pvecthermo[th.index_th_exp_m_kappa];
      
      /* Print values of time variables */
      fprintf (stderr, "%+16e ", tau);
      fprintf (stderr, "%+16e ", a);
      fprintf (stderr, "%+16e ", log10(a/ba.a_eq));
    }
  
  
    /* Columns from 4 to tp2_size+3 are the sources */
    for (int index_type = 0; index_type < tp2_size; ++index_type) {
      
      if (WRITE_ZERO == _FALSE_) {
        /* We only computed the sources for the case k1 >= k2 because they are symmetric with respect to k1 <-> k2 */
        if (index_k1 >= index_k2)
          source = pt2.sources[index_type][index_k1][index_k2][index_tau*k3_size + INDEX_K3];
        else
          source = pt2.sources[index_type][index_k2][index_k1][index_tau*k3_size + INDEX_K3];
      }
      else {
        source = 0;
      }
        
      fprintf (stderr, "%+16e ", source);
    }
  
    fprintf (stderr, "%+16e ", kappa_dot);
    fprintf (stderr, "%+16e ", exp_minus_kappa);
  
    fprintf (stderr, "\n");
  
  } // end of for(running_index)
  
  
  
  
  // ================================================
  // =              Free structures                 =
  // ================================================
  
  /* Free the memory associated with the line-of-sight sources for the considered k1 */
  if ((pr2.load_sources_from_disk == _TRUE_) || (pr2.store_sources_to_disk == _TRUE_)) {
    if (perturb2_free_k1_level(&pt2, index_k1) == _FAILURE_) {
      printf("\n\nError in perturb2_free_k1_level \n=>%s\n",pt2.error_message);
      return _FAILURE_;    
    }
  }
  
  if (perturb2_free(&pr2,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
    return _FAILURE_;
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
