/** @file print_transfers2.c 
 *
 * Print to screen the second-order transfer functions. The definition
 * of second-order transfer function can be found in eq. 3.69 of
 * http://arxiv.org/abs/1405.2280.
 * 
 * The transfer functions will be printed for a fixed (k1,k2) slice with varying k. 
 *
 *  usage:     ./print_transfers2 <ini file> <pre file> <index_k1> <index_k2> <format string>
 *             ./print_transfers2 <run_directory> <index_k1> <index_k2> <format string>
 *
 * The format string has the form X_L_M, where X=T,E,B is the field (T for temperature,
 * E for E-mode polarisation, B for B-mode polarisation), L=2,3,4... is the angular scale
 * and M=0,1,2... is the azimuthal mode. Note that index_k1 must be larger than index_k2.
 *
 *  example:   ./print_transfers2 params.ini params.pre 1 0 T_2_1
 *
 * Created by Guido W. Pettinari on 03.07.2012
 * Last edited by Guido W. Pettinari on 24.05.2015
 */
 
#include "song.h"

int main (int argc, char **argv) {

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

  // ===============================================================================
  // =                               Parse arguments                               =
  // ===============================================================================

  /* We introduce the n_args variable to differentiate between Song arguments and the arguments
  for this function */
  int n_args = 3;

  /* The arguments are:
  - index_k1: index of k1 wavemode inside ppt2->k
  - index_k2: index of k2 wavemode inside ppt2->k
  - transfer_string: a string of the format X_l_m, where X is the desired type (T for
    temperature, E for E-mode polarization, B for B-mode polarization) and (l,m) is
    the desired angular multipole */

	int index_k1, index_k2, l, m;
  char type;

	/* SONG can accept either one argument (the run directory)... */
  if (argc == 2 + n_args) {
    index_k1 = atoi(argv[2]);
    index_k2 = atoi(argv[3]);
    sscanf (argv[4], "%c_%d_%d", &type, &l, &m);
  }
  /* ... or two arguments (ini and pre files) */
  else if (argc == 3 + n_args) {
    index_k1 = atoi(argv[3]);
    index_k2 = atoi(argv[4]);
    sscanf (argv[5], "%c_%d_%d", &type, &l, &m);
  }
  else {
    printf ("usage:     %s <ini file> <pre file> <index_k1> <index_k2> <format string>\n", argv[0]);
    printf ("           %s <run_directory> <index_k1> <index_k2> <format string>\n", argv[0]);
    printf ("\n");
    printf ("example:   %s params.ini params.pre 1 0 T_2_1\n", argv[0]);
    return _FAILURE_;
  }

  /* Check that index_k2 is correct */
  if (index_k1 < index_k2) {
    printf ("Error: index_k1=%d must be larger than index_k2=%d.\n", index_k1, index_k2);
    return _FAILURE_;
  }


  // ===============================================================================
  // =                           Compute perturbations                             =
  // ===============================================================================

  /* Decrease the argument counter.  The reason is that CLASS should be fed only its
  default arguments, that is the parameter files and, optionally, the run directory */
  argc -= n_args;

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,
    &sp,&nl,&le,&bs,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&tr,&bs,&bs2,
    &tr2,&pm, &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (pt2.has_perturbations2 == _FALSE_) {
    printf ("\n\nSpecify an output requiring a second-order computation (e.g. out=tCl2)\n");
    return _FAILURE_;
  }

  if (pt2.stop_at_perturbations2 == _TRUE_) {
    printf ("\n\nCannot compute late transfer function if pt2.stop_at_perturbations2 is on");
    return _FAILURE_;
  }

  /* Check that we are going to compute the needed transfer type */
  if (type == 'T') {
    if (pt2.has_cmb_temperature==_FALSE_) {
      printf ("ERROR: you cannot ask for CMB temperature if you don't compute it first");
      return _FAILURE_;
    }
  }
  else if (type == 'E') {
    if (pt2.has_cmb_polarization_e==_FALSE_) {
      printf ("ERROR: you cannot ask for E-modes if you don't compute them first");
      return _FAILURE_;
    }
  }
  else if (type == 'B') {
    if (pt2.has_cmb_polarization_b==_FALSE_) {
      printf ("ERROR: you cannot ask for B-modes if you don't compute them first");
      return _FAILURE_;
    }
  }
  
  /* Print the transfer functions without the artificial rescaling factor sin(theta)^m,
  so that they match the standard definition in Eq. 3.69 of http://arxiv.org/abs/1405.2280.
  The transfer functions produced in this way cannot be used to obtain bispectra,
  because the bispectrum module relies on the presence of the rescaling (see eq. 6.26).
  Note that transfer functions loaded from disk are not affected by this option. */
  pt2.rescale_cmb_sources = _FALSE_;

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
  
  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }
  
  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&ba,&th,&tr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (bessel2_init(&pr,&pr2,&pt2,&bs,&bs2) == _FAILURE_) {
    printf("\n\nError in bessel2_init \n =>%s\n",bs2.error_message);
    return _FAILURE_;
  }

  if (transfer2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2) == _FAILURE_) {
    printf("\n\nError in transfer2_init \n=>%s\n",tr2.error_message);
    return _FAILURE_;
  }


  // ===============================================================================
  // =                               What to print?                                =
  // ===============================================================================

  /* Find the index corresponding to the desired transfer-type monopole */
  int index_tt_monopole;
  
  if (type == 'T')
    index_tt_monopole = tr2.index_tt2_T;

  else if (type == 'E')
    index_tt_monopole = tr2.index_tt2_E;

  else if (type == 'B')
    index_tt_monopole = tr2.index_tt2_B;


  /* Find the indices associated to l and m */
  int index_l = tr.index_l[l];
  int index_m = pr2.index_m[m];
  if (index_l<0) {
    printf("\n\nERROR: l=%d is not contained the l-list, which is given by:\n", l);

    for (index_l=0; index_l < (bs.l_size-1); ++index_l)
      printf("%d,", bs.l[index_l]);

    printf("%d\n", bs.l[index_l]);
    return _FAILURE_;
  }
  if (index_m<0) {
    printf("\n\nERROR: m=%d is not contained in the m-list, which is given by:\n", m);

    for (index_m=0; index_m < pr2.m_size; ++index_m)
      printf("%d,", pr2.m[index_m]);

    printf("\n");
    return _FAILURE_;
  }
  if (fabs(m)>l) {
    printf("\n\nERROR: specify l>=m\n");
    return _FAILURE_;
  }
  
  
  /* Account for overshooting of the inputs */
  
  int k1_size = pt2.k_size;
  int k2_size = k1_size - index_k1;
  
  if(index_k1 > k1_size-1) index_k1 = k1_size-1;
  if(index_k1 < 0) index_k1 = 0;
  if(index_k2 > k1_size-1) index_k2 = k1_size-1;  // We use k1_size because index_k2 refers to pt2.k
  if(index_k2 < 0) index_k2 = 0;
  
  double k1 = pt2.k[index_k1];
  double k2 = pt2.k[index_k2];
  
  
  // ===============================================================================
  // =                                  Find k-axis                                =
  // ===============================================================================
  
  /* Extract the first column that will be printed, that is, the k3 axis
  of the transfer function T_lm (k1,k2,k3) */
  double * k;
  class_alloc(k, tr2.k_size_k1k2[index_k1][index_k2]*sizeof(double), tr2.error_message);
  
  int dump;
  if (transfer2_get_k3_list(&pr,&pr2,&pt2,&bs,&bs2,&tr2,index_k1,index_k2,k,&dump)==_FAILURE_) {
    printf("\n\nError running transfer2_get_k3_list\n=>%s\n",tr2.error_message); 
    return _FAILURE_;
  }
  
  
  
  // ====================================================================================
  // =                                Print debug info                                  =
  // ====================================================================================
  
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
  
  /* Find the index_tt corresponding to the (l,m) multipole */
  int index_tt = index_tt_monopole + tr2.lm_array[index_l][index_m];
  // int index_tt = index_tt_monopole
  //   + multipole2offset_indexl_indexm (l, m, tr2.l, tr2.l_size, tr2.m, tr2.m_size);
  
  /* Number of rows to be printed */
  int k_size = tr2.k_size_k1k2[index_k1][index_k2];
  
  /* Print information to stdout */
  printf("Printing (index_k1,index_k2)=(%d,%d) out of %d\n",
    index_k1, index_k2, pt2.k_size);
  printf("corresponding to (l, m, k1, k2, index_tt) = (%d, %d, %g, %g, %d)\n",
    tr2.l[index_l], tr2.m[index_m], k1, k2, index_tt);

  /* Print information to stderr */
  fprintf(stderr, "# k1 = %g, k2 = %g, k3_size = %d, l_size = %d\n",
    k1, k2, k_size, tr2.l_size);
    
  fprintf(stderr, "# Time-sampling of quadsources with %d points from tau=%g to %g\n",
    pt.tau_size_quadsources, pt.tau_sampling_quadsources[0],
    pt.tau_sampling_quadsources[pt.tau_size_quadsources-1]);
  fprintf(stderr, "# Time-sampling of sources with %d points from tau=%g to %g\n",
    pt2.tau_size, pt2.tau_sampling[0], pt2.tau_sampling[pt2.tau_size-1]);
  fprintf(stderr, "# k-sampling of transfer functions with %d points from k=%g to %g\n",
    k_size, k[0], k[k_size-1]);
  if (pr.bispectra_k3_extrapolation != no_k3_extrapolation) {
    int start = tr2.k_physical_start_k1k2[index_k1][index_k2];
    int end = tr2.k_physical_start_k1k2[index_k1][index_k2]
      + tr2.k_physical_size_k1k2[index_k1][index_k2] - 1;
    fprintf(stderr, "# extrapolation for k<%g (row<%d, excluding comments) and k>%g (row>%d)\n",
      k[start], start, k[end], end);
  }
  else {
    fprintf(stderr, "# no extrapolation\n");
  }
  fprintf(stderr, "# l-sampling of transfer functions with %d multipoles from l=%d to %d\n",
    tr2.l_size, tr2.l[0], tr2.l[tr2.l_size - 1]);
  fprintf(stderr, "# Printing transfer function for (index_tt,index_k1,index_k2) = (%d,%d,%d)\n",
    index_tt, index_k1, index_k2);
  fprintf(stderr, "# corresponding to (l,m,k1,k2) = (%d,%d,%g,%g)\n",
    tr2.l[index_l], tr2.m[index_m], k1, k2);
  
  

  // ====================================================================================
  // =                                    Load from disk                                =
  // ====================================================================================
  
  /* Load the transfer functions if we stored them to disk previously, either during a separate
  run or during this run. */
  if ( (pr2.load_transfers_from_disk == _TRUE_) || (pr2.store_transfers_to_disk == _TRUE_) ) {
    if (transfer2_load_transfers_from_disk(&pt2, &tr2, index_tt) == _FAILURE_) {
      printf("\n\nError in transfer2_load_transfers_from_disk \n=>%s\n", tr2.error_message);
      return _FAILURE_;
    }
  }
  else if (pr.load_run == _TRUE_) {
    printf("Error: you are trying to load the transfers from a run that has not them.\n");
    return _FAILURE_;
  }
  

  
  // ==================================================================================
  // =                                  Print transfers                               =
  // ==================================================================================

  int index_print = 1; /* Running index used to number the columns */
  fprintf (stderr, "%25s(%03d) ", "k",index_print++); /* label of x axis */
  fprintf (stderr, "%10s(%03d) ", tr2.tt2_labels[index_tt], index_print++); /* transfers labels */
  fprintf(stderr, "\n");
    
  for (int index_k = 0; index_k < k_size; ++index_k) {
    
    fprintf (stderr, "%+30.17e ", k[index_k]);
          
    double var = tr2.transfer [index_tt]
                              [index_k1]
                              [index_k2]
                              [index_k];
    
    fprintf (stderr, "%+15e ", var);
  
    fprintf (stderr, "\n");
  
  } // end of for(index_k)
  
  
  
  // ==================================================================================
  // =                                 Free structures                                =
  // ==================================================================================
  
  /* Free the memory associated with the line-of-sight transfers for the considered k1 */
  if ( (pr2.load_transfers_from_disk == _TRUE_) || (pr2.store_transfers_to_disk == _TRUE_) ) {
    if (transfer2_free_type_level(&pt2, &tr2, index_tt) == _FAILURE_) {
      printf("\n\nError in transfer2_free_type_level \n=>%s\n",tr2.error_message);
      return _FAILURE_;    
    }
  }
  
  if (transfer2_free(&pr2,&pt2,&tr2) == _FAILURE_) {
    printf("\n\nError in transfer2_free \n=>%s\n",tr2.error_message);
    return _FAILURE_;
  }

  if (bessel2_free(&pr,&pr2,&bs,&bs2) == _FAILURE_)  {
    printf("\n\nError in bessel2_free \n=>%s\n",bs2.error_message);
    return _FAILURE_;
  }
  
  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
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


