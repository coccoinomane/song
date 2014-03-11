/** @file print_bispectra.c 
 * Created by Guido W. Pettinari on 11.08.2012
 * Last edited by Guido W. Pettinari on 11.08.2012
 *
 * Print to screen bispectra stored in bi.bispectra.
 *
 * usage:     print_bispectra <bispectrum type> <print mode> <options . . .> <ini file> <pre file> 
 *            print_bispectra <bispectrum type> <options . . .> <run_directory> 
 * example:   print_bispectra local_ttt equilateral params.ini params.pre
 *
 */
 
#include "song.h"
#include "input2.h"

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

  // =================================================================================
  // =                              Parse arguments                                  =
  // =================================================================================
  
  /* Indices that address the levels of bi.bispectra[index_bt][index_bp][index_l1][index_l2][index_l3] */
  int index_bt, index_bp, index_l1, index_l2, index_l3;

  /* Names of the selected bispectra (X,Y,Z=T,E...)*/
  char bt_label[128], X_label[128], Y_label[128], Z_label[128];

  /* What to print, read from command line: "triangular", "equilateral", "squeezed_pitrou". These are the meanings
    of the arguments according to the chosen mode:
      "triangular"              :  two arguments, that is the fixed values of l1 and l2
      "equilateral"             :  no argument needed, the output will be a b_L_L_L versus L
      "squeezed_pitrou"         :  1 argument 'ratio'. Tale l2=l3=ratio*l1 and let l1 vary. Requires interpolation.
      "squeezed_large_scale"    :  1 argument 'short_mode'. Take l1=l2=short_mode and let l3 vary.
      "squeezed_small_scale"    :  1 argument 'long_mode'. Take l1=long_mode and let l2=l3 vary.
  */
  char print_mode[32];

  /* Arguments specific for the triangular case */
  int l1_fixed;
  int l2_fixed;
  
  /* Argument specific for the squeezed_pitrou case */
  int short_to_long_ratio;

  /* Arguments specific for the squeezed_small_scale and squeezed_large_scale mode */
  int l_long;
  int l_short;


  /* Minimum number of command line arguments that this function will accept. 
  This is given by the minimum number of arguments accepted by CLASS (one) plus that
  needed by the function itself (two: bispectrum type and configurations to print) */
  int n_min_args = 1 + 2;
  int current_arg = 1;

  /* Usage for the function */
  char usage[32768];
  sprintf (usage, "  *  usage:     %s <bispectra_type> <print mode> <print mode args> <ini file> \n", argv[0]);
  sprintf (usage, "%s  *             %s <bispectra_type> <print mode> <print mode args> <ini file> <pre file> \n", usage, argv[0]);
  sprintf (usage, "%s  *             %s <bispectra_type> <print mode> <print mode args> <run_directory>\n", usage, argv[0]);
  sprintf (usage, "%s  *  The type of bispectra should be in the form 'model_probe', e.g. 'local_ttt', 'intrinsic_tee', 'orthogonal_eee'",
    usage, argv[0]);
  sprintf (usage, "%s  *  \n", usage);
  sprintf (usage, "%s  *  The allowed print modes are:\n", usage);
  sprintf (usage, "%s  *              - triangular <l1> <l2>\n", usage);
  sprintf (usage, "%s  *              - equilateral\n", usage);
  sprintf (usage, "%s  *              - squeezed_pitrou <short_to_long_ratio>\n", usage);
  sprintf (usage, "%s  *              - squeezed_small_scale <long mode>\n", usage);
  sprintf (usage, "%s  *              - squeezed_large_scale <short mode>\n", usage);
  sprintf (usage, "%s  *  Example:\n", usage);
  sprintf (usage, "%s  *              %s local_ttt squeezed_small_scale 6 ini/primordial.ini ini/primordial.pre\n", usage, argv[0]);
  
  /* Parse compulsory arguments */
  if (argc > n_min_args) {
    sscanf(argv[current_arg++], "%[^_]_%c%c%c", bt_label, X_label, Y_label, Z_label);
    strcpy(print_mode, argv[current_arg++]);
  }
  else {
    printf("ERROR: Insufficient number of arguments.\n");
    printf("%s", usage);
    return _FAILURE_;
  }
  
  /* Number of extra arguments needed for the various configurations */
  int n_args_triangular = 2;
  int n_args_equilateral = 0;
  int n_args_squeezed_pitrou = 1;
  int n_args_squeezed_small_scale = 1;
  int n_args_squeezed_large_scale = 1;

  /* Initialize the flags to _FALSE_ */
  short PRINT_TRIANGULAR = _FALSE_;
  short PRINT_EQUILATERAL = _FALSE_;
  short PRINT_SQUEEZED_PITROU = _FALSE_;
  short PRINT_SQUEEZED_SMALL_SCALE = _FALSE_;
  short PRINT_SQUEEZED_LARGE_SCALE = _FALSE_;


  // *** Parse the arguments given according to the chosen 'print_mode'
  
  int remaining_args = argc - current_arg;

  if (strcmp(print_mode, "triangular") == 0) {
    
    if (remaining_args < n_args_triangular) {
      printf("ERROR: The mode 'triangular' needs two arguments.\n");
      printf ("usage:     %s triangular <l1> <l2> <CLASS/SONG arguments ...> \n", argv[0]);
      return _FAILURE_;
    }

    PRINT_TRIANGULAR = _TRUE_;
    if (strlen(argv[current_arg]) > 4) printf ("WARNING: Mixing arguments between SONG and print mode?");
    l1_fixed = atoi(argv[current_arg++]);
    if (strlen(argv[current_arg]) > 4) printf ("WARNING: Mixing arguments between SONG and print mode?");
    l2_fixed = atoi(argv[current_arg++]);
    
  } // end of if(triangular)
  else if (strcmp(print_mode, "equilateral") == 0) {

    if (remaining_args < n_args_equilateral) {
      printf("ERROR: The mode 'equilateral' does not need any argument.\n");
      printf ("usage:     %s ... equilateral <CLASS/SONG arguments ...> \n", argv[0]);
      return _FAILURE_;
    }
  
    PRINT_EQUILATERAL = _TRUE_;
      
  } // end of if(equilateral)
  else if (strcmp(print_mode, "squeezed_pitrou") == 0) {

    if (remaining_args < n_args_squeezed_pitrou) {
      printf("ERROR: The mode 'squeezed_pitrou' needs one argument.\n");
      printf ("usage:     %s ... squeezed_pitrou <short to long mode ratio> <CLASS/SONG arguments ...> \n", argv[0]);
      return _FAILURE_;
    }

    PRINT_SQUEEZED_PITROU = _TRUE_;
    if (strlen(argv[current_arg]) > 4) printf ("WARNING: Mixing arguments between SONG and print mode?");
    short_to_long_ratio = atoi(argv[current_arg++]);

  } // end of if(squeezed_pitrou)
  else if (strcmp(print_mode, "squeezed_small_scale") == 0) {

    if (remaining_args < n_args_squeezed_small_scale) {
      printf("ERROR: The mode 'squeezed_small_scale' needs one argument.\n");
      printf ("usage:     %s ... squeezed_small_scale <long mode> <CLASS/SONG arguments ...> \n", argv[0]);
      return _FAILURE_;
    }

    PRINT_SQUEEZED_SMALL_SCALE = _TRUE_;
    if (strlen(argv[current_arg]) > 4) printf ("WARNING: Mixing arguments between SONG and print mode?");
    l_long = atoi(argv[current_arg++]);

  } // end of if(squeezed_small_scale)
  else if ((strcmp(print_mode, "squeezed_large_scale") == 0) || strcmp(print_mode, "squeezed_large_scales") == 0) {

    if (remaining_args < n_args_squeezed_large_scale) {
      printf("ERROR: The mode 'squeezed_large_scale' needs one argument.\n");
      printf ("usage:     %s ... squeezed_large_scale <short mode> <CLASS/SONG arguments ...> \n", argv[0]);
      return _FAILURE_;
    }

    PRINT_SQUEEZED_LARGE_SCALE = _TRUE_;
    if (strlen(argv[current_arg]) > 4) printf ("WARNING: Mixing arguments between SONG and print mode?");
    l_short = atoi(argv[current_arg++]);

  } // end of if(squeezed_large_scale)

  else {
    printf ("ERROR: Wrong 'print_mode'. Correct syntax follows:\n");
    printf ("%s", usage);
    return _FAILURE_;

  } // end of if(print_mode)


  // =================================================================================
  // =                                  Calculate stuff                              =
  // =================================================================================

  /* Determine which arguments to pass to CLASS/SONG. The '+1' accounts for argv[0] that should be the name of 
  the executable */
  int argc_for_SONG = argc - current_arg + 1;
  char ** argv_for_SONG = malloc (sizeof(char *)*argc_for_SONG);
  int i; for (i=0; i<argc_for_SONG; ++i) argv_for_SONG[i] = malloc(sizeof(char) * 1024);

  /* The first argument is always the name of the binary file being executed */
  strcpy (argv_for_SONG[0], argv[0]);
  
  for (i=1; i<argc_for_SONG; ++i) {
    if (strlen(argv[current_arg]) < 5) {
      printf ("WARNING: the %dth argument, which will be passed to SONG ('%s'), is probably wrong. Correct syntax follows:\n",
      current_arg, argv[current_arg]);
      printf ("%s", usage);
    }    
    strcpy (argv_for_SONG[i], argv[current_arg++]);
  }
  
  /* Some debug */
  // for (i=0; i<argc_for_SONG; ++i)
  //   printf("argv_for_SONG[%d] = %s\n", i, argv_for_SONG[i]);

  if (input_init_from_arguments(argc_for_SONG,argv_for_SONG,&pr,&ba,&th,&pt,&bs,
    &tr,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (input2_init_from_arguments(argc_for_SONG,argv_for_SONG,&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,
    &tr,&tr2,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  /* Flags needed to access the quantities that make up the analytical approximation for the squeezed limit */
  bi.has_intrinsic_squeezed = _TRUE_;
  pt.has_cl_cmb_zeta = _TRUE_;
  sp.compute_cl_derivative = _TRUE_;
  
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  
  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }
  
  if (pt2.has_perturbations2 == _FALSE_)
    if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
      printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
      return _FAILURE_;
    }
  
  if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }
  
  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }
  
  if (bessel2_init(&pr,&pr2,&pt2,&bs,&bs2) == _FAILURE_) {
    printf("\n\nError in bessel2_init \n =>%s\n",bs2.error_message);
    return _FAILURE_;
  }
  
  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  
  if (transfer2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2) == _FAILURE_) {
    printf("\n\nError in transfer2_init \n=>%s\n",tr2.error_message);
    return _FAILURE_;
  }
  
  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
  
  if (spectra_init(&pr,&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
    
  if (bispectra_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra_init \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }
  
  if (bispectra2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra2_init \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }

  /* Uncomment to compute the Fisher matrix too */
  // if (fisher_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi,&fi) == _FAILURE_) {
  //   printf("\n\nError in fisher_init \n=>%s\n",fi.error_message);
  //   return _FAILURE_;
  // }
  
  // =================================================================================
  // =                          Which bispectra to print?                            =
  // =================================================================================
  
  /* Find out the indices corresponding to the requested bispectrum */
  index_bt = 0;
  while (strcmp(bi.bt_labels[index_bt], bt_label)!=0) {
    index_bt++;
    if (index_bt >= bi.bt_size) {
      printf ("ERROR: bispectrum type '%s' does not exist or it wasn't computed. Try one between: ", bt_label);
      for (i=0; i < bi.bt_size; ++i) printf ("%s ", bi.bt_labels[i]); printf ("\n");
      return _FAILURE_;
    }
  }

  /* Find the required fields */
  int X=0, Y=0, Z=0;
  while (strcmp(bi.bf_labels[X], X_label)!=0) {
    X++;
    if (X >= bi.bf_size) {
      printf ("ERROR: field '%s' does not exist or it wasn't computed. Try with: ", X_label);
      for (i=0; i < bi.bf_size; ++i) printf ("%s ", bi.bf_labels[i]); printf ("\n");
      return _FAILURE_;
    }
  }
  while (strcmp(bi.bf_labels[Y], Y_label)!=0) {
    Y++;
    if (Y >= bi.bf_size) {
      printf ("ERROR: field '%s' does not exist or it wasn't computed. Try with: ", Y_label);
      for (i=0; i < bi.bf_size; ++i) printf ("%c ", bi.bf_labels[i]); printf ("\n");
      return _FAILURE_;
    }
  }
  while (strcmp(bi.bf_labels[Z], Z_label)!=0) {
    Z++;
    if (Z >= bi.bf_size) {
      printf ("ERROR: field '%s' does not exist or it wasn't computed. Try with: ", Z_label);
      for (i=0; i < bi.bf_size; ++i) printf ("%c ", bi.bf_labels[i]); printf ("\n");
      return _FAILURE_;
    }
  }

  /* The warning below is useless now, as we just print the configurations where l1>=l2>=l3 when considering
  the analytical approximation */
  //   if (((bi.has_intrinsic_squeezed == _TRUE_) && (index_bt == bi.index_bt_intrinsic_squeezed)))
  //     printf ("WARNING: while printing the squeezed-limit bispectrum, keep in mind that we write as \
  // in eq. 4.1 of Lewis 2012 (http://arxiv.org/abs/1204.5018), which is not symmetrised. Hence, you \
  // should ignore those configurations where the condition l1<=l2<=l3 is not met.\n");


  // =================================================================================
  // =                           Load bispectra from disk                            =
  // =================================================================================
  
  /* Load bispectrum from disk, if it was not already computed. Note that we load the bispectrum from disk
  only if it a secondary bispectrum, since the primary ones are always recomputed as they are very quick
  to obtain. */
  if ((bi.bispectrum_type[index_bt] == non_separable_bispectrum)
   || (bi.bispectrum_type[index_bt] == intrinsic_bispectrum)) {
    if ( (bi.load_bispectra_from_disk == _TRUE_) || (bi.store_bispectra_to_disk == _TRUE_) ) {
      if (bispectra_load_from_disk(&bi, index_bt) == _FAILURE_) {
        printf("\n\nError in bispectra_load_from_disk \n=>%s\n", bi.error_message);
        return _FAILURE_;
      }
    }
  }


  // =================================================================================
  // =                                  Quick access                                 =
  // =================================================================================

  /* Uncomment to quickly access the bispectra without all the complicated machinery that follows. */
  // for (index_l1=0; index_l1 < bi.l_size; ++index_l1) {
  //   
  //   int l1 = bi.l[index_l1];
  //   
  //   for (index_l2=0; index_l2 <= index_l1; ++index_l2) {
  //     
  //     int l2 = bi.l[index_l2];
  // 
  //     /* Skip those configurations that are forbidden by the triangular condition (optimization) */
  //     if (l2 < l1/2)
  //       continue;
  // 
  //     /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
  //     int index_l3_min = bi.index_l_triangular_min[index_l1][index_l2];
  //     int index_l3_max = MIN (index_l2, bi.index_l_triangular_max[index_l1][index_l2]);
  //     
  //     for (index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  
  // 
  //       int l3 = bi.l[index_l3];
  // 
  //       /* Index of the current (l1,l2,l3) configuration */
  //       long int index_l1_l2_l3 = bi.index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  // 
  //       double bispectrum = bi.bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
  // 
  //       /* Print a squeezed configuration */
  //       if ((l3==6) && (l1==l2))
  //         fprintf (stderr, "%5d %17.7g\n", l2, bispectrum);
  //       
  //       /* Print an equilateral configuration */
  //       // if ((l3==l2) && (l1==l2))
  //       //   fprintf (stderr, "%5d %17.7g\n", l1, bispectrum);  
  //     }
  //   }
  // } 
  // 
  // return _SUCCESS_;
  
  /* Uncomment to quickly access the UNSYMMETRISED bispectrum without all the complicated machinery that
  follows. */
  // int index_l_k2, index_l_k3, index_l_k1;
  // 
  // /* EQUILATERAL */
  // for (index_l_k2 = 0; index_l_k2 < bi.l_size; ++index_l_k2) {
  //   int l2 = bi.l[index_l_k2];
  //   for (index_l_k3 = 0; index_l_k3 < bi.l_size; ++index_l_k3) {
  //     int l3 = bi.l[index_l_k3];
  //     int index_l_k1_min = bi.index_l_triangular_min[index_l_k2][index_l_k3];
  //     int index_l_k1_max = bi.index_l_triangular_max[index_l_k2][index_l_k3];
  //     for (index_l_k1=index_l_k1_min; index_l_k1<=index_l_k1_max; ++index_l_k1) {  
  //       int l1 = bi.l[index_l_k1];
  //       double b = bi.unsymmetrised_bispectrum[index_l_k2][index_l_k3][index_l_k1-index_l_k1_min];
  //       
  //       /* The unsymmetrised bispectrum contains <T_l1 E_l2 B_l3> */
  //       if ((l1==l2) && (l2==l3))
  //         fprintf (stderr, "%5d %17.7g\n", l2, b);
  // 
  //     }
  //   }
  // }
  // 
  // /* SLICE */
  // for (index_l_k2 = 0; index_l_k2 < bi.l_size; ++index_l_k2) {
  //   int l2 = bi.l[index_l_k2];
  //   for (index_l_k3 = 0; index_l_k3 < bi.l_size; ++index_l_k3) {
  //     int l3 = bi.l[index_l_k3];
  //     int index_l_k1_min = bi.index_l_triangular_min[index_l_k2][index_l_k3];
  //     int index_l_k1_max = bi.index_l_triangular_max[index_l_k2][index_l_k3];
  //     for (index_l_k1=index_l_k1_min; index_l_k1<=index_l_k1_max; ++index_l_k1) {  
  //       int l1 = bi.l[index_l_k1];
  //       double b = bi.unsymmetrised_bispectrum[index_l_k2][index_l_k3][index_l_k1-index_l_k1_min];
  //       
  //       /* The unsymmetrised bispectrum contains <T_l1 E_l2 B_l3> */
  //       if ((l1==105) && (l2==l3))
  //         fprintf (stderr, "%5d %17.7g\n", l2, b);
  // 
  //     }
  //   }
  // }
  // 
  // return _SUCCESS_;
  
  
  // =================================================================================
  // =                                  Perform checks                               =
  // =================================================================================
  
  /* Maximum and minimum scales considered */
  int l_min = bi.l[0];
  int l_max = bi.l[bi.l_size-1];
    
  /* Check that the l's specified by the user correspond to multipoles for which we computed the bispectrum */
  int index_l;
  
  int index_l1_fixed;
  int index_l2_fixed;
  
  short l1_in_list = _FALSE_;
  short l2_in_list = _FALSE_;
  
  int index_l_short;
  int index_l_long;
  
  if (PRINT_TRIANGULAR == _TRUE_) {
  
    index_l1_fixed = bs.index_l[l1_fixed];
    index_l2_fixed = bs.index_l[l2_fixed];
  
    if ((index_l1_fixed<0) || (index_l2_fixed<0)) {
  
      printf("\n\nERROR: l1 and l2 should be both chosen from the following list:\n");
      for (index_l=0; index_l < (bi.l_size-1); ++index_l)
        printf("%d,", bi.l[index_l]);
      printf("%d\n", bi.l[index_l]);

      return _FAILURE_;
    }
  }
  else if ( PRINT_SQUEEZED_SMALL_SCALE == _TRUE_ ) {
  
    index_l_long = bs.index_l[l_long];
  
    if (index_l_long < 0) {
  
      printf("\n\nERROR: l_long should be chosen from the following list:\n");  
      for (index_l=0; index_l < (bi.l_size-1); ++index_l)
        printf("%d,", bi.l[index_l]);
      printf("%d\n", bi.l[index_l]);
  
      return _FAILURE_;
    }
  }
  else if ( PRINT_SQUEEZED_LARGE_SCALE == _TRUE_ ) {
  
    index_l_short = bs.index_l[l_short];
  
    if (index_l_short < 0) {
  
      printf("\n\nERROR: l_short should be chosen from the following list:\n");  
      for (index_l=0; index_l < (bi.l_size-1); ++index_l)
        printf("%d,", bi.l[index_l]);
      printf("%d\n", bi.l[index_l]);
  
      return _FAILURE_;
    }
  }
  
  
  // =================================================================================
  // =                                 Print debug info                              =
  // =================================================================================
  
  /* Print Information about the cosmological model */
  double h = ba.h;
  
  fprintf(stderr, "# Cosmological parameters:\n");
  fprintf(stderr, "# tau0 = %g, a_equality = %g, Omega_b = %g, Tcmb = %g, Omega_cdm = %g, omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.conformal_age, ba.a_eq, ba.Omega0_b, ba.T_cmb, ba.Omega0_cdm, ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  fprintf(stderr, "# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);
  
  
  /* Info about the used gauge */
  fprintf(stderr, "# gauge = ");
  if (pt.gauge == newtonian)
    fprintf(stderr, "Newtonian gauge\n");
  if (pt.gauge == synchronous)
    fprintf(stderr, "synchronous gauge\n");
  
  /* Which bispectrum are we printing? */
  fprintf(stderr, "# Showing the '%s_%s%s%s' bispectrum (index_bt=%d)\n", bt_label, X_label, Y_label, Z_label, index_bt);
  
  /* Print information on what we are going to print */
  if (PRINT_TRIANGULAR == _TRUE_) {
    
    int l3_size = bi.l_triangular_size[index_l1_fixed][index_l2_fixed];
    int index_l3_min = bi.index_l_triangular_min[index_l1_fixed][index_l2_fixed]; 
    int l3_min = bi.l[index_l3_min];
  
    fprintf (stderr, "# Tabulating TRIANGULAR configurations of b(l1,l2,l3) with l1=%d, l2=%d, l3=L free to vary\n", l1_fixed, l2_fixed);
    fprintf (stderr, "# L spans %d values from %d to %d\n", l3_size, l3_min, bi.l[index_l3_min + l3_size - 1]);
  }
  else if (PRINT_EQUILATERAL == _TRUE_) {
  
    fprintf (stderr, "# Tabulating EQUILATERAL configurations of b(l1,l2,l3) with l1=l2=l3=L free to vary\n");
    fprintf (stderr, "# L spans %d values from %d to %d\n", bi.l_size, bi.l[0], bi.l[bi.l_size-1]);
  } 
  else if (PRINT_SQUEEZED_PITROU == _TRUE_) {
  
    fprintf (stderr, "# Tabulating SQUEEZED_PITROU configurations of b(l1,l2,l3) with %d*l1=l2=l3, with l1=L free to vary\n", short_to_long_ratio);
  }
  else if (PRINT_SQUEEZED_SMALL_SCALE == _TRUE_) {
  
    fprintf (stderr, "# Tabulating SQUEEZED_SMALL_SCALE configurations of b(l1,l2,l3) with l1=%d, l2=l3=L free to vary\n", l_long);
    fprintf (stderr, "# L spans %d values from %d to %d\n", bi.l_size, bi.l[0], bi.l[bi.l_size-1]);
  }
  else if (PRINT_SQUEEZED_LARGE_SCALE == _TRUE_) {
  
    fprintf (stderr, "# Tabulating SQUEEZED_LARGE_SCALE configurations of b(l1,l2,l3) with l1=l2=%d, l3=L free to vary\n", l_short);
    fprintf (stderr, "# L spans %d values from %d to %d\n", bi.l_size, bi.l[0], bi.l[bi.l_size-1]);
  }

  
  // =================================================================================
  // =                                    Print labels                               =
  // =================================================================================
    
  /* Running index used to number the columns */
  int index_print = 1;
  
  /* Multipole L */
  fprintf (stderr, "%15s(%03d) ", "L", index_print++);

  /* Reduced bispectrum of the brightness delta (momentum integrated distribution function) */
  fprintf (stderr, "%20s(%03d) ", "brightness_D", index_print++);

  /* Reduced bispectrum of the brightness temperature */
  fprintf (stderr, "%20s(%03d) ", "brightness_T", index_print++);

  /* Reduced bispectrum of the bolometric temperature */
  fprintf (stderr, "%20s(%03d) ", "bolometric_T", index_print++);

  /* Correction needed to convert between the brightness temperature bispectrum to the bolometric temperature one
  (bolometric = brightness + correction) */
  fprintf (stderr, "%20s(%03d) ", "T_correction", index_print++);

  /* Bispectra with the factor used by Komatsu (see for example Fig. 4.5 of Komatsu's thesis, or Pitrou et al. 2010 in Fig. 3) */
  fprintf (stderr, "%20s(%03d) ", "brightness_T_K", index_print++);
  fprintf (stderr, "%20s(%03d) ", "bolometric_T_K", index_print++);
  

  /* Analytical approximation in Lewis 2012 (L2012).
  We normalize all the quantities in these papers by -1/12 C_long * C_short */
  if ( (PRINT_SQUEEZED_SMALL_SCALE == _TRUE_) || (PRINT_SQUEEZED_LARGE_SCALE == _TRUE_) ) {

    /* Normalized brightness temperature reduced bispectrum */
    fprintf (stderr, "%25s(%03d) ", "brightness_T_norm", index_print++);

    /* Normalized bolometric temperature reduced bispectrum */
    fprintf (stderr, "%25s(%03d) ", "bolometric_T_norm", index_print++);

    /* Sum of the eq. 2.5 and 2.6 in Lewis 2012 */
    fprintf (stderr, "%25s(%03d) ", "bolometric_T_lewis", index_print++);

    /* Approximation in Lewis 2012 (eq. 2.5) */
    fprintf (stderr, "%25s(%03d) ", "bolometric_T_ricci_lewis", index_print++);

    /* Approximation in Lewis 2012 (eq. 2.6) */
    fprintf (stderr, "%25s(%03d) ", "bolometric_T_redsh_lewis", index_print++);

    /* Approximation in CPV2011 (eq 4.3) */
    fprintf (stderr, "%25s(%03d) ", "bolometric_T_lensing_cpv", index_print++);

    /* Normalized temperature correction */
    fprintf (stderr, "%25s(%03d) ", "T_correction_norm", index_print++);

    /* Normalization used to compare the analytical solution with the numerical one */
    fprintf (stderr, "%25s(%03d) ", "normalization", index_print++);
      
    /* C_l for zeta and temperature. Should be equal to the magenta line in Fig. 3 of L2012 */
    fprintf (stderr, "%25s(%03d) ", "cl_Xz_short", index_print++);
    fprintf (stderr, "%25s(%03d) ", "cl_Xz_long", index_print++);
    fprintf (stderr, "%25s(%03d) ", "cl_XX_short", index_print++);
    fprintf (stderr, "%25s(%03d) ", "cl_XX_long", index_print++);
    
    /* In the squeezed limit, the local bispectrum is equal to -5*(Cl1*Cl2+Cl1*Cl3+Cl2*Cl3) */
    // fprintf (stderr, "%25s(%03d) ", "-5*(Cl1*Cl2+Cl1*Cl3+Cl2*Cl3)", index_print++);  

  } // end of(PRINT_SQUEEZED_SMALL_OR_LARGE_SCALE)
  
  fprintf(stderr, "\n");
  
  
  // ===================================================================================
  // =                                    Print columns                                =
  // ===================================================================================
  
  
  /* We now loop over the (l1,l2,l3) configurations that satisfy the triangular condition. Since we
  computed the bispectra only for l1>=l2>=l3, we obtain the values outside that range using
  the symmetry of the bispectrum wrt to permutations of (l1,l2,l3) */
  for (index_l1=0; index_l1 < bi.l_size; ++index_l1) {
    
    int l1 = bi.l[index_l1];
    
    for (index_l2=0; index_l2 < bi.l_size; ++index_l2) {
      
      int l2 = bi.l[index_l2];
  
      for (index_l3=bi.index_l_triangular_min[index_l1][index_l2]; index_l3<=bi.index_l_triangular_max[index_l1][index_l2]; ++index_l3) {  
  
        int l3 = bi.l[index_l3];
  
        /* We have computed the bispectrum only for those configurations where l1>=l2>=l3. However, we
        print it for all (l1,l2,l3) configurations by permuting the field indices XYZ, as all physical
        bispectra of the form <X_l1 Y_l2 Z_l3> are symmetric with respec to the exchanges
        (l1,X) <-> (l2,Y), (l1,X) <-> (l3,Z) and (l2,Y) <-> (l3,Z). This cannot be done
        for the bispectra that are inherently asymmetric, like it is the case for the analytical approximation
        for the squeezed configurations. In that case, we only print the configurations where l1>=l2>=l3, as
        the other ones would not make sense. (In practice, we enforce l3>=l2>=l1 instead of l1>=l2>=l3
        because we have assigned l1 to be l_long, that is, the smallest of the three.) */
        if ((((bi.has_intrinsic_squeezed == _TRUE_) && (index_bt == bi.index_bt_intrinsic_squeezed)))
           ||((bi.has_local_squeezed == _TRUE_) && (index_bt == bi.index_bt_local_squeezed)))
          if ((l1>l2) || (l1>l3) || (l2>l3))
            continue;
  
        /* Skip unphysical configurations. For intensity, they are those where l1+l2+l3 is odd.
        However, when only the scalar bispectrum was computed, there is no need to do so because
        in that case we have extracted the Gaunt structure analytically. */
        /* TODO: update for B-modes */
        // if (pr2.m_max_2nd_order > 0)
        //   if ((bi.l[index_l1]+bi.l[index_l2]+bi.l[index_l3])%2!=0)
        //     continue;

        // ---------------------------------------------------------------------------
        // -                        Sort the l1,l2,l3 multipoles                     -
        // ---------------------------------------------------------------------------
  
        /* Sort the 3 indices */
        int index_largest, index_mid, index_smallest;
  
        if (l1 >= l2) {
          if (l2 >= l3) {
            index_largest = index_l1;
            index_mid = index_l2;
            index_smallest = index_l3;
          }
          else /* if l3 > l2 */ {
            if (l1 >= l3) {
              index_largest = index_l1;
              index_mid = index_l3;
              index_smallest = index_l2;
            }
            else /* if l3 > l1 */ {
              index_largest = index_l3;
              index_mid = index_l1;
              index_smallest = index_l2;
            }
          }
        }
        else /* if l2 > l1 */ {
          if (l1 >= l3) {
            index_largest = index_l2;
            index_mid = index_l1;
            index_smallest = index_l3;
          }
          else /* if l3 > l1 */ {
            if (l2 >= l3) {
              index_largest = index_l2;
              index_mid = index_l3;
              index_smallest = index_l1;
            }
            else /* if l3 > l2 */ {
              index_largest = index_l3;
              index_mid = index_l2;
              index_smallest = index_l1;
            }
          }
        }

        /* Double check that the sorting worked */
        int l_largest = bi.l[index_largest];
        int l_mid = bi.l[index_mid];
        int l_smallest = bi.l[index_smallest];

        if ((l_smallest>l_mid) || (l_smallest>l_largest) || (l_mid>l_largest)) {
         printf ("ERROR in sorting l1,l2,l3\n");
         return _FAILURE_;
        }
  
        // -------------------------------------------------------------------------------------------------------
        // -                                   Extract the computed bispectra                                    -
        // -------------------------------------------------------------------------------------------------------

        /* Index of the current (l1,l2,l3) configuration */
        int index_max_for_smallest = MIN (index_mid, bi.index_l_triangular_max[index_largest][index_mid]);
        long int index_l1_l2_l3 = bi.index_l1_l2_l3[index_largest][index_largest-index_mid][index_max_for_smallest-index_smallest];
        
        /* Value of the bolometric temperature bispectrum in (l1,l2,l3) */
        double bolometric_T;

        /* Since we computed the bispectra only for l1>=l2>=l3, we obtain the values outside that range using
        the symmetry of the bispectrum with respect to permutations of (l1,l2,l3). This approach won't work
        for the squeezed-limit approximations of the i_squeezed and l_squeezed bispectra (see comment above at
        the beginning of l3 loop), but in that case the values of l1,l2 and l3 in this loop are constrained
        to satisfy l1<=l2<=l3. */
             if ((l1>=l2) && (l2>=l3))
          bolometric_T = bi.bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
        else if ((l1>=l3) && (l3>=l2))
          bolometric_T = bi.bispectra[index_bt][X][Z][Y][index_l1_l2_l3];
        else if ((l2>=l1) && (l1>=l3))
          bolometric_T = bi.bispectra[index_bt][Y][X][Z][index_l1_l2_l3];
        else if ((l2>=l3) && (l3>=l1))
          bolometric_T = bi.bispectra[index_bt][Y][Z][X][index_l1_l2_l3];
        else if ((l3>=l1) && (l1>=l2))
          bolometric_T = bi.bispectra[index_bt][Z][X][Y][index_l1_l2_l3];
        else if ((l3>=l2) && (l2>=l1))
          bolometric_T = bi.bispectra[index_bt][Z][Y][X][index_l1_l2_l3];
        
        /* ADD BY HAND THE QUADRATIC CORRECTIONS (temp + redshift) */
        // bolometric_T += (4-3)/8. * bi.bispectra[bi.index_bt_quadratic][X][Y][Z][index_l1_l2_l3];
        
        /* Corrective factor to convert the BOLOMETRIC TEMPERATURE bispectrum to the BRIGHTNESS TEMPERATURE bispectrum.
        These are related by: BRIGHTNESS = BOLOMETRIC + 3*(Cl1*Cl2 + Cl2*Cl3 + Cl1*Cl3). See, e.g., Pitrou et al. 2010
        or eq. 2.14 Nitta et al. 2009.

        IMPORTANT: subtracting the temperature correction produces numerical noise if the bispectrum
        is small, as it is the case when m!=0. This happens because in the bispectrum module we summed
        the (big) temperature correction to the (small) bispectrum and here we are subtracting it.
        However, this is only a visualisation issue relevant to this function that does not affect Song. */
        /* TODO: update for polarization. This means loading a full bispectrum rather than computing it here
        on the fly */
        double Cl1 = sp.cl[pt.index_md_scalars][index_l1 * sp.ct_size + bi.index_ct_of_bf[X][X]];
        double Cl2 = sp.cl[pt.index_md_scalars][index_l2 * sp.ct_size + bi.index_ct_of_bf[X][X]];
        double Cl3 = sp.cl[pt.index_md_scalars][index_l3 * sp.ct_size + bi.index_ct_of_bf[X][X]];
        double temperature_correction = 3*(Cl1*Cl2 + Cl1*Cl3 + Cl2*Cl3);
  
        /* Value of the brightness temperature bispectrum */
        double brightness_T = bolometric_T + temperature_correction;

        /* Only for the intrinsic bispectra there is distinction between bolometric and brightness temperature.
        We also kill the quadratic contributions when pt2.has_quadratic_sources=_FALSE_  */
        if (((strstr(bt_label,"intrinsic")==NULL) && (strstr(bt_label,"squeezed")==NULL))
          ||((strstr(bt_label,"intrinsic")!=NULL) && (pt2.has_quadratic_sources==_FALSE_)))
          brightness_T = bolometric_T;
          
        /* Value of the bispectrum of the brightness (momentum-integrated distribution function) */
        double brightness_temperature_factor = (2*l1+1.)*(2*l2+1.)*(2*l3+1.) * pow(2.*_PI_,3);
        double brightness_D = brightness_T * brightness_temperature_factor;

  
        // ----------------------------------------------------------------------------------------------
        // -                                 Find out which l's to print                                -
        // ----------------------------------------------------------------------------------------------
  
        /* First column will be 'l' */
        int l;
          
        // ****  TRIANGULAR CASE  ****
        if ((PRINT_TRIANGULAR == _TRUE_) && (l1==l1_fixed) && (l2==l2_fixed)) {
          
          l = l3;
          
        } // end of if(triangular)
          
        // ****  EQUILATERAL CASE  ****        
        else if ((PRINT_EQUILATERAL == _TRUE_) && (l1==l2) && (l2==l3)) {
          
          l = l3;
          
        } // end of if(equilateral)
            
        // ****  SQUEEZED_PITROU CASE  ****        
        else if ((PRINT_SQUEEZED_PITROU == _TRUE_) && (l2==l3) && ((short_to_long_ratio*l1)==l2)) {
          
          l = l1;
          
        } // end of if(squeezed_pitrou)
          
        // ****  SQUEEZED_SMALL_SCALE CASE  ****        

        /* We choose the long wavemode to be l1, which is associated to the field X. The short wavemodes
        will be l2 and l3, associated respectively to Y and Z */
        
        else if ((PRINT_SQUEEZED_SMALL_SCALE == _TRUE_) && (l1==l_long) && (l2==l3)) {
            
          l = l3;
          l_short = l3;
          index_l_short = index_l3;
          
        } // end of if(squeezed_small_scale)
          
        // ****  SQUEEZED_LARGE_SCALE CASE  ****        
        else if ((PRINT_SQUEEZED_LARGE_SCALE == _TRUE_) && (l1==l_short) && (l2==l_short)) {
          
          l = l3;
          l_long = l3;
          index_l_long = index_l3;
          
        } // end of if(squeezed_large_scale)
          
        /* If nothing has to be printed, then continue */
        else {
          continue;
        }
        
        /* Coefficient adopted by Komatsu and by Pitrou et et al. 2010 in Fig. 3 */
        double komatsu_coefficient = 1e16 * l*l * (l+1)*(l+1) / ((2*_PI_)*(2*_PI_));
          
        fprintf(stderr, "%20d %25.7g %25.7g %25.7g %25.7g %25.7g %25.7g ",
          l,
          brightness_D, 
          brightness_T, 
          bolometric_T, 
          temperature_correction, 
          komatsu_coefficient*brightness_T, 
          komatsu_coefficient*bolometric_T
        );
                    
        // ========================================================================================================
        // =                                      Squeezed analytical approximation                               =
        // ========================================================================================================
          
        if ((PRINT_SQUEEZED_SMALL_SCALE==_TRUE_) || (PRINT_SQUEEZED_LARGE_SCALE==_TRUE_)) {
          
          int mode = pt.index_md_scalars;
                    
          /* C_l's for the short and long modes */
          double cl_YZ_short    = bi.cls[bi.index_ct_of_bf[Y][Z]][l_short-2];
          double cl_XX_long     = bi.cls[bi.index_ct_of_bf[X][X]][l_long-2];
          double cl_XX_short    = bi.cls[bi.index_ct_of_bf[X][X]][l_short-2];
          double cl_XY_long     = bi.cls[bi.index_ct_of_bf[X][Y]][l_long-2];

          /* Determine which field has to be correlated with the comoving curvature perturbation zeta.
          The resulting C_l^zeta always gets the largest scale multipole */
          int index_ct_Xz, index_ct_Yz, index_ct_Zz;

          if ((bi.has_bispectra_t == _TRUE_) && (X == bi.index_bf_t)) {
            index_ct_Xz = sp.index_ct_tz;
          }
          else if ((bi.has_bispectra_e == _TRUE_) && (X == bi.index_bf_e)) {
            index_ct_Xz = sp.index_ct_ez;
          }
          else {
            index_ct_Xz = 0;
            printf ("WARNING: C_l's for <X * zeta> where X=%s not found. Do not trust approximations.\n", bi.bf_labels[X]);
          }
          double cl_Xz_long = bi.cls[index_ct_Xz][l_long-2];
          double cl_Xz_short = bi.cls[index_ct_Xz][l_short-2];

          /* We include an l-dependent normalization factor so that the comparison between numerical and analytical 
          result won't depend on the adopted conventions for the Cl's or the primordial power spectrum.
          A good candidate would have been the approximation for the local squeezed limit, given by
          - 5 / (12 * cl_Xz_long * cl_YZ_short). However, this normalization function crosses the zero, and would
          result in infinities. Therefore we choose simply - 1 / (12 * cl_TT_long * cl_TT_short) for the TTT and mixed
          bispectra, and - 1 / (12 * cl_EE_long * cl_EE_short) for the EEE bispectrum. */
          double normalization;

          if (bi.has_bispectra_t == _TRUE_)
            normalization = - 1 / (12 * bi.cls[sp.index_ct_tt][l_long-2] * bi.cls[sp.index_ct_tt][l_short-2]);
          else if (bi.bf_size == 1)
            normalization = - 1 / (12 * bi.cls[bi.index_ct_of_bf[X][X]][l_long-2] * bi.cls[bi.index_ct_of_bf[X][X]][l_short-2]);
          else {
            printf ("ERROR: case of 2 non-temperature fields not implemented yet.\n");
          }
        
          /* Normalized version of the brightness temperature bispectrum */
          double brightness_T_normalized = normalization * brightness_T;
          
          /* Normalized version of the bolometric temperature bispectrum */
          double bolometric_T_normalized = normalization * bolometric_T;
          
          /* Normalized temperature correction. Note that in this squeezed_small_scale limit, it is also given by
          -3 * cl_long*cl_short * (2 + cl_short/cl_long) * normalization */
          double temperature_correction_normalized = normalization * temperature_correction;
          
          // -----------------------------------------------------------------
          // -                         Lewis 2012                            -
          // -----------------------------------------------------------------

          /* Here we compute the approximations in eq. 4.1 and 4.2 of Lewis 2012. This is the general
          formula that includes polarisation. With respect to Lewis' formula, i->X, j->Y, k->Z and
          greek zeta -> z. */
          double dcl_YZ_short = bi.d_lsq_cls[bi.index_ct_of_bf[Y][Z]][l_short-2];
          
          /* Ricci focussing in Lewis 2012 (eq. 4.1) */
          double bolometric_T_lewis_ricci = normalization * (- cl_Xz_long * dcl_YZ_short/l_short);
          
          /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
          double bolometric_T_lewis_redshift = 0;
          
          if (bi.has_bispectra_t == _TRUE_) {
            double cl_Xt_long = 0; double cl_Yt_short = 0; double cl_Zt_short = 0;
            cl_Xt_long = bi.cls[bi.index_ct_of_bf[X][bi.index_bf_t]][l_long-2];
            if (Z == bi.index_bf_t) cl_Yt_short = bi.cls[bi.index_ct_of_bf[Y][bi.index_bf_t]][l_short-2];
            if (Y == bi.index_bf_t) cl_Zt_short = bi.cls[bi.index_ct_of_bf[Z][bi.index_bf_t]][l_short-2];
            bolometric_T_lewis_redshift = normalization * (cl_Xt_long * (cl_Yt_short + cl_Zt_short));
          }
          
          /* Sum of Ricci focussing and redshift modulation */
          double bolometric_T_lewis = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

          // -------------------------------------------------------
          // -         Creminelli, Pitrou & Vernizzi 2011          -
          // -------------------------------------------------------
          
          /* Lensing contribution, that is eq. 4.3 in CPV2011 */
          double cos_theta = (l*l - l_long*l_long - l_short*l_short)/(2.*l_long*l_short);
          double cos_2_theta = cos_theta*cos_theta - (1-cos_theta*cos_theta);
          double bolometric_T_cpv_lensing = normalization * 6*cl_XX_long*(cl_XX_short*2*cos_2_theta - (1 + cos_2_theta)*dcl_YZ_short/l_short);
          
          // -------------------------------------------------------
          // -                 Print extra columns                 -
          // -------------------------------------------------------
          
          fprintf(stderr, "%30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g %30.7g",
            brightness_T_normalized,
            bolometric_T_normalized,
            bolometric_T_lewis,
            bolometric_T_lewis_ricci,
            bolometric_T_lewis_redshift,
            bolometric_T_cpv_lensing,
            temperature_correction_normalized,
            normalization,
            l*(l+1.)/(2.*_PI_)*ba.T_cmb*1e6*cl_Xz_short,
            l*(l+1.)/(2.*_PI_)*ba.T_cmb*1e6*cl_Xz_long,
            l*(l+1.)/(2.*_PI_)*pow(ba.T_cmb*1e6,2)*cl_XX_short,
            l*(l+1.)/(2.*_PI_)*pow(ba.T_cmb*1e6,2)*cl_XX_long
            // -5*komatsu_coefficient*(Cl1*Cl2 + Cl1*Cl3 + Cl2*Cl3)
            // -5 * cl_short/cl_zeta_short              /* Difference between lewis_ricci and analytical_cpv_no_lensing */
            );
          
        } // end of if(PRINT_SQUEEZED_SMART_SCALE || PRINT_SQUEEZED_LARGE_SCALE)
          
        fprintf(stderr, "\n");
  
      } // end of for(index_l1)
    } // end of for(index_l2)
  } // end of for(index_l3)
  
  
  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================
  
  for (i=0; i<argc_for_SONG; ++i) free (argv_for_SONG[i]);
  free (argv_for_SONG);
  
  
  if (bispectra_free(&pt,&sp,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra_free \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }
  
  if (pt.has_cls == _TRUE_) {
    if (spectra_free(&sp) == _FAILURE_) {
      printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
      return _FAILURE_;
    }
  }
  
  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
  
  if (transfer2_free(&pt2,&tr2) == _FAILURE_) {
    printf("\n\nError in transfer2_free \n=>%s\n",tr2.error_message);
    return _FAILURE_;
  }
  
  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  
  if (bessel2_free(&pr,&pr2,&bs,&bs2) == _FAILURE_)  {
    printf("\n\nError in bessel2_free \n=>%s\n",bs2.error_message);
    return _FAILURE_;
  }
  
  if (bessel_free(&pr,&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
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




