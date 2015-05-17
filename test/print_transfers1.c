/** @file print_transfers2.c 
 * Created by Guido W. Pettinari on 17.07.2012
 * Last edited by Guido W. Pettinari on 17.07.2012
 *
 * Print to screen the transfer array inside the transfer structure, with k varying.
 *
 * usage:     print_transfers <ini file> [<pre file>] <n_columns>
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
	int n_columns;
	
	/* CLASS can accept either one or two arguments */
  if (argc == 2 + n_args) {
    n_columns = atoi(argv[2]);
  }
  else if (argc == 3 + n_args) {
    n_columns = atoi(argv[3]);
  }
  else {
    printf ("usage:     %s <ini file> <pre file> <n_columns>\n", argv[0]);
    printf ("           %s <run_directory> <n_columns>\n", argv[0]);
    return _FAILURE_;
  }


  // ===========================
  // = Calculate perturbations =
  // ===========================

  /* Decrease the argument counter. The reason is that CLASS should be fed only its default arguments,
  that is the parameter files and, optionally, the run directory */
  argc -= n_args;

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  /* Compute the second-order sources and transfer functions no matter what is written in params.ini */
  pt2.has_perturbations2 = _TRUE_;
	pt.has_cls = pt2.has_cls = tr2.has_cls = _TRUE_;
	pt.has_cl_cmb_temperature = pt2.has_cmb_temperature = _TRUE_;
  bi.has_bispectra = _FALSE_;

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (pt.has_perturbations2 == _TRUE_) {
    if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
      printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
      return _FAILURE_;
    }
  }
  else {
    if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
      printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
      return _FAILURE_;
    } 
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (pt.has_perturbations2 == _TRUE_) {
    if (bessel2_init(&pr,&pr2,&pt2,&bs,&bs2) == _FAILURE_) {
      printf("\n\nError in bessel2_init \n =>%s\n",bs2.error_message);
      return _FAILURE_;
    }
  }

  if (pt2.k_size > 1)
    if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
      printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }

  if (pt.has_perturbations2 == _TRUE_) {
    if (transfer2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2) == _FAILURE_) {
      printf("\n\nError in transfer2_init \n=>%s\n",tr2.error_message);
      return _FAILURE_;
    }
  }
  else {
    if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
      printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }
  }
  
  /* Information about the cosmological model */
  double h = ba.h;
  double a_equality = ba.a_eq;  
  fprintf(stderr,"# Cosmological parameters:\n");
  fprintf(stderr,"# tau0 = %g, a_equality = %g, Omega_b = %g, Tcmb = %g, Omega_cdm = %g, omega_lambda = %g, Omega_ur = %g, Omega_fld = %g, h = %g, tau0 = %g\n",
    ba.conformal_age, ba.a_eq, ba.Omega0_b, ba.T_cmb, ba.Omega0_cdm, ba.Omega0_lambda, ba.Omega0_ur, ba.Omega0_fld, ba.h, ba.conformal_age);
  fprintf(stderr,"# omega_b = %g, omega_cdm = %g, omega_lambda = %g, omega_ur = %g, omega_fld = %g\n",
    ba.Omega0_b*h*h, ba.Omega0_cdm*h*h, ba.Omega0_lambda*h*h, ba.Omega0_ur*h*h, ba.Omega0_fld*h*h);

  /* Shortcuts */
  int k_size = tr.k_size[index_mode];
  int tt_size = tr.tt_size[index_mode];
  int l_size = tr.l_size[index_mode];  
  
  /* Account for overshooting of the inputs */
  if((n_columns>(tt_size*l_size)) || (n_columns<=0)) n_columns = tt_size*l_size;
  
  /* Some debug info */
  fprintf(stderr,"# gauge = ");
  if (pt.gauge == newtonian)
    fprintf(stderr,"Newtonian gauge\n");
  if (pt.gauge == synchronous)
    fprintf(stderr,"synchronous gauge\n");
  
  fprintf(stderr,"# Time-sampling of quadsources with %d points from tau=%g to %g\n",
    pt.tau_size, pt.tau_sampling[0], pt.tau_sampling[pt.tau_size-1]);
  fprintf(stderr,"# K-sampling of 1st-order transfer functions with %d points from k=%g to %g\n",
    k_size, tr.k[index_mode][0], tr.k[index_mode][k_size-1]);
  fprintf(stderr,"# Number of transfer types, multipoles, columns printed x type: %d, %d, %d\n",
    tt_size, l_size, n_columns);
  


  /* Running index used to number the columns */
  int index_print=1;

  /* First row contains the labels of the different types */
  fprintf (stderr,"%10s(%03d) ", "k",index_print++);
  char label[32];

  /* Print labels of transfer types */
  for (int index_tt = 0; index_tt < MIN(tt_size,n_columns); ++index_tt) {
    for (int index_l = 0; index_l < MIN(l_size,n_columns); ++index_l) {
      if (index_tt==tr.index_tt_t)
        sprintf (label, "T_%d", tr.l[index_l]);
      else if (index_tt==tr.index_tt_e)
        sprintf (label, "E_%d", tr.l[index_l]);
      fprintf (stderr,"%10s(%03d) ", label, index_print++);
    }
  }
  fprintf(stderr,"\n");
  

  // **************    Loop on k    *****************//
  for (int index_k = 0; index_k < k_size; ++index_k) {
  
    double k = tr.k[index_mode][index_k];
  
    /* Column containing k */
    fprintf (stderr,"%+15e ", k);
      
    /* Columns from 3 to tt2_size+2 are the transfer functions */
    for (int index_tt = 0; index_tt < MIN(tt_size,n_columns); ++index_tt) {
      for (int index_l = 0; index_l < MIN(l_size,n_columns); ++index_l) {
  
        double var = tr.transfer[index_mode] 
          [((index_ic * tt_size + index_tt) * l_size + index_l) * k_size + index_k];
        
        /* Compensate extra factor added by CLASS */
        // if (index_tt==tr.index_tt_e) {
        //   int l = tr.l[index_l];
        //   double factor = sqrt((l+2.) * (l+1.) * l * (l-1.));
        //   var /= factor;
        // }
  
        fprintf (stderr,"%+15e ", var);
      }
    }
    fprintf (stderr,"\n");
  
  } // end of for(index_k)
  
  
  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  
  if (bessel_free(&pr,&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
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