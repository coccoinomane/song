/** @file song.c
 *
 * Main executable for SONG.
 *
 * Last modified by Guido W. Pettinari, 26.04.2015
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

  /* Read parameters from input files */
  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,
    &pt,&tr,&pm,&sp,&nl,&le,&bs,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,
    &pt,&pt2,&tr,&bs,&bs2,&tr2,&pm, &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  /* This file is meant only for computations that involve second-order perturbations */
  if (pt2.has_perturbations2 == _FALSE_) {
    printf ("\nThe computation you requested is linear. Use 'class' rather than 'song'.\n");
    return _FAILURE_;
  }

  /* Compute background quantities */
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  /* Compute recombination and reionisation quantities */
  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  /* Compute the first-order C_l */
  if (pt.has_cls || pt.has_pks) {
    if (compute_cls (&pr,&ba,&th,&pt,&sp,&le,errmsg) == _FAILURE_) {
      printf("\n\nError in compute_cls \n=>%s\n",errmsg);
      return _FAILURE_;
    }
  }

  /* Compute first and second-order perturbations */
  if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  /* Compute primordial power spectrum from inflation */
  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  /* Compute nonlinear corrections */
  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  /* Compute first-order transfer functions using the line of sight formalism */
  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  
  /* Compute geometrical factors needed for the bispectrum integration */
  if (bessel_init(&pr,&ba,&th,&tr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }
  
  /* Compute geometrical factors needed for the line of sight integration
  at second order */
  if (bessel2_init(&pr,&pr2,&pt2,&bs,&bs2) == _FAILURE_) {
    printf("\n\nError in bessel2_init \n =>%s\n",bs2.error_message);
    return _FAILURE_;
  }
  
  /* Compute second-order transfer functions using the line of sight formalism */
  if (transfer2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2) == _FAILURE_) {
    printf("\n\nError in transfer2_init \n=>%s\n",tr2.error_message);
    return _FAILURE_;
  }
  
  /* Compute bispectra */
  if (bispectra_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra_init \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }
  
  /* Compute the intrinsic bispectrum */
  if (bispectra2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra2_init \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }
  
  /* Compute the intrinsic C_l */
  if (spectra2_init(&pr,&pr2,&ba,&th,&pt,&pt2,&bs,&bs2,&tr,&tr2,&pm,&le,&bi,&sp) == _FAILURE_) {
    printf("\n\nError in spectra2_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  /* Compute Fisher matrix */
  if (fisher_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi,&fi) == _FAILURE_) {
    printf("\n\nError in fisher_init \n=>%s\n",fi.error_message);
    return _FAILURE_;
  }
  
  /* Write output files */
  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&bi,&fi,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }


  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================
  
  if (fisher_free(&bi,&fi) == _FAILURE_) {
    printf("\n\nError in fisher_free \n=>%s\n",fi.error_message);
    return _FAILURE_;
  }

  if (bispectra_free(&pr,&pt,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra_free \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (pt.has_cls && spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
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

  if (input2_free(&pr2) == _FAILURE_) {
    printf("\n\nError in input2_free \n=>%s\n",pr2.error_message);
    return _FAILURE_;
  }

  parser_free(pr.parameter_files_content);
  free (pr.parameter_files_content);

  return _SUCCESS_;

}
