#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions (1st-order) */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions (1st-order) */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra (1st-order) */
  struct bispectra bi;        /* for bispectra */
  struct fisher fi;           /* for fisher matrix */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&bi,&fi,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }
  
  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
  
  /* Compute C_l's (lensed and unlensed). If we don't need the lensed C_l's
  all the way to l_max, then execute the standard CLASS moduels. Otherwise
  call the function 'compute_cls' which extends l_max to l_max + delta_l_max. */
  if (pr.use_lensed_cls_in_fisher == _FALSE_) {
  
    if (spectra_init(&pr,&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
      printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
      return _FAILURE_;
    }
  
    if (nonlinear_init(&pr,&ba,&th,&pm,&sp,&nl) == _FAILURE_) {
      printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }
  
    if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
      printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
      return _FAILURE_;
    }
  }
  else {

    if (compute_cls (&pr,&ba,&th,&sp,&nl,&le,errmsg) == _FAILURE_) {
      printf("\n\nError in compute_cls \n=>%s\n",errmsg);
      return _FAILURE_;
    }
  }  
  
  if (bispectra_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi) == _FAILURE_) {
    printf("\n\nError in bispectra_init \n=>%s\n",bi.error_message);
    return _FAILURE_;
  }
     
  if (fisher_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&bi,&fi) == _FAILURE_) {
    printf("\n\nError in fisher_init \n=>%s\n",fi.error_message);
    return _FAILURE_;
  }
  
  if (output_init(&ba,&pt,&sp,&nl,&le,&bi,&fi,&op) == _FAILURE_) {
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
  
  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }
  
  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
  
  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
    
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
