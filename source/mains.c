#include "mains.h"

int compute_cls(
     struct precision * ppr,
     struct background * pba,
     struct thermo * pth,
     // struct perturbs * ppt,
     // struct bessels * pbs,
     // struct transfers * ptr,
     // struct primordial * ppm,
     struct spectra * psp,
     struct nonlinear * pnl,
     struct lensing * ple,
     // struct bispectra * pbi,
     // struct fisher * pfi,
     // struct output * pop,
     ErrorMsg error_message
     )
{

  /* Structures that are not of interest for SONG */
  struct precision pr;
  struct perturbs pt;
  struct bessels bs;
  struct transfers tr;
  struct primordial pm;
  struct bispectra bi;
  struct fisher fi;
  struct output op;


  // ======================================================================================
  // =                                Modify input parameters                             =
  // ======================================================================================
  
  /* Structure that contains the two input parameter files */
  struct file_content * pfc = ppr->input_file_content;

  if (ppr->use_lensed_cls_in_fisher == _TRUE_) {

    /* Read old value of l_max */
    int found = _FALSE_;
    int old_l_max_scalars;
    class_call (parser_read_int (pfc,"l_max_scalars", &old_l_max_scalars, &found, error_message),
	       error_message,
	       error_message);
    class_test (found == _FALSE_,
      error_message,
      "make sure that the parameter 'l_max_scalar' is set in the input file");
         
    /* Increase l_max so that the lensing C_l's are computed up to l_max */
    int new_l_max_scalars = old_l_max_scalars + ppr->delta_l_max;
    char buffer[_ARGUMENT_LENGTH_MAX_];
    sprintf (buffer, "%d", new_l_max_scalars);    
    class_call (parser_overwrite_entry (pfc, "l_max_scalars", buffer, error_message),
      error_message,
      error_message);
      
    if ((psp->spectra_verbose>0) || (ple->lensing_verbose>0))
      printf ("Will increase l_max for lensed C_l's\n");
  }

  /* Do not create the run directory twice */ 
  class_call (parser_overwrite_entry (pfc, "store_run", "no", error_message),
    error_message,
    error_message);
  
  /* Re-read the parameters, this time from the modified 'file_content' structure. */
  class_call (input_init (pfc,&pr,pba,pth,&pt,&bs,&tr,&pm,psp,&bi,&fi,pnl,ple,&op,error_message),
    error_message,
    error_message);

  /* Turn off the second-order perturbations by setting the flag 'pt.has_perturbations2' to
  _FALSE_. In this way, SONG is ignored and the standard CLASS will compute the linear C_l's
  for the parameters specified in the input files. */    
  int had_perturbations2 = pt.has_perturbations2;
  pt.has_perturbations2 = _FALSE_;
  
  
  // ======================================================================================
  // =                                      Run CLASS                                     =
  // ======================================================================================
    
  pt.perturbations_verbose = 0;
  class_call (perturb_init(&pr,pba,pth,&pt),
    error_message,
    error_message);

  bs.bessels_verbose = 0;
  class_call (bessel_init(&pr,&bs),
    error_message,
    error_message);

  tr.transfer_verbose = 0;
  class_call (transfer_init(&pr,pba,pth,&pt,&bs,&tr),
    error_message,
    error_message);
  
  pm.primordial_verbose = 0;
  class_call (primordial_init(&pr,&pt,&pm),
    error_message,
    error_message);

  class_call (spectra_init(&pr,pba,&pt,&tr,&pm,psp),
    error_message,
    error_message);
  
  class_call (nonlinear_init(&pr,pba,pth,&pm,psp,pnl),
    error_message,
    error_message);
      
  class_call (lensing_init(&pr,&pt,psp,pnl,ple),
    error_message,
    error_message);


  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================
  
  class_call (primordial_free(&pm),
    error_message,
    error_message);
    
  class_call (transfer_free(&tr),
    error_message,
    error_message);
  
  class_call (bessel_free(&pr,&bs),
    error_message,
    error_message);
  
  class_call (perturb_free(&pr,&pt),
    error_message,
    error_message);
         
  return _SUCCESS_;
  
}