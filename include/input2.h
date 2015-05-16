/** @file input.h Documented includes for input module */

#ifndef __INPUT2__
#define __INPUT2__

#include "input.h"
#include "common2.h"
#include "bessel2.h"
#include "perturbations2.h"
#include "transfer2.h"


/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input2_init_from_arguments(
		 int argc, 
		 char **argv,
		 struct precision * ppr,
		 struct precision2 * ppr2,
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
     struct perturbs2 *ppt2,
		 struct transfers *ptr,
		 struct bessels *pbs,
		 struct bessels2 *pbs2,
     struct transfers2 *ptr2,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct bispectra *pbi,
     struct fisher *pfi,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input2_init(
		 struct file_content * pfc,
		 struct precision * ppr,
		 struct precision2 * ppr2,
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
     struct perturbs2 *ppt2,
		 struct transfers *ptr,
		 struct bessels * pbs,
		 struct bessels2 * pbs2,
     struct transfers2 *ptr2,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct bispectra *pbi,
     struct fisher *pfi,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input2_default_params(
			   struct background *pba,
			   struct thermo *pth,
			   struct perturbs *ppt,  
         struct perturbs2 *ppt2,
			   struct transfers *ptr,
			   struct bessels * pbs,
			   struct bessels2 * pbs2,
         struct transfers2 *ptr2,
			   struct primordial *ppm,
			   struct spectra *psp,
			   struct nonlinear *pnl,
			   struct lensing *ple,
			   struct bispectra *pbi,
         struct fisher *pfi,
			   struct output *pop
			   );
  
  int input2_default_precision(
			      struct precision2 * ppr2
			      );

  int input2_free(
            struct precision2 * ppr2
            );

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
