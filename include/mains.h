/** @file mains.h */

#ifndef __MAINS__
#define __MAINS__

/* class modules */
#include "common.h"
#include "input.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "perturbations2.h"
#include "threej.h"
#include "bessel.h"
#include "transfer.h"
#include "primordial.h"
#include "spectra.h"
#include "nonlinear.h"
#include "lensing.h"
#include "fisher.h"
#include "output.h"


/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
  extern "C" {
#endif

  int compute_cls(
       struct precision * ppr,
       struct background * pba,
       struct thermo * pth,
       struct spectra * psp,
       struct nonlinear * pnl,
       struct lensing * ple,
       ErrorMsg error_message
  );

#ifdef __cplusplus
  }
#endif

/**************************************************************/

#endif
