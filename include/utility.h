/** @file utility.h Header file for utility.c */

#ifndef __UTILITY__
#define __UTILITY__

#include "song.h"

/****************************************************************/

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
       struct lensing * ple,
       ErrorMsg error_message
  );

#ifdef __cplusplus
  }
#endif

/**************************************************************/

#endif
