#ifndef __CRYSTAL__
#define __CRYSTAL__

#include <bispectra.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_coupling.h>


/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int mesh_int (
      int type,
      struct bispectra * pbi,
      int index_bt,
      double l1,
      double l2,
      double l3,
      double * interpolated_value
      );

  int mesh_sort (
      int type,
      struct bispectra * pbi,
      int index_bt,
      double ** vals
      );

  double distance (double * vec1, double * vec2);

#ifdef __cplusplus
}
#endif

#endif