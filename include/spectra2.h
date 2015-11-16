/** @file spectra2.h Documented includes for intrinsic spectra module */

#ifndef __SPECTRA2__
#define __SPECTRA2__

#include "transfer2.h"
#include "spectra.h"

/**
 * Macro to access the second-order sources
 */
#undef sources
#define sources(INDEX_TAU,INDEX_K_TRIANGULAR) \
  ppt2->sources[index_tp2]\
               [index_k1]\
               [index_k2]\
               [(INDEX_TAU)*k_pt_size + (INDEX_K_TRIANGULAR)]


/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int spectra2_init(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct background * pba,
        struct thermo * pth,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers * ptr,
        struct transfers2 * ptr2,
        struct primordial * ppm,
        struct lensing * ple,
        struct bispectra * pbi,
        struct spectra * psp
        );

  int spectra2_cls (
       struct precision * ppr,
       struct precision2 * ppr2,
       struct background * pba,
       struct thermo * pth,
       struct perturbs * ppt,
       struct perturbs2 * ppt2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       struct transfers * ptr,
       struct transfers2 * ptr2,
       struct primordial * ppm,
       struct lensing * ple,
       struct bispectra * pbi,
       struct spectra * psp
       );

  int spectra2_pks (
       struct precision * ppr,
       struct precision2 * ppr2,
       struct background * pba,
       struct thermo * pth,
       struct perturbs * ppt,
       struct perturbs2 * ppt2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       struct transfers * ptr,
       struct transfers2 * ptr2,
       struct primordial * ppm,
       struct lensing * ple,
       struct bispectra * pbi,
       struct spectra * psp
       );

  int spectra2_integrate_fourier_sym(
  			struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct primordial * ppm,
        struct perturbs2 * ppt2,
        struct spectra * psp
  			);
			
  int spectra2_integrate_fourier(
  			struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct primordial * ppm,
        struct perturbs2 * ppt2,
        struct spectra * psp
  			);			

  int spectra2_interpolate_sources_in_k(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct spectra * psp,
        int index_k1,
        int index_k2,
        int index_tp2,
        double * k_grid,
        double * sources_k_spline,
        double * interpolated_sources_in_k
        );
      

  int spectra2_interpolate_sources_in_k2(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct spectra * psp,
        int index_kt1,
        int index_kt3,
        int index_tp2,
        double * k_grid,
        double * sources_k_spline,
        double * interpolated_sources_in_k
        );

  int spectra2_get_k3_size (
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct spectra * psp
        );
      
#ifdef __cplusplus
}
#endif

#endif
