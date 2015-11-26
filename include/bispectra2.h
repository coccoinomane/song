/** @file bispectra2.h Documented header file for the intrinsic bispectra module */

#ifndef __BISPECTRA2__
#define __BISPECTRA2__

#include "bispectra.h"
#include "bessel2.h"
#include "perturbations2.h"
#include "transfer2.h"


/**
 * Workspace that contains the intermediate results for the integration of an intrinsic
 * bispectrum.
 *
 */

struct bispectra_workspace_intrinsic {

  /* Considered bispectrum XYZ (e.g. TTT, TEE, ...) */
  int X;
  int Y;
  int Z;

  /* Array that relates the bispectrum field indices (T,E,B) to the index of the second order
  transfer functions */
  int index_tt2_of_bf[_MAX_NUM_FIELDS_];

  /* What is the parity of the currently considered bispectrum? Even or odd? */
  int bispectrum_parity;

  /* Support points for the second-order transfer function in the smooth directions k1 and k2 */
  double * k_smooth_grid;
  int k_smooth_size;

  /* Temporary pointer to store arrays that have pwb->k_smooth_size elements (one per thread) */
  double ** f;

  /* Window function for the interpolation of the second-order transfer function. It is indexed as
  pwb->k_window[index_k] where 'index_k' indexes pwb->k_smooth_grid. Its size is pwb->k_smooth_size. */
  double * k_window;

  /* Inverse window function for the interpolation of the second-order transfer function. It is indexed as
  pwb->k_window[index_k] where 'index_k' indexes ptr->k. Its size is ptr->q_size. */
  double * k_window_inverse;

  /* Grid in the integration variable 'r'.  This is the parameter that stems from the Rayleigh expansion 
  of the Dirac Delta \delta(\vec{k1}+\vec{k2}+\vec{k3}) */
  double r_min;
  double r_max;
  int r_size;
  double * r;


  /* Array to contain the integral over k3:
  
                             /
     INT_l3 (r,k1,k2)   =   |  dk3  k3^2 * j_l3(r*k3) * T_l3(k3,k1,k2)
                            /

    The array is indexed as pbi->integral_over_k3[index_l3][index_r][index_k1][index_k2]  */
  double **** integral_over_k3;
                           
  /* Integration grid in k3 for a given k1 and k3, one for each thread: k3_grid[thread][index_k3] */
  double ** k3_grid;

  /* Factor by which the second-order transfer function has to be rescaled in order to obtain the
  m!=0 contributions. It is indexed as pwb->T_rescaling_factor[thread][index_k3]. */
  double ** T_rescaling_factor;
  double M3_coefficient[_MAX_NUM_AZIMUTHAL_]; /* Coefficient entering the rescaling factor */


  /* Array to contain the integral over k2:  
  
                             /
     INT_l2_l3 (r,k1)   =   |  dk2  k2^2 * j_l2(r*k2) * INT_l3 (r,k1,k2) * T_l2(k2)
                            /

    The array is indexed as pbi->integral_over_k2[index_l2][index_l3][index_r][index_k1] */
  double **** integral_over_k2;



  /* Array to contain the integral over k1:

                             /
     INT_l1_l2_l3 (r)   =   |  dk1  k1^2 * j_l1(r*k1) * INT_l2_l3 (r,k1) * T_l1(k1)
                            /

    The array is indexed as pbi->integral_over_k1[index_l1][index_l2][index_l3-index_l_triangular_min][index_r] */
  double **** integral_over_k1;



  /* Array to contain the integral over r:
  
                          /
     INT_l1_l2_l3    =   |  dr  r^2 * INT_l1_l2_l3(r)
                         /

    The array is indexed as pbi->integral_over_r[index_l1][index_l2][index_l3-index_l_triangular_min] */
  double *** integral_over_r;


  /* Array to contain the unsymmetrised bispectrum. This is basically the integral over r times messy
  geometrical factors summed over all M3,L3,L1 configurations. 
  Indexed as pbi->unsymmetrised_bispectrum[index_l1][index_l2][index_l3-index_l_triangular_min], 
  The ordering of its 6 levels is such that:
  < X^(2)_l1 Y^(1)_l2 Z^(1)_l3 > = pwb->unsymmetrised_bispectrum[X][Y][Z][l1][l2][l3]
  that is, the second-order transfer function always corresponds to the first field and to the
  the first multipole index of the unsymmetrised bispectrum array. */
  double ****** unsymmetrised_bispectrum;
  

  /* Array that contains the interpolated values of the above integrals in ptr->k. Each thread has one.
  Indexed as integral_splines[thread][index_k] and interpolated_integral[thread][index_k], where
  index_k belongs to ptr->k. */
  double ** integral_splines;
  double ** interpolated_integral;
  
  /* Same as above, but for the k3 integration grid (one per thread) */
  double ** delta_k3;

  /* Array that contains the pwb->r[i+1] - pwb->r[i-1] values needed for the trapezoidal rule by 
  the function 'bispectra_smooth_integration_over_r' */
  double * delta_r;

  
  /* Keep track of the summations over M3,L3,L1 */
  int M3;
  int abs_M3;
  int offset_L3;
  int offset_L1;


  /* Memory counters for the various arrays */
  long int count_allocated_for_unsymmetrised_bispectrum;
  long int count_allocated_for_integral_over_k1;
  long int count_allocated_for_integral_over_k2;
  long int count_allocated_for_integral_over_k3;
  long int count_allocated_for_integral_over_r;

  long int count_memorised_for_unsymmetrised_bispectrum;
  long int count_memorised_for_integral_over_k1;
  long int count_memorised_for_integral_over_k2;
  long int count_memorised_for_integral_over_k3;
  long int count_memorised_for_integral_over_r;



};










/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif


  int bispectra2_init(
       struct precision * ppr,
       struct precision2 * ppr2,
       struct background * pba,
       struct thermo *pth,
       struct perturbs * ppt,
       struct perturbs2 * ppt2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       struct transfers * ptr,
       struct transfers2 * ptr2,
       struct primordial * ppm,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi
       );

  int bispectra2_harmonic(
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
      struct spectra * psp,
      struct bispectra * pbi
      );
  
  int bispectra2_intrinsic_init(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_intrinsic_workspace_init(
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
      struct spectra * psp,
      struct bispectra * pbi,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_intrinsic_workspace_free(
      struct perturbs2 * ppt2,
      struct transfers2 * ptr2,
      struct bispectra * pbi,
      struct bispectra_workspace_intrinsic * pwb
      );

  
  int bispectra2_intrinsic_integrate_over_k3(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      int index_tt2_k3,
      int index_M3,
      int offset_L3,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_interpolate_over_k2(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      int index_r,
      int index_k1,
      int index_l3,
      double * integral_splines,
      double * interpolated_integral,
      double * f,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_intrinsic_integrate_over_k2(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      int index_tt_k2,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_interpolate_over_k1(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      int index_r,
      int index_l3,
      int index_l2,
      double * integral_splines,
      double * interpolated_integral,
      double * f,
      struct bispectra_workspace_intrinsic * pwb
      );
  

  int bispectra2_intrinsic_integrate_over_k1(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      int index_tt_k1,
      int offset_L1,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_intrinsic_integrate_over_r(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi,
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_intrinsic_geometrical_factors(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,      
      struct bispectra * pbi,
      int index_bt,
      int index_M3,
      int offset_L3,
      int offset_L1,
      double *** unsymmetrised_bispectrum, /* out */
      struct bispectra_workspace_intrinsic * pwb
      );

  int bispectra2_add_quadratic_corrections (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct background * pba,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct primordial * ppm,
      struct spectra * psp,
      struct bispectra * pbi
      );

  
#ifdef __cplusplus
}
#endif

#endif
