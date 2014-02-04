/** @file spectra.h Documented includes for spectra module */

#ifndef __FISHER__
#define __FISHER__

#include "bispectra.h"
#include "lensing.h"
#include "spectra.h"
#include "mesh_interpolation.h"


/* Which kind of interpolation for the bispectra?  */
enum bispectra_interpolation_method {
  smart_interpolation,
  trilinear_interpolation,
  mesh_interpolation,
  mesh_interpolation_2d,
  sum_over_all_multipoles
};


/* Maximum number of frequency bands of the experiment, for the purpose of Fisher matrix computation. */
#define _N_FREQUENCY_CHANNELS_MAX_ 100



/**
 *
 * Note that the f_NL we use is the one for the gravitational potential during matter domination, which is
 * related to the comoving curvature perturbation by a factor -3/5 (fnl_R = -3/5 fnl_psi). 
 *
 */

struct fisher {


  // =============================================================================
  // =                             Flags and indices                             =
  // =============================================================================

  short has_fisher;

  // ===============================================================================
  // =                                 Arrays                                      =
  // ===============================================================================
  
  int l_min;                 /* min value where the bispectrum is known (pbi->l[0]) */
  int l_max;                 /* max value where the bispectrum is known (pbi->l[pbi->l_size-1]) */
  int full_l_size;           /* equal to l_max - l_min + 1 */
  
  int l_min_estimator;       /* minimum l in the estimator sum, default is pbi->l[0] */
  int l_max_estimator;       /* maximum l in the estimator sum, default is pbi->l[pbi->l_size-1] */

  /* Debug variables which, by default, are set using l_min_estimator and l_max_estimator */
  int l1_min_global;
  int l2_min_global;
  int l3_min_global;
  int l1_max_global;
  int l2_max_global;
  int l3_max_global;

  /* Arrays over which the Fisher sum will run */
  int * l1;
  int l1_size;
  int * l2;
  int l2_size;
  int * l3;
  int l3_size;
  
  /* Contribution to the Fisher matrix coming from a given l1 and for a given XYZ bispectrum,
  where XYZ=TTT,TTE,TET, etc. This is the sum over l2, l3, A, B, C of
  b^XYZ(l1,l2,l3) * b^ABC(l1,l2,l3) * cov^XYZABC(l1,l2,l3), with l1>=l2>=l3.
  Indexed as pfi->fisher_matrix_XYZ_l1[X][Y][Z][index_l1][index_bt_1][index_bt_2],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ****** fisher_matrix_XYZ_l1;
  
  /* Same as above, but for l3, the smallest multipole, and l3 belonging to pfi->l3[index_l3] */
  double ****** fisher_matrix_XYZ_l3;
  
  /* Same as fisher_matrix_XYZ_l1, but summed over XYZ */
  double *** fisher_matrix_l1;

  /* Same as fisher_matrix_XYZ_l3, but summed over XYZ */
  double *** fisher_matrix_l3;

  /* Fisher matrix for the considered experiment, as a function of the angular resolution
  and for a given bispectrum XYZ. This is obtained as sum_{lmin<=l1<=lmax} fisher_matrix_XYZ_l1,
  with lmin fixed (=2) and lmax varying.
  Indexed as pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l1][index_bt_1][index_bt_2]
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ****** fisher_matrix_XYZ_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double ****** fisher_matrix_XYZ_lmin; 

  /* Fisher matrix for the considered experiment, as a function of the angular resolution.
  It is obtained as sum_{lmin<=l1<=l_max,XYZ} fisher_matrix_XYZ_l1.
  Indexed as pbi->fisher_matrix_lmax[index_l1][index_bt_1][index_bt_2],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double *** fisher_matrix_lmax;
  double *** inverse_fisher_matrix_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double *** fisher_matrix_lmin;        
  double *** inverse_fisher_matrix_lmin;

  /* Array that contains 1/sqrt(F^ii), with i=1,..,pbi->bt_size. For a given bispectrum type,
  it corresponds to the value of f_NL that could be detected by an experiment with a
  resolution of lmin<=l<=l_max, with lmin fixed (=2) and lmax varying.
  Indexed in the same way as fisher_matrix_lmax. 
  Indexed as pfi->sigma_fnl_lmax[index_l1][index_bt],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ** sigma_fnl_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double ** sigma_fnl_lmin;
    
    
  /* Array containing the quantity I_l1_l2_l3 in eq. 13 of Komatsu, Spergel & Wandelt (2005):
  
       I_l1_l2_l3 = sqrt( (2L1+1)(2L2+1)(2L3+1)/4*pi ) *  (L1 L2 L3)
                                                          (0  0  0 )
  
  which is needed to compute the Fisher matrix. (This is the factor that converts a reduced
  bispectrum to an angular averaged one.)  */
  double * I_l1_l2_l3;


  /* Cross-power spectrum of the C_l's, and its inverse; it is needed to compute the covariance
  matrix between the various bispectra. For example, if we consider temperature and polarisation
  bispectra, the full covariance matrix is a 8x8 matrix and the cross-power spectrum is given by
  C = ( C_l^TT C_l^TE )
      ( C_l^TE C_l^EE ).
  The C and inverse_C arrays are indexed as pfi->C[l-2][pfi->index_fp_X][pfi->index_fp_Y],
  where X and Y are the considered probes (T and E for temperature and polarisation) */
  double *** C;
  double *** inverse_C;


  // ==========================================================================================
  // =                                       Noise model                                      =
  // ==========================================================================================

  /* Beam for each frequency band of the considered experiment With respect to Table I of astro-ph/0506396v2,
    'beam[index_channel]' is theta_fwhm in radians for that frequency channel. */
  int n_channels;
  double beam[_N_FREQUENCY_CHANNELS_MAX_];
  
  /* Amplitude of the noise. With respect to Table I of astro-ph/0506396v2, 'noise' is 'sigma' */
  double noise_t[_N_FREQUENCY_CHANNELS_MAX_];
  double noise_e[_N_FREQUENCY_CHANNELS_MAX_];
  double noise_r[_N_FREQUENCY_CHANNELS_MAX_];

  /* Total noise as a function of l, defined by C_l_experiment = C_l_theory + N_l. This includes co-added
  contributions from all frequency channels, as explained in eq. 29 of astro-ph/0506396v2.
  It is indexed as pfi->N_l[pbi->index_bf][l-2], where pbi->index_bt=T,E... */
  double ** N_l;

  /* Sky coverage of the experiment. Equal to 1 for a full-sky experiment. */
  double f_sky;


  // ============================================================================================
  // =                                   Bispectra interpolation                                =
  // ============================================================================================
  
  /* Variable set to the type of desired interpolation (trilinear, mesh or sum) */
  enum bispectra_interpolation_method bispectra_interpolation;
    
  /* Number of meshes in which to partition the 3D l-space */
  int n_meshes;
  double * link_lengths;
  double * group_lengths;
  double * soft_coeffs;

  /* l-multipole at which we shall switch between one mesh and another. Indexed as pbi->l_turnover[index_mesh],
    where index_mesh runs from 0 to n_mesh_grids-1 */
  int * l_turnover;
  
  /* Array of pointers to interpolation workspaces, indexed as mesh_workspaces[index_bt][i][j][k][index_mesh],
  where i,j,k are the field indices (eg. TET). */
  struct mesh_interpolation_workspace ****** mesh_workspaces;

  /* C_l power spectrum to use for the bispectrum interpolation window function */
  int index_ct_window;


  // ===========================================================================================
  // =                                    Technical parameter                                  =
  // ===========================================================================================
  
  short fisher_verbose;           /* Flag regulating the amount of information sent to standard output (none if set to zero) */                                                    
  char info[8192];                /* Store Fisher matrix information to be printed out to screen and saved to file */
  ErrorMsg error_message;         /* Zone for writing error messages */

};




/**
 * Variables and arrays needed to compute the f_NL estimator given a pair of bispectra.
 *
 */
struct fisher_workspace {
  
  /* Weights for the linear interpolation. Used for the l1 and l2 sums in the Fisher matrix
  estimator. */
  double * delta_l;
  
  /* Weights for the linear interpolation along the triangular direction (l3). Indexed as
  pfi->delta_l3[thread][index_l3] */
  double ** delta_l3;

};








// /*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif


  int fisher_init(
       struct precision * ppr,
       struct background * pba,
       struct thermo *pth,
       struct perturbs * ppt,
       struct bessels * pbs,
       struct transfers * ptr,
       struct primordial * ppm,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi
       );

  int fisher_free(
       struct bispectra * pbi,
       struct fisher * pfi
       );

  int fisher_indices(
          struct precision * ppr,
          struct background * pba,
          struct perturbs * ppt,
          struct bessels * pbs,
          struct transfers * ptr,
          struct primordial * ppm,
          struct spectra * psp,
          struct lensing * ple,
          struct bispectra * pbi,
          struct fisher * pfi
          );

  int fisher_cross_cls(
          struct precision * ppr,
          struct background * pba,
          struct perturbs * ppt,
          struct bessels * pbs,
          struct transfers * ptr,
          struct primordial * ppm,
          struct spectra * psp,
          struct lensing * ple,
          struct bispectra * pbi,
          struct fisher * pfi
          );

  int fisher_noise(
          struct precision * ppr,
          struct background * pba,
          struct perturbs * ppt,
          struct bessels * pbs,
          struct transfers * ptr,
          struct primordial * ppm,
          struct spectra * psp,
          struct lensing * ple,
          struct bispectra * pbi,
          struct fisher * pfi
          );

  int fisher_create_interpolation_mesh(
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        );

  int fisher_interpolation_mesh(
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        );


  int fisher_interpolate_bispectrum (
        struct bispectra * pbi,
        struct fisher * pfi,
        int index_bt,
        int i,
        int j,
        int k,
        double l1,
        double l2,
        double l3,
        double * interpolated_value
        );

  int fisher_compute(
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        );

  int fisher_cross_correlate_mesh (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        struct fisher_workspace * pw
        );

  int fisher_cross_correlate_nodes (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        struct fisher_workspace * pw
        );
        
  int fisher_interpolation_weights (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        int index_l1,
        int index_l2,
        double * delta_l3,
        struct fisher_workspace * pw
        );
        

#ifdef __cplusplus
}
#endif

#endif
