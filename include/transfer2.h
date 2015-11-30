/** @file transfer2.h Documented header file for the second-order transfer module. */

#ifndef __TRANSFER2__
#define __TRANSFER2__

#include "transfer.h"
#include "perturbations.h"
#include "perturbations2.h"
#include "bessel.h"
#include "bessel2.h"

/**
 * Method to determine the sampling of the k3 direction of the second-order transfer
 * functions. This is the direction that is integrated in the line of sight integral
 * (k in 5.95 of http://arxiv.org/abs/1405.2280). All methods use the same grid in
 * ppt2->k3 up to a certain k3, and afterwards they use a finer grid.
 */
enum transfer2_k3_sampling {
  bessel_k3_sampling, /**< add points to grid based on x-sampling of the Bessel functions */
  class_transfer2_k3_sampling   /**< add points to grid based on input from user */
};

/**
 * Method to determine the integration grid in time for line of sight integral. All
 * methods use the grid in ppt2->tau_sampling plus some extra points.
 */
enum transfer2_tau_sampling {
  sources_tau_sampling, /**< Use the same sampling as the line-of-sight sources (ppt2->tau_sampling). This is the only
                        option that does not requires interpolation of the sources. */ 
  custom_tau_sampling,  /**< Add extra points to the sources sampling based on input from user */
  bessel_tau_sampling,  /**< Add extra points to the sources sampling based on the x-sampling of the Bessel functions */
};


/** 
 * Macro used to index the first level ptr2->transfer.
 */
#define lm_cls(index_l,index_m) ptr2->lm_array[index_l][index_m]


/**
 * In order to access the sources for the line of sight integration,
 * we define a preprocessor macro that takes as arguments the time and
 * k3 indices.
 */
#undef sources
#define sources(INDEX_TAU,INDEX_K_TRIANGULAR) \
  ppt2->sources[index_tp2]\
               [index_k1]\
               [index_k2]\
               [(INDEX_TAU)*k_pt_size + (INDEX_K_TRIANGULAR)]


/**
 * Maximum number of transfer function types computed in this module.
 *
 * Feel free to increase it, it is just a memory parameter.
 */
#define _MAX_NUM_TRANSFERS_ 4096


/**
 * Structure containing everything about the second-order transfer functions in
 * harmonic space \f$ \T_l^{X} (k1,k2,k3) \f$ that other modules need to know.
 *
 * Once initialized by transfer2_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested azimuthal modes m, type (temperature, polarization, etc), multipole l
 * and wavenumber (k1,k2,k3).
 * 
 * The content of this structure is entirely computed in the transfer2.c module,
 * given the content of the 'precision', 'background', 'thermodynamics',
 * 'perturbation', 'transfer' and 'bessel' structures.
 */

struct transfers2 {


  // ==========================================================================================
  // =                                    Transfer functions                                  =
  // ==========================================================================================

  /** 
   * Array containing the transfer function T_lm(k1,k2,k3) for all required types (T,E,B).
   *
   * The ptr2->transfer array should be indexed as follows:
   * 
   *     ptr2->transfer [index_tt2]
   *                    [index_k1]
   *                    [index_k2]
   *                    [index_k]
   * 
   * - index_tt2 is a composite index that includes both the field (x=T,E,B...) and the multipole
   *   (l,m); it is expanded as index_tt2 = ptr2->index_tt2_X + ptr2->lm_array[index_l][index_m].
   * - index_k1 goes from 0 to ppt2->k_size-1.    
   * - index_k2 goes from 0 to ppt2->k_size-index_k1-1 due to symmetry reasons.
   * - index_k goes from 0 to ptr2->k_size_k1k2[index_k1][index_k2]-1.
   */

  double **** transfer; 
  
  short * transfers_available;  /**< If transfers_available[index_tt] is true, then
                                ptr2->transfer[index_tt] has been filled and is ready
                                for use */

  short * transfers_allocated;  /**< If transfers_allocated[index_tt] is true, then
                                ptr2->transfer[index_tt] is fully allocated and can
                                be used to store the transfer function. */

  
  int index_tt2_T;              /* Index for transfer type = temperature */
  int index_tt2_E;              /* Index for transfer type = E-polarization */
  int index_tt2_B;              /* Index for transfer type = B-polarization */

  int tt2_size;                 /* Number of requested transfer types */

  /* Array of strings that contain the labels of the various transfer types
  For example,  tt2_labels[index_tt2_T] is equal to "T_00" */
  char (*tt2_labels)[_MAX_LENGTH_LABEL_];

  /* True if the module has been executed. Useful to free memory only if needed. */
  short has_cls;



  // ==============================================================================
  // =                                  Multipoles                                =
  // ==============================================================================

  // *** Number and list of multipoles

  int l_size;        /* number of multipole values */
  int * l;           /* list of multipole values l[index_l] */
  int * m;           /* list of azimuthal multipole values m[index_m] */
  int m_size;        /* number of of azimuthal multipole values m[index_m] */
  int n_transfers;   /* number of possible (l,m) combinations attainable in ptr2->l and ptr2-m */

  int n_nonzero_transfers_E;   /**< Number of nonzero transfers to be computed for photon E-polarization; this is basically
                                    ptr2->n_transfers minus the l=0 and l=1 modes */
  int n_nonzero_transfers_B;   /**< Number of nonzero transfers to be computed for photon B-polarization; this is basically
                                    ptr2->n_transfers minus the l=0, l=1 and m=0 modes */

  // *** Correspondance between (type,l,m) and index_tt

  int * tt2_to_index_l;    /**< Array with the correspondence between ptr2->index_tt2_XXX
                           index and l multipole index */
  int * tt2_to_index_m;    /**< Array with the correspondence between ptr2->index_tt2_XXX
                           index and m multipole index */
  int * tt2_start;         /**< Array with the correspondence between ptr2->index_tt2_XXX
                           index and the start of the hierarchy to which the index belongs */
  int * tp2_start;         /**< Array with the correspondence between ptr2->index_tt2_XXX
                           index and the start of the hierarchy to which the index belongs 
                           in the source array */


  /* 'lm_array[index_l][index_m]' contains the index associated with a given (l,m) couple.
  To extract the (l,m) multipole of the temperature transfer function, do this:
  ptr2->transfer[index_monopole_t + lm_array[index_l][index_m]].
  Note that this array is different than ppt2->lm_array. In the perturbations module we
  needed only all the l's up to a certain value (usually l_max ~ 10) for the Boltzmann hierarchy,
  while here we need many sparsely sampled l's to compute the Cl's or the bispetrcum
  (usually l_max ~ 2000). */
  int ** lm_array;



  // ===============================================================================
  // =                                 k-sampling                                  =
  // ===============================================================================

  /* Which k-sampling should we use for the second-order transfer functions? */  
  enum transfer2_k3_sampling k_sampling;

  /* Which time-sampling should we use for the second-order transfer functions? */  
  enum transfer2_tau_sampling tau_sampling;

  /* For a given (k1,k2), number of considered k3 values. This number includes points outside the physical
  (i.e. triangular) regime when extrapolation is requested. */
  int ** k_size_k1k2;

  /* For a given (k1,k2), index of the first value of k3 that satisfies the triangular condition. All
  entries must be equal zero when no extrapolation is used */
  int ** k_physical_start_k1k2;

  /* For a given (k1,k2), number of k3 values that satisfy the triangular condition */
  int ** k_physical_size_k1k2;

  /* For a given (k1,k2), extrema of the k3 grid. */
  double ** k_min_k1k2;
  double ** k_max_k1k2;
  
  /* Maximum extent of the k3 grid for all possible pairs of (k1,k2) */
  int k3_size_max;
  


  // ====================================================================================
  // =                                 Disk storage                                     =
  // ====================================================================================

  /** @ingroup StorageFiles
   * Parameters related to the reading and writing of storage files.
   *
   * Refer to the documentation in perturbations2.h (\ref StorageFiles) for details.
   */
  //@{
  char storage_dir[_FILENAMESIZE_]; /**< Directory containing the transfer function storage files. If it already
                                    exists, and load_transfers is true, the transfers will be read from this folder
                                    to the array ptr2->transfer. If it does not exist, and store_transfers is true,
                                    the transfers will be first computed and then written to this folder. Either way,
                                    the directory contains one binary file for transfer type, for a total of
                                    tt2_size files. */

  char ** storage_paths; /**< storage_paths[index_tt2] is the path to the file with the transfer function
                         of type index_tt2. Used only if either of the flags store_transfers or
                         load_transfers are true. */

  FILE ** storage_files; /**< storage_files[index_tt2] is the pointer to the file with the transfer function
                         of type index_tt2. Used only if either of the flags store_transfers or
                         load_transfers are true. */
  //@}


  // ==================================================================================
  // =                            Cosmological quantities                             =
  // ==================================================================================

  double tau0;                /* Conformal age */  
  double tau_rec;             /* Conformal time at recombination */
  double rs_rec;              /* Comoving sound horizon at recombination */



  // ====================================================================================
  // =                                 Debug parameters                                 =
  // ====================================================================================

  short transfer2_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  /* Count the number of values in ppt2->sources as we fill the array (debug only) */
  long int count_allocated_transfers;
  long int count_memorised_transfers;

  /* Logical array. If the index_k1 position is true, then ppt2->transfers[index_k1] is allocated */
  short * has_allocated_transfers;

  short stop_at_transfers2; /**< If _TRUE_, SONG will stop execution after having run the transfer2.c
                            module. Useful to debug the second-order transfer functions today. */
  

};








/**
 * Just a collection of often-used, temporary parameters that are passed through
 * various functions in the transfer2 module.  Each set of (l,m,k1,k2,k) for which
 * we compute the transfer functions has a reserved workspace.
 */
struct transfer2_workspace {

  /* Value of the transfer function for a given set of (l,m,k1,k2,k) */
  double transfer;
  
  /* Index in ptr2->transfer of the transfer type we are working on */
  int index_tt;
  
  /* Indexes and values of the set (l,m,k1,k2,k) where we are computing the transfer function */
  int index_k;
  int index_k1;
  int index_k2;
  double k;
  double k1;
  double k2;
  double cosk1k2;
  
  /* Integration grid. All the following arrays have size tau_grid_size */
  int tau_grid_size;
  double * tau_grid;
  double * tau0_minus_tau;        /* List of tau0-tau values, tau0_minus_tau[index_tau_grid] */
  double * delta_tau;             /* List of delta_tau values for trapezoidal rule, delta_tau[index_tau_grid] */

  /* Sampling in k where we shall compute the transfer function */
  double * k_grid;

  /* Array that will contain the interpolated value of the sources at
  the times contained in the integration grid. */
  double ** interpolated_sources_in_time;

  /* We shall interpolate the sources in ptr2->k and in pw->tau_grid. In order to speed up the table look-up for the
  time interpolation, we define pw->index_tau_left[index_tau_grid] which gives the index of the element of ppt2->tau_sampling
  closest (to the left) to pw->tau_grid[index_tau_grid] */
  int * index_tau_left;

  /* How many line of sight sources should be kept? This is updated according to whether we consider
  temperature (ppr2->l_max_los_t) or polarization (ppr->l_max_los_p). */
  int L_max;

  /* Number of the thread that is currently using this workspace */
  int thread;

};


// ------------------------------------------------------------------------------------
// -                                 Debug shortcuts                                  -
// ------------------------------------------------------------------------------------


#undef fprintf_k_debug
#define fprintf_k_debug(args...) {                                  \
  if((pw->index_k1==ppt2->index_k1_debug) &&                      \
     (pw->index_k2==ppt2->index_k2_debug) &&                      \
     (pw->index_k==ppt2->index_k3_debug)) {                      \
    fprintf (args);                                                 \
  }                                                                 \
}

#undef printf_k_debug
#define printf_k_debug(args...) {                                   \
  fprintf_k_debug (stdout, args)                                    \
}


/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int transfer2_init(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct background * pba,
        struct thermo * pth,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,    
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers * ptr,
        struct transfers2 * ptr2
        );
    
  int transfer2_free(
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct transfers2 * ptr2
          );
  
  int transfer2_indices_of_transfers(
              struct precision * ppr,
              struct precision2 * ppr2,
              struct perturbs2 * ppt2,
              struct bessels * pbs,
              struct bessels2 * pbs2,
              struct transfers * ptr,
              struct transfers2 * ptr2
              );
  
  int transfer2_get_l_list(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2
        );
  
  
  int transfer2_get_lm_lists(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct bessels * pbs,
          struct bessels2 * pbs2,
          struct transfers2 * ptr2
          );



  int transfer2_get_k3_sizes(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct bessels * pbs,
          struct bessels2 * pbs2,
          struct transfers * ptr,
          struct transfers2 * ptr2
          );

  
  int transfer2_get_k3_sizes(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct bessels * pbs,
          struct bessels2 * pbs2,
          struct transfers * ptr,
          struct transfers2 * ptr2
          );

  int transfer2_interpolate_sources_in_k(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1,
        int index_k2,
        int index_tp,
        double * k_grid,
        double * interpolated_sources
        );


  int transfer2_get_time_grid(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1,
        int index_k2,
        int index_k,
        struct transfer2_workspace * pw
        );

  int transfer2_interpolate_sources_in_time(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs * ppt,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_tp2,                          
        double * interpolated_sources_in_k,
        double * interpolated_sources_in_time,
        struct transfer2_workspace * pw
        );
  
  int transfer2_get_k3_size(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1,
        int index_k2
        );

  
  int transfer2_get_k3_list(
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1,
        int index_k2,
        double * k3,
        int * last_used_index_pt
        );

  int transfer2_compute (
          struct precision * ppr,
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct bessels * pbs,
          struct bessels2 * pbs2,
          struct transfers2 * ptr2,
          int index_k1,
          int index_k2,
          int index_k,
          int index_l,  
          int index_m,
          int index_tt,
          double ** interpolated_sources_in_time,
          struct transfer2_workspace * pw
          );

  int transfer2_integrate (
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1,
        int index_k2,
        int index_k,
        int index_l,  
        int index_m,
        double ** interpolated_sources_in_time,
        int index_J,
        int index_source_monopole,
        struct transfer2_workspace * pw,
        double * integral
        );


  int transfer2_store(
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct transfers2 * ptr2,
        int index_k1
        );

  int transfer2_load(
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct transfers2 * ptr2,
        int index_tt
        );


  int transfer2_allocate_type_level(
       struct perturbs2 * ppt2,
       struct transfers2 * ptr2,
       int index_tt
       );

  int transfer2_free_type_level(
       struct perturbs2 * ppt2,
       struct transfers2 * ptr2,
       int index_tt
       );

  int transfer2_allocate_k1_level(
       struct perturbs2 * ppt2,
       struct transfers2 * ptr2,
       int index_k1
       );

  int transfer2_free_k1_level(
       struct perturbs2 * ppt2,
       struct transfers2 * ptr2,
       int index_k1
       );


#ifdef __cplusplus
}
#endif

#endif
