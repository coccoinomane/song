/** @file transfer.h Documented includes for transfer module. */

#ifndef __TRANSFER2__
#define __TRANSFER2__

#include "transfer.h"
#include "perturbations.h"
#include "perturbations2.h"
#include "bessel.h"
#include "bessel2.h"


/* Which k-sampling should we use for the second-order transfer functions? */
enum transfer2_k_sampling {bessel_k_sampling, class_transfer2_k_sampling};

/* Which time-sampling should we use for the second-order transfer functions? */
enum transfer2_tau_sampling {bessel_tau_sampling, custom_transfer2_tau_sampling};


/**
 * Structure containing everything about transfer functions in harmonic space \f$ \Delta_l^{X} (k) \f$ that other modules need to know.
 *
 * Once initialized by transfer2_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, type
 * (temperature, polarization, etc), multipole l and wavenumber k.
 * 
 * The content of this structure is entirely computed in this module,
 * given the content of the 'precision', 'bessels', 'background',
 * 'thermodynamics' and 'perturbation' structures.
 */

struct transfers2 {




  // ==========================================================================================
  // =                                    Transfer functions                                  =
  // ==========================================================================================

  /* Table of transfer functions for each mode, initial condition, type, multipole and wavenumber.
  Note that the k-grid of ptr2->transfer is much finer than the k-grid of ppt2->sources, as the
  projection functions in the line-of-sight integral (basically spherical Bessel functions) make
  the transfer functions to oscillate wildly in the k direction.
    
    The ptr->transfer array should be indexed as follows:

      ptr->transfer [index_tt2_X + lm_cls(index_l,index_m)]
                    [index_k1]
                    [index_k2]
                    [index_k]

  Where X can be T, E, B, etc. for temperature, e-modes, b-modes, etc.
  Important: as in ppt2->sources, index_k2 goes from 0 to index_k1. */
  double **** transfer; 
  
  int index_tt2_T;              /* Index for transfer type = temperature */
  int index_tt2_E;              /* Index for transfer type = E-polarization */
  int index_tt2_B;              /* Index for transfer type = B-polarization */

  int tt2_size;                 /* Number of requested transfer types */


  /* Array of strings that contain the labels of the various transfer types
  For example,  tt2_labels[index_tt2_T] is equal to "T_00" */
  char ** tt2_labels;

  /* True if the module has been executed. Useful to free memory only if needed. */
  short has_cls;






  // ==============================================================================
  // =                                  Multipoles                                =
  // ==============================================================================


  // *** Number and list of multipoles

  int l_size;        /* number of multipole values */
  int l_size_max;    /* greatest of all l_size */
  int * l;           /* list of multipole values l[index_l] */
  int * m;           /* list of azimuthal multipole values m[index_m] */
  int m_size;        /* number of of azimuthal multipole values m[index_m] */
  int n_transfers;   /* number of possible (l,m) combinations attainable in ptr2->l and ptr2-m */



  // *** Correspondance between (type,l,m) and index_tt

  int * index_tt2_monopole;      /* index_tt2_monopole[index_tt] is the index in ptr2->transfer of the monopole corresponding to index_tt */
  int * index_pt2_monopole;      /* index_pt2_monopole[index_tt] is the index in ppt2->sources of the monopole corresponding to index_tt */
  int * corresponding_index_l;   /* corresponding_index_l[index_tt] is the index in ptr2->l corresponding to index_tt */
  int * corresponding_index_m;   /* corresponding_index_m[index_tt] is the index in ptr2->m corresponding to index_tt */

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
  enum transfer2_k_sampling k_sampling;

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
  






  // =================================================================================
  // =                        Storage of intermediate results                        =
  // =================================================================================

  /* Should we store the content of ptr2->transfer to disk? */
  short store_transfers_to_disk;

  /* Should we bother computing the transfer functions, or they will be loaded from disk? This flag is
    on only if both ppr->load_run and 'store_transfers_to_disk' are on. */
  short load_transfers_from_disk;



  /* Files where the transfer functions will be stored (one file for each transfer type) */
  char transfers_run_directory[_FILENAMESIZE_];
  FILE ** transfers_run_files;
  char ** transfers_run_paths;

  /* File that will keep track how how many transfer files have been succesfully written */
  FILE * transfers_status_file;
  char transfers_status_path[_FILENAMESIZE_];




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

  /* If true, compute only the 2nd-order transfer functions today, and do not care about
  other flags invoking the subsequent modules */
  short has_transfers2_only;




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

  /* Arrays that will contain the second derivatives and the interpolated value of the sources, respectively, at
    the times contained in the integration grid. */
  double ** sources_time_spline;
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
        double * splines,
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
        double * sources_time_spline,
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


  int transfer2_save_transfers_to_disk(
          struct perturbs2 * ppt2,
          struct transfers2 * ptr2,
          int index_k1
          );

  int transfer2_load_transfers_from_disk(
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
