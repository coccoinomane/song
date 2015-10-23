/** @file transfer2.c
 *
 * Module to compute and store today's value of the second-order transfer
 * functions.
 *
 * The computation consists of solving the line of sight integral, a convolution
 * in conformal time of the line of sight sources with the Bessel projection functions
 * (computed in perturbations2.c and bessel2.c, respectively). This line of sight formalism
 * allows us to infer the current value of the transfer functions today without having to
 * evolve them all the way to today.
 *
 * Every passage of the computation is documented in the source code. For
 * a more detailed explanation of the physics behind this module,
 * please refer to chapter 5 of my thesis (http://arxiv.org/abs/1405.2280).
 *
 * The main output of this module, the second-order transfer functions evaluated
 * today, is stored in the ptr2->transfer array. These are used to compute 
 * the intrinsic bispectrum in bispectra2.c. See chapter 6 of my thesis
 * (link above) for details on the bispectrum computation.
 *
 * The main functions that can be called externally are:
 * -# transfer2_init() to run the module; requires background_init(), thermodynamics_init()
 *    and perturb2_init().
 * -# transfer2_free() to free all the memory associated to the module.
 * 
 * If the user specified 'store_transfers_to_disk=yes', the module will save the 
 * transfer functions to disk after computing them, and then free the associated 
 * memory. To reload them from disk, use transfer2_load_transfers_from_disk().
 * To free again the memory associated to the sources, call
 * transfer2_free_type_level().
 * 
 * Created by Guido W. Pettinari on 04.06.2012 based on transfer.c by the CLASS
 * team (http://class-code.net/).
 * Last modified by Guido W. Pettinari on 04.06.2015.
 */

#include "transfer2.h"


/**
 * Fill all the fields in the transfers2 structure, especially the ptr2->transfer
 * array.
 * 
 * This function calls the other transfer2_XXX functions in the order needed to 
 * solve the line of sight integral at second order. In the process, it fills the
 * ptr2 structure so that it can be used in the subsequent modules
 * (spectra2.c, bispectra2.c ...).
 *
 * Before calling this function, make sure that background_init(),
 * thermodynamics_init() and perturb2_init() have already been executed.
 * 
 * Details on the physics and on the adopted method can be found in my thesis
 * (http://arxiv.org/abs/1405.2280) in chapter 5. The code itself is
 * extensively documented and hopefully will give you further insight.
 *
 * In detail, this function does:
 *
 * -# Based on the line of sight sources computed in the perturbations2.c module,
 *    determine which transfer functions need to be computed and their k3-sampling via
 *    transfer2_indices_of_transfers().
 *
 * -# Allocate one workspace per thread to store temporary values related
 *    to the line of sight integration.
 *
 * -# For each combination of k1 and k2, interpolate the line of sight sources at the 
 *    k3 values determined above, using transfer2_interpolate_sources_in_k().
 * 
 * -# For each combination of k1, k2 and k3, determine the best integration grid in
 *    time using transfer2_get_time_grid(), and interpolate the line of sight sources
 *    in such grid using transfer2_interpolate_sources_in_time().
 *
 * -# For each combination of k1, k2, k3 and type index, solve the line of sight
 *    integral using transfer2_compute(); this will fill the ptr2->transfer array.
 *
 * -# Free the workspaces and, if requested, store to disk the content of ptr2->transfer
 *    and free the array.
 */

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
      )
{

  // =================================================================================
  // =                              Preliminary checks                               =
  // =================================================================================


  /* Check whether any spectrum in harmonic space (i.e., any C_l's) is actually requested */

  if ((ppt2->has_perturbations2 == _FALSE_)
     || ((ppt2->has_cls == _FALSE_)
     && (ppt2->has_bispectra == _FALSE_)
     && (ptr2->has_transfers2_only == _FALSE_)) ) {

    ptr2->has_cls = _FALSE_;
    if (ptr2->transfer2_verbose > 0)
      printf("Second-order transfer module skipped.\n");

    return _SUCCESS_;

  }
  else
    ptr2->has_cls = _TRUE_;

  if (ptr2->transfer2_verbose > 0)
    printf("Computing second-order transfer functions\n");

  /* Get conformal age, recombination time and comoving sound horizon at recombination
  from background / thermodynamics structures (only place where these structures
  are used in this module) */
  ptr2->tau0 = pba->conformal_age;
  ptr2->tau_rec = pth->tau_rec;
  ptr2->rs_rec = pth->rs_rec;



  // ====================================================================================
  // =                             Indices and k-samplings                              =
  // ===================================================================================

  /* Initialize all indices in the transfers structure and allocate all its arrays. This
  function also calls transfer2_get_lm_lists() that fills ptr2->l, and transfer2_get_k3_sizes()
  that finds for all (k1,k2) pairs the allowed k3 range. */
  class_call (transfer2_indices_of_transfers (ppr,ppr2,ppt2,pbs,pbs2,ptr,ptr2),
    ptr2->error_message,
    ptr2->error_message);



  // ====================================================================================
  // =                         1st-order transfer function                              =
  // ====================================================================================

  /* Check that we have enough k1 values.  If, for example, we are testing the code
  with just one or two values of k1, it does not make sense to continue as the transfer.c
  module relies on the spline interpolation on the k1 grid. */
  if ((ppt2->k_size < 4) && (ptr2->has_transfers2_only == _FALSE_)) {
    printf("# WARNING: cannot do cubic interpolation in k1 with less than 4 values. ");
    printf("Increase k_size, or do not trust results for 1st-order transfer functions.\n");
  }


  /* Apart from ptr2->transfer, all the arrays needed by the subsequent modules have been filled.
  If the user requested to load the transfer functions from disk, we can stop the execution of
  this module now without regrets. */
  if (ppr2->load_transfers_from_disk == _TRUE_) {
    
    if (ptr2->transfer2_verbose > 0)
      printf(" -> leaving transfer2 module; transfer functions will be read from disk\n");
    
    return _SUCCESS_;
  }

  if (ptr2->transfer2_verbose > 0)
    printf(" -> starting actual computation of second-order transfer functions\n");


  // ==================================================================================
  // =                               Allocate workspaces                              =
  // ==================================================================================

  // *** Parallelization variables

  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;

  #ifdef _OPENMP
  #pragma omp parallel
  number_of_threads = omp_get_num_threads();
  #endif
  
  // *** Allocate a workspace per thread
    
  struct transfer2_workspace ** ppw;
  class_alloc(ppw, number_of_threads*sizeof(struct transfer2_workspace*), ptr2->error_message);
    
  /* We shall allocate the time arrays with the maximum possible number of integration steps. The
  actual steps depend on the wavemode k considered, with the largest k's having more time steps.
  This is because k acts as the frequency of the highly oscillatory projection functions J(k(tau0-tau)). */
  int tau_size_max;

  if (ptr2->tau_sampling == bessel_tau_sampling) {
    tau_size_max = ppt2->tau_size + pbs2->xx_size;
  }

  else if (ptr2->tau_sampling == custom_transfer2_tau_sampling) {
    double k_max = ptr2->k_max_k1k2[ppt2->k_size-1][ppt2->k_size-1];
    double tau_step_min = ppr2->tau_step_trans_song/k_max;
    double tau_max = ppt2->tau_sampling[ppt2->tau_size-1];
    tau_size_max = ppt2->tau_size + tau_max/tau_step_min + 1;
  }
  
  /* Print some information on the finest time grid that will be used */
  if (ptr2->transfer2_verbose > 0)
    printf (" -> maximum number of time steps in the LOS integration = %d\n", tau_size_max);
  
  /* Allocate arrays in the workspace */
  abort = _FALSE_;
  #pragma omp parallel shared(ppw,ptr2) private(thread)
  {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
    
    /* Allocate workspace array */
    class_alloc_parallel(
      ppw[thread],
      sizeof(struct transfer2_workspace),
      ptr2->error_message);

    /* Allocate the k-grid with the maximum possible number of k3-values.  */
    class_alloc_parallel(
      ppw[thread]->k_grid,
      ptr2->k3_size_max*sizeof(double),
      ptr2->error_message);

    /* Allocate the integration grid array. */
    class_alloc_parallel(
      ppw[thread]->tau_grid,
      tau_size_max*sizeof(double),
      ptr2->error_message);
    
    /* Allocate tau0_minus_tau */
    class_alloc_parallel(
      ppw[thread]->tau0_minus_tau,
      tau_size_max*sizeof(double),
      ptr2->error_message);
    
    /* Allocate delta_tau, which is needed for the trapezoidal approximation of the line-of-sight
      integration. This array is defined as tau(i+1)-tau(i-1) except for the first and last elements,
      which are, respectively, tau(1)-tau(0) and tau(N)-tau(N-1). */
    class_alloc_parallel(
      ppw[thread]->delta_tau,
      tau_size_max*sizeof(double),
      ptr2->error_message);
    
    /* Allocate index_tau_left, an array useful for the time interpolation of the sources (see header file) */
    class_alloc_parallel(
      ppw[thread]->index_tau_left,
      tau_size_max*sizeof(double),
      ptr2->error_message);

    /* Allocate the array that will contain the second derivatives of the sources with respect to time  */
    class_alloc_parallel (ppw[thread]->sources_time_spline, ppt2->tp2_size*sizeof(double *), ptr2->error_message);
      
    int index_tp;
    for (index_tp=0; index_tp<ppt2->tp2_size; ++index_tp)
      class_alloc_parallel(
        ppw[thread]->sources_time_spline[index_tp],
        tau_size_max*sizeof(double),
        ptr2->error_message);


    /* Allocate the array that will contain the sources interpolated at the points of the integration grid */
    class_alloc_parallel (ppw[thread]->interpolated_sources_in_time, ppt2->tp2_size*sizeof(double *), ptr2->error_message);

    for (index_tp=0; index_tp<ppt2->tp2_size; ++index_tp)
      class_alloc_parallel(
        ppw[thread]->interpolated_sources_in_time[index_tp],
        tau_size_max*sizeof(double),
        ptr2->error_message);

  } // end of parallel region
  
  if (abort == _TRUE_) return _FAILURE_;

  #ifdef _OPENMP
  if (ptr2->transfer2_verbose > 3)
    printf("In %s: Split the computation of the 2nd-order transfer functions between %d threads\n",
     __func__,number_of_threads);
  #endif

  /* Vector that will contain the second derivative of the sources for a given k1,k2 and source type */
  double ** sources_k_spline;
  class_alloc (sources_k_spline, ppt2->tp2_size*sizeof(double *), ptr2->error_message);
  
  /* Vector that will contain the interpolated sources for a given k1,k2 and source type */
  double ** interpolated_sources_in_k;
  class_alloc (interpolated_sources_in_k, ppt2->tp2_size*sizeof(double *), ptr2->error_message);
  

  // =====================================================================================
  // =                            Main loop on (l,m,k1,k2,k3)                            =
  // =====================================================================================

  /* Four loops over k1, k2, transfer type and k follow */
  for (int index_k1 = 0; index_k1 < ppt2->k_size; ++index_k1) {

    if (ptr2->transfer2_verbose > 1)
      printf ("     * computing transfer functions today for index_k1=%d of %d, k1=%g\n",
        index_k1, ppt2->k_size, ppt2->k[index_k1]);

    // ---------------------------------------------------------------------
    // -                       Load second-order sources                   -
    // ---------------------------------------------------------------------
    
    /* Allocate the k1 level of ptr2->transfer[index_type][index_k1][index_k2][index_k]
    the array used to store the results of the computations made in this module.  */
    class_call(transfer2_allocate_k1_level(ppt2, ptr2, index_k1),
      ptr2->error_message, ptr2->error_message);

    /* Load sources from disk if they were previously stored.  This can be true either because we are loading
    them from a precomputed run, or because we stored them in this run. */
    if ((ppr2->load_sources_from_disk == _TRUE_) || (ppr2->store_sources_to_disk == _TRUE_))
      class_call(perturb2_load_sources_from_disk(ppt2, index_k1),
          ppt2->error_message,
          ptr2->error_message);

    /* We only need to consider those k2's that are equal to or larger than k1,
    as the quadratic sources were symmetrised in the perturbation2 module */
    for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {

      /* Print some info */
      if (ptr2->transfer2_verbose > 2)
        printf(" -> computing transfer function for (k1,k2) = (%.3g,%.3g)\n", ppt2->k[index_k1], ppt2->k[index_k2]);

      // ---------------------------------------------------------------------
      // -                    Interpolate sources in k3                      -
      // ---------------------------------------------------------------------      

      /* Find the grid in k3 for the current (k1,k2) pair */
      int last_used_index_pt;
      double * k_grid_temp;
      class_alloc(k_grid_temp, ptr2->k3_size_max*sizeof(double), ptr2->error_message);

      class_call_parallel (transfer2_get_k3_list(
                    ppr,
                    ppr2,
                    ppt2,
                    pbs,
                    pbs2,
                    ptr2,
                    index_k1,
                    index_k2,
                    k_grid_temp,  /* output */
                    &last_used_index_pt
                    ),
        ptr2->error_message,
        ptr2->error_message);

      /* Copy the k3 grid in the pwb workspace. We do so in a parallel region because
      we have as many pwb[thread] workspaces as the number of threads */
      for (int thread=0; thread < number_of_threads; ++thread)
        for (int index_k=0; index_k < ptr2->k_size_k1k2[index_k1][index_k2]; ++index_k)
          ppw[thread]->k_grid[index_k] = k_grid_temp[index_k];
      
      /* Print some information */
      if (ptr2->transfer2_verbose > 3)
        printf("     * (k1,k2)=(%.3g,%.3g): the k3-grid comprises sources+transfer+left+right=%d+%d+%d+%d points from %g to %g\n",
          ppt2->k[index_k1], ppt2->k[index_k2],
          last_used_index_pt,
          ptr2->k_physical_size_k1k2[index_k1][index_k2] - last_used_index_pt,
          ptr2->k_physical_start_k1k2[index_k1][index_k2],
          ptr2->k_size_k1k2[index_k1][index_k2] - ptr2->k_physical_size_k1k2[index_k1][index_k2] - ptr2->k_physical_start_k1k2[index_k1][index_k2],
          ppw[thread]->k_grid[0], ppw[thread]->k_grid[ptr2->k_size_k1k2[index_k1][index_k2]-1]);

      free (k_grid_temp);
      

      // -----------------------------------------------------------------------
      // -                      Interpolate sources in k                       -
      // -----------------------------------------------------------------------

      for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
      
        class_alloc(
          sources_k_spline[index_tp],
          ppt2->k3_size[index_k1][index_k2]*ppt2->tau_size*sizeof(double),
          ptr2->error_message);
      
        class_alloc(
          interpolated_sources_in_k[index_tp],
          ptr2->k_size_k1k2[index_k1][index_k2]*ppt2->tau_size*sizeof(double),
          ptr2->error_message);
      
        class_call (transfer2_interpolate_sources_in_k(
                      ppr,
                      ppr2,
                      ppt,
                      ppt2,
                      pbs,
                      pbs2,
                      ptr2,
                      index_k1,
                      index_k2,
                      index_tp,
                      ppw[0]->k_grid,                       /* All the ppw->[thread]->k_grid are filled with the same values */
                      sources_k_spline[index_tp],           /* Will be filled with second-order derivatives */
                      interpolated_sources_in_k[index_tp]   /* Will be filled with interpolated values in ptr2->k(k1,k2) */
                      ),
          ptr2->error_message,
          ptr2->error_message);
      
      } // end of for (index_tp)

                    
      /* Beginning of parallel region */
      abort = _FALSE_;    
      #pragma omp parallel shared (ppw,ppr,ppt2,pbs,ptr2,abort) private (thread)
      {

        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        /* Update workspace */
        ppw[thread]->thread = thread;
        ppw[thread]->index_k1 = index_k1;
        ppw[thread]->k1 = ppt2->k[index_k1];
        ppw[thread]->index_k2 = index_k2;
        ppw[thread]->k2 = ppt2->k[index_k2];

        #pragma omp for schedule (static)
        for (int index_k = 0; index_k < ptr2->k_size_k1k2[index_k1][index_k2]; ++index_k) { 

          /* Update workspace */
          ppw[thread]->index_k = index_k;
          ppw[thread]->k = ppw[thread]->k_grid[index_k];


          // -----------------------------------------------------------------------
          // -                    Interpolate sources in time                      -
          // -----------------------------------------------------------------------

          /* Get the integration grid in time for the given k-mode */
          class_call_parallel (transfer2_get_time_grid(
                                 ppr,
                                 ppr2,
                                 ppt,
                                 ppt2,
                                 pbs,
                                 pbs2,
                                 ptr2,
                                 index_k1,
                                 index_k2,
                                 index_k,
                                 ppw[thread]
                                 ),
            ptr2->error_message,
            ptr2->error_message);

          /* Interpolate the sources at the right value of time */
            
          for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
      
            class_call_parallel (transfer2_interpolate_sources_in_time(
                          ppr,
                          ppr2,
                          ppt,
                          ppt2,
                          pbs,
                          pbs2,
                          ptr2,
                          index_tp,
                          interpolated_sources_in_k[index_tp],                 /* Must be already filled */
                          ppw[thread]->sources_time_spline[index_tp],          /* Will be filled with second-order derivatives */
                          ppw[thread]->interpolated_sources_in_time[index_tp], /* Will be filled with interpolated values in pw->tau_grid */
                          ppw[thread]
                          ),
              ptr2->error_message,
              ptr2->error_message);
      
          } // end of for (index_tp)

          
          for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {


            // ----------------------------------------------------------------------
            // -                    Compute transfer functions                      -
            // ----------------------------------------------------------------------

            /* Compute the transfer function for this (k1,k2,k). Note that index_tt is needed by the function
            for (i) find the type index corresponding to index_tt in ppt2->sources, (ii) store the result in
            ptr2->tranfer[index_tt] */
            class_call_parallel (transfer2_compute (
                                   ppr,
                                   ppr2,
                                   ppt2,
                                   pbs,
                                   pbs2,
                                   ptr2,
                                   index_k1,
                                   index_k2,
                                   index_k,
                                   ptr2->corresponding_index_l[index_tt],
                                   ptr2->corresponding_index_m[index_tt],
                                   index_tt,
                                   interpolated_sources_in_k,
                                   ppw[thread]
                                   ),
              ptr2->error_message,
              ptr2->error_message);

          } // end of for(index_tt)

          #pragma omp flush(abort)

        } // end of for(index_k) 

      } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */

      /* Free the memory for the interpolated sources */
      for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
        free(sources_k_spline[index_tp]);
        free(interpolated_sources_in_k[index_tp]);
      }

    } // end of for(index_k2)

    /* Free the memory associated with the line-of-sight sources for the considered k1.
    We won't need them anymore because the different k1 modes are independent. Note that
    this memory was either allocate at the beginning of the k1 loop, in this module, or
    at the beginning of the k1 loop in perturb2_init. */
    class_call (perturb2_free_k1_level (ppt2, index_k1),
      ppt2->error_message, ppt2->error_message);


    /* Save all transfer functions for the given k1, and free the memory associated with them.
    The next time we'll need them, we shall load them from disk.  */
    if (ppr2->store_transfers_to_disk == _TRUE_) {
      
      class_call (transfer2_store_transfers_to_disk (ppt2, ptr2, index_k1),
          ptr2->error_message,
          ptr2->error_message);

      class_call (transfer2_free_k1_level (ppt2, ptr2, index_k1),
        ptr2->error_message, ptr2->error_message);
    }

  } // end of for(index_k1)

  free (sources_k_spline);
  free (interpolated_sources_in_k);
  
  #pragma omp parallel shared(ppw) private(thread)
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    free(ppw[thread]->k_grid);
    free(ppw[thread]->tau_grid);
    free(ppw[thread]->tau0_minus_tau);
    free(ppw[thread]->delta_tau);
    free(ppw[thread]->index_tau_left);    
    int index_tp;
    for (index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
      free(ppw[thread]->sources_time_spline[index_tp]);
      free(ppw[thread]->interpolated_sources_in_time[index_tp]);
    }
    free(ppw[thread]->sources_time_spline);
    free(ppw[thread]->interpolated_sources_in_time);
    free(ppw[thread]);
  }
    
  free(ppw);
  
  /* We are finished filling the transfer function files, so close them */
  if (ppr2->store_transfers_to_disk == _TRUE_)
    for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++)
      fclose (ptr2->transfers_files[index_tt]);

  if (ptr2->transfer2_verbose > 1)
    printf (" -> filled ptr2->transfer with %ld values (%g MB)\n",
      ptr2->count_memorised_transfers, ptr2->count_memorised_transfers/1e6*8);
  
  /* Check that the number of filled values corresponds to the number of allocated space */
  if (ppr2->load_transfers_from_disk == _FALSE_)
    class_test (ptr2->count_allocated_transfers != ptr2->count_memorised_transfers,
      ptr2->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!", ptr2->count_allocated_transfers, ptr2->count_memorised_transfers);
  
  /* Do not evaluate the subsequent modules if ppt2->has_early_transfers2_only == _TRUE_ */
  if (ptr2->has_transfers2_only == _TRUE_) {
    ppt->has_cls = _FALSE_;
    ppt2->has_cls = _FALSE_;
    ppt2->has_bispectra = _FALSE_;
  }
  
  return _SUCCESS_;
}


/**
 * Allocate the type level of the transfer functions array.
 *
 * Allocate the transfer type (tt) level of the array that will contain the second-order
 * transfer functions: ptr2->transfer[index_tt][index_k1][index_k2][index_k].
 *
 * This function makes space for the transfer functions; it is called before loading
 * them from disk in transfer2_load_transfers_from_disk(). It relies on values that are
 * computed by transfer2_indices_of_perturbs() and transfer2_get_k3_sizes(), so make
 * sure to call them beforehand.
 *
 */
int transfer2_allocate_type_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_tt
     )
{

  long int count=0;
  int k1_size = ppt2->k_size;

  class_alloc(
    ptr2->transfer[index_tt],
    k1_size * sizeof(double **),
    ptr2->error_message);

  for (int index_k1=0; index_k1<k1_size; ++index_k1) {

    /* Allocate k2 level.  Note that, as for ppt2->sources, the size of this level is smaller
    than ppt2->k_size, and it depends on k1.  The reason is that we only need to compute
    the transfer functions for those k2's that are smaller than k1 (our equations are symmetrised
    wrt to k1<->k2) */
    int k2_size = index_k1 + 1;
  
    class_alloc(
      ptr2->transfer[index_tt][index_k1],
      k2_size * sizeof(double *),
      ptr2->error_message);
  
    for (int index_k2=0; index_k2<=index_k1; ++index_k2) {

      /* Allocate k level. Note that we are using ptr2->k_size here instead of ppt2->k_size.  The
      reason is that ptr2->k is sampled much more finely than ppt2->k in order to catch the wild
      oscillation of the Bessel functions in k.  Furthermore, we only allocate memory for the k's
      that are compatible with the values of k1 and k2, i.e. we impose fabs(cosk1k2) <= 1.  See
      transfer2_get_k3_sizes() for further details. */
      class_alloc(
        ptr2->transfer[index_tt][index_k1][index_k2],
        ptr2->k_size_k1k2[index_k1][index_k2] * sizeof(double),
        ptr2->error_message);

      #pragma omp atomic
      ptr2->count_allocated_transfers += ptr2->k_size_k1k2[index_k1][index_k2];
      count += ptr2->k_size_k1k2[index_k1][index_k2];

    } // end of for(index_k2)
  } // end of for(index_k1)
  
  /* Print some debug information on memory consumption */
  if (ptr2->transfer2_verbose > 2) {
    printf("     * allocated ~ %.2f MB in ptr2->transfer (%ld doubles) for index_tt=%d; it's size is now ~ %.3g MB;\n",
      count*sizeof(double)/1e6, count, index_tt, ptr2->count_allocated_transfers*sizeof(double)/1e6);
  }

  return _SUCCESS_;

}


/**
 * Load the transfer functions from disk for a given transfer type.
 *
 * The transfer functions will be read from the file given in ptr2->transfers_paths[index_type]
 * and stored in the array ptr2->transfer. Before running this function, make sure to
 * allocate the corresponding type level of ptr2->transfer using transfer2_allocate_type_level().
 *
 * This function is used in the bispectra2.c module.
 */
int transfer2_load_transfers_from_disk(
        struct perturbs2 * ppt2,
        struct transfers2 * ptr2,
        int index_tt
        )
{

  /* Running indexes */
  int index_k1, index_k2;

  /* Allocate memory to keep the transfer functions */
  class_call (transfer2_allocate_type_level(ppt2, ptr2, index_tt),
    ptr2->error_message, ptr2->error_message);

  /* Print some debug */
  if (ptr2->transfer2_verbose > 2)
    printf("     * transfer2_load_transfers_from_disk: reading results for index_tt=%d from '%s' ...",
      index_tt, ptr2->transfers_paths[index_tt]);
  
  /* Open file for reading */
  class_open (ptr2->transfers_files[index_tt], ptr2->transfers_paths[index_tt], "rb", ptr2->error_message);
  
  /* Two loops follow to read the file */
  for (index_k1 = 0; index_k1 < ppt2->k_size; ++index_k1) {
  
    for (index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
  
      int n_to_read = ptr2->k_size_k1k2[index_k1][index_k2];
  
      /* Some debug */
      // printf ("reading %d entries for index_tt=%d, index_k1=%d, index_k2=%d\n", n_to_read, index_tt, index_k1, index_k2);
  
      /* Read a chunk with all the k-values for this set of (type,k1,k2) */
      int n_read = fread(
              ptr2->transfer[index_tt][index_k1][index_k2],
              sizeof(double),
              n_to_read,
              ptr2->transfers_files[index_tt]);
  
      class_test(n_read != n_to_read,
        ptr2->error_message,
        "Could not read in '%s' file, read %d entries but expected %d (index_tt=%d,index_k1=%d,index_k2=%d)",
          ptr2->transfers_paths[index_tt], n_read, n_to_read, index_tt, index_k1, index_k2);
  
      /* Update the counter for the values stored in ptr2->transfers */
      #pragma omp atomic
      ptr2->count_memorised_transfers += ptr2->k_size_k1k2[index_k1][index_k2];
  
    } // end of for(index_k2)
    
  } // end of for(index_k1)
  
  /* Close file */
  fclose(ptr2->transfers_files[index_tt]);

  if (ptr2->transfer2_verbose > 2)
    printf ("Done.\n");

  return _SUCCESS_; 
  
}


/**
 * This routine frees all the memory space allocated by transfer2_init().
 *
 * To be called at the end of each run, only when no further calls to
 * transfer_functions_at_k() are needed.
 */ 
int transfer2_free(
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct transfers2 * ptr2
      )
{

  if (ptr2->has_cls == _TRUE_) {
  
    int k1_size = ppt2->k_size;

    /* Free the k1 level of ptr2->transfer only if we are neither loading nor storing the transfers to disk.
    The memory management in those two cases is handled separately via the transfer_store and transfer_load
    functions  */
    if ((ppr2->store_transfers_to_disk==_FALSE_) && (ppr2->load_transfers_from_disk==_FALSE_)) {
      for (int index_k1 = 0; index_k1 < k1_size; ++index_k1)
        if (ptr2->has_allocated_transfers[index_k1] == _TRUE_)
          class_call(transfer2_free_k1_level(ppt2, ptr2, index_k1), ptr2->error_message, ptr2->error_message);
      for (int index_tt = 0; index_tt < ptr2->tt2_size; ++index_tt)
        free (ptr2->transfer[index_tt]);
    }

    free (ptr2->transfer);

    free (ptr2->has_allocated_transfers);

    for(int index_k1=0; index_k1<ppt2->k_size; ++index_k1) {
      free (ptr2->k_size_k1k2[index_k1]);
      free (ptr2->k_physical_start_k1k2[index_k1]);
      free (ptr2->k_physical_size_k1k2[index_k1]);
      free (ptr2->k_min_k1k2[index_k1]);
      free (ptr2->k_max_k1k2[index_k1]);
    }
    free (ptr2->k_size_k1k2);
    free (ptr2->k_physical_start_k1k2);
    free (ptr2->k_physical_size_k1k2);
    free (ptr2->k_min_k1k2);
    free (ptr2->k_max_k1k2);
      
    /* Free labels */
    for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++)
      free(ptr2->tt2_labels[index_tt]);

    free (ptr2->l);
    free (ptr2->m);
    free (ptr2->tt2_labels);

    free (ptr2->corresponding_index_l);
    free (ptr2->corresponding_index_m);
    free (ptr2->index_tt2_monopole);
    free (ptr2->index_pt2_monopole);
  
    for (int index_l=0; index_l<ptr2->l_size; ++index_l)
      free (ptr2->lm_array[index_l]);
    free (ptr2->lm_array);

    /* Free file arrays */
    if ((ppr2->store_transfers_to_disk == _TRUE_) || (ppr2->load_transfers_from_disk == _TRUE_)) {

      // fclose(ptr2->transfers_status_file);

      for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt)
        free (ptr2->transfers_paths[index_tt]);

      free (ptr2->transfers_files);
      free (ptr2->transfers_paths);
    }

  } // end of if(has_cls)

  return _SUCCESS_;
  
}


/**
 *
 * Initialize indices and arrays in the second-order transfer functions structure.
 *
 * In detail, this function does:
 *
 *  -# Computes the list of multipoles where to compute the transfer functions via
 *     transfer2_get_l_list().    
 *
 *  -# Determine which transfer functions to compute based on the previous modules, and
 *     assign them the ptr2->index_tt2_XXX indices.
 *
 *  -# Define the indexing strategy of the transfer functions array (ptr2->transfer)
 *     and establish a correspondence between line of sight sources and transfer 
 *     functions, via transfer2_get_lm_lists().
 *
 *  -# In transfer2_get_k3_sizes(), fill the array ptr2->k_size_k1k2 which, for each
 *     (k1,k2) pair, determines the size of the k3 grid of the transfer functions, based
 *     on the triangular condition.
 *
 *  -# Allocate the type level of the transfer array, ptr2->transfer[index_tt].
 *
 *  -# Open the files where we will store the transfer functions at the end of the computation.
 *
 */
int transfer2_indices_of_transfers(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct perturbs2 * ppt2,
          struct bessels * pbs,
          struct bessels2 * pbs2,
          struct transfers * ptr,
          struct transfers2 * ptr2
          )
{

  /* Get l values using transfer2_get_l_list() (same for all modes) */
  class_call (transfer2_get_l_list(ppr,ppr2,ppt2,pbs,pbs2,ptr2),
    ptr2->error_message,
    ptr2->error_message);
  
  /* Copy the m array from ppr->m */
  ptr2->m_size = ppr2->m_size;
  class_alloc (ptr2->m, ptr2->m_size*sizeof(int), ptr2->error_message);
  for (int index_m=0; index_m<ptr2->m_size; ++index_m)
    ptr2->m[index_m] = ppr2->m[index_m];


  // ======================================================================================
  // =                            What transfers to compute?                              =
  // ======================================================================================

  /* Index that will keep track of the trasfer-function type to compute */
  int index_tt = 0;

  /* Number of transfer functions that need to be computed for each source type (temperature,
  polarisation...). This is given by all the possible (l,m) combinations contained in ptr2->l
  and ppr2->m */
  ptr2->n_transfers = size_indexl_indexm (ptr2->l, ptr2->l_size_max, ppr2->m, ppr2->m_size);
  
  /* Number of non-vanishing E-mode transfer functions, for debug purposes */
  ptr2->n_nonzero_transfers_E = ptr2->n_transfers;
  if (ppr2->compute_m[0])
    ptr2->n_nonzero_transfers_E -= 2;
  if (ppr2->compute_m[1])
    ptr2->n_nonzero_transfers_E -= 1;

  /* Number of non-vanishing B-mode transfer functions, for debug purposes */
  ptr2->n_nonzero_transfers_B = ptr2->n_transfers;
  if (ppr2->compute_m[0])
    ptr2->n_nonzero_transfers_B -= ptr2->l_size_max;
  if (ppr2->compute_m[1])
    ptr2->n_nonzero_transfers_B -= 1;


  /* Photon temperature transfer functions */
  if (ppt2->has_cmb_temperature == _TRUE_) {
  
    ptr2->index_tt2_T = index_tt;
    index_tt += ptr2->n_transfers; 
  }

  /* Photon E-mode polarisation transfer functions */
  if (ppt2->has_cmb_polarization_e == _TRUE_) {
  
    ptr2->index_tt2_E = index_tt;
    index_tt += ptr2->n_transfers;    
  }

  /* Photon B-mode polarisation transfer functions. To be computed
  only for non scalar modes, otherwise they just vanish. */
  if (ppt2->has_cmb_polarization_b == _TRUE_) {
  
    if (ppr2->m_max_2nd_order>0) {

      ptr2->index_tt2_B = index_tt;
      index_tt += ptr2->n_transfers;
    }
        
  }

  /* Total number of transfer functions to compute */
  ptr2->tt2_size = index_tt;
 

  /* Allocate memory for the labels of the transfer types */
  class_alloc(ptr2->tt2_labels, ptr2->tt2_size*sizeof(char *), ptr2->error_message);
  for (int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt)
    class_alloc(ptr2->tt2_labels[index_tt], 64*sizeof(char), ptr2->error_message);


  if (ptr2->transfer2_verbose > 1) {
    printf (" -> will compute tt2_size=%d transfer functions: ", ptr2->tt2_size);
    if (ppt2->has_cmb_temperature == _TRUE_) printf ("T=%d ", ptr2->n_transfers);
    if (ppt2->has_cmb_polarization_e == _TRUE_) printf ("E=%d (%d non-zero) ",
      ptr2->n_transfers, ptr2->n_nonzero_transfers_E);
    if (ppt2->has_cmb_polarization_b == _TRUE_) printf ("B=%d (%d non-zero) ",
      ptr2->n_transfers, ptr2->n_nonzero_transfers_B);
    printf ("\n");
  }  

  
  // ==============================================================================
  // =                            Fill (l,m) arrays                               =
  // ==============================================================================

  class_call (transfer2_get_lm_lists (
                ppr,
                ppr2,
                ppt2,
                pbs,
                pbs2,
                ptr2),
    ptr2->error_message,
    ptr2->error_message);



  // ==================================================================================
  // =                       Determine range of k3(k1,k2)                             =
  // ==================================================================================

  class_call (transfer2_get_k3_sizes (
                ppr,
                ppr2,
                ppt2,
                pbs,
                pbs2,
                ptr,
                ptr2),
    ptr2->error_message,
    ptr2->error_message);



  // =======================================================================================
  // =                      Allocate first levels of ptr2->transfer                        =
  // =======================================================================================
  
  /* ptr2->transfer has five levels that should be indexed in the following way:
        ptr2->transfer [index_tt]
                       [index_k1]
                       [index_k2]
                       [index_k]
    Important: as in ppt2->sources, index_k2 goes from 0 to ppt2->k_size-index_k1.  */

  /* Counter to keep track of the memory usage of ptr2->transfer */
  ptr2->count_memorised_transfers = 0;
  ptr2->count_allocated_transfers = 0;
    
  /* Allocate transfer-type (tt2) level */  
  class_alloc (
        ptr2->transfer,
        ptr2->tt2_size * sizeof(double ***),
        ptr2->error_message);


  /* Allocate level that will address k1.  The ptr2->transfer array has also a k and a k2
  levels, which we do not allocate here.  They will be allocated later in transfer2_init,
  but only if the user did not ask to load the transfer functions from a run directory */
  for (int index_tt = 0; index_tt < ptr2->tt2_size; ++index_tt) {
  
    int k1_size = ppt2->k_size;
    
    class_alloc (
      ptr2->transfer[index_tt],
      k1_size * sizeof(double **),
      ptr2->error_message);
  
  } // end of for(index_type)
  

  /* Allocate and initialize the logical array that keeps track of the memory state of ptr2->transfers */
  class_calloc(ptr2->has_allocated_transfers, ppt2->k_size, sizeof(short), ptr2->error_message);


  // ==================================================================================
  // =                             Create transfers files                             =
  // ==================================================================================
  

  /* Create the files to store the transfer functions in */
  if ((ppr2->store_transfers_to_disk == _TRUE_) || (ppr2->load_transfers_from_disk == _TRUE_)) {

    /* We are going to store the transfers in n=k_size files, one for each requested k1 */
    class_alloc (ptr2->transfers_files, ptr2->tt2_size*sizeof(FILE *), ptr2->error_message);
    class_alloc (ptr2->transfers_paths, ptr2->tt2_size*sizeof(char *), ptr2->error_message);

    for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt) {
      /* The name of each transfers file will have the tt index in it */
      class_alloc (ptr2->transfers_paths[index_tt], _FILENAMESIZE_*sizeof(char), ptr2->error_message);
      sprintf (ptr2->transfers_paths[index_tt], "%s/transfers_%03d.dat", ptr2->transfers_dir, index_tt);
      if (ppr2->store_transfers_to_disk == _TRUE_)
        class_open (ptr2->transfers_files[index_tt], ptr2->transfers_paths[index_tt], "wb", ptr2->error_message);
      
    } // end of loop on index_tt

    if (ptr2->transfer2_verbose > 2)
      printf ("     * created %d files to store transfer functions\n", ptr2->tt2_size);

  } // end of if(ppr2->store_transfers_to_disk)
  

  return _SUCCESS_;

}










/**
 * This routine defines the number and values of mutipoles l for all modes.
 * It basically copies the pbs->l into ptr->l.  Note that pbs->l was computed
 * in the Bessel module in 'bessel_get_l_list', ans it is just a log + linear
 * sampling.
 *
 * @param ppr  Input : pointer to precision structure
 * @param ppt  Input : pointer to perturbation structure
 * @param pbs  Input : pointer to bessels structure
 * @param ptr  Input/Output : pointer to transfers structure containing l's
 * @return the error status
 */

int transfer2_get_l_list (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers2 * ptr2
      )
{


  /* Just consider all l's in pbs->l regardless of the mode */
  ptr2->l_size = pbs->l_size; 
  ptr2->l_size_max = pbs->l_size;
  
  /* Copy the l's from pbs->l into ptr2->l up to ptr2->l_size_max (original CLASS) */
  class_alloc(ptr2->l, ptr2->l_size_max*sizeof(int), ptr2->error_message);

  for (int index_l=0; index_l < ptr2->l_size_max; index_l++) {
    ptr2->l[index_l] = pbs->l[index_l];

    // *** Some debug
    // printf("ptr2->l[%d] = %d\n", index_l, ptr2->l[index_l]);
  }
  
  return _SUCCESS_;

}








/**
 * Allocate and fill all the arrays that deal with the (l,m) indexing of the hierarchies
 * in the ptr2->transfer array.
 */
int transfer2_get_lm_lists (
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2
        )
{ 

  // =========================================================================================
  // =                                (l,m) indexing of arrays                               =
  // =========================================================================================

  /* Fill ptr2->lm_array, that contains the index associated with a given (l,m) couple.
  lm_array is useful every time you need a specific (l,m) value from the transfer
  function array: ptr2->transfer[index_type + lm_array[index_l][index_m].
  Note that this array is used only by the preprocessor macros lm_cls(index_l,index_m),
  defined in common2.h. */
  class_alloc(ptr2->lm_array, ptr2->l_size*sizeof(int*), ptr2->error_message);    
  
  for (int index_l=0; index_l<ptr2->l_size; ++index_l) {
  
    int l = ptr2->l[index_l];
  
    class_calloc (ptr2->lm_array[index_l],
                  ppr2->index_m_max[l]+1,
                  sizeof(int),
                  ptr2->error_message);

    /* We enter the loop only if index_m_max >= 0. This is false if l<ptr2->m[0],
    for example for l=2 and m_min=3 */
    for (int index_m=0; index_m <= ppr2->index_m_max[l]; ++index_m) {

      ptr2->lm_array[index_l][index_m] = multipole2offset_indexl_indexm (
                                           ptr2->l[index_l],
                                           ptr2->m[index_m],
                                           ptr2->l, ptr2->l_size,
                                           ptr2->m, ptr2->m_size
                                           );
                                                           
      /* Debug lm_array */
      // printf ("(l,m) = (%d,%d) corresponds to an offset of %d\n",
      //   ptr2->l[index_l], ptr2->m[index_m], ptr2->lm_array[index_l][index_m]);

    } // end of for(index_m)
  } // end of for(index_l)



  // ======================================================================================
  // =                          Transfers - Sources correspondence                        =
  // ======================================================================================

  /* In the transfer2 module a transfer function type (index_tt) encloses both the
  multipole (l,m) and the source type (temperature, polarization, etc). Here, we
  connect index_tt to the corresponding source types index_tp and multipoles
  by defining ad hoc arrays */
  
  /* For a given transfer type index index_tt, these arrays contain the corresponding
  indices in ptr2->l and ptr2->m */
  class_alloc (ptr2->corresponding_index_l, ptr2->tt2_size*sizeof(int), ptr2->error_message);
  class_alloc (ptr2->corresponding_index_m, ptr2->tt2_size*sizeof(int), ptr2->error_message);

  /* For a given transfer type index index_tt, these arrays contain the indices corresponding
  to the monopole of index_tt in ptr2->transfers and ppt2->sources, respectively. */
  class_alloc (ptr2->index_tt2_monopole, ptr2->tt2_size*sizeof(int), ptr2->error_message);
  class_alloc (ptr2->index_pt2_monopole, ptr2->tt2_size*sizeof(int), ptr2->error_message);

  for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    
    /* Initialise the vectors to -1 */
    int index_l = -1;
    int index_m = -1;

    // *** Photon temperature ***
    if ((ppt2->has_source_T==_TRUE_)
    && (index_tt >= ptr2->index_tt2_T) && (index_tt < ptr2->index_tt2_T+ptr2->n_transfers)) {

      /* Find the position of the monopole of the same type as index_tp */
      ptr2->index_tt2_monopole[index_tt] = ptr2->index_tt2_T;
      ptr2->index_pt2_monopole[index_tt] = ppt2->index_tp2_T;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->index_tt2_monopole[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "T_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("T, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->index_tt2_monopole[index_tt], ptr2->index_pt2_monopole[index_tt]);

    }

    // *** Photon E-mode polarization ***
    else if ((ppt2->has_source_E==_TRUE_)
    && (index_tt >= ptr2->index_tt2_E) && (index_tt < ptr2->index_tt2_E+ptr2->n_transfers)) {

      /* Find the position of the monopole of the same type as index_tp */
      ptr2->index_tt2_monopole[index_tt] = ptr2->index_tt2_E;
      ptr2->index_pt2_monopole[index_tt] = ppt2->index_tp2_E;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->index_tt2_monopole[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "E_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("E, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->index_tt2_monopole[index_tt], ptr2->index_pt2_monopole[index_tt]);

    }

    // *** Photon B-mode polarization ***
    else if ((ppt2->has_source_B==_TRUE_)
    && (index_tt >= ptr2->index_tt2_B) && (index_tt < ptr2->index_tt2_B+ptr2->n_transfers)) {

      /* Find the position of the monopole of the same type as index_tp */
      ptr2->index_tt2_monopole[index_tt] = ptr2->index_tt2_B;
      ptr2->index_pt2_monopole[index_tt] = ppt2->index_tp2_B;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->index_tt2_monopole[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "B_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("B, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->index_tt2_monopole[index_tt], ptr2->index_pt2_monopole[index_tt]);

    }

    /* Check the result */
    class_test ((index_l>=ptr2->l_size) || (index_m>=ptr2->m_size) || (index_l<0) || (index_m<0),
      ptr2->error_message,
      "index_tt=%d: result (index_l,index_m)=(%d,%d) is out of bounds index_l=[%d,%d], index_m=[%d,%d]\n",
      index_tt, index_l, index_m, 0, ptr2->l_size-1, 0, ptr2->m_size-1);

    /* Write the result */
    ptr2->corresponding_index_l[index_tt] = index_l;
    ptr2->corresponding_index_m[index_tt] = index_m;
    
  } // end of for(index_tt)


  return _SUCCESS_;

}
    







/**
 * For a given (k1,k2) pair, find the size of the corresponding k3 grid.
 *
 * At second order, the transfer functions are sampled in a 3D Fourier space (k1,k2,k3).
 * The k3 direction is special because in the line of sight integral (see eq. 5.95 of
 * http://arxiv.org/abs/1405.2280) it is convolved with a high frequency projection function
 * J_Llm(k3*tau). Therefore, the k3 direction requires a special sampling that can catch the
 * oscillations of J_Llm. Here we find the size of such sampling for a given (k1,k2) pair and
 * store it in the array ptr2->k_size_k1k2[index_k1][index_k2].
 * 
 * Extrapolation: to solve the bispectrum integral (eq. 6.36) it is useful to extrapolate
 * the transfer functions in the k3 direction beyond their physical limits, that is, beyond
 * k3_min = |k1-k2| and k3_max = k1+k2. The reason is purely numerical: eventually, the
 * contributions from the extrapolated regions will cancel. The extrapolation has the purpose
 * of stabilizing an problematic integration.
 *
 * The following arrays are filled in this function:
 *  - ptr2->k_min_k1k2[index_k1][index_k2]
 *  - ptr2->k_max_k1k2[index_k1][index_k2]
 *  - ptr2->k_size_k1k2[index_k1][index_k2]
 *  - ptr2->k_physical_start_k1k2[index_k1][index_k2]
 *  - ptr2->k_physical_size_k1k2[index_k1][index_k2]
 */
int transfer2_get_k3_size (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers2 * ptr2,
      int index_k1,
      int index_k2
      )
{

  class_test (index_k2 > index_k1,
              ptr2->error_message,
              "stop to avoid segmentation fault, as we always require k1 >= k2.");

  /* Limits of k3 as determined by the triangular condition \vec{k3} = \ve{k1} + \vec{k2} */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  double k_min_pt = ppt2->k3[index_k1][index_k2][0];
  double k_max_pt = ppt2->k3[index_k1][index_k2][k_pt_size-1];

  /* Some debug on the limits of k3 */
  // if ((index_k1==85) && (index_k2==63)) {
  //   
  //   printf ("(k_min_pt,|k1-k2|) = (%g,%g)\n", k_min_pt, fabs(ppt2->k[index_k1]-ppt2->k[index_k2]));
  //   printf ("(k_max_pt, k1+k2 ) = (%g,%g)\n", k_max_pt, ppt2->k[index_k1]+ppt2->k[index_k2]);
  //   
  // }

  /* By default, for the transfer functions we take the same k-limits used for the sources */
  double k_min_tr = k_min_pt;
  double k_max_tr = k_max_pt;

  /* Extend the k-grid beyond the physical limits to make the bispectrum integral numerically stable. */
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {

    /* The physical range is the one dictated by the triangular condition: k1 + k2 = k3 */
    double physical_k3_range = k_max_pt - k_min_pt;
    
    /* The extended range allows for the Bessel function j(k3*(tau0 - tau_rec)) to develop at least
    ppr->extra_k3_oscillations oscillations in k3 in order to stabilize the bispectrum integral. */
    double extended_k3_range_right = 2.*_PI_*ppr->extra_k3_oscillations_right/(ptr2->tau0 - ptr2->tau_rec);
    double extended_k3_range_left  = 2.*_PI_*ppr->extra_k3_oscillations_left/(ptr2->tau0 - ptr2->tau_rec);

    if (physical_k3_range < extended_k3_range_right)
      k_max_tr += (extended_k3_range_right - physical_k3_range);    

    /* Never go beyond the maximum allowed value for k1 and k2 */
    k_max_tr = MIN (k_max_tr, ppt2->k[ppt2->k_size-1]);

    if (physical_k3_range < extended_k3_range_left)
      k_min_tr -= (extended_k3_range_left - physical_k3_range);

    /* Never go below the smallest between the usual k-range and the physical bound. */
    k_min_tr = MAX (k_min_tr, MIN(ppt2->k[0], k_min_pt));

    /* If the user specifies a negative number of oscillations, then extend the integration
    range all the way to the minimum and maximum k's considered in the code */    
    if (ppr->extra_k3_oscillations_left < 0)
      k_min_tr = ppt2->k[0];

    if (ppr->extra_k3_oscillations_right < 0)
      k_max_tr = ppt2->k[ppt2->k_size-1];

    /* Some debug */
    // printf("physical_k3_range=%g, extended_k3_range_right=%g\n", physical_k3_range, extended_k3_range_right);
    // printf("PRE:  k_max_tr = %g\n", k_max_pt);
    // printf("POST: k_max_tr = %g\n", k_max_tr);

  }

  /* Update the transfer structure with the k3 limits */
  ptr2->k_min_k1k2[index_k1][index_k2] = k_min_tr;
  ptr2->k_max_k1k2[index_k1][index_k2] = k_max_tr;
  
  /* Some debug */
  // printf ("k_max_tr(%d,%d) = %g\n", index_k1, index_k2, k_max_tr);
  

  /* Check that we the x-grid of the projection functions J is large enough. We do not need to test the lower
  limit as our x-grid always starts from zero. */
  class_test (k_max_tr > pbs2->xx_max/ptr2->tau0,
    ptr2->error_message,
    "xx_max is not large enough. The projection functions J(k*(tau0-tau)) do not cover the needed k-range \
(k_max=%g,tau0=%g,xx_max=%g). Maybe you are using too much extrapolation.", k_max_tr, ptr2->tau0, pbs2->xx_max);

  /* Determine type of k-sampling. */
  /* TODO: For later times in the line-of-sight integration, where the frequency tau0-tau is small,
  one can think of using a less dense k-grid */
  double k_step_max;

  if (ptr2->k_sampling == bessel_k3_sampling) {
    /* We require the maximum step between two points to be determined by the sampling of the projection
    functions J. The argument of J is x=k*(tau0-tau) and is linearly sampled in x. Hence, the maximum
    step in k is given by the x/tau0. */
    k_step_max = pbs2->xx_step/ptr2->tau0;
  }
  else if (ptr2->k_sampling == class_transfer2_k3_sampling) {

    k_step_max = 2.*_PI_/(ptr2->tau0-ptr2->tau_rec)*ppr2->q_linstep_song;
    
    /* Old style sampling */
    if ((ppr->load_run == _TRUE_) && (ppr2->old_run == _TRUE_))
      k_step_max = 2.*_PI_/ptr2->rs_rec*ppr2->q_linstep_song;
  }

  // *** Count the number of necessary values

  /* First point */
  
  int index_k_pt = 0;
  int index_k_tr = 0;  
  double k = k_min_tr;
  index_k_tr++;

  /* TODO: Magic number */
  double n = 1e99;

  /* Add unphysical points to the left of k_min_tr in order to stabilize the bispectrum integral */
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {

    /* Count points outside the physical regime, i.e. k3 values smaller than k_min_pt. Only matters when
    k_min_tr != k_min_pt, which happens only if extrapolation is turned on */
    double first_physical_step = ppt2->k3[index_k1][index_k2][1] - ppt2->k3[index_k1][index_k2][0];

    class_test (first_physical_step < ppr->smallest_allowed_variation,
      ptr2->error_message,
      "stopping to avoid segmentation fault, step=%g", first_physical_step);

    /* Take a linear step in the extrapolation regime */
    while (k < k_min_pt) {
      k += MIN (k_step_max, n*first_physical_step);
      index_k_tr++;
    }

    if (k != k_min_pt)
      k = k_min_pt;
  }
  
  index_k_pt++;

  /* The regime where the triangular condition is satisfied starts here */
  ptr2->k_physical_start_k1k2[index_k1][index_k2] = index_k_tr-1;

  /* Points taken from perturbation module if the step is small enough */
  while ((index_k_pt < k_pt_size) && ((ppt2->k3[index_k1][index_k2][index_k_pt] - k) < k_step_max)) {
    k = ppt2->k3[index_k1][index_k2][index_k_pt];
    index_k_pt++;
    index_k_tr++;
  }

  /* All the points that we added in the above while-loop satisfy the triangular condition */
  ptr2->k_physical_size_k1k2[index_k1][index_k2] = index_k_pt;

  /* Make sure that we do not introduce a jump in the extrapolation */
  double last_triangular_step = ppt2->k3[index_k1][index_k2][k_pt_size-1] - ppt2->k3[index_k1][index_k2][k_pt_size-2];

  class_test(k_step_max == 0.,
    ptr2->error_message,
    "stop to avoid infinite loop");

  /* Then, points spaced linearily with step k_step_max. Note that when using k3-extrapolation, we
  add points beyond the physical limit dictated by the triangular condition. */
  while (k < k_max_tr) {

    k += MIN (k_step_max, n*last_triangular_step);
    index_k_tr++;
    
    /* Take note of where the physical regime ends */
    if ((k <= k_max_pt) && (k >= k_min_pt))
      ptr2->k_physical_size_k1k2[index_k1][index_k2]++;
    
  }
  
  /* Get number of points and allocate list */

  if (k > k_max_tr)
    ptr2->k_size_k1k2[index_k1][index_k2] = index_k_tr-1;

  else
    ptr2->k_size_k1k2[index_k1][index_k2] = index_k_tr;

  
  class_test (ptr2->k_size_k1k2[index_k1][index_k2] < 2,
    ptr2->error_message,
    "found less than 2 valid values for the k-grid of modes k1=%f, k2=%f. This should not happen \
unless you set k3_size_min to a value lower than 2, which is not allowed as you need at least 2 points\
to do linear interpolation.\n",
    ppt2->k[index_k1], ppt2->k[index_k2])
  
  class_test ((ptr2->k_physical_size_k1k2[index_k1][index_k2] != ptr2->k_size_k1k2[index_k1][index_k2])
    && (ppr->bispectra_k3_extrapolation == no_k3_extrapolation),
    ptr2->error_message,
    "when no extrapolation is requested, all k's should be physical, i.e. satisfy the triangular condition.");
  
  class_test ((ptr2->k_physical_start_k1k2[index_k1][index_k2] != 0)
    && (ppr->bispectra_k3_extrapolation == no_k3_extrapolation),
    ptr2->error_message,
    "when no extrapolation is requested, all k's should be physical, i.e. satisfy the triangular condition.");
        
  return _SUCCESS_;
  
  
} // end of transfer2_get_k3_size




/**
 * For a given (k1,k2) pair, find the size of the corresponding k3 grid.
 * 
 * This function relies on transfer2_get_k3_sizes() being called beforehand to compute
 * the array ptr2->k_size_k1k2[index_k1][index_k2].
 *
 * Please refer to the documentation of transfer2_get_k3_size() for further detail.
 */
int transfer2_get_k3_list (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers2 * ptr2,
      int index_k1,
      int index_k2,
      double * k3, /**< output, array of size ptr2->k_size_k1k2[index_k1][index_k2] with
                      the k3 sampling of the transfer function in (k1,k2) */
      int * last_used_index_pt
      )
{

  /* MY MODIFICATIONS */
  // if (ppt2->k3_size[index_k1][index_k2] > 0)
  //   /* Implement the case when k3_size < 0 */
  
  /* Limits of k3 as determined by the triangular condition \vec{k3} = \ve{k1} + \vec{k2} */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];  
  double k_min_pt = ppt2->k3[index_k1][index_k2][0];
  double k_max_pt = ppt2->k3[index_k1][index_k2][k_pt_size-1];

  /* We assume that the limits for the current k3-grid were already computed in
  transfer2_get_k3_sizes() */
  int k_tr_size = ptr2->k_size_k1k2[index_k1][index_k2];
  double k_min_tr = ptr2->k_min_k1k2[index_k1][index_k2];
  double k_max_tr = ptr2->k_max_k1k2[index_k1][index_k2];
  
  /* Determine type of k-sampling. */
  /* TODO: For later times in the line-of-sight integration, where the frequency tau0-tau is small,
  one can think of using a less dense k-grid */
  double k_step_max;

  if (ptr2->k_sampling == bessel_k3_sampling) {
    /* We require the maximum step between two points to be determined by the sampling of the projection
    functions J. The argument of J is x=k*(tau0-tau) and is linearly sampled in x. Hence, the maximum
    step in k is given by the x/tau0. */
    k_step_max = pbs2->xx_step/ptr2->tau0;
  }
  else if (ptr2->k_sampling == class_transfer2_k3_sampling) {

    k_step_max = 2.*_PI_/(ptr2->tau0-ptr2->tau_rec)*ppr2->q_linstep_song;
    
    /* Old style sampling */
    if ((ppr->load_run == _TRUE_) && (ppr2->old_run == _TRUE_))
      k_step_max = 2.*_PI_/ptr2->rs_rec*ppr2->q_linstep_song;
  }

  /* - first point */

  int index_k_pt = 0;
  int index_k_tr = 0;
  double k = k_min_tr;
  k3[0] = k_min_tr;
  index_k_tr++;

  /* TODO: Magic number */
  double n = 1e99;

  /* Add unphysical points to the left of k_min_tr in order to stabilize the bispeptrum integral */
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {

    /* Count points outside the physical regime, i.e. k3 values smaller than k_min_pt. Only matters when
    k_min_tr != k_min_pt, which happens only if extrapolation is turned on */
    double first_physical_step = ppt2->k3[index_k1][index_k2][1] - ppt2->k3[index_k1][index_k2][0];
  
    while (k < k_min_pt) {

      k += MIN (k_step_max, n*first_physical_step);
      k3[index_k_tr] = k;
      index_k_tr++;    
    }

    if (k != k_min_pt) {
      k = k_min_pt;
      k3[index_k_tr - 1] = k;
    }
  }
  
  index_k_pt++;

  /* Take the next point from the sources sampling (ppt->k3) if the step is small enough.
  This assumes that ppt->k3 is a growing array. */
  while ((index_k_pt < k_pt_size) && ((ppt2->k3[index_k1][index_k2][index_k_pt] - k) < k_step_max)) {
      k = ppt2->k3[index_k1][index_k2][index_k_pt];
      k3[index_k_tr] = k;
      index_k_pt++;
      index_k_tr++;
  }

  *last_used_index_pt = index_k_pt;

  /* Make sure that we do not introduce a jump in the extrapolation */
  double last_triangular_step = ppt2->k3[index_k1][index_k2][k_pt_size-1] 
    - ppt2->k3[index_k1][index_k2][k_pt_size-2];

  /* Then, points spaced linearily with step k_step_max. Note that when using k3-extrapolation,
  we add points beyond the physical limit dictated by the triangular condition. */
  while ((index_k_tr < k_tr_size) && (k < k_max_tr)) {
    k += MIN (k_step_max, n*last_triangular_step);
    k3[index_k_tr] = k;
    index_k_tr++;
  }

  /* Consistency check on the maximum value of k3 */
  class_test (k3[k_tr_size-1] > k_max_tr,
    ptr2->error_message,
    "bug in k list calculation, k_max_transfers2=%.17f is larger than k_max_tr=%.17f, should never happen. This can happen for a race condition.",
    k3[k_tr_size-1], k_max_tr);

  /* Some debug */
  // int index_k1_debug = 0;
  // int index_k2_debug = 0;
  //
  // if ((index_k1==index_k1_debug) && (index_k2==index_k2_debug)) {
  //
  //   fprintf (stderr, "# ~~~~ (index_k1,index_k2)=(%d,%d), (k1,k2)=(%g,%g), k_tr_size=%d ~~~~~\n",
  //     index_k1, index_k2, ppt2->k[index_k1], ppt2->k[index_k2], ptr2->k_size_k1k2[index_k1][index_k2]);
  //
  //   int first_k_phys = ptr2->k_physical_start_k1k2[index_k1][index_k2];
  //   int last_k_phys = ptr2->k_physical_start_k1k2[index_k1][index_k2] + ptr2->k_physical_size_k1k2[index_k1][index_k2] - 1;
  //
  //   fprintf (stderr, "# ~~~~  K-SAMPLING OF SOURCES (k1=%g, k2=%g, %d k's in [%g,%g], of which used: %d) ~~~~~\n",
  //     ppt2->k[index_k1], ppt2->k[index_k2], k_pt_size, k_min_pt, k_max_pt, *last_used_index_pt);
  //   for (index_k_pt=0; index_k_pt < k_pt_size; ++index_k_pt)
  //     fprintf(stderr, "# %12d %26.17f\n", index_k_pt, ppt2->k3[index_k1][index_k2][index_k_pt]);
  //
  //   fprintf (stderr, "# ~~~~  LEFT-EXTRAPOLATION (%d k's in (%g,%g)) ~~~~~\n",
  //     first_k_phys, k3[0], k3[first_k_phys]);
  //   for (index_k_tr=0; index_k_tr < first_k_phys; ++index_k_tr)
  //     fprintf(stderr, "# %12d %26.17f\n", index_k_tr, k3[index_k_tr]);
  //
  //   fprintf (stderr, "# ~~~~  PHYSICAL K's (%d k's in [%g,%g]) ~~~~~\n",
  //     last_k_phys-first_k_phys+1, k3[first_k_phys], k3[last_k_phys]);
  //   for (index_k_tr=first_k_phys; index_k_tr < last_k_phys+1; ++index_k_tr)
  //     fprintf(stderr, "# %12d %26.17f\n", index_k_tr, k3[index_k_tr]);
  //
  //   fprintf (stderr, "# ~~~~  RIGHT-EXTRAPOLATION (%d k's in (%g,%g)) ~~~~~\n",
  //     k_tr_size-last_k_phys-1, k3[last_k_phys], k3[k_tr_size-1]);
  //   for (index_k_tr=last_k_phys+1; index_k_tr < k_tr_size; ++index_k_tr)
  //     fprintf(stderr, "# %12d %26.17f\n", index_k_tr, k3[index_k_tr]);
  //
  // } // end of debug

  /* Check that the computed k3-grid always satisfies the triangular condition in the physical regime */
  int physical_size = ptr2->k_physical_size_k1k2[index_k1][index_k2];
  int first_physical_index = ptr2->k_physical_start_k1k2[index_k1][index_k2];
  int last_physical_index = first_physical_index + physical_size - 1;

  for (int index_k3=first_physical_index; index_k3<=last_physical_index; ++index_k3) {
    class_test ((k3[index_k3] <= fabs(ppt2->k[index_k1]-ppt2->k[index_k2]))
      || (k3[index_k3] >= (ppt2->k[index_k1]+ppt2->k[index_k2])),
      ptr2->error_message,
      "triangular condition not satisfied in the physical regime for (index_k1,index_k2)=(%d,%d) bug!",
      index_k1, index_k2);
  }

  /* Check that k3 is strictly ascending */
  for (int index_k3=0; index_k3 < ptr2->k_size_k1k2[index_k1][index_k2]-1; ++index_k3)
    class_test (k3[index_k3] >= k3[index_k3+1],
      ptr2->error_message,
      "bug in k list calculation, k3-grid is not strictly ascending: k[%d]=%g > k[%d]=%g.",
      index_k3, index_k3+1, k3[index_k3], k3[index_k3 + 1]);

  return _SUCCESS_;

} // end of transfer2_get_k3_list




/**
 * Find the size of the k3 grid for all possible (k1,k2) pairs.
 *
 * This function allocates the arrays in the transfer2 structure indexed by
 * index_k1 and index_k2, and fills them calling transfer2_get_k3_size() in a loop.
 * It also computes ptr2->k3_size_max, which is later used to allocate several
 * arrays.
 * 
 * Please refer to the documentation of transfer2_get_k3_size() for further details.
 */
int transfer2_get_k3_sizes (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2
      )
{

  /* Check that the sources were computed for enough k3-values. Because the transfer2.c module relies
  on the linear interpolation of the source on the k3 grid, it does not make sense to proceed further
  if less than 2 points were computed. */
  for(int index_k1=0; index_k1<ppt2->k_size; ++index_k1)
    for(int index_k2=0; index_k2<index_k1+1; ++index_k2)
      class_test (
        ((ppt2->k3_size[index_k1][index_k2] < 2) && (ppr2->sources_k3_interpolation==linear_interpolation))
        || ((ppt2->k3_size[index_k1][index_k2] < 4) && (ppr2->sources_k3_interpolation==cubic_interpolation)),
        ptr2->error_message,
        "index_k1=%d, index_k2=%d: cannot do interpolation in k3 with just k3_size=%d values. Increase k3_size.",
        index_k1, index_k2, ppt2->k3_size[index_k1][index_k2]);
  
  /* The k3-dependence of the second-order transfer function T(k1,k2,k3) sets the frequency
  of oscillation of the Bessel function in the line-of-sight integral. This means that T(k1,k2,k3)
  oscillates wildly in k.  On the other hand, the k1 and k2 directions of T(k1,k2,k3) are smooth
  because the (k1,k2) dependence comes directly from the smooth line-of-sight sources.
            
  Note that, as for the sampling of the sources, the grid in k is a function of the considered
  k1 and k2 according to the triangular relation \vec{k} = \vec{k1} + \vec{k2}. */
  
  /* Allocate k1 level */
  int k1_size = ppt2->k_size;
  class_alloc(ptr2->k_size_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_physical_start_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_physical_size_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_min_k1k2, k1_size*sizeof(double *), ptr2->error_message);
  class_alloc(ptr2->k_max_k1k2, k1_size*sizeof(double *), ptr2->error_message);

  for(int index_k1=0; index_k1<k1_size; ++index_k1) {
  
    double k1 = ppt2->k[index_k1];
  
    /* Allocate k2 level */
    int k2_size = index_k1 + 1;

    class_alloc(ptr2->k_size_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_physical_start_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_physical_size_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_min_k1k2[index_k1], k2_size*sizeof(double), ptr2->error_message);
    class_alloc(ptr2->k_max_k1k2[index_k1], k2_size*sizeof(double), ptr2->error_message);
  
    /* Fill k_size_k1k2, k_min and k_max */
    for(int index_k2=0; index_k2<=index_k1; ++index_k2) {
  
      double k2 = ppt2->k[index_k2];

      /* Compute the number of points in the k3_grid for the given pair of (k1,k2) */
      class_call (transfer2_get_k3_size (
                    ppr,
                    ppr2,
                    ppt2,
                    pbs,
                    pbs2,
                    ptr2,
                    index_k1,
                    index_k2
                    ),
        ptr2->error_message,
        ptr2->error_message);

    } // end of for(index_k2)
  } // end of for(index_k1)
  
  /* We shall allocate the arrays that depend on k3 with the largest possible number of k3-values */
  ptr2->k3_size_max = 0;
  long int k_size_sum = 0;
  long int k1_k2_pairs = 0;

  for (int index_k1 = 0; index_k1 < ppt2->k_size; ++index_k1) {
    for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
      ptr2->k3_size_max = MAX (ptr2->k3_size_max, ptr2->k_size_k1k2[index_k1][index_k2]);
      k_size_sum += ptr2->k_size_k1k2[index_k1][index_k2];
      k1_k2_pairs++;
   }
  }
  
  if (ptr2->transfer2_verbose > 1)
    printf (" -> maximum size of k3-grid = %d, average size = %g\n",
      ptr2->k3_size_max, (double)k_size_sum/k1_k2_pairs);

  return _SUCCESS_;    

}





/**
 * This routine computes the transfer functions \Delta_l^{X} (k)
 * as a function of wavenumber k for each mode, initial condition,
 * type and multipole l passed in input. 
 *
 * For a given value of k, the transfer function is infered from 
 * the source function (passed in input in the array interpolated_sources)
 * and from Bessel functions (passed in input in the bessels structure),
 * either by convolving them along tau, or by a Limber appoximation.
 * This elementary task is distributed either to transfer2_integrate()
 * or to transfer2_limber(). The task of this routine is mainly to
 * loop over k values, and to decide at which k_max the calculation can
 * be stopped, according to some approximation scheme designed to find a 
 * compromise between execution time and precision. The approximation scheme
 * is defined by parameters in the precision structure.
 * 
 * @param ppr                   Input : pointer to precision structure 
 * @param ppt                   Input : pointer to perturbation structure
 * @param pbs                   Input : pointer to bessels structure 
 * @param ptr                   Input/output : pointer to transfers structure (result stored there)
 * @param index_tt              Input : index of type of transfer
 * @param index_l               Input : index of multipole
 * @param interpolated_sources  Input : array containing the sources
 * @param ptw                   Input : pointer to transfer2_workspace structure (allocated in transfer2_init() to avoid numerous reallocation) 
 * @return the error status
 */

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
        double ** interpolated_sources_in_k,
        struct transfer2_workspace * pw
        )
{
    
  /* Compute the transfer function for this set of parameters (l,m,k1,k2,k) and the type
  of transfer function given by index_tt2 */
  int transfer_type = ptr2->index_tt2_monopole[index_tt];

  /* Print some info */
  if (ptr2->transfer2_verbose > 4)
    printf("     * computing transfer function for (l,m) = (%d,%d)\n", ptr2->l[index_l], ptr2->m[index_m]);



  // =====================================================================================
  // =                         Temperature transfer function                             =
  // =====================================================================================
  
  if ((ppt2->has_cmb_temperature==_TRUE_) && (transfer_type==ptr2->index_tt2_T)) {
  
    /* Number of multipole sources to consider */
    pw->L_max = ppr2->l_max_los_t;
  
    /* Intensity does not mix with anything else, so it only has the T->T contribution */
    class_call (transfer2_integrate(
                  ppr,
                  ppr2,
                  ppt2,
                  pbs,
                  pbs2,
                  ptr2,
                  index_k1,
                  index_k2,
                  index_k,
                  index_l,
                  index_m,
                  pw->interpolated_sources_in_time,
                  pbs2->index_J_TT,   /* Temperature projection function */
                  ppt2->index_tp2_T,  /* Temperature source function */
                  pw,
                  &(pw->transfer)),
      ptr2->error_message,
      ptr2->error_message);
      
  } // end of temperature



  // =====================================================================================
  // =                     E-mode polarisation transfer function                         =
  // =====================================================================================
  
  else if ((ppt2->has_cmb_polarization_e==_TRUE_) && (transfer_type==ptr2->index_tt2_E)) {

    /* Number of multipole sources to consider */
    pw->L_max = ppr2->l_max_los_p;

    /* E-modes and B-modes mix while they propagate. There are two contributions: a
    direct one (E->E) and a mixing one (B->E) */
    double direct_contribution=0, mixing_contribution=0;
  
    /* Direct contribution (E -> E) */
    class_call (transfer2_integrate(
                  ppr,
                  ppr2,
                  ppt2,
                  pbs,
                  pbs2,
                  ptr2,
                  index_k1,
                  index_k2,
                  index_k,
                  index_l,
                  index_m,
                  pw->interpolated_sources_in_time,
                  pbs2->index_J_EE,   /* E-mode projection function */
                  ppt2->index_tp2_E,  /* E-mode source function */
                  pw,
                  &(direct_contribution)),
      ptr2->error_message,
      ptr2->error_message);

    /* Mixing contribution (B -> E), only for non-scalar modes */

    if (ptr2->m[index_m] != 0) {

      class_call (transfer2_integrate(
                    ppr,
                    ppr2,
                    ppt2,
                    pbs,
                    pbs2,
                    ptr2,
                    index_k1,
                    index_k2,
                    index_k,
                    index_l,
                    index_m,
                    pw->interpolated_sources_in_time,
                    pbs2->index_J_EB,   /* EB mixing projection function */
                    ppt2->index_tp2_B,  /* B-mode source function */
                    pw,
                    &(mixing_contribution)),
        ptr2->error_message,
        ptr2->error_message);
      
      }
      
      /* The integral is given by the sum of the E->E and B->E contributions. */
      pw->transfer = direct_contribution + mixing_contribution;
      
  } // end of E-modes



  // =====================================================================================
  // =                     B-mode polarisation transfer function                         =
  // =====================================================================================
  
  /* We consider the B-modes only for non-scalar modes. Even if didn't do so, both the 
  direct and mixed contributions would vanish. The direct contribution vanishes because
  it involves the B-mode sources at m=0, which vanish. The mixed contribution vanishes
  because the projection function vanishes for m=0 (see odd function in eq. 5.104 of
  http://arxiv.org/abs/1405.2280) */

  if ((ppt2->has_cmb_polarization_b==_TRUE_) && (transfer_type==ptr2->index_tt2_B)) {
  
    if (ptr2->m[index_m] == 0) {
    
      pw->transfer = 0;
    
    }
  
    else {

      /* Number of multipole sources to consider */
      pw->L_max = ppr2->l_max_los_p;

      /* E-modes and B-modes mix while they propagate. There are two contributions: a
      direct one B->B and a mixing one E->B */
      double direct_contribution=0, mixing_contribution=0;
  
      /* Direct contribution (B -> B). The projection function for this contribution (J_BB)
      is equal to the one for E->E (J_EE). */
      class_call (transfer2_integrate(
                    ppr,
                    ppr2,
                    ppt2,
                    pbs,
                    pbs2,
                    ptr2,
                    index_k1,
                    index_k2,
                    index_k,
                    index_l,
                    index_m,
                    pw->interpolated_sources_in_time,
                    pbs2->index_J_EE,   /* E-mode projection function (equal to B's) */
                    ppt2->index_tp2_B,  /* B-mode source function, vanishes for m=0 */
                    pw,
                    &(direct_contribution)),
        ptr2->error_message,
        ptr2->error_message);

      /* Mixing contribution (E->B). The projection function for this contribution (J_BE)
      is equal to minues the one for B->E (-J_EB). */
      class_call (transfer2_integrate(
                    ppr,
                    ppr2,
                    ppt2,
                    pbs,
                    pbs2,
                    ptr2,
                    index_k1,
                    index_k2,
                    index_k,
                    index_l,
                    index_m,
                    pw->interpolated_sources_in_time,
                    pbs2->index_J_EB,   /* EB mixing projection function (equal to minus BE),
                                           vanishes for m=0 (see eq. 5.104 of http://arxiv.org/abs/1405.2280)*/
                    ppt2->index_tp2_E,  /* E-mode source function */
                    pw,
                    &(mixing_contribution)),
        ptr2->error_message,
        ptr2->error_message);
      
        /* The integral is given by the sum of the B->B and E->B contributions.
        IMPORTANT: the mixing projection function J_BE is given by -J_EB. This follows from 
        the properties of the H matrix in appendix B of http://arxiv.org/abs/1102.1524. Hence,
        here we multiply the mixing contribution by a minus sign. */
        pw->transfer = direct_contribution - mixing_contribution;
      
    } // end of if m!=0
    
  } // end of B-modes



  // ====================================================
  // =           Rescale transfer functions             =
  // ====================================================

  /**
   * This is the best place to apply several factors on our second-order transfer functions.
   *
   * The second-order transfer functions that we just computed are the (l,m) spherical harmonic
   * expansion of the momentum-integrated distribution function (aka brightness). However, the
   * first-order transfer functions outputted by CLASS are the Legendre multipoles T_l(k,tau0),
   * which are related to the spherical harmonic  counterpart by a (2l+1) factor:
   * T_l0(k,tau0) = (2l+1) * T_l(k,tau0).  Here, we divide our second-order transfer functions
   * by a (2l+1) factor in order to match the first-order definition.
   *
   * We also divide our second-order transfer functions by an overall 4 factor, in order to
   * convert our brightness multipoles into brightness TEMPERATURE multipoles:
   * Theta_brightness = Delta/4 (see sec. 4.1 of http://arxiv.org/abs/1405.2280, especially
   * eq. 4.69, or sec. 4.1 of Pitrou et al. 2010). In the bispectra module, we shall further
   * transfrom this brightness temperature into a bolometric temperature by adding an integrated
   * correction which is quadratic at first-order.
   *   
   * Finally, we divide the transfer function by a factor 2. The reason is that we defined
   * our perturbations as X ~ X^(1) + X^(2)/2, but to compute the bispectrum we need the
   * full second-order part, that is < X^(2)/2  X^(1)  X^(1) >
   */     
          
  /* Brightness -> Brightness temperature */
  pw->transfer /= 4.;
        
  /* Y_lm expansion -> Legendre expansion */
  pw->transfer /= (2.*ptr2->l[index_l] + 1);
    
  /* Take the full second-order part of the temperature perturbation */
  pw->transfer /= 2.;



  // ===============================================
  // =           Store transfer function           =
  // ===============================================

  /* Store transfer function in transfer structure */
  ptr2->transfer[index_tt][index_k1][index_k2][index_k] = pw->transfer;


  /* Counter to keep track of the number of values fit into ptr2->transfer */
  #pragma omp atomic
  ++ptr2->count_memorised_transfers;


  return _SUCCESS_;

}





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
      )
{

  /* Shortcuts */
  int l = ptr2->l[index_l];
  int m = ptr2->m[index_m];
  double k = pw->k;
  
  /* Set the value of the transfer function to zero */
  *integral = 0;


  // =====================================================================================
  // =                            Check integration bounds                               =
  // =====================================================================================

  /* Check that we computed enough of the J's. As pbs2->xx_max is proportional to l_max,
  this is basically a test on whether we are computing enough l's for the chosen
  k-grid. If we chose to use ptr2->k to be sampled as pbs2->xx/tau0, with the option
  ptr2->k_sampling=bessels, this test will always succeed as pbs2->xx_max is by construction
  larger than pw->k*ptr2->tau0_minus_tau[0]. */
  class_test (k*pw->tau0_minus_tau[0] > pbs2->xx_max,
              ptr2->error_message,
              "not enough J's.  Increase l_max to %g or decrease k_max to %.3g.",
              ceil(k*pw->tau0_minus_tau[0]/ppr->k_max_tau0_over_l_max),
              pbs2->xx_max/pw->tau0_minus_tau[0]);
  
  /* Minimum value of x=k*(tau0-tau) above which the J_Llm(x) start to be non-negligible.  This
  always occurs for L = L_max because then J_Llm has contributions from the largest number
  of Bessels functions (J_Llm = sum over l1 of all Bessels j_l1(x) with |l-L| <= l1 <= l+L). */
  double x_min_bessel = pbs2->x_min_J[index_J][pw->L_max][index_l][index_m];
      
  /* Maximum value of x=k*(tau0-tau) needed for the line-of-sight integration.  This is linked to
  the particular k-mode we are computing. Note that ptr2->tau0_minus_tau[index_tau] has its
  largest value for index_tau=0. */
  double tau_max_needed = pw->tau0_minus_tau[0];
  double x_max_needed = k*tau_max_needed;
  
  /* The projection functions J_Llm(m) are computed only for the values of x where they are non-negligible.
  If there is no overlap between the sources and the region whereJ_Llm(x) is non-zero, return zero.
  This usually happens for high l's and small k's, as the large scales do not affect the high multipoles. */
  if (x_max_needed < x_min_bessel) {
  
    if (ptr2->transfer2_verbose > 4)
      printf("     \\ set transfer function to zero because k*tau0 << L (x_min_bessel=%g, k*tau0~%g)\n",
      x_min_bessel, x_max_needed);
  
    *integral = 0;
  
    return _SUCCESS_;
  }
      
  /* We now want to see how much overlap there is between the sources and the projection functions.
  We define an 'index_tau_max' as the time-index of the sources after which the J's become
  negligible.  If 'index_tau_max' is close to zero, then the J's are almost always negligible
  during all the sources regime.  If it is close to ppt2->tau_size, then the J's are almost never
  negligible. */
  int index_tau_max = pw->tau_grid_size-1;
  while (k*pw->tau0_minus_tau[index_tau_max] < x_min_bessel)
    index_tau_max--;
      
  if ((ptr2->transfer2_verbose > 4) && (index_tau_max!=pw->tau_grid_size-1))
    printf("     \\ adjusted integration range from tau=[%g,%g] to tau=[%g,%g]\n",
      pw->tau_grid[0], pw->tau_grid[pw->tau_grid_size-1], pw->tau_grid[0], pw->tau_grid[index_tau_max]);
  
    
  
  // =====================================================================================
  // =                             Perform the integration                               =
  // =====================================================================================
  
  /* We perform the line of sight integration starting from the initial time where the
  sources are sampled (that is ppt2->tau_sampling[0]) up to pw->tau_grid[index_tau_max]
  (see above for details). We implement the integration following the positive direction
  of time, which is the negative direction of x=k(tau0-tau) */
  
  /* Loop over time. We do not include the last point as we are cycling through integration
  intervals rather than through grid points */
  for (int index_tau=0; index_tau<index_tau_max; ++index_tau) {
  
    /* Value of x at the required time */
    double x = k * pw->tau0_minus_tau[index_tau];
      
    /* Position at the left of x in the array pbs2->xx where we sampled the J (assumes
    pbs2->xx starts from zero) */
    int index_x = (int)(x/pbs2->xx_step);
  
    /* Interpolation weight of the left point for J(x) */
    double a_J = (pbs2->xx[index_x+1] - x)/pbs2->xx_step;
    
    /* The integrand function is the sum over L of the product between the source
    S_Lm(tau) and the projection functions J_Llm(k(tau0-tau)).*/
    double integrand = 0;

  
    // **********    Sum over L    ***********

    /* The upper limit for L is controlled by the various ppr2->l_max_los_xxx, which is O(few)
    when considering scattering sources only */
    for(int index_L=0; index_L<=pw->L_max; ++index_L) {
  
      /* The 3j symbols in the definition of J force the azimuthal number 'm' to be smaller than L */
      int L = pbs2->L[index_L];
      if (abs(m) > MIN(L,l)) continue;
      
      // *** Skip if J is negligible
      
      /* The projection function J_Llm is negligible when x is too small. If this is the case, the integral
      won't give any contribution, and we can skip this L. This is just a refinement of the check we made
      above, where instead of the current x and L, we used 'x_max_needed' and 'L_max'. Note that if you do
      not include this check, you will get random segmentation faults, because you would end up addressing
      J_Llm_x with a negative x index. */

      /* Minimum value of x where J_Llm(x) is non-negligible */
      int index_x_min = pbs2->index_xmin_J[index_J][index_L][index_l][index_m];

      /* Skip if we would add only negligible contributions. In this way we also skip the case
      when we are dealing with polarisation and L<2. */
      if (index_x < index_x_min)
        continue;    
  
      /* Position of x in the pbs2->J_Llm_x array. This is now safe as we excluded negative values with
      the above check. */
      int index_x_in_J = index_x - index_x_min;
      
      /* Before accessing the projection functions, check that we won't overshoot the maximum value
      of x for which they were computed. This turns out to be an expensive check, uncomment only
      if needed. */
      // int x_size = pbs2->x_size_J[index_J][index_L][index_l][index_m];
      // class_test (index_x_in_J >= x_size,
      //   ptr2->error_message,
      //   "indexing pbs2->J_Llm_x with index_x=%d while its size is x_size=%d", index_x_in_J, x_size);

  
      // *** Interpolate the projection function in x

      /* Interpolate J in the right points for this Llm. We use the fact that the two needed values of J are contiguous in memory  */
      double * J = &(pbs2->J_Llm_x[index_J][index_L][index_l][index_m][index_x_in_J]);
      double J_Llm_left = *J;
      double J_Llm_right = *(J+1);
      
      double J_Llm;
      if (ppr->bessels_interpolation == linear_interpolation) {
        J_Llm = a_J*J_Llm_left + (1-a_J)*J_Llm_right;
      }
      else if (ppr->bessels_interpolation == cubic_interpolation) {
        double * ddJ = &(pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m][index_x_in_J]);
        double ddJ_Llm_left = *ddJ;
        double ddJ_Llm_right = *(ddJ+1);
        J_Llm = (a_J * J_Llm_left +                                
                 (1.-a_J) * (J_Llm_right
               - a_J * ((a_J+1.) * ddJ_Llm_left
                +(2.-a_J) * ddJ_Llm_right) 
               * pbs2->xx_step * pbs2->xx_step / 6.0) );
      }

      // *** Increment the integral estimate

      /* Pre-interpolated source function in tau and k3 for the desired source type */
      double source_Lm = interpolated_sources_in_time[index_source_monopole + lm(L,m)][index_tau];

      /* Increment the sum between J_Llm and S_Lm */
      integrand += J_Llm * source_Lm;
  
      /* Some debug - print the transfer functions to screen */
      // if ((index_tau == 50) && (pw->index_k1==1) && (pw->index_k2==1)) {
      //
      //   int index_xmin_J = pbs2->index_xmin_J[index_J][index_L][index_l][index_m];
      //   double x_left = pbs2->xx[index_x-index_xmin_J];
      //   double x_right = pbs2->xx[index_x-index_xmin_J+1];
      //
      //   if ((L == 2) && (l==4) && (m==0)) {
      //     if (x_left > -1) {
      //       double J_Llm_linear = a_J*J_Llm_left + (1-a_J)*J_Llm_right;
      //       fprintf(stderr, "%.17f %.17f %.17f %.17f %.17f\n",
      //         x, J_Llm, J_Llm_linear, J_Llm_left, J_Llm_right);
      //       fprintf(stderr, "JL(%d,%d,%d,x=%10.6g[%d]) = %10.6g\n",
      //         L, l, m, x_left, index_x-index_xmin_J, J_Llm_left);
      //       fprintf(stderr, "J(%d,%d,%d,x=%10.6g) = %10.6g\n",
      //         L, l, m, x, J_Llm);
      //       fprintf(stderr, "JR(%d,%d,%d,x=%10.6g[%d]) = %10.6g\n",
      //         L, l, m, x_right, index_x-index_xmin_J+1, J_Llm_right);
      //       fprintf(stderr, "\n");
      //     }
      //   }
      // }

      /* Some debug - print the integrand function for a given set of (l,m,k1,k2,k) */
      // if ( (l==100) && (m==0) ) {
      //   if ( (pw->index_k1==0) && (pw->index_k2==1) && (pw->index_k==2500) ) {
      //     fprintf (stderr, "%15f %15f\n", ptr2->tau0_minus_tau[index_tau], integrand);
      //   }
      // }

    } // end of for(index_L)    
    
    /* Increment the result with the contribution from the considered time-step */
    *integral += integrand * pw->delta_tau[index_tau];
  
    /* Some debug */
    // if (index_tau == index_tau_max-2)
    //   printf("transfer = %g\n", *integral);
  
  } // end of for(index_tau)
    
  /* Correct for factor 1/2 from the trapezoidal rule */
  *integral *= 0.5;
  
  /* Test for nans */
  class_test (isnan(*integral),
    ptr2->error_message,
    "found nan in second-order transfer function");

  return _SUCCESS_;
}



/**
 * This function interpolates sources S(k1, k2, k, tau) at the needed
 * values of k and tau.  Note that the time sampling for the integration
 * grid could be computed once and for all, as the k grid is fixed (the
 * nodes do not depend on k1 and k2).
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ppt2                  Input : pointer to 2nd-order perturbation structure
 * @param pbs                   Input : pointer to Bessel structure
 * @param ptr2                  Input : pointer to 2nd-order transfer structure
 * @param tau0                  Input : conformal time today
 * @return the error status
 */

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
      )
{

  /* Shortcuts */
  double k = pw->k;
  
  
  // **************                    Determine size of the time integration grid                ****************

  
  /* The line-of-sight integral convolves the projection functions J_Llm(k(tau0-tau)) with the
  source functions S(k1,k2,k,tau). The integration grid in time should be fine enough to catch
  the variations of both S and J. The sampling for J is linear with step pbs2->xx_step, while the
  sampling for the sources is ppt2->tau_sampling and is arbitrary.  */


  /* The limits for the time-integration-grid for the line-of-sight integral are completely determined by the sources */
  int tau_size_pt = ppt2->tau_size;
  double tau_min_pt = ppt2->tau_sampling[0];
  double tau_max_pt = ppt2->tau_sampling[tau_size_pt-1];
  
  /* Determine the maximum linear step in time */
  double tau_step_max;

  if (ptr2->tau_sampling == bessel_tau_sampling) {
    /* Merge the J and S samplings for the given k-mode in such a way that we do not oversample the J's when
    the sources grid is very fine. One could just add the two grids, but our treatment is better in the
    periods where the J and S grids have a similar stride. */
    tau_step_max = pbs2->xx_step/k;
  }
  else if (ptr2->tau_sampling == custom_transfer2_tau_sampling) {
    tau_step_max = ppr2->tau_step_trans_song/k;
  }
  
  class_test (tau_step_max == 0,
              ptr2->error_message,
              "stopping to prevent segmentation fault");
  
  /* Some debug */
  // printf("tau_step_max = %g\n", tau_step_max);

  /* We expand the grid for the sources by including extra points when the J would not be sampled well enough */
  int index_tau_tr=0;
  int index_tau_pt=0;

  /* The first point in the grid is always determined by the sources sampling */
  double tau = tau_min_pt;
  index_tau_tr++;
  index_tau_pt++;
  
  while (index_tau_pt < tau_size_pt) {
    
    while ((ppt2->tau_sampling[index_tau_pt] - tau) > tau_step_max) {
      tau += tau_step_max;
      index_tau_tr++;
    }
    
    tau = ppt2->tau_sampling[index_tau_pt];  
    ++index_tau_tr;
    ++index_tau_pt;    
  
  } // end of while(index_tau < tau_size_pt)

  pw->tau_grid_size = index_tau_tr;
  
  

  // **************                    Fill the time integration grid                ****************
  
  /* Repeat exactly the same steps as above, but this time filling the time arrays */
  index_tau_tr=0;
  index_tau_pt=0;
  
  /* We shall take note of the position of pw->tau_grid[index_tau_tr] inside ppt2->sampling by means of the index_tau_left
    array. This is very useful to spare an interpolation table look-up afterwards.. */
  pw->index_tau_left[0] = 0;
  
  /* The first point in the grid is always determined by the sources sampling */
  tau = tau_min_pt;
  pw->tau_grid[0] = tau;
  index_tau_tr++;
  index_tau_pt++;

  while (index_tau_pt < tau_size_pt) {

    while ((ppt2->tau_sampling[index_tau_pt] - tau) > tau_step_max) {

      tau += tau_step_max;
      pw->tau_grid[index_tau_tr] = tau;
      pw->index_tau_left[index_tau_tr] = index_tau_pt;
      index_tau_tr++;

      // printf("*** ADDING POINT: k=%g, index_k=%d, tau=%.7f, index_tau_pt=%d of %d, index_tau_tr=%d, index_tau_left=%d, delta_tau_pt=%g, tau_step_max=%g\n",
      //   k, index_k, tau, index_tau_pt, tau_size_pt, index_tau_tr, pw->index_tau_left[index_tau_tr], ppt2->tau_sampling[index_tau_pt] - tau, tau_step_max);
  
    }
    
    tau = ppt2->tau_sampling[index_tau_pt];
    pw->tau_grid[index_tau_tr] = tau;
    pw->index_tau_left[index_tau_tr] = index_tau_pt;
    index_tau_tr++;
    index_tau_pt++;    

    // printf("***             : k=%g, index_k=%d, tau=%.7f, index_tau_pt=%d of %d, index_tau_tr=%d, index_tau_left=%d, delta_tau_pt=%g, tau_step_max=%g\n",
    //   k, index_k, tau, index_tau_pt, tau_size_pt, index_tau_tr, pw->index_tau_left[index_tau_tr], ppt2->tau_sampling[index_tau_pt] - tau, tau_step_max);

  } // end of while(index_tau < tau_size_pt)

    
  /* The left index cannot be the last index (otherwise you'd get segmentation faults)  */
  for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size; ++index_tau_tr)
    if (pw->index_tau_left[index_tau_tr] == (ppt2->tau_size-1))
      pw->index_tau_left[index_tau_tr] = ppt2->tau_size-2;  
  
  /* Print some info */
  pw->cosk1k2 = (pw->k*pw->k - pw->k1*pw->k1 - pw->k2*pw->k2)/(2*pw->k1*pw->k2);
  if (ptr2->transfer2_verbose > 4)
    printf("     * (k,cosk1k2)=(%.3g,%.3g): the time-integration grid comprises %d+%d points from %g to %g\n",
      pw->k, pw->cosk1k2, ppt2->tau_size, pw->tau_grid_size-ppt2->tau_size, pw->tau_grid[0], pw->tau_grid[pw->tau_grid_size-1]);
  
  /* Check that the grid is sorted in ascending order */
  for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size-1; ++index_tau_tr)
    class_test (pw->tau_grid[index_tau_tr] >= pw->tau_grid[index_tau_tr+1],
              ptr2->error_message,
              "integration grid is not sorted, bug in the code");
  
  /* Check that the grid has at least as many points as in the sources grid */
  class_test (pw->tau_grid_size < ppt2->tau_size,
              ptr2->error_message,
              "integration grid is not fine enough, bug in the code");
  
  /* Check that the last and first points of the grid are equal to those of the sources sampling */
  class_test ((pw->tau_grid[0] != tau_min_pt) || (pw->tau_grid[pw->tau_grid_size-1] != tau_max_pt),
              ptr2->error_message,
              "the limits of the time integration grid do not correspond to the time limits of the sources");

    
  /* Some debug - print the integration grid */
  // int index_k_debug = 0, index_k1_debug = 0, index_k2_debug = 1;
  // 
  // fprintf (stderr, "# index_k=%d, k=%g, tau_step_max=%g, tau_min_pt=%g, tau0-tau_min_pt=%g\n",
  //   pw->index_k, pw->k, tau_step_max, tau_min_pt, ptr2->tau0-tau_min_pt);
  // 
  // if ((pw->index_k == index_k_debug) && (pw->index_k1 == index_k1_debug) && (pw->index_k2 == index_k2_debug))
  //   for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size-1; ++index_tau_tr)
  //     fprintf (stderr, "%6d %15.7g %10d\n",
  //       index_tau_tr, pw->tau_grid[index_tau_tr], pw->index_tau_left[index_tau_tr]);



  /* Fill tau0_minus_tau */
  for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size; ++index_tau_tr) 
    pw->tau0_minus_tau[index_tau_tr] = ptr2->tau0 - pw->tau_grid[index_tau_tr]; 

  /* Fill delta_tau, that is the measure for the trapezoidal integral over time */
  pw->delta_tau[0] = pw->tau_grid[1]-pw->tau_grid[0];
      
  for (index_tau_tr=1; index_tau_tr < pw->tau_grid_size-1; ++index_tau_tr)
    pw->delta_tau[index_tau_tr] = pw->tau_grid[index_tau_tr+1]-pw->tau_grid[index_tau_tr-1];
      
  pw->delta_tau[pw->tau_grid_size-1] = pw->tau_grid[pw->tau_grid_size-1]-pw->tau_grid[pw->tau_grid_size-2];

  return _SUCCESS_;
  
} // end of transfer2_get_tau_list








/**
 * This function interpolates sources S(k1, k2, k, tau) at the needed
 * values of k and tau.  Note that the time sampling for the integration
 * grid could be computed once and for all, as the k grid is fixed (the
 * nodes do not depend on k1 and k2).
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ppt2                  Input : pointer to 2nd-order perturbation structure
 * @param pbs                   Input : pointer to Bessel structure
 * @param ptr2                  Input : pointer to 2nd-order transfer structure
 * @param tau0                  Input : conformal time today
 * @return the error status
 */

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
      int index_tp2,
      double * k_grid,
      double * sources_k_spline,
      double * interpolated_sources_in_k
      )
{


  /* Shortcuts */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  double * k_pt = ppt2->k3[index_k1][index_k2];
  int k_tr_size = ptr2->k_size_k1k2[index_k1][index_k2];
  double * k_tr = k_grid;
  
  /* Cycle index */
  int index_tau;
  
  /* Find second derivative of original sources with respect to k in view of spline interpolation */
  if (ppr2->sources_k3_interpolation == cubic_interpolation) {

    class_call (array_spline_table_columns (
                  ppt2->k3[index_k1][index_k2],
                  ppt2->k3_size[index_k1][index_k2],
                  ppt2->sources[index_tp2][index_k1][index_k2],
                  ppt2->tau_size,
                  sources_k_spline,
                  _SPLINE_EST_DERIV_,
                  ptr2->error_message),
         ptr2->error_message,
         ptr2->error_message);
  }

  // =======================================================
  // =                    Interpolation                    =
  // =======================================================

  /* Limits where for which we shall interpolate the sources */
  int physical_size = ptr2->k_physical_size_k1k2[index_k1][index_k2];
  int first_physical_index = ptr2->k_physical_start_k1k2[index_k1][index_k2];
  int last_physical_index = first_physical_index + physical_size - 1;

  /* Interpolate at each k value using the usual spline interpolation algorithm */
  int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
  
  int index_k_tr;
    
  for (index_k_tr = first_physical_index; index_k_tr <= last_physical_index; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., ptr2->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;
      
    /* Interpolate for each value of conformal time */
    if (ppr2->sources_k3_interpolation == linear_interpolation) {
      for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1);
    }
    else if (ppr2->sources_k3_interpolation == cubic_interpolation) {
      for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1)
          + ((a*a*a-a) * sources_k_spline[index_tau*k_pt_size + index_k]
          +(b*b*b-b) * sources_k_spline[index_tau*k_pt_size + index_k+1])*h*h/6.0;
    }

  } // end of for (index_k_tr)


  // =======================================================
  // =                    Extrapolation                    =
  // =======================================================

  /* Extrapolate the source functions in the non-physical regime to stabilize the bispectrum integral */
  /* For m!=0, taking the first and last physical points is equivalent to taking a value which is
  equal to zero because there sin(theta_1)=sin(theta_2)=0. Hence we end up not extrapolating anything.
  We have fixed this by analytically including a 1/sin(theta)^m factor in the perturbations2.c module in 
  all the sources, so that the amplitude at the borders of the k3 range is non-vanishing. This is also 
  the factor needed to perform the bispetrum integration, ter */
    
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {
    
    /* Extrapolation on the right (k > k_max_pt) */
    for (index_k_tr = last_physical_index+1; index_k_tr < k_tr_size; ++index_k_tr) {
      
      if (ppr->bispectra_k3_extrapolation == flat_k3_extrapolation)
        for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
          interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
            interpolated_sources_in_k[last_physical_index*ppt2->tau_size + index_tau];
    }
    /* Extrapolation on the left (k < k_max_pt) */
    for (index_k_tr = 0; index_k_tr < first_physical_index; ++index_k_tr) {
      
      if (ppr->bispectra_k3_extrapolation == flat_k3_extrapolation)
        for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
          interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
            interpolated_sources_in_k[first_physical_index*ppt2->tau_size + index_tau];
    }
  }
  

  return _SUCCESS_;
  
} // end of transfer_interpolate_sources_in_k




int transfer2_interpolate_sources_in_time (
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
      )
{

  /* Shortcuts */
  int tau_size_pt = ppt2->tau_size;
  double * tau_pt = ppt2->tau_sampling;
  int tau_size_tr = pw->tau_grid_size;
  double * tau_tr = pw->tau_grid;
  
  /* Find second derivative of original sources with respect to k in view of spline interpolation */
  if (ppr2->sources_time_interpolation == cubic_interpolation) {

    class_call (array_spline_table_columns (
                  tau_pt,
                  tau_size_pt,
                  interpolated_sources_in_k + pw->index_k * tau_size_pt, /* Start from index_k */
                  1,  /* How many columns (k-indices) to consider (desired size of the slow index) */
                  sources_time_spline,
                  _SPLINE_EST_DERIV_,
                  ptr2->error_message),
         ptr2->error_message,
         ptr2->error_message);

  }

  
  /* Interpolate at each k value using the usual spline interpolation algorithm */
  int index_tau = 0;
  double h = tau_pt[index_tau+1] - tau_pt[index_tau];
  
  int index_tau_tr;
    
  for (index_tau_tr = 0; index_tau_tr < tau_size_tr; ++index_tau_tr) {
    
    while (((index_tau+1) < tau_size_pt) && (tau_pt[index_tau+1] < tau_tr[index_tau_tr])) {
      index_tau++;
      h = tau_pt[index_tau+1] - tau_pt[index_tau];
    }
    
    class_test(h==0., ptr2->error_message, "stop to avoid division by zero");
    
    double b = (tau_tr[index_tau_tr] - tau_pt[index_tau])/h;
    double a = 1.-b;

    /* Actual interpolation */
    if (ppr2->sources_time_interpolation == linear_interpolation) {

      interpolated_sources_in_time[index_tau_tr] = 
        a * interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau]
        + b * interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau+1];
    }
    else if (ppr2->sources_time_interpolation == cubic_interpolation) {

      interpolated_sources_in_time[index_tau_tr] = 
        a * interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau]
        + b * interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau+1]
        + ((a*a*a-a) * sources_time_spline[index_tau]
        +(b*b*b-b) * sources_time_spline[index_tau+1])*h*h/6.0;
    }

    /* Some debug */
    // if ((index_tp2==0) && (pw->index_k1==1) && (pw->index_k2==1) && (pw->index_k==1000)) {
    //   printf("index_tau_tr = %d, index_tau = %d\n", index_tau_tr, index_tau);
    //   printf("tau_tr = %g\n", tau_tr[index_tau_tr]);
    //   printf("interpolated_sources_in_time[index_tau_tr] = %g\n", interpolated_sources_in_time[index_tau_tr]);
    //   printf("a = %g, b = %g\n", a, b);
    //   printf("tau_pt = %g, source = %g, dd = %g\n",
    //     tau_pt[index_tau], interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau],
    //       (a*a*a-a)*sources_time_spline[index_tau]*h*h/6.0);
    //   printf("tau_pt[+1] = %g, source = %g, dd = %g\n",
    //     tau_pt[index_tau+1], interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau+1],
    //       (b*b*b-b)*sources_time_spline[index_tau+1]*h*h/6.0);
    //   printf("\n");
    // }
  
  } // end of for (index_tau_tr)
  
  /* Some debug - print the original array and the interpolation */
  // if ((index_tp2==0) && (pw->index_k1==1) && (pw->index_k2==1) && (pw->index_k==1000)) {
  // 
  //   fprintf (stderr, "\n\n");
  //   
  //   for (index_tau=0; index_tau < tau_size_pt; ++index_tau)
  //     fprintf (stderr, "%17.7g %17.7g\n", tau_pt[index_tau], interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau]);
  // 
  //   fprintf (stderr, "\n");
  // 
  //   for (index_tau_tr = 0; index_tau_tr < tau_size_tr; ++index_tau_tr)
  //     fprintf (stderr, "%17.7g %17.7g\n", tau_tr[index_tau_tr], interpolated_sources_in_time[index_tau_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  // } 

  
  return _SUCCESS_;
  
} // end of transfer_interpolate_sources_in_time





/**
 * Save the transfer functions to disk for a given transfer type.
 * 
 * The transfer functions will be saved to the file in ptr2->transfers_paths[index_tt].
 */
int transfer2_store_transfers_to_disk(
        struct perturbs2 * ppt2,
        struct transfers2 * ptr2,
        int index_k1
        )
{

  /* We shall loop over transfer types and k2 */
  int index_tt, index_k2;

  /* Print some info */
  if (ptr2->transfer2_verbose > 1)
    printf("     \\ writing transfer function for index_k1=%d ...\n", index_k1);
          
  for (index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    
    /* Open file for writing */
    // class_open (ptr2->transfers_files[index_tt], ptr2->transfers_paths[index_tt], "a+b", ptr2->error_message);

    /* Print some info */
    if (ptr2->transfer2_verbose > 3)
      printf("     * writing transfer function for (index_tt,index_k1)=(%d,%d) on '%s' ...\n",
        index_tt, index_k1, ptr2->transfers_paths[index_tt]);

    for (index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
  
      /* Write a chunk with all the k-values for this set of (type,k1,k2) */
      fwrite(
            ptr2->transfer[index_tt][index_k1][index_k2],
            sizeof(double),
            ptr2->k_size_k1k2[index_k1][index_k2],
            ptr2->transfers_files[index_tt]
            );

    } // end of for(index_k2)

    // if (ptr2->transfer2_verbose > 3)
    //   printf(" Done.\n");

    /* Close file */
    // fclose(ptr2->transfers_files[index_tt]);
      
  } // end of for(index_tt)

  // if (ptr2->transfer2_verbose == 3)
  //   printf(" Done.\n");
    
  return _SUCCESS_; 
  
}
  
  
  
  
  
/**
 * Allocate the k1 level of the transfer functions array.
 *
 * This function relies on values that are computed by transfer2_indices_of_perturbs()
 * and transfer2_get_k3_sizes(), so make sure to call them beforehand.
 *
 */
int transfer2_allocate_k1_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_k1
     )
{

  /* Issue an error if ptr2->transfers[index_k1] has already been allocated */
  class_test (ptr2->has_allocated_transfers[index_k1] == _TRUE_,
    ptr2->error_message,
    "the index_k1=%d level of ptr2->transfers is already allocated, stop to prevent error",
    index_k1);

  long int count=0;

  for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt) {

    int k1_size = ppt2->k_size;

    /* Allocate k2 level.  Note that, as for ppt2->sources, the size of this level is smaller
      than ppt2->k_size, and it depends on k1.  The reason is that we only need to compute
      the transfer functions for those k2's that are smaller than k1 (our equations are symmetrised
      wrt to k1<->k2) */
    int k2_size = index_k1 + 1;
  
    class_alloc(
      ptr2->transfer[index_tt][index_k1],
      k2_size * sizeof(double *),
      ptr2->error_message);
  
    for(int index_k2=0; index_k2<=index_k1; ++index_k2) {      

      /* Allocate k level. Note that we are using ptr2->k_size here instead of ppt2->k_size.  The
        reason is that ptr2->k is sampled much more finely than ppt2->k in order to catch the wild
        oscillation of the Bessel functions in k.  Furthermore, we only allocate memory for the k's
        that are compatible with the values of k1 and k2, i.e. we impose fabs(cosk1k2) <= 1.  See
        transfer2_get_k3_size() for further details. */
      class_alloc(
        ptr2->transfer[index_tt][index_k1][index_k2],
        ptr2->k_size_k1k2[index_k1][index_k2] * sizeof(double),
        ptr2->error_message);

      #pragma omp atomic
      ptr2->count_allocated_transfers += ptr2->k_size_k1k2[index_k1][index_k2];
      count += ptr2->k_size_k1k2[index_k1][index_k2];

    } // end of for(index_k2)
  } // end of for(index_tt)
  
  /* Print some debug information on memory consumption */
  if (ptr2->transfer2_verbose > 2) {
    printf("     * allocated ~ %.2f MB in ptr2->transfer (%ld doubles) for index_k1=%d; its size is now ~ %.3g MB;\n",
      count*sizeof(double)/1e6, count, index_k1, ptr2->count_allocated_transfers*sizeof(double)/1e6);
  }

  /* We succesfully allocated the k1 level of ptr2->transfer */
  ptr2->has_allocated_transfers[index_k1] = _TRUE_;

  return _SUCCESS_;

}



int transfer2_free_type_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_tt
     )
{

  int k1_size = ppt2->k_size;
  
  for (int index_k1=0; index_k1<k1_size; ++index_k1) {

    for (int index_k2=0; index_k2<=index_k1; ++index_k2)
      free(ptr2->transfer[index_tt][index_k1][index_k2]);

    free(ptr2->transfer[index_tt][index_k1]);
  
  } // end of for(index_k1)
  
  free(ptr2->transfer[index_tt]);

  return _SUCCESS_;

}





int transfer2_free_k1_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_k1
     )
{

  /* Issue an error if ptr2->transfers[index_k1] has already been freed */
  class_test (ptr2->has_allocated_transfers[index_k1] == _FALSE_,
    ptr2->error_message,
    "the index_k1=%d level of ptr2->transfers is already free, stop to prevent error", index_k1);

  int index_tt, index_k2;
  int k1_size = ppt2->k_size;

  long int count = 0;

  for (index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    for(index_k2=0; index_k2<=index_k1; ++index_k2) {
      free(ptr2->transfer[index_tt][index_k1][index_k2]);
      count += ptr2->k_size_k1k2[index_k1][index_k2];
    }

    free(ptr2->transfer[index_tt][index_k1]);
  }

  /* Print some debug information on memory consumption */
  if (ptr2->transfer2_verbose > 2)
    printf("     * freed ~ %.2f MB from ptr2->transfer (%ld doubles) for index_k1=%d\n",
      count*sizeof(double)/1e6, count, index_k1);

  /* We succesfully freed the k1 level of ptr2->transfer */
  ptr2->has_allocated_transfers[index_k1] = _FALSE_;

  return _SUCCESS_;

}
