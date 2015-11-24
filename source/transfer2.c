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
 * If the user specified 'store_transfers=yes', the module will save the 
 * transfer functions to disk after computing them, and then free the associated 
 * memory. To reload them from disk, use transfer2_load().
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

  if (!ppt2->has_perturbations2 ||
     (!ppt2->has_cls && !ppt2->has_cmb_bispectra && !ptr2->stop_at_transfers2)) {

    ptr2->has_cls = _FALSE_;
    if (ptr2->transfer2_verbose > 0)
      printf("Second-order transfer module skipped.\n");

    return _SUCCESS_;

  }
  else
    ptr2->has_cls = _TRUE_;

  if (ptr2->transfer2_verbose > 0)
    printf("Computing second-order transfer functions\n");


  /* Check that we have enough k1 values. Unless we are testing the code with just one
  or two values of k1, it does not make sense to continue with less than 4 values, because
  the modules that follow rely on the spline interpolation on the k1 grid. */

  if ((ppt2->k_size < 4) && (ptr2->stop_at_transfers2 == _FALSE_)) {
    printf("# WARNING: cannot do cubic interpolation in k1 with less than 4 values. ");
    printf("Increase k_size, or do not trust results for 1st-order transfer functions.\n");
  }

  /* Get conformal age, recombination time and comoving sound horizon at recombination
  from the background and thermodynamics structures (only place where these structures
  are used in this module) */
  ptr2->tau0 = pba->conformal_age;
  ptr2->tau_rec = pth->tau_rec;
  ptr2->rs_rec = pth->rs_rec;



  // ====================================================================================
  // =                             Indices and k-samplings                              =
  // ===================================================================================

  /* Initialize all indices in the transfers structure and allocate all its arrays. This
  function calls transfer2_get_lm_lists() that fills ptr2->l, and transfer2_get_k3_sizes()
  that finds for all (k1,k2) pairs the allowed k3 range. */

  class_call (transfer2_indices_of_transfers (ppr,ppr2,ppt2,pbs,pbs2,ptr,ptr2),
    ptr2->error_message,
    ptr2->error_message);


  /* Apart from ptr2->transfer, all the arrays needed by the subsequent modules have been filled.
  If the user requested to load the transfer functions from disk, we can stop the execution of
  this module now without regrets. */

  if (ppr2->load_transfers == _TRUE_) {
    
    if (ptr2->transfer2_verbose > 0)
      printf(" -> leaving transfer2 module; transfer functions will be read from disk\n");
    
    return _SUCCESS_;
  }


  // ==================================================================================
  // =                               Allocate workspaces                              =
  // ==================================================================================

  /* - Parallelization variables */

  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;

  #ifdef _OPENMP
  #pragma omp parallel
  number_of_threads = omp_get_num_threads();
  #endif
  

  /* - Allocate a workspace per thread */
    
  struct transfer2_workspace ** ppw;
  class_alloc(ppw, number_of_threads*sizeof(struct transfer2_workspace*), ptr2->error_message);
    
  /* We shall allocate the time arrays with the maximum possible number of integration steps */
  int tau_size_max;

  /* In the sources time sampling, the integration grid matches the sources time sampling */
  if (ptr2->tau_sampling == sources_tau_sampling) {
    tau_size_max = ppt2->tau_size;
  }

  /* In the custom time sampling, the number of steps depend on the considered wavemode; the
  largest k will have more time steps because the projection function, J(k(tau0-tau)),
  oscillates faster in time. */
  else if (ptr2->tau_sampling == custom_tau_sampling) {
    double k_max = ptr2->k_max_k1k2[ppt2->k_size-1][ppt2->k_size-1];
    double tau_step_min = 2*_PI_/k_max*ppr2->tau_linstep_song;
    double tau_max = ppt2->tau_sampling[ppt2->tau_size-1];
    tau_size_max = ppt2->tau_size + ceil(tau_max/tau_step_min) + 1;
  }

  /* In the bessel sampling, we take as many steps as the points where the Bessel functions
  are sampled */
  else if (ptr2->tau_sampling == bessel_tau_sampling) {
    tau_size_max = ppt2->tau_size + pbs2->xx_size;
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

    /* Allocate the k-grid with the maximum possible number of k3-values  */
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
    
    /* Allocate delta_tau, the trapezoidal measure of the line-of-sight integral. This array
    is defined as tau(i+1)-tau(i-1) except for the first and last elements, which are,
    respectively, tau(1)-tau(0) and tau(N)-tau(N-1). */
    class_alloc_parallel(
      ppw[thread]->delta_tau,
      tau_size_max*sizeof(double),
      ptr2->error_message);
    
    /* Allocate index_tau_left, an array useful for the time interpolation of the sources (see header file) */
    class_alloc_parallel(
      ppw[thread]->index_tau_left,
      tau_size_max*sizeof(double),
      ptr2->error_message);

    /* Allocate the array that will contain the second derivatives of the sources with respect to time,
    in view of spline interpolation */
    class_alloc_parallel (
      ppw[thread]->sources_time_spline,
      ppt2->tp2_size*sizeof(double *),
      ptr2->error_message);

    for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp)
      class_alloc_parallel(
        ppw[thread]->sources_time_spline[index_tp],
        tau_size_max*sizeof(double),
        ptr2->error_message);

    /* Allocate the array that will contain the sources interpolated at the exact time-values in the integration grid */
    class_alloc_parallel (
      ppw[thread]->interpolated_sources_in_time,
      ppt2->tp2_size*sizeof(double *),
      ptr2->error_message);

    for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp)
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

  /* Array that will contain the second derivative of the sources with respect to k3,
  in view of spline interpolation */
  double ** sources_k_spline;
  class_alloc (sources_k_spline, ppt2->tp2_size*sizeof(double *), ptr2->error_message);
  
  /* Array that will contain the interpolated sources at the exact k3-values needed
  of the intergration grid */
  double ** interpolated_sources_in_k;
  class_alloc (interpolated_sources_in_k, ppt2->tp2_size*sizeof(double *), ptr2->error_message);
  


  // =====================================================================================
  // =                             Main loop on (k1,k2,tt2)                             =
  // =====================================================================================

  /* We shall now compute the transfer function array (ptr2->transfer) by calling the
  transfer2_compute() function in a loop over its levels. The order of the loops is
  index_k1, index_k2 and index_tt2. The latter includes the transfer type (T,E,B) and
  the multipole indices (l,m). */

  if (ptr2->transfer2_verbose > 0)
    printf(" -> starting actual computation of second-order transfer functions\n");

  for (int index_k1 = 0; index_k1 < ppt2->k_size; ++index_k1) {

    if (ptr2->transfer2_verbose > 1)
      printf ("     * computing transfer functions today for index_k1=%d of %d, k1=%g\n",
        index_k1, ppt2->k_size, ppt2->k[index_k1]);

    /* Allocate the remaining levels of ptr2->transfer */
    class_call(transfer2_allocate_k1_level(ppt2, ptr2, index_k1),
      ptr2->error_message,
      ptr2->error_message);


    // -----------------------------------------------------------------------------
    // -                            Load sources from disk                         -
    // -----------------------------------------------------------------------------
    
    /* Load sources from disk if needed */
    class_call(perturb2_load(ppt2, index_k1),
        ppt2->error_message,
        ptr2->error_message);
          

    /* We only need to consider those k2 that are equal to or larger than k1,
    as the quadratic sources were symmetrised in the perturbation2.c module */
    for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {

      if (ptr2->transfer2_verbose > 2)
        printf(" -> computing transfer function for (k1,k2) = (%.3g,%.3g)\n", ppt2->k[index_k1], ppt2->k[index_k2]);


      // -----------------------------------------------------------------------------
      // -                        Interpolate sources in k3                          -
      // -----------------------------------------------------------------------------

      /* Find the integration grid in k3 for the current (k1,k2) pair */
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
          ptr2->k_size_k1k2[index_k1][index_k2]
            - ptr2->k_physical_size_k1k2[index_k1][index_k2] - ptr2->k_physical_start_k1k2[index_k1][index_k2],
          ppw[thread]->k_grid[0], ppw[thread]->k_grid[ptr2->k_size_k1k2[index_k1][index_k2]-1]);

      free (k_grid_temp);
      

      // -----------------------------------------------------------------------------
      // -                          Interpolate sources in k                         -
      // -----------------------------------------------------------------------------

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
                      ppw[0]->k_grid, /* Grid of desired k-values for integration; all ppw->[thread]->k_grid are filled with the same values */
                      sources_k_spline[index_tp], /* Will be filled with second-order derivatives */
                      interpolated_sources_in_k[index_tp] /* Will be filled with interpolated values in ptr2->k(k1,k2) */
                      ),
          ptr2->error_message,
          ptr2->error_message);
      
      } // for (index_tp)

                    
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
                          interpolated_sources_in_k[index_tp], /* Must be already filled by transfer2_interpolate_sources_in_k() */
                          ppw[thread]->sources_time_spline[index_tp], /* Will be filled with second-order derivatives */
                          ppw[thread]->interpolated_sources_in_time[index_tp], /* Will be filled with interpolated values in pw->tau_grid */
                          ppw[thread]
                          ),
              ptr2->error_message,
              ptr2->error_message);
      
          } // for (index_tp)

          
          // ----------------------------------------------------------------------------
          // -                         Compute transfer functions                       -
          // ----------------------------------------------------------------------------

          /* Now that we have interpolated the source function as a function of time in
          the right point of k, we have all the ingredients to compute the second-order
          transfer functions. We do so by looping over index_tt, the composite index that
          includes both the field (T,E,B) and multipole (l,m) dependences. */

          for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
  
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
                                   ptr2->tt2_to_index_l[index_tt],
                                   ptr2->tt2_to_index_m[index_tt],
                                   index_tt,
                                   ppw[thread]->interpolated_sources_in_time,
                                   ppw[thread]
                                   ),
              ptr2->error_message,
              ptr2->error_message);

          } // for(index_tt)

          #pragma omp flush(abort)

        } // for(index_k) 

      } if (abort) return _FAILURE_; /* end of parallel region */

      /* Free the memory for the interpolated sources */
      for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
        free(sources_k_spline[index_tp]);
        free(interpolated_sources_in_k[index_tp]);
      }

    } // for(index_k2)

    /* Free the memory associated with the source function for the considered k1.
    We won't need them anymore because the different k1 modes are independent. */
    class_call (perturb2_free_k1_level (ppt2, index_k1),
      ppt2->error_message, ppt2->error_message);

    /* Save all transfer functions for the given k1, and free the memory associated with
    them. The next time we'll need them, we shall load them from disk.  */
    if (ppr2->store_transfers) {
      
      class_call (transfer2_store (ppt2, ptr2, index_k1),
          ptr2->error_message,
          ptr2->error_message);

      class_call (transfer2_free_k1_level (ppt2, ptr2, index_k1),
        ptr2->error_message, ptr2->error_message);
    }

  } // for(index_k1)

  /* Mark the transfer functions as ready to be used */
  for (int index_tt=0; index_tt < ptr2->tt2_size; ++index_tt)
    ptr2->transfers_available[index_tt] = _TRUE_;



  // ====================================================================================
  // =                                  Clean & exit                                    =
  // ====================================================================================

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
  
  if (ptr2->transfer2_verbose > 1)
    printf (" -> filled ptr2->transfer with %ld values (%g MB)\n",
      ptr2->count_memorised_transfers, ptr2->count_memorised_transfers/1e6*8);
  
  /* Check that the number of filled values corresponds to the number of allocated space */
  if (ppr2->load_transfers)
    class_test (ptr2->count_allocated_transfers != ptr2->count_memorised_transfers,
      ptr2->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!", ptr2->count_allocated_transfers, ptr2->count_memorised_transfers);

  /* Do not evaluate the subsequent modules if ppt2->has_transfers2_only is true */
  if (ptr2->stop_at_transfers2) {
    ppt->has_cls = _FALSE_;
    ppt->has_cmb_bispectra = _FALSE_;
    ppt2->has_cls = _FALSE_;
    ppt2->has_pks = _FALSE_;
    ppt2->has_cmb_bispectra = _FALSE_;
  }
    
  return _SUCCESS_;
}


/**
 * Allocate all levels beyond the transfer-type level of the transfer functions array.
 *
 * Allocate the transfer type (tt) level of the array that will contain the second-order
 * transfer functions: ptr2->transfer[index_tt][index_k1][index_k2][index_k].
 *
 * This function makes space for the transfer functions; it is called before loading
 * them from disk in transfer2_load(). It relies on values that are
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

  /* Allocate memory only if needed */
  if (ptr2->transfers_allocated[index_tt])
    return _SUCCESS_;

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

    } // for(index_k2)
  } // for(index_k1)
  
  /* Print some debug information on memory consumption */
  if (ptr2->transfer2_verbose > 2) {
    printf("     * allocated ~ %.2f MB in ptr2->transfer (%ld doubles) for index_tt=%d; it's size is now ~ %.3g MB;\n",
      count*sizeof(double)/1e6, count, index_tt, ptr2->count_allocated_transfers*sizeof(double)/1e6);
  }

  /* We succesfully allocated the k1 level of ppt2->sources */
  ptr2->transfers_allocated[index_tt] = _TRUE_;

  return _SUCCESS_;

}



/**
 * Load the source function T^X_lm(k1,k2,k3) from disk for a given (X,l,m) value.
 *
 * See the documentation in perturbations2.h (\ref StorageFiles) for more
 * details.
 */

int transfer2_load(
      struct perturbs2 * ppt2,
      struct transfers2 * ptr2,
      int index_tt
      )
{

  /* Load only if needed */
  if (ptr2->transfers_available[index_tt])
    return _SUCCESS_;

  /* Print some info */
  if (ptr2->transfer2_verbose > 2)
    printf("     * transfer2_load: reading results for index_tt=%d from '%s' ...",
      index_tt, ptr2->storage_paths[index_tt]);
  
  /* Complain if there is no file to load */
  struct stat st;
  class_test (stat (ptr2->storage_paths[index_tt], &st) != 0,
		ptr2->error_message,
		"cannot load transfer functions for index_tt=%d, file '%s' does not exist",
    index_tt, ptr2->storage_paths[index_tt]);

  /* Make space */
  class_call (transfer2_allocate_type_level(ppt2, ptr2, index_tt),
    ptr2->error_message,
    ptr2->error_message);

  /* Open file for reading */
  class_open (ptr2->storage_files[index_tt],
    ptr2->storage_paths[index_tt],
    "rb", ptr2->error_message);
  
  /* Read from file */
  for (int index_k1 = 0; index_k1 < ppt2->k_size; ++index_k1) {  
    for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
  
      int n_to_read = ptr2->k_size_k1k2[index_k1][index_k2];
  
      /* Debug: print some info */
      // printf ("reading %d entries for index_tt=%d, index_k1=%d, index_k2=%d\n",
      //   n_to_read, index_tt, index_k1, index_k2);
  
      /* Read a chunk with all the k-values for this set of (type,k1,k2) */
      int n_read = fread(
              ptr2->transfer[index_tt][index_k1][index_k2],
              sizeof(double),
              n_to_read,
              ptr2->storage_files[index_tt]);
  
      class_test(n_read != n_to_read,
        ptr2->error_message,
        "Could not read in '%s' file, read %d entries but expected %d (index_tt=%d,index_k1=%d,index_k2=%d)",
          ptr2->storage_paths[index_tt], n_read, n_to_read, index_tt, index_k1, index_k2);
  
      /* Update the counter for the values stored in ptr2->transfers */
      #pragma omp atomic
      ptr2->count_memorised_transfers += ptr2->k_size_k1k2[index_k1][index_k2];
  
    } 
  }
  
  /* Close file */
  fclose (ptr2->storage_files[index_tt]);
  
  /* Transfers are now ready to use */
  ptr2->transfers_available[index_tt] = _TRUE_;

  return _SUCCESS_; 
  
}






/**
 * Free all the memory space allocated by transfer2_init().
 */ 

int transfer2_free(
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct transfers2 * ptr2
      )
{

  if (ptr2->has_cls == _TRUE_) {
  
    free (ptr2->tt2_labels);
  
    int k1_size = ppt2->k_size;

    for (int index_tt = 0; index_tt < ptr2->tt2_size; ++index_tt) {

      class_call (transfer2_free_type_level(ppt2, ptr2, index_tt),
          ptr2->error_message,
          ptr2->error_message);

    }

    free (ptr2->transfer);
    free (ptr2->transfers_allocated);
    free (ptr2->transfers_available);

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
      
    free (ptr2->l);
    free (ptr2->m);

    free (ptr2->tt2_to_index_l);
    free (ptr2->tt2_to_index_m);
    free (ptr2->tt2_start);
    free (ptr2->tp2_start);
  
    for (int index_l=0; index_l<ptr2->l_size; ++index_l)
      free (ptr2->lm_array[index_l]);
    free (ptr2->lm_array);

    /* Free file arrays */
    if (ppr2->store_transfers || ppr2->load_transfers) {

      for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt)
        free (ptr2->storage_paths[index_tt]);

      free (ptr2->storage_files);
      free (ptr2->storage_paths);
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

  /* Get l values where to compute the transfer functions */
  class_call (transfer2_get_l_list(ppr,ppr2,ppt2,pbs,pbs2,ptr2),
    ptr2->error_message,
    ptr2->error_message);
  
  /* Get the m values where to compute the transfer functions, copying them
  from ppr2->m. */
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
  ptr2->n_transfers = size_indexl_indexm (ptr2->l, ptr2->l_size, ppr2->m, ppr2->m_size);
  
  /* Number of non-vanishing E-mode transfer functions, for debug purposes */
  ptr2->n_nonzero_transfers_E = ptr2->n_transfers;
  if (ppr2->compute_m[0])
    ptr2->n_nonzero_transfers_E -= 2;
  if (ppr2->compute_m[1])
    ptr2->n_nonzero_transfers_E -= 1;

  /* Number of non-vanishing B-mode transfer functions, for debug purposes */
  ptr2->n_nonzero_transfers_B = ptr2->n_transfers;
  if (ppr2->compute_m[0])
    ptr2->n_nonzero_transfers_B -= ptr2->l_size;
  if (ppr2->compute_m[1])
    ptr2->n_nonzero_transfers_B -= 1;

  /* Photon temperature transfer functions */
  if (ppt2->has_source_T) {

    ptr2->index_tt2_T = index_tt;
    index_tt += ptr2->n_transfers; 
  }

  /* Photon E-mode polarisation transfer functions */
  if (ppt2->has_source_E) {
  
    ptr2->index_tt2_E = index_tt;
    index_tt += ptr2->n_transfers;    
  }

  /* Photon B-mode polarisation transfer functions. They will be computed
  only for non scalar modes, otherwise they just vanish. */
  if (ppt2->has_source_B) {

    ptr2->index_tt2_B = index_tt;
    index_tt += ptr2->n_transfers;
  }


  /* Total number of transfer functions to compute */
  ptr2->tt2_size = index_tt;

  /* Initialise the labels of the transfer types */
  class_calloc (ptr2->tt2_labels,
    ptr2->tt2_size*_MAX_LENGTH_LABEL_,
    sizeof(char),
    ptr2->error_message);

  if (ptr2->transfer2_verbose > 1) {
    printf (" -> will compute tt2_size=%d transfer functions: ", ptr2->tt2_size);
    if (ppt2->has_source_B == _TRUE_) printf ("T=%d ", ptr2->n_transfers);
    if (ppt2->has_source_E == _TRUE_) printf ("E=%d (%d non-zero) ",
      ptr2->n_transfers, ptr2->n_nonzero_transfers_E);
    if (ppt2->has_source_B == _TRUE_) printf ("B=%d (%d non-zero) ",
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
  
  /* The array where we shall store the second-order transfer function, ptr2->transfer, has
  five levels that should be indexed in the following way:
    
        ptr2->transfer [index_tt]
                       [index_k1]
                       [index_k2]
                       [index_k]
    
  index_tt includes both the transfer type (T,E,B) and the multipole indices (l,m).
  index_k1 goes from 0 to ppt2->k_size.    
  index_k2 goes from 0 to ppt2->k_size-index_k1.
  index_k3 goes from 0 to ppt2->k3_size[index_k1][index_k2]. */

  /* Counter to keep track of the memory usage of ptr2->transfer */
  ptr2->count_memorised_transfers = 0;
  ptr2->count_allocated_transfers = 0;
    
  /* Allocate transfer-type (tt2) level */  
  class_alloc (
    ptr2->transfer,
    ptr2->tt2_size * sizeof(double ***),
    ptr2->error_message);

  /* Allocate the index_tt level. The remaining k1, k2 and k levels will be allocated later
  using the function transfer2_allocate_k1_level(). */
  for (int index_tt = 0; index_tt < ptr2->tt2_size; ++index_tt) {
  
    int k1_size = ppt2->k_size;
    
    class_alloc (
      ptr2->transfer[index_tt],
      k1_size * sizeof(double **),
      ptr2->error_message);
  
  }
  
  /* Allocate and initialize the logical arrays keeping track of the state of ptr2->transfer */
  class_calloc(ptr2->transfers_allocated, ptr2->tt2_size, sizeof(short), ptr2->error_message);
  class_calloc(ptr2->transfers_available, ptr2->tt2_size, sizeof(short), ptr2->error_message);



  // ==================================================================================
  // =                             Create transfers files                             =
  // ==================================================================================  

  /* Create the files to store the transfer functions in */
  
  if ((ppr2->store_transfers == _TRUE_) || (ppr2->load_transfers == _TRUE_)) {

    /* We are going to store the transfers in n=k_size files, one for each requested k1 */
    class_alloc (ptr2->storage_files, ptr2->tt2_size*sizeof(FILE *), ptr2->error_message);
    class_alloc (ptr2->storage_paths, ptr2->tt2_size*sizeof(char *), ptr2->error_message);

    for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt) {
      
      /* The name of each transfers file will have the tt index in it */
      class_alloc (ptr2->storage_paths[index_tt], _FILENAMESIZE_*sizeof(char), ptr2->error_message);
      sprintf (ptr2->storage_paths[index_tt], "%s/transfers_%03d.dat", ptr2->storage_dir, index_tt);

      /* Uncomment to open the transfer storage files immediately, rather than in the store
      function. This might be an improvement in speed, becuase in this way we would avoid
      opening and closing all the transfer storage files for each k1 value. */
      // if (ppr2->store_transfers)
      //   class_open (ptr2->storage_files[index_tt],
      //     ptr2->storage_paths[index_tt],
      //     "wb", ptr2->error_message);
      
    }

    if (ptr2->transfer2_verbose > 2)
      printf ("     * created %d files to store transfer functions\n", ptr2->tt2_size);

  }
  

  return _SUCCESS_;

}




/**
 * Determine the l-values where the transfer functions will be computed.
 * 
 * This function basically copies the l-sampling of the first-order transfer
 * functions (ptr->l, which coincides with pbs->l) into ptr2->l.  This original
 * sampling was computed in transfer_get_l_list(), and it is just a log + linear
 * sampling.
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
  
  /* Copy the l's from pbs->l into ptr2->l */
  class_alloc(ptr2->l, ptr2->l_size*sizeof(int), ptr2->error_message);

  for (int index_l=0; index_l < ptr2->l_size; index_l++) {
    ptr2->l[index_l] = pbs->l[index_l];
  }
  
  return _SUCCESS_;

}




/**
 * Fill the arrays that deal with the (l,m) indexing of the ptr2->transfer
 * array, and create a correspondence between source and transfer functions.
 * 
 * The following arrays will be filled:
 *
 * - ptr2->lm_array[index_l][index_m], which contains the index associated with a given (l,m) couple
 *   and is used to access the ptr2->transfer array.
 * - ptr2->tt2_to_index_l[index_tt2], which gives the index in ptr2->l for the 
 *   considered transfer type.
 * - ptr2->tt2_to_index_m[index_tt2], which gives the index in ptr2->m for the 
 *   considered transfer type.
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

  /* Fill ptr2->lm_array, which contains the index associated with a given (l,m) couple.
  lm_array is useful every time you need a specific (l,m) value from the transfer
  function array: ptr2->transfer[index_type + ptr2->lm_array[index_l][index_m]]. */

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
                                                           
      /* Debug - Print the lm_array */
      // printf ("(l,m) = (%d,%d) corresponds to an offset of %d\n",
      //   ptr2->l[index_l], ptr2->m[index_m], ptr2->lm_array[index_l][index_m]);

    } // end of for(index_m)
  } // end of for(index_l)



  // ======================================================================================
  // =                            Transfer-source correspondence                          =
  // ======================================================================================

  /* In the transfer2.c module, a transfer type (index_tt2) encloses both the multipole
  (l,m) and the source type (temperature, polarization, etc). The lm_array we
  have defined above serve the purpose of generating an univocal index_tt2 given the
  (l,m) multipole. Here, we do the inverse operation by univocally associating to index_tt2
  the multipole indices (l,m) */

  /* For a given transfer type index index_tt, these arrays contain the corresponding
  indices in ptr2->l and ptr2->m */
  class_alloc (ptr2->tt2_to_index_l, ptr2->tt2_size*sizeof(int), ptr2->error_message);
  class_alloc (ptr2->tt2_to_index_m, ptr2->tt2_size*sizeof(int), ptr2->error_message);

  /* For a given transfer type index index_tt, these arrays contain the indices corresponding
  to the monopole of index_tt in ptr2->transfers and ppt2->sources, respectively. */
  class_alloc (ptr2->tt2_start, ptr2->tt2_size*sizeof(int), ptr2->error_message);
  class_alloc (ptr2->tp2_start, ptr2->tt2_size*sizeof(int), ptr2->error_message);

  for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    
    /* Initialise the current indices to -1 */
    int index_l = -1;
    int index_m = -1;


    /* - Photon temperature */
    if ((ppt2->has_source_T==_TRUE_)
    && (index_tt >= ptr2->index_tt2_T) && (index_tt < ptr2->index_tt2_T+ptr2->n_transfers)) {

      /* Store the position of the temperature monopole in ppt2->sources and ptr2->transfer */
      ptr2->tp2_start[index_tt] = ppt2->index_tp2_T;
      ptr2->tt2_start[index_tt] = ptr2->index_tt2_T;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->tt2_start[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "T_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("T, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->tt2_start[index_tt], ptr2->tp2_start[index_tt]);
    }


    /* - Photon E-mode polarization */
    
    else if ((ppt2->has_source_E==_TRUE_)
    && (index_tt >= ptr2->index_tt2_E) && (index_tt < ptr2->index_tt2_E+ptr2->n_transfers)) {

      /* Store the position of the E-modes monopole in ppt2->sources and ptr2->transfer */
      ptr2->tp2_start[index_tt] = ppt2->index_tp2_E;
      ptr2->tt2_start[index_tt] = ptr2->index_tt2_E;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->tt2_start[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "E_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("E, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->tt2_start[index_tt], ptr2->tp2_start[index_tt]);
    }


    /* - Photon B-mode polarization */

    else if ((ppt2->has_source_B==_TRUE_)
    && (index_tt >= ptr2->index_tt2_B) && (index_tt < ptr2->index_tt2_B+ptr2->n_transfers)) {

      /* Store the position of the B-modes monopole in ppt2->sources and ptr2->transfer */
      ptr2->tp2_start[index_tt] = ppt2->index_tp2_B;
      ptr2->tt2_start[index_tt] = ptr2->index_tt2_B;

      /* Find (l,m) associated with index_tt */
      int lm_offset = index_tt - ptr2->tt2_start[index_tt];
      offset2multipole_indexl_indexm (lm_offset, ptr2->l, ptr2->l_size, ptr2->m, ptr2->m_size,
        &index_l, &index_m);

      /* Set the labels of the transfer types */
      sprintf(ptr2->tt2_labels[index_tt], "B_%d_%d",ptr2->l[index_l],ptr2->m[index_m]);

      /* Some debug */
      // printf("B, index_tt=%d: lm_offset=%d -> (%d,%d), label=%s, monopole_tr=%d, monopole_pt=%d\n",
      //   index_tt, lm_offset, ptr2->l[index_l], ptr2->m[index_m],
      //   ptr2->tt2_labels[index_tt],
      //   ptr2->tt2_start[index_tt], ptr2->tp2_start[index_tt]);
    }

    /* Check the result */
    class_test ((index_l>=ptr2->l_size) || (index_m>=ptr2->m_size) || (index_l<0) || (index_m<0),
      ptr2->error_message,
      "index_tt=%d: result (index_l,index_m)=(%d,%d) is out of bounds index_l=[%d,%d], index_m=[%d,%d]\n",
      index_tt, index_l, index_m, 0, ptr2->l_size-1, 0, ptr2->m_size-1);

    /* Write the result */
    ptr2->tt2_to_index_l[index_tt] = index_l;
    ptr2->tt2_to_index_m[index_tt] = index_m;
    
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
 * More in detail, the oscillation frequency of the transfer function T(k1,k2,k3) in the k3
 * direction is dictated by the Bessel functions in the line of sight integral. These have
 * argument k3*(tau_0-tau); since most of the sources are localised at recombination, it results
 * that k3 oscillates with a frequency of tau_0-tau_rec~tau_0. The other directions, k1 and
 * k2, inherit the oscillation frequency of the sources at recombination, which is roughly
 * tau_rec/sqrt(3) because of tight coupling. Therefore, the k3 direction has a frequency
 * tau_0/tau_rec~80 times larger than the k1 and k2 directions. On the other hand, for the
 * k1 and k2 directions we keep the same sampling as in the sources, that is, ppt2->k.
 * 
 * Extrapolation: to solve the bispectrum integral (eq. 6.36) it is useful to extrapolate
 * the transfer functions in the k3 direction beyond their physical limits, that is, beyond
 * k3_min = |k1-k2| and k3_max = k1+k2. The reason is purely numerical: eventually, the
 * contributions from the extrapolated regions will cancel. The extrapolation has the purpose
 * of stabilizing an problematic integration.
 *
 * TODO: For later times in the line-of-sight integration, where the frequency tau0-tau is
 * small, one can think of using a less dense k3-grid.
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


  // ====================================================================================
  // =                               Determine k3 limits                                =
  // ====================================================================================

  /* Extract the minimum and maximum values of k3 computed in the perturbations2.c module.
  They are computed taking into consideration the triangular condition, whereby
  \vec{k3} = \vec{k1} + \vec{k2}, and the overall minimum and maximum k-values
  (ppt2->k[0] and ppt2->k[ppt2->k_size-1], respectively). */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  double k_min_pt = ppt2->k3[index_k1][index_k2][0];
  double k_max_pt = ppt2->k3[index_k1][index_k2][k_pt_size-1];

  /* By default, for the transfer functions we take the same k-limits used for the sources */
  double k_min_tr = k_min_pt;
  double k_max_tr = k_max_pt;

  /* When computing spectra and bispectra, it is convenient to extend the k-grid beyond the
  physical limits dictated by the triangular condition in order to make the successive
  integrations numerically stabler. Here we do such extension of k_min_tr and k_max_tr. */
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

    /* Debug - Print the change in the k3-range */
    // printf("physical_k3_range=%g, extended_k3_range_right=%g\n", physical_k3_range, extended_k3_range_right);
    // printf("PRE:  k_max_tr = %g\n", k_max_pt);
    // printf("POST: k_max_tr = %g\n", k_max_tr);

  }

  /* Update the transfer structure with the k3 limits */
  ptr2->k_min_k1k2[index_k1][index_k2] = k_min_tr;
  ptr2->k_max_k1k2[index_k1][index_k2] = k_max_tr;
  
  /* Check that we the grid where we sampled the the projection functions J can accomodate for the required
  values of k3. We do not need to test the lower limit as the Bessel grid always starts from zero. */
  class_test (k_max_tr > pbs2->xx_max/ptr2->tau0,
    ptr2->error_message,
    "xx_max is not large enough. The projection functions J(k*(tau0-tau)) do not cover the needed k-range\
 (k_max=%g,tau0=%g,xx_max=%g). Maybe you are using too much extrapolation.", k_max_tr, ptr2->tau0, pbs2->xx_max);



  // ====================================================================================
  // =                            Determine linear step                                 =
  // ====================================================================================

  /* Find the linear step to use for the k3 sampling, based on the requested type of
  k-sampling. */
  double k_step_max;

  /* In the CLASS-type sampling, we reason in terms of oscillations. We know that, for sources
  peaked at recombination, the transfer function oscillates in k3 with a frequency of 
  roughly tau0-tau_rec. Here we choose a linear step in k3 which is proportional to such
  frequency, via the parameter q_linstep_song. */

  if (ptr2->k_sampling == class_transfer2_k3_sampling) {

    k_step_max = 2.*_PI_/(ptr2->tau0-ptr2->tau_rec)*ppr2->q_linstep_song;
    
    /* In older versions, the step was determined with respect to the comoving sound horizon
    at recombination. */
    if ((ppr->load_run == _TRUE_) && (ppr2->old_run == _TRUE_))
      k_step_max = 2.*_PI_/ptr2->rs_rec*ppr2->q_linstep_song;

  }

  /* In the Bessel-type sampling, we require the maximum step between two points to be
  determined by the sampling of the projection functions J, which in turn is controlled
  by the user via the parameter bessel_x_step_song. The argument of J is x=k*(tau0-tau)
  and is linearly sampled in x. Hence, the maximum step in k is given by the x/tau0. */

  else if (ptr2->k_sampling == bessel_k3_sampling) {
    
    k_step_max = pbs2->xx_step/ptr2->tau0;

  }

  class_test(k_step_max == 0,
    ptr2->error_message,
    "stop to avoid infinite loop");


  // ====================================================================================
  // =                              Determine k3 size                                   =
  // ====================================================================================

  /* We shall sample k3 in all the points that were used for the line of sight sources,
  plus all the points from a linear sampling with a step of k_step_max. */
  
  int index_k_pt = 0;
  int index_k_tr = 0;  
  double k = k_min_tr;
  index_k_tr++;


  /* If the user requested for extrapolation, add unphysical points to the left of the
  physical limit */
  
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {

    /* We space the extra points linearly with the same step used for the physical regime */
    double k_step = k_step_max;
    
    /* Uncomment to use the same step separating the first two points in the sources, instead.
    This step can be very small, so we use an extra parameter to increase it. */
    // double n = 1;
    // double first_physical_step = ppt2->k3[index_k1][index_k2][1] - ppt2->k3[index_k1][index_k2][0];
    // k_step = MIN (k_step_max, n*first_physical_step);
    
    class_test (k_step < ppr->smallest_allowed_variation,
      ptr2->error_message,
      "stopping to avoid segmentation fault, step=%g", k_step);

    while (k < k_min_pt) {
      k += k_step;
      index_k_tr++;
    }

    if (k != k_min_pt)
      k = k_min_pt;
  }
  
  index_k_pt++;

  /* The regime where the triangular condition is satisfied starts here */
  ptr2->k_physical_start_k1k2[index_k1][index_k2] = index_k_tr-1;

  /* Take the next points from the sources sampling as long as the step is small enough.
  This assumes that ppt->k3 is a growing array. */
  while ((index_k_pt < k_pt_size) && ((ppt2->k3[index_k1][index_k2][index_k_pt] - k) < k_step_max)) {
    k = ppt2->k3[index_k1][index_k2][index_k_pt];
    index_k_pt++;
    index_k_tr++;
  }

  /* All the points that we added in the above while-loop satisfy the triangular condition */
  ptr2->k_physical_size_k1k2[index_k1][index_k2] = index_k_pt;

  /* Then, points spaced linearily with step k_step_max. Note that if the user asked for k3-extrapolation,
  here we also add points beyond the physical limit dictated by the triangular condition. */
  while (k < k_max_tr) {

    k += k_step_max;
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
    "found less than 2 valid values for the k-grid of modes k1=%f, k2=%f. This should not happen\
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
 * To understand what this function does, please refer to the documentation of
 * transfer2_get_k3_size(), of which this function is largely a copy.
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
                   the k3 sampling-points of the transfer function in (k1,k2) */
      int * last_used_index_pt /**< output, last point in ppt2->k3 that appears in the
                               transfer function k3 sampling; for debug purposes only. */
      )
{

  /* Limits of k3 from the perturbations2.c module */

  int k_pt_size = ppt2->k3_size[index_k1][index_k2];  
  double k_min_pt = ppt2->k3[index_k1][index_k2][0];
  double k_max_pt = ppt2->k3[index_k1][index_k2][k_pt_size-1];

  /* The limits for the current k3-grid were already computed in transfer2_get_k3_sizes() */

  int k_tr_size = ptr2->k_size_k1k2[index_k1][index_k2];
  double k_min_tr = ptr2->k_min_k1k2[index_k1][index_k2];
  double k_max_tr = ptr2->k_max_k1k2[index_k1][index_k2];
  
  /* Find the linear step to use for the k3 sampling */

  double k_step_max;

  if (ptr2->k_sampling == class_transfer2_k3_sampling) {

    k_step_max = 2.*_PI_/(ptr2->tau0-ptr2->tau_rec)*ppr2->q_linstep_song;
    
    if ((ppr->load_run == _TRUE_) && (ppr2->old_run == _TRUE_))
      k_step_max = 2.*_PI_/ptr2->rs_rec*ppr2->q_linstep_song;
  }

  else if (ptr2->k_sampling == bessel_k3_sampling) {
    k_step_max = pbs2->xx_step/ptr2->tau0;
  }

  class_test(k_step_max == 0,
    ptr2->error_message,
    "stop to avoid infinite loop");

  /* First point in the k-sampling */
    
  int index_k_pt = 0;
  int index_k_tr = 0;
  double k = k_min_tr;
  k3[0] = k_min_tr;
  index_k_tr++;

  /* Add points to the left */
  
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {

    double k_step = k_step_max;

    /* Uncomment to use the same step separating the first two points in the sources, instead.
    This step can be very small, so we use an extra parameter to increase it. */
    // double n = 1;
    // double first_physical_step = ppt2->k3[index_k1][index_k2][1] - ppt2->k3[index_k1][index_k2][0];
    // k_step = MIN (k_step_max, n*first_physical_step);
    
    class_test (k_step < ppr->smallest_allowed_variation,
      ptr2->error_message,
      "stopping to avoid segmentation fault, step=%g", k_step);

    while (k < k_min_pt) {
      k += k_step;
      k3[index_k_tr] = k;
      index_k_tr++;    
    }

    if (k != k_min_pt) {
      k = k_min_pt;
      k3[index_k_tr - 1] = k;
    }
  }
  
  index_k_pt++;

  /* Take the next points from the sources sampling as long as the step is small enough */
  while ((index_k_pt < k_pt_size) && ((ppt2->k3[index_k1][index_k2][index_k_pt] - k) < k_step_max)) {
      k = ppt2->k3[index_k1][index_k2][index_k_pt];
      k3[index_k_tr] = k;
      index_k_pt++;
      index_k_tr++;
  }

  *last_used_index_pt = index_k_pt;

  /* Then, points spaced linearily with step k_step_max */
  while ((index_k_tr < k_tr_size) && (k < k_max_tr)) {
    k += k_step_max;
    k3[index_k_tr] = k;
    index_k_tr++;
  }

  /* Consistency check on the maximum value of k3 */
  class_test (k3[k_tr_size-1] > k_max_tr,
    ptr2->error_message,
    "bug in k list calculation, k_max_transfers2=%.17f is larger than k_max_tr=%.17f,\
 should never happen. This can happen for a race condition.",
    k3[k_tr_size-1], k_max_tr);

  /* Check that the number of k3 points matches the precomputed one */
  class_test (index_k_tr != k_tr_size,
    ptr2->error_message,
    "mismatch between k3 grid sizes: index_k_tr-1=%d, k_tr_size=%d",
    index_k_tr-1, k_tr_size);

  /* Debug - Print the k3 list */
  // int index_k1_debug = index_k1;
  // int index_k2_debug = index_k2;
  //
  // if ((index_k1==index_k1_debug) && (index_k2==index_k2_debug)) {
  //
  //   fprintf (stderr, "# ~~~~ (index_k1,index_k2)=(%d,%d), (k1,k2)=(%g,%g), k_tr_size=%d ~~~~~\n",
  //     index_k1, index_k2, ppt2->k[index_k1], ppt2->k[index_k2], ptr2->k_size_k1k2[index_k1][index_k2]);
  //
  //   int first_k_phys = ptr2->k_physical_start_k1k2[index_k1][index_k2];
  //   int last_k_phys = ptr2->k_physical_start_k1k2[index_k1][index_k2]
  //     + ptr2->k_physical_size_k1k2[index_k1][index_k2] - 1;
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

#ifdef DEBUG

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

#endif // DEBUG

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

  /* Check that the sources were computed for enough k3-values to allow a meaningful interpolation */
  for(int index_k1=0; index_k1<ppt2->k_size; ++index_k1)
    for(int index_k2=0; index_k2<index_k1+1; ++index_k2)
      class_test (
        ((ppt2->k3_size[index_k1][index_k2] < 2) && (ppr2->sources_k3_interpolation==linear_interpolation))
        || ((ppt2->k3_size[index_k1][index_k2] < 4) && (ppr2->sources_k3_interpolation==cubic_interpolation)),
        ptr2->error_message,
        "index_k1=%d, index_k2=%d: cannot do interpolation in k3 with just k3_size=%d values. Increase k3_size.",
        index_k1, index_k2, ppt2->k3_size[index_k1][index_k2]);
  
  /* Allocate k1 level of arrays */
  int k1_size = ppt2->k_size;
  class_alloc(ptr2->k_size_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_physical_start_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_physical_size_k1k2, k1_size*sizeof(int *), ptr2->error_message);
  class_alloc(ptr2->k_min_k1k2, k1_size*sizeof(double *), ptr2->error_message);
  class_alloc(ptr2->k_max_k1k2, k1_size*sizeof(double *), ptr2->error_message);

  for(int index_k1=0; index_k1<k1_size; ++index_k1) {
  
    double k1 = ppt2->k[index_k1];
  
    /* Allocate k2 level of arrays */
    int k2_size = ppt2->k_size;
    class_alloc(ptr2->k_size_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_physical_start_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_physical_size_k1k2[index_k1], k2_size*sizeof(int), ptr2->error_message);
    class_alloc(ptr2->k_min_k1k2[index_k1], k2_size*sizeof(double), ptr2->error_message);
    class_alloc(ptr2->k_max_k1k2[index_k1], k2_size*sizeof(double), ptr2->error_message);
  
    for(int index_k2=0; index_k2<k2_size; ++index_k2) {
  
      double k2 = ppt2->k[index_k2];

      /* Compute the number of points in the k3_grid for the given pair of (k1,k2), and 
      store it in ptr2->k3_size[index_k1][index_k2] */
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
 * Compute the second-order transfer function T_lm(k1,k2,k) for the 
 * transfer type given by index_tt.
 *
 * This function is basically a wrapper to transfer2_integrate(), which
 * solves the line of sight integral by convolving the projection function J
 * (computed in the bessel2.c module) with the source function S
 * (passed with the array interpolated_sources_in_time). Details on the line
 * of sight formalism can be found in Sec. 5.5.
 *
 * The temperature tranfer function is obtained by a single call
 * to transfer2_integrate(), while the polarised transfer functions
 * require a second call in order to compute the extra contribution from
 * the mixing of the E and B-modes. Please refer to Sec. 5.5.1.4 of
 * http://arxiv.org/abs/1405.2280 for further details on the polarised 
 * transfer functions.
 *
 * The output is just a number and will be written both in pw->transfer 
 * and in ptr2->transfer[index_tt][index_k1][index_k2][index_k].
 */

int transfer2_compute (
        struct precision * ppr,
        struct precision2 * ppr2,
        struct perturbs2 * ppt2,
        struct bessels * pbs,
        struct bessels2 * pbs2,
        struct transfers2 * ptr2,
        int index_k1, /**< input, ppt2->k[index_k1] is the k1-value for which we shall compute T_lm(k1,k2,k) */
        int index_k2, /**< input, ppt2->k[index_k2] is the k2-value for which we shall compute T_lm(k1,k2,k) */
        int index_k,  /**< input, pw->k_grid[index_k] is the k-value for which we shall compute T_lm(k1,k2,k) */
        int index_l,  /**< input, ptr2->l[index_l] is the l-value for which we shall compute T_lm(k1,k2,k) */
        int index_m,  /**< input, ptr2->m[index_m] is the m-value for which we shall compute T_lm(k1,k2,k) */
        int index_tt, /**< input, index that identifies which field (T,E,B) we are computing T_lm(k1,k2,k) for */
        double ** interpolated_sources_in_time, /**< input, value of the source functions S_lm(k1,k2,k) at all times in pw->tau_grid */
        struct transfer2_workspace * pw /**< input and output, workspace containing the time grid (pw->tau_grid) and the value of k (pw->index_k)  */
        )
{
    

  if (ptr2->transfer2_verbose > 4)
    printf("     * computing transfer function for (l,m) = (%d,%d)\n", ptr2->l[index_l], ptr2->m[index_m]);


  /* Determine what we are deasling with: temperature, E-modes or B-modes? */
  int transfer_type = ptr2->tt2_start[index_tt];


  // =====================================================================================
  // =                                   Temperature                                     =
  // =====================================================================================
  
  if (ppt2->has_source_T && transfer_type==ptr2->index_tt2_T) {
  
    /* Number of multipole sources to consider */
    pw->L_max = ppr2->l_max_los_t;

    /* The transfer function for photon intensity requires a single evaluation of the
    line of sight integral, as in eq. 5.95 of http://arxiv.org/abs/1405.2280. */
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
                  interpolated_sources_in_time,
                  pbs2->index_J_TT,   /* Temperature projection function */
                  ppt2->index_tp2_T,  /* Temperature source function */
                  pw,
                  &(pw->transfer)),
      ptr2->error_message,
      ptr2->error_message);
      
  } // end of temperature



  // =====================================================================================
  // =                                     E-modes                                       =
  // =====================================================================================
  
  else if (ppt2->has_source_E && transfer_type==ptr2->index_tt2_E) {

    /* Number of multipole sources to consider */
    pw->L_max = ppr2->l_max_los_p;

    /* E-modes and B-modes mix while they propagate. Therefore, we have both a direct
    contribution (E->E), as for temperature, and a mixing one (B->E) that encodes the
    conversion of B-modes into E-modes as photons propagate. Each contribution has 
    a different source and projection function. See Sec. 5.5.1.4 of
    http://arxiv.org/abs/1405.2280. */
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
                  interpolated_sources_in_time,
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
                    interpolated_sources_in_time,
                    pbs2->index_J_EB,   /* EB mixing projection function, vanishes for m=0 */
                    ppt2->index_tp2_B,  /* B-mode source function, vanishes for m=0 */
                    pw,
                    &(mixing_contribution)),
        ptr2->error_message,
        ptr2->error_message);
      
      }
      
      /* The integral is given by the sum of the E->E and B->E contributions. */
      pw->transfer = direct_contribution + mixing_contribution;
      
  } // end of E-modes



  // =====================================================================================
  // =                                     B-modes                                       =
  // =====================================================================================
  
  /* We consider the B-modes only for non-scalar modes (ppt2->has_source_B is false 
  for m=0). Even if didn't do so, both the direct and mixed contributions would vanish.
  The direct contribution vanishes because it involves the B-mode sources at m=0, which
  vanish. The mixed contribution vanishes because the projection function vanishes for
  m=0 (see odd function in eq. 5.104 of http://arxiv.org/abs/1405.2280) */

  else if (ppt2->has_source_B && transfer_type==ptr2->index_tt2_B) {
  
    /* Number of multipole sources to consider */
    pw->L_max = ppr2->l_max_los_p;

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
                  interpolated_sources_in_time,
                  pbs2->index_J_EE,   /* E-mode projection function (same as J_BB) */
                  ppt2->index_tp2_B,  /* B-mode source function, vanishes for m=0 */
                  pw,
                  &(direct_contribution)),
      ptr2->error_message,
      ptr2->error_message);

    /* Mixing contribution (E->B). The projection function for this contribution (J_BE)
    is equal to minus the one for B->E (-J_EB). We shall adjust for this sign below,
    when we sum the direct and mixed contributions. */
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
                  interpolated_sources_in_time,
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

  } // end of B-modes



  // =====================================================================================
  // =                             Rescale transfer functions                            =
  // =====================================================================================

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



  // =====================================================================================
  // =                              Store transfer function                              =
  // =====================================================================================

  /* Store transfer function in transfer structure */
  ptr2->transfer[index_tt][index_k1][index_k2][index_k] = pw->transfer;

  /* Counter to keep track of the number of values fit into ptr2->transfer */
  #pragma omp atomic
  ++ptr2->count_memorised_transfers;


  return _SUCCESS_;

}




/**
 * Solve the line of sight integral for the set of parameters (k1,k2,k,l,m)
 * using the trapezoidal method.
 *
 * The line of sight integral is a convolution over time between the source
 * function S_Lm(k1,k2,k,tau) and the projection function J_Llm(x), with
 * x=k*(tau_0-tau) and L summed from 0 to L_max.
 *
 * The time grid for the integration was computed in the transfer2_get_time_grid()
 * function, and is stored in pw->time_grid.
 *
 * The second-order source function is passed as an array of values 
 * at the grid points, interpolated from the ppt2->sources array; the interpolation
 * was performed in the functions transfer2_interpolate_sources_in_k() and
 * transfer2_interpolate_sources_in_time().
 *
 * The projection function will be interpolated in x=k*(tau_0-tau) from the precomputed
 * table stored in the bessel2 structure.
 */
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
  
  /* Initialise the output value of the transfer function */
  *integral = 0;

#ifdef DEBUG
  /* Check that we computed the projection functions J(x) for all the needed values of x.
  Since x appears in the line of sight integral as x=k*(tau0-tau), here we check that x_max
  is larger than k*(tau0-tau_min). */
  class_test (k*pw->tau0_minus_tau[0] > pbs2->xx_max,
    ptr2->error_message,
    "not enough J's.  Increase l_max to %g or decrease k_max to %.3g.",
    ceil(k*pw->tau0_minus_tau[0]/ppr->k_max_tau0_over_l_max),
    pbs2->xx_max/pw->tau0_minus_tau[0]);
#endif // DEBUG
    


  // =====================================================================================
  // =                             Perform the integration                               =
  // =====================================================================================
  
  /* Solve the line of sight integral by looping over time. The integral is in eq. 5.95
  of http://arxiv.org/abs/1405.2280. The time grid is built in transfer2_get_time_grid()
  to match the sampling of the line of sight sources (ppt2->tau_sampling), with the addition
  of extra points to follow the oscillations of the projection functions J. We implement
  the integration following the positive direction of time, which is the negative direction
  of x=k(tau0-tau). */
  
  for (int index_tau=0; index_tau < pw->tau_grid_size; ++index_tau) {
  
    /* Argument of the projection function at the considered time */
    double x = k * pw->tau0_minus_tau[index_tau];
      
    /* Position of x inside pbs2->xx, the array that we used to sample the projection
    functions */
    int index_x = (int)(x/pbs2->xx_step);
  
    /* Interpolation weight assigned to x */
    double a_J = (pbs2->xx[index_x+1] - x)/pbs2->xx_step;
    
    /* The integrand function is the sum over L of the product between the source
    S_Lm(k1,k2,k,tau) and the projection functions J_Llm(k(tau0-tau)); therefore,
    we shall increment the integrand function inside a loop over L from 0 to 
    pw->L_max. In principle, pw->L_max should be of order O(2000), but in practice
    it is O(few) at recombination due to tight-coupling suppression of higher-order
    multipoles. */
    double integrand = 0;

    for(int index_L=0; index_L<=pw->L_max; ++index_L) {
  
      /* The 3j symbol in the definition of J forces the azimuthal number m to be smaller
      than both l and L (see eq. 5.97 of http://arxiv.org/abs/1405.2280) */
      int L = pbs2->L[index_L];
      if (abs(m) > MIN(L,l))
        continue;
            
      /* The projection function J_Llm is negligible when x is too small; index_x_min is the first
      value in pbs2->xx where J_Llm(x) is non-negligible */
      int index_x_min = pbs2->index_xmin_J[index_J][index_L][index_l][index_m];

      /* Skip the contribution to the integral from this L if the projection function is
      negligible. The projection function J_Llm(x) is basically a Bessel function of order l,
      so it vanishes when x<<l. Since its argument is x=k*(tau_0-tau), this means that J
      vanishes when tau is much larger than tau_0-l/k; that is, high-l multipoles do not get
      contributions from low-k modes. If you do not include this check, you will get random
      segmentation faults, because you would end up addressing J_Llm_x with a negative x index.
      Note that this check also ensures that we skip L<2 configurations for polarisation,
      as it should be. */
      if (index_x < index_x_min)
        continue;
  
      /* Index needed to address the x-level of the projection function array, pbs2->J_Llm_x.
      The above check ensures that the index is non-negative. */
      int index_x_in_J = index_x - index_x_min;
      
#ifdef DEBUG
      /* Check that index_x is within the limits for which we have computed the projection functions */
      int x_size = pbs2->x_size_J[index_J][index_L][index_l][index_m];
      class_test ((index_x_in_J >= x_size) || (index_x_in_J < 0),
        ptr2->error_message,
        "L=%d, l=%d, m=%d, x=%g: wrong indexing for pbs2->J_Llm_x; index_x_in_J=%d, x_size=%d, index_x_min=%d",
        L, l, m, x, index_x_in_J, x_size, index_x_min);
#endif // DEBUG
  
  
      // ---------------------------------------------------------------------------
      // -                              Interpolate J(x)                           -
      // ---------------------------------------------------------------------------

      /* Interpolate J in x=k*(tau0-tau) */
      double J_Llm;
      double * J = &(pbs2->J_Llm_x[index_J][index_L][index_l][index_m][index_x_in_J]);
      double J_Llm_left = *J;
      double J_Llm_right = *(J+1);
      
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

      /* Debug - Check the interpolation of the projection function */
      // if ((index_tau == ppt2->index_tau_rec) && (pw->index_k1==1) && (pw->index_k2==0)) {
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


      // ---------------------------------------------------------------------------
      // -                           Build the integrand                           -
      // ---------------------------------------------------------------------------

      /* Pre-interpolated source function in tau and k3 for the desired source type */
      double source_Lm = interpolated_sources_in_time[index_source_monopole + lm(L,m)][index_tau];

      /* Increment the integrand function by adding another L-multipole  */
      integrand += J_Llm * source_Lm;
  
      /* Debug - Print the integrand function for a given set of (l,m,k1,k2,k) */
      // if ( (l==100) && (m==0) ) {
      //   if ( (pw->index_k1==0) && (pw->index_k2==1) && (pw->index_k==2500) ) {
      //     fprintf (stderr, "%15f %15f\n", ptr2->tau0_minus_tau[index_tau], integrand);
      //   }
      // }

    } // end of for(index_L)    
    
    /* Increment the result with the contribution from the considered time-step */
    *integral += integrand * pw->delta_tau[index_tau];
  
    /* Debug - Output the integrand as a function of time */
    // if (l == 2) {
    //   if (index_tau==0)
    //     fprintf_k_debug (stderr, "# k1=%g, k2=%g, k=%g\n", ppt2->k[index_k1], ppt2->k[index_k2], k);
    //   p2 (pw->tau_grid[index_tau], integrand);
    // }
  
  } // end of for(index_tau)

  /* Correct for factor 1/2 from the trapezoidal rule */
  *integral *= 0.5;


#ifdef DEBUG
  /* Test for nans */
  class_test (isnan(*integral),
    ptr2->error_message,
    "found nan in second-order transfer function");
#endif // DEBUG


  return _SUCCESS_;
}



/**
 * Determine the integration grid in time for the line of sight integral.
 * 
 * The line-of-sight integral convolves the projection functions J(k*(tau0-tau))
 * with the source functions S(k1,k2,k,tau). The integration grid in time should be
 * fine enough to catch the variations of both S and J.
 *
 * Our strategy is to start from the time sampling of the source functions
 * (ppt2->tau_sampling) and to add points where it is not dense enough
 * to follow the oscillations in the projection functions.
 *
 * This function will fill the following arrays in the transfer workspace:
 *
 * - pw->tau_grid[index_tau]
 * - pw->tau_grid_size
 * - pw->index_tau_left[index_tau]
 * - tau0_minus_tau[index_tau]
 * - pw->delta_tau[index_tau]
 *
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
  
  double k = pw->k;
  
  /* The limits for the time-integration-grid for the line-of-sight integral are completely
  determined by the sources */
  int tau_size_pt = ppt2->tau_size;
  double tau_min_pt = ppt2->tau_sampling[0];
  double tau_max_pt = ppt2->tau_sampling[tau_size_pt-1];


  // ====================================================================================
  // =                             Determine linear step                                =
  // ====================================================================================
  
  /* We are going to keep all points in the sources sampling (ppt2->tau_sampling), and
  complement them with points linearly sampled according to tau_step_max. Therefore,
  tau_step_max should reflect the oscillation period of the projection functions, rather
  than that of the source functions. */
  double tau_step_max;

  /* In the sources time sampling, we take the time sampling of the line of sight source
  as our integration grid. In this way, we will have a sparser time sampling than with
  the other options, but the execution will be faster because the sources will not need to
  be interpolated. */

  if (ptr2->tau_sampling == sources_tau_sampling) {
    
      pw->tau_grid_size = tau_size_pt;
      
      for (int index_tau=0; index_tau < pw->tau_grid_size; ++index_tau)
        pw->tau_grid[index_tau] = ppt2->tau_sampling[index_tau];
      
      goto fill_delta_tau;
  }

  /* In the custom time sampling, we reason in terms of oscillations. In the line of sight
  integral, the Bessel functions appear with argument k*(tau_0-tau), meaning that they oscillate
  in time with a frequency of roughly 2*pi/k. Here we choose a linear step in tau which is
  proportional to such frequency, via the parameter tau_linstep_song. */

  else if (ptr2->tau_sampling == custom_tau_sampling) {

    tau_step_max = 2*_PI_/k*ppr2->tau_linstep_song;

  }

  /* In the Bessel-type sampling, we require the maximum step between two points to be
  determined by the sampling of the projection functions J, which in turn is controlled
  by the user via the parameter bessel_x_step_song. */

  else if (ptr2->tau_sampling == bessel_tau_sampling) {
    
    tau_step_max = pbs2->xx_step/k;

  }
    
  class_test (tau_step_max == 0,
    ptr2->error_message,
    "stopping to prevent segmentation fault");


  // ====================================================================================
  // =                              Determine grid size                                 =
  // ====================================================================================
  
  int index_tau_tr=0;
  int index_tau_pt=0;

  /* The first point in the grid is the first point in the time sampling of the sources */
  double tau = tau_min_pt;
  index_tau_tr++;
  index_tau_pt++;
  
  /* The following points are taken from the source sampling. If the time step in such sampling
  is larger than tau_step_max, then use tau_step_max as a step, but never throw away the points
  in the source sampling. */
  while (index_tau_pt < tau_size_pt) {
    
    while ((ppt2->tau_sampling[index_tau_pt] - tau) > tau_step_max) {
      tau += tau_step_max;
      index_tau_tr++;
    }
    
    tau = ppt2->tau_sampling[index_tau_pt];  
    ++index_tau_tr;
    ++index_tau_pt;    
  }

  pw->tau_grid_size = index_tau_tr;


  // ====================================================================================
  // =                                 Fill time grid                                   =
  // ====================================================================================
  
  /* - Repeat exactly the same steps as above, but this time filling the time arrays */
  
  index_tau_tr=0;
  index_tau_pt=0;
  
  /* First point in the grid */
  tau = tau_min_pt;
  pw->tau_grid[0] = tau;
  index_tau_tr++;
  index_tau_pt++;

  /* We shall take note of the position of each grid–point inside ppt2->tau_sampling,
  in view of table look-up for interpolating the sources. The first grid point
  coincides with the first point in ppt2->tau_sampling, so we give it index 0. */
  pw->index_tau_left[0] = 0;

  while (index_tau_pt < tau_size_pt) {

    while ((ppt2->tau_sampling[index_tau_pt] - tau) > tau_step_max) {

      tau += tau_step_max;
      pw->tau_grid[index_tau_tr] = tau;
      pw->index_tau_left[index_tau_tr] = index_tau_pt;
      index_tau_tr++;

      /* Debug - Show the time points as we add them to the integration grid */
      // printf("* k=%g,index_k=%d,tau=%.7f,index_tau_pt=%d/%d,index_tau_tr=%d,index_tau_left=%d,delta_tau_pt=%g,tau_step_max=%g\n",
      //   k, index_k, tau, index_tau_pt, tau_size_pt, index_tau_tr, pw->index_tau_left[index_tau_tr],
      //   ppt2->tau_sampling[index_tau_pt] - tau, tau_step_max);
    }
    
    tau = ppt2->tau_sampling[index_tau_pt];
    pw->tau_grid[index_tau_tr] = tau;
    pw->index_tau_left[index_tau_tr] = index_tau_pt;
    index_tau_tr++;
    index_tau_pt++;    

    // printf("k=%g,index_k=%d,tau=%.7f,index_tau_pt=%dof%d,index_tau_tr=%d,index_tau_left=%d,delta_tau_pt=%g,tau_step_max=%g\n",
    //   k, index_k, tau, index_tau_pt, tau_size_pt, index_tau_tr, pw->index_tau_left[index_tau_tr],
    //   ppt2->tau_sampling[index_tau_pt] - tau, tau_step_max);
  }
    
    
  /* The left index cannot be the last index (otherwise you'd get segmentation faults)  */
  for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size; ++index_tau_tr)
    if (pw->index_tau_left[index_tau_tr] == (ppt2->tau_size-1))
      pw->index_tau_left[index_tau_tr] = ppt2->tau_size-2;  
  
  /* Print some info */
  pw->cosk1k2 = (pw->k*pw->k - pw->k1*pw->k1 - pw->k2*pw->k2)/(2*pw->k1*pw->k2);
  if (ptr2->transfer2_verbose > 4)
    printf("     * (k,cosk1k2)=(%.3g,%.3g): the time-integration grid comprises %d+%d points from %g to %g\n",
      pw->k, pw->cosk1k2, ppt2->tau_size, pw->tau_grid_size-ppt2->tau_size, pw->tau_grid[0], pw->tau_grid[pw->tau_grid_size-1]);
  
  /* Debug - Print the integration grid */
  // int index_k_debug = 0, index_k1_debug = 0, index_k2_debug = 1;
  // 
  // fprintf (stderr, "# index_k=%d, k=%g, tau_step_max=%g, tau_min_pt=%g, tau0-tau_min_pt=%g\n",
  //   pw->index_k, pw->k, tau_step_max, tau_min_pt, ptr2->tau0-tau_min_pt);
  // 
  // if ((pw->index_k == index_k_debug) && (pw->index_k1 == index_k1_debug) && (pw->index_k2 == index_k2_debug))
  //   for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size-1; ++index_tau_tr)
  //     fprintf (stderr, "%6d %15.7g %10d\n",
  //       index_tau_tr, pw->tau_grid[index_tau_tr], pw->index_tau_left[index_tau_tr]);

  
#ifdef DEBUG

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

#endif // DEBUG



  // ====================================================================================
  // =                                 Fill delta_tau                                   =
  // ====================================================================================

  fill_delta_tau:

  /* Fill tau0_minus_tau, the argument of the projection function in the line
  of sight integral */
  for (index_tau_tr=0; index_tau_tr < pw->tau_grid_size; ++index_tau_tr) 
    pw->tau0_minus_tau[index_tau_tr] = ptr2->tau0 - pw->tau_grid[index_tau_tr]; 

  /* Fill delta_tau, the measure for the trapezoidal integral */
  pw->delta_tau[0] = pw->tau_grid[1]-pw->tau_grid[0];
      
  for (index_tau_tr=1; index_tau_tr < pw->tau_grid_size-1; ++index_tau_tr)
    pw->delta_tau[index_tau_tr] = pw->tau_grid[index_tau_tr+1]-pw->tau_grid[index_tau_tr-1];
      
  pw->delta_tau[pw->tau_grid_size-1] = pw->tau_grid[pw->tau_grid_size-1]-pw->tau_grid[pw->tau_grid_size-2];


  return _SUCCESS_;

} // end of transfer2_get_time_grid








/**
 * Interpolate the source function S(k1, k2, k, tau) at the desired values of k, for
 * all values of time.
 *
 * This function takes a source function for fixed values of k1 and k2, and returns
 * an array containing the same source function interpolated in the desired values
 * of k, for all times where the source function is sampled.
 *
 * The output arrays, sources_k_spline and interpolated_sources_in_k, are matrices of
 * size k3_size*tau_size, where k3_size = ptr2->k3_size[index_k1][index_k2] and
 * tau_size = ppt2->tau_sampling. They are addressed as:
 * 
 *   sources_k_spline[index_k3*ppt2->tau_size + index_tau]
 *   interpolated_sources_in_k[index_k3*ppt2->tau_size + index_tau]
 * 
 * and must be preallocated.
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
      int index_tp2, /**< input, type of the requested source function */
      double * k_grid, /**< input, list of k-values where to interpolate the source function */
      double * sources_k_spline, /**< output, second derivative of the source function with respect to k, for all time values */
      double * interpolated_sources_in_k /**< output, source function at the desired k-values, for all time values */
      )
{


  /* Shortcuts */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  double * k_pt = ppt2->k3[index_k1][index_k2];
  int k_tr_size = ptr2->k_size_k1k2[index_k1][index_k2];
  double * k_tr = k_grid;
    
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


  // ====================================================================================
  // =                                   Interpolation                                  =
  // ====================================================================================

  /* Limits for which we shall interpolate the sources */
  int physical_size = ptr2->k_physical_size_k1k2[index_k1][index_k2];
  int first_physical_index = ptr2->k_physical_start_k1k2[index_k1][index_k2];
  int last_physical_index = first_physical_index + physical_size - 1;

  /* Interpolate the source function at each k value contained in k_grid, using the usual
  spline interpolation algorithm */
  int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];

  /* Debug - Print the sources as a function of time */
  // index_K = 50;
  // if ((index_k1 == 1) && (index_k2 == 0))
  //   if (index_tp2 == (ppt2->index_tp2_T + lm(2,0)))
  //     for (int index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
  //         printf ("%12g %12g\n", ppt2->tau_sampling[index_tau], sources(index_tau,index_K));

  /* Debug - Print the sources as a function of k3 */
  // int index_tau = ppt2->index_tau_rec;
  // if ((index_k1 == 1) && (index_k2 == 0))
  //   if (index_tp2 == (ppt2->index_tp2_T + lm(1,0)))
  //     for (int index_k=0; index_k < ppt2->k3_size[index_k1][index_k2]; ++index_k)
  //       printf ("%12g %12g\n", ppt2->k3[index_k1][index_k2][index_k], sources(index_tau,index_k));


  for (int index_k_tr = first_physical_index; index_k_tr <= last_physical_index; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0, ptr2->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1-b;

    /* We shall interpolate for each value of conformal time, hence the loop
    on index_tau */
    if (ppr2->sources_k3_interpolation == linear_interpolation) {
      for (int index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1);
    }
    else if (ppr2->sources_k3_interpolation == cubic_interpolation) {
      for (int index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1)
          + ((a*a*a-a) * sources_k_spline[index_tau*k_pt_size + index_k]
          +(b*b*b-b) * sources_k_spline[index_tau*k_pt_size + index_k+1])*h*h/6.0;
    }

  } // end of for (index_k_tr)


  /* Debug - Print the sources as a function of the new grid */
  // int index_tau = ppt2->index_tau_rec;
  // if ((index_k1 == 1) && (index_k2 == 0))
  //   if (index_tp2 == (ppt2->index_tp2_T + lm(1,0)))
  //     for (int index_k_tr = first_physical_index; index_k_tr <= last_physical_index; ++index_k_tr)
  //       printf ("%12g %12g\n", k_tr[index_k_tr], interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau]);


  // ====================================================================================
  // =                                   Extrapolation                                  =
  // ====================================================================================

  /* In the perturbations2.c module, we have computed the source functions assuming the
  triangular condition: |k1-k2|<=k3<=k1+k2. In the transfer2.c module, we need to 
  relax this assumption in order to stabilise numerically the bispectrum integral. This
  extrapolation requires adding extra non-physical points at the left of k_min_pt and
  at the right of k_max_pt. For these points, we assume that the source functions are
  constant and equal to the closest physical point. In this way, the resulting transfer 
  function will only pick the oscillations coming from the convolution with the Bessel
  projection function.
  
  Note that for non scalar modes (m!=0), the transfer functions vanish on the edges of
  the triangular limits (k3=k1+k2 or k3=|k1-k2|), because these correspond
  to configurations where all three wavemodes are aligned (sin(theta_1)=sin(theta_2)=0).
  Hence, we would end up extrapolating zero. We have fixed this by analytically including
  a sin(theta)^(-m) factor in the perturbations2.c module in all the sources, so that the
  amplitude at the borders of the k3 range is non-vanishing. This is also  the factor needed
  to perform the bispetrum integration in bispectra2.c. */
    
  if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {
    
    /* Extrapolation on the right (k > k_max_pt) */
    for (int index_k_tr = last_physical_index+1; index_k_tr < k_tr_size; ++index_k_tr) {
      
      if (ppr->bispectra_k3_extrapolation == flat_k3_extrapolation)
        for (int index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
          interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
            interpolated_sources_in_k[last_physical_index*ppt2->tau_size + index_tau];
    }
    /* Extrapolation on the left (k < k_max_pt) */
    for (int index_k_tr = 0; index_k_tr < first_physical_index; ++index_k_tr) {
      
      if (ppr->bispectra_k3_extrapolation == flat_k3_extrapolation)
        for (int index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
          interpolated_sources_in_k[index_k_tr*ppt2->tau_size + index_tau] = 
            interpolated_sources_in_k[first_physical_index*ppt2->tau_size + index_tau];
    }
  }
  

  return _SUCCESS_;
  
} // end of transfer_interpolate_sources_in_k



/**
 * Interpolate the source function S(k1, k2, k, tau) at the desired values of time.
 *
 * This function takes a source function for fixed values of k1, k2 and k, and returns
 * an array containing the same source function interpolated in the desired values
 * of time.
 *
 * The output arrays, sources_time_spline and interpolated_sources_in_time, are vectors
 * of size pw->tau_grid_size; they are addressed as:
 * 
 *   sources_time_spline[index_tau]
 *   interpolated_sources_in_time[index_tau]
 * 
 * and must be preallocated.
 */

int transfer2_interpolate_sources_in_time (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2,
      struct transfers2 * ptr2,
      int index_tp2, /**< input, type of the requested source function */
      double * interpolated_sources_in_k, /**< input, source function at the desired k-values, for all time values;
                                          this is the output of transfer2_interpolate_sources_in_k(). */
      double * sources_time_spline, /**< output, second derivative of the source function with respect to time */
      double * interpolated_sources_in_time, /**< output, source function at the desired time-value */
      struct transfer2_workspace * pw /**< input, workspace containing the grid of required times (pw->tau_grid) and the value of k (pw->index_k) */
      )
{

  /* Shortcuts */
  int tau_size_pt = ppt2->tau_size;
  double * tau_pt = ppt2->tau_sampling;
  int tau_size_tr = pw->tau_grid_size;
  double * tau_tr = pw->tau_grid;
    
  /* If the integration grid matches the time sampling of the sources, there is no need for
  interpolation */
  if (ptr2->tau_sampling == sources_tau_sampling) {
    
    for (int index_tau = 0; index_tau < pw->tau_grid_size; ++index_tau)     
      interpolated_sources_in_time[index_tau]
        = interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau];
    
  }
  
  else {

    /* Find second derivative of original sources with respect to k in view of spline
    interpolation */
    if (ppr2->sources_time_interpolation == cubic_interpolation) {

      class_call (array_spline_table_columns (
                    tau_pt,
                    tau_size_pt,
                    interpolated_sources_in_k + pw->index_k * tau_size_pt, /* start from index_k */
                    1,  /* we interpolate in time for only 1 k-value (pw->index_k) */
                    sources_time_spline,
                    _SPLINE_EST_DERIV_,
                    ptr2->error_message),
           ptr2->error_message,
           ptr2->error_message);
    }
  
    /* Interpolate the source function at each tau value contained in pw->tau_grid, using the
    usual spline interpolation algorithm */

    for (int index_tau_tr = 0; index_tau_tr < tau_size_tr; ++index_tau_tr) {
    
      /* Determine the index to the left of tau_tr in the sources time sampling, using
      the array pw->index_tau_left, precomputed  in transfer2_get_time_grid() */
      int index_tau = pw->index_tau_left[index_tau_tr];
      double h = tau_pt[index_tau+1] - tau_pt[index_tau];

      class_test(h==0., ptr2->error_message, "stop to avoid division by zero");
    
      double b = (tau_tr[index_tau_tr] - tau_pt[index_tau])/h;
      double a = 1-b;

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

      /* Debug - Print information for the current time of interpolation */
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
    
  } // end of if not sources time sampling

  
  /* Debug - Print the source function at the node points together with its interpolated value */
  // if ((pw->index_k1==1) && (pw->index_k2==0) && (pw->index_k==50)) {
  //
  //   if (index_tp2 == (ppt2->index_tp2_T + lm(2,0))) {
  //
  //     fprintf (stderr, "\n\n");
  //
  //     for (int index_tau=0; index_tau < tau_size_pt; ++index_tau)
  //       fprintf (stderr, "%17.7g %17.7g\n", tau_pt[index_tau], interpolated_sources_in_k[pw->index_k*tau_size_pt + index_tau]);
  //
  //     fprintf (stderr, "\n");
  //
  //     for (int index_tau_tr = 0; index_tau_tr < tau_size_tr; ++index_tau_tr)
  //       fprintf (stderr, "%17.7g %17.7g\n", tau_tr[index_tau_tr], interpolated_sources_in_time[index_tau_tr]);
  //
  //     fprintf (stderr, "\n\n");
  //
  //   }
  // }

  
  return _SUCCESS_;
  
} // end of transfer_interpolate_sources_in_time





/**
 * Save the transfer function T^X_lm(k1,k2,k3) to disk for a given k1 value.
 *
 * This function will write on as many files as the number of transfer types
 * (ptr2->tt2_size) the transfer function value in k1.
 *
 * See the documentation in perturbations2.h (\ref StorageFiles) for more
 * details.
 */

int transfer2_store(
        struct perturbs2 * ppt2,
        struct transfers2 * ptr2,
        int index_k1
        )
{

  /* Print some info */
  if (ptr2->transfer2_verbose > 1)
    printf("     \\ writing transfer function for index_k1=%d ...\n", index_k1);

  /* Since the next block will write on many files, we enclose it in a 
  critical openmp directive so that it is executed serially */
  
  int abort = _FALSE_;

  #pragma omp critical
  {

    for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    
      /* Open file */
      class_open_parallel (ptr2->storage_files[index_tt],
        ptr2->storage_paths[index_tt],
        "ab", ptr2->error_message);
    
      /* Print some info */
      if (ptr2->transfer2_verbose > 3)
        printf("     * writing transfer function for (index_tt,index_k1)=(%d,%d) on '%s' ...\n",
          index_tt, index_k1, ptr2->storage_paths[index_tt]);

      for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
  
        /* Write a chunk with all the k-values for this set of (type,k1,k2) */
        fwrite(
              ptr2->transfer[index_tt][index_k1][index_k2],
              sizeof(double),
              ptr2->k_size_k1k2[index_k1][index_k2],
              ptr2->storage_files[index_tt]
              );

      }

      /* Close file */
      fclose (ptr2->storage_files[index_tt]);
    
    }
    
  } // parallel region

  if (abort)
    return _FAILURE_;

  return _SUCCESS_; 
  
}
  
  
  
  
  
/**
 * Allocate all levels beyond the k1 level of the transfer functions array.
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

  long int count=0;

  for(int index_tt=0; index_tt<ptr2->tt2_size; ++index_tt) {

    int k1_size = ppt2->k_size;

    /* Allocate k2 level.  Note that, as for ppt2->sources, the size of this level is smaller
    than the k1 level, and it depends on k1. The reason is that we only need to compute
    the transfer functions for those k2's that are smaller than k1, because our equations
    are symmetrised with respect to k1<->k2. */
    int k2_size = index_k1 + 1;
  
    class_alloc(
      ptr2->transfer[index_tt][index_k1],
      k2_size * sizeof(double *),
      ptr2->error_message);
  
    for(int index_k2=0; index_k2<=index_k1; ++index_k2) {      

      /* Allocate k level */
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

  return _SUCCESS_;

}

/**
 * Free all levels beyond the type level of the transfer functions array.
 */

int transfer2_free_type_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_tt
     )
{

  /* Free memory only if needed */
  if (!ptr2->transfers_allocated[index_tt])
    return _SUCCESS_;

  for (int index_k1=0; index_k1 < ppt2->k_size; ++index_k1) {

    for (int index_k2=0; index_k2 <= index_k1; ++index_k2)
      free(ptr2->transfer[index_tt][index_k1][index_k2]);

    free(ptr2->transfer[index_tt][index_k1]);

  }
  
  free(ptr2->transfer[index_tt]);

  ptr2->transfers_allocated[index_tt] = _FALSE_;
  ptr2->transfers_available[index_tt] = _FALSE_;

  return _SUCCESS_;

}




/**
 * Free all levels beyond the k1 level of the transfer functions array.
 */

int transfer2_free_k1_level(
     struct perturbs2 * ppt2,
     struct transfers2 * ptr2,
     int index_k1
     )
{

  int k1_size = ppt2->k_size;

  long int count = 0;

  for (int index_tt = 0; index_tt < ptr2->tt2_size; index_tt++) {
    for(int index_k2=0; index_k2<=index_k1; ++index_k2) {
      free(ptr2->transfer[index_tt][index_k1][index_k2]);
      count += ptr2->k_size_k1k2[index_k1][index_k2];
    }

    free(ptr2->transfer[index_tt][index_k1]);
  }

  /* Print some debug information on memory consumption */
  if (ptr2->transfer2_verbose > 2)
    printf("     * freed ~ %.2f MB from ptr2->transfer (%ld doubles) for index_k1=%d\n",
      count*sizeof(double)/1e6, count, index_k1);

  return _SUCCESS_;

}
