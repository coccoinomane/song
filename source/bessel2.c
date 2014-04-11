/** @file bessel2.c Documented Bessel module.
 *
 * Guido W. Pettinari, 17.03.2013
 *
 * This module loads spherical Bessel functions
 * (either read from file or computed from scratch).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# bessel_init() at the beginning (anytime after input_init() and before transfer_init())
 * -# bessel_at_x() at any time for computing a value j_l(x) at any x by interpolation
 * -# bessel_free() at the end
 */

#include "bessel2.h"




/**
 * Get spherical Bessel functions (either read from file or compute
 * from scratch).
 *
 * Each table of spherical Bessel functions \f$ j_l(x) \f$ corresponds
 * to a set of values for: 
 *
 * -# pbs->l[index_l]: list of l values l of size pbs->l_size
 * -# pbs->x_step: step dx for sampling Bessel functions \f$ j_l(x) \f$
 * -# pbs->x_max: last value of x (always a multiple of x_step!)
 * -# pbs->j_cut: value of \f$ j_l \f$ below which it is approximated by zero (in the region x << l)
 *
 * This function checks whether there is alread a file "bessels.dat"
 * with the same l's, x_step, x_max, j_cut). 
 * If yes, it fills the table of bessel functions (and
 * their second derivatives, needed for spline interpolation) with the
 * values read from the file. If not, it computes all values using
 * bessel_j_for_l(), and stores them both in the bessels
 * stucture pbs, and in a file "bessels.dat" (in view of the next
 * runs).
 *
 * @param ppr Input : pointer to precision strucutre
 * @param pbs Output: initialized bessel structure 
 * @return the error status
 */

int bessel2_init(
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2
    )
{

  // *** Do we need to compute the 2nd-order projection functions?
  /* TODO: include a flag so that we don't compute them if we are loading the transfer
  functions from disk */
  if ((pbs->l_max==0) || ((ppt2->has_cls == _FALSE_) && (ppt2->has_bispectra == _FALSE_))
  || ((ppr2->load_transfers_from_disk == _TRUE_) && (ppr->load_bispectra_from_disk == _TRUE_))) {

    if (pbs2->bessels2_verbose > 0)
      printf("No second-order projection functions needed. Second-order Bessel module skipped.\n");

    return _SUCCESS_;
  }
  else {
    if (pbs2->bessels2_verbose > 0)
      printf("Computing 2nd-order projection functions\n");
  }



  // ==============================================================================
  // =                               Preparations                                 =
  // ==============================================================================

  /* Initialize counter for the memory allocated in the projection functions */
  pbs2->count_allocated_Js = 0;

  /* Determine minimum allowed values of the Bessels and of the J's */
  pbs2->j_l1_cut  = ppr2->bessel_j_cut_2nd_order;
  pbs2->J_Llm_cut = ppr2->bessel_J_cut_2nd_order;

  /* Determine the grid in x where J_Llm(x) will be sampled */
  class_call (bessel2_get_xx_list (ppr, ppr2, pbs, pbs2),
    pbs2->error_message,
    pbs2->error_message);
    



  
  // ================================================================================
  // =                          Count projection functions                          =
  // ================================================================================
  
  int index_J = 0;
  
  if (ppt2->has_cmb_temperature == _TRUE_) {

    /* T->T projection function */
    pbs2->has_J_TT = _TRUE_;
    pbs2->index_J_TT = index_J++;
    
  }

  /* We compute the projection functions for both the E and B-modes, regardless of
  which polarisation type is requested, because free streaming makes the two types
  of polarisation mix in the line of sight integral */
  if ((ppt2->has_cmb_polarization_e == _TRUE_) || (ppt2->has_cmb_polarization_b == _TRUE_)) {

    /* E->E projection function (equal to B->B) */
    pbs2->has_J_EE = _TRUE_;
    pbs2->index_J_EE = index_J++;

    /* B->E projection function (equal to minus E->B) */
    pbs2->has_J_EB = _TRUE_;
    pbs2->index_J_EB = index_J++;
    
  }

  pbs2->J_size = index_J;

  if (pbs2->bessels2_verbose > 0)
    printf(" -> will compute size_J=%d projection functions\n", pbs2->J_size);





  // ==========================================================================
  // =                            Compute j_l1(x)                             =
  // ==========================================================================
      
  /* In order to compute the J_Llm_x functions, we need to evaluate the spherical
  Bessel functions j_l1(x) for a wide range of l1's and x's.  The range of l1's
  is larger than the standard pbs2->l_list:  for a given L and l, we will need all
  the l1's between abs(l-L) and l+L.  Here, we precompute j_l1(x) for the extended
  range contained in pbs2->l1.  The process is equivalent to what it is done in
  in the first part of bessel_init, where pbs2->j was filled. */

  /* Get the list where we will need to compute Bessel functions */
  class_call (bessel2_get_l1_list (ppr,ppr2,pbs,pbs2),
    pbs2->error_message,
    pbs2->error_message);


  /* Allocate the l1-level of the arrays */
  class_alloc (pbs2->j_l1,pbs2->l1_size*sizeof(double*),pbs2->error_message);
  if (ppr->bessels_interpolation == cubic_interpolation)
    class_alloc (pbs2->ddj_l1,pbs2->l1_size*sizeof(double*),pbs2->error_message);
  class_alloc (pbs2->index_xmin_l1,pbs2->l1_size*sizeof(int*),pbs2->error_message);
  class_alloc (pbs2->x_size_l1,pbs2->l1_size*sizeof(int*),pbs2->error_message);
  class_alloc (pbs2->x_min_l1,pbs2->l1_size*sizeof(double*),pbs2->error_message);

  if (pbs2->bessels2_verbose > 1)
    printf (" -> computing spherical Bessel functions in ell_1\n");

  /* Beginning of parallel region */
  int abort = _FALSE_;
  #pragma omp parallel shared(ppr,pbs,abort)
  {

    /* Loop over l1 and x values, compute j_l1(x) for each of them */
    int index_l1;
    
    #pragma omp for schedule (dynamic)
    for (index_l1 = 0; index_l1 < pbs2->l1_size; index_l1++) {

      class_call_parallel (bessel2_j_for_l1 (ppr,ppr2,pbs,pbs2,index_l1),
        pbs2->error_message,
        pbs2->error_message);

      #pragma omp flush(abort)

    } // end of for(index_l1)

  } // end of parallel region
  if (abort == _TRUE_) return _FAILURE_;


  /* Debug the computation of the Bessel functions */
  // {
  //   double j;
  //   int index_l = 120;
  //   double x = 233.24142;
  //   bessel2_l1_at_x (pbs2, x, index_l, &j);
  //   printf ("j_%d(%g) = %g\n", pbs2->l1[index_l], x, j);
  //   bessel2_l1_at_x_linear (pbs2, x, index_l, &j);
  //   printf ("j_%d(%g) = %g\n", pbs2->l1[index_l], x, j);
  // }



  // ========================================================================================
  // =                           Compute the projection functions                           =
  // ========================================================================================

  /* Compute the J's, that is the combination of Bessels and 3j-symbols needed to compute
  the line of sight integral at second order. The J's are 3-index objects depending on
  an argument x.  We shall denote them as J_Llm(x). We shall store the projection function
  for temperature (TT), polarisation (EE) and for the mixing of polarisation (EB) */

  /* Create the m array by copying it from ppr2->m */
  pbs2->m_size = ppr2->m_size;
  class_alloc (pbs2->m, pbs2->m_size*sizeof(int), pbs2->error_message);
  for (int index_m=0; index_m<pbs2->m_size; ++index_m)
    pbs2->m[index_m] = ppr2->m[index_m];
  

  /* Create the L array. L is the multipole index that will be summed over in the 
  line of sight integral to obtain the second-order transfer function. */
  pbs2->L_size = pbs2->L_max + 1;
  class_alloc (pbs2->L, pbs2->L_size*sizeof(int), pbs2->error_message);
  for (int index_L=0; index_L<pbs2->L_size; ++index_L)
    pbs2->L[index_L] = index_L;
    

  // ---------------------------------------------------------------------
  // -                     Allocate all levels but x                     -
  // ---------------------------------------------------------------------
  
  /* The following arrays have the (J,L,l,m) levels:
     - pbs2->x_size_J
     - pbs2->index_xmin_J
     - pbs2->x_min_J
     - pbs2->J_Llm_x
     - pbs2->ddJ_Llm_x
  pbs2->J_Llm_x and pbs2->ddJ_Llm_x also have an extra x-level. */

  /* Level of the projection function type */
  class_alloc (pbs2->index_xmin_J,  pbs2->J_size*sizeof(int***), pbs2->error_message);
  class_alloc (pbs2->x_min_J,  pbs2->J_size*sizeof(double***), pbs2->error_message);
  class_alloc (pbs2->x_size_J, pbs2->J_size*sizeof(int***), pbs2->error_message);
  class_alloc (pbs2->J_Llm_x,  pbs2->J_size*sizeof(double****), pbs2->error_message);
  if (ppr->bessels_interpolation == cubic_interpolation)
    class_alloc (pbs2->ddJ_Llm_x,  pbs2->J_size*sizeof(double****), pbs2->error_message);

  for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {

    /* L-level */
    class_alloc (pbs2->index_xmin_J[index_J],  pbs2->L_size*sizeof(int**), pbs2->error_message);
    class_alloc (pbs2->x_min_J[index_J],  pbs2->L_size*sizeof(double**), pbs2->error_message);
    class_alloc (pbs2->x_size_J[index_J], pbs2->L_size*sizeof(int**), pbs2->error_message);
    class_alloc (pbs2->J_Llm_x[index_J],  pbs2->L_size*sizeof(double***), pbs2->error_message);
    if (ppr->bessels_interpolation == cubic_interpolation)
      class_alloc (pbs2->ddJ_Llm_x[index_J],  pbs2->L_size*sizeof(double***), pbs2->error_message);
  
    /* l-level */
    for (int index_L=0; index_L<pbs2->L_size; ++index_L) {
  
      class_alloc (pbs2->index_xmin_J[index_J][index_L], pbs->l_size*sizeof(int*), pbs2->error_message);
      class_alloc (pbs2->x_min_J[index_J][index_L], pbs->l_size*sizeof(double*), pbs2->error_message);
      class_alloc (pbs2->x_size_J[index_J][index_L], pbs->l_size*sizeof(int*), pbs2->error_message);
      class_alloc (pbs2->J_Llm_x[index_J][index_L], pbs->l_size*sizeof(double**), pbs2->error_message);
      if (ppr->bessels_interpolation == cubic_interpolation)
        class_alloc (pbs2->ddJ_Llm_x[index_J][index_L], pbs->l_size*sizeof(double**), pbs2->error_message);
  
      /* m-level */
      for (int index_l=0; index_l<pbs->l_size; ++index_l) {

        /* Only allocate those m's that satisfy the condition m<=MIN(L,l) */
        // int m_size = ppr2->index_m_max[pbs2->L[index_L]] + 1;
        int m_size = MIN (ppr2->index_m_max[pbs2->L[index_L]], ppr2->index_m_max[pbs->l[index_l]]) + 1;

        class_alloc (pbs2->index_xmin_J[index_J][index_L][index_l], m_size*sizeof(int), pbs2->error_message);
        class_alloc (pbs2->x_min_J[index_J][index_L][index_l], m_size*sizeof(double), pbs2->error_message);
        class_alloc (pbs2->x_size_J[index_J][index_L][index_l], m_size*sizeof(int), pbs2->error_message);
        class_alloc (pbs2->J_Llm_x[index_J][index_L][index_l], m_size*sizeof(double*), pbs2->error_message);
        if (ppr->bessels_interpolation == cubic_interpolation)
          class_alloc (pbs2->ddJ_Llm_x[index_J][index_L][index_l], m_size*sizeof(double*), pbs2->error_message);

      } // end of for(index_l)
    } // end of for(index_L)
  } // end of loop on projection functions type


  // ---------------------------------------------------------------------------------
  // -                  Loop on the type of projection functions                     -
  // ---------------------------------------------------------------------------------
  
  for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {
  
    if (pbs2->bessels2_verbose > 1)
      printf (" -> computing the projection function #%d\n", index_J);
  
    /* Loop on the source index */
    for (int index_L = 0; index_L < pbs2->L_size; ++index_L) {
      
      /* Beginning of parallel region */
      abort = _FALSE_;
      #pragma omp parallel shared(index_L,ppr,pbs,abort)
      {
  
        /* Loop on the external l-index */
        #pragma omp for schedule (dynamic)
        for (int index_l = 0; index_l < pbs->l_size; ++index_l) {
        
          /* There are two important consideration to take into account when dealing with
          the m!=0 case.  First, while (L,l) are independent, not all m-values are allowed.
          In fact, abs(m) should be smaller than MIN(L,l).  Secondly, there is no need to
          compute J_Llm(x) for negative m's. The projection functions for negative m's
          can be obtained by flipping the sign of the second line of one of the 3j's, thus
          yielding a (-1)^(l+l1+L) sign. The factor l+l1+L is even for the direct projection
          functions (TT, EE and BB), while it is odd for the mixing ones (EB and BE). This 
          property comes from the the definition of H_X_X' in Beneke & Fidler 2010 (below eq. B.9).
          As a result:

          J_TT_Ll-m = J_TT_Llm
          J_EE_Ll-m = J_EE_Llm
          J_BB_Ll-m = J_BB_Llm
          J_EB_Ll-m = - J_EB_Llm
          J_BE_Ll-m = - J_BE_Llm.

          It follows from the above relations that for m=0, both J_EB and J_BE vanish.  */
          int index_m_max = MIN (ppr2->index_m_max[pbs2->L[index_L]], ppr2->index_m_max[pbs->l[index_l]]);
        
          /* Loop on the external m-index */
          for (int index_m = 0; index_m <= index_m_max; ++index_m) {
    
            int L = pbs2->L[index_L];
            int l = pbs->l[index_l];
            int m = pbs2->m[index_m];
  
            class_test_parallel (abs(m) > MIN(L,l),
              pbs2->error_message,
              "abs(m) should always be smaller that MIN(L,l), check ppr2->index_m_max");
                       
            class_call_parallel (bessel2_J_for_Llm(
                      ppr,
                      ppr2,
                      pbs,
                      pbs2,
                      index_J,
                      index_L,
                      index_l,
                      index_m),
                 pbs2->error_message,
                 pbs2->error_message);
  
            #pragma omp flush(abort)    
       
          } // end of for(index_m)
        } // end of for(index_l)
      } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
    } // end of for(index_L)
  } // end of loop on type of projection functions
  
  /* Determine the maximum size of the x-level in pbs2->J_Llm_x */
  pbs2->x_size_max_J = 0;
  for (int index_J = 0; index_J < pbs2->J_size; ++index_J)
    for (int index_L = 0; index_L < pbs2->L_size; ++index_L)
      for (int index_l = 0; index_l < pbs->l_size; ++index_l)
        for (int index_m = 0; index_m < MIN(ppr2->index_m_max[pbs2->L[index_L]],ppr2->index_m_max[pbs->l[index_l]])+1; ++index_m)   
          if (pbs2->x_size_J[index_J][index_L][index_l][index_m] > pbs2->x_size_max_J)
            pbs2->x_size_max_J = pbs2->x_size_J[index_J][index_L][index_l][index_m];
  

  /* Print information on the used memory */
  if (pbs2->bessels2_verbose > 1)
    printf (" -> memory in use to store the 2nd-order projection functions: ~ %3g MB\n",
    8*pbs2->count_allocated_Js/1e6);






  // ===================================================================================
  // =                               Spline interpolation                              =
  // ===================================================================================
  
  if (ppr->bessels_interpolation == cubic_interpolation) {
    
    /* Compute second derivatives of the J_Llm's in view of the spline interpolation needed
    in the second-order transfer module */
    
    for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {

      for (int index_L = 0; index_L < pbs2->L_size; ++index_L) {

        /* Beginning of parallel region */
        int abort = _FALSE_;    

        #pragma omp parallel shared(pbs2,abort)
        {

          #pragma omp for schedule (dynamic)
          for (int index_l = 0; index_l < pbs->l_size; ++index_l) {

            int index_m_max = MIN (ppr2->index_m_max[pbs2->L[index_L]], ppr2->index_m_max[pbs->l[index_l]]);
        
            for (int index_m = 0; index_m <= index_m_max; ++index_m) {

              /* TODO: Here I am doing a bit of free style, find a criterion */
              int SPLINE_METHOD = (pbs->l[index_l]==2 ? _SPLINE_EST_DERIV_:_SPLINE_NATURAL_);

              class_call_parallel (array_spline_table_one_column(
                                     pbs2->xx + pbs2->index_xmin_J[index_J][index_L][index_l][index_m],
                                     pbs2->x_size_J[index_J][index_L][index_l][index_m],
                                     pbs2->J_Llm_x[index_J][index_L][index_l][index_m],
                                     1,         /* Not used */    
                                     0,         /* We need to spline just one function */   
                                     pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m],
                                     SPLINE_METHOD,
                                     pbs2->error_message
                                     ),
                pbs2->error_message,
                pbs2->error_message);

            } // end of for(index_m)

            #pragma omp flush(abort)
          
          } // end of for(index_l)
        } if (abort == _TRUE_) return _FAILURE_;
      } // end of for(index_L)
    } // end of loop of projection function type


    /* Compute second derivatives of the j1_l's in view of the spline interpolation needed
    in the intrinsic bispectrum module */
    int abort = _FALSE_;    
    #pragma omp parallel for shared(pbs,abort)
    for (int index_l1 = 0; index_l1 < pbs2->l1_size; ++index_l1) {

        int SPLINE_METHOD = _SPLINE_NATURAL_;

        class_call_parallel (array_spline_table_one_column(
                               pbs2->xx + pbs2->index_xmin_l1[index_l1],
                               pbs2->x_size_l1[index_l1],
                               pbs2->j_l1[index_l1],
                               1,                                                 /* Not used */    
                               0,                                                 /* We need to spline just one function */   
                               pbs2->ddj_l1[index_l1],
                               SPLINE_METHOD,
                               pbs2->error_message
                               ),
          pbs2->error_message,
          pbs2->error_message);

      #pragma omp flush(abort)
    } // end of for(index_l1)
    if (abort == _TRUE_) return _FAILURE_;



  } // end of spline calculation


  return _SUCCESS_;

}












/**
 * Free all memory space allocated by bessel_init().
 *
 * To be called at the end of each run.
 *
 * @param pbs Input : Initialized bessel structure 
 * @return the error status
 */

int bessel2_free( 
    struct precision * ppr,
    struct precision2 * ppr2,
    struct bessels * pbs,
    struct bessels2 * pbs2
    )
{

  int index_l, index_J, index_L, index_m, index_l1;

  for (index_J = 0; index_J < pbs2->J_size; ++index_J) {
  
    for (index_L = 0; index_L < pbs2->L_size; index_L++) {
  
      for (index_l = 0; index_l < pbs->l_size; index_l++) {
  
        int index_m_max = MIN (ppr2->index_m_max[pbs2->L[index_L]], ppr2->index_m_max[pbs->l[index_l]]);
  
        for (index_m = 0; index_m <= index_m_max; ++index_m) {
        
          free(pbs2->J_Llm_x[index_J][index_L][index_l][index_m]);
          if (ppr->bessels_interpolation == cubic_interpolation)
            free(pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m]);
        
        }  // end of for(index_m)
      
        free(pbs2->index_xmin_J[index_J][index_L][index_l]);
        free(pbs2->x_size_J[index_J][index_L][index_l]);
        free(pbs2->x_min_J[index_J][index_L][index_l]);
        free(pbs2->J_Llm_x[index_J][index_L][index_l]);
        if (ppr->bessels_interpolation == cubic_interpolation)
          free(pbs2->ddJ_Llm_x[index_J][index_L][index_l]);
      
      }  // end of for(index_l)
    
      free(pbs2->index_xmin_J[index_J][index_L]);
      free(pbs2->x_size_J[index_J][index_L]);
      free(pbs2->x_min_J[index_J][index_L]);
      free(pbs2->J_Llm_x[index_J][index_L]);
      if (ppr->bessels_interpolation == cubic_interpolation)
        free(pbs2->ddJ_Llm_x[index_J][index_L]);
    
    }  // end of for(index_L)
  
    free(pbs2->index_xmin_J[index_J]);
    free(pbs2->x_size_J[index_J]);
    free(pbs2->x_min_J[index_J]);
    free(pbs2->J_Llm_x[index_J]);
    if (ppr->bessels_interpolation == cubic_interpolation)
      free(pbs2->ddJ_Llm_x[index_J]);
  
  } // end of for(index_J)
  
  free(pbs2->index_xmin_J);
  free(pbs2->x_size_J);
  free(pbs2->x_min_J);
  free(pbs2->J_Llm_x);
  if (ppr->bessels_interpolation == cubic_interpolation)
    free(pbs2->ddJ_Llm_x);
  
  
  free(pbs2->l1);
  
  free(pbs2->xx);
  
  free(pbs2->x_size_l1);
  free(pbs2->x_min_l1);
  for(index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
    free(pbs2->j_l1[index_l1]);
  free(pbs2->j_l1);
  if (ppr->bessels_interpolation == cubic_interpolation) {
    for(index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
      free(pbs2->ddj_l1[index_l1]);
    free(pbs2->ddj_l1);
  }
  
  
  free(pbs2->L);
  free(pbs2->m);

  return _SUCCESS_; 
}


/**
 * Define number and values of mutipoles l. This is crucial since not
 * only the Bessel functions, but also the transfer functions and
 * anisotropy spectra C_l will automatically be sampled at the same
 * values (there would be no logic in having various l lists differing
 * from each other).
 *
 *
 * @param ppr Input : pointer to precision structure
 * @param pbs Input/Output : pointer to besseld structure (result stored here) 
 * @return the error status
 */

int bessel2_get_l1_list(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct bessels * pbs,
          struct bessels2 * pbs2
          )
{

  /* The line of sight integration at second-order requires to compute a linear combination of
  spherical Bessel functions with 3j-symbols.  In particular, the Bessel function of order l1
  is required several times.  However, l1 goes from |l-L| to l+L, where the upper limit of L
  is given by pbs2->L_max.  This means that we need to compute the spherical Bessels in an
  extended range that includes the l's that are within a range of +/- pbs2->L_max from 
  each member of pbs2->l.
  
  Similarly, to compute the contribution to the bispectrum integral from a certain azimuthal number
  m, we will need to compute Bessel functions with order between |l-|m|| and l+|m|. Hence, we need
  to compute the spherical Bessels in an extended range that includes the l's that are within a range
  of +/- MAX(pbs2->L_max, ppr2->m_max_second_order) from each member of pbs2->l.*/
    
  if (pbs2->bessels2_verbose > 2) {
    printf (" -> will compute %d l's in the range l=(%d,%d)\n", pbs->l_size, pbs->l[0], pbs->l[pbs->l_size-1]);    
    printf ("    * l-list:\n       ");

    for (int index_l=0; index_l < pbs->l_size; ++index_l)
      printf ("%d /\\ ", pbs->l[index_l]);
      
    printf ("\n");
  }

  
  // *** Determine the number of elements in the extended list pbs2->l1

  /* The maximum number of elements pbs2->l1 can have is given by the maximum l in
  pbs2->l plus L_max=MAX(pbs2->L_max, ppr2->m_max_2nd_order). See the long comment
  above for details. We use a logical array (i.e. an array of 1s and 0s) to keep track
  of which l1's are to be kept. */
  int L_max = pbs2->L_max;
  int l1_max_size = pbs->l[pbs->l_size - 1] + L_max;
  short * l1_logical = (short *)calloc(l1_max_size+1, sizeof(int));

  /* Each of the l1's is to be kept only if it is in the range  abs(l-L_max) <= l1 <= l+L_max  */
  pbs2->l1_size = 0;

  for(int l1=0; l1<=l1_max_size; ++l1) {

    l1_logical[l1] = _FALSE_;                         // Initialize flag

    for(int index_l=0; index_l<pbs->l_size; ++index_l) {

      int l = pbs->l[index_l];

      if ( (l1>=(l-L_max)) && (l1<=l+L_max) )         // It is important to omit abs()
        l1_logical[l1] = _TRUE_;                      // Flag the l1 to be kept
                  
    }                                                 // end of for(index_l)
    
    if (l1_logical[l1] == _TRUE_)                     // Increment the counter
      pbs2->l1_size++;

  }                                                   // end of for(index_l1)

  // *** Some debug
  // for(index_l1=0; index_l1<l1_max_size; ++index_l1)
  //   if (l1_logical[index_l1] == _FALSE_)
  //     printf("We shall not compute the l1=%d bessel.\n", index_l1);
  // printf("pbs->l_size = %d, pbs2->l1_size = %d\n", pbs->l_size, pbs2->l1_size);
  
  
  // *** Write pbs2->l1 list using l1_logical
  class_alloc (pbs2->l1, pbs2->l1_size*sizeof(int), pbs2->error_message);
  int l1 = 0;
  for(int index_l1=0; index_l1<pbs2->l1_size; ++index_l1) {

    while (l1_logical[l1] == _FALSE_) l1++;             // Look for the first l1 to keep

    pbs2->l1[index_l1] = l1++;                          // And store it in pbs2->l1
    
    // Some debug
    // printf("pbs2->l1[%d] = %d\n", index_l1, pbs2->l1[index_l1]);
  }

  free(l1_logical);

  /* Find out the index in pbs2->l1 corresponding to a given l. */
  class_alloc (pbs2->index_l1, (pbs2->l1[pbs2->l1_size-1]+1)*sizeof(int), pbs2->error_message);
  for(l1=0; l1<=pbs2->l1[pbs2->l1_size-1]; ++l1) {
  
    pbs2->index_l1[l1] = -1;
  
    for (int index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
      if (l1==pbs2->l1[index_l1]) pbs2->index_l1[l1] = index_l1;
    
    /* Some debug */
    // printf("pbs2->index_l1[%d] = %d\n", l1, pbs2->index_l1[l1]);
  }

  return _SUCCESS_;

}





/** 
 * Fill the array pbs2->xx.  This is the grid in x where J_Llm(x) will be sampled.
 *
 * For each (L,l,m), the sampling will start at a different index_x because 
 * J_Llm(x) is negligible if x is close enough to zero (unless L=l).
 * 
 * 
 * @param ppr       Input: pointer to precision structure
 * @param pbs       Input/Output: pointer to bessels structure
 * @return the error status
 */
int bessel2_get_xx_list(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      )
{
  
  /* We shall use linear sampling starting at zero */
  pbs2->xx_size = pbs2->xx_max/pbs2->xx_step + 1;

  class_alloc (pbs2->xx, pbs2->xx_size*sizeof(double), pbs2->error_message);
  lin_space (pbs2->xx, 0, pbs2->xx_max, pbs2->xx_size);
  
  return _SUCCESS_;

}





/** 
 * Bessel factors J_Llm for arbitrary argument x. 
 *
 * Evaluates the J_Llm factor needed for the 2nd-order line of sight integrations at a given
 * value of x by interpolating in the pre-computed table pbs2->J_Llm_x.  This function can be
 * called from whatever module at whatever time, provided that bessel_init() has been called
 * before, and bessel_free() has not been called yet.
 *
 * @param pbs       Input: pointer to bessels structure
 * @param x         Input: argument x
 * @param index_L   Input: index defining L = pbs2->L[index_L]
 * @param index_l   Input: index defining l = pbs->l[index_l]
 * @param index_m   Input: index defining m = pbs2->m[index_m]
 * @param J_Llm_x   Ouput: J_Llm(x)
 * @return the error status
 */
int bessel2_J_Llm_at_x (
    struct bessels * pbs,
    struct bessels2 * pbs2,
    double x,
    int index_L,
    int index_l,
    int index_m,
    double * J_Llm_x
    )
{


  return _SUCCESS_;

}





/**
 * Compute the J_Llm(x) coefficients for all the allowed values of L.  The result
 * will be stored in pbs2->J_Llm_x[index_J][index_L][index_l][index_m][index_x] where:
 *
 * - the number of 'l' indices is determined in bessel_get_l_list
 * - the number of 'm' indices depends on the modes (scalar, vector, tensor) asked by the user
 * - the number of 'L' indices is determined by the option pbs2->L_max
 * - the 'x' values are x=x_min_J+xx_step*index_x, where x_min_J is determined
 *   in this function through bisection, xx_step comes from the parameter file (default 0.3),
 *   and the number of steps is simply xx_size=xx_max-xx_min/(xx_step+1), with xx_max computed in
 *   input.c as k_max*MAX(tau-tau0) ~ k_max*tau0.
 *
 *
 * @param ppr Input : pointer to precision structure
 * @param pbs Input/Output : pointer to bessel structure (store result in pbs2->J_Llm_x) 
 * @param l   Input: l value
 * @param m   Input: m value
 * @return the error status
 *
 * TODO:
 *    - Create xx_max with second-order k_max*tau_0
 *
 */

int bessel2_J_for_Llm (
       struct precision * ppr,
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       int index_J,
       int index_L,
       int index_l,
       int index_m
       )
{

  // ==============================================================
  // =                        Preparations                        =
  // ==============================================================

  /* Values of L,l,m */
  int L = pbs2->L[index_L];
  int l = pbs->l[index_l];
  int m = pbs2->m[index_m];
  int S;

  /* Determine the type of projection function to compute. Each of them has a spin associated,
  according to appendix B of http://arxiv.org/abs/1102.1524 (S = -F). */
  enum projection_function_types projection_function;
  
  if ((pbs2->has_J_TT==_TRUE_) && (index_J==pbs2->index_J_TT)) {
    projection_function = J_TT;
    S = 0;
  }
  else if ((pbs2->has_J_EE==_TRUE_) && (index_J==pbs2->index_J_EE)) {
    projection_function = J_EE;
    S = 2;
  }  
  else if ((pbs2->has_J_EB==_TRUE_) && (index_J==pbs2->index_J_EB)) {
    projection_function = J_EB;
    S = 2;
  }


  /* The spin S appears in the definition of the projection functions in the second line of the
  3j-symbol, below both l and L. Here we skip those values of l and L that are smaller than S. */
  short constraint_spin = ((abs(S)>l) || (abs(S)>L));

  /* If m=0, the projection function for the E-B mixing vanishes (for details, refer
  to the long comment in init function). */
  short constraint_EB = ((projection_function==J_EB) && (m==0));

  /* If you want to skip the forbidden combinations of (S,L,l,m), you should put the relevant
  condition in the following if block. It is not sufficient to just write "if(condition) continue";
  we also need to assign values to the x_min, x_size, J_Llm_x arrays in order to avoid segmentation
  faults. */  
  if (constraint_spin || constraint_EB) {

    pbs2->x_size_J[index_J][index_L][index_l][index_m] = 1;
    pbs2->index_xmin_J[index_J][index_L][index_l][index_m] = pbs2->xx_size - 1;
    pbs2->x_min_J[index_J][index_L][index_l][index_m] = pbs2->xx_max;
    class_calloc (pbs2->J_Llm_x[index_J][index_L][index_l][index_m],1,sizeof(double),pbs2->error_message);
    pbs2->J_Llm_x[index_J][index_L][index_l][index_m][0] = 0.;
    if (ppr->bessels_interpolation == cubic_interpolation) {
      class_calloc (pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m],1,sizeof(double),pbs2->error_message);
      pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m][0] = 0.;
    }

    if (pbs2->bessels2_verbose > 2)
      printf("     \\ skipping J(%d,%d,%d,%d) due to the geometrical constraints\n",
        index_J, L, l, m);

    return _SUCCESS_;
  }


  /* Some debug */
  if (pbs2->bessels2_verbose > 2)
    printf("    * bessel_J_for_Llm: processing (index_J,index_L,index_l,index_m) = (%d,%d,%d,%d), (L,l,m) = (%d,%d,%d)\n",
    index_J, index_L, index_l, index_m, L, l, m);

  /* Index for x and value x=x_min_J(L,l,m)+xx_step*index_x */
  int index_x = 0;
  double x;

  /* Value of J_Llm(x) */
  double J = 0.;

  /* For computing x_min */
  double x_min_up, x_min_down, x_min_J;
  int index_xmin_up, index_xmin_down, index_xmin_J;
  
  /* Flags used to determine whether J_Llm(x) is negligible for the x-domain we are
  interested in, i.e. for all x's in pbs2->xx */
  short might_be_negligible = _FALSE_;
  short is_negligible = _FALSE_;
  

  // ==============================================================
  // =                        Compute 3j's                        =
  // ==============================================================


  /* We shall memorize the results of the 3j computation inside the simple
  ad-hoc structure bessel_3j_data-> */
  struct J_Llm_data * bessel_3j_data;
  class_alloc (bessel_3j_data, sizeof(struct J_Llm_data), pbs2->error_message);

  /* Allocate temporary arrays to store spherical Bessel functions & 3j's */
  class_alloc (bessel_3j_data->bessels, (2*pbs2->L_max+1)*sizeof(double), pbs2->error_message);
  class_alloc (bessel_3j_data->first_3j, (2*pbs2->L_max+1)*sizeof(double), pbs2->error_message);
  class_alloc (bessel_3j_data->second_3j, (2*pbs2->L_max+1)*sizeof(double), pbs2->error_message);
  double l1_min_D, l1_max_D;

  /* Compute the 3j symbol
    (    l     l1      L   )
    (   -S      0      S   )
  for all allowed values of 'l1'. Since l1 has to be first, we rearrange it as
    (   l1      L      l   )
    (    0      S     -S   )
  */  
  class_call (drc3jj(
                L, l, S, -S,
                &l1_min_D, &l1_max_D,
                bessel_3j_data->first_3j,
                (2*pbs2->L_max+1),
                pbs2->error_message       
                ),
    pbs2->error_message,
    pbs2->error_message);

  bessel_3j_data->l1_min = (int)(l1_min_D + _EPS_);
  bessel_3j_data->l1_max = (int)(l1_max_D + _EPS_);
  bessel_3j_data->l1_size = bessel_3j_data->l1_max - bessel_3j_data->l1_min + 1;

  /*   Compute the 3j symbol
     (    l     l1      L   )
     (    m      0     -m   )
  for all allowed values of 'l1'. Since l1 has to be first, we rearrange it as
    (   l1      L      l   )
    (    0     -m      m   )
  Careful here: we compute the second 3j only if it is different from the first one
  */
  if ((S!=0) || (m!=0)) {
    class_call (drc3jj(
                  L, l, -m, m,
                  &l1_min_D, &l1_max_D,
                  bessel_3j_data->second_3j,
                  (2*pbs2->L_max+1),
                  pbs2->error_message       
                  ),
      pbs2->error_message,
      pbs2->error_message);
  }
  else {
    for(int index_l1=0; index_l1<bessel_3j_data->l1_size; ++index_l1)
      bessel_3j_data->second_3j[index_l1] = bessel_3j_data->first_3j[index_l1];
  }

  /* Print information on the sum over l1 */
  if (pbs2->bessels2_verbose>2)
    printf("    * bessel_J_for_Llm: the sum over l1 will run from %d to %d.\n",
    bessel_3j_data->l1_min, bessel_3j_data->l1_max);      

  /* Compute the position of l1_min in the array pbs2->l1 */
  bessel_3j_data->index_l1_min = pbs2->index_l1[bessel_3j_data->l1_min];
  class_test (bessel_3j_data->index_l1_min < 0,
    pbs2->error_message,
    "error in ell_1 list");



  // ========================================================================
  // =                    Find where J isn't negligible                     =
  // ========================================================================
  
  /* When L==l, we have a contribution to J_Llm(x) coming from the j_0(x) = sin(x)/x.
  Unlike all the other spherical Bessels functions, j_0(x) does not vanish for x->0
  but it asymptotes to 1.  Hence, there is no regime close to x=0 where J_llm(x) is
  negligible, and when L=l we just take x_min=0 */
  if(L==l) {
    index_xmin_J = 0;  
    x_min_J = 0.;
    is_negligible = _FALSE_;
  }
  // *** If L!=l, do the bisection ***
  else {
  
    /* 'x_min_up' is the starting point for the bisection.  We shall look for the first x<=x_min_up
    that makes J_Llm(x) smaller than pbs2->J_Llm_cut. Since J_Llm(x) involves a sum over spherical Bessel
    functions with order  |l-L| <= l1 <= l+L, we shall set the initial guess in such a way that
    j_l1(x_min_up) is always non-negligible, and not into the oscillatory regime.  This is attained
    by choosing x_min_up = minimum allowed vale of l1 + 0.5 */
    x_min_up = bessel_3j_data->l1_min + 0.5;
    x_min_down = 0.;
    index_xmin_down = 0;

    // *** Determine index_xmin_up
    if (x_min_up >= pbs2->xx_max) {

      /* When our guess is larger than the maximum value of x in pbs2->xx, we scale it down to
      the maximum value contained in pbs2->xx.  Most likely, this means that the J_Llm(x) is
      negligible for all values of x in pbs2->xx, as J_Llm is just a collection of spherical
      Bessel functions. */
      index_xmin_up = pbs2->xx_size-1;
      x_min_up = pbs2->xx_max;
      might_be_negligible = _TRUE_;
    }
    else {

      /* Find the index in pbs2->xx following x_min_up */
      index_xmin_up = 0;
      while(pbs2->xx[index_xmin_up] < x_min_up) index_xmin_up++;
      
      /* Update the value of x_min_up to reflect an actual entry of the xx array */
      x_min_up = pbs2->xx[index_xmin_up];
    }
    
  
    /* Value of J_Llm in our upper limit */
    class_call (bessel2_J_Llm(
                  ppr2,
                  pbs,
                  pbs2,
                  projection_function,
                  L,
                  l,
                  m,
                  index_xmin_up,
                  bessel_3j_data,
                  &J),
      pbs2->error_message,
      pbs2->error_message);
  
    // *** Some debug
    // if ((L==4) && (l==10))
    //   printf("bisection: index_xmin_down=%d, index_xmin_up=%d\n", index_xmin_down, index_xmin_up);
    // if ((L==4) && (l==10))
    //   printf("bisection: x_min_down=%g, x_min_up=%g, J=%g, J_Llm_cut=%g\n", x_min_down, x_min_up, J, pbs2->J_Llm_cut);

    /* If the value of J_Llm in x_min_up is already below the target cut, then it means
    that we hit a point where the J_Llm is negligible with our first attempt.  If xx_max
    is smaller than l1_min+0.5, then this means that J_Llm(x) is negligible for all the
    needed values of x, and the bisection cycle is not needed.  On the other hand, if
    J_Llm(x_min_up) is negligible but xx_max > x_min_up = l1_min+0.5, it means that either our
    initial guess was way too high and the Bessels are all very damped (x_min_up >> l1_min+0.5),
    or that we entered into the oscillatory regime of J_Llm and hit an x-value corresponding
    to a zero crossing.  In both cases, we return an error.  */
    if (fabs(J) < pbs2->J_Llm_cut)                // if cut is already satisfied
      if (might_be_negligible)                    // if pbs2->xx_max < l1_min+0.5
        is_negligible == _TRUE_;
      else                                        // if pbs2->xx_max > l1_min+0.5
        class_test(fabs(J) < pbs2->J_Llm_cut,      
           pbs2->error_message,
           "in bisection, wrong initial guess for x_min_up.");

    // --------------------------------------------------------------
    // -                       Actual bisection                     -
    // --------------------------------------------------------------
  
    /* We use the same tolerance as for the normal Bessels */  
    if (!is_negligible) while (index_xmin_down!=(index_xmin_up-1)) {

      int index_xmid = (index_xmin_up+index_xmin_down)/2;
      
      class_call (bessel2_J_Llm (
                    ppr2,
                    pbs,
                    pbs2,
                    projection_function,
                    L,
                    l,
                    m,
                    index_xmid,
                    bessel_3j_data,
                    &J),
        pbs2->error_message,
        pbs2->error_message);  

      if (fabs(J) >= pbs2->J_Llm_cut) {
        index_xmin_up = index_xmid;
        x_min_up = pbs2->xx[index_xmid];
      }
      else {
        index_xmin_down = index_xmid;
        x_min_down = pbs2->xx[index_xmid];
      }
  
      // *** Some debug
      // if ((L==4) && (l==10))
      //   printf("bisection: index_xmin_down=%d, index_xmin_up=%d\n", index_xmin_down, index_xmin_up);
      // if ((L==4) && (l==10))
      //   printf("bisection: x_min_down=%g, x_min_up=%g, J=%g, J_Llm_cut=%g\n", x_min_down, x_min_up, J, pbs2->J_Llm_cut);

    } // end of if(!is_negligible)
  
    /* End of bisection */
    index_xmin_J = index_xmin_down;
    x_min_J = x_min_down;

  } // end of if (l!=L)


  /* Define number of x values to be stored (one if all values of j_l(x) were negligible for this l) */  
  if (is_negligible)
    pbs2->x_size_J[index_J][index_L][index_l][index_m] = 1;
  else
    pbs2->x_size_J[index_J][index_L][index_l][index_m] = pbs2->xx_size - index_xmin_J;


  // *** Some debug
  // printf("pbs2->x_size_J[%d][%d][%d][%d] = %d, pbs->x_max = %g\n",
  //   index_J, index_L, index_l, index_m, pbs2->x_size_J[index_J][index_L][index_l][index_m], pbs->x_max);


  // =========================================================================
  // =                    Compute the projection function                    =
  // =========================================================================

  /* Allocate the x-level of J_Llm_x, taking care that it is initialised to zero */
  int x_size_J = pbs2->x_size_J[index_J][index_L][index_l][index_m];

  class_calloc (pbs2->J_Llm_x[index_J][index_L][index_l][index_m], x_size_J, sizeof(double), pbs2->error_message);
  #pragma omp atomic
  pbs2->count_allocated_Js += x_size_J;
  
  if (ppr->bessels_interpolation == cubic_interpolation) {
    class_calloc (pbs2->ddJ_Llm_x[index_J][index_L][index_l][index_m], x_size_J, sizeof(double), pbs2->error_message);
    #pragma omp atomic
    pbs2->count_allocated_Js += x_size_J;
  }
  
  /* Shortcut for J_Llm_x */
  double * J_Llm_x = pbs2->J_Llm_x[index_J][index_L][index_l][index_m];
  
  
  // *** Case when all values of J_Llm(x) were negligible for this (L,l,m)
  if (is_negligible) {
  
    if(pbs2->bessels2_verbose > 0)
      printf("Setting the J(L=%d,l=%d,m=%d)=0 for xx_max=%g\n", L,l,m,pbs2->xx_max);
      
    pbs2->index_xmin_J[index_J][index_L][index_l][index_m] = pbs2->xx_size - 1;
    pbs2->x_min_J[index_J][index_L][index_l][index_m] = pbs2->xx_max;
    J_Llm_x[0] = 0.;
    
  }
  // *** Otherwise, write first non-negligible value and then loop over x
  else {
  
    pbs2->index_xmin_J[index_J][index_L][index_l][index_m] = index_xmin_J;
    pbs2->x_min_J[index_J][index_L][index_l][index_m] = x_min_J;
    
    /* Loop over all other non-negligible values */
    for (index_x=index_xmin_J; index_x < pbs2->xx_size; index_x++) {
    
      class_call (bessel2_J_Llm(
                    ppr2,
                    pbs,
                    pbs2,
                    projection_function,
                    L,
                    l,
                    m,
                    index_x,
                    bessel_3j_data,
                    &(J_Llm_x[index_x - index_xmin_J])
                    ),
        pbs2->error_message,
        pbs2->error_message);
    
      /* Some debug: print the J_Llm(x) to screen */
      // if ((L==0) && (l==100) && (m==0))
      //   fprintf(stderr, "%g %g\n", pbs2->xx[index_x], J_Llm_x[index_x - index_xmin_J]);
    
    } // end of for(index_x)
  } // end of if(x_min_J >= xx_max)


  // *** Free 3j-symbol arrays
  free(bessel_3j_data->bessels);
  free(bessel_3j_data->first_3j);
  free(bessel_3j_data->second_3j);
  free(bessel_3j_data);
    
  return _SUCCESS_;
  
}






/**
 * Compute the J coefficients using spherical Bessel function and the 3j symbols.
 * The result will depend on the given values of L, l, m and x, where:
 * 
 * - 'L' is the index that will be summed with the sources in the line-of sight
 *   integration.
 *
 * - 'l' and 'm' are the indices associated to the (l,m)-th multipole for which we are
 *   computing the line of sight integral.
 *
 * - 'x' is the argument of the Bessel function, that will be equal to k(tau-tau0)
 *   in the line of sight integration.
 * 
 * @param pbs Input : pointer to bessel structure (used only for error message)
 * @param L   Input: L value
 * @param l   Input: l value
 * @param m   Input: m value
 * @param x   Input: x value
 * @param J_Llm_x  Output: J_Llm(x) value
 * @return the error status
 */

int bessel2_J_Llm (
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       enum projection_function_types projection_function,
       int L,
       int l,
       int m,
       int index_x,
       struct J_Llm_data * bessel_3j_data,
       double * J_Llm_x
       )
{

  /* Considered x */
  double x = pbs2->xx[index_x];
  
  /* We are going to sum over l1 */
  int l1, index_l1, index_x_in_jl1, index_l1_in_jl1;
  double j_l1_x = 0;         // Bessel function j_l1(x)
  double i_prefactor;        // Coefficient determining the sign of each l1 contribution
  (*J_Llm_x) = 0;            // We shall accumulate J over l1

  /* Summing up to 'l1_size' is equivalent to sum over all the allowed l1's.
  Note that this index_l1 IS NOT the one that should be used to access
  pbs2->l1.  To do so, you need to use bessel_3j_data->index_l1_min + index_l1. */
  for(index_l1=0; index_l1<bessel_3j_data->l1_size; ++index_l1) {


    // -------------------------------------------
    // -       Skip vanishing contributions      -
    // -------------------------------------------

    /* Skip the contribution when the three-j vanishes */
    if (bessel_3j_data->first_3j[index_l1]==0.) continue;

    /* Obtain the value of l1 */
    index_l1_in_jl1 = bessel_3j_data->index_l1_min + index_l1;
    index_x_in_jl1  = index_x - pbs2->index_xmin_l1[index_l1_in_jl1];
    l1 = pbs2->l1[index_l1_in_jl1];

    /* If the Bessel function that is being summed is smaller than the treshold set in 
    ppr2->bessel_j_cut_2nd_order, then skip it */
    if (index_x_in_jl1 < 0) {
      if (pbs2->bessels2_verbose > 2)
        printf("     \\ Assumed that the l1=%d contribution to J_%d_%d_%d(%g) vanishes (i.e. j_%d(%g) < %g)\n",
          l1, L, l, m, x, l1, x, pbs2->j_l1_cut);      
      continue;
    }

    /* Use the parity properties of the projection function to skip vanishing contributions.
    J_TT, J_EE and J_BB are different from zero only for even l+l1+L. For J_TT, this is
    enforced by the three-j symbol. For J_EE and J_BB, this is enforced by the parity of
    the mixing matrix H (see eq. 76 of Beneke & Fidler 2010). The mixing matrix also enforces
    that J_EB and J_BE are different from zero only for odd l+l1+L.
    
    The prefactor is i^(l-l1-L) for the even parity projection functions (J_TT, J_EE, J_BB), and
    either i^(l-l1-L-1) or i^(l-l1-L+1) for the odd parity ones (J_EB and J_BE respectively). In both
    cases, the parity forces the prefactors to be just (real) alternating signs. */
    short is_even = ((l-l1-L)%2==0);

    if ((projection_function==J_TT)||(projection_function==J_EE))
      if (is_even)
        i_prefactor = alternating_sign((l-l1-L)/2);
      else
        continue;

    else if (projection_function==J_EB)
      if (!is_even)
        i_prefactor = alternating_sign((l-l1-L-1)/2);
      else
        continue;



    // -------------------------------------------
    // -          Build l1 contribution          -
    // -------------------------------------------
    
    /*  Value of the spherical Bessel function j_l(x), precomputed in pbs2->j_l1[index_l1][index_x].
    Note that the l1 level of pbs2->j_l1 should be addressed the same way as pbs2->l1, while
    the x level should be addressed as index_x - pbs2->index_xmin_l1[index_l1], where index_x
    is the index of pbs2->xx. */
    j_l1_x = pbs2->j_l1[index_l1_in_jl1][index_x_in_jl1];

    // *** Some debug
    // if ((L==4) && (l==1172)) printf("j(%d, %10.6g) = %10.6g\n", l1, x, j_l1_x);
    // if ((L==4) && (l==1172)) printf("3j(%d,%d,%d)(%d,0,%d) = %10.6g\n", l1, L, l, -m, m, bessel_3j_data->second_3j[index_l1]);    
    
    /* Increment the result */
    (*J_Llm_x) +=   i_prefactor * (2*l1+1)
                  * bessel_3j_data->first_3j[index_l1]
                  * bessel_3j_data->second_3j[index_l1]
                  * j_l1_x;
    
  } // end of for(index_l1)


  /* Apply constant factors */
  (*J_Llm_x) *= alternating_sign(m) * (2*l+1);

  /* Some debug */
  // if ((projection_function==J_EB) && (L==2) && (l==102) && (m==1))
  //   fprintf(stderr, "%10.6g %10.6g\n", x, *J_Llm_x);
  //   // fprintf(stderr, "J(%d,%d,%d,%10.6g[%d]) = %10.6g\n", L, l, m, x, index_x, *J_Llm_x);
  
  /* Mathematica batch debug */
  // double val = *J_Llm_x;
  // if (projection_function == J_TT) {
  //   if (val==0.) {
  //     fprintf (stderr, "Abs[J[%d,%d,%d][%g]]<10^%g,\n",
  //     L,l,m,x,log10(ppr2->bessel_J_cut_2nd_order));
  //   }
  // }
  // double val = *J_Llm_x;
  // if (projection_function == J_TT) {
  //   if (val!=0) {
  //     if ((l<250) && (l<250)) {
  //       if (m==1) {
  //         if ((index_x<3) || (index_x%128==0)) {
  //           fprintf (stderr, "Abs[1 - J[%d,%d,%d][%g]/((%g)*10^(%f))] < 10^-2,\n",
  //             L,l,m,x,sign(val),log10(fabs(val)));
  //         }
  //       }
  //     }
  //   }
  // }
  
  return _SUCCESS_;

}



/**
 * Compute J_Llm(x) without relying on the precomputed spherical Bessel
 * in pbs2->j_l1.  To do so, use the SLATEC subroutine BESJ.  This function
 * might be useful to debug bessel_J_Llm, but it is never used in the
 * code.
 */
int bessel2_J_Llm_at_x_exact(
       struct bessels * pbs,
       struct bessels2 * pbs2,
       int L,
       int l,
       int m,
       double x,
       struct J_Llm_data * bessel_3j_data,
       double * J_Llm_x
       )
{
  
  
  /* Spherical Bessel functions in x for all the needed values of l1. We
  do not include yet the sqrt(pi/(2x)) factor that turns the Bessel functions J_l(x) into
  the spherical ones j_l(x).  We shall do it at a later time in order to save time on the
  computation of sqrt(x).  Also note that the function below is a single precision one.
  This gives ~7 digits of precision, which is enough for any cosmological application.
  Comment the following line if you want to use CLASS version of the spherical Bessels. */
  besselJ_l1(
     (float)bessel_3j_data->l1_min + 0.5,    // j_l ~ J_(l+1/2)
     (float)x,                          
     bessel_3j_data->l1_size,
     bessel_3j_data->bessels,
     pbs2->error_message          
     );

  // *********      Sum over l1     ***********
  int l1, index_l1;
  double j_l1_x = 0;         // Bessel function j_l(x)
  double i_coefficient;      // Coefficient i^(l-l1-L)
  (*J_Llm_x) = 0;            // We shall accumulate J over l1

  /* Summing up to 'l1_size' is equivalent to sum over all the allowed l1's */
  for(index_l1=0; index_l1<bessel_3j_data->l1_size; ++index_l1) {
  
    /* The first 3j-symbol vanishes for odd l-l1-L */
    if (bessel_3j_data->first_3j[index_l1]==0.) continue;
    
    /* Value of l1 */
    l1 = bessel_3j_data->l1_min + index_l1;
    
    /*  Due to the symmetries of the first 3j-symbol, we need to consider only
      even l-l1-L. Hence, the coefficient i^(l-l1-L) depends on the parity of
      of (l-l1-L)/2.  */
    i_coefficient = alternating_sign((l-l1-L)/2);
    
    /* Value of the spherical Bessel function j_l(x), comment if you want
    to use SLATEC version. */
    // class_call (bessel_j(pbs,
    //            l1,
    //            x, 
    //            &j_l1_x),
    //         pbs2->error_message,
    //         pbs2->error_message);

    /* Value of the spherical Bessel function j_l(x), comment if you want
    to use CLASS version. */    
    j_l1_x = bessel_3j_data->bessels[index_l1];
    
    /* Increment the result */
    (*J_Llm_x) += i_coefficient * (2*l1+1) * bessel_3j_data->first_3j[index_l1] * bessel_3j_data->second_3j[index_l1] * j_l1_x;
    
  } // end of for(index_l1)

  /*  Include the sqrt(pi/(2x)) factor that converts the J_(l+1/2) into j_l.
    Comment if you are using CLASS spherical Bessel functions. */
  (*J_Llm_x) *= sqrt_pi_over_2/sqrt(x);

  /* Apply constant factors */
  (*J_Llm_x) *= alternating_sign(m) * (2*l+1);

  // *** Some debug
  // printf("J(%d,%d,%d,%20.15g) = %.15g\n", L, l, m, x, *J_Llm_x);
  
  return _SUCCESS_;
  
}







/** 
 * Spherical Bessel function for arbitrary argument x, with order l1. 
 *
 * Evaluates the spherical Bessel function x at a given value of x by
 * interpolating in the pre-computed table.  This function can be
 * called from whatever module at whatever time, provided that
 * bessel_init() has been called before, and bessel_free() has not
 * been called yet.
 *
 * The index for l1 refers to the array pbs2->l1, which is an extension
 * of pbs2->l.
 *
 * @param pbs      Input: pointer to bessels structure
 * @param x        Input: argument x
 * @param index_l1 Input: index defining l1 = pbs2->l1[index_l1]
 * @param j        Ouput: j_l1(x)
 * @return the error status
 */
int bessel2_l1_at_x(
    struct bessels2 * pbs2,
    double x,
    int index_l1,
    double * j_l1
    )
{
  
  /* If index_l1 is too large to be in the interpolation table, return  an error */
  class_test(index_l1 > pbs2->l1_size,
             pbs2->error_message,
             "index_l1=%d>l1_size=%d; increase l_max.", index_l1, pbs2->l1_size);


  /* If x is too small to be in the interpolation table, return 0 */
  if (x < pbs2->x_min_l1[index_l1]) {
    
    *j_l1 = 0;
    return _SUCCESS_;
    
  }
  else {

    /* If x is too large to be in the interpolation table, return an error  */
    class_test(x > pbs2->xx_max,
               pbs2->error_message,
               "x=%e>xx_max=%e in bessel structure",x,pbs2->xx_max);


    /* Find index_x, i.e. the position of x in the table; no complicated algorithm needed,
      since values are linearly spaced with a known step and known first value. */

    int index_x = (int)((x-pbs2->x_min_l1[index_l1])/pbs2->xx_step);
    double a = (pbs2->x_min_l1[index_l1]+pbs2->xx_step*(index_x+1) - x)/pbs2->xx_step;

    /* Find result with spline interpolation */

    *j_l1 = a * pbs2->j_l1[index_l1][index_x] 
      + (1.-a) * ( pbs2->j_l1[index_l1][index_x+1]
          - a * ((a+1.) * pbs2->ddj_l1[index_l1][index_x]
           +(2.-a) * pbs2->ddj_l1[index_l1][index_x+1]) 
          * pbs2->xx_step * pbs2->xx_step / 6.0);

  }
  

  return _SUCCESS_;

}






int bessel2_l1_at_x_linear(
    struct bessels2 * pbs2,
    double x,
    int index_l1,
    double * j_l1
    )
{
  

  /* Define local variables */
  int index_x;
  int l1 = pbs2->l1[index_l1];


  /* If index_l1 is too large to be in the interpolation table, return  an error */
  class_test(index_l1 > pbs2->l1_size,
             pbs2->error_message,
             "index_l1=%d>l1_size=%d; increase l_max.", index_l1, pbs2->l1_size);


  /* If x is too small to be in the interpolation table, return 0 */
  if (x < pbs2->x_min_l1[index_l1]) {
    
    *j_l1 = 0;
    return _SUCCESS_;
    
  }
  else {

    /* If x is too large to be in the interpolation table, return an error  */
    class_test(x > pbs2->xx_max,
               pbs2->error_message,
               "x=%e>xx_max=%e in bessel structure",x,pbs2->xx_max);


    /* Find index_x, i.e. the position of x in the table; no complicated algorithm needed,
      since values are linearly spaced with a known step and known first value. */

    index_x = (int)((x-pbs2->x_min_l1[index_l1])/pbs2->xx_step);

    /* Find result with linear interpolation */

    // double delta_x, a, y_1, y_2, h, b, x_1, x_2;
    // delta_x = pbs2->x_min_l1[index_l1] + pbs2->xx_step * (index_x+1) - x;
    // a = delta_x / pbs2->xx_step;
    double a, y_1, y_2, h, b, x_1, x_2;
    y_1 = pbs2->j_l1[index_l1][index_x];
    y_2 = pbs2->j_l1[index_l1][index_x+1];

    x_1 = pbs2->x_min_l1[index_l1] + (index_x)*pbs2->xx_step;
    x_2 = pbs2->x_min_l1[index_l1] + (index_x+1)*pbs2->xx_step;

    h = x_2 - x_1;
    b = (x-x_1)/h;
    a = 1-b;

    *j_l1 = a * y_1 + b * y_2;

    // *** Some debug
    // if (0==0) {
    //   printf("from table: j(%d, %20.15g) = %g\n", l1, pbs2->x_min_l1[index_l1] + (index_x)*pbs2->xx_step, pbs2->j_l1[index_l1][index_x]);
    //   printf("from table: j(%d, %20.15g) = %g\n", l1, pbs2->x_min_l1[index_l1] + (index_x+1)*pbs2->xx_step, pbs2->j_l1[index_l1][index_x+1]);
    //   printf("interp    : j(%d, %20.15g) = %g\n", l1, x, *j_l1);
    // }
    
  }

  return _SUCCESS_;

}




/**
 * Get spherical Bessel functions for given value of l1.
 *
 * This function it is equivalent to bessel_j_for_l. The only difference
 * is that it computes the spherical Bessels for the extended list pbs2->l1
 * instead that pbs2->l.
 *
 * @param ppr Input : pointer to precision structure
 * @param pbs Input/Output : pointer to bessel structure (store result here) 
 * @param index_l1 Input : index of the multipole
 * @return the error status
 */

int bessel2_j_for_l1(
       struct precision * ppr,
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       int index_l1
       )
{
  

  /** - Define local variables */

  /* Index for x and value x=x_min[index_l1]+xx_step*index_x */
  int index_x;
  double x;

  /* Value of j_l1(x) */
  double j;

  /* For computing x_min by bisection */
  double x_min_up;
  double x_min_down;
  double x_min;

  index_x = 0;
  j = 0.;
  int l1 = pbs2->l1[index_l1];


  /*   j_0(x) = sin(x)/x does not vanish for x->0, but it asymptotes to 1.
     Hence, there is no regime close to x=0 where J_llm(x) is negligible,
     and we just take x_min=0 when l1=0.  */
  if(l1==0) {
    x_min = 0.;
    pbs2->index_xmin_l1[index_l1] = 0;
  }
  else {

    /** - Find x_min_l1[index_l1] by bisection */
    x_min_up=(double)pbs2->l1[index_l1]+0.5;
    x_min_down=0.;

    class_call (bessel_j (pbs,
            pbs2->l1[index_l1],   // l1
            x_min_up,            // x
            &j),                 // j_l1(x)
         pbs->error_message,
         pbs2->error_message);
  
    class_test(j < pbs2->j_l1_cut,
         pbs2->error_message,
         "in bisection, wrong initial guess for x_min_up.");
 
    while ((x_min_up-x_min_down)/x_min_down > ppr->bessel_tol_x_min) {

      /* When we get very close to the lower limit, just start sampling the Bessel from there */
      if ((x_min_up-x_min_down) < ppr->smallest_allowed_variation) {
        x_min = x_min_up = x_min_down;
        break;
      }
      
      class_call (bessel_j(pbs,
        pbs2->l1[index_l1],
        0.5 * (x_min_up+x_min_down),
        &j),
           pbs->error_message,
           pbs2->error_message);
    
      if (j >= pbs2->j_l1_cut) 
        x_min_up=0.5 * (x_min_up+x_min_down);
      else
        x_min_down=0.5 * (x_min_up+x_min_down);

    }
  
    x_min = x_min_down;


    /* Find the index in pbs2->xx following x_min */
    if (x_min < pbs2->xx_max) {
    
      pbs2->index_xmin_l1[index_l1] = 0;
      while (pbs2->xx[pbs2->index_xmin_l1[index_l1]] < x_min) pbs2->index_xmin_l1[index_l1]++;
    
      /* We shall store the spherical Bessels from the index preceding x_min */
      if (pbs2->index_xmin_l1[index_l1] != 0) --pbs2->index_xmin_l1[index_l1];
    
      /* Consistency check */
      class_test (((pbs2->index_xmin_l1[index_l1] < 0) || (pbs2->index_xmin_l1[index_l1] > pbs2->xx_size-1)),
                  pbs2->error_message,
                  "stopping to prevent segmentation fault");
    } 
      
  } // end of xmin block
    



  /** Define number of x values to be stored (one if all values of j_l1(x) were negligible for this l1) */  
  if (x_min >= pbs2->xx_max) {
    pbs2->index_xmin_l1[index_l1] = pbs2->xx_size - 1;
    pbs2->x_size_l1[index_l1] = 1;
  }
  else {
    pbs2->x_size_l1[index_l1] = pbs2->xx_size - pbs2->index_xmin_l1[index_l1];
  }


  class_test (((pbs2->index_xmin_l1[index_l1] < 0) || (pbs2->index_xmin_l1[index_l1] > pbs2->xx_size-1)),
              pbs2->error_message,
              "stopping to prevent segmentation fault");


  /** Allocate memory for j_l1[index_l1] */
  class_alloc (pbs2->j_l1[index_l1], sizeof(double)*pbs2->x_size_l1[index_l1], pbs2->error_message);
  #pragma omp atomic
  pbs2->count_allocated_Js += pbs2->x_size_l1[index_l1];

  if (ppr->bessels_interpolation == cubic_interpolation) {
    class_alloc (pbs2->ddj_l1[index_l1], sizeof(double)*pbs2->x_size_l1[index_l1], pbs2->error_message);
    #pragma omp atomic
    pbs2->count_allocated_Js += pbs2->x_size_l1[index_l1];
  }  


  /** Case when all values of j_l(x) were negligible for this l*/
  if (x_min >= pbs2->xx_max) {
    
    pbs2->x_min_l1[index_l1] = pbs2->xx_max;
    pbs2->j_l1[index_l1][0] = 0.;

  }
  /** Otherwise, write first non-negligible value and then loop over x */
  else {

    pbs2->x_min_l1[index_l1] = pbs2->xx[pbs2->index_xmin_l1[index_l1]];

    /* Loop over other non-negligible values */
    for (index_x=0; index_x < pbs2->x_size_l1[index_l1]; index_x++) {

      class_call (bessel_j(pbs,
                           l1,
                           pbs2->xx[pbs2->index_xmin_l1[index_l1] + index_x],
                           &pbs2->j_l1[index_l1][index_x]),
        pbs->error_message,
        pbs2->error_message);
         

      // *** Some debug
      // if (l1>2500) printf("l1=%d, index_x=%d, x=%10.5g, j=%10.5g\n", l1, index_x, pbs2->xx[pbs2->index_xmin_l1[index_l1] + index_x], pbs2->j_l1[index_l1][index_x]);

    }
  }
  
  return _SUCCESS_;

}




/**
 * Compute the convolution integral between a spherical Bessel function and up two two arrays.
 * This is the same function as bessel_convolution in the bessel module, except that the order
 * of the Bessel function has to belong to pbs2->l1, which is an extended version of pbs2->l. 
 */
int bessel2_convolution (
    struct precision * ppr,
    struct bessels2 * pbs2,
    double * kk,
    double * delta_kk,
    int k_size,
    double * f,
    double * g,
    int index_l,
    double r,
    double * integral,
    ErrorMsg error_message
    )
{

  /* Loop variable */
  int index_k;

  /* Initialize the integral */
  *integral = 0;

  /* Find the value of l from the Bessel structure */
  int l = pbs2->l1[index_l];
     
  /* We shall store the value of j_l2(r*k2) in here */
  double j;

  
  // *** Actual integration ***
    
  for (index_k = 0; index_k < k_size; ++index_k) {

    /* Value of the function f in k */
    double f_in_k = f[index_k];

    /* If the function f vanishes, do not bother computing the Bessel function,
    and jump to the next iteration without incrementing the integral.  Note that this
    is important because CLASS transfer functions at first-order are set to zero above
    a certain value of k. */
    if (f_in_k == 0.)
      continue;

    /* Same for the g function, if it is defined */
    double g_in_k = 1.;
    
    if (g != NULL) {

      g_in_k = g[index_k];

      if (g_in_k == 0.)
        continue;
    }

    /* Value of the considered k */
    double k = kk[index_k];

    
    // *** Bessel interpolation ***

    double x = k*r;

    /* j_l(x) vanishers for x < x_min(l) */
    if (x < pbs2->x_min_l1[index_l])
      continue;
    
    int index_x = (int)((x-pbs2->x_min_l1[index_l])/pbs2->xx_step);
    double a = (pbs2->x_min_l1[index_l]+pbs2->xx_step*(index_x+1) - x)/pbs2->xx_step;

    /* Store in 'j' the value of j_l(r*k) */
    if (ppr->bessels_interpolation == linear_interpolation) {
      
      j = a * pbs2->j_l1[index_l][index_x] 
          + (1.-a) * pbs2->j_l1[index_l][index_x+1];

    }
    else if (ppr->bessels_interpolation == cubic_interpolation) {

      j = a * pbs2->j_l1[index_l][index_x] 
          + (1.-a) * ( pbs2->j_l1[index_l][index_x+1]
            - a * ((a+1.) * pbs2->ddj_l1[index_l][index_x]
            +(2.-a) * pbs2->ddj_l1[index_l][index_x+1]) 
            * pbs2->xx_step * pbs2->xx_step / 6.0);
 
    }


    /* Value of the integrand function */
    double integrand = k * k * j * f_in_k * g_in_k;

    /* Increment the estimate of the integral */
    *integral += integrand * delta_kk[index_k];
    
  } // end of for(index_k)
   
   
  /* Divide the integral by a factor 1/2 to account for the trapezoidal rule */
  (*integral) *= 0.5;
   
  return _SUCCESS_;
  
} // end of bessel_convolution








