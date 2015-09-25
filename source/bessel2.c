/** @file bessel2.c
 * 
 * Module to compute and store the projection functions for SONG.
 *
 * The projection functions J are purely geometrical object needed to convert
 * the line of sight sources at recombination, which were computed by the
 * perturbations2.c module, into the transfer functions today, which will be
 * computed by the transfer2.c module.
 *
 * The procedure is detailed in sec. 5.5 of http://arxiv.org/abs/1405.2280,
 * and is a generalisation of the classical line of sight formalism initially
 * proposed in 1996 by Seljak and Zaldarriaga.
 *
 * The projection functions are linear combinations of spherical Bessel functions.
 * The temperature projection function (J_TT) is given in eq. 5.97 of
 * http://arxiv.org/abs/1405.2280, while the polarised ones (J_EE, J_EB)
 * are in eq. 5.103 and 5.104, respectively.
 *
 * The following functions can be called from other modules:
 *
 * -# bessel2_init() to run the module.
 * -# bessel2_free() to free the memory allocated by the module.
 * -# bessel2_convolution() to convolve the projection functions computed
 *    in this module with arbitrary functions.
 *
 * Created by Guido W. Pettinari on 17.03.2013 based on bessel.c by the
 * CLASS team (http://class-code.net/).
 * Last modified by Guido W. Pettinari on 29.08.2015
 */

#include "bessel2.h"



/**
 * Compute and store the J projection functions in pbs2->J_Llm_x.
 *
 * In detail, this function does:
 *
 * -# Determine the x-grid where J_Llm(x) will be sampled and store it in pbs2->xx.
 *    The list is the same for all (L,l,m), but we will sample each J_Llm(x) at a
 *    different point in the list (pbs2->index_xmin_J), for optimization purposes.
 * 
 * -# Determine which projection functions are needed based on the perturbations2.c
 *    module; the options are J_TT, J_EE and J_EB, and are shown and explained in eqs. 
 *    5.97, 5.101 and 5.102 of http://arxiv.org/abs/1405.2280.
 * 
 * -# Compute the spherical Bessel functions needed to build the projection functions,
 *    via bessel2_j_for_l1(). We label these as j_l1(x), to differentiate them from
 *    the j_l(x) computed in the bessel.c module.
 *
 * -# Allocate memory for the array of the projection functions pbs2->J_Llm_x and its
 *    second derivative pbs->ddJ_Llm_x(x), in view of spline interpolation.
 *
 * -# For each multipole value in pbs->l, call bessel_j_for_l() to compute and store
 *    the spherical Bessel functions.
 *
 * The computation of the Bessel function is determined by the following parameters:
 *
 * -# pbs->l[index_l]: list of size pbs->l_size with the multipole values l where we shall
 *    compute the projection functions; it is determined in the transfer.c module.
 * -# pbs2->xx_step: step dx for sampling the projection functions J(x), given by the
 *    user via the parameter file.
 * -# pbs2->xx_max: maximum value of x where to compute the projection functions; determined
 *    in the input2.c module
 * -# ppr2->bessel_J_cut_song: value of J_Llm(x) below which it is approximated by zero in the
 *    region x<<l
 *
 * This function requires the following modules:
 *
 * -# perturbations2_init(), to determine which projection functions are needed.
 * -# bessel_init(), to access the list of multipoles in pbs->l.
 */

int bessel2_init(
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2
    )
{

  // ==============================================================================
  // =                               Preparations                                 =
  // ==============================================================================

  /* Do we need to compute the 2nd-order projection functions? */
  if ((pbs->l_max==0) || ((ppt2->has_cls==_FALSE_) && (ppt2->has_cmb_bispectra==_FALSE_))
    || ((ppr2->load_transfers_from_disk==_TRUE_) && (ppr->load_bispectra_from_disk==_TRUE_))) {

    if (pbs2->bessels2_verbose > 0)
      printf("Second-order Bessel module skipped.\n");

    return _SUCCESS_;
  }
  else {
    if (pbs2->bessels2_verbose > 0)
      printf("Computing 2nd-order projection functions\n");
  }

  /* Initialize counter for the memory allocated in the projection functions */
  pbs2->count_allocated_Js = 0;

  /* Determine minimum allowed values of the Bessels and of the J's */
  pbs2->j_l1_cut  = ppr2->bessel_j_cut_song;
  pbs2->J_Llm_cut = ppr2->bessel_J_cut_song;

  /* Determine the grid in x where J_Llm(x) will be sampled */
  class_call (bessel2_get_xx_list (ppr, ppr2, ppt2, pbs, pbs2),
    pbs2->error_message,
    pbs2->error_message);

  /* Copy the array with the m values (azimuthal number of the spherical harmonics)
  from the ppr2 structure */
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


  
  // ====================================================================================
  // =                          Count projection functions                              =
  // ====================================================================================
  
  /* Determine which projection functions to compute, based on the required observables.
  In SONG, we consider three types of projection functions:
   - the temperature projection function (J_TT), which propagates the temperature
     perturbations in time;
   - the polarisation projection function (J_EE=J_BB), which propagates the E and B-modes
     perturbations in time;
   - the mixing projection function (J_EB=-J_BE), which turns E-modes into B-modes and
     viceversa; it exist only for non-scalar modes.
  For details on how these functions are defined and why they're important, refer to
  Sec. 5.5.1 of http://arxiv.org/abs/1405.2280. */
  
  int index_J = 0;
  
  if (ppt2->has_cmb_temperature == _TRUE_) {

    /* T->T projection function (eq. 5.97 of http://arxiv.org/abs/1405.2280) */
    pbs2->has_J_TT = _TRUE_;
    pbs2->index_J_TT = index_J++;
  }

  /* We compute the projection functions for both the E and B-modes, regardless of
  which polarisation type is requested, because free streaming makes the two types
  of polarisation mix in the line of sight integral */
  if ((ppt2->has_cmb_polarization_e == _TRUE_) || (ppt2->has_cmb_polarization_b == _TRUE_)) {

    /* E->E projection function (equal to B->B). See eq. 5.103 of
    http://arxiv.org/abs/1405.2280 */
    pbs2->has_J_EE = _TRUE_;
    pbs2->index_J_EE = index_J++;

    /* B->E projection function (equal to minus E->B). See eq. 5.103 of
    http://arxiv.org/abs/1405.2280. Since they vanish for scalar modes, we
    compute them only when we have at least a non-scalar mode */
    if (ppr2->m_max_2nd_order>0) {
      pbs2->has_J_EB = _TRUE_;
      pbs2->index_J_EB = index_J++;
    }
  }

  pbs2->J_size = index_J;

  if ((pbs2->bessels2_verbose > 0) && (ppr2->load_transfers_from_disk == _FALSE_))
    printf(" -> will compute size_J=%d projection functions\n", pbs2->J_size);



  // ====================================================================================
  // =                                  Compute j_l1(x)                                 =
  // ====================================================================================
      
  /* In order to compute the projection functions J_Llm(x), we need the spherical Bessel
  functions j_l1(x). For a given L and l, we will need j_l1(x) for all l1 between abs(l-L)
  and l+L. Here we build a list of required l1 values and store it in pbs2->. */
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
    printf (" -> computing spherical Bessel functions in l1\n");

  /* Compute j_l1(x) in a parallel loop over l1 */
  int abort = _FALSE_;
  #pragma omp parallel for schedule (dynamic)
  for (int index_l1 = 0; index_l1 < pbs2->l1_size; index_l1++) {

    class_call_parallel (bessel2_j_for_l1 (ppr,ppr2,pbs,pbs2,index_l1),
      pbs2->error_message,
      pbs2->error_message);

    #pragma omp flush(abort)

  }
  if (abort == _TRUE_) return _FAILURE_;

  /* Debug - Check the computation of the Bessel functions */
  // {
  //   double j;
  //   int index_l = 120;
  //   double x = 233.24142;
  //   bessel2_l1_at_x (pbs2, x, index_l, &j);
  //   printf ("j_%d(%g) = %g\n", pbs2->l1[index_l], x, j);
  //   bessel2_l1_at_x_linear (pbs2, x, index_l, &j);
  //   printf ("j_%d(%g) = %g\n", pbs2->l1[index_l], x, j);
  // }


  /* The projection functions are needed only to compute the transfer functions. If the
  latter are loaded from disk, there is no need for the former */
  if (ppr2->load_transfers_from_disk == _TRUE_) {
    if (pbs2->bessels2_verbose > 0)
      printf (" -> No second-order projection functions needed.\n");

    return _SUCCESS_;
  }



  // ====================================================================================
  // =                             Allocate result arrays                               =
  // ====================================================================================
  
  /* Allocate the (J,L,l,m) levels of the following arrays:

     - pbs2->x_size_J
     - pbs2->index_xmin_J
     - pbs2->x_min_J
     - pbs2->J_Llm_x
     - pbs2->ddJ_Llm_x
  
  We shall allocate the last level (x) of pbs2->J_Llm_x and pbs2->ddJ_Llm_x later,
  in bessel2_J_for_Llm(), because only then we will know for each (J,L,l,m) the domain
  where J_Llm(x) does not vanish. */

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



  // ====================================================================================
  // =                                    Compute J                                     =
  // ====================================================================================

  /* Compute the projection function J. We remember that J is a combination of spherical
  Bessel functions and 3j-symbols, and they are needed to compute the line of sight
  integral at second order. J is a three-index objects depending on the argument x.
  We shall denote it as J_Llm(x). There are three types of projection functions:
  temperature (J_TT), polarisation (J_EE) and polarisation mixing (J_EB). For more
  detail, see Sec. 5.5.1.4 of http://arxiv.org/abs/1405.2280. */
  
  /* Loop on the type of projection function */
  for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {
  
    if (pbs2->bessels2_verbose > 1)
      printf (" -> computing the projection function #%d\n", index_J);
  
    /* Loop on the source index L */
    for (int index_L = 0; index_L < pbs2->L_size; ++index_L) {
      
      /* Parallel loop on the multipole index l */
      abort = _FALSE_;
      #pragma omp parallel for schedule (dynamic)
      for (int index_l = 0; index_l < pbs->l_size; ++index_l) {
      
        /* There are two important consideration to take into account when dealing with
        the m!=0 case.  First, not all m-values are allowed: abs(m) should be smaller than
        MIN(L,l).  Secondly, there is no need to compute J_Llm(x) for negative m, because 
        it can be obtained by flipping the sign of the second line of one of the 3j, thus
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
      
        /* Loop on the azimuthal index m */
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

        } // end of for(index_m)
        #pragma omp flush(abort)    
      } // end of for(index_l)
      if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
    } // end of for(index_L)
  } // end of loop on type of projection functions
  
  /* Determine the maximum size of the x-level in pbs2->J_Llm_x */
  pbs2->x_size_max_J = 0;
  for (int index_J = 0; index_J < pbs2->J_size; ++index_J)
    for (int index_L = 0; index_L < pbs2->L_size; ++index_L)
      for (int index_l = 0; index_l < pbs->l_size; ++index_l)
        for (int index_m = 0; index_m < MIN(ppr2->index_m_max[pbs2->L[index_L]],ppr2->index_m_max[pbs->l[index_l]])+1; ++index_m)                 
          pbs2->x_size_max_J = MAX (pbs2->x_size_max_J, pbs2->x_size_J[index_J][index_L][index_l][index_m]);
  
  /* Print information on the used memory */
  if (pbs2->bessels2_verbose > 1)
    printf (" -> memory in use to store the 2nd-order projection functions: ~ %3g MB\n",
    8*pbs2->count_allocated_Js/1e6);



  // ====================================================================================
  // =                               Spline interpolation                               =
  // ====================================================================================
  
  if (ppr->bessels_interpolation == cubic_interpolation) {
    
    /* Compute second derivatives of the J_Llm in view of the spline interpolation. 
    The J will be used in the second-order transfer module to solve the line of sight
    integral. */
    
    for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {

      for (int index_L = 0; index_L < pbs2->L_size; ++index_L) {

        /* Beginning of parallel region */
        int abort = _FALSE_;    
        #pragma omp parallel for schedule (dynamic)
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
        if (abort == _TRUE_) return _FAILURE_;
      } // end of for(index_L)
    } // end of loop of projection function type


    /* Compute second derivatives of the j1_l in view of the spline interpolation. */

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
    }
    if (abort == _TRUE_) return _FAILURE_;


  } // end of spline calculation


  return _SUCCESS_;

}




/**
 * Compute the convolution integral between a spherical Bessel function and up two
 * two arrays.
 * 
 * Given the integration domain kk and the arrays f[index_kk] and g[index_kk], compute
 * the following integral using the trapezoidal rule:
 * 
 *     /
 *    |  dk k^2 f[k] * g[k] * j_l1[k*r]
 *    /
 *
 * The spherical Bessel functions are interpolated from the pbs2->j_l1 table in the
 * bessel2 structure, which contains more l1 values than the similar table in the
 * bessel structure (pbs->j_l).
 *
 * This is the same function as bessel_convolution() in the bessel module, except
 * that the order of the Bessel function has to belong to pbs2->l1, which is an extended
 * version of pbs2->l. 
 * 
 * The delta_k array should contain the trapezoidal measure around a given k[i],
 * that is delta_k[i]=k[i+1]-k[i-1], with delta_k[0]=k[1]-k[0] and
 * delta_k[k_size-1]=k[k_size-1]-k[k_size-2].
 *
 * If you give a NULL pointer for the g function, then it is not used at all.
 */

int bessel2_convolution (
    struct precision * ppr, /**< pointer to precision2 structure */
    struct bessels2 * pbs2, /**< pointer to Bessel2 structure, should be already initiated
                            with bessel2_init() */
    double * kk, /**< array with the integration grid in k */
    double * delta_kk, /**< trapezoidal measure, compute as delta_k[i]=k[i+1]-k[i-1],
                       and delta_k[0]=k[1]-k[0], delta_k[k_size-1]=k[k_size-1]-k[k_size-2] */
    int k_size, /**< size of the integration grid in k */
    double * f, /**< array with the integrand function f, of size k_size */
    double * g, /**< array with the integrand function g, of size k_size; pass NULL
                to automatically set it to unity */
    int index_l, /**< order of the Bessel function, taken from the multipole array pbs2->l1 */
    double r, /**< frequency of the Bessel function */
    double * integral, /**< output, estimate of the integral */
    ErrorMsg error_message /**< string to write error message */
    )
{
  
#ifdef DEBUG
  /* Test that the Bessel functions have been computed for the requested
  multipole index (index_l) and argument (x=k*r) */
  class_test ((index_l<0) || (index_l>=pbs2->l1_size),
    error_message,
    "index_l=%d out of bounds (l1_size=%d)", index_l, pbs2->l1_size);

  class_test ((kk[k_size-1]*r)>pbs2->xx[pbs2->xx_size-1],
    error_message,
    "r*k_max=%g is larger than x_max=%g (index_l1=%d,r=%g,k_max=%g)",
    kk[k_size-1]*r, pbs2->xx[pbs2->xx_size-1], index_l,r,kk[k_size-1]);
#endif // DEBUG

  /* Initialize the integral */
  *integral = 0;

  /* Find the value of l from the Bessel structure */
  int l = pbs2->l1[index_l];

  /* Loop over the integration grid */ 
  for (int index_k = 0; index_k < k_size; ++index_k) {

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

    /* - Bessel interpolation j_l(k*r) */
    double j;
    double x = k*r;

    /* j_l(x) vanishes for x < x_min(l) */
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



/** 
 * Compute the spherical Bessel function j_l1(x) using spline interpolation.
 *
 * Evaluates the spherical Bessel function x at a given value of x by interpolating
 * in the pre-computed table pbs2->j_l1. This function can be called from whatever
 * module at whatever time, provided that bessel2_init() has been called before,
 * and bessel2_free() has not been called yet.
 *
 * The index for l1 refers to the array pbs2->l1, which is an extension of ptr->l.
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





/** 
 * Compute the spherical Bessel function j_l1(x) using linear interpolation.
 *
 * Evaluates the spherical Bessel function x at a given value of x by interpolating
 * in the pre-computed table pbs2->j_l1. This function can be called from whatever
 * module at whatever time, provided that bessel2_init() has been called before,
 * and bessel2_free() has not been called yet.
 *
 * The index for l1 refers to the array pbs2->l1, which is an extension of ptr->l.
 */

int bessel2_l1_at_x_linear(
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

    /* Find result with linear interpolation (same as spline, but keep only
    the terms linear in a) */

    *j_l1 = a * pbs2->j_l1[index_l1][index_x] 
          + (1.-a) * pbs2->j_l1[index_l1][index_x+1];

  }
  
  return _SUCCESS_;

}



/**
 * Free all memory space allocated by bessel2_init().
 *
 * To be called at the end of each run.
 */

int bessel2_free( 
    struct precision * ppr,
    struct precision2 * ppr2,
    struct bessels * pbs,
    struct bessels2 * pbs2
    )
{

  if (ppr2->load_transfers_from_disk == _FALSE_) {

    for (int index_J = 0; index_J < pbs2->J_size; ++index_J) {
  
      for (int index_L = 0; index_L < pbs2->L_size; index_L++) {
  
        for (int index_l = 0; index_l < pbs->l_size; index_l++) {
  
          int index_m_max = MIN (ppr2->index_m_max[pbs2->L[index_L]], ppr2->index_m_max[pbs->l[index_l]]);
  
          for (int index_m = 0; index_m <= index_m_max; ++index_m) {
        
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

  }
  
  free(pbs2->l1);
  free(pbs2->xx);
  
  free(pbs2->x_size_l1);
  free(pbs2->x_min_l1);
  for (int index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
    free(pbs2->j_l1[index_l1]);
  free(pbs2->j_l1);
  if (ppr->bessels_interpolation == cubic_interpolation) {
    for (int index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
      free(pbs2->ddj_l1[index_l1]);
    free(pbs2->ddj_l1);
  }
  
  free(pbs2->L);
  free(pbs2->m);
  
  return _SUCCESS_; 
}


/**
 * Determine list of l1 multipoles where to compute the Bessel functions
 * j_l1(x).
 *
 * The j_l1(x) Bessels are needed to build the projection function J_Llm(x),
 * according to eq. 5.97 of http://arxiv.org/abs/1405.2280.
 *
 * The argument l1 ranges between l1_min=|l-L| and l1_max=l+L. Here we build
 * the list of all l1 values that will be ever needed. This is given by the
 * list of l (pbs->l) with the addition of extra points.
 *
 * This function will fill pbs2->l1 and pbs2->l1_size.
 */

int bessel2_get_l1_list(
          struct precision * ppr,
          struct precision2 * ppr2,
          struct bessels * pbs,
          struct bessels2 * pbs2
          )
{

  if (pbs2->bessels2_verbose > 2) {
    printf (" -> will compute %d l's in the range l=(%d,%d)\n", pbs->l_size, pbs->l[0], pbs->l[pbs->l_size-1]);    
    printf ("    * l-list:\n       ");

    for (int index_l=0; index_l < pbs->l_size; ++index_l)
      printf ("%d /\\ ", pbs->l[index_l]);
      
    printf ("\n");
  }

  
  // ====================================================================================
  // =                                 Count l1 values                                  =
  // ====================================================================================

  /* The line of sight integration at second-order requires to compute a linear
  combination of spherical Bessel functions with 3j-symbols (see eq. 5.97 of
  http://arxiv.org/abs/1405.2280). In particular, the Bessel function of order l1
  is required several times.  However, l1 goes from |l-L| to l+L, where the upper
  limit of L is given by pbs2->L_max. This means that we need to compute the spherical
  Bessels in an extended range that includes the l that are within a range of
  +/- pbs2->L_max from  each member of pbs2->l.
  
  Similarly, to compute the contribution to the bispectrum integral from a certain
  azimuthal number m, we will need to compute Bessel functions with order between
  |l-|m|| and l+|m| (see eq. 6.36 of http://arxiv.org/abs/1405.2280). Hence, l1
  needs to include all l that are within a range of
  +/- MAX(pbs2->L_max, ppr2->m_max_second_order) from each member of pbs2->l. */
    

  /* The maximum number of elements of pbs2->l1 is given by the maximum l in pbs->l
  plus MAX(pbs2->L_max, ppr2->m_max_2nd_order). We use a logical array to keep track
  of which l1 are to be kept. */
  int L_max = pbs2->L_max + ppr2->m_max_2nd_order; // +1;
  int l1_max = pbs->l[pbs->l_size - 1] + L_max;
  short * l1_logical;
  class_calloc (l1_logical, l1_max+1, sizeof(int), pbs2->error_message);

  /* Each of the l1's is to be kept only if it is in the range abs(l-L_max)<=l1<=l+L_max */
  pbs2->l1_size = 0;

  for(int l1=0; l1<(l1_max+1); ++l1) {

    l1_logical[l1] = _FALSE_;                   // Initialize flag

    for(int index_l=0; index_l<pbs->l_size; ++index_l) {

      int l = pbs->l[index_l];

      if ((l1>=(l-L_max)) && (l1<=l+L_max))     // It is important to omit abs()
        l1_logical[l1] = _TRUE_;                // Flag the l1 to be kept
                  
    }
    
    if (l1_logical[l1] == _TRUE_)               // Increment the counter
      pbs2->l1_size++;

  }

  /* Debug - Print which l1 shall be computed */
  // for(index_l1=0; index_l1<(l1_max+1); ++index_l1)
  //   if (l1_logical[index_l1] == _FALSE_)
  //     printf("We shall not compute the l1=%d bessel.\n", index_l1);
  // printf("pbs->l_size = %d, pbs2->l1_size = %d\n", pbs->l_size, pbs2->l1_size);
  
  

  // ====================================================================================
  // =                                 Write l1 values                                  =
  // ====================================================================================

  class_alloc (pbs2->l1, pbs2->l1_size*sizeof(int), pbs2->error_message);

  int l1 = 0;

  for(int index_l1=0; index_l1<pbs2->l1_size; ++index_l1) {

    while (l1_logical[l1] == _FALSE_)
      l1++;                         // Look for the first l1 to keep

    pbs2->l1[index_l1] = l1++;      // and store it in pbs2->l1
    
  }

  free(l1_logical);


  /* It is useful to define an array that contains for each value of l1
  the index that value is assigned in pbs->l1. */

  class_alloc (pbs2->index_l1, (pbs2->l1[pbs2->l1_size-1]+1)*sizeof(int), pbs2->error_message);

  for(l1=0; l1<=pbs2->l1[pbs2->l1_size-1]; ++l1) {
  
    pbs2->index_l1[l1] = -1;       // Assign a negative value if l1 is not included in pbs2->l1
  
    for (int index_l1=0; index_l1<pbs2->l1_size; ++index_l1)
      if (l1==pbs2->l1[index_l1]) pbs2->index_l1[l1] = index_l1;

  }

  return _SUCCESS_;

}





/** 
 * Determine the x-grid where J_Llm(x) will be sampled and store it in pbs2->xx.
 *
 * Since the projection functions J_Llm(x) are a linear combination of spherical Bessel
 * functions j_l1(x), te sample them linearly.
 *
 * The step and maximum value of the linear sampling are determined in the input2.c
 * module based on the user-given parameters ppr2->k_max_tau0_over_l_max,
 * ppt->l_scalar_max and ppr2->bessel_x_step_song.
 * 
 */
int bessel2_get_xx_list(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct bessels * pbs,
      struct bessels2 * pbs2
      )
{
    
  /* The argument of the projection function in the line of sight integral is x=k*(tau_0-tau).
  Therefore, we choose the maximum x-value where to compute J, pbs2->xx_max, to be as large
  as the maximum k used in SONG, times the tau_0, the conformal age of the Universe. This is 
  done in the input2.c module, where we have set
    pbs2->xx_max = ppr2->k_max_tau0_over_l_max * ppt->l_scalar_max.
  Here we set pbs2->xx_max to be at least as large as the maximum x used to compute the Bessel
  functions of the bessel.c module. */
  pbs2->xx_max = MAX (pbs2->xx_max, pbs->x_max);
  
  /* Sample xx linearly between 0 and pbs2->xx_max with step pbs2->xx_step */
  pbs2->xx_size = pbs2->xx_max/pbs2->xx_step + 1;
  class_alloc (pbs2->xx, pbs2->xx_size*sizeof(double), pbs2->error_message);
  lin_space (pbs2->xx, 0, pbs2->xx_max, pbs2->xx_size);
  
  return _SUCCESS_;

}





/**
 * Compute the projection function J_Llm(x) for all values of x where it is
 * non-negligible, given the source index L, the multipole index l and the
 * azimuthal index m.
 * 
 * The projection functions will be convolved with the second-order source
 * functions in the transfer2.c module, as described in Sec. 5.5.1 of
 * http://arxiv.org/abs/1405.2280. For reference, this is the physical meaning
 * of the involved parameters, assuming all sources are located at recombination;
 * 
 * - L is the multipole index of the source function that will be convolved with
 *   the projection function. The most important contribution comes from the
 *   L=0,1,2 multipoles, because the others are tight-coupling suppressed.
 *
 * - l and m are the angular scale and distribution, respectively, of the
 *   photon perturbations as seen today.
 *
 * - x is the argument of the projection function; in the line of sight integral
 *   it is equal to k(tau-tau0). Since the projection function shares the same
 *   geometrical properties of a Bessel function, this means that the largest
 *   contribution today at the angular scale l comes from those modes that have
 *   k~l/tau_0.
 * 
 * In detail, this function does:
 *
 * -# Compute the 3j symbols that will be needed for the computation of J(x).
 * -# Using bisection, find the first point in pbs2->xx where J(x) is non-negligible,
 *    and the number of such points.
 * -# Compute J(x) only in the points where it is not negligible, using the function
 *    bessel2_J_Llm().
 * 
 * The following arrays will be filled:
 * - pbs2->J_Llm_x[index_J][index_L][index_l][index_m][index_x], containing the
 *   projection function J_Llm(x).
 * - pbs2->index_xmin_J[index_J][index_L][index_l][index_m], the first index in pbs2->xx
 *   where J(x) is non-negligible.
 * - pbs2->xmin_J[index_J][index_L][index_l][index_m], the first point in pbs2->xx where
 *   J(x) is non-negligible.
 * - pbs2->x_size_J[index_J][index_L][index_l][index_m], the number of points in pbs2->xx
 *   where J(x) is non-negligible.
 *
 * J(x) will be computed only in the values of pbs2->xx where it is non-negligible. The
 * number of such values depends on the considered (L,l,m) configuration. In general,
 *  x = x_min_J(L,l,m) + xx_step*index_x
 * where x_min_J(L,l,m) is determined in this function through bisection, xx_step comes
 * from the parameter file (usually 0.2), and index_x is the index of the last level of
 * pbs2->J_Llm_x.
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

  /* Shortcuts */
  int L = pbs2->L[index_L];
  int l = pbs->l[index_l];
  int m = pbs2->m[index_m];

  /* Index in our x domain (pbs2->xx) from which we shall start to sample J(x).
  If set to zero, then we shall consider all points, if set to pbs2->xx_size-1,
  we shall set J(x) to vanish for all x */
  int index_xmin_J = 0;


  // ====================================================================================
  // =                               Determine which J                                  =
  // ====================================================================================

  /* Determine the requested projection function, and its spin. Temperature has spin zero,
  while E and B-modes have spin 2. See Sec. 5.5.1.4 of http://arxiv.org/abs/1405.2280, or
  appendix B of http://arxiv.org/abs/1102.1524, for further details. */
  enum projection_function_types projection_function;
  int S;
  
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

  /* The spin S appears in the definition of the projection functions in the second line
  of the 3j-symbol, below both l and L. This means that S should be smaller of both l and
  L, a constraint that we enforce here. */
  short spin_constraint = ((abs(S)>l) || (abs(S)>L));

  /* If m=0, the projection function for the E-B mixing vanishes */
  short EB_constraint = ((projection_function==J_EB) && (m==0));

  /* If any of the above condition is met, the projection functions J vanish. We enforce
  this by telling SONG that J(x) is negligible for all x-values in pbs2->xx. */  
  if (spin_constraint || EB_constraint) {
    index_xmin_J = pbs2->xx_size - 1;
    goto end_of_bisection; 
  }

  if (pbs2->bessels2_verbose > 2)
    printf("    * processing (index_J,index_L,index_l,index_m) = (%d,%d,%d,%d), (L,l,m) = (%d,%d,%d)\n",
    index_J, index_L, index_l, index_m, L, l, m);



  // ====================================================================================
  // =                               Compute 3j symbols                                 =
  // ====================================================================================

  /* The projection function requires the computation of two Wigner 3j symbols. We shall
  store the results of the 3j computation inside the simple ad-hoc structure
  bessel_3j_data */
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
  for all allowed values of l1. Since l1 has to be first, we rearrange it as
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
  for all allowed values of l1. Since l1 has to be first, we rearrange it as
    (   l1      L      l   )
    (    0     -m      m   )
  We compute the second 3j only if it is different from the first one.
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



  // ====================================================================================
  // =                            Exclude negligible values                             =
  // ====================================================================================
  
  /* The projection function J_Llm(x) is basically a Bessel function of order l, and as
  such it vanishes when x is much smaller than l. Here we determine the value of x after
  which J start to be non-negligible, using bisection; we denote it as x_min_J, and we
  denote the position of x inside pbs2->xx as index_xmin_J. The parameter we use to
  determine the non-negligible threshold is given directly from the parameter file:
  pbs2->J_Llm_cut=ppr2->bessel_J_cut_song. It is usually of order 1e-6. */

  double x_min_J;

  
  /* Special case: when L==l, we have a contribution to J_Llm(x) coming from
  j_0(x) = sin(x)/x. Unlike all the other spherical Bessels functions, j_0(x) does
  not vanish for x->0 but it asymptotes to 1. Ther Hence, there is no regime close to
  x=0 where J_llm(x) is negligible, and when L=l we just take x_min=0 */

  if (L==l) {
    index_xmin_J = 0;
    x_min_J = pbs2->xx[index_xmin_J];
    goto end_of_bisection;
  }

  // -------------------------------------------------------------------------------
  // -                          Determine bisection limits                         -
  // -------------------------------------------------------------------------------

  /* The lower limit for the bisection is the lowest value in the x domain (pbs2->xx) */
  int index_xmin_down = 0;
  double x_min_down = pbs2->xx[index_xmin_down];    

  /* We determine the upper limit for the bisection using the properties of the Bessel
  functions. The projection function J_Llm(x) involves a sum over spherical Bessel
  functions with order |l-L|<=l1<=l+L. The spherical Bessel function j_l1(x) grows as x^l1
  until it peaks at x~l1+0.5. Afterwards it starts to oscillate with an envelope
  that goes as 1/x. We shall set the upper limit for the bisection, x_min_up, at the first peak,
  that is, we shall set x_min_up = |l-L|+0.5. In this way, we ensure that we start after the
  regime where J_Llm is strongly suppressed (x<l1_min) but before it starts oscillating 
  (x>l1_min). */
  
  double x_min_up = bessel_3j_data->l1_min + 0.5;
  int index_xmin_up;

  /* Find the index in pbs2->xx corresponding to the bisection upper limit */
  
  if (x_min_up >= pbs2->xx_max) {

    /* When our guess is larger than the maximum value of x in pbs2->xx, it means that 
    none of the terms in J have their peak in our x domain. It is then likely that the
    whole projection function J is negligible. However, we still perform the bisection
    out of caution, setting the upper limit to the largest point in the x-domain */
    index_xmin_up = pbs2->xx_size-1;
    x_min_up = pbs2->xx_max;
  }
  else {

    /* Find the index in pbs2->xx that follows x_min_up */
    index_xmin_up = 0;
    while(pbs2->xx[index_xmin_up] < x_min_up)
       index_xmin_up++;
    
    /* Update the value of x_min_up to reflect an actual entry of the xx array */
    x_min_up = pbs2->xx[index_xmin_up];
  }
  
  /* Compute the value of J_Llm in the upper limit */  
  double J = 0;
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

  /* Debug - Print this first step of the bisection */
  // if ((L==4) && (l==10))
  //   printf("bisection: index_xmin_down=%d, index_xmin_up=%d\n", index_xmin_down, index_xmin_up);
  // if ((L==4) && (l==10))
  //   printf("bisection: x_min_down=%g, x_min_up=%g, J=%g, J_Llm_cut=%g\n", x_min_down, x_min_up, J, pbs2->J_Llm_cut);

  /* If the value of J_Llm in x_min_up is already below the target cut, then it means
  that we already hit a point where the J_Llm is negligible. We can then safely flag
  J_Llm(x) to be negligible for all the values in pbs2->xx. */

  if (fabs(J) < pbs2->J_Llm_cut) {

    /* J is negligible in all points of xx */
    index_xmin_J = pbs2->xx_size - 1;

    /* Check on our Bessel functions. If J_Llm(x_min_up) is negligible but
    xx_max>x_min_up=l_min+0.5, it means that either our initial guess was way too high
    and the Bessels are all very damped (x_min_up >> l1_min+0.5), or that we entered into
    the oscillatory regime of J_Llm and hit an x-value close to a zero crossing. In both
    cases, we return an error. */

    class_test(x_min_up < pbs2->xx_max,
      pbs2->error_message,
      "in bisection, wrong initial guess for x_min_up");
      
    goto end_of_bisection;
      
  }


  // -----------------------------------------------------------------------------
  // -                              Actual bisection                             -
  // -----------------------------------------------------------------------------

  /* Stop the bisection only when we have unambiguously bracketed the first point in
  pbs2->xx whereby |J_Llm(x)| < pbs2->J_Llm_cut  */

  while (index_xmin_down != (index_xmin_up-1)) {

    /* Find the middle point between the lower and upper limits on x, by using the
    fact that pbs2->xx is sampled linearly */
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

    /* Debug - Check the current bisection step */
    // if ((L==4) && (l==10))
    //   printf("bisection: index_xmin_down=%d, index_xmin_up=%d\n", index_xmin_down, index_xmin_up);
    // if ((L==4) && (l==10))
    //   printf("bisection: x_min_down=%g, x_min_up=%g, J=%g, J_Llm_cut=%g\n", x_min_down, x_min_up, J, pbs2->J_Llm_cut);

  }

  /* Store the bisection results */
  index_xmin_J = index_xmin_down;
  x_min_J = x_min_down;

  /* Debug - Print the number of non-negligible x for J */
  // printf("pbs2->x_size_J[%d][%d][%d][%d] = %d, pbs->x_max = %g\n",
  //   index_J, index_L, index_l, index_m, pbs2->x_size_J[index_J][index_L][index_l][index_m], pbs->x_max);

  end_of_bisection:



  // ====================================================================================
  // =                                 Compute and store J                              =
  // ====================================================================================

  /* Define number of x values to be stored (one if all values of j_l(x) were negligible for this l) */  
  pbs2->x_size_J[index_J][index_L][index_l][index_m] = pbs2->xx_size - index_xmin_J;

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
  
  /* Define a shorthand for J_Llm_x */
  double * J_Llm_x = pbs2->J_Llm_x[index_J][index_L][index_l][index_m];
  
  /* In case J(x) was found to be negligible for all the x values in our domain,
  set the projection function to zero */
  
  if (index_xmin_J == (pbs2->xx_size-1)) {
  
    pbs2->index_xmin_J[index_J][index_L][index_l][index_m] = pbs2->xx_size - 1;
    pbs2->x_min_J[index_J][index_L][index_l][index_m] = pbs2->xx_max;
    J_Llm_x[0] = 0.;

    if (pbs2->bessels2_verbose > 3)
      printf("      \\ Set J(L=%d,l=%d,m=%d)=0 for all x (xx_max=%g)\n",
        L, l, m, pbs2->xx_max);
          
  }

  /* Otherwise, compute the projection function for all the x-values where it is
  not negligible */

  else {
  
    pbs2->index_xmin_J[index_J][index_L][index_l][index_m] = index_xmin_J;
    pbs2->x_min_J[index_J][index_L][index_l][index_m] = x_min_J;
    
    /* Loop over all other non-negligible values */
    
    for (int index_x=index_xmin_J; index_x < pbs2->xx_size; index_x++) {
    
      /* Note that the first index in J_Llm_x corresponds to the first non-negligible value,
      rather than to the first value of x inside pbs2->xx. */
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
    
      /* Debug - Print J_Llm(x) to screen */
      // if ((L==0) && (l==100) && (m==0))
      //   fprintf(stderr, "%g %g\n", pbs2->xx[index_x], J_Llm_x[index_x - index_xmin_J]);
    
    }
  }


  /* Free 3j-symbol arrays */
  if (!spin_constraint && !EB_constraint) {
    free(bessel_3j_data->bessels);
    free(bessel_3j_data->first_3j);
    free(bessel_3j_data->second_3j);
    free(bessel_3j_data);
  }
  
  return _SUCCESS_;
  
}






/**
 * Compute the projection function J_Llm(x) using the precomputed table of spherical
 * Bessel functions (pbs2->j_l1) and output the result in J_Llm_x.
 *
 * The projection functions are linear combinations of spherical Bessel functions.
 * The temperature projection function (J_TT) is given in eq. 5.97 of
 * http://arxiv.org/abs/1405.2280, while the polarised ones (J_EE, J_EB)
 * are in eq. 5.103 and 5.104, respectively.
 * 
 * This function relies on the 3j symbols appearing in the definition of J(x)
 * to be precomputed and ready in the bessel_3j_data structure.
 * 
 * Make sure that the input l belongs to pbs->l, otherwise you will get segmentation
 * faults.
 */

int bessel2_J_Llm (
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       enum projection_function_types projection_function,
       int L, /**< input, source multipole */
       int l, /**< input, angular scale today */
       int m, /**< input, angular distribution today */
       int index_x, /**< input, argument of the projection function J(x) */
       struct J_Llm_data * bessel_3j_data, /**< input, data structure with 3j symbols ready */
       double * J_Llm_x /**< output, the computed value of J_Llm(x) */
       )
{

  /* Initialise the output value; we shall accumulate it in a loop over l1 */
  (*J_Llm_x) = 0;

  /* Loop over the l1. We are considering only the l1 values allowed by the 3j symbol
  symmetries, that is |l-L|<=l1<=l+L. Note that this index_l1 IS NOT the one that should
  be used to access pbs2->l1, which instead requires bessel_3j_data->index_l1_min+index_l1. */

  for(int index_l1=0; index_l1<bessel_3j_data->l1_size; ++index_l1) {


    // -----------------------------------------------------------------------------
    // -                        Skip vanishing contributions                       -
    // -----------------------------------------------------------------------------

    /* Skip the contribution when the three-j vanishes */
    if (bessel_3j_data->first_3j[index_l1]==0)
      continue;

    /* Obtain the value of l1 */
    int index_l1_in_jl1 = bessel_3j_data->index_l1_min + index_l1;
    int index_x_in_jl1  = index_x - pbs2->index_xmin_l1[index_l1_in_jl1];
    int l1 = pbs2->l1[index_l1_in_jl1];
    
    /* If the Bessel function that is being summed is smaller than the treshold set in 
    ppr2->bessel_j_cut_song, then skip it */
    if (index_x_in_jl1 < 0) {
      if (pbs2->bessels2_verbose > 2)
        printf("     \\ Assumed that the l1=%d contribution to J_%d_%d_%d(%g) vanishes (i.e. j_%d(%g) < %g)\n",
          l1, L, l, m, pbs2->xx[index_x], l1, pbs2->xx[index_x], pbs2->j_l1_cut);      
      continue;
    }

    /* Use the parity properties of the projection function to skip vanishing contributions.
    J_TT, J_EE and J_BB are different from zero only for even l+l1+L. For J_TT, this is
    enforced by the three-j symbol. For J_EE and J_BB, this is enforced by the parity of
    the mixing matrix H (see eq. 76 of Beneke & Fidler 2010). The mixing matrix also enforces
    that J_EB and J_BE are different from zero only for odd l+l1+L.
    
    The prefactor is i^(l-l1-L) for the even parity projection functions (J_TT, J_EE, J_BB),
    and either i^(l-l1-L-1) or i^(l-l1-L+1) for the odd parity ones (J_EB and J_BE respectively).
    In both cases, the parity forces the prefactors to be just (real) alternating signs. */

    short is_even = ((l-l1-L)%2==0);
    double i_prefactor;

    if ((projection_function==J_TT)||(projection_function==J_EE)) {
      if (is_even)
        i_prefactor = ALTERNATING_SIGN((l-l1-L)/2);
      else
        continue;
    }
    else if (projection_function==J_EB) {
      if (!is_even)
        i_prefactor = ALTERNATING_SIGN((l-l1-L-1)/2);
      else
        continue;
    }


    // -----------------------------------------------------------------------------
    // -                           Build l1 contribution                           -
    // -----------------------------------------------------------------------------
    
    /*  Value of the spherical Bessel function j_l(x), precomputed in pbs2->j_l1[index_l1][index_x].
    The l1 level of pbs2->j_l1 should be addressed the same way as pbs2->l1, while the x level should
    be addressed as index_x - pbs2->index_xmin_l1[index_l1], where index_x is the index of pbs2->xx. */
    double j_l1_x = pbs2->j_l1[index_l1_in_jl1][index_x_in_jl1];
    
    /* Increment the result */
    (*J_Llm_x) += i_prefactor * (2*l1+1)
                  * bessel_3j_data->first_3j[index_l1]
                  * bessel_3j_data->second_3j[index_l1]
                  * j_l1_x;

    /* Debug - Print Bessel function and 3j symbol */
    // if ((L==4) && (l==1172)) printf("j(%d, %10.6g) = %10.6g\n", l1, x, j_l1_x);
    // if ((L==4) && (l==1172)) printf("3j(%d,%d,%d)(%d,0,%d) = %10.6g\n", l1, L, l, -m, m, bessel_3j_data->second_3j[index_l1]);    
    
  } // end of for(index_l1)

  /* Apply l1-independent factors */
  (*J_Llm_x) *= ALTERNATING_SIGN(m) * (2*l+1);


  /* Debug - Print projection function */
  // if ((projection_function==J_EB) && (L==2) && (l==102) && (m==1))
  //   fprintf(stderr, "%10.6g %10.6g\n", x, *J_Llm_x);
  //   fprintf(stderr, "J(%d,%d,%d,%10.6g[%d]) = %10.6g\n", L, l, m, x, index_x, *J_Llm_x);
  

  /* Mathematica batch debug */
  // double val = *J_Llm_x;
  // if (projection_function == J_TT) {
  //   if (val==0.) {
  //     fprintf (stderr, "Abs[J[%d,%d,%d][%g]]<10^%g,\n",
  //     L,l,m,x,log10(ppr2->bessel_J_cut_song));
  //   }
  // }
  // double val = *J_Llm_x;
  // if (projection_function == J_TT) {
  //   if (val!=0) {
  //     if ((l<250) && (l<250)) {
  //       if (m==1) {
  //         if ((index_x<3) || (index_x%128==0)) {
  //           fprintf (stderr, "Abs[1 - J[%d,%d,%d][%g]/((%g)*10^(%f))] < 10^-2,\n",
  //             L,l,m,x,SIGN(val),log10(fabs(val)));
  //         }
  //       }
  //     }
  //   }
  // }
  
  return _SUCCESS_;

}



/**
 * Compute the projection function J_Llm(x) without relying on the precomputed
 * spherical Bessel in pbs2->j_l1.
 *
 * To do so, use the numerical recipes subroutine bessel_j. This function might
 * be useful to debug bessel_J_Llm, but it is never used in SONG.
 * 
 * This function only supports the temperature projection function (J_TT);
 * extending it to J_EE and J_EB is straightforward.
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
 
  double j_l1_x = 0;         // Bessel function j_l(x) 
  (*J_Llm_x) = 0;            // We shall accumulate J over l1

  /* Summing up to l1_size is equivalent to sum over all the allowed l1 */
  for(int index_l1=0; index_l1<bessel_3j_data->l1_size; ++index_l1) {
  
    /* The first 3j-symbol vanishes for odd l-l1-L */
    if (bessel_3j_data->first_3j[index_l1]==0.) continue;
    
    /* Value of l1 */
    int l1 = bessel_3j_data->l1_min + index_l1;
    
    /*  Due to the symmetries of the first 3j-symbol, we need to consider only
    even l-l1-L. Hence, the coefficient i^(l-l1-L) depends on the parity of
    of (l-l1-L)/2.  */
    double i_coefficient = ALTERNATING_SIGN((l-l1-L)/2);
    
    /* Value of the spherical Bessel function j_l(x) */
    class_call (bessel_j(pbs,
                  l1,
                  x,
                  &j_l1_x),
      pbs2->error_message,
      pbs2->error_message);

    /* Increment the result */
    (*J_Llm_x) += i_coefficient * (2*l1+1) * bessel_3j_data->first_3j[index_l1]
                  * bessel_3j_data->second_3j[index_l1] * j_l1_x;
    
  } // end of for(index_l1)

  /* Apply constant factors */
  (*J_Llm_x) *= ALTERNATING_SIGN(m) * (2*l+1);

  /* Debug - print the projection function */
  // printf("J(%d,%d,%d,%20.15g) = %.15g\n", L, l, m, x, *J_Llm_x);
  
  return _SUCCESS_;
  
}




/**
 * Compute the spherical Bessel functions for the given value of l1, which must belong
 * to the array pbs2->l1, and store the result in pbs2->j_l1[index_l1].
 *
 * This function it is equivalent to bessel_j_for_l(), from CLASS. The only
 * differences is that it computes the spherical Bessels for the extended list
 * pbs2->l1 rather than for pbs->l; see bessel2_get_l1_list() for more details
 * on the l1 list.
 */

int bessel2_j_for_l1(
       struct precision * ppr,
       struct precision2 * ppr2,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       int index_l1 /**< input, index inside pbs2->l1 for which to compute the Bessel function */
       )
{
  
  int l1 = pbs2->l1[index_l1];

  // ====================================================================================
  // =                            Exclude negligible values                             =
  // ====================================================================================

  /* First value of x whereby the j_l1(x) is not negligible. Will be determined
  by bisection below. */
  double x_min;

  /* j_0(x) = sin(x)/x does not vanish for x->0, but it asymptotes to 1.
  Hence, there is no regime close to x=0 where J_llm(x) is negligible,
  and we just take x_min=0 when l1=0. */
  if (l1 == 0) {
    x_min = 0;
    pbs2->index_xmin_l1[index_l1] = 0;
  }
  else {

    /* Find x_min by bisection. As lower limit, we set the smallest point in the x domain,
    x=0, while for the the upper limit we choose the peak of the spherical Bessel function,
    x=l1+0.5. In this way, we are sure that we start after the area where j_l1(x) 
    is suppressed (x<l1) but before it starts oscillating (x>l1) */
    
    double x_min_down=0;    
    double x_min_up=(double)pbs2->l1[index_l1]+0.5;
    double j;

    class_call (bessel_j (pbs,
            l1,
            x_min_up,
            &j),
         pbs->error_message,
         pbs2->error_message);
  
    class_test(j < pbs2->j_l1_cut,
      pbs2->error_message,
      "in bisection, wrong initial guess for x_min_up.");
 
    while ((x_min_up-x_min_down)/x_min_down > ppr->bessel_tol_x_min) {

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
    } 
      
  } // end of if (l1!=0)
    

  /* Define number of x values to be stored; set it to unity if all values of j_l1(x) were
  negligible for this l1. */  
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


  /* Allocate memory for j_l1[index_l1] */

  class_alloc (pbs2->j_l1[index_l1], sizeof(double)*pbs2->x_size_l1[index_l1], pbs2->error_message);
  #pragma omp atomic
  pbs2->count_allocated_Js += pbs2->x_size_l1[index_l1];

  if (ppr->bessels_interpolation == cubic_interpolation) {
    class_alloc (pbs2->ddj_l1[index_l1], sizeof(double)*pbs2->x_size_l1[index_l1], pbs2->error_message);
    #pragma omp atomic
    pbs2->count_allocated_Js += pbs2->x_size_l1[index_l1];
  }  


  // ====================================================================================
  // =                             Compute and store j_l1                               =
  // ====================================================================================

  /* Case when all values of j_l(x) were negligible for this l*/

  if (x_min >= pbs2->xx_max) {
    
    pbs2->x_min_l1[index_l1] = pbs2->xx_max;
    pbs2->j_l1[index_l1][0] = 0.;

  }

  /* Otherwise, write first non-negligible value and then loop over x */

  else {

    pbs2->x_min_l1[index_l1] = pbs2->xx[pbs2->index_xmin_l1[index_l1]];

    /* Loop over other non-negligible values */
    for (int index_x=0; index_x < pbs2->x_size_l1[index_l1]; index_x++) {

      class_call (bessel_j(pbs,
                           l1,
                           pbs2->xx[pbs2->index_xmin_l1[index_l1] + index_x],
                           &pbs2->j_l1[index_l1][index_x]),
        pbs->error_message,
        pbs2->error_message);
    }
  }
  
  return _SUCCESS_;

}









