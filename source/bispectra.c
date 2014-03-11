/** @file bispectra.c documented spectra module for second-order perturbations
 *
 * Guido W Pettinari, 19.07.2012.
 */

#include "bispectra.h"


/**
 * This routine initializes the spectra structure (in particular, 
 * computes table of anisotropy and Fourier spectra \f$ C_l^{X}, P(k), ... \f$)
 * 
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfer structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Output: pointer to initialized spectra structure
 * @return the error status
 */

int bispectra_init (
     struct precision * ppr,
     struct background * pba,
     struct thermo * pth,
     struct perturbs * ppt,
     struct bessels * pbs,
     struct transfers * ptr,
     struct primordial * ppm,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi
     )
{

  /* Check whether we need to compute bispectra at all */  
  if (pbi->has_bispectra == _FALSE_) {

    if (pbi->bispectra_verbose > 0)
      printf("No bispectra requested. Bispectra module skipped.\n");

    pbi->bt_size=0;

    return _SUCCESS_;
  }
  else {
    if (pbi->bispectra_verbose > 0)
      printf("Computing bispectra\n");
  }


  // =========================================================
  // =                    Preparations                       =
  // =========================================================

  /* Initialize indices & arrays in the bispectra structure */

  class_call (bispectra_indices (ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);





  // =======================================================
  // =                  Compute bispectra                  =
  // =======================================================
  
  class_call (bispectra_harmonic (ppr,pba,ppt,pbs,ptr,ppm,psp,pbi),
    pbi->error_message,
    pbi->error_message);



	/* Debug - Print some bispectra configuration (make sure you're not loading the bispectrum from disk, though) */
  // for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {
  //   for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
  //     for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  //       int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
  //       int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  //       for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
  //         long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  //         int l1 = pbi->l[index_l1]; int l2 = pbi->l[index_l2]; int l3 = pbi->l[index_l3];
  //         double cl_1 = pbi->cls[psp->index_ct_tt][l1-2];;
  //         double cl_2 = pbi->cls[psp->index_ct_tt][l2-2];;
  //         double cl_3 = pbi->cls[psp->index_ct_tt][l3-2];;
  //         double squeezed_normalisation = 1/(12*cl_1*cl_3);
  //         double equilateral_normalisation = 1e16 * l1*l1 * (l1+1)*(l1+1) / ((2*_PI_)*(2*_PI_));
  //         if (strstr(pbi->bt_labels[index_bt], "equilateral")!=NULL) {
  //           /* Squeezed configuration */
  //           if ((l3==6) && (l1==l2)) {
  //             fprintf (stderr, "%12d %12g %12g %12g %12g %12g %12g %12g %12g\n",
  //               l1,
  //               pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][1][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][1][0][index_l1_l2_l3] * equilateral_normalisation
  //             );
  //           }
  //           /* Equilateral configuration */
  //           // if ((l1==l2) && (l2==l3)) {
  //           //   fprintf (stderr, "%12d %12g %12g %12g\n",
  //           //   l1, pbi->bispectra[index_bt][index_l1_l2_l3], equilateral_normalisation,
  //           //   pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3]*equilateral_normalisation);
  //           // }
  //         }
  //       } // end of for(index_l3)
  //     } // end of for(index_l2)
  //   } // end of for(index_l1)
  // } // end of for(index_bt)

  return _SUCCESS_;

}






/**
 * This routine frees all the memory space allocated by bispectra_init().
 *
 * @param pbi Input: pointer to bispectra structure (which fields must be freed)
 * @return the error status
 */

int bispectra_free(
     struct perturbs * ppt,
     struct spectra * psp,
     struct bispectra * pbi
     )
{

  if (pbi->has_bispectra == _TRUE_) {

    free(pbi->l);
    free(pbi->pk);

    for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {

      free(pbi->l_triangular_size[index_l1]);
      free(pbi->index_l_triangular_min[index_l1]);
      free(pbi->index_l_triangular_max[index_l1]);
    
    } // end of for(index_l1)

    free(pbi->l_triangular_size);
    free(pbi->index_l_triangular_min);
    free(pbi->index_l_triangular_max);


    /* Free pbi->bispectra */
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      for (int X = 0; X < pbi->bf_size; ++X) {
        for (int Y = 0; Y < pbi->bf_size; ++Y) {
          for (int Z = 0; Z < pbi->bf_size; ++Z)
            free (pbi->bispectra[index_bt][X][Y][Z]);
          free (pbi->bispectra[index_bt][X][Y]);
        }
        free (pbi->bispectra[index_bt][X]);
      }
      free (pbi->bispectra[index_bt]);
    }
    free (pbi->bispectra);


    /* Arrays specific to the primordial models */
  
    if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
        
      for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
        for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
          free (pbi->alpha[index_bf][index_l]);
          free (pbi->beta[index_bf][index_l]);
        }
        free (pbi->alpha[index_bf]);
        free (pbi->beta[index_bf]);
      }
      free (pbi->alpha);
      free (pbi->beta);
    } // end of local model

    if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
  
      for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
        for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
          free (pbi->gamma[index_bf][index_l]);
          free (pbi->delta[index_bf][index_l]);
        }
        free (pbi->gamma[index_bf]);
        free (pbi->delta[index_bf]);
      }
      free (pbi->gamma);
      free (pbi->delta);
    } // end of equilateral and orthogonal models  
    
    free (pbi->delta_k);
    
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->cls[index_ct]);
    free (pbi->cls);
    if (pbi->has_intrinsic_squeezed == _TRUE_) {
      for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
        free (pbi->d_lsq_cls[index_ct]);
      free (pbi->d_lsq_cls);
    }
      
    
  } // end of if(has_bispectra)

  
  return _SUCCESS_;
 
}








/**
 * This routine defines indices and allocates tables in the bispectra structure 
 *
 */

int bispectra_indices (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi
        )
{ 

  // ================================================================================================
  // =                                        Count bispectra                                       =
  // ================================================================================================

  /* Find out which kind of bispectra to compute and assign them indices and labels */


  // ------------------------------------------------------------------
  // -                        Bispectra fields                        -
  // ------------------------------------------------------------------  

  /* Generate indices for the probes (T for temperature, E for E-mode polarisation,
  R for Rayleigh...) that we will use to build the bispectra and the Fisher matrix
  elements. */
  
  int index_bf = 0;
  
  pbi->has_bispectra_t = _FALSE_;
  pbi->has_bispectra_e = _FALSE_;
  pbi->has_bispectra_b = _FALSE_;
  pbi->has_bispectra_r = _FALSE_;
    
  if (ppt->has_bi_cmb_temperature == _TRUE_) {
    pbi->has_bispectra_t = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "t");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 0;
    pbi->index_bf_t = index_bf++;
  }
  
  if (ppt->has_bi_cmb_polarization == _TRUE_) {
    pbi->has_bispectra_e = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "e");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 2;
    pbi->index_bf_e = index_bf++;
  }

  if (ppt->has_bi_cmb_rayleigh == _TRUE_) {
    pbi->has_bispectra_r = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "r");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 0;
    pbi->index_bf_r = index_bf++;
  }

  pbi->bf_size = index_bf;
  pbi->n_probes = pow(pbi->bf_size, 3);

  class_test (pbi->bf_size > _MAX_NUM_FIELDS_,
    pbi->error_message,
    "we cannot compute the bispectrum for more than %d fields (e.g. T and E), reduce your expectations :-)");

  class_test (pbi->bf_size < 1,
    pbi->error_message,
    "no probes requested");

  /* Create labels for the full bispectra */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        sprintf (pbi->bfff_labels[X][Y][Z], "%s%s%s",
        pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);


  /* Associate to each field (T,E,...) its transfer function, which was computed in the transfer.c module,
  and to each possible pair of fields (TT,EE,TE,...) their power spectra, which were computed
  in the spectra.c module. */
  for (int X = 0; X < pbi->bf_size; ++X) {
    
    if ((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_t;
      pbi->index_ct_of_bf[X][X] = psp->index_ct_tt;
    }

    if ((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_e;
      pbi->index_ct_of_bf[X][X] = psp->index_ct_ee; 
    }

    if ((pbi->has_bispectra_r == _TRUE_) && (X == pbi->index_bf_r)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_r;
      pbi->index_ct_of_bf[X][X] = psp->index_ct_rr; 
    }

    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      if (((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_e == _TRUE_) && (Y == pbi->index_bf_e)))
        pbi->index_ct_of_bf[X][Y] = pbi->index_ct_of_bf[Y][X] = psp->index_ct_te;

      if (((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r)))
        pbi->index_ct_of_bf[X][Y] = pbi->index_ct_of_bf[Y][X] = psp->index_ct_tr;

      if (((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r)))
        class_test (1==1, pbi->error_message, "polarization-Rayleigh bispectrum not implemented yet");
        // pbi->index_ct_of_bf[X][Y] = pbi->index_ct_of_bf[Y][X] = psp->index_ct_er;

      if (((pbi->has_bispectra_b == _TRUE_) && (X == pbi->index_bf_b))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r)))
        class_test (1==1, pbi->error_message, "polarization-Rayleigh bispectrum not implemented yet");
        // pbi->index_ct_of_bf[X][Y] = pbi->index_ct_of_bf[Y][X] = psp->index_ct_br;
    }
  }
  

  // ------------------------------------------------------------------
  // -                        Bispectra types                         -
  // ------------------------------------------------------------------  
  
  int index_bt = 0;
  
  pbi->n[separable_bispectrum] = 0;
  pbi->n[non_separable_bispectrum] = 0;
  pbi->n[analytical_bispectrum] = 0;
  pbi->n[intrinsic_bispectrum] = 0;

  // *** Separable bispectra
  
  if (pbi->has_local_model) {
    pbi->index_bt_local = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_equilateral_model) {
    pbi->index_bt_equilateral = index_bt;
    strcpy (pbi->bt_labels[index_bt], "equilateral");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_orthogonal_model) {
    pbi->index_bt_orthogonal = index_bt;
    strcpy (pbi->bt_labels[index_bt], "orthogonal");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  // *** Non-separable bispectra

  if (pbi->has_galileon_model) {

    /* Bispectrum induced by pi_dot*pi_grad^2 */
    pbi->index_bt_galileon_gradient = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_grad");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    index_bt++;

    /* Bispectrum induced by pi_dot^3 */
    pbi->index_bt_galileon_time = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_time");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    index_bt++;
  }

  // *** Analytical bispectra

  if (pbi->has_local_squeezed == _TRUE_) {
    pbi->index_bt_local_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "l-squeezed");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_intrinsic_squeezed == _TRUE_) {
    pbi->index_bt_intrinsic_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "i-squeezed");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_cosine_shape == _TRUE_) {
    pbi->index_bt_cosine = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cosine");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_isw_lensing == _TRUE_) {
    pbi->index_bt_isw_lensing = index_bt;
    strcpy (pbi->bt_labels[index_bt], "isw-lensing");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }
  
  if (pbi->has_quadratic_correction == _TRUE_) {
    pbi->index_bt_quadratic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "quadratic");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  // *** Intrinsic (i.e. second-order) bispectra

  if (pbi->has_intrinsic == _TRUE_) {
    pbi->index_bt_intrinsic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic");
    pbi->bispectrum_type[index_bt] = intrinsic_bispectrum;
    pbi->n[intrinsic_bispectrum]++;
    index_bt++;
  }

  pbi->bt_size = index_bt;

  /* Perform some checks */
  class_test_permissive (
    (pbi->has_isw_lensing==_TRUE_) && ((pbi->bf_size>1) || ((pbi->bf_size==1 && pbi->has_bispectra_t==_FALSE_))),
    pbi->error_message,
    "the isw_lensing bispectrum is implemented only for temperature.");

  // =================================================================================================
  // =                                      Determine l-sampling                                     =
  // =================================================================================================
  
  /* We compute the bispectrum on a mesh where l1>=l2>=l3, with l3 determined by the triangular
  condition. All the multipoles are drawn from pbi->l, which is a copy of pbs->l.  */

  pbi->l_size = pbs->l_size;
  class_alloc (pbi->l, pbi->l_size*sizeof(int), pbi->error_message);
  for(int index_l=0; index_l<pbi->l_size; ++index_l)
    pbi->l[index_l] = pbs->l[index_l];

  /* Maximum value in pbi->l */
  pbi->l_max = pbi->l[pbi->l_size-1];
  pbi->full_l_size = pbi->l_max - 2 + 1;

  // *** Allocate & fill pb1->l_triangular_size and pb1->index_l_triangular_min
  
  /* The three multipole indexes l1, l2, l3 satisfy the triangular condition |l1-l2| <= l3 <= l1+l2.
  We choose to impose the condition on l3, for a fixed couple (l1,l2).  As a consequence, the
  allowed values for l3 will only be a subset of those contained in the vector pbs->l. The subset
  is determined by pbi->index_l_triangular_min[index_l1][index_l2] and pbi->l_triangular_size[index_l1][index_l2],
  which we fill below. */

  class_alloc (pbi->l_triangular_size, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l_triangular_min, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l_triangular_max, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l1_l2_l3, pbi->l_size*sizeof(long int **), pbi->error_message);

  /* Number of multipole configurations (l1,l2,l3) that satisfy the triangular condition and l1>=l2>=l3 */
  pbi->n_independent_configurations = 0;

  /* Number of multipole configurations (l1,l2,l3) that satisfy the triangular condition */
  pbi->n_total_configurations = 0;

  /* Initialise the index for a given (l1,l2,l3) triplet */
  long int index_l1_l2_l3 = 0;

  for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {

    int l1 = pbi->l[index_l1];

    /* Allocate l1 level */
    class_alloc (pbi->l_triangular_size[index_l1], pbi->l_size*sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_min[index_l1], pbi->l_size*sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_max[index_l1], pbi->l_size*sizeof(int), pbi->error_message);

    /* We consider only configurations whereby l1>=l2 */
    class_alloc (pbi->index_l1_l2_l3[index_l1], (index_l1+1)*sizeof(long int *), pbi->error_message);

    /* Fill pbi->l_triangular_size and pbi->index_l_triangular_min */
    for(int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
      
      int l2 = pbi->l[index_l2];

      /* Limits imposed on l3 by the triangular condition */
      int l_triangular_min = abs(l1-l2);
      int l_triangular_max = l1+l2;

      /* Find the index corresponding to l_triangular_min inside pbi->l */
      class_test (l_triangular_min >= pbi->l[pbs->l_size-1],
        pbi->error_message,
        "could not fulfill triangular condition for l3 using the multipoles in pbi->l");
      
      int index_l_triangular_min = 0;
      while (pbi->l[index_l_triangular_min] < l_triangular_min) ++index_l_triangular_min;


      /* Find the index corresponding to l_triangular_max inside pbi->l */
      class_test (l_triangular_max <= pbi->l[0],
        pbi->error_message,
        "could not fulfill triangular condition for l3 using the multipoles in pbi->l");
      
      int index_l_triangular_max = pbi->l_size-1;
      while (pbi->l[index_l_triangular_max] > l_triangular_max) --index_l_triangular_max;
    

      /* Fill pbi->index_l_triangular_min and pbi->l_triangular_size */
      pbi->index_l_triangular_min[index_l1][index_l2] = index_l_triangular_min;
      pbi->index_l_triangular_max[index_l1][index_l2] = index_l_triangular_max;
      pbi->l_triangular_size[index_l1][index_l2] = index_l_triangular_max - index_l_triangular_min + 1;

      /* Update counter of triangular configurations */
      pbi->n_total_configurations += pbi->l_triangular_size[index_l1][index_l2];

      /* We shall store the bispectra only for those configurations that contemporaneously satisfy 
      the triangular condition and the l1>=l2>=l3 condition. We use the pbi->index_l1_l2_l3 array
      to keep track of the index assigned to a given allowed configuration. */
      if (index_l2<=index_l1) {

        /* When the triangular condition is not compatible with index_l3<=index_l2, then index_l3_max < index_l3_min+1 will
        be either zero or negative. */ 
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
        int l3_size = MAX (0, index_l3_max-index_l3_min+1);
        class_alloc (pbi->index_l1_l2_l3[index_l1][index_l1-index_l2], l3_size*sizeof(long int), pbi->error_message);
        
        /* The indexing of pbi->index_l1_l2_l3 reflects the l1>=l2>=l3 constraint */
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
          pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3] = index_l1_l2_l3++;

      }
      
      /* Print some debug */
      // printf("l1=%d, l2=%d, l_triangular_min[%d]=%d, l_triangular_max[%d]=%d, l[index_l_triangular_min]=%d, l=[index_l_triangular_max]=%d, l_triangular_size=%d\n", 
      //   l1, l2, index_l_triangular_min, l_triangular_min, index_l_triangular_max, l_triangular_max, pbs->l[index_l_triangular_min], pbs->l[index_l_triangular_max], pbi->l_triangular_size[index_l1][index_l2]);
      
    } // end of for(index_l2)
  } // end of for(index_l1)

  /* Each bispectrum will be computed for the following number of configurations */
  pbi->n_independent_configurations = index_l1_l2_l3;
  

  /* Inform the user on how much her machine will have to suffer */
  if (pbi->bispectra_verbose > 1) {
    printf(" -> we shall compute %dx%d=%d bispectr%s for %ld configurations of (l1,l2,l3)\n",
      pbi->bt_size, pbi->n_probes, pbi->bt_size*pbi->n_probes,
      ((pbi->bt_size*pbi->n_probes)!=1?"a":"um"), pbi->n_independent_configurations);
    // printf("    with (L1,L2,L3) ranging from l=%d to %d (l_size=%d)\n", pbi->l[0], pbi->l[pbi->l_size-1], pbi->l_size);
  }
  
  
  
  
  
  // ============================================================================================
  // =                               Integration grid in k1 and k2                              =
  // ============================================================================================
  
  
  /* Sampling of the first-order transfer functions  */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];


  /* Allocate & fill delta_k, which is needed for the trapezoidal integration of the bispectrum.
  This array is has ptr->k_size elements, and is defined as k(i+1)-k(i-1) except for the first
  k(1)-k(0) and last k(N)-k(N-1) elements.  Note that when dealing with non-separable shapes,
  we shall not use this array for the integration over k3, as in that case the grid in k3 is
  not fixed but it depends on k1 and k2. */
  class_alloc (pbi->delta_k, k_tr_size * sizeof(double), pbi->error_message);
  
  /* Fill pbi->delta_k */
  pbi->delta_k[0] = k_tr[1] - k_tr[0];
      
  for (int index_k=1; index_k < k_tr_size-1; ++index_k)
    pbi->delta_k[index_k] = k_tr[index_k+1] - k_tr[index_k-1];
      
  pbi->delta_k[k_tr_size-1] = k_tr[k_tr_size-1] - k_tr[k_tr_size-2];
  
  
  
  
  
  
  // ==============================================================================================
  // =                             Allocate memory for bispectra                                  =
  // ==============================================================================================
  
  class_alloc (pbi->bispectra, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
  for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {

    class_alloc (pbi->bispectra[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

    for (int X = 0; X < pbi->bf_size; ++X) {
      
      class_alloc (pbi->bispectra[index_bt][X], pbi->bf_size*sizeof(double **), pbi->error_message);

      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        
        class_alloc (pbi->bispectra[index_bt][X][Y], pbi->bf_size*sizeof(double *), pbi->error_message);
        
        for (int Z = 0; Z < pbi->bf_size; ++Z)
          class_calloc (pbi->bispectra[index_bt][X][Y][Z], pbi->n_independent_configurations, sizeof(double), pbi->error_message);
      }
    }
  }
  
  pbi->count_allocated_for_bispectra = pbi->bt_size*pbi->n_probes*pbi->n_independent_configurations;
  
  if (pbi->bispectra_verbose > 2)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the bispectra array\n",
      pbi->count_allocated_for_bispectra*sizeof(double)/1e6, pbi->count_allocated_for_bispectra);




  
  
  // ==============================================================================================
  // =                            Create files for the bispectra                                  =
  // ==============================================================================================
  
  
  /* If it does not exist, create the directory for the bispectra.  Also create/open
  a status file.  The status file is an ASCII file used to determine the bispectra type that
  we already succesfully stored the bispectra. */
  if (pbi->store_bispectra_to_disk == _TRUE_) {

    sprintf(pbi->bispectra_run_directory, "%s/bispectra", ppr->run_directory);
    
    /* If we are not in a load run, just create the bispectra directory */
    if (ppr->load_run == _FALSE_) {
      class_test (mkdir (pbi->bispectra_run_directory, 0777)!=0,
                  pbi->error_message,
                  "could not create directory '%s', maybe it already exists?", pbi->bispectra_run_directory);
    }
    /* Otherwise, determine whether to skip the bispectrum module based on the content of the existence of the bispectra directory */
    else {

      struct stat st;
      short bispectra_dir_exists = (stat(pbi->bispectra_run_directory, &st)==0);

      /* If the bispectra folder already exists, assume that the bispectra have already been computed */
      if (bispectra_dir_exists) {

        if (pbi->bispectra_verbose > 1)
          printf (" -> found bispectra folder.\n");

        pbi->load_bispectra_from_disk = _TRUE_;
      }
      
      /* If the bispectra folder does not exist, add it to the load run and compute the bispectra */
      else {

        if (pbi->bispectra_verbose > 1)
          printf ("     * bispectra folder not found, will create it.\n");

        class_test (mkdir (pbi->bispectra_run_directory, 0777)!=0,
                    pbi->error_message,
                    "could not create directory '%s', maybe it already exists?", pbi->bispectra_run_directory);
        
        pbi->load_bispectra_from_disk = _FALSE_;

      } // end of if(bispectra_dir_exists)
      
    } // end of if(load_run)
    
  
    /* Create/open the status file. The 'a+' mode means that if the file does not exist it will be created,
    but if it exist it won't be erased (append mode) */
    sprintf(pbi->bispectra_status_path, "%s/bispectra_status_file.txt", ppr->run_directory);
    class_open(pbi->bispectra_status_file, pbi->bispectra_status_path, "a+", pbi->error_message);
    
  
    // *** Allocate bispectra files
    int index_bt;
  
    if (pbi->bispectra_verbose > 1)
      if (ppr->load_run == _FALSE_)
        printf ("     * will create %d files for the bispectra\n", pbi->bt_size);
    
    /* We are going to store the bispectra in n=k_size files, one for each requested k1 */
    class_alloc (pbi->bispectra_run_files, pbi->bt_size*sizeof(FILE *), pbi->error_message);
    class_alloc (pbi->bispectra_run_paths, pbi->bt_size*sizeof(char *), pbi->error_message);
  
    for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      
      /* Include the name of the bispectrum in its file */
      class_alloc (pbi->bispectra_run_paths[index_bt], _FILENAMESIZE_*sizeof(char), pbi->error_message);
      sprintf (pbi->bispectra_run_paths[index_bt], "%s/bispectra_%s.dat", pbi->bispectra_run_directory, pbi->bt_labels[index_bt]);
      
    } // end of for(index_bt)
  
  } // end of if(pbi->store_bispectra_to_disk)
  
  


  // ===========================================================================================
  // =                                     Interpolate P(k)                                    =
  // ===========================================================================================
  
  /* The primordial power spectrum for the comoving curvature perturbatio R was already
  computed when we initialized the ppm structure, but we store it locally for the k-values
  of ptr-k to access it faster */
  class_call (bispectra_primordial_power_spectrum(
                pba,
                ppt,
                ptr,
                ppm,
                pbi),
    pbi->error_message,
    pbi->error_message);





  // =============================================================================================
  // =                                   Interpolate the C_l's                                   =
  // =============================================================================================

  /* Interpolate the Cl's in all l-values */
  class_call (bispectra_cls(
                ppt,
                psp,
                ple,
                pbi),
    pbi->error_message,
    pbi->error_message);





  return _SUCCESS_;

}







/**
 * Compute and store in pbi->pk[index_k] the primordial power spectrum of the Newtonian
 * time potential psi. The power spectrum is computed in the points contained in ptr->k, as
 * these are the points where the transfer functions are computed.
 *
 * The purpose of this function is twofold. First, it stores the primordial power spectrum
 * into memory for faster access by the bispectra module, as calling primordial_spectrum_at_k
 * is fairly expensive. Secondly, we convert the dimensionless spectrum of the curvature
 * perturbation R outputted by the primordial module of CLASS, Delta_R(k), into the power
 * spectrum for the Newtonian curvature potential, P_Phi(k). The
 * two spectra are related by:
 *  
 *  P_Phi(k) = 2*Pi^2/k^3 * Delta_R(k)
 *
 * where
 * 
 * Delta_R(k) = A_s * (k/k_pivot)^(n_s-1)
 *
 */  
int bispectra_primordial_power_spectrum (
    struct background * pba,
    struct perturbs * ppt,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi
    )
{

  /* Allocate the pbi->pk vector so that it contains ptr->k_size values */
  int k_size = ptr->k_size[ppt->index_md_scalars];
  class_alloc (pbi->pk, k_size*sizeof(double), pbi->error_message);
  
  /* Fill pk with the values of the primordial power spectrum, as obtained in the ppm module */

  for (int index_k=0; index_k<k_size; ++index_k) {
    
    double k = ptr->k[ppt->index_md_scalars][index_k];

    class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k,
                  &(pbi->pk[index_k])),
      ppm->error_message,
      pbi->error_message);

    /* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    pbi->pk[index_k] = 2*_PI_*_PI_/(k*k*k) * pbi->pk[index_k];
    
  } // end of for(index_k)
  
  

  /* Do the same, but with ppt->k */
  int k_pt_size = ppt->k_size[ppt->index_md_scalars];
  class_alloc (pbi->pk_pt, k_pt_size*sizeof(double), pbi->error_message);
  
  /* Fill pk with the values of the primordial power spectrum, as obtained in the ppm module */

  for (int index_k_pt=0; index_k_pt<k_pt_size; ++index_k_pt) {
    
    double k_pt = ppt->k[ppt->index_md_scalars][index_k_pt];

    class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k_pt,
                  &(pbi->pk_pt[index_k_pt])),
      ppm->error_message,
      pbi->error_message);

    pbi->pk_pt[index_k_pt] = 2*_PI_*_PI_/(k_pt*k_pt*k_pt) * pbi->pk_pt[index_k_pt];
    
  } // end of for(index_k_pt)
  
  
  return _SUCCESS_;
  
}



int bispectra_cls (
    struct perturbs * ppt,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi
    )
{
  

  // ==============================================================================================
  // =                                          Allocate arrays                                   =
  // ==============================================================================================

  /* Allocate the array that will contain the C_l's for all types and all l's. */
  class_alloc (pbi->cls, psp->ct_size*sizeof(double*), pbi->error_message);
  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    class_alloc (pbi->cls[index_ct], pbi->full_l_size*sizeof(double), pbi->error_message);
  
  /* Do the same for the logarithmic derivative of the C_l's */
  if (pbi->has_intrinsic_squeezed == _TRUE_) {
    class_alloc (pbi->d_lsq_cls, psp->ct_size*sizeof(double*), pbi->error_message);
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      class_alloc (pbi->d_lsq_cls[index_ct], pbi->full_l_size*sizeof(double), pbi->error_message);
  }
  
  
  /* We shall now call the CLASS function 'spectra_cl_at_l'. This gives three outputs:
  the total Cl's for each probe (T, E, B...); the Cl's divided by probe and mode
  (scalar, vector, tensor); the Cl's divided by probe, mode, and initial condition
  (adiabatic, isocurvature...). We have to allocate these three arrays before being
  able to call 'spectra_cl_at_l'. We do so copying what is done in the function
  'output_total_cl_at_l' in the output module */
  double * cl;        /* cl_md_ic[index_ct] */
  double ** cl_md;    /* cl_md[index_mode][index_ct] */
  double ** cl_md_ic; /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */
  class_alloc(cl, psp->ct_size*sizeof(double), pbi->error_message);	
  class_alloc(cl_md_ic, psp->md_size*sizeof(double *), pbi->error_message);
  class_alloc(cl_md, psp->md_size*sizeof(double *), pbi->error_message);
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {
    if (psp->md_size > 1)
      class_alloc(cl_md[index_mode], psp->ct_size*sizeof(double), pbi->error_message);	
    if (psp->ic_size[index_mode] > 1)
      class_alloc(cl_md_ic[index_mode], psp->ic_ic_size[index_mode]*psp->ct_size*sizeof(double), pbi->error_message);
  }
  
  
  
  // ==========================================================================================================
  // =                                                Store C_l's                                             =
  // ==========================================================================================================
  
  for (int l=2; l<=pbi->l_max; ++l) {
    
    class_call(spectra_cl_at_l(
                 psp,
                 (double)l,
                 cl,
                 cl_md,
                 cl_md_ic),
      psp->error_message,
      pbi->error_message);
      
    /* Store the total Cl's into an array as a function of l and probe. By 'total' we mean the power
    spectrum summed over all the modes and initial conditions */
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      pbi->cls[index_ct][l-2] = cl[index_ct];

    /* To compute the squeezed limit approximation case, we need the derivative of 
    l*l*C_l. The best way to obtain these is to take the derivative of the
    spline-interpolated C_l's that we have just stored in pbi->cls. The reason is
    that the C_l's are smooth and the cubic interpolation does a very good job at
    filling the gaps in psp->cl, which contains only the C_l's computed by the spectra.c
    module. On the other hand, directly taking the derivative of the sparsely sampled C_l's
    in psp->cl will give a bad numerical convergence. We have tested that the first
    approach (derivative of spline-interpolated C_l's) converges much better than
    the second one (derivative of computed C_l's), with respect to the number of points
    taken in the l-grid. Uncomment the following lines to use the first approach. */

    // if (pbi->has_intrinsic_squeezed == _TRUE_) {
    // 
    //   class_call(spectra_dcl_at_l(
    //                psp,
    //                (double)l,
    //                cl,
    //                cl_md,
    //                cl_md_ic),
    //     psp->error_message,
    //     pbi->error_message);
    //   
    //   /* Store the total Cl's into an array as a function of l and probe. By 'total' we mean the power
    //   spectrum summed over all the modes and initial conditions */
    //   for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    //     pbi->d_lsq_cls[index_ct][l-2] = cl[index_ct];
    // 
    //   /* Some debug - print the cls */
    //   // printf ("%d %g %g %g\n", l,
    //   //   pbi->cls[psp->index_ct_tt][l-2],
    //   //   pbi->cls[psp->index_ct_tz][l-2],
    //   //   pbi->d_lsq_cls[psp->index_ct_tt][l-2]);
    // }
    
  } // end of for loop on the l's

  /* Compute the derivative of l*l*C_l from the spline-interpolated C_l's. This corresponds
  to the second method mentioned in the long comment above. */
  if (pbi->has_intrinsic_squeezed == _TRUE_) {
  
    double * l_array;
    class_alloc (l_array, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);
    for (int l=2; l<=pbi->l_max; ++l)
      l_array[l-2] = l;
  
    double * lsq_cl, * dd_lsq_cl;
    class_alloc (lsq_cl, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);
    class_alloc (dd_lsq_cl, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);
  
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct) {
  
      for (int l=2; l<=pbi->l_max; ++l)
        lsq_cl[l-2] = l*l*pbi->cls[index_ct][l-2];
    
      /* Compute the second derivatives of l^2*C_l */
      class_call(array_spline_table_lines(
                   l_array,
                   pbi->l_max-2+1,
                   lsq_cl,
                   1,
                   dd_lsq_cl,
                   _SPLINE_EST_DERIV_,
                   pbi->error_message),
        pbi->error_message,
        pbi->error_message);
        
      /* Compute the first derivative using the above information */
      class_call (array_spline_derive_table_lines(
                    l_array,
                    pbi->l_max-2+1,
                    lsq_cl,
                    dd_lsq_cl,
                    1,
                    pbi->d_lsq_cls[index_ct],
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
  
    } // end of loop on index_ct
    
    free (l_array);
    free(lsq_cl);
    free(dd_lsq_cl);
    
  } // end of if (has_intrinsic_squeezed)

  /* Free memory */
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {    
    if (psp->md_size > 1) free(cl_md[index_mode]);  
    if (psp->ic_size[index_mode] > 1) free(cl_md_ic[index_mode]);
  }  
  free(cl_md_ic);
  free(cl_md);
      
  return _SUCCESS_;
  
}



/**
 * This routine computes a table of values for all harmonic spectra C_l's,
 * given the transfer functions and primordial spectra.
 * 
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfers structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int bispectra_harmonic (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
    struct bispectra * pbi
    )
{

  /* Initialize counter for the values that will be memorised in pbi->bispectra */
  pbi->count_memorised_for_bispectra = 0;

  // =============================================================================================
  // =                             Compute separable bispectra                                   =
  // =============================================================================================

  /* We first obtain the bispectra whose primordial shape function is separable in (k1,k2,k3) because they
  are quick to compute. */
      
  if (pbi->n[separable_bispectrum] > 0) {

    struct bispectra_workspace_separable * pwb_sep;
    class_alloc (pwb_sep, sizeof(struct bispectra_workspace_separable), pbi->error_message);

    /* Compute the separable integrals, and integrate them together in r according to the chosen models */
    class_call (bispectra_separable_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  pwb_sep),
      pbi->error_message,
      pbi->error_message);
  

    /* Free the 'pwb_sep' workspace */
    class_call (bispectra_separable_workspace_free(
                  pbi,
                  pwb_sep),
      pbi->error_message,
      pbi->error_message);

    /* IMPORTANT: No symmetrization is needed, as the primary bispectrum is automatically symmetrised if we feed a
    primordial bispectrum B(k1,k2,k3) symmetrised in k, as we do. */

  }



  // ======================================================================================================
  // =                                    Compute analytic bispectra                                      =
  // ======================================================================================================
  
  
  /* Compute the bispectra obtained from simple analytical formulas, such as the lensing one */
  
  if (pbi->n[analytical_bispectrum] > 0) {
  
    class_call (bispectra_analytical_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  psp,
                  pbi),
      pbi->error_message,
      pbi->error_message);
      
  }
  
  /* Apart from the secondary bispectra in pbi->bispectra, all the arrays needed by the subsequent modules
  have been filled. If the user requested to load the bispectra from disk, we can stop the execution of
  this module now without regrets. */
  if (pbi->load_bispectra_from_disk == _TRUE_) {
    
    if (pbi->bispectra_verbose > 0)
      printf(" -> the intrinsic and non-separable bispectra will be read from disk\n");
    
    return _SUCCESS_;
  }
  
  
  // ===========================================================================================================
  // =                                   Compute non-separable bispectra                                       =
  // ===========================================================================================================
  
  if (pbi->n[non_separable_bispectrum] > 0) {
  
    struct bispectra_workspace_non_separable * pwb_nonsep;
    class_alloc (pwb_nonsep, sizeof(struct bispectra_workspace_non_separable), pbi->error_message);
  
    /* Compute the non-separable bispectra */
    class_call (bispectra_non_separable_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  pwb_nonsep),
      pbi->error_message,
      pbi->error_message);
  
  
    /* Free the 'pwb_nonsep' workspace */
    class_call (bispectra_non_separable_workspace_free(
                  pbi,
                  pwb_nonsep),
      pbi->error_message,
      pbi->error_message);
  }
  
  /* Print information on memory usage */
  if ((pbi->bispectra_verbose > 1) && (pbi->n[intrinsic_bispectrum]) < 1)
    printf(" -> memorised ~ %.3g MB (%ld doubles) in the bispectra array\n",
      pbi->count_memorised_for_bispectra*sizeof(double)/1e6, pbi->count_memorised_for_bispectra);
  
  
   
  // =======================================================================================================
  // =                                       Save bispectra to disk                                        =
  // =======================================================================================================
  
  if (pbi->store_bispectra_to_disk == _TRUE_) {

    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

      /* Save the non-separable bispectra to disk. We do not save the other ones because
      they take no time to recompute. */
      if (pbi->bispectrum_type[index_bt] == non_separable_bispectrum)
        class_call (bispectra_save_to_disk (
                      pbi,
                      index_bt),
          pbi->error_message,
          pbi->error_message);
    }  
  }
  
  /* Check that we correctly filled the bispectra array (but only if there are no
  intrinsic bispectra left to be computed)*/
  if (pbi->n[intrinsic_bispectrum] < 1)
    class_test_permissive (pbi->count_allocated_for_bispectra != pbi->count_memorised_for_bispectra,
      pbi->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!",
      pbi->count_allocated_for_bispectra, pbi->count_memorised_for_bispectra);
  
  // ============================================================================
  // =                      Check bispectra against nan's                       =
  // ============================================================================

  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* We have not computed the intrinsic bispectra yet, so we skip them */
    if (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {
    
          for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
            for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
     
              /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
              int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
              int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
     
              for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

                /* Index of the current (l1,l2,l3) configuration */
                long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
          
                double bispectrum = pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
   
                if (isnan(bispectrum))
                  printf ("@@@ WARNING: b(%d,%d,%d) = %g for bispectrum '%s_%s'.\n",
                  pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], bispectrum,
                  pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

              } // end of for(index_l3)
            } // end of for(index_l2)
          } // end of for(index_l1)
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)
  } // end of for(index_bt)
  
  // ============================================================================
  // =                      Check bispectra against nan's                       =
  // ============================================================================

  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* We have not computed the intrinsic bispectra yet, so we skip them */
    if (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {
    
          for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
            for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
     
              /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
              int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
              int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
     
              for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

                /* Index of the current (l1,l2,l3) configuration */
                long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
          
                double bispectrum = pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
   
                if (isnan(bispectrum))
                  printf ("@@@ WARNING: b(%d,%d,%d) = %g for bispectrum '%s_%s'.\n",
                  pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], bispectrum,
                  pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

              } // end of for(index_l3)
            } // end of for(index_l2)
          } // end of for(index_l1)
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)
  } // end of for(index_bt)
  
  return _SUCCESS_;

}








int bispectra_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{

  // ===================================================================================================
  // =                                    Prepare integration grid                                     =
  // ===================================================================================================
  
  /* We set the r-sampling as if it were a time sampling.  We do so because 'r' has the right dimensions, and
  it always appears in the argument of a Bessel function multiplying a wavemode, just as it was for conformal
  time in the line-of-sight integral.  This is the only place in the module where the background structure
  is accessed. */
  pwb->r_min = pbi->r_min;
  pwb->r_max = pbi->r_max;
  pwb->r_size = pbi->r_size;
    
  /* We decide to sample r linearly */
  class_alloc (pwb->r, pwb->r_size*sizeof(double), pbi->error_message);
  lin_space (pwb->r, pwb->r_min, pwb->r_max, pwb->r_size);
  
  
  /* Allocate & fill delta_r, the measure for the trapezoidal integration over r */
  class_alloc (pwb->delta_r, pwb->r_size * sizeof(double), pbi->error_message);

  /* Fill pwb->delta_r */
  pwb->delta_r[0] = pwb->r[1] - pwb->r[0];
      
  for (int index_r=1; index_r < pwb->r_size-1; ++index_r)
    pwb->delta_r[index_r] = pwb->r[index_r+1] - pwb->r[index_r-1];
      
  pwb->delta_r[pwb->r_size-1] = pwb->r[pwb->r_size-1] - pwb->r[pwb->r_size-2];






  // =====================================================================================================
  // =                                  Allocate memory for filter functions                             =
  // =====================================================================================================
  
  /* Cycle variables */
  int index_l;
  

  /* Allocate the separable integrals needed for each requested model of primordial non-Gaussianity.
  Refer to the header file for details on the models. The integrand arrays need to be written
  by different threads at the same time, hence we allocate one of them for each thread. */
      
  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel
  #ifdef _OPENMP
  number_of_threads = omp_get_num_threads();
  #endif

  // -------------------------------------------------------
  // -            Needed by all bispectra models           -
  // -------------------------------------------------------
  
  /* The alpha and beta functions need to be computed for all models */      
  if ((pbi->has_local_model == _TRUE_) ||
      (pbi->has_equilateral_model == _TRUE_) ||
      (pbi->has_orthogonal_model == _TRUE_)) {

    class_alloc (pbi->alpha, pbi->bf_size*sizeof(double **), pbi->error_message);
    class_alloc (pbi->beta, pbi->bf_size*sizeof(double **), pbi->error_message);

    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

      class_alloc (pbi->alpha[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
      class_alloc (pbi->beta[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
  
      /* Allocate r-level of the integral arrays */
      for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
        class_calloc (pbi->alpha[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
        class_calloc (pbi->beta[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
      }
    }
    
    /* Allocate memory for the integrand functions */
    class_alloc (pwb->alpha_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    class_alloc (pwb->beta_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    
    #pragma omp parallel private (thread)
    {
      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
      
      class_calloc_parallel (pwb->alpha_integrand[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);
      class_calloc_parallel (pwb->beta_integrand[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);      
    }
    
  } // end of if all models



  // ---------------------------------------------------------
  // -          Specific to equilateral & orthogonal         -
  // ---------------------------------------------------------

  if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {

    class_alloc (pbi->gamma, pbi->bf_size*sizeof(double **), pbi->error_message);
    class_alloc (pbi->delta, pbi->bf_size*sizeof(double **), pbi->error_message);

    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

      class_alloc (pbi->gamma[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
      class_alloc (pbi->delta[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
  
      /* Allocate r-level of the integral arrays */
      for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
        class_calloc (pbi->gamma[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
        class_calloc (pbi->delta[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
      }
    }
    
    /* Allocate memory for the integrand functions */
    class_alloc (pwb->gamma_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    class_alloc (pwb->delta_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    
    #pragma omp parallel private (thread)
    {
      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
      
      class_calloc_parallel (pwb->gamma_integrand[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);
      class_calloc_parallel (pwb->delta_integrand[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);      
    }
    
  } // end of if equilateral || orthogonal

  return _SUCCESS_;
  
}




int bispectra_separable_filter_functions (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bf,
    struct bispectra_workspace_separable * pwb
    )
{
  
  
  /* We shall integrate over the same k-points where we computed the first-order transfer functions */
  int k_size = ptr->k_size[ppt->index_md_scalars];

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel                                             \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)                         \
    private (thread)
  {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    /* Cycle variables, allocate inside parallel loop */
    int index_l, index_k, index_r;
      
    #pragma omp for schedule (dynamic)
    for (int index_l = 0; index_l < pbi->l_size; ++index_l) {

      if (pbi->bispectra_verbose > 2)
        printf("      \\ computing filter functions for l=%d, index_l=%d\n", pbi->l[index_l], index_l);

      /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
      int tt_size = ptr->tt_size[ppt->index_md_scalars];
      int l_size = ptr->l_size[ppt->index_md_scalars];

      double * transfer = &(ptr->transfer
        [ppt->index_md_scalars]
        [((ppt->index_ic_ad * tt_size + pbi->index_tt_of_bf[index_bf]) * l_size + index_l) * k_size]);


      // ==============================================
      // =          Build integrand functions         =
      // ==============================================

      /* Build the integrand functions, leaving out the k^2 * j_l(k*r) bit and any other numerical factor. */
    
      // *** Local and equilateral model ***
      if ((pbi->has_local_model == _TRUE_) ||
          (pbi->has_orthogonal_model == _TRUE_) ||
          (pbi->has_equilateral_model == _TRUE_)) {

        for (int index_k=0; index_k < ptr->k_size[ppt->index_md_scalars]; ++index_k) {
        
          double pk = pbi->pk[index_k];
          double tr = transfer[index_k];

          pwb->alpha_integrand[thread][index_k] = tr;
          pwb->beta_integrand[thread][index_k] = tr * pk;
          
          /* Print out the transfer function for a given field */
          // if (ptr->l[index_l]==284)
          //   if ((ppt->has_cl_cmb_polarization == _TRUE_) && (pbi->index_tt_of_bf[index_bf] == ptr->index_tt_e))
          //     fprintf (stderr, "%12g %12g\n", ptr->k[ppt->index_md_scalars][index_k], tr);
          
        }

      } // end of all models
    
      // *** Equilateral and orthogonal models ***
      if ((pbi->has_equilateral_model == _TRUE_) ||
          (pbi->has_orthogonal_model == _TRUE_)) {
     
        /* Here we basically copy eqs. 15-18 in Creminelli et al. 2006, keeping out the 2/pi factor and the
        Bessel function, and multiplying the rest by k (we use the dimensional power spectrum and we already
        factored out a k^2 factor) */
        for (int index_k=0; index_k < ptr->k_size[ppt->index_md_scalars]; ++index_k) {      

          double pk_one_third = pow(pbi->pk[index_k], 1/3.);
          double pk_two_thirds = pk_one_third*pk_one_third;
          double tr = transfer[index_k];

          pwb->gamma_integrand[thread][index_k] = tr * pk_one_third;
          pwb->delta_integrand[thread][index_k] = tr * pk_two_thirds;

        }
     
      } // end of equilateral and orthogonal model





      // =========================================================
      // =         Convolve with the Bessel function             =
      // =========================================================
      
      for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
        double r = pwb->r[index_r];
      
        if (pbi->bispectra_verbose > 3)
          printf("       \\ r=%g, index_r=%d\n", r, index_r);
        
        
        // *** All models ***
        if ((pbi->has_local_model == _TRUE_) ||
            (pbi->has_equilateral_model == _TRUE_) ||
            (pbi->has_orthogonal_model == _TRUE_)) {
                      
          /* alpha integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->k[ppt->index_md_scalars],
                        pbi->delta_k,
                        ptr->k_size[ppt->index_md_scalars],
                        pwb->alpha_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->alpha[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);

          /* beta integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->k[ppt->index_md_scalars],
                        pbi->delta_k,
                        ptr->k_size[ppt->index_md_scalars],
                        pwb->beta_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->beta[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);
                
        } // end of local model
        
        // /* Some debug */
        // if (index_r == 70)
        //   fprintf (stderr, "%10d %17.7g %17.7g\n",
        //   pbi->l[index_l], pbi->alpha[index_bf][index_l][index_r], pbi->beta[index_bf][index_l][index_r]);

        
        // *** Equilateral and orthogonal models ***
        if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {

          /* gamma integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->k[ppt->index_md_scalars],
                        pbi->delta_k,
                        ptr->k_size[ppt->index_md_scalars],
                        pwb->gamma_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->gamma[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);

          /* delta integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->k[ppt->index_md_scalars],
                        pbi->delta_k,
                        ptr->k_size[ppt->index_md_scalars],
                        pwb->delta_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->delta[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);
      
        } // end of equilateral and orthogonal models
        
        #pragma omp flush(abort)
      
      } // end of for(index_r)
  
    } // end of for(index_l)
    
  } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */

    
  /* Output the filter functions */
  // int index_l = 60;
  // 
  // if ((pbi->has_bispectra_t == _TRUE_) && (index_bf == pbi->index_bf_t) ) {
  // 
  //   /* Output beta = b_l */ 
  //   fprintf (stderr, "# Printing BETA filter function for l=%d\n", pbi->l[index_l]);
  //   fprintf (stderr, "%17s %17s\n", "r", "beta(l,r)");
  //   for (int index_r = 0; index_r < pwb->r_size; ++index_r)
  //     fprintf (stderr, "%17.7g %17.7g\n", pwb->r[index_r], pbi->beta[index_bf][index_l][index_r]);
  // 
  //   /* Output alpha = b_nl */  
  //   fprintf (stderr, "# Printing ALPHA filter function for l=%d\n", pbi->l[index_l]);
  //   fprintf (stderr, "%17s %17s\n", "r", "alpha(l,r)");
  //   for (int index_r = 0; index_r < pwb->r_size; ++index_r)
  //     fprintf (stderr, "%17.7g %17.7g\n", pwb->r[index_r], pbi->alpha[index_bf][index_l][index_r]);
  //   }
  // }

  return _SUCCESS_;
  
}









int bispectra_separable_integrate_over_r (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int X, // index for the first field (T,B,E)
    int Y, // index for the second field (T,B,E)
    int Z, // index for the third field (T,B,E)
    struct bispectra_workspace_separable * pwb
    )

{

  if (pbi->bispectra_verbose > 2)
    printf("     * computing the r-integral for the bispectrum %s_%s%s%s\n",
    pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  /* Shortcuts */
  double *** alpha = pbi->alpha;
  double *** beta = pbi->beta;
  double *** gamma = pbi->gamma;
  double *** delta = pbi->delta;
  
  // =======================================================================================
  // =                              Compute bispectrum r-integral                          =
  // =======================================================================================

  /* We now proceed to the final step of the bispectrum computation with the smooth_separable technique. This time
  we shall loop over all the l1-l2-l3 configurations, compute the simple integral over r, and store the bispectrum
  in pbi->bispectra.

    When computing the primordial bispectra, we shall assume a primordial f_NL of unity for phi. So far, we
  derived transfer functions for R. The curvature perturbation R and the potential during matter domination
  phi are related by R = -5/3 phi. The bispectrum integral has two power spectra (equivalent to 4 R's)
  and 3 transfer functions.  When converting to phi, each R brings a (-5/3) while each transfer function
  brings a (-3/5). We are then left with an overall factor of (-5/3), which means that a bispectrum
  with fnl_phi=1 is equivalent to a bispectrum with fnl_R=-3/5.
  
     We write down the r-integral using the same notation as in Yadav & Wandelt 2010 (eq. 32,
  35 and 38 of http://arxiv.org/abs/1006.0275v3), so that the i,j,k indices there correspond 
  to those in this function (T, E, ...). See also eq. 12 of Yadav & Wandelt 2007. */

  double fnl_R = -3/5.;

  /* We parallelize the outer loop over 'l1'. */
  #pragma omp parallel                                     \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)       \
    private (thread)
  {
  
    /* Cycle variables, declared inside parallel loop */ 
    int index_l1, index_l2, index_l3, index_r;
    
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      if (pbi->bispectra_verbose > 2)
        printf("      \\ computing r-integral for l1=%d, index_l1=%d\n", pbi->l[index_l1], index_l1);
    
      for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (pbi->l[index_l2] < pbi->l[index_l1]/2)
          continue;
  
        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

          /* Index of the current (l1,l2,l3) configuration */
          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  
          double integral = 0;
  
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
            /* Value of the considered r */
            double r = pwb->r[index_r];
            double integrand;
    
            // ------------------------------------------------------------------------
            // -                               Local model                            -
            // ------------------------------------------------------------------------

            if ((pbi->has_local_model == _TRUE_) && (index_bt == pbi->index_bt_local)) {
  
              /* The primordial bispectrum for the local model has only two extra contributions due to the symmetrization
              P(k1)*P(k2) + P(k1)*P(k3) + P(k2)*P(k3). It is easy to see that each contribution has the same value, but
              we sum them nonetheless to be consistent. */
              integrand = 2 * fnl_R * r*r * (
                  alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                );
            
            } // end of local model

            // -------------------------------------------------------------------------
            // -                           Equilateral model                           -
            // -------------------------------------------------------------------------
            
            else if ((pbi->has_equilateral_model == _TRUE_) && (index_bt == pbi->index_bt_equilateral)) {
  
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                - beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                - beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                    
                /* Symmetrical part */
                - 2 * delta[X][index_l1][index_r] * delta[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + beta[X][index_l1][index_r]  * gamma[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * delta[Y][index_l2][index_r] * gamma[Z][index_l3][index_r]
                + gamma[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * delta[Z][index_l3][index_r]
                + gamma[X][index_l1][index_r] * delta[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                + delta[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * gamma[Z][index_l3][index_r]
                + delta[X][index_l1][index_r] * gamma[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
              
              );

            } // end of equilateral model

            // -------------------------------------------------------------------------------
            // -                              Orthogonal model                               -
            // -------------------------------------------------------------------------------
            
            else if ((pbi->has_orthogonal_model == _TRUE_) && (index_bt == pbi->index_bt_orthogonal)) {
  
              /* We take the formula from Senatore et al. 2010, also shown in Komatsu et al. 2011 (WMAP7 paper, eq. 64) */
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - 3 * alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                - 3 * beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                - 3 * beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                    
                /* Symmetrical part. We found what we think is a typo in eq. 38 of http://arxiv.org/abs/1006.0275v3,
                where the coefficient is -2/3*3 = -2 instead of -8. We think -8 is the correct coefficient, as it can
                be verified from eq. 3.2 of Senatore, Smith & Zaldarriaga 2010.  */
                - 8 * delta[X][index_l1][index_r] * delta[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + 3 * beta[X][index_l1][index_r]  * gamma[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                + 3 * beta[X][index_l1][index_r]  * delta[Y][index_l2][index_r] * gamma[Z][index_l3][index_r]
                + 3 * gamma[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * delta[Z][index_l3][index_r]
                + 3 * gamma[X][index_l1][index_r] * delta[Y][index_l2][index_r] * beta[Z][index_l3][index_r]             
                + 3 * delta[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * gamma[Z][index_l3][index_r]
                + 3 * delta[X][index_l1][index_r] * gamma[Y][index_l2][index_r] * beta[Z][index_l3][index_r]             
              
              );

            } // end of orthogonal model

            
            /* Increment the estimate of the integral */
            integral += integrand * pwb->delta_r[index_r];

            /* Print integrand as a function of r */
            // if ((index_l1 == pbi->l_size-1) && (index_l2 == pbi->l_size-1) && (index_l3 == pbi->l_size-1))
            // if ((index_l1 == 0) && (index_l2 == 0) && (index_l3 == 0))
            // if ((index_l1 == 0) && (index_l2 == 0) && (index_l3 == 0))
            //   fprintf (stderr, "%15.7g %15.7g\n", r, integrand);

  
          } // end of for(index_r)
  

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3), including the factor 1/2
          from trapezoidal rule */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = 0.5 * integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. In KSW2005, this factor
          was split between the alpha and beta integrals, but from the numerical point of view it is preferable
          to include it at the end of the computation (as in eq. 17 of Fergusson & Shellard 2007). */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] *= pow(2./_PI_,3);
  
          /* Some debug */
          // if ((index_l1==index_l2) && (index_l2==index_l3))
          //   fprintf (stderr, "%10d %17.7g\n", pbi->l[index_l1], pbi->bispectra[index_bt][index_l1][index_l2][index_l3-index_l3_min]);
    
          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_for_bispectra++;

        } // end of for(index_l3)
      } // end of for(index_l2)
      
      #pragma omp flush(abort)
  
    } // end of for(index_l1)
  
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  return _SUCCESS_; 
  
}
  


/**
 * Compute the angular bispectrum for the separable primordial bispectra.
 *
*/
int bispectra_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{


  // ===================================================================================
  // =                                Initialise workspace                             =
  // ===================================================================================

  /* Allocate arrays inside the integration workspace */
  class_call (bispectra_separable_workspace_init(
                ppr,
                pba,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);


  if (pbi->bispectra_verbose > 1)
    printf (" -> computing separable bispectra; r sampled %d times in [%g,%g]\n", pwb->r_size, pwb->r_min, pwb->r_max);
    




  // ===================================================================================
  // =                             Compute filter functions                            =
  // ===================================================================================

  /* Compute the filter functions: alpha(l,r), beta(l,r), gamma(l,r), delta(l,r) for each field 
  (T,E,B). The filter functions are convolution of P(k)^A*T(k)^B with a Bessel function. */

  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

    if (pbi->bispectra_verbose > 1)
        printf("     * computing %s filter functions for the separable bispectra ...\n", pbi->bf_labels[index_bf]);

    class_call (bispectra_separable_filter_functions(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  index_bf,
                  pwb),
      pbi->error_message,
      pbi->error_message);

  }


  // ==========================================================================================
  // =                               Perform final integration                                =
  // ==========================================================================================

  /* Compute the bispectra by integrating the filter functions in r. This is where the primordial shape
  functions are used. */
  
  /* Loop on the type of bispectrum (local, equilateral, orthogonal...) */
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {
  
    if (pbi->bispectrum_type[index_bt] != separable_bispectrum)
      continue;

    /* Loop on the considered probe (TTT, TTE, TEE, EEE...) */
    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {
  
          class_call (bispectra_separable_integrate_over_r(
                        ppr,
                        pba,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_bt,
                        X,
                        Y,
                        Z,
                        pwb),
            pbi->error_message,
            pbi->error_message);
  
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)
    
  } // end of for(index_bt)
  
  
  return _SUCCESS_; 
  
}






int bispectra_separable_workspace_free (
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{

  free (pwb->r);
  free (pwb->delta_r);
 
  /* Arrays specific to the primordial models */

  int thread = 0;
  
  #pragma omp parallel private (thread)
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
      free (pwb->alpha_integrand[thread]);
      free (pwb->beta_integrand[thread]);
    }

    if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
      free (pwb->gamma_integrand[thread]);
      free (pwb->delta_integrand[thread]);
    }
  }

  if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
    free (pwb->alpha_integrand);
    free (pwb->beta_integrand);
  }

  if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
    free (pwb->gamma_integrand);
    free (pwb->delta_integrand);
  }
 
  return _SUCCESS_; 
  
}






/**
 * Compute the angular bispectrum for the analytical primordial bispectra.
 *
 */
int bispectra_analytical_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
    struct bispectra * pbi
    )
{  

  if (pbi->bispectra_verbose > 0)
    printf (" -> computing analytic bispectra\n");

  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel
  {
    #ifdef _OPENMP
    number_of_threads = omp_get_num_threads();
    #endif
  }

  // ===================================================================================
  // =                           Allocate three-j symbols                              =
  // ===================================================================================

  /* Allocate arrays needed to compute the three-j symbols. These are only needed to obtain 
  the quadratic correction bispectrum. For details on why we compute these three-j symbols,
  refer to the very long comment in the innermost loop, below. We will need 7 different
  combination of 3j symbols */
  enum {
    l1_l2_l3_F_X3_mF_X3,
    l1_l2_l3_0_mF_X3,
    l2_l3_l1_F_X1_mF_X1,
    l2_l3_l1_0_mF_X1,
    l3_l1_l2_F_X2_mF_X2,
    l3_l1_l2_0_mF_X2,
    l1_l2_l3_0_0,
  };  
  int n_geometrical_factors = 7;
  int size[number_of_threads][n_geometrical_factors];
  int min[number_of_threads][n_geometrical_factors];
  int max[number_of_threads][n_geometrical_factors];
  double *** value;

  if (pbi->has_quadratic_correction == _TRUE_) {  
    /* Temporary arrays and values needed to store the results of the 3j and 6j computations */
    class_alloc (value, number_of_threads*sizeof(double **), pbi->error_message);
    for (int thread=0; thread < number_of_threads; ++thread) {
      class_alloc (value[thread], n_geometrical_factors*sizeof(double *), pbi->error_message);
      for (int i=0; i < n_geometrical_factors; ++i)
        class_alloc (value[thread][i], (2*pbi->l_max+1)*sizeof(double), pbi->error_message);
    }
  }

  // ===================================================================================
  // =                                   Main loop                                     =
  // ===================================================================================

  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
    
    if (pbi->bispectrum_type[index_bt] != analytical_bispectrum)
      continue;

      for (int X = 0; X < pbi->bf_size; ++X) {
        for (int Y = 0; Y < pbi->bf_size; ++Y) {
          for (int Z = 0; Z < pbi->bf_size; ++Z) {
    
          if (pbi->bispectra_verbose > 2)
            printf(" -> computing the bispectrum (%s_%s)\n",
              pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);
    
          /* We parallelize the outer loop over 'l1'. */
          #pragma omp parallel for shared (ppt,pbs,ptr,ppm,pbi,abort) private (thread)
          for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
      
            #ifdef _OPENMP
            thread = omp_get_thread_num();
            #endif
      
            int l1 = pbi->l[index_l1];
            double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
        
            for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
    
              /* Skip those configurations that are forbidden by the triangular condition (optimization) */
              if (pbi->l[index_l2] < pbi->l[index_l1]/2)
                continue;
    
              int l2 = pbi->l[index_l2];
              double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
    
              /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
              int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
              int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);

              for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  
    
                int l3 = pbi->l[index_l3];
                double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
                
                /* Parity of the considered configuration */
                int L = l3-l1-l2;
    
                /* Index of the current (l1,l2,l3) configuration */
                long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
    
                /* Initialize bispectrum */
                pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = 0;
        
                // ---------------------------------------------------------------------------------------------
                // -                          Squeezed limit of the intrinsic bispectrum                       -
                // ---------------------------------------------------------------------------------------------
        
                /* Here we compute the approximations in eq. 4.1 and 4.2 of Lewis 2012. This is the general
                formula that includes polarisation. With respect to Lewis' formula, (i,l1)->(Z,l3), (j,l2)->(X,l1),
                (k,l3)->(Y,l2) and \zeta -> z. It is crucial that l3 is associated with the C_l that correlates
                with the comoving curvature perturbation (zeta), because l3 in this loop is the smallest multipole,
                which describes the long wavelength mode. Also, l3 must be associated with Z because our convention
                for the bispectrum is <X_l1 Y_l2 Z_l3>. If you associate l3 with another field, then the Fisher matrix
                estimator will associate to that field the wrong covariance matrix, and the result will change
                drastically. */
                if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed)) {
    
                  /* Determine which field has to be correlated with the comoving curvature perturbation 'z'.
                  The resulting C_l^Zz always gets the largest scale multipole */
                  int index_ct_Zz;
                      
                  if ((pbi->has_bispectra_t == _TRUE_) && (Z == pbi->index_bf_t))
                    index_ct_Zz = psp->index_ct_tz;
                  else if ((pbi->has_bispectra_e == _TRUE_) && (Z == pbi->index_bf_e))
                    index_ct_Zz = psp->index_ct_ez;
                  else {
                    index_ct_Zz = 0;
                    printf ("WARNING: Could not find C_l's for <Z * zeta> with Z=%s not found. Do not trust squeezed approximation.\n",
                    pbi->bf_labels[Z]);
                  }
                      
                  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
                  choice as the l-loop we are into is constructed to have l1>=l2>=l3  */
                  double cl3_Zz = pbi->cls[index_ct_Zz][l3-2];
                  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf[X][Y]][l1-2];
                  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf[X][Y]][l2-2];
                            
                  /* Ricci focussing in Lewis 2012 (eq. 4.1). It should be noted that this expression
                  is NOT symmetric with respect to a simultaneous exchange of (l1,X) <-> (l2,Y),
                  (l1,X) <-> (l3,Z), contrary to the other bispectra, which are computed as <X_l1 Y_l2 Z_l3>.
                  Therefore, for this bispectrum one cannot obtain the configurations outside l1>=l2>=l3 by
                  permuting the XYZ indices, as it is done, for example, in the print_bispectra function. */
                  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);
                  
                  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
                  double bolometric_T_lewis_redshift = 0;
                            
                  if (pbi->has_bispectra_t == _TRUE_) {
                    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
                    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf[Z][pbi->index_bf_t]][l3-2];
                    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf[X][pbi->index_bf_t]][l1-2];
                    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf[Y][pbi->index_bf_t]][l2-2];
                    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
                  }
                            
                  class_test_parallel ((pbi->has_bispectra_e == _TRUE_) &&
                    ((X == pbi->index_bf_e) && (Y == pbi->index_bf_e))
                    && (bolometric_T_lewis_redshift!=0),
                    pbi->error_message,
                    "anisotropic redshifting should vanish (eq. 4.2 of Lewis 2012)", bolometric_T_lewis_redshift);
                            
                  /* Sum of Ricci focussing and redshift modulation */
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

                  /* At some point you might want to symmetrise the analytic approximation with respect to the
                  exchanges (l1,X) <-> (l2,Y) and (l1,X) <-> (l3,Z). Not sure it would make
                  a lot of sense, because one would include in a given (l1,l2,l3) configurations contributions
                  where the analytical approximation is not valid, i.e. one would have the short mode in the
                  C_l with the zeta correlation and the long mode in the C_l derivative.  */
                  /* Determine which field has to be correlated with the comoving curvature perturbation zeta.
                  The resulting C_l^Xzeta always gets the largest scale multipole */
                  // int index_ct_Xz;
                  // int index_ct_Yz;
                  // int index_ct_Zz;
                  //     
                  // if (pbi->has_bispectra_t == _TRUE_) {
                  //   if (X == pbi->index_bf_t) index_ct_Xz = psp->index_ct_tz;
                  //   if (Y == pbi->index_bf_t) index_ct_Yz = psp->index_ct_tz;
                  //   if (Z == pbi->index_bf_t) index_ct_Zz = psp->index_ct_tz;
                  // }
                  // else if (pbi->has_bispectra_e == _TRUE_) {
                  //   if (X == pbi->index_bf_e) index_ct_Xz = psp->index_ct_ez;
                  //   if (Y == pbi->index_bf_e) index_ct_Yz = psp->index_ct_ez;
                  //   if (Z == pbi->index_bf_e) index_ct_Zz = psp->index_ct_ez;
                  // }
                  // else {
                  //   index_ct_Xz = 0;
                  //   printf ("WARNING: C_l's for <X * zeta> where X=%s not found. Do not trust squeezed approximation.\n",
                  //   pbi->bf_labels[X]);
                  // }
                  // 
                  // double cl1_Xz = pbi->cls[index_ct_Xz][l1-2];
                  // double cl2_Yz = pbi->cls[index_ct_Yz][l2-2];
                  // double cl3_Zz = pbi->cls[index_ct_Zz][l3-2];
                  // 
                  // double dcl1_XZ = pbi->d_lsq_cls[pbi->index_ct_of_bf[X][Z]][l1-2];
                  // double dcl1_YX = pbi->d_lsq_cls[pbi->index_ct_of_bf[Y][X]][l1-2];
                  // double dcl2_YZ = pbi->d_lsq_cls[pbi->index_ct_of_bf[Y][Z]][l2-2];
                  // double dcl2_YX = pbi->d_lsq_cls[pbi->index_ct_of_bf[Y][X]][l2-2];
                  // double dcl3_YZ = pbi->d_lsq_cls[pbi->index_ct_of_bf[Y][Z]][l3-2];
                  // double dcl3_XZ = pbi->d_lsq_cls[pbi->index_ct_of_bf[X][Z]][l3-2];
                  // 
                  // /* The following formula is redundant (it can be reduced to three terms), but in this
                  // way it is clearer */
                  // double bolometric_T_lewis_ricci =
                  //   (- 0.5 * cl1_Xz * (dcl2_YZ/l2 + dcl3_YZ/l3)      /* <X_l1 Y_l2 Z_l3> */
                  //    - 0.5 * cl1_Xz * (dcl3_YZ/l3 + dcl2_YZ/l2)      /* <X_l1 Z_l3 Y_l2> */
                  //                                                    
                  //    - 0.5 * cl2_Yz * (dcl1_XZ/l1 + dcl3_XZ/l3)      /* <Y_l2 X_l1 Z_l3> */
                  //    - 0.5 * cl2_Yz * (dcl3_XZ/l3 + dcl1_XZ/l1)      /* <Y_l2 X_l1 Z_l3> */
                  //                                                    
                  //    - 0.5 * cl3_Zz * (dcl1_YX/l1 + dcl2_YX/l2)      /* <Z_l3 X_l1 Y_l2> */
                  //    - 0.5 * cl3_Zz * (dcl2_YX/l2 + dcl1_YX/l1))/6;  /* <Z_l3 Y_l2 X_l1> */
                  // 
                  // pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = bolometric_T_lewis_ricci;

                }
    
                // ---------------------------------------------------------------------------------------------
                // -                          Squeezed limit of the local bispectrum                           -
                // ---------------------------------------------------------------------------------------------
    
                if ((pbi->has_local_squeezed == _TRUE_) && (index_bt == pbi->index_bt_local_squeezed)) {
            
                  /* TODO: This is easily coded. Just copy bits from print_bispectra and code
                  - 5 / (12 * cl_Xz_long * cl_YZ_short) */
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = 6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
    
                }
    
                // ---------------------------------------------------------------------------------------------
                // -                          Squeezed limit of the local bispectrum                           -
                // ---------------------------------------------------------------------------------------------
    
                if ((pbi->has_cosine_shape == _TRUE_) && (index_bt == pbi->index_bt_cosine)) {
            
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] =
                    6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) * cos((l1+l2+l3)/50.);
    
                }
    
                // ----------------------------------------------------------------------------------------------
                // -                                   Lensing-ISW bispectrum                                   -
                // ----------------------------------------------------------------------------------------------
          
                /* TODO: This is not the correct formula, it is just a test. For the correct formula, refer
                to eq. 4.5 of Lewis, Challinor & Hanson 2011. */
                if ((pbi->has_isw_lensing == _TRUE_) && (index_bt == pbi->index_bt_isw_lensing)) {
            
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] +=
                    + (l1*(l1+1) + l2*(l2+1) - l3*(l3+1))
                    + l1*(l1+1) - l2*(l2+1) + l3*(l3+1)
                    - l1*(l1+1) + l2*(l2+1) + l3*(l3+1);
              
                }
            
                /* Update the counter */
                #pragma omp atomic
                pbi->count_memorised_for_bispectra++;
    
              } // end of for(index_l3)
              
              // ---------------------------------------------------------------------------------------------
              // -                                   Quadratic correction                                    -
              // ---------------------------------------------------------------------------------------------
            
              /* The CMB temperature is obtained from the second-order distribution function by adding to it
              a quadratic term (see Eq. 3.12 of arXiv:1401.3296).
              At the bispectrum level this extra term translates to a four-point function in the first-order
              perturbations, which in turn is expressed in terms of products of C_l's.
              The same correction term appears for the delta_tilde transformation (see Eq. 3.5 of arXiv:1401.3296).
              The correction in both cases has this form (see Eq. 3.6, 3.7 and 3.9 of same paper):
              
                QC_{l1 l2 l3} = S_X3 * i^(L+V_X3) * ( l'  l''  |    l3 ) * (  l'  l'' | l3 )  
                                                    ( 0   F_X3 | -F_X3 )   (  m'  m'' | m3  ) 
                                * <a^X1_l1m1 * a^X2_l2m2 * a^I_l'm' * a^T_X3_l''m''> + 2 permutations (1->2->3)
                
              where L=l3-l'-l'', all the a_lm's are first-order and, in this loop, X1=X, X2=Y and X3=Z. The
              permutations go over 1->2->3 and refer to the position of the second-order perturbation in
              the bispectrum; in the formula above, the positioning is <a^(1)X1_l1m1 * a^(1)X2_l2m2 * a^(2)X3_l3m3>.
              The coefficients take different values according to which field is X3:

                X3=I -> F=0, S=1, V_X3=0,  T_X3=I
                X3=E -> F=2, S=2, V_X3=0,  T_X3=E
                X3=B -> F=2, S=2, V_X3=-1, T_X3=E
              
              For X3=E-modes, the sum over l' and l'' only includes EVEN values of l3-l'-l'', for X3=B-modes
              only includes ODD values of l3-l'-l''.

              By employing Wick's theorem, the above can be expressed in terms of the angular power spectrum
              of the CMB, C_l = <a_lm a^*_lm> (note that the a_lm's already include the 1/4 factor):
              
                QC_{l1 l2 l3} = 4 * B_{l1 l2 l3} * i^V_X3 * G^m1m2m3_l1l2l3 * S_X3 * [ 
                                                   ( l1    l2     l3 ) * C^{X1,I}_l1 * C^{X2,X3}_l2
                                                   ( 0   F_X3  -F_X3 )
                                                 + (   l1  l2     l3 ) * C^{X1,X3}_l1 * C^{X2,I}_l2
                                                   ( F_X3   0  -F_X3 )
                                                                          ] + 2 permutations (1->2->3) ,
              
              where
              
                B_{l1 l2 l3} = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*_PI_)) * ( l1    l2   l3 )
                                                                           ( m1    m2   m3 ) .
              
              Below we compute the quadratic correction term using this general formula and store it in
              pbi->bispectra[index_bt_quadratic].
              
              For intensity (X1=X2=X3=I), the formula reduces to 8 * [ C_l1*C_l2 + C_l1*C_l2 + C_l1*C_l2 ], while             
              any contribution form the second-order B-modes has the V_X=-1 and is therefore purely imaginary.
              Furthermore, the B-mode contribution only includes one term, theone where a^B_lm is second order,
              as we ignore the first-order B-modes. In any case, we do not compute the quadratic term for a B-mode
              spectrum.
              
              */

              if ((pbi->has_quadratic_correction == _TRUE_) && (index_bt == pbi->index_bt_quadratic)) {
              
                // ---------------------------------------------------------------
                // -                  Determine field coefficients               -
                // ---------------------------------------------------------------
                
                /* Relabel X,Y,Z to match the notation in the long comment above */
                int X1=X, X2=Y, X3=Z;
                
                /* The spin of the fields is either 0 or 2, and we have stored in in the 'indices' function */
                int F_X1 = pbi->field_spin[X1];
                int F_X2 = pbi->field_spin[X2];
                int F_X3 = pbi->field_spin[X3];
                
                /* Determine the amplitude S_X of the term in QC that involves the second-order part of X,
                Also determine the field that goes in the C_l's, T_X, which is computed as T_I=I, T_E=E, T_B=E.  */
                double S_X1, S_X2, S_X3;
                int T_X1, T_X2, T_X3;
                
                /* Index of the cross power spectrum between X1, X2, X3 and I. It has to be set by hand,
                rather than using the array pbi->index_ct_of_bf because pbi->index_ct_of_bf only contains
                information on the fields appearing in one of the requested bispectrum. For example,
                if you only request EEE, then pbi->index_ct_of_bf does not contain information about <ET> */
                int index_ct_X1_I, index_ct_X2_I, index_ct_X3_I;
                
                if (pbi->has_bispectra_t == _TRUE_) {
                  if (X1 == pbi->index_bf_t) {S_X1 = 1; T_X1 = pbi->index_bf_t; index_ct_X1_I = psp->index_ct_tt;}
                  if (X2 == pbi->index_bf_t) {S_X2 = 1; T_X2 = pbi->index_bf_t; index_ct_X2_I = psp->index_ct_tt;}
                  if (X3 == pbi->index_bf_t) {S_X3 = 1; T_X3 = pbi->index_bf_t; index_ct_X3_I = psp->index_ct_tt;}
                }
                if (pbi->has_bispectra_r == _TRUE_) {
                  if (X1 == pbi->index_bf_r) {S_X1 = 1; T_X1 = pbi->index_bf_r; index_ct_X1_I = psp->index_ct_tr;}
                  if (X2 == pbi->index_bf_r) {S_X2 = 1; T_X2 = pbi->index_bf_r; index_ct_X2_I = psp->index_ct_tr;}
                  if (X3 == pbi->index_bf_r) {S_X3 = 1; T_X3 = pbi->index_bf_r; index_ct_X3_I = psp->index_ct_tr;}
                }
                if (pbi->has_bispectra_e == _TRUE_) {
                  if (X1 == pbi->index_bf_e) {S_X1 = 2; T_X1 = pbi->index_bf_e; index_ct_X1_I = psp->index_ct_te;}
                  if (X2 == pbi->index_bf_e) {S_X2 = 2; T_X2 = pbi->index_bf_e; index_ct_X2_I = psp->index_ct_te;}
                  if (X3 == pbi->index_bf_e) {S_X3 = 2; T_X3 = pbi->index_bf_e; index_ct_X3_I = psp->index_ct_te;}
                }
                if (pbi->has_bispectra_b == _TRUE_) { /* Note that <TB> vanishes, hence the negative values */
                  if (X1 == pbi->index_bf_b) {S_X1 = 2; T_X1 = pbi->index_bf_e; index_ct_X1_I = -1;}
                  if (X2 == pbi->index_bf_b) {S_X2 = 2; T_X2 = pbi->index_bf_e; index_ct_X2_I = -1;}
                  if (X3 == pbi->index_bf_b) {S_X3 = 2; T_X3 = pbi->index_bf_e; index_ct_X3_I = -1;}
                }
              
                // -----------------------------------------------------------------------------------
                // -                             Compute three-j symbols                             -
                // -----------------------------------------------------------------------------------
                
                /* Compute the three-j symbols. We want to do it outside the loop on l3 as our
                function for the 3j's computes them for all allowed values of l3 */
              
                double min_D, max_D;
              
                //                        l1,      l2,    l3,
                //                        /*0,*/ F_X3, -F_X3,
                class_call_parallel (drc3jj (
                                       l1, l2, 0, F_X3,
                                       &min_D, &max_D,
                                       value[thread][l1_l2_l3_F_X3_mF_X3],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l1_l2_l3_F_X3_mF_X3] = (int)(min_D + _EPS_);
                max[thread][l1_l2_l3_F_X3_mF_X3] = (int)(max_D + _EPS_);
                size[thread][l1_l2_l3_F_X3_mF_X3] = max[thread][l1_l2_l3_F_X3_mF_X3] - min[thread][l1_l2_l3_F_X3_mF_X3] + 1;
                
                //                        l1,      l2,    l3,
                //                        /*F_X3*/, 0, -F_X3,
                class_call_parallel (drc3jj (
                                       l1, l2, F_X3, 0,
                                       &min_D, &max_D,
                                       value[thread][l1_l2_l3_0_mF_X3],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l1_l2_l3_0_mF_X3] = (int)(min_D + _EPS_);
                max[thread][l1_l2_l3_0_mF_X3] = (int)(max_D + _EPS_);
                size[thread][l1_l2_l3_0_mF_X3] = max[thread][l1_l2_l3_0_mF_X3] - min[thread][l1_l2_l3_0_mF_X3] + 1;
                
                //                        l2,      l3,    l1,
                //                        /*0*/, F_X1, -F_X1,
                class_call_parallel (drc3jj (
                                       l1, l2, -F_X1, 0,
                                       &min_D, &max_D,
                                       value[thread][l2_l3_l1_F_X1_mF_X1],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l2_l3_l1_F_X1_mF_X1] = (int)(min_D + _EPS_);
                max[thread][l2_l3_l1_F_X1_mF_X1] = (int)(max_D + _EPS_);
                size[thread][l2_l3_l1_F_X1_mF_X1] = max[thread][l2_l3_l1_F_X1_mF_X1] - min[thread][l2_l3_l1_F_X1_mF_X1] + 1;
              
                //                        l2,      l3,    l1,
                //                        /*F_X1*/, 0, -F_X1,
                class_call_parallel (drc3jj (
                                       l1, l2, -F_X1, F_X1,
                                       &min_D, &max_D,
                                       value[thread][l2_l3_l1_0_mF_X1],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l2_l3_l1_0_mF_X1] = (int)(min_D + _EPS_);
                max[thread][l2_l3_l1_0_mF_X1] = (int)(max_D + _EPS_);
                size[thread][l2_l3_l1_0_mF_X1] = max[thread][l2_l3_l1_0_mF_X1] - min[thread][l2_l3_l1_0_mF_X1] + 1;  
              
                //                        l3,      l1,    l2,
                //                        /*0*/, F_X2, -F_X2,
                class_call_parallel (drc3jj (
                                       l1, l2, F_X2, -F_X2,
                                       &min_D, &max_D,
                                       value[thread][l3_l1_l2_F_X2_mF_X2],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l3_l1_l2_F_X2_mF_X2] = (int)(min_D + _EPS_);
                max[thread][l3_l1_l2_F_X2_mF_X2] = (int)(max_D + _EPS_);
                size[thread][l3_l1_l2_F_X2_mF_X2] = max[thread][l3_l1_l2_F_X2_mF_X2] - min[thread][l3_l1_l2_F_X2_mF_X2] + 1;
                
                //                        l3,      l1,    l2,
                //                        /*F_X2*/, 0, -F_X2,
                class_call_parallel (drc3jj (
                                       l1, l2, 0, -F_X2,
                                       &min_D, &max_D,
                                       value[thread][l3_l1_l2_0_mF_X2],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l3_l1_l2_0_mF_X2] = (int)(min_D + _EPS_);
                max[thread][l3_l1_l2_0_mF_X2] = (int)(max_D + _EPS_);
                size[thread][l3_l1_l2_0_mF_X2] = max[thread][l3_l1_l2_0_mF_X2] - min[thread][l3_l1_l2_0_mF_X2] + 1;
              
                //                        l1,   l2, l3,
                //                        /*0*/, 0,  0,
                class_call_parallel (drc3jj (
                                       l1, l2, 0, 0,
                                       &min_D, &max_D,
                                       value[thread][l1_l2_l3_0_0],
                                       (2*pbi->l_max+1),
                                       pbi->error_message       
                                       ),
                  pbi->error_message,
                  pbi->error_message);
              
                min[thread][l1_l2_l3_0_0] = (int)(min_D + _EPS_);
                max[thread][l1_l2_l3_0_0] = (int)(max_D + _EPS_);
                size[thread][l1_l2_l3_0_0] = max[thread][l1_l2_l3_0_0] - min[thread][l1_l2_l3_0_0] + 1;
              
                // -----------------------------------------------------------------------------------
                // -                                    Loop over l3                                 -
                // -----------------------------------------------------------------------------------
                
                for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  
                  
                  int l3 = pbi->l[index_l3];
                  double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
                
                  /* Parity of the considered configuration */
                  int L = l3-l1-l2;
                  
                  /* Index of the current (l1,l2,l3) configuration */
                  long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
                
                  /* If only TTT is needed, compute the quadratic (reduced) bispectrum in the simple way */
                  if ((pbi->has_bispectra_t == _TRUE_) && (pbi->bf_size == 1)) {
                    pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
                    goto end_of_quadratic_bispectrum; 
                  }
              
                  /* Take into account that the contribution from E^(2) only exists when l3-l1-l2 is even,
                  while that from B^(2) only when l3-l1-l2 is odd. */
                  if ((pbi->has_bispectra_e == _TRUE_) && (L%2!=0)) {
                    if (X1 == pbi->index_bf_e) S_X1 = 0;
                    if (X2 == pbi->index_bf_e) S_X2 = 0;
                    if (X3 == pbi->index_bf_e) S_X3 = 0;
                  }
                  if ((pbi->has_bispectra_b == _TRUE_) && (L%2==0)) {
                    if (X1 == pbi->index_bf_b) S_X1 = 0;
                    if (X2 == pbi->index_bf_b) S_X2 = 0;
                    if (X3 == pbi->index_bf_b) S_X3 = 0;
                  }
                          
                  /* B-mode bispectrum not implemented yet */
                  class_test_parallel (pbi->has_bispectra_b == _TRUE_,
                    pbi->error_message,
                    "quadratic correction for B-mode bispectrum not implemented yet.");
                  
                  /* Get the C_l's. When implementing the B-modes, remember to 
                  set the C_l's to zero, so that the only contribution comes from the
                  bispectrum with the second-order a^B_lm.  */
                  int I = pbi->index_bf_t;
                  double C_l1_X1_I   = pbi->cls[index_ct_X1_I][l1-2];
                  double C_l1_X1_TX2 = pbi->cls[pbi->index_ct_of_bf[ X1 ][ T_X2 ]][l1-2];
                  double C_l1_X1_TX3 = pbi->cls[pbi->index_ct_of_bf[ X1 ][ T_X3 ]][l1-2];
              
                  double C_l2_X2_I   = pbi->cls[index_ct_X2_I][l2-2];
                  double C_l2_X2_TX1 = pbi->cls[pbi->index_ct_of_bf[ X2 ][ T_X1 ]][l2-2];
                  double C_l2_X2_TX3 = pbi->cls[pbi->index_ct_of_bf[ X2 ][ T_X3 ]][l2-2];
                  
                  double C_l3_X3_I   = pbi->cls[index_ct_X3_I][l3-2];                  
                  double C_l3_X3_TX1 = pbi->cls[pbi->index_ct_of_bf[ X3 ][ T_X1 ]][l3-2];
                  double C_l3_X3_TX2 = pbi->cls[pbi->index_ct_of_bf[ X3 ][ T_X2 ]][l3-2];
                    
                  /* The sum includes three terms, corresponding to the three possible types
                  of second-order perturbations. The T_X1_term, for example, refers to the
                  term where a^X1_lm is second order. */
              
                  double threej_F_X3_mF_X3 = value[thread][l1_l2_l3_F_X3_mF_X3][l3 - min[thread][l1_l2_l3_F_X3_mF_X3]];
                  double threej_0_mF_X3 = value[thread][l1_l2_l3_0_mF_X3][l3 - min[thread][l1_l2_l3_0_mF_X3]];
                  double T_X3_term = 4 * S_X3 * (threej_F_X3_mF_X3*C_l1_X1_I*C_l2_X2_TX3 + threej_0_mF_X3*C_l1_X1_TX3*C_l2_X2_I);
              
                  double threej_F_X1_mF_X1 = value[thread][l2_l3_l1_F_X1_mF_X1][l3 - min[thread][l2_l3_l1_F_X1_mF_X1]];
                  double threej_0_mF_X1 = value[thread][l2_l3_l1_0_mF_X1][l3 - min[thread][l2_l3_l1_0_mF_X1]];
                  double T_X1_term = 4 * S_X1 * (threej_F_X1_mF_X1*C_l2_X2_I*C_l3_X3_TX1 + threej_0_mF_X1*C_l2_X2_TX1*C_l3_X3_I);
              
                  double three_j_F_X2_mF_X2 = value[thread][l3_l1_l2_F_X2_mF_X2][l3 - min[thread][l3_l1_l2_F_X2_mF_X2]];
                  double three_j_0_mF_X2 = value[thread][l3_l1_l2_0_mF_X2][l3 - min[thread][l3_l1_l2_0_mF_X2]];
                  double T_X2_term = 4 * S_X2 * (three_j_F_X2_mF_X2*C_l3_X3_I*C_l1_X1_TX2 + three_j_0_mF_X2*C_l3_X3_TX2*C_l1_X1_I);
              
                  /* Perform the sum over the three possible positions of the second-order perturbation */
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = (T_X3_term + T_X2_term + T_X1_term);
                  
                  /* Comupte the reduced bispectrum */
                  pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] /= value[thread][l1_l2_l3_0_0][l3 - min[thread][l1_l2_l3_0_0]];
                
                  /* Check that for <TTT> the correction is equal to 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) */
                  if ((pbi->has_bispectra_t==_TRUE_) && (X==pbi->index_bf_t) && (X==Y) && (X==Z)) {
                    double exact = 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
                    double diff = fabs (1-pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3]/exact);
                    class_test_parallel (diff > _SMALL_,
                     pbi->error_message,
                     "quadratic bispectrum correction does not simplify to 8 (C_l1*C_l2 + perm); b=%g, exact=%g, diff=%g",
                     pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3], exact, diff);
                  }
                  
                  end_of_quadratic_bispectrum: ;
              
                  /* Update the counter */
                  #pragma omp atomic
                  pbi->count_memorised_for_bispectra++;
              
                  /* Debug. Print the quadratic contribution */
                  // printf ("B_quadratic_%3s[%3d,%3d,%3d] = %g\n",
                  //   pbi->bfff_labels[X][Y][Z], l1, l2, l3,
                  //   pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3]);
                
                } // end of for(index_l3)
                
              } // end of quadratic correction
              
            } // end of for(index_l2)
    
            #pragma omp flush(abort)
    
          } // end of for(index_l1)
      
          if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
    
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)

  } // end of for(index_bt)

  /* Free 3j values array */
  if (pbi->has_quadratic_correction == _TRUE_) {
    for (thread=0; thread < number_of_threads; ++thread) {
      for (int i=0; i < n_geometrical_factors; ++i)
        free (value[thread][i]);
      free (value[thread]);
    }
    free (value);
  }

  return _SUCCESS_;

}







int bispectra_non_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Allocate arrays inside the integration workspace */
  class_call (bispectra_non_separable_workspace_init(
                ppr,
                pba,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);
  

  if (pbi->bispectra_verbose > 1)
    printf (" -> computing non-separable bispectra; r sampled %d times in [%g,%g]\n", pwb->r_size, pwb->r_min, pwb->r_max);
  
  
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* Skip the bispectrum if it not of the non-separable type */
    if (pbi->bispectrum_type[index_bt] != non_separable_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {

      pwb->X = X;

      // ==================================================================================================
      // =                              Determine the bispectrum to compute                               =
      // ==================================================================================================
    
      /* Which primordial shape is needed? */
      if ((pbi->has_galileon_model==_TRUE_) && (index_bt==pbi->index_bt_galileon_gradient))
        pwb->shape_function = bispectra_galileon_gradient;

      else if ((pbi->has_galileon_model==_TRUE_) && (index_bt==pbi->index_bt_galileon_time))
        pwb->shape_function = bispectra_galileon_time;

      /* Debug - Uncomment to compute a custom bispectrum; useful to test against the separable shapes */
      // pwb->shape_function = bispectra_orthogonal_model;

      if (pbi->bispectra_verbose > 0)
        printf(" -> computing %s bispectrum involving %s^(2)\n",
        pbi->bt_labels[index_bt], pbi->bf_labels[X]);

      /* Compute fist integral over k3 */
      class_call (bispectra_non_separable_integrate_over_k3(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    index_bt,
                    pbi->index_tt_of_bf[X],
                    pwb),
        pbi->error_message,
        pbi->error_message);
        
      for (int Y = 0; Y < pbi->bf_size; ++Y) {

        pwb->Y = Y;

        /* Compute second integral over k2 */
        class_call (bispectra_non_separable_integrate_over_k2(
                      ppr,
                      pba,
                      ppt,
                      pbs,
                      ptr,
                      ppm,
                      pbi,
                      index_bt,
                      pbi->index_tt_of_bf[Y],
                      pwb),
          pbi->error_message,
          pbi->error_message);

        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          pwb->Z = Z;

          if (pbi->bispectra_verbose > 1)
            printf("   \\ computing bispectrum %s_%s%s%s\n",
            pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);

          /* Compute the third integral over k1 */
          class_call (bispectra_non_separable_integrate_over_k1(
                        ppr,
                        pba,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_bt,
                        pbi->index_tt_of_bf[Z],
                        pwb),
            pbi->error_message,
            pbi->error_message);
      
          /* Compute the fourth and last integral over r */
          class_call (bispectra_non_separable_integrate_over_r(
                        ppr,
                        pba,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        pbi->bispectra[index_bt][X][Y][Z],
                        pwb),
            pbi->error_message,
            pbi->error_message);

        } // end of for(Z)
      } // end of for(Y)
    } // end of for(X)
  } // end of for(index_bt)
  
  return _SUCCESS_;
  
}






int bispectra_non_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{

  // ===================================================================================================
  // =                                    Prepare integration grid                                     =
  // ===================================================================================================


  // ------------------------------------------------------------
  // -                         Grid in r                        -
  // ------------------------------------------------------------
  
  /* We set the r-sampling as if it were a time sampling.  We do so because 'r' has the right dimensions, and
  it always appears in the argument of a Bessel function multiplying a wavemode, just as it was for conformal
  time in the line-of-sight integral.  This is the only place in the module where the background structure
  is accessed. */
  pwb->r_min = pbi->r_min;
  pwb->r_max = pbi->r_max;
  pwb->r_size = pbi->r_size;
    

  /* We decide to sample r linearly */
  class_alloc (pwb->r, pwb->r_size*sizeof(double), pbi->error_message);
  lin_space (pwb->r, pwb->r_min, pwb->r_max, pwb->r_size);
    
  /* Allocate & fill delta_r, the measure for the trapezoidal integration over r */
  class_alloc (pwb->delta_r, pwb->r_size * sizeof(double), pbi->error_message);

  /* Fill pwb->delta_r */
  pwb->delta_r[0] = pwb->r[1] - pwb->r[0];
      
  for (int index_r=1; index_r < pwb->r_size-1; ++index_r)
    pwb->delta_r[index_r] = pwb->r[index_r+1] - pwb->r[index_r-1];
      
  pwb->delta_r[pwb->r_size-1] = pwb->r[pwb->r_size-1] - pwb->r[pwb->r_size-2];



  // -----------------------------------------------------------------------
  // -                       Grid for the shape function                   -
  // -----------------------------------------------------------------------
  
  /* We shall sample the primordial shape function only in the k-points determined in the perturbation
  module, which are much sparsely sampled than those of the transfer functions. Hence, we are assuming
  that the shape function is a smooth function of (k1,k2,k3) */
  
  pwb->k_smooth_size = ppt->k_size[ppt->index_md_scalars];
  
  class_alloc (pwb->k_smooth_grid, pwb->k_smooth_size*sizeof(double), pbi->error_message);
  
  for (int index_k=0; index_k < pwb->k_smooth_size; ++index_k)
    pwb->k_smooth_grid[index_k] = ppt->k[ppt->index_md_scalars][index_k];
  
    
  
  
  // -----------------------------------------------------------------------
  // -                       Grid for the shape function                   -
  // -----------------------------------------------------------------------

  /* Here we set the integration limits on k3. These can be made equal to those of k1 and k2 (that is equal to
  the range where we computed the transfer functions) but we can do better than that. In fact, 
  the r-integral enforces the triangular condition |k1-k2| <= k3 <= k1+k2 which means that we can restrict
  our range to those configurations. However, we should not get too close to the triangular limits otherwise
  the integral becomes numerically unstable. */

  /* Sampling of the first-order transfer functions.  */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];
  pwb->k3_size_max = k_tr_size;
  
  class_alloc (pwb->index_k3_lower, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  class_alloc (pwb->index_k3_upper, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  class_alloc (pwb->k3_grid_size, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  
  for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
  
    double k1 = pwb->k_smooth_grid[index_k1];
    class_alloc (pwb->index_k3_lower[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
    class_alloc (pwb->index_k3_upper[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
    class_alloc (pwb->k3_grid_size[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
        
    for (int index_k2=0; index_k2 <= index_k1; ++index_k2) {

      double k2 = pwb->k_smooth_grid[index_k2];
      
      /* Uncomment to use the triangular condition. This will give an imprecise result, one
      needs to extend the range a bit, same as we do for the second-order bispectrum. */
      // double k_tr_min = fabs(k1-k2);
      // double k_tr_max = k1+k2;

      /* Uncomment to take the whole k3 range for the integration (safer) */
      double k_tr_min = k_tr[0];
      double k_tr_max = k_tr[k_tr_size-1];
      
      /* Find the index corresponding to k3_lower inside k_tr */
      int index_k3_lower = 0;
      while (k_tr[index_k3_lower] < k_tr_min) ++index_k3_lower;
      pwb->index_k3_lower[index_k1][index_k2] = index_k3_lower;

      /* Find the index corresponding to k3_upper inside ptr->k */
      int index_k3_upper = k_tr_size - 1;
      while (k_tr[index_k3_upper] > k_tr_max) --index_k3_upper;
      pwb->index_k3_upper[index_k1][index_k2] = index_k3_upper;
      
      /* Number of points in the k_tr grid between 'k3_lower' and 'k3_upper' */
      pwb->k3_grid_size[index_k1][index_k2] = index_k3_upper - index_k3_lower + 1;

      /* Some debug - print out the k3_grid list for a special configuration */      
      // if ((index_k1==100) && (index_k2>-1)) {
      //   fprintf (stderr, "k1[%d]=%.5e, k2[%d]=%.5e, k3_grid_size=%d, k3_min=%.5e, k3_max=%.5e\n",
      //     index_k1, k_tr[index_k1], index_k2, k_tr[index_k2], pwb->k3_grid_size[index_k1][index_k2], k_tr_min, k_tr_max);
      //   for (int index_k3=0; index_k3 < pwb->k3_grid_size[index_k1][index_k2]; ++index_k3)
      //     fprintf(stderr, "%d %.17f /\\ ", index_k3, k_tr[index_k3]);
      // 
      //   fprintf (stderr, "\n\n");
      // }

    
    } // end of for(k2)
  } // end of for(k1)



  // ========================================================================================
  // =                               Determine window function                              =
  // ========================================================================================
  
  /* Determine the window function for the interpolation of the k-space bispectrum in k1 and k2. We need
  a window function because the primordial bispectrum will usually contain products of power spectra 
  that diverge as k^-3 for k -> 0. */
    
  /* Window function will have pbi->k_smooth_size elements */
  class_alloc (pwb->k_window, pwb->k_smooth_size*sizeof(double), pbi->error_message);

  for (int index_k=0; index_k < pwb->k_smooth_size; ++index_k) {
    
    double k = pwb->k_smooth_grid[index_k];
    double pk = pbi->pk_pt[index_k];
    
    pwb->k_window[index_k] = pow(k,2);

  }

  /* Inverse window function will have ptr->k_size[ppt->index_md_scalars] elements */
  class_alloc (pwb->k_window_inverse, ptr->k_size[ppt->index_md_scalars]*sizeof(double), pbi->error_message);
  for (int index_k=0; index_k < ptr->k_size[ppt->index_md_scalars]; ++index_k) {

    double k = ptr->k[ppt->index_md_scalars][index_k];
    double pk = pbi->pk[index_k];

    pwb->k_window_inverse[index_k] = 1./pow(k,2);

  }



  // =========================================================================================
  // =                                    Allocate arrays                                    =
  // =========================================================================================
  
  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel private (thread)
  #ifdef _OPENMP
  number_of_threads = omp_get_num_threads();
  #endif
  
  /* We need a k3_grid per thread because it varies with k1 and k2 due to the triangular condition */
  class_alloc (pwb->k3_grid, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->delta_k3, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->integral_splines, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->interpolated_integral, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->f, number_of_threads*sizeof(double*), pbi->error_message);
    
  /* Beginning of parallel region */
  abort = _FALSE_;
  #pragma omp parallel shared(pwb,pbi) private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    /* Allocate the integration grid in k3 with the maximum possible number of k3-values. This is given by the
    number of k-values in ptr->k */
    class_alloc_parallel(pwb->k3_grid[thread], pwb->k3_size_max*sizeof(double), pbi->error_message);
    class_alloc_parallel(pwb->delta_k3[thread], pwb->k3_size_max*sizeof(double), pbi->error_message);
  
    /* Allocate memory for the interpolation arrays (used only for the k2 and k3 integrations) */
    class_alloc_parallel (pwb->integral_splines[thread], ptr->k_size[ppt->index_md_scalars]*sizeof(double), pbi->error_message);
    class_alloc_parallel (pwb->interpolated_integral[thread], ptr->k_size[ppt->index_md_scalars]*sizeof(double), pbi->error_message);
    class_alloc_parallel (pwb->f[thread], pwb->k_smooth_size*sizeof(double), pbi->error_message);
  
  } // end of parallel region
  
  if (abort == _TRUE_) return _FAILURE_;

  return _SUCCESS_;
  
} // end of bispectra_non_separable_workspace_init







int bispectra_non_separable_integrate_over_k3 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k3,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Integration grid */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];  

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  // =========================================================================================
  // =                             Allocate memory for I_l3(r,k1,k2)                         =
  // =========================================================================================
  
  /* Initialize counter */
  pwb->count_allocated_for_integral_over_k3 = 0;
  
  /* Allocate l3-level.  Note that, even if l3 must satisfy the triangular
  inequality, we allocate this level for all the possible l3 values.  We do so because this
  array is going to be used by all (l1,l2) computations that follow, which means that l3 will
  eventually cover all the allowed range. */
  class_alloc (pwb->integral_over_k3, pbi->l_size*sizeof(double ***), pbi->error_message);
  
  for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
    
    /* Allocate r-level */
    class_alloc (pwb->integral_over_k3[index_l3], pwb->r_size*sizeof(double **), pbi->error_message);
  
    /* Allocate 'k1' level */
    for (int index_r=0; index_r < pwb->r_size; ++index_r) {
  
      int k1_size = pwb->k_smooth_size;
      class_alloc (pwb->integral_over_k3[index_l3][index_r], k1_size*sizeof(double *), pbi->error_message);

      /* Allocate 'k2' level */
      for (int index_k1=0; index_k1<k1_size; ++index_k1) {
  
        int k2_size = index_k1 + 1;
        class_calloc (pwb->integral_over_k3[index_l3][index_r][index_k1], k2_size, sizeof(double), pbi->error_message);
        
        /* Increase memory counter */
        pwb->count_allocated_for_integral_over_k3 += k2_size;
      } // end of for(index_k1)
    } // end of for(index_r)
  } // end of for(index_l3)
    
  if (pbi->bispectra_verbose > 2)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the k3-integral array (k_size=%d)\n",
      pwb->count_allocated_for_integral_over_k3*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k3, pwb->k_smooth_size);
  



  // ===================================================================================================
  // =                               Compute the INT_l3(r, k1, k2)  integral                           =
  // ===================================================================================================
  
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k3 = 0;

  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
  
    #pragma omp for schedule (dynamic)
    for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {

      double k1 = pwb->k_smooth_grid[index_k1];
      double pk_1 = pbi->pk_pt[index_k1];
      
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k3 integral for k1=%g, index_k1=%d\n", pwb->k_smooth_grid[index_k1], index_k1);
  
      /* We only need to consider those k2's that are equal to or larger than k1,
        as the shape function is assumed to be symmetric woth respect to k1<->k2<->k3 */      
      for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
    
        double k2 = pwb->k_smooth_grid[index_k2];
        double pk_2 = pbi->pk_pt[index_k2];
    
        /* Get the size of the integration grid. Note that when extrapolation is turned on, the k3-grid will also
        include values that do not satisfty the triangular condition k1 + k2 = k3. */
        int k3_size = pwb->k3_grid_size[index_k1][index_k2];
        int index_k3_lower = pwb->index_k3_lower[index_k1][index_k2];
        int index_k3_upper = pwb->index_k3_upper[index_k1][index_k2];
  
        /* Determine the integration grid. This is given by the portion of the transfer function grid
        starting at 'index_k3_lower' and ending at 'index_k3_upper' */
        for (int index_k3=index_k3_lower; index_k3<=index_k3_upper; ++index_k3)
          pwb->k3_grid[thread][index_k3-index_k3_lower] = k_tr[index_k3];

        /* If there are no points in k_tr that satisfy the triangular condition for the current (k1,k2), then
          there is no contribution to the integral */
        class_test_parallel (k3_size<=0, pbi->error_message, "include when triangular condition does not fit with ptr->k");
          
        /* Determine the measure for the trapezoidal rule for k3 */  
        pwb->delta_k3[thread][0] = pwb->k3_grid[thread][1] - pwb->k3_grid[thread][0];
  
        for (int index_k3=1; index_k3<(k3_size-1); ++index_k3)
          pwb->delta_k3[thread][index_k3] = pwb->k3_grid[thread][index_k3 + 1] - pwb->k3_grid[thread][index_k3 - 1];
  
        pwb->delta_k3[thread][k3_size-1] = pwb->k3_grid[thread][k3_size - 1] - pwb->k3_grid[thread][k3_size - 2];
        
        /* Shape function for this (k1,k2) slice */
        for (int index_k3=index_k3_lower; index_k3 <= index_k3_upper; ++index_k3) {
          
          double k3 = k_tr[index_k3];
          double pk_3 = pbi->pk[index_k3];

          class_call_parallel (pwb->shape_function (
                        ppm,
                        pbi,
                        k1, k2, k3,
                        pk_1, pk_2, pk_3,
                        &pwb->interpolated_integral[thread][index_k3]),
            pbi->error_message,
            pbi->error_message);
          
        }
        
  
        /* We compute the integral over k3 for all possible l-values */
        for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

          /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
          int tt_size = ptr->tt_size[ppt->index_md_scalars];
          int l_size = ptr->l_size[ppt->index_md_scalars];

          double * transfer = &(ptr->transfer
            [ppt->index_md_scalars]
            [((ppt->index_ic_ad * tt_size + index_tt_k3) * l_size + index_l3) * k_tr_size]);

          /* Some debug - print transfer function */
          // for (index_k3=0; index_k3 < k3_size; ++index_k3)
          //   if ((index_k1==3) && (index_k2==2) && (index_l3==0))
          //     fprintf(stderr, "%g %g\n", pwb->k3_grid[thread][index_k3], transfer[index_k3]);
    
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr + index_k3_lower,
                          pwb->delta_k3[thread],
                          k3_size,
                          transfer + index_k3_lower,
                          pwb->interpolated_integral[thread] + index_k3_lower,
                          index_l3,
                          pwb->r[index_r],
                          &(pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
                

            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k3;
  
            #pragma omp flush(abort)
  
          } // end of for(index_r)          
        } // end of for(index_l3)
      } // end of for(index_k2)
    } // end of for(index_k1)
  } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */    
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k3-integral array\n",
      pwb->count_memorised_for_integral_over_k3*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k3);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k3 != pwb->count_allocated_for_integral_over_k3,
              pbi->error_message,
              "there is a mismatch between allocated (%ld) and used (%ld) space!",
                pwb->count_allocated_for_integral_over_k3, pwb->count_memorised_for_integral_over_k3);

  return _SUCCESS_;
  
}








int bispectra_non_separable_integrate_over_k2 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k2,
    struct bispectra_workspace_non_separable * pwb
    )
{

  /* Integration grid */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];
    
  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  // ======================================================================================================
  // =                                   Allocate memory for INT_l3_l2(r,k1)                              =
  // ======================================================================================================
  
  /* Because we shall recycle the array for more than one fields (T,E...) we make sure we allocate it
  only once. This is achieved by performing the allocation only at the beginning of the loop over
  the field (that is, when pwb->Y==0) */

  if (pwb->Y == 0) {
  
    /* Initialize counter */
    pwb->count_allocated_for_integral_over_k2 = 0;
  
    /* Allocate l3-level */
    class_alloc (pwb->integral_over_k2, pbi->l_size*sizeof(double ***), pbi->error_message);
  
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
  
      /* Allocate l2-level. We only need l2<=l3 because of the k2<->k3 symmetry of the shape function */
      class_alloc (pwb->integral_over_k2[index_l3], (index_l3+1)*sizeof(double **), pbi->error_message);
  
      for (int index_l2=0; index_l2<=index_l3; ++index_l2) {
    
        /* Allocate r-level */
        class_alloc (pwb->integral_over_k2[index_l3][index_l2], pwb->r_size*sizeof(double *), pbi->error_message);
  
        /* Allocate 'k1' level */
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {
  
          int k1_size = pwb->k_smooth_size;
          class_alloc (pwb->integral_over_k2[index_l3][index_l2][index_r], k1_size*sizeof(double), pbi->error_message);
  
          /* Increase memory counter */
          pwb->count_allocated_for_integral_over_k2 += k1_size;
  
        } // end of for(index_k1)
      } // end of for(index_r)
    } // end of for(index_l2)
    
    if (pbi->bispectra_verbose > 2)
      printf("     * allocated ~ %.3g MB (%ld doubles) for the k2-integral array\n",
        pwb->count_allocated_for_integral_over_k2*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k2);
  
  } // end of if (pwb->Y==0)


  // ==============================================================================================================
  // =                                   Compute the INT_l3_l2(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k2 = 0;
    
  abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k2 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);
  
      for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {
          
        for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
      
          /* Interpolate the integral I_l2(k1,k2,r) that we computed above in the integration grid of k2.
          Note that we pass the integral_splines, interpolated_integral and pwb->f arrays separately rather
          than accessing them from pwb, because they are thread dependent (while pwb isn't). */
          class_call_parallel (bispectra_non_separable_interpolate_over_k2(
                      ppr,
                      ppt,
                      pbs,
                      ptr,
                      ppm,
                      pbi,
                      index_r,
                      index_k1,
                      index_l3,
                      pwb->integral_splines[thread],
                      pwb->interpolated_integral[thread],
                      pwb->f[thread],
                      pwb),
            pbi->error_message,
            pbi->error_message);
      
          for (int index_l2 = 0; index_l2 <= index_l3; ++index_l2) {  
  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k2) * l_size + index_l2) * k_tr_size]);

            /* Some debug - print transfer function */
            // if (index_r==0)
            //   for (index_k3=0; index_k3 < k_tr_size; ++index_k3)
            //     if ((index_k1==3) && (index_l3==2) && (index_l2==0))
            //       fprintf(stderr, "%g %g\n", k_tr[index_k3], transfer[index_k3]);

            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l2,
                          pwb->r[index_r],
                          &(pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k2;
  
            #pragma omp flush(abort)
  
          } // end of for(index_l2)
        } // end of for(index_l3)
      } // end of for(index_k1)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k2-integral array\n",
      pwb->count_memorised_for_integral_over_k2*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k2);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k2 != pwb->count_allocated_for_integral_over_k2,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k2, pwb->count_memorised_for_integral_over_k2);
  
  /* Free the memory that was allocated for the I_l_k3 integral, but only if we have already computed it
  for all the required probes */
  if (pwb->Y == (pbi->bf_size-1)) {
    for (int index_l2=0; index_l2 < pbi->l_size; ++index_l2) {
      for (int index_r=0; index_r < pwb->r_size; ++index_r) {
        for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
          free (pwb->integral_over_k3[index_l2][index_r][index_k1]);
        } // end of for(index_k1)
        free (pwb->integral_over_k3[index_l2][index_r]);
      } // end of for(index_r)
      free (pwb->integral_over_k3[index_l2]);
    } // end of for(index_l2)    
    free (pwb->integral_over_k3);
  }
    
  return _SUCCESS_;
}








int bispectra_non_separable_interpolate_over_k2 (
    struct precision * ppr,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_k1,
    int index_l3,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  int index_k, index_k_tr, index_k2;

  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];


  /* So far, we always assumed that k1>=k2 because the shape function is symmetric wrt k1<->k2.
  Interpolating an array with such property is complicated, hence we build a temporary array
  where f(index_k1, index_k2) = f(index_k2, index_k1) when index_k1 < index_k2. */
  for (index_k2=0; index_k2 < k_pt_size; ++index_k2) {

    f[index_k2] = (index_k1 > index_k2 ? pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]:
                                         pwb->integral_over_k3[index_l3][index_r][index_k2][index_k1]);

    /* Multiply by window function */
    f[index_k2] *= pwb->k_window[index_k2];

  }
  
  
  if (ppr->transfers_k2_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider (desired size of the slow index) */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;

    /* Interpolate for each value of l3, r, k1 */
    if (ppr->transfers_k2_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k2_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1] + ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Revert the effect of the window function */
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  } // end of for (index_k_tr)


  /* Some debug - print the original array and the interpolation */
  // if ((index_k1==1) && (index_l3==0) && (index_r==0)) {
  // 
  //   fprintf (stderr, "\n\n");
  // 
  //   /* Node points before window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]/pwb->k_window[index_k]);
  // 
  //   fprintf (stderr, "\n");
  //   
  //   /* Node points after window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]);
  // 
  //   fprintf (stderr, "\n");
  // 
  //   /* Interpolation after inverse window */  
  //   for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  //  }

  return _SUCCESS_;

}






int bispectra_non_separable_integrate_over_k1 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k1,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Integration grid */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  // ==========================================================================================================
  // =                                    Allocate memory for INT_l3_l2_l1(r)                                 =
  // ==========================================================================================================
  
  /* We shall allocate the array pwb->integral_over_k1[index_l_k2][index_l_k3][index_l_k1-index_l_k1_min][index_r]
  so that the l_k1 level is the one satisfying the triangular inequality (|l_k2-l_k3| <= l_k1 <= l_k2+l_k3). 
  pwb->integral_over_k1 is recycled by the different iterations in Z-field loop, hence we allocate
  it only at the first iteration (pwb->Z==0) */
  
  if (pwb->Z==0) {

    /* The integral over k1 yields a function of (l3,l2,l1) that is symmetric under permutations of (l3,l2,l1).
    Hence, we compute and store it only for l3>=l2>=l1 configurations (that satisfy the triangular condition) */
    pwb->count_allocated_for_integral_over_k1 = pwb->r_size * pbi->n_independent_configurations;
  
    /* Allocate (l3,l2,l1)-level */
    class_alloc (pwb->integral_over_k1, pbi->n_independent_configurations*sizeof(double *), pbi->error_message);
  
    for (long int index_l3_l2_l1=0; index_l3_l2_l1 < pbi->n_independent_configurations; ++index_l3_l2_l1)
      class_alloc (pwb->integral_over_k1[index_l3_l2_l1], pwb->r_size*sizeof(double), pbi->error_message);
    
    if (pbi->bispectra_verbose > 2)
      printf("     * allocated ~ %.3g MB (%ld doubles) for the k1-integral array\n",
        pwb->count_allocated_for_integral_over_k1*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k1);
  }
  
  // ==========================================================================================================
  // =                                    Compute  INT_l3_l2_l1(r)  integral                                  =
  // ==========================================================================================================

  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k1 = 0;
  
  /* As for the other integrals, we parallelize the loop over 'r'. */
  abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k1 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);

      for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
    
        for (int index_l2 = 0; index_l2 <= index_l3; ++index_l2) {
  
          /* Skip those configurations that are forbidden by the triangular condition (optimization) */
          if (pbi->l[index_l2] < pbi->l[index_l3]/2)
            continue;

          /* Interpolate the integral I_l3_l2(k1,r) that we computed above in the integration grid of k1 */
          class_call_parallel (bispectra_non_separable_interpolate_over_k1 (
                        ppr,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_r,
                        index_l3,
                        index_l2,
                        pwb->integral_splines[thread],
                        pwb->interpolated_integral[thread],
                        pwb->f[thread],
                        pwb),
            pbi->error_message,
            pbi->error_message);      
  
          /* Determine the limits for l1, which come from the triangular inequality |l3-l2| <= l1 <= l3+l2 */
          int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
          int index_l1_max = MIN (index_l2, pbi->index_l_triangular_max[index_l3][index_l2]);
  
          for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {  

            /* Index of the current (l3,l2,l1) configuration */
            long int index_l3_l2_l1 = pbi->index_l1_l2_l3[index_l3][index_l3-index_l2][index_l1_max-index_l1];
  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k1) * l_size + index_l1) * k_tr_size]);

            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l1,
                          pwb->r[index_r],
                          &(pwb->integral_over_k1[index_l3_l2_l1][index_r]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k1;
              
            #pragma omp flush(abort)
  
          } // end of for(index_l1)
        } // end of for(index_l2)
      } // end of for(index_l3)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k1-integral array\n",
      pwb->count_memorised_for_integral_over_k1*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k1);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k1 != pwb->count_allocated_for_integral_over_k1,
              pbi->error_message,
              "there is a mismatch between allocated (%ld) and used (%ld) space!",
                pwb->count_allocated_for_integral_over_k1, pwb->count_memorised_for_integral_over_k1);
  
  /* Free the memory that was allocated for the integral over k2, but do that only when we are at
  the last iteration of the Z and Y loops. */
  if ((pwb->Y == (pbi->bf_size-1)) && (pwb->Z == (pbi->bf_size-1))) {
    for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
      for (int index_l1=0; index_l1<=index_l2; ++index_l1) {
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {      
          free (pwb->integral_over_k2[index_l2][index_l1][index_r]);
        } // end of for(index_r)
        free (pwb->integral_over_k2[index_l2][index_l1]);
      } // end of for(index_l1)
      free (pwb->integral_over_k2[index_l2]);
    } // end of for(index_l2)
    free (pwb->integral_over_k2);
  }
  
  return _SUCCESS_;
  
}








int bispectra_non_separable_interpolate_over_k1 (
    struct precision * ppr,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_l3,
    int index_l2,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  int index_k, index_k_tr, index_k1;
  
  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];
  
  /* Define the function to be interpolated, and multiply it by a window function */
  for (index_k1=0; index_k1 < k_pt_size; ++index_k1) {
    
    f[index_k1] = pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1];
    f[index_k1] *= pwb->k_window[index_k1];
  }
  
  if (ppr->transfers_k1_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider (desired size of the slow index) */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;
      
    /* Interpolate for each value of l3, l2, r */
    if (ppr->transfers_k1_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k1_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1] + ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Revert the effect of the window function */
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  } // end of for (index_k_tr)


  /* Some debug - print the original array and the interpolation */
  // if ((index_l3==48) && (index_l2==35) && (index_r==50)) {
  // 
  //   fprintf (stderr, "\n\n");
  // 
  //   /* Node points before window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]/pwb->k_window[index_k]);
  // 
  //   fprintf (stderr, "\n");
  //   
  //   /* Node points after window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]);
  // 
  //   fprintf (stderr, "\n");
  // 
  //   /* Interpolation after inverse window */  
  //   for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  //  }

  return _SUCCESS_;
  
}









int bispectra_non_separable_integrate_over_r (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    double * bispectrum,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* We now proceed to the final step of the bispectrum computation where we integrate over
  the r-dependence and obtain the bispectrum. We also multiply the result by the appropriate
  coefficients.  */  
  
  /* We parallelize the outer loop over 'l1'. */
  int abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)
  {
  
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      if (pbi->bispectra_verbose > 2)
        printf("     * computing the r-integral for l1=%d, index_l1=%d\n", pbi->l[index_l1], index_l1);
    
      for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (pbi->l[index_l2] < pbi->l[index_l1]/2)
          continue;

        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];    
          double * I = pwb->integral_over_k1[index_l1_l2_l3];
          double integral = 0;

          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
            double r = pwb->r[index_r];
            double integrand = r*r * I[index_r];
      
            if (integrand == 0.)
              continue;
      
            integral += integrand * pwb->delta_r[index_r];

            /* Some debug - output intermediate results on stderr for a custom (l1,l2,l3) configuration */
            // if ( (l1==100) && (l2==200) && (l3==300) ) {
            //   if (index_r==0) {
            //     fprintf(stderr, "##########    l1 = %d, l2 = %d, l3 = %d, n_rows = %d    ##########\n\n",
            //       l1, l2, l3, r, pwb->r_size);
            //     fprintf(stderr, "%12s %17s %17s %17s\n", "r", "r^2*I_l1_l2_l3(r)", "integral","delta_r");
            //   }
            //   else {
            //     fprintf(stderr, "%12.7g %17.7g %17.7g %17.7g\n", r, integrand, integral, pwb->delta_r[index_r]);
            //   }
            // }
   
          } // end of for(index_r)

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3) with 1/2 from trapezoidal rule */
          bispectrum[index_l1_l2_l3] = 0.5 * integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. This factor is seen
          (see, for instance, eq. 17 of Fergusson & Shellard 2007). */
          bispectrum[index_l1_l2_l3] *= pow(2./_PI_,3);

          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_for_bispectra++;

          /* Some debug - output the integral as a function of r on stderr for a custom (l1,l2,l3) */
          // if ( (l1==l2) && (l2==l3) ) {
          //   fprintf(stderr, "%12d %17.7g\n", l1, pwb->integral_over_r[index_l1][index_l2][index_l3-index_l3_min]);
          // }

        } // end of for(index_l3)
      } // end of for(index_l2)
      
      #pragma omp flush(abort)
      
    } // end of for(index_l1)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  /* Free the memory that was allocated for the integral over k1, but do that only when we are at
  the last iteration of the Z and Y . */
  if ((pwb->Z == (pbi->bf_size-1))) {
    for (long int index_l1_l2_l3 = 0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
      free (pwb->integral_over_k1[index_l1_l2_l3]);
    free (pwb->integral_over_k1);
  }
    
  return _SUCCESS_;
  
}




int bispectra_non_separable_workspace_free (
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  free (pwb->r);
  free (pwb->delta_r);  
  free (pwb->k_smooth_grid);

  for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
    free (pwb->index_k3_lower[index_k1]);
    free (pwb->index_k3_upper[index_k1]);
    free (pwb->k3_grid_size[index_k1]);    
  }
  free (pwb->index_k3_lower);
  free (pwb->index_k3_upper);
  free (pwb->k3_grid_size);

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  #pragma omp parallel private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    free(pwb->k3_grid[thread]);
    free(pwb->delta_k3[thread]);
    free(pwb->integral_splines[thread]);
    free(pwb->interpolated_integral[thread]);
    free(pwb->f[thread]);
    
  }  if (abort == _TRUE_) return _FAILURE_;
  
  free(pwb->k3_grid);
  free(pwb->delta_k3);
  free(pwb->integral_splines);
  free(pwb->interpolated_integral);
  free(pwb->f);
  
  return _SUCCESS_;
    
}




/**
 * Save a bispectrum to disk. It is important to keep the index_bt dependence in the argument list,
 * as we might want to selectively save only the secondary bispectra.
 */
int bispectra_save_to_disk (
    struct bispectra * pbi,
    int index_bt
    )
{

  /* Open file for writing */
  class_open (pbi->bispectra_run_files[index_bt], pbi->bispectra_run_paths[index_bt], "a+b", pbi->error_message);

  /* Print some debug */
  if (pbi->bispectra_verbose > 2)
    printf("     * writing bispectra to disk for index_bt=%d on '%s'\n",
      index_bt, pbi->bispectra_run_paths[index_bt]);

  /* Write all the independent (l1,l2,l3) triplets for this bispectrum */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        fwrite(
              pbi->bispectra[index_bt][X][Y][Z],
              sizeof(double),
              pbi->n_independent_configurations,
              pbi->bispectra_run_files[index_bt]
              );

  /* Close file */
  fclose(pbi->bispectra_run_files[index_bt]);
  
  return _SUCCESS_;
  
}






/**
 * Load a bispectrum from disk. It is important to keep the index_bt dependence in the argument list,
 * as we might want to selectively load only the secondary bispectra.
 */
int bispectra_load_from_disk(
    struct bispectra * pbi,
    int index_bt
    )
{

  /* Open file for reading */
  class_open (pbi->bispectra_run_files[index_bt], pbi->bispectra_run_paths[index_bt], "rb", pbi->error_message);

  /* Print some debug */
  if (pbi->bispectra_verbose > 2)
    printf("     * reading bispectra from disk for index_bt=%d on'%s'\n", index_bt, pbi->bispectra_run_paths[index_bt]);

  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {

        int n_to_read = pbi->n_independent_configurations;
  
        /* Read a chunk with all the independent (l1,l2,l3) triplets to pbi->bispectra[index_bt] */
        int n_read = fread(
                pbi->bispectra[index_bt][X][Y][Z],
                sizeof(double),
                n_to_read,
                pbi->bispectra_run_files[index_bt]);

        class_test(n_read != n_to_read,
          pbi->error_message,
          "Could not read in '%s' file, read %d entries but expected %d",
            pbi->bispectra_run_paths[index_bt], n_read, n_to_read);        
      }
    }
  }
  
  /* Close file */
  fclose(pbi->bispectra_run_files[index_bt]); 

  return _SUCCESS_;
  
}



/**
 * Bispectrum produced by pi_dot * grad_pi^2 in the Galileon Lagrangian, from eq. 22 of arXiv:0905.3746v3
 *
 */
int bispectra_galileon_gradient (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3,
  double pk_1, double pk_2, double pk_3,
  double * out
  )
{ 
  
  double K1 = k1 + k2 + k3;
  double K1_squared = K1*K1;
  double K1_cubed = K1_squared*K1;
  double K1_fourth = K1_squared*K1_squared;
  double K1_sixth = K1_fourth*K1_squared;
  
  double K2_squared = k1*k2 + k2*k3 + k3*k1;
  double K2_fourth = K2_squared*K2_squared;
  
  double K3_cubed = k1*k2*k3;
  double K3_sixth = K3_cubed*K3_cubed;
  double K3_ninth = K3_sixth*K3_cubed;
  
  *out = -27/17. *

      (

        + 24*K3_sixth - 8*K2_squared*K3_cubed*K1 - 8*K2_fourth*K1_squared
        + 22*K3_cubed*K1_cubed - 6*K2_squared*K1_fourth + 2*K1_sixth

      ) / (K3_ninth*K1_cubed);
  
  /* Multiply by the amplitude of primordial fluctuations. The factor -3/5 derives by the
  fact that we defined the transfer functions wrt R, and that a bispectrum with fnl_phi=1
  (which is what we want) is equivalent to a bispectrum with fnl_R=-3/5.  */
  double fnl_R = -3/5.; 
  *out *= fnl_R * (2*_PI_*_PI_*ppm->A_s) * (2*_PI_*_PI_*ppm->A_s);
  
  return _SUCCESS_; 
  
}



/**
 * Bispectrum produced by pi_dot^3 in the Galileon Lagrangian, from eq. 22 of arXiv:0905.3746v3
 *
 */
int bispectra_galileon_time (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double K1 = k1+k2+k3;
  double K1_cubed = K1*K1*K1;
  double K3_cubed = k1*k2*k3;
  
  *out = 162 / (K3_cubed * K1_cubed);
  
  /* Multiply by the amplitude of primordial fluctuations */
  double fnl_R = -3/5.; 
  *out *= fnl_R * (2*_PI_*_PI_*ppm->A_s) * (2*_PI_*_PI_*ppm->A_s);
  
  return _SUCCESS_; 
  
}



/**
 * Local template bispectrum. Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_local_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  *out = 2 * (-3/5.) * ( pk_1*pk_2 + pk_1*pk_3 + pk_2*pk_3 );
  
  return _SUCCESS_; 
  
}


/**
 * Equilateral template bispectrum. Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_equilateral_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double pk_one_third_1 = pow(pk_1, one_third);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, one_third);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, one_third);
  double pk_two_thirds_3 = pk_one_third_3*pk_one_third_3;

  
  *out = 6 * (-3/5.) * (
      
      - pk_1 * pk_2
      - pk_1 * pk_3
      - pk_2 * pk_3
        
      - 2 * pk_two_thirds_1 * pk_two_thirds_2 * pk_two_thirds_3
        
      + pk_1 * pk_one_third_2 * pk_two_thirds_3
      + pk_1 * pk_one_third_3 * pk_two_thirds_2
      + pk_2 * pk_one_third_1 * pk_two_thirds_3
      + pk_2 * pk_one_third_3 * pk_two_thirds_1
      + pk_3 * pk_one_third_1 * pk_two_thirds_2
      + pk_3 * pk_one_third_2 * pk_two_thirds_1
        
  );
  

  return _SUCCESS_; 
  
}
 
 
   
/**
 * Orthogonal template bispectrum (arXiv:0905.3746). Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_orthogonal_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double pk_one_third_1 = pow(pk_1, one_third);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, one_third);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, one_third);
  double pk_two_thirds_3 = pk_one_third_3*pk_one_third_3;

  
  *out = 6 * (-3/5.) * (
      
      - 3 * pk_1 * pk_2
      - 3 * pk_1 * pk_3
      - 3 * pk_2 * pk_3
        
      - 8 * pk_two_thirds_1 * pk_two_thirds_2 * pk_two_thirds_3
        
      + 3 * pk_1 * pk_one_third_2 * pk_two_thirds_3
      + 3 * pk_1 * pk_one_third_3 * pk_two_thirds_2
      + 3 * pk_2 * pk_one_third_1 * pk_two_thirds_3
      + 3 * pk_2 * pk_one_third_3 * pk_two_thirds_1
      + 3 * pk_3 * pk_one_third_1 * pk_two_thirds_2
      + 3 * pk_3 * pk_one_third_2 * pk_two_thirds_1
        
  );
  

  return _SUCCESS_; 
  
}








