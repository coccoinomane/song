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


  // =====================================================================================
  // =                                    Preparations                                   =
  // =====================================================================================

  /* Initialize indices & arrays in the bispectra structure */

  class_call (bispectra_indices (ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);





  // =====================================================================================
  // =                                 Compute bispectra                                 =
  // =====================================================================================
  
  class_call (bispectra_harmonic (ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi),
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
     struct precision * ppr,
     struct perturbs * ppt,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi
     )
{

  if (pbi->has_bispectra == _TRUE_) {

    free(pbi->l);
    free(pbi->pk);
    free(pbi->pk_pt);

    for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {

      free(pbi->l_triangular_size[index_l1]);
      free(pbi->index_l_triangular_min[index_l1]);
      free(pbi->index_l_triangular_max[index_l1]);
    
    } // end of for(index_l1)

    free(pbi->l_triangular_size);
    free(pbi->index_l_triangular_min);
    free(pbi->index_l_triangular_max);

    for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {
      for(int index_l2=0; index_l2<=index_l1; ++index_l2)
          free (pbi->index_l1_l2_l3[index_l1][index_l1-index_l2]);
      free (pbi->index_l1_l2_l3[index_l1]);
    } free (pbi->index_l1_l2_l3);

    /* Free pbi->bispectra */
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      for (int X = 0; X < pbi->bf_size; ++X) {
        for (int Y = 0; Y < pbi->bf_size; ++Y) {
          for (int Z = 0; Z < pbi->bf_size; ++Z)
            free (pbi->bispectra[index_bt][X][Y][Z]);
          free (pbi->bispectra[index_bt][X][Y]);
        } free (pbi->bispectra[index_bt][X]);
      } free (pbi->bispectra[index_bt]);
    } free (pbi->bispectra);

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
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->d_lsq_cls[index_ct]);
    free (pbi->d_lsq_cls);
    if (pbi->include_lensing_effects == _TRUE_) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_d_lsq_cls[index_lt]);
      free (pbi->lensed_d_lsq_cls);
    }
    
    if (pbi->include_lensing_effects == _TRUE_) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_cls[index_lt]);
      free (pbi->lensed_cls);
    }
      
    /* Free file arrays */
    if ((ppr->store_bispectra_to_disk == _TRUE_) || (ppr->load_bispectra_from_disk == _TRUE_)) {
    
      // fclose(pbi->bispectra_status_file);
    
      for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt)
        free (pbi->bispectra_paths[index_bt]);
    
      free (pbi->bispectra_files);
      free (pbi->bispectra_paths);
    
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
  // =                                   Count bispectra fields                                     =
  // ================================================================================================

  /* Find out which kind of bispectra to compute and assign them indices and labels
  Generate indices for the probes (T for temperature, E for E-mode polarisation,
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

  class_test (pbi->bf_size > 2,
    pbi->error_message,
    "we cannot compute the bispectrum for more than %d fields (e.g. T and E), reduce your expectations :-)", pbi->bf_size);

  class_test (pbi->bf_size < 1,
    pbi->error_message,
    "no probes requested");

  /* Create labels for the full bispectra */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        sprintf (pbi->bfff_labels[X][Y][Z], "%s%s%s",
        pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);


  /* Associate to each field X=T,E,... its transfer function, which was computed in the transfer.c
  module, the C_l correlation <X phi> with the lensing potential, the C_l correlation <X zeta> with
  the curvature perturbation and, to each possible pair of fields (TT,EE,TE,...), their power spectra,
  which were computed in the spectra.c module. */
  for (int X = 0; X < pbi->bf_size; ++X) {
    
    if ((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_t;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_tp;
      pbi->index_ct_of_zeta_bf[X] = psp->index_ct_tz;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_tt;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_tt;
      if (pbi->include_lensing_effects == _TRUE_)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_tt;
    }

    if ((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_e;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_ep;
      pbi->index_ct_of_zeta_bf[X] = psp->index_ct_ez;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_te;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_ee;
      if (pbi->include_lensing_effects == _TRUE_)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_ee;
    }

    if ((pbi->has_bispectra_r == _TRUE_) && (X == pbi->index_bf_r)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_r;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_tr;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_rr;
      /* TODO: lensed Rayleigh not implemented yet */
      // pbi->index_ct_of_phi_bf[X] = psp->index_ct_rp;
      // pbi->index_ct_of_zeta_bf[X] = psp->index_ct_rz;
      // if (pbi->include_lensing_effects == _TRUE_)
      //   pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_rr;
    }

    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      if (((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_e == _TRUE_) && (Y == pbi->index_bf_e))) {

        pbi->index_ct_of_bf_bf[X][Y] = pbi->index_ct_of_bf_bf[Y][X] = psp->index_ct_te;
        if (pbi->include_lensing_effects == _TRUE_)
          pbi->index_lt_of_bf_bf[X][Y] = pbi->index_lt_of_bf_bf[Y][X] = ple->index_lt_te;
      }

      if (((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r))) {

        pbi->index_ct_of_bf_bf[X][Y] = pbi->index_ct_of_bf_bf[Y][X] = psp->index_ct_tr;
        /* TODO: lensed Rayleigh not implemented yet */
        // if (pbi->include_lensing_effects == _TRUE_)
        //   pbi->index_lt_of_bf_bf[X][Y] = pbi->index_lt_of_bf_bf[Y][X] = ple->index_lt_tr;
      }

      if (((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r)))
        class_test (1==1, pbi->error_message, "polarization-Rayleigh bispectrum not implemented yet");

      if (((pbi->has_bispectra_b == _TRUE_) && (X == pbi->index_bf_b))
       && ((pbi->has_bispectra_r == _TRUE_) && (Y == pbi->index_bf_r)))
        class_test (1==1, pbi->error_message, "polarization-Rayleigh bispectrum not implemented yet");
    }
  }

  // ================================================================================================
  // =                                    Count bispectra types                                     =
  // ================================================================================================
  
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
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_equilateral_model) {
    pbi->index_bt_equilateral = index_bt;
    strcpy (pbi->bt_labels[index_bt], "equilateral");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_orthogonal_model) {
    pbi->index_bt_orthogonal = index_bt;
    strcpy (pbi->bt_labels[index_bt], "orthogonal");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  // *** Non-separable bispectra

  if (pbi->has_galileon_model) {

    /* Bispectrum induced by pi_dot*pi_grad^2 */
    pbi->index_bt_galileon_gradient = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_grad");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;

    /* Bispectrum induced by pi_dot^3 */
    pbi->index_bt_galileon_time = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_time");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  // *** Analytical bispectra

  if (pbi->has_local_squeezed == _TRUE_) {
    pbi->index_bt_local_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local(sqz)");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_intrinsic_squeezed == _TRUE_) {
    pbi->index_bt_intrinsic_squeezed = index_bt;
    if (pbi->lensed_intrinsic == _TRUE_)
      strcpy (pbi->bt_labels[index_bt], "intrinsic(sqz,lens)");
    else
      strcpy (pbi->bt_labels[index_bt], "intrinsic(sqz,unl)");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_intrinsic_squeezed_unlensed == _TRUE_) {
    pbi->index_bt_intrinsic_squeezed_unlensed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic(sqz,unl)");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_cosine_shape == _TRUE_) {
    pbi->index_bt_cosine = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cosine");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_cmb_lensing == _TRUE_) {
    pbi->index_bt_cmb_lensing = index_bt;
    strcpy (pbi->bt_labels[index_bt], "CMB-lensing");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }
  
  if (pbi->has_cmb_lensing_squeezed == _TRUE_) {
    pbi->index_bt_cmb_lensing_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "CMB-lensing(sqz)");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }
  
  /* The kernel for the squeezed CMB-lensing bispectrum is needed to compute
  the lensing contribution to the variance. In the final Fisher matrix, the kernel
  will be multiplied by C_l^{X\phi} to give the actual squeezed bispectrum (see
  eq. 5.20 of http://uk.arxiv.org/abs/1101.2234); it will therefore show up as
  'CMB-lensing(sqz)' */
  if (pbi->has_cmb_lensing_kernel == _TRUE_) {
    pbi->index_bt_cmb_lensing_kernel = index_bt;
    strcpy (pbi->bt_labels[index_bt], "CMB-lensing(sqz)");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }
  
  if (pbi->has_quadratic_correction == _TRUE_) {
    pbi->index_bt_quadratic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "quadratic");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    index_bt++;
  }

  // *** Intrinsic (i.e. second-order) bispectra

  if (pbi->has_intrinsic == _TRUE_) {
    pbi->index_bt_intrinsic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic");
    pbi->bispectrum_type[index_bt] = intrinsic_bispectrum;
    pbi->n[intrinsic_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_; /* will be overwritten in bispectra2.c */
    index_bt++;
  }

  pbi->bt_size = index_bt;

  class_test (pbi->bt_size > _MAX_NUM_BISPECTRA_,
   "exceeded maximum number of allowed bispectra, increase _MAX_NUM_BISPECTRA_ in common.h",
   pbi->error_message);

  /* Are the Wigner 3j-symbols needed to compute the requested bispectra? */
  pbi->need_3j_symbols = ((pbi->has_bispectra_e) &&
  ((pbi->has_quadratic_correction == _TRUE_)
  || (pbi->has_cmb_lensing == _TRUE_)
  || (pbi->has_cmb_lensing_squeezed == _TRUE_)
  || (pbi->has_cmb_lensing_kernel == _TRUE_)));


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
  // =                                  Assign window functions                                 =
  // ============================================================================================
  
  /* Assign to each bispectrum a window function suitable for its interpolation. See header file
  for mode details. Comment out a bispectrum to have it interpolated without a window function.
  
  It seems that using window functions that cross the zero might generate some issues. This
  has yet to be verified rigorously, though.
  
  */
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    pbi->window_function[index_bt] = NULL;
    
    /* Previously, SONG used to interpolate the bispectra in the two smallest l-multipoles. That
    was not a good idea, and we needed to use window functions. Now we interpolate the two 
    largest l's instead, and we don't need window functions anymore. */
    continue;

    /* In absence of reionisation and for polarisation, better results are obtained without a window
    function. The reason is that the C_l's for polarisation are tiny at small l's without
    reionisation, so multiplying and dividing by the window functions creates numerical
    instabilities. The patology is stronger for bispectra that peak in the squeezed limit,
    as in that case the large scales are those where most of the signal comes from. */
    /* TODO: the above statement has to be corrected in terms of the new interpolation,
    which relies on interpolating the two largest multipoles */
    if ((ppr->has_reionization == _FALSE_) && (pbi->has_bispectra_e == _TRUE_) && (pbi->bf_size > 1)) {
      continue;
    }
    else {
  
      if ((pbi->has_local_model == _TRUE_) && (index_bt == pbi->index_bt_local))
        pbi->window_function[index_bt] = bispectra_local_window_function;
          
      else if ((pbi->has_intrinsic == _TRUE_) && (index_bt == pbi->index_bt_intrinsic))
        pbi->window_function[index_bt] = bispectra_local_window_function;
          
      else if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed))
        pbi->window_function[index_bt] = bispectra_local_window_function;
    }
  
    /* For the non-squeezed bispectra, we always use the window function. */
    if ((pbi->has_equilateral_model == _TRUE_) && (index_bt == pbi->index_bt_equilateral))
      pbi->window_function[index_bt] = bispectra_local_window_function;
    
    else if ((pbi->has_orthogonal_model == _TRUE_) && (index_bt == pbi->index_bt_orthogonal))
      pbi->window_function[index_bt] = bispectra_local_window_function;
    
    else if ((pbi->has_galileon_model==_TRUE_) && (index_bt == pbi->index_bt_galileon_gradient))
      pbi->window_function[index_bt] = bispectra_local_window_function;
    
    else if ((pbi->has_galileon_model==_TRUE_) && (index_bt == pbi->index_bt_galileon_time))
      pbi->window_function[index_bt] = bispectra_local_window_function;
  
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
  
  
  /* Create the files to store the bispectra in */
  if ((ppr->store_bispectra_to_disk == _TRUE_) || (ppr->load_bispectra_from_disk == _TRUE_)) {

    /* We are going to store the bispectra in n=bt_size files, one for each requested type of bispectrum */
    class_alloc (pbi->bispectra_files, pbi->bt_size*sizeof(FILE *), pbi->error_message);
    class_alloc (pbi->bispectra_paths, pbi->bt_size*sizeof(char *), pbi->error_message);
  
    for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      
      /* Include the name of the bispectrum in its file */
      class_alloc (pbi->bispectra_paths[index_bt], _FILENAMESIZE_*sizeof(char), pbi->error_message);
      sprintf (pbi->bispectra_paths[index_bt], "%s/bispectra_%s.dat", pbi->bispectra_dir, pbi->bt_labels[index_bt]);
      
    } // end of for(index_bt)

    if (ppr->store_bispectra_to_disk == _TRUE_)
      if (pbi->bispectra_verbose > 1)
        printf ("     * will create %d files for the bispectra\n", pbi->bt_size);
      
  } // end of if(ppr->store_bispectra_to_disk)
  
  


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
                ppr,
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
    struct precision * ppr,
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
    class_calloc (pbi->cls[index_ct], pbi->full_l_size, sizeof(double), pbi->error_message);
  
  /* Do the same for the logarithmic derivative of the C_l's */
  class_alloc (pbi->d_lsq_cls, psp->ct_size*sizeof(double*), pbi->error_message);
  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    class_calloc (pbi->d_lsq_cls[index_ct], pbi->full_l_size, sizeof(double), pbi->error_message);

  /* If the the effect of lensing is to be included in the bispectra, interpolate also
  the lensed C_l's */ 
  if (pbi->include_lensing_effects == _TRUE_) {

    class_alloc (pbi->lensed_cls, ple->lt_size*sizeof(double*), pbi->error_message);
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
      class_calloc (pbi->lensed_cls[index_lt], pbi->full_l_size, sizeof(double), pbi->error_message);
    
    /* The squeezed approximation of the intrinsic bispectrum requires the computation of the
    C_l derivatives */
    class_alloc (pbi->lensed_d_lsq_cls, ple->lt_size*sizeof(double*), pbi->error_message);
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
      class_calloc (pbi->lensed_d_lsq_cls[index_lt], pbi->full_l_size, sizeof(double), pbi->error_message);
  }
  
  
  /* We shall call the CLASS function 'spectra_cl_at_l'. This gives three outputs:
  the total Cl's for each probe (T, E, B...); the Cl's divided by probe and mode
  (scalar, vector, tensor); the Cl's divided by probe, mode, and initial condition
  (adiabatic, isocurvature...). We have to allocate these three arrays before being
  able to call 'spectra_cl_at_l'. We do so copying what is done in the function
  'output_total_cl_at_l' in the output module */
  double * cl;        /* cl_md_ic[index_ct] */
  double ** cl_md;    /* cl_md[index_mode][index_ct] */
  double ** cl_md_ic; /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */
  class_alloc(cl, MAX(psp->ct_size,ple->lt_size)*sizeof(double), pbi->error_message);	
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
    

    /* Store the lensed C_l's */
    if (pbi->include_lensing_effects == _TRUE_) {

        /* The lensed C_l's must have been computed up to the maximum required multipole,
        that is, up to pbi->l_max. This check is made inside 'lensing_cl_at_l' */
        class_call(lensing_cl_at_l(
                     ple,
                     l,
                     cl),
          ple->error_message,
          pbi->error_message);
    
        for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
          pbi->lensed_cls[index_lt][l-2] = cl[index_lt];
        
        /* Debug - print temperature-lensing potential C_l's */
        // double factor = l*(l+1.)/(2*_PI_);
        // fprintf (stderr, "%4d %16g %16g %16g %16g\n",
        //   l, factor*pbi->cls[psp->index_ct_tt][l-2], factor*pbi->lensed_cls[ple->index_lt_tt][l-2],
        //   factor*sqrt(l*(l+1))*pbi->cls[psp->index_ct_tp][l-2], factor*sqrt(l*(l+1))*pbi->lensed_cls[ple->index_lt_tp][l-2]);
        
    }
    
    /* Uncomment to turn the CMB-lensing C_l to zero on small scales, where we cannot trust them.
    This won't change the result because these C_l's are very small for large l's. In CAMB, Antony
    sets C_l^TP=0 for l>300 and C_l^EP=0 for l>40. */
    if ((pbi->has_cmb_lensing == _TRUE_)
    || (pbi->has_cmb_lensing_squeezed == _TRUE_)
    || (pbi->has_cmb_lensing_kernel == _TRUE_)) {
      
      pbi->lmax_lensing_corrT = 300;
      pbi->lmax_lensing_corrE = 300;
      
      if ((l > pbi->lmax_lensing_corrT) && (pbi->has_bispectra_t == _TRUE_))
        pbi->cls[psp->index_ct_tp][l-2] = 0;

      if ((l > pbi->lmax_lensing_corrE) && (pbi->has_bispectra_e == _TRUE_))
        pbi->cls[psp->index_ct_ep][l-2] = 0;

    }

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
    
  } // end of for loop on the l's

  /* Compute the derivative of l*l*C_l from the spline-interpolated C_l's. This corresponds
  to the second method mentioned in the long comment above. */
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
  
  /* Do the same for the lensed C_l's */
  if (pbi->include_lensing_effects == _TRUE_) {
  
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt) {

      for (int l=2; l<=pbi->l_max; ++l)
        lsq_cl[l-2] = l*l*pbi->lensed_cls[index_lt][l-2];
  
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
                    pbi->lensed_d_lsq_cls[index_lt],
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);

    } // end of loop on index_lt
  } // end of if 
  
  free (l_array);
  free(lsq_cl);
  free(dd_lsq_cl);

  /* Free memory */
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {    
    if (psp->md_size > 1) free(cl_md[index_mode]);  
    if (psp->ic_size[index_mode] > 1) free(cl_md_ic[index_mode]);
  }  
  free(cl_md_ic);
  free(cl_md);
  free(cl);
  
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
    struct lensing * ple,
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
                  ple,
                  pbi),
      pbi->error_message,
      pbi->error_message);
      
  }
  
  /* Apart from the secondary bispectra in pbi->bispectra, all the arrays needed by the subsequent modules
  have been filled. If the user requested to load the bispectra from disk, we can stop the execution of
  this module now without regrets. */
  if (ppr->load_bispectra_from_disk == _TRUE_) {
    
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
  
  if (ppr->store_bispectra_to_disk == _TRUE_) {

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
    
          #pragma omp parallel for
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
  pwb->r_min = ppr->r_min;
  pwb->r_max = ppr->r_max;
  pwb->r_size = ppr->r_size;
    
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

    if (pbi->bispectra_verbose > 1)
        printf("     * integrating the filter functions over r ...\n");

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
 
  free (pwb);
 
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
    struct lensing * ple,
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
  // =                           Define analytic bispectra                             =
  // ===================================================================================

  /* Associate to each analytic bispectrum a function. In order to add your custom
  function, just add a line here, */
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    pbi->bispectrum_function[index_bt] = NULL;

    if ((pbi->has_cmb_lensing == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_bispectrum;

    else if ((pbi->has_cmb_lensing_squeezed == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_bispectrum;

    else if ((pbi->has_cmb_lensing_kernel == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing_kernel))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_kernel;

    else if ((pbi->has_local_squeezed == _TRUE_) && (index_bt == pbi->index_bt_local_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_local_squeezed_bispectrum;
    
    else if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_intrinsic_squeezed_bispectrum;
     
    else if ((pbi->has_intrinsic_squeezed_unlensed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed_unlensed))
      pbi->bispectrum_function[index_bt] = bispectra_intrinsic_squeezed_unlensed_bispectrum;
     
    else if ((pbi->has_quadratic_correction == _TRUE_) && (index_bt == pbi->index_bt_quadratic))
      pbi->bispectrum_function[index_bt] = bispectra_quadratic_bispectrum;
    
    else if ((pbi->has_cosine_shape == _TRUE_) && (index_bt == pbi->index_bt_cosine))
      pbi->bispectrum_function[index_bt] = bispectra_cosine_bispectrum;

  }
  

  // ===================================================================================
  // =                                   Main loop                                     =
  // ===================================================================================
    
  /* We parallelize the outer loop over 'l1'. */
  #pragma omp parallel for shared (ppt,pbs,ptr,ppm,pbi,abort) private (thread)
  for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    int l1 = pbi->l[index_l1];

    for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {

      int l2 = pbi->l[index_l2];

      /* Skip those configurations that are forbidden by the triangular condition (optimization) */
      if (l2 < l1/2)
        continue;

      /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);

      /* Uncomment test the accuracy of threej_ratio_M_recursive */        
      // int M=4;
      // double threej_num[2*pbi->l_max+1], threej_den[2*pbi->l_max+1];
      // int l3_min_num, l3_min_den;
      // double min_D, max_D;
      // class_call_parallel (drc3jj (
      //                        MAX(l1,M), MAX(l2,M), 0, -M,
      //                        &min_D, &max_D,
      //                        &(threej_num[0]),
      //                        (2*pbi->l_max+1),
      //                        pbi->error_message),
      //   pbi->error_message,
      //   pbi->error_message);
      // l3_min_num = (int)(min_D + _EPS_);          
      // class_call_parallel (drc3jj (
      //                        MAX(l1,M), MAX(l2,M), 0, 0,
      //                        &min_D, &max_D,
      //                        &(threej_den[0]),
      //                        (2*pbi->l_max+1),
      //                        pbi->error_message),
      //   pbi->error_message,
      //   pbi->error_message);
      // l3_min_den = (int)(min_D + _EPS_);          
      // 
      // for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
      //   int l3 = pbi->l[index_l3];
      //   if ((l1+l2+l3)%2==0) {
      //     if ((l1<M) || (l2<M) || (l3<M)) continue;
      //     double * ratio = malloc(sizeof(double)*(M+1));
      //     class_call_parallel (threej_ratio_M_recursive(l1, l2, l3, M, ratio, pbi->error_message),
      //       pbi->error_message, pbi->error_message);
      //     double res_1 = threej_num[l3-l3_min_num]/threej_den[l3-l3_min_den];
      //     double res_2 = ratio[M];
      //     double frac = 1-res_1/res_2;
      //     class_test_parallel (fabs(frac) > _SMALL_,
      //       pbi->error_message,
      //       "(%3d,%3d,%3d,M=%d), res_1=%14.6g, res_2=%14.6g, diff=%14.6g\n",
      //       l1, l2, l3, M, res_1, res_2, frac);
      //   }
      // }
          

      // -----------------------------------------------------------------------------------
      // -                               Compute bispectra                                 -
      // -----------------------------------------------------------------------------------

      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

        /* Index of the current (l1,l2,l3) configuration */
        int l3 = pbi->l[index_l3];
        long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
        
        /* Compute 3J ratios, needed only for polarised bispectra such as CMB-lensing
        and the quadratic correction. */
        double threej_ratio_20m2, threej_ratio_m220, threej_ratio_0m22;

        if (pbi->need_3j_symbols == _TRUE_) {          

          class_call_parallel (threej_ratio_M (l2, l1, l3, 2, &threej_ratio_20m2, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l3, l1, l2, 2, &threej_ratio_m220, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l1, l2, l3, 2, &threej_ratio_0m22, pbi->error_message),
            pbi->error_message, pbi->error_message);

        } // end of 3j computation

        for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

          if (pbi->bispectrum_type[index_bt] != analytical_bispectrum)
            continue;
      
          /* Check that the current bispectrum has a function associated to it */
          class_test_parallel (pbi->bispectrum_function[index_bt]==NULL,
            pbi->error_message,
            "no function associated for the bispectrum '%s'. Maybe it's not analytical?",
            pbi->bt_labels[index_bt]);

          for (int X = 0; X < pbi->bf_size; ++X) {
            for (int Y = 0; Y < pbi->bf_size; ++Y) {
              for (int Z = 0; Z < pbi->bf_size; ++Z) {

                  /* Compute the bispectrum using the function associated to index_bt */
                  class_call_parallel ((*pbi->bispectrum_function[index_bt]) (
                                ppr, psp, ple, pbi,
                                l1, l2, l3,
                                X, Y, Z,
                                threej_ratio_20m2,
                                threej_ratio_m220,
                                threej_ratio_0m22,
                                &(pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3])),
                    pbi->error_message,
                    pbi->error_message);
      
                  /* Update the counter */
                  increase_counter:;
                  #pragma omp atomic
                  pbi->count_memorised_for_bispectra++;

              } // end of for(Z)
            } // end of for(Y)
          } // end of for(X)
        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_bt)
    #pragma omp flush(abort)
  } // end of for(index_l1)
  if (abort == _TRUE_) return _FAILURE_;  // end of parallel region

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
  pwb->r_min = ppr->r_min;
  pwb->r_max = ppr->r_max;
  pwb->r_size = ppr->r_size;
    

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
  class_open (pbi->bispectra_files[index_bt], pbi->bispectra_paths[index_bt], "a+b", pbi->error_message);

  /* Print some debug */
  if (pbi->bispectra_verbose > 2)
    printf("     * writing bispectra to disk for index_bt=%d on '%s'\n",
      index_bt, pbi->bispectra_paths[index_bt]);

  /* Write all the independent (l1,l2,l3) triplets for this bispectrum */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        fwrite(
              pbi->bispectra[index_bt][X][Y][Z],
              sizeof(double),
              pbi->n_independent_configurations,
              pbi->bispectra_files[index_bt]
              );

  /* Close file */
  fclose(pbi->bispectra_files[index_bt]);
  
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
  class_open (pbi->bispectra_files[index_bt], pbi->bispectra_paths[index_bt], "rb", pbi->error_message);

  /* Print some debug */
  if (pbi->bispectra_verbose > 2)
    printf("     * reading bispectra from disk for index_bt=%d on'%s'\n", index_bt, pbi->bispectra_paths[index_bt]);

  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {

        int n_to_read = pbi->n_independent_configurations;
  
        /* Read a chunk with all the independent (l1,l2,l3) triplets to pbi->bispectra[index_bt] */
        int n_read = fread(
                pbi->bispectra[index_bt][X][Y][Z],
                sizeof(double),
                n_to_read,
                pbi->bispectra_files[index_bt]);

        class_test(n_read != n_to_read,
          pbi->error_message,
          "Could not read in '%s' file, read %d entries but expected %d",
            pbi->bispectra_paths[index_bt], n_read, n_to_read);        
      }
    }
  }
  
  /* Close file */
  fclose(pbi->bispectra_files[index_bt]); 

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


/*
 * Implement the formula for the CMB lensing bispectrum including polarisation, in Eq. 4.5 of
 * Lewis, Challinor & Hanson 2011 (http://uk.arxiv.org/abs/1101.2234). With respect to that
 * formula, in SONG we have i->X1, j->X2, k->X3, and l1<->l3. We have not included the imaginary
 * term in square brackets which is only needed when including B-mode polarisation (not
 * supported yet, TODO).
 *
 * This formula is non-perturbative (in the sense of Sec. 3.2, ibidem) and is valid for all
 * (l1,l2,l3) configurations, including non-squeezed ones.
 */

int bispectra_cmb_lensing_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------
  
  /* When only including temperature, the CMB-lensing bispectrum reduces to
    b^TTT_l1l2l3 = 1/2 * [l2(l2+1) + l3(l3+1) - l1(l1+1)] C_l2^{T\psi} C_l3^{TT} + 5 permutations,
  where \psi is the lensing potential. */
  
  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    /* Temperature-lensing potential correlation */
    double C_l1_Tp = pbi->cls[psp->index_ct_tp][l1-2];
    double C_l2_Tp = pbi->cls[psp->index_ct_tp][l2-2];
    double C_l3_Tp = pbi->cls[psp->index_ct_tp][l3-2];
                        
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
    double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
                  
    /* Use lensed temperature C_l's if available */
    if (pbi->include_lensing_effects == _TRUE_) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
      C_l3 = pbi->lensed_cls[ple->index_lt_tt][l3-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l2*(l2+1) + l3*(l3+1) - l1*(l1+1) ) * C_l2_Tp * C_l3
      + ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * C_l3_Tp * C_l2
      + ( l1*(l1+1) + l3*(l3+1) - l2*(l2+1) ) * C_l1_Tp * C_l3
      + ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * C_l3_Tp * C_l1
      + ( l1*(l1+1) + l2*(l2+1) - l3*(l3+1) ) * C_l1_Tp * C_l2
      + ( l2*(l2+1) + l1*(l1+1) - l3*(l3+1) ) * C_l2_Tp * C_l1
    );
  
    /* Debug - print temperature-lensing potential C_l's */
    // if ((l1==pbi->l_max) && (l2==pbi->l_max))
    //   fprintf (stderr, "%4d %10g %10g %10g\n", l3, C_l3_TT, C_l3, C_l3_Tp);
  
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // end of temperature only
  
  // -------------------------------------------------------------------------------
  // -                         Determine field coefficients                        -
  // -------------------------------------------------------------------------------
  
  /* Set the amplitude of the geometrical factor in Eq. 4.4 of Lewis et al. 2011. This is either
  2 (if the two 3j-symbols add up) or 0 (if they cancel). Whether they add up or cancel depends
  on the parity of the considered field (even for T and E and odd for B). */
  
  double S_X1=0, S_X2=0, S_X3=0;
  int L = l3-l1-l2;
  
  if (pbi->has_bispectra_t == _TRUE_) {
    if (X1 == pbi->index_bf_t) S_X1 = 2;
    if (X2 == pbi->index_bf_t) S_X2 = 2;
    if (X3 == pbi->index_bf_t) S_X3 = 2;
  }
  /* TODO: rayleigh-phi correlation not implemented yet */
  // if (pbi->has_bispectra_r == _TRUE_) {
  //   if (X1 == pbi->index_bf_r) S_X1 = 2;
  //   if (X2 == pbi->index_bf_r) S_X2 = 2;
  //   if (X3 == pbi->index_bf_r) S_X3 = 2;
  // }
  if (pbi->has_bispectra_e == _TRUE_) {
    if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
    if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
    if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
  }
  if (pbi->has_bispectra_b == _TRUE_) { /* Note that <TB> vanishes, hence the negative values */
    if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
    if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
    if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
  }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "CMB-lensing bispectrum for B-modes not implemented yet.");
  

  // ----------------------------------------------------------------------------------
  // -                               Obtain the C_l's                                 -
  // ----------------------------------------------------------------------------------
  
  /* Get the C_l's involving the lensing potential \phi. When implementing the B-modes, remember
  to set them to zero. */
  double C_l1_X1_p = pbi->cls[pbi->index_ct_of_phi_bf[ X1 ]][l1-2];
  double C_l2_X2_p = pbi->cls[pbi->index_ct_of_phi_bf[ X2 ]][l2-2];
  double C_l3_X3_p = pbi->cls[pbi->index_ct_of_phi_bf[ X3 ]][l3-2];
  
  /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
  double C_l3_X1_X3 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X3 ]][l3-2];
  double C_l2_X1_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X2 ]][l2-2];
  double C_l3_X2_X3 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X3 ]][l3-2];
  double C_l1_X2_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X1 ]][l1-2];
  double C_l2_X3_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ X2 ]][l2-2];
  double C_l1_X3_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ X1 ]][l1-2];
    
  /* Use lensed temperature C_l's if available */
  if (pbi->include_lensing_effects == _TRUE_) {
    C_l3_X1_X3 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X3 ]][l3-2];
    C_l2_X1_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X2 ]][l2-2];
    C_l3_X2_X3 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X3 ]][l3-2];
    C_l1_X2_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X1 ]][l1-2];
    C_l2_X3_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X3 ][ X2 ]][l2-2];
    C_l1_X3_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X3 ][ X1 ]][l1-2];
  }
  
  // ----------------------------------------------------------------------------------
  // -                               Get the right 3j's                               -
  // ----------------------------------------------------------------------------------
    
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];
  int F_X3 = pbi->field_spin[X3];
  
  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l1_l2_l3_FX1_0_mFX1 = 1, threej_l1_l3_l2_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    // class_call (threej_ratio_recursive (l2, l1, l3, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l1_l2_l3_FX1_0_mFX1 = ratio[M];
    // class_call (threej_ratio_recursive (l3, l1, l2, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l1_l3_l2_FX1_0_mFX1 = ratio[M];
    threej_l1_l2_l3_FX1_0_mFX1 = threej_ratio_20m2;
    threej_l1_l3_l2_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l2_l3_l1_FX2_0_mFX2 = 1, threej_l2_l1_l3_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    // class_call (threej_ratio_recursive (l3, l2, l1, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l2_l3_l1_FX2_0_mFX2 = ratio[M];
    // class_call (threej_ratio_recursive (l1, l2, l3, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l2_l1_l3_FX2_0_mFX2 = ratio[M];
    threej_l2_l3_l1_FX2_0_mFX2 = threej_ratio_m220;
    threej_l2_l1_l3_FX2_0_mFX2 = threej_ratio_0m22;
  }
    
  double threej_l3_l1_l2_FX3_0_mFX3 = 1, threej_l3_l2_l1_FX3_0_mFX3 = 1;
  if (F_X3==2) {
    // class_call (threej_ratio_recursive (l1, l3, l2, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l3_l1_l2_FX3_0_mFX3 = ratio[M];
    // class_call (threej_ratio_recursive (l2, l3, l1, M, ratio, pbi->error_message),
    //   pbi->error_message, pbi->error_message);
    // threej_l3_l2_l1_FX3_0_mFX3 = ratio[M];
    threej_l3_l1_l2_FX3_0_mFX3 = threej_ratio_0m22;
    threej_l3_l2_l1_FX3_0_mFX3 = threej_ratio_20m2;
  }

  /* Obtain the geometric factor F^+s_l1l2l3 */
  double F_l1_l2_l3_X1 = 0.25 * ( l2*(l2+1) + l3*(l3+1) - l1*(l1+1) ) * S_X1 * threej_l1_l2_l3_FX1_0_mFX1; /* 1-2-3 */
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l1_l3_X2 = 0.25 * ( l1*(l1+1) + l3*(l3+1) - l2*(l2+1) ) * S_X2 * threej_l2_l1_l3_FX2_0_mFX2; /* 2-1-3 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
  double F_l3_l1_l2_X3 = 0.25 * ( l1*(l1+1) + l2*(l2+1) - l3*(l3+1) ) * S_X3 * threej_l3_l1_l2_FX3_0_mFX3; /* 3-1-2 */
  double F_l3_l2_l1_X3 = 0.25 * ( l2*(l2+1) + l1*(l1+1) - l3*(l3+1) ) * S_X3 * threej_l3_l2_l1_FX3_0_mFX3; /* 3-2-1 */
  
  // ---------------------------------------------------------------------------------------
  // -                                  Bispectrum formula                                 -
  // ---------------------------------------------------------------------------------------
  
  /* CMB-lensing bispectrum formula, from Eq. 4.5 of Lewis, Challinor & Hanson 2011. */
  *result = 
      C_l2_X2_p * C_l3_X1_X3 * F_l1_l2_l3_X1   /* 1-2-3 */
    + C_l3_X3_p * C_l2_X1_X2 * F_l1_l3_l2_X1   /* 1-3-2 */
    + C_l1_X1_p * C_l3_X2_X3 * F_l2_l1_l3_X2   /* 2-1-3 */
    + C_l3_X3_p * C_l1_X2_X1 * F_l2_l3_l1_X2   /* 2-3-1 */
    + C_l1_X1_p * C_l2_X3_X2 * F_l3_l1_l2_X3   /* 3-1-2 */
    + C_l2_X2_p * C_l1_X3_X1 * F_l3_l2_l1_X3;  /* 3-2-1 */
                  
  /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "CMB-lensing bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
     l1, l2, l3, *result, exact, diff);
  }

  /* Debug - Print the 3j-symbols */
  // printf ("threej_l3_l2_l1_FX3_0_mFX3 = %g\n", threej_l3_l2_l1_FX3_0_mFX3);
  // printf ("threej_l3_l1_l2_FX3_0_mFX3 = %g\n", threej_l3_l1_l2_FX3_0_mFX3);
  // printf ("threej_l2_l1_l3_FX2_0_mFX2 = %g\n", threej_l2_l1_l3_FX2_0_mFX2);
  // printf ("threej_l2_l3_l1_FX2_0_mFX2 = %g\n", threej_l2_l3_l1_FX2_0_mFX2);
  // printf ("threej_l1_l3_l2_FX1_0_mFX1 = %g\n", threej_l1_l3_l2_FX1_0_mFX1);
  // printf ("threej_l1_l2_l3_FX1_0_mFX1 = %g\n", threej_l1_l2_l3_FX1_0_mFX1);
  // printf("\n");
  
  /* Debug - Print the C_l's */
  // printf ("C_l1_X1_p = %g\n", C_l1_X1_p);
  // printf ("C_l2_X2_p = %g\n", C_l2_X2_p);
  // printf ("C_l3_X3_p = %g\n", C_l3_X3_p);
  // printf ("\n");

  return _SUCCESS_;
  
}


/*
 * Compute the CMB lensing bispectrum kernel including polarisation in the squeezed limit, valid
 * when l3<<l1 and l3<<l2. Here we code just the kernel in Eq. 5.20 of Lewis, Challinor & Hanson
 * 2011 (http://uk.arxiv.org/abs/1101.2234); the full squeezed bispectrum is given by the
 * product between the kernel and C_l3^{X3\phi}.
 *
 * With respect to Eq. 5.20 (ibidem) in SONG we use i->X1, j->X2, k->X3 and then we perform
 * the substitution (X1,l1)<->(X3,l3), as in our convention l3 rather than l1 is the smallest
 * multipole. Therefore, this is what we code here:
 *
 * A^{X1,X2}_{l1l2l3} = \tilde{C}^{X2,X1}_{l1} * F^{X1}_{l1l3l2}
 *                    + \tilde{C}^{X1,X2}_{l2} * F^{X2}_{l2l3l1}
 *
 */

int bispectra_cmb_lensing_squeezed_kernel (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------
  
  /* When only including temperature, the CMB-lensing bispectrum reduces to
    b^TTT_l1l2l3 = 1/2 * [l2(l2+1) + l3(l3+1) - l1(l1+1)] C_l2^{T\psi} C_l3^{TT} + 5 permutations,
  where \psi is the lensing potential. */
  
  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
                  
    /* Use lensed temperature C_l's if available */
    if (pbi->include_lensing_effects == _TRUE_) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * C_l2
      + ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * C_l1
    );

    /* Uncomment to compute the bispectrum rather than the kernel */
    // ttt *= pbi->cls[psp->index_ct_tp][l3-2];                      
    
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // end of temperature only
  
  // -------------------------------------------------------------------------------
  // -                         Determine field coefficients                        -
  // -------------------------------------------------------------------------------
  
  /* Set the amplitude of the geometrical factor in Eq. 4.4 of Lewis et al. 2011. This is either
  2 (if the two 3j-symbols add up) or 0 (if they cancel). Whether they add up or cancel depends
  on the parity of the considered field (even for T and E and odd for B). */
  
  double S_X1=0, S_X2=0;
  int L = l3-l1-l2;
  
  if (pbi->has_bispectra_t == _TRUE_) {
    if (X1 == pbi->index_bf_t) S_X1 = 2;
    if (X2 == pbi->index_bf_t) S_X2 = 2;
  }
  /* TODO: rayleigh-phi correlation not implemented yet */
  // if (pbi->has_bispectra_r == _TRUE_) {
  //   if (X1 == pbi->index_bf_r) S_X1 = 2;
  //   if (X2 == pbi->index_bf_r) S_X2 = 2;
  // }
  if (pbi->has_bispectra_e == _TRUE_) {
    if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
    if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  }
  if (pbi->has_bispectra_b == _TRUE_) { /* Note that <TB> vanishes, hence the negative values */
    if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
    if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "CMB-lensing squeezed bispectrum for B-modes not implemented yet.");
  
  // ----------------------------------------------------------------------------------
  // -                               Obtain the C_l's                                 -
  // ----------------------------------------------------------------------------------
  
  /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
  double C_l2_X1_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X2 ]][l2-2];
  double C_l1_X2_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X1 ]][l1-2];
    
  /* Use lensed temperature C_l's if available */
  if (pbi->include_lensing_effects == _TRUE_) {
    C_l2_X1_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X2 ]][l2-2];
    C_l1_X2_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X1 ]][l1-2];
  }
  
  // ----------------------------------------------------------------------------------
  // -                               Get the right 3j's                               -
  // ----------------------------------------------------------------------------------
  
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];

  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l1_l2_l3_FX1_0_mFX1 = 1, threej_l1_l3_l2_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    threej_l1_l2_l3_FX1_0_mFX1 = threej_ratio_20m2;
    threej_l1_l3_l2_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l2_l3_l1_FX2_0_mFX2 = 1, threej_l2_l1_l3_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    threej_l2_l3_l1_FX2_0_mFX2 = threej_ratio_m220;
    threej_l2_l1_l3_FX2_0_mFX2 = threej_ratio_0m22;
  }

  /* Obtain the geometric factor F^+s_l1l2l3 */
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
      
  // ---------------------------------------------------------------------------------------
  // -                                  Bispectrum formula                                 -
  // ---------------------------------------------------------------------------------------
  
  /* Kernel of the CMB-lensing bispectrum in the squeezed limit, from Eq. 5.20 of Lewis,
  Challinor & Hanson 2011. This is simply the general formula with C_l1_X1_p=C_l2_X2_p=0 
  (see bispectra_cmb_lensing_bispectrum) */
  *result = C_l2_X1_X2 * F_l1_l3_l2_X1
          + C_l1_X2_X1 * F_l2_l3_l1_X2;
                  
  /* Uncomment to compute the bispectrum rather than the kernel */
  // *result *= pbi->cls[pbi->index_ct_of_phi_bf[ X3 ]][l3-2];
                  
  /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "CMB-lensing squeezed bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
     l1, l2, l3, *result, exact, diff);
  }

  return _SUCCESS_;
  
}

/*
 * Compute the CMB lensing bispectrum including polarisation in the squeezed limit, valid
 * when l3<<l1 and l3<<l2. This is given by the kernel computed in 'bispectra_cmb_lensing_squeezed_kernel'
 * times the cross-correlation with the lensing potential, C_l3^{X3\phi} (see Eq. 5.20 of Lewis,
 * Challinor & Hanson 2011 (http://uk.arxiv.org/abs/1101.2234).
 *
 * This formula is non-perturbative (in the sense of Sec. 3.2, ibidem) and is valid only
 * for squeezed configurations, where l3<<l1 and l3<<l2; it is obtained from the general
 * formula (Eq. 4.5, ibidem) by setting C_l1_X1_p=C_l2_X2_p=0. It is an excellent approximation
 * for the CMB-lensing bispectrum, as most of its signal is in these squeezed configurations.
 *
 * Look at the comment above 'bispectra_intrinsic_squeezed_bispectrum' for details
 * about the symmetry properties of this special squeezed bispectrum.
 * 
 */
int bispectra_cmb_lensing_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Compute the kernel */
  class_call (bispectra_cmb_lensing_squeezed_kernel (
                ppr, psp, ple, pbi,
                l1, l2, l3,
                X1, X2, X3,
                threej_ratio_20m2,
                threej_ratio_m220,
                threej_ratio_0m22,
                result),
    pbi->error_message,
    pbi->error_message);

  /* Obtain the bispectrum in the squeezed limit by multiplication with C_l3^{X3\phi} */
  *result *= pbi->cls[pbi->index_ct_of_phi_bf[ X3 ]][l3-2];

  return _SUCCESS_;

}


/** 
 * Squeezed approximation for the local bispectrum. This is the generalisation to polarisation
 * of the temperature approximation first described in Gangui et al. 1994, Komatsu & Spergel 2001.
 *
 */
int bispectra_local_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
  choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double cl1_XY = pbi->cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double cl2_XY = pbi->cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Use lensed temperature C_l's if available */
  if (pbi->include_lensing_effects == _TRUE_) {
    cl1_XY = pbi->lensed_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];    
    cl2_XY = pbi->lensed_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];    
  }

  /* Obvious generalisation of the TTT result in Eq. 2.16 of http://arxiv.org/abs/1201.1010 */
  *result = 6/5. * cl3_Zz * (cl1_XY + cl2_XY);

  /* TTT result */
  // *result = 6/5. * (
  //     pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
  //   + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
  //   + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2]
  // );      

  /* EEE result */
  // *result = 6/5. * (
  //     pbi->cls[psp->index_ct_ee][l1-2] * pbi->cls[psp->index_ct_ee][l2-2]
  //   + pbi->cls[psp->index_ct_ee][l2-2] * pbi->cls[psp->index_ct_ee][l3-2]
  //   + pbi->cls[psp->index_ct_ee][l3-2] * pbi->cls[psp->index_ct_ee][l1-2]
  // );      

  return _SUCCESS_;

}



/** 
 * Compute the squeezed-limit approximation for the intrinsic bispectrum, as reported in eq. 4.1 and 4.2 of
 * Lewis 2012 (see also Creminelli et al. 2004, Creminelli et al. 2011, Bartolo et al. 2012). This is the
 * general formula that includes polarisation. With respect to Lewis' formula, (i,l1)->(Z,l3), (j,l2)->(X,l1),
 * (k,l3)->(Y,l2) and \zeta -> z.
 * 
 * ~~~ CONSIDERATIONS THAT APPLY TO ALL "SQUEEZED" BISPECTRA ~~~
 *
 * In this and in the other "squeezed" approximation functions, l3 is taken to be the long-wavelength
 * mode, l1 and l2 the short ones. This choice is preferred because in SONG we loop over (l1,l2,l3)
 * configurations that satisfy the condition l1>=l2>=l3.
 *
 * Therefore, it is crucial that that, in the formulas below, the multipole l3 is associated with the
 * C_l that correlates with the comoving curvature perturbation zeta (for the CMB lensing bispectrum
 * it would be the lensing potential phi) because l3 is the smallest multipole, which describes the
 * long wavelength mode. Also, l3 must be associated with Z (or X3) because our convention
 * for the bispectrum is <X_l1 Y_l2 Z_l3>. If you give l3 to another field, then the Fisher matrix
 * estimator will associate to that field the wrong covariance matrix, and the result will change
 * drastically.
 *
 * Because of the special role played by l3, the bispectrum computed in this and the other 'squeezed'
 * approximation functions is NOT symmetric with respect to an exchange of (l1,X) <-> (l3,Z) or of
 * or (l2,Y) <-> (l3,Z), contrary to the other bispectra, which are computed as <X_l1 Y_l2 Z_l3>.
 * Therefore, for this bispectrum one cannot obtain the configurations outside l1>=l2>=l3 by
 * permuting the XYZ indices, as it is done, for example, in the 'print_bispectra' function.
 * 
 */

int bispectra_intrinsic_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* Uncomment to exclude the polarisation when it is assigned the large-scale mode */
  // if ((pbi->has_bispectra_e==_TRUE_) && (Z==pbi->index_bf_e)) {
  //   *result = 0;
  //   return _SUCCESS_;
  // }

  /* Uncomment to restrict the approximation to squeezed configurations, as in Sec. 6 of Creminelli,
  Pitrou & Vernizzi 2011. Note that in our case the smallest mode is l3, while in that paper it
  is l1. IMPORTANT: setting a bispectrum to zero might screw up some matrix inversions done in the
  Fisher module, especially when 'pfi->include_lensing_effects' is _TRUE_. */
  // if (!((l3<=100) && (l2>=10*l3))) {
  //   *result = 0;
  //   return _SUCCESS_;
  // }

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
  choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Use lensed temperature C_l's if available */
  if ((pbi->include_lensing_effects == _TRUE_) && (pbi->lensed_intrinsic == _TRUE_)) {
    dcl1_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l1-2];    
    dcl2_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l2-2];    
  }
          
  /* Ricci focussing in Lewis 2012 (eq. 4.1) */
  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);

  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
  double bolometric_T_lewis_redshift = 0;
          
  if (pbi->has_bispectra_t == _TRUE_) {
    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf_bf[Z][pbi->index_bf_t]][l3-2];
    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf_bf[X][pbi->index_bf_t]][l1-2];
    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    if ((pbi->include_lensing_effects == _TRUE_) && (pbi->lensed_intrinsic == _TRUE_)) {
      if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[X][pbi->index_bf_t]][l1-2];
      if (X == pbi->index_bf_t) cl2_Yt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    }
    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
  }
          
  class_test ((pbi->has_bispectra_e == _TRUE_) &&
    ((X == pbi->index_bf_e) && (Y == pbi->index_bf_e))
    && (bolometric_T_lewis_redshift!=0),
    pbi->error_message,
    "anisotropic redshifting should vanish (eq. 4.2 of Lewis 2012)", bolometric_T_lewis_redshift);
          
  /* Sum of Ricci focussing and redshift modulation */
  *result = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

  return _SUCCESS_;

}


/**
 * Same as above, but using the unlensed power spectrum, which gives more signal as the
 * lensed C_l's oscillate less. See 'bispectra_intrinsic_squeezed_bispectrum' for details.
 */

int bispectra_intrinsic_squeezed_unlensed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
  choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Ricci focussing in Lewis 2012 (eq. 4.1) */
  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);

  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
  double bolometric_T_lewis_redshift = 0;
          
  if (pbi->has_bispectra_t == _TRUE_) {
    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf_bf[Z][pbi->index_bf_t]][l3-2];
    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf_bf[X][pbi->index_bf_t]][l1-2];
    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
  }
          
  /* Sum of Ricci focussing and redshift modulation */
  *result = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

  return _SUCCESS_;

}


/** 
 * 
 * The CMB temperature is obtained from the second-order distribution function by adding to it
 * a quadratic term (see Eq. 3.12 of arXiv:1401.3296).
 * At the bispectrum level this extra term translates to a four-point function in the first-order
 * perturbations, which in turn is expressed in terms of products of C_l's.
 * The same correction term appears for the delta_tilde transformation (see Eq. 3.5 of arXiv:1401.3296).
 * The correction in both cases has this form (see Eq. 3.6, 3.7 and 3.9 of same paper):
 * 
 *   QC_{l1 l2 l3} = S_X3 * i^(L+V_X3) * ( l'  l''  |    l3 ) * (  l'  l'' | l3 )  
 *                                       ( 0   F_X3 | -F_X3 )   (  m'  m'' | m3  ) 
 *                   * <a^X1_l1m1 * a^X2_l2m2 * a^I_l'm' * a^T_X3_l''m''> + 2 permutations (1->2->3)
 *   
 * where L=l3-l'-l'', all the a_lm's are first-order and, in this loop, X1=X, X2=Y and X3=Z. The
 * permutations go over 1->2->3 and refer to the position of the second-order perturbation in
 * the bispectrum; in the formula above, the positioning is <a^(1)X1_l1m1 * a^(1)X2_l2m2 * a^(2)X3_l3m3>.
 * The coefficients take different values according to which field is X3:
 * 
 *   X3=I -> F=0, S=1, V_X3=0,  T_X3=I
 *   X3=E -> F=2, S=2, V_X3=0,  T_X3=E
 *   X3=B -> F=2, S=2, V_X3=-1, T_X3=E
 * 
 * For X3=E-modes, the sum over l' and l'' only includes EVEN values of l3-l'-l'', for X3=B-modes
 * only includes ODD values of l3-l'-l''.
 * 
 * By employing Wick's theorem, the above can be expressed in terms of the angular power spectrum
 * of the CMB, C_l = <a_lm a^*_lm> (note that the a_lm's already include the 1/4 factor):
 * 
 *   QC_{l1 l2 l3} = 4 * B_{l1 l2 l3} * i^V_X3 * G^m1m2m3_l1l2l3 * S_X3 * [ 
 *                                      ( l1    l2     l3 ) * C^{X1,I}_l1 * C^{X2,X3}_l2
 *                                      ( 0   F_X3  -F_X3 )
 *                                    + (   l1  l2     l3 ) * C^{X1,X3}_l1 * C^{X2,I}_l2
 *                                      ( F_X3   0  -F_X3 )
 *                                                             ] + 2 permutations (1->2->3) ,
 * 
 * where
 * 
 *   B_{l1 l2 l3} = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*_PI_)) * ( l1    l2   l3 )
 *                                                              ( m1    m2   m3 ) .
 * 
 * Below we compute the quadratic correction term using this general formula and store it in
 * pbi->bispectra[index_bt_quadratic].
 * 
 * For intensity (X1=X2=X3=I), the formula reduces to 8 * [ C_l1*C_l2 + C_l1*C_l2 + C_l1*C_l2 ], while             
 * any contribution form the second-order B-modes has the V_X=-1 and is therefore purely imaginary.
 * Furthermore, the B-mode contribution only includes one term, theone where a^B_lm is second order,
 * as we ignore the first-order B-modes. In any case, we do not compute the quadratic term for a B-mode
 * spectrum.
 */
int bispectra_quadratic_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------

  /* The quadratic (reduced) bispectrum for TTT is very simple: 
    b^TTT_l1l2l3 = 8 * (C_l1*C_l2 + 2 perms) */

  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    double C_l1_TT = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2_TT = pbi->cls[psp->index_ct_tt][l2-2];
    double C_l3_TT = pbi->cls[psp->index_ct_tt][l3-2];
  
    ttt = 8 * (C_l1_TT*C_l2_TT + C_l1_TT*C_l3_TT + C_l2_TT*C_l3_TT);
    
    /* If only temperature is requested, skip what follows exit from the 'quadratic' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      goto end_quadratic; 
    }
  }

  // ---------------------------------------------------------------
  // -                  Determine field coefficients               -
  // ---------------------------------------------------------------
                
  /* Determine various indices and coefficients:
  1) The amplitude S_X of the term in QC that involves the second-order part of X. When doing so,
     take into account that the contribution from E^(2) only exists when l3-l1-l2 is even,
     while that from B^(2) only when l3-l1-l2 is odd.
  2) The field that goes in the C_l's, T_X, which is computed as T_I=I, T_E=E, T_B=E.
  3) The indices of the cross power spectra between X1, X2, X3 and I. These have
     to be set by hand, rather than using the array pbi->index_ct_of_bf_bf because pbi->index_ct_of_bf_bf
     only contains information on the fields appearing in one of the requested bispectrum. For example,
     if you only request EEE, then pbi->index_ct_of_bf_bf does not contain information about <ET>,
     because no bispectrum containing T is requested */
  
  double S_X1=0, S_X2=0, S_X3=0;
  int T_X1=0, T_X2=0, T_X3=0;
  int index_ct_X1_I, index_ct_X2_I, index_ct_X3_I;
  int L = l3-l1-l2;
  
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
    if (X1 == pbi->index_bf_e) {S_X1 = (L%2==0)?2:0; T_X1 = pbi->index_bf_e; index_ct_X1_I = psp->index_ct_te;}
    if (X2 == pbi->index_bf_e) {S_X2 = (L%2==0)?2:0; T_X2 = pbi->index_bf_e; index_ct_X2_I = psp->index_ct_te;}
    if (X3 == pbi->index_bf_e) {S_X3 = (L%2==0)?2:0; T_X3 = pbi->index_bf_e; index_ct_X3_I = psp->index_ct_te;}
  }
  if (pbi->has_bispectra_b == _TRUE_) { /* Note that <TB> vanishes, hence the negative values */
    if (X1 == pbi->index_bf_b) {S_X1 = (L%2!=0)?2:0; T_X1 = pbi->index_bf_e; index_ct_X1_I = -1;}
    if (X2 == pbi->index_bf_b) {S_X2 = (L%2!=0)?2:0; T_X2 = pbi->index_bf_e; index_ct_X2_I = -1;}
    if (X3 == pbi->index_bf_b) {S_X3 = (L%2!=0)?2:0; T_X3 = pbi->index_bf_e; index_ct_X3_I = -1;}
  }

  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "quadratic correction for B-mode bispectrum not implemented yet.");
  
  // ----------------------------------------------------------------------------------
  // -                               Obtain the C_l's                                 -
  // ----------------------------------------------------------------------------------
  
  /* Get the C_l's. When implementing the B-modes, remember to 
  set the C_l's to zero, so that the only contribution comes from the
  bispectrum with the second-order a^B_lm.  */
  /* TODO: do we need the lensed C_l's? */
  /* TODO: CHECK FACTORS 4!!!! */
    
  double C_l1_X1_I   = pbi->cls[index_ct_X1_I][l1-2];
  double C_l1_X1_TX2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ T_X2 ]][l1-2];
  double C_l1_X1_TX3 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ T_X3 ]][l1-2];
              
  double C_l2_X2_I   = pbi->cls[index_ct_X2_I][l2-2];
  double C_l2_X2_TX1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ T_X1 ]][l2-2];
  double C_l2_X2_TX3 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ T_X3 ]][l2-2];
  
  double C_l3_X3_I   = pbi->cls[index_ct_X3_I][l3-2];                  
  double C_l3_X3_TX1 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ T_X1 ]][l3-2];
  double C_l3_X3_TX2 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ T_X2 ]][l3-2];
    
  
  // ----------------------------------------------------------------------------------
  // -                               Get the right 3j's                               -
  // ----------------------------------------------------------------------------------
  
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];
  int F_X3 = pbi->field_spin[X3];
  
  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l2_l3_l1_0_FX1_mFX1 = 1, threej_l2_l3_l1_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    threej_l2_l3_l1_0_FX1_mFX1 = threej_ratio_20m2;
    threej_l2_l3_l1_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l3_l1_l2_0_FX2_mFX2 = 1, threej_l3_l1_l2_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    threej_l3_l1_l2_0_FX2_mFX2 = threej_ratio_m220;
    threej_l3_l1_l2_FX2_0_mFX2 = threej_ratio_0m22;
  }
  
  double threej_l1_l2_l3_0_FX3_mFX3 = 1, threej_l1_l2_l3_FX3_0_mFX3 = 1;
  if (F_X3==2) {
    threej_l1_l2_l3_0_FX3_mFX3 = threej_ratio_0m22;
    threej_l1_l2_l3_FX3_0_mFX3 = threej_ratio_20m2;
  }

  // ---------------------------------------------------------------------------------------
  // -                                  Bispectrum formula                                 -
  // ---------------------------------------------------------------------------------------

  /* The sum includes three terms, corresponding to the three possible types
  of second-order perturbations. The first term, involving T_X1, for example,
  refers to the term where a^X1_lm is second order. */
  *result = 4 * S_X1 * (threej_l2_l3_l1_0_FX1_mFX1*C_l2_X2_I*C_l3_X3_TX1 + threej_l2_l3_l1_FX1_0_mFX1*C_l2_X2_TX1*C_l3_X3_I)
          + 4 * S_X2 * (threej_l3_l1_l2_0_FX2_mFX2*C_l3_X3_I*C_l1_X1_TX2 + threej_l3_l1_l2_FX2_0_mFX2*C_l3_X3_TX2*C_l1_X1_I)
          + 4 * S_X3 * (threej_l1_l2_l3_0_FX3_mFX3*C_l1_X1_I*C_l2_X2_TX3 + threej_l1_l2_l3_FX3_0_mFX3*C_l1_X1_TX3*C_l2_X2_I);

  /* Check that for <TTT> the correction is equal to 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "quadratic bispectrum correction does not simplify to 8 (C_l1*C_l2 + perm); b=%g, exact=%g, diff=%g",
     *result, exact, diff);
  }

  /* Debug. Print the quadratic contribution */
  // printf ("B_quadratic_%3s[%3d,%3d,%3d] = %g\n",
  //   pbi->bfff_labels[X][Y][Z], l1, l2, l3, *result);

  end_quadratic: ;
  
  return _SUCCESS_;

}


/** 
 * Simple bispectrum used for testing purposes:
 *  b_l1l2l3 = 6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) * cos((l1+l2+l3)/50).
 * It is only defined for temperature.
 */
int bispectra_cosine_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{                   
  
  /* Get the squeezed approximations for the local bispectrum */
  class_call (bispectra_local_squeezed_bispectrum (
                ppr, psp, ple, pbi,
                l1, l2, l3,
                X1, X2, X3,
                threej_ratio_20m2,
                threej_ratio_m220,
                threej_ratio_0m22,
                result),
    pbi->error_message,
    pbi->error_message);

  *result *= cos((l1+l2+l3)/50.);

  return _SUCCESS_;
  
}




/** 
 * Window function for the local (l1,l2,l3) bispectrum.
 *
 * The natural scaling for the bispectrum (both the templates and the second-order one) is given by a product 
 * of two power spectra. When available, we always use the C_l's for the temperature as, when multiplied by l^2,
 * they are approximately the same order of magnitude for all l's, contrary to those for polarisation which are
 * very small for l<200.
 *
 */
int bispectra_local_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Interpolate all bispectra with temperature's C_l1*C_l2 */
  if (pbi->has_bispectra_t == _TRUE_) {
    *result = 
        pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
      + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
      + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  }
  else if (pbi->bf_size == 1) {
    *result =
        pbi->cls[0][l1-2] * pbi->cls[0][l2-2]
      + pbi->cls[0][l2-2] * pbi->cls[0][l3-2]
      + pbi->cls[0][l3-2] * pbi->cls[0][l1-2];
  }
  else {
    class_stop (pbi->error_message, "bispectra with two non-temperature fields are not implemented yet.\n");
  }

  /* Interpolate only temperature */
  // if ((pbi->has_bispectra_t == _TRUE_) && (X==pbi->index_bf_t) && (X==Y) && (X==Z)) {
  //   *result = 
  //       pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
  //     + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
  //     + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  // }
  // else {
  //   *result = 1;
  // }

  return _SUCCESS_;

}



/** 
 * Window function for the intrinsic (l1,l2,l3) bispectrum.
 *
 */
int bispectra_intrinsic_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  if (pbi->has_bispectra_t == _TRUE_) {
    *result = 
        pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
      + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
      + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  }
  else if (pbi->bf_size == 1) {
    *result =
        pbi->cls[0][l1-2] * pbi->cls[0][l2-2]
      + pbi->cls[0][l2-2] * pbi->cls[0][l3-2]
      + pbi->cls[0][l3-2] * pbi->cls[0][l1-2];
  }
  else {
    class_stop (pbi->error_message, "bispectra with two non-temperature fields are not implemented yet.\n");
  }

  return _SUCCESS_;

}



