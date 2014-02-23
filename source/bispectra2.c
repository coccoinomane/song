/** @file bispectra.c documented spectra module for second-order perturbations
 *
 * Guido W Pettinari, 19.07.2012
 */

#include "bispectra2.h"


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

int bispectra2_init (
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
     struct lensing * ple,
     struct bispectra * pbi
     )
{

  /* Check whether we need to compute intrinsic spectra at all */  
  if ((ppt2->has_bispectra == _FALSE_)       /* Have we computed the needed 2nd-order perturbations? */
     || (pbi->has_bispectra == _FALSE_)      /* Did we load the 1st-order bispectrum module? */
     || (pbi->n[intrinsic_bispectrum]<1)) {  /* Did we request at least a second-order bispectrum? */

    if (pbi->bispectra_verbose > 0)
      printf("No intrinsic bispectra requested. Second-order bispectra module skipped.\n");

    return _SUCCESS_;
  }
  
  

  // =======================================================
  // =                  Compute bispectra                  =
  // =======================================================
  
  class_call (bispectra2_harmonic (ppr,ppr2,pba,ppt,ppt2,pbs,pbs2,ptr,ptr2,ppm,psp,pbi),
    pbi->error_message,
    pbi->error_message);


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

int bispectra2_harmonic (
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
    )
{


  /* If the user requested to load the bispectra from disk, we can stop the execution of this function here */
  if (pbi->load_bispectra_from_disk == _TRUE_) {
    return _SUCCESS_;
  }


  // ================================================================================
  // =                          Compute intrinsic bispectra                         =
  // ================================================================================
  
  if (pbi->n[intrinsic_bispectrum] > 0) {
  
    struct bispectra_workspace_intrinsic * pwb;
    class_alloc (pwb, sizeof(struct bispectra_workspace_intrinsic), pbi->error_message);

    /* Allocate arrays inside the integration workspace */
    class_call (bispectra2_intrinsic_workspace_init(
                  ppr,
                  ppr2,
                  ppt,
                  ppt2,
                  pbs,
                  pbs2,
                  ptr,
                  ptr2,
                  ppm,
                  pbi,
                  pwb),
      pbi->error_message,
      pbi->error_message);

    /* Perform the bispectrum integration */
    class_call (bispectra2_intrinsic_init(
                  ppr,
                  ppr2,
                  ppt,
                  ppt2,
                  pbs,
                  pbs2,
                  ptr,
                  ptr2,
                  ppm,
                  pbi,
                  pwb),
      pbi->error_message,
      pbi->error_message);
  
    /* Free the 'pwb' workspace */
    class_call (bispectra2_intrinsic_workspace_free(
                  ppt2,
                  ptr2,
                  pbi,
                  pwb),
      pbi->error_message,
      pbi->error_message);
  
  }
  
  
  
  
  
    
  // ===================================================================================
  // =                             Add temperature correction                          =
  // ===================================================================================

  /* In the following loop, we add the correction that turns the brightness bispectrum, which we
  computed above, into the bolometric temperature bispectrum (see eq. 5.8 of http://arxiv.org/abs/1302.0832).
  The loop will run only for the intrinsic bispectra, as for the primary ones there is no difference between
  brightness and bolometric temperature. 
  
  Comment the loop if you want to use the brightness bispectrum rather than the bolometric temperature
  one. */

  /* TODO: update this for polarization, both E and B-modes. The best way to do this is to define
  another non-Fisher bispectrum type index_bt_temp_correction */

  if ((pbi->bf_size == 1) && (ppt->has_cl_cmb_temperature == _TRUE_)) {

    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {
  
      /* Skip the bispectrum if it not of the non-separable type */
      if (pbi->bispectrum_type[index_bt] != intrinsic_bispectrum)
        continue;
  
      for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
        
        double C_l1 = pbi->cls[psp->index_ct_tt][pbi->l[index_l1]-2];
        
        for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
     
          double C_l2 = pbi->cls[psp->index_ct_tt][pbi->l[index_l2]-2];
     
          /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
          int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
          int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
     
          for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
      
            double C_l3 = pbi->cls[psp->index_ct_tt][pbi->l[index_l3]-2];
  
            /* Index of the current (l1,l2,l3) configuration */
            long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  
            /* By default, we assume the standard correction to convert brightness temperature
            in bolometric temperature */              
            double bolometric_correction = - 3 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
          
            /* We assume that when the quadratic sources are not included, the user is trying to run the
            code as a first-order one, and hence we turn off the bolometric temperature correction. */
            if (ppt2->has_quadratic_sources == _FALSE_)
              bolometric_correction = 0;
  
            /* If the transfer functions are for delta_tilde, and not delta, we need a different correction to obtain the bolometric
            temperature (see Huang & Vernizzi 2012) */
            else if ((ppt2->use_delta_tilde_in_los == _TRUE_) || (ppt2->use_delta_tilde_approx_in_los == _TRUE_))
              bolometric_correction = C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3;
  
            /* Add the correction */
            pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3] += bolometric_correction;
  
            /* Update the counter */
            #pragma omp atomic
            pbi->count_memorised_for_bispectra++;
  
            /* Check against nan's */
            if (isnan(pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3]))
              printf ("\n### WARNING ###\nb(%d,%d,%d) = %g!!!\n\n",
              pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3]);
  
          } // end of for(index_l3)
        } // end of for(index_l2)
      } // end of for(index_l1)
    } // end of for(index_bt)

    /* Check that we correctly filled the bispectra array */
    class_test_permissive (pbi->count_allocated_for_bispectra != pbi->count_memorised_for_bispectra,
      pbi->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!",
      pbi->count_allocated_for_bispectra, pbi->count_memorised_for_bispectra);
  
    /* Print information on memory usage */
    if (pbi->bispectra_verbose > 1)
      printf(" -> memorised ~ %.3g MB (%ld doubles) in the bispectra array\n",
        pbi->count_memorised_for_bispectra*sizeof(double)/1e6, pbi->count_memorised_for_bispectra);

  } // end of if temperature only 
  
  
  // ============================================================================================
  // =                                   Save bispectra to disk                                 =
  // ============================================================================================
  
  if (pbi->store_bispectra_to_disk == _TRUE_)
    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt)
      if (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum)
        class_call (bispectra_save_to_disk (
                      pbi,
                      index_bt),
          pbi->error_message,
          pbi->error_message);
  
  
  return _SUCCESS_;

}






int bispectra2_intrinsic_init (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_intrinsic * pwb
    )
{
 
  // ==============================================================================================
  // =                                Documentation (loops)                                       =
  // ==============================================================================================
  
  /**
   *   The intrinsic bispectrum consists in a 3D sum of a 4D integral. Our strategy is to tackle the
   * integrations separately from the sums. 
   *
   *   The integration goes over the magnitude of the three wavemodes (k1,k2,k3) and a numerical
   * parameter (r) of dimension time, which arises from the Rayleigh expansion of the Dirac delta
   * k1+k2+k3. This is analogous to the primordial reduced bispectra computed in the 'bispectra' module,
   * with the difference that here the second-order transfer function (which plays the role
   * of the primordial shape function) is always non-separable in (k1,k3,k3). We detail the integral
   * in the long comment inside the three loops, below.
   *
   *   The indices of the 3 summations are M3, L3 and L1. The azimuthal number M3 can assume any value
   * from -inf to +inf. In practice, only the first few values contribute, with M3=0 (scalar modes)
   * usually being the dominant contribution. For M3=0, the two other sums over L3 and L1 collapse
   * into Dirac deltas. The resulting 4D integral has almost the same form as the ones treated
   * in the bispectra module.
   *
   *   The sums over L3 and L1 have the following ranges:
   *     |l3-|M3|| <= L3 <= l3 + |M3|
   *     |l1-|M3|| <= L1 <= l1 + |M3|
   * and involve a set of four 3j-symbols and one 6j-symbol.
   * 
   *   We compute the 4D integral for all possible (l1,l2,l3) values and a single set of (M3,L3,L1).
   * We do so in the four functions bispectra2_intrinsic_integrate_over_XX where XX=k3,k2,k1,r. We
   * wnclose the integration in three loops over (M3,L3,L1). Since the L3 and L1 domains depend on
   * l3 and l1, we have to loop over offsets given by:
   * L3 = |l3-|M3|| + offset_L3
   * L1 = |l1-|M1|| + offset_L1
   *
   * This is a schematic of what we do:
   * 
   * for (index_M3)
   *   for (offset_L3)
   *     for (offset_L1)
   *   
   *       obtain the 4D integral for all (l1,l2,l3) and for the considered (M3,offset_L3,offset_L1)
   *       by integrating a function of (M3,L3,L1,l1,l2,l3,k3,k2,k1,r) over k3, k2, k1, r.
   *
   *       for (index_l1)
   *         for (index_l2)
   *           for (index_l3)
   *             
   *              bispectrum_l1_l2_l3 (M3,L3,L1) =
   *                integral_l1_l2_l3 (M3,offset_L3,offset_L1) * geometrical_factors (M3,offset_L3,offset_L1)
   *              
   *              bispectrum_l1_l2_l3 += bispectrum_l1_l2_l3 (M3,L3,L1)
   *
   *
   *    It is important to realise that this complicated sum reduces to a simple formula for scalar modes
   * (M3=0). The M3=0 formula is very similar to the usual primordial bispectrum in, for instance,
   * eq. 17 of Fergusson & Shellard 2007. The only difference is that here we have a non-separable
   * second-order transfer function (with angular dependence in k1,k2,k3) instead of a primordial
   * shape function.
   *
   */

  // ==============================================================================================
  // =                              Documentation (integration)                                   =
  // ==============================================================================================

  /**
   *   The formula for the 4D integral involves integrations in the following 4 variables:
   *
   *    - k3, where the second-order transfer functions T_l3(k3,k1,k2) is evaluated
   *    - k2, where the first-order transfer function T_l2(k2) is evaluated
   *    - k1, where the first-order transfer function T_l1(k1) is evaluated
   *    - r, that sets the frequency of the product of spherical Bessel functions j_l1(r*k1) * j_l2(r*k2) * j_l3(r*k3)
   *
   *   We estimate the reduced bispectrum in a straightforward way.  We first perform the integration on k3, involving
   * T_l3(k3,k1,k2) * j_l3(r*k3), on a grid in k1, k2, r.  The grid in k1 and k2 is coarse grained, as k1 and k2 are
   * the smooth directions of the second-order transfer function.  We also assume that the 'r' dependence of the 
   * integral is smooth, an assumption that is confirmed by the results.  We call the resulting integral
   * I_l3 (r,k1,k2), and store its values in pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2].
   *
   *   We then integrate I_l3(r,k1,k2) * T_l2(k2) * j_l2(r*k2) over k2.  As T_l2(k2) oscillates in k2 with Bessel
   * frequency, we need to interpolate I_l3 in k2.  This is not an issue as we know that it is smooth in k2.
   * We call the resulting integral I_l2_l3 (r,k1), and store its values in
   * pwb->integral_over_k2[index_l2][index_l3][index_r][index_k1].
   *
   *   Similarly, we convolve I_l2_l3 (r,k1) with T_l1(k1) * j_l1(r*k1) to obtain I_l1_l2_l3 (r), and store its values in
   * pwb->integral_over_k1[index_l1][index_l2][index_l3-index_l_triangular_min][index_r].
   *
   *   Finally, we integrate I_l1_l2_l3 (r) to obtain the bispectrum.
   *
   *   As we computed all the transfer functions and Bessel functions in a fixed grid, every integration is performed
   * using the simple trapezoidal rule.  Adopting more elaborated techniques would not necessarily give an increase
   * in precision, but they are likely to be slower.  As hinted in Fergusson & Shellard (2009), it could be interesting
   * to compute the geometrical integral involving the 3 Bessel functions separately, possibly using the analytical 
   * formulas presented in Mehrem, Londergan & Macfarlane (1990) (same as in Whelan, 1993).  Note however that, as
   * opposed to Fergusson & Shellard, we cannot identify a separable shape function as there is no such thing as
   * scale invariance for our second-order transfer functions.
   *
   */

  // =====================================================================================
  // =                               Bispectrum computation                              =
  // =====================================================================================

  /* Loop over the various intrinsic bispectra (so far, only one type) */
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* Skip the bispectrum if it not of the non-separable type */
    if (pbi->bispectrum_type[index_bt] != intrinsic_bispectrum)
      continue;        

    /* The XYZ indices refer to the considered field (T,E...). The first, X, refers to the
    second order perturbation, X=T^(2),E^(2)..., while Y and Z refer to the first-order ones:
    X -> second-order, Y -> first-order, Z -> first-order. This association ceases to be valid
    after we add the two extra bispectrum permutations, which effectively symmetrise the 
    bispectrum with respect to the position of the second-order field. 
    Nota also that, in this function, each field has a fixed wavemode and multipole associated:
    X -> (l3,k3), Y -> (l2,k2), Z -> (l1,k1). */

    for (int X=0; X < pbi->bf_size; ++X) {

      /* Do not consider the second-order part of Rayleigh scattering */
      if ((pbi->has_bispectra_r == _TRUE_) && (X == pbi->index_bf_r)) {
        if (pbi->bispectra_verbose > 0)
          printf(" -> will set to zero the %s bispectrum involving %s^(2)\n",
          pbi->bt_labels[index_bt], pbi->bf_labels[X]);
        continue;
      }

      pwb->X = X;

      if (pbi->bispectra_verbose > 0)
        printf(" -> computing %s bispectrum involving %s^(2)\n",
        pbi->bt_labels[index_bt], pbi->bf_labels[X]);

      // -------------------------------------------------------------------------------------
      // -                              Cycle over M3, L3, L1                                -
      // -------------------------------------------------------------------------------------
    
      for (int index_M3=0; index_M3 < ppr2->m_size; ++index_M3) {

        /* Update the structure with the current value of M3 */
        pwb->M3 = ppr2->m[index_M3];
        pwb->abs_M3 = abs(ppr2->m[index_M3]);

        /* Cycle on offset_L3, which is related to L3 by L3 = |l3-|M3|| + offset_L3 */
        for (int offset_L3=0; offset_L3 < (2*pwb->abs_M3+1); ++offset_L3) {

          pwb->offset_L3 = offset_L3;
        
          if (pbi->bispectra_verbose > 2)
            printf ("   \\ computing (M3=%d, offset_L3=%d) contribution to the bispectrum\n",
              pwb->M3, offset_L3);

          /* The geometrical factors in the bispectrum formula include the 3j symbol (l3,L3,|M3|)(M3,0,-M3).
          By itself, this 3j does not constrain offset_L3. However, for a given M3 we have a sum
          over |M3| and -|M3|. The only quantities that depend on the sign of M3 are:

          (l3,L3,|M3|)(M3,0,-M3) * \bar{T}_l3_M3
        
          For intensity and E-mode polarisation, \bar{T}_l3_-|M3| = \bar{T}_l3_|M3| which means that the bispectrum
          is proportional to (l3,L3,|M3|)(|M3|,0,-|M3|) + (l3,L3,|M3|)(-|M3|,0,|M3|). Switching the sign of the
          second line of a 3j introduces a factor (-1)^(l3+L3+|M3|) which is equal to (-1)^offset_L3. Hence,
          for intensity offset_L3 has to be even.
        
          For B-modes, \bar{T}_l3_-|M3| = - \bar{T}_l3_|M3| which means that the bispectrum is
          proportional to (l3,L3,|M3|)(|M3|,0,-|M3|) - (l3,L3,|M3|)(-|M3|,0,|M3|). Hence,
          for B-modes, offset_L3 has to be odd. */
          short skip_condition;

          if (pwb->bispectrum_parity == _EVEN_)
            skip_condition = (offset_L3%2!=0);
          else
            skip_condition = (offset_L3%2==0);

          if (skip_condition) {
            if (pbi->bispectra_verbose > 2)
              printf ("      \\ skipping offset_L3=%d for symmetry reasons\n", offset_L3);
            continue;
          }

          // -------------------------------------------------------------------
          // -                        Compute 4D integral                      -
          // -------------------------------------------------------------------

          /* Compute fist integral over k3 */
          class_call (bispectra2_intrinsic_integrate_over_k3(
                        ppr,
                        ppr2,
                        ppt,
                        ppt2,
                        pbs,
                        pbs2,
                        ptr,
                        ptr2,
                        ppm,
                        pbi,
                        pwb->index_tt2_of_bf[X],
                        index_M3,
                        offset_L3,
                        pwb),
            pbi->error_message,
            pbi->error_message);


          for (int Y=0; Y < pbi->bf_size; ++Y) {

            pwb->Y = Y;
        
            /* Compute second integral over k2 */
            class_call (bispectra2_intrinsic_integrate_over_k2(
                          ppr,
                          ppr2,
                          ppt,
                          ppt2,
                          pbs,
                          pbs2,
                          ptr,
                          ptr2,
                          ppm,
                          pbi,
                          pbi->index_tt_of_bf[Y],
                          pwb),
              pbi->error_message,
              pbi->error_message);

            /* Cycle on offset_L1, which is related to L1 by L1 = |l1-|M1|| + offset_L1 */
            for (int offset_L1=0; offset_L1 < (2*pwb->abs_M3+1); ++offset_L1) {
            
              pwb->offset_L1 = offset_L1;
            
              if (pbi->bispectra_verbose > 3)
                printf ("     * computing (M3=%d, offset_L3=%d, offset_L1=%d) contribution to the bispectrum\n",
                  pwb->M3, offset_L3, offset_L1);
            
              /* The geometrical factors involve a 3j symbol (l1,L1,|M3|)(0,0,0) which enforces
              that l1+L1+M3 is even. Hence offset_L1 =L1-|l1-|M3|| has to be even. This
              is valid for the B-modes bispectrum, too. */
              if (offset_L1%2!=0) { 
                if (pbi->bispectra_verbose > 2)
                  printf ("      \\ skipping offset_L1=%d for symmetry reasons\n", offset_L1);
                continue;
              }  
            
              /* Loop on the last first-order perturbation, k=T,E,... */
              for (int Z=0; Z < pbi->bf_size; ++Z) {

                pwb->Z = Z;
                          
                if ((pbi->bispectra_verbose > 1) && (ppr2->m_max_2nd_order==0))
                  printf("   \\ computing bispectrum %s_%s%s%s for m=%d\n",
                  pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z], pwb->M3);

                else if ((pbi->bispectra_verbose > 1) && (ppr2->m_max_2nd_order>0))
                  printf("   \\ computing bispectrum %s_%s%s%s for (m,offset_l3,offset_l1)=(%d,%d,%d)\n",
                  pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z],
                  pwb->M3, offset_L3, offset_L1);

                          
                /* Compute the third integral over k1 */
                class_call (bispectra2_intrinsic_integrate_over_k1 (
                              ppr,
                              ppr2,
                              ppt,
                              ppt2,
                              pbs,
                              pbs2,
                              ptr,
                              ptr2,
                              ppm,
                              pbi,
                              pbi->index_tt_of_bf[Z],
                              offset_L1,
                              pwb),
                  pbi->error_message,
                  pbi->error_message);
                      
                      
                /* Compute the fourth and last integral over r */
                class_call (bispectra2_intrinsic_integrate_over_r (
                              ppr,
                              ppr2,
                              ppt,
                              ppt2,
                              pbs,
                              pbs2,
                              ptr,
                              ptr2,
                              ppm,
                              pbi,
                              pwb),
                  pbi->error_message,
                  pbi->error_message);
                      
                // ===========================================================================
                // =                           Deal with geometry                            =
                // ===========================================================================
                          
                /* Compute the geometrical factors (3j's and a 6j) appearing outside of the integral.
                Store the result in pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2][index_l3] */
                class_call (bispectra2_intrinsic_geometrical_factors (
                              ppr,
                              ppr2,
                              ppt,
                              ppt2,
                              pbs,
                              pbs2,
                              ptr,
                              ptr2,
                              ppm,
                              pbi,
                              index_bt,
                              index_M3,
                              offset_L3,
                              offset_L1,
                              pwb->unsymmetrised_bispectrum[X][Y][Z], /* out */
                              pwb),
                  pbi->error_message,
                  pbi->error_message);
                          
              } // end of loop on field 'Z'
            } // end of loop on L1
              
          } // end of loop on field 'Y'
        } // end of loop on L3
      } // end of loop on M3  
    } // end of loop on field 'X'

    // ====================================================================================
    // =                               Consistency checks                                 =
    // ====================================================================================
  
    /* Variables useful for performing the checks below */
    int all_counter = 0;
    int non_negligible_counter = 0;
    int correct_counter = 0;
    int wrong_counter = 0;
    double incremental_diff_correct = 0.;
    double incremental_diff_wrong = 0.;
  
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
      for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
    
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
    
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
    
          for (int X=0; X < pbi->bf_size; ++X) {
            for (int Y=0; Y < pbi->bf_size; ++Y) {
              for (int Z=0; Z < pbi->bf_size; ++Z) {
          
                // -------------------------------------------------------------------
                // -                           Parity check                          -
                // -------------------------------------------------------------------
                
                /* For temperature and E-mode polarisation, check that when L=l1+l2+l3 is odd the bispectrum vanishes.
                For B-mode polarisation, check that when L=l1+l2+l3 is even the bispectrum vanishes. */
                
                double bispectrum = pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2][index_l3-index_l3_min];
                short parity = (pbi->field_parity[X] + pbi->field_parity[Y] + pbi->field_parity[Z])%2;
                short is_even_configuration = ((pbi->l[index_l1]+pbi->l[index_l2]+pbi->l[index_l3])%2==0);
    
                if (ppr2->m_max_2nd_order != 0) {
    
                  if (parity == _EVEN_) {
                    class_test ((!is_even_configuration) && (fabs(bispectrum)>_MINUSCULE_),
                      pbi->error_message,
                      "the computed bispectrum isn't parity invariant! Possible bug in the sum of 3js and 6js.");
                  }
                  else {
                    class_test ((is_even_configuration) && (fabs(bispectrum)>_MINUSCULE_),
                      pbi->error_message,
                      "the computed bispectrum isn't parity invariant! Possible bug in the sum of 3js and 6js.");
                  }
                }

                // -------------------------------------------------------------------
                // -                          l1 <-> l2 check                        -
                // -------------------------------------------------------------------

                /* Check that, when the two first-order transfer functions (Y and Z) are equal, like in TEE
                or ETT, the bispectrum does not change. */
                /* TODO: This kind of works to the 3% level when I include extrapolation. Can it be improved? */
                if (Y==Z) {

                  all_counter++;
                  int index_l2_min = pbi->index_l_triangular_min[index_l1][index_l3];
                  double b_1 = pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2][index_l3-index_l3_min];
                  double b_2 = pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l3][index_l2-index_l2_min];
    
                  if (fabs(b_2) > _MINUSCULE_) {
    
                    non_negligible_counter++;
                    double diff = 1-b_1/b_2;
    
                    /* Increment the counter of mismatches. We are happy with a 10% level match because
                    this is just an internal check. */
                    if (fabs(diff) > 0.05) {
                      wrong_counter++;
                      incremental_diff_wrong += fabs(diff);
                    }
                    else {
                      correct_counter++;
                      incremental_diff_correct += fabs(diff);
                    }
    
                    /* This is a very rough test. Generally the precision is of ~1e-6, but since the bispectrum
                    crosses the zero several times, here we just set the threshold to O(1) */
                    // if (fabs(diff) > 2)
                    //   printf ("b(%d,%d,%d)=%.17g is different from b(%d,%d,%d)=%.17g, diff=%g\n",
                    //   pbi->l[index_l2],pbi->l[index_l3],pbi->l[index_l1], b_1,
                    //   pbi->l[index_l1],pbi->l[index_l3],pbi->l[index_l2], b_2,
                    //   diff);
    
                  } // end of if bispectrum > _MINUSCULE_
                } // end of if Y==Z

              } // end of for(X)
            } // end of for(Y)
          } // end of for(Z)

        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_l1)
    
    /* We worry only if we have more than 1% of mismatches. We keep the treshold so high because
    the bispectrum crosses the zero many times. */
    class_test_permissive ((wrong_counter/(double)all_counter)>0.01,
      pbi->error_message,
      "l1<->l2 symmetry violated: n_configurations=%d, non_negligible=%d, wrong=%d, wrong/non_negl=%g, <diff> of matches=%g, <diff> of wrongs=%g\n",
      all_counter, non_negligible_counter, wrong_counter, wrong_counter/(double)non_negligible_counter,
      incremental_diff_correct/(double)correct_counter, incremental_diff_wrong/(double)wrong_counter);
  
  
    // ================================================================================================
    // =                                 Add bispectrum permutations                                  =
    // ================================================================================================

    /* At second-order, the bispectrum < X_l1 Y_l2 Z_l3 > is approximated by:
      < X^(2)_l1 Y^(1)_l2 Z^(1)_l3 >
    + < X^(1)_l1 Y^(2)_l2 Z^(1)_l3 >
    + < X^(1)_l1 Y^(1)_l2 Z^(2)_l3 >. 
    Here, we build this object from pwb->unsymmetrised_bispectrum. The ordering of its 6 levels is
    such that:
    < X^(2)_l1 Y^(1)_l2 Z^(1)_l3 > = pwb->unsymmetrised_bispectrum[X][Y][Z][l1][l2][l3]
    that is, the second-order transfer function always corresponds to the first field and to the
    the first multipole index of the unsymmetrised bispectrum array.
    We store the bispectrum for each combination of XYZ and only for those configurations where
    l1>=l2>=l3 and the triangular condition is satisfied. */

    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
    
      for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
    
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
    
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
    
          all_counter++;
    
          /* Index of the current (l1,l2,l3) configuration */
          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

          int index_l1_min = pbi->index_l_triangular_min[index_l2][index_l3];          
          int index_l2_min = pbi->index_l_triangular_min[index_l1][index_l3];
          
          for (int X=0; X < pbi->bf_size; ++X) {
            for (int Y=0; Y < pbi->bf_size; ++Y) {
              for (int Z=0; Z < pbi->bf_size; ++Z) {

                pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] =
                  pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2][index_l3-index_l3_min]
                + pwb->unsymmetrised_bispectrum[Y][X][Z][index_l2][index_l1][index_l3-index_l3_min]
                + pwb->unsymmetrised_bispectrum[Z][Y][X][index_l3][index_l2][index_l1-index_l1_min];

              } // end of for(Z)
            } // end of for(Y)
          } // end of for(X)

        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_l1)

  } // end of loop on index_bt
  
  return _SUCCESS_; 
  
}
  
  
  
  
  
  

int bispectra2_intrinsic_workspace_init (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_intrinsic * pwb
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

  /* Print the r-grid */
  // printf ("# ~~~ r-sampling for the bispectrum integration ~~~\n");
  // for (index_r=0; index_r < pwb->r_size; ++index_r)
  //   printf ("%d %g\n", index_r, pwb->r[index_r]);


  // -----------------------------------------------------------------------
  // -                           Grid in k1 and k2                         -
  // -----------------------------------------------------------------------

  /* In the second-order transfer functions T(k1,k2,k), k1 and k2 are the smooth directions */
  pwb->k_smooth_grid = ppt2->k;
  pwb->k_smooth_size = ppt2->k_size;
  
  
  
  // ========================================================================================
  // =                               Determine window function                              =
  // ========================================================================================

  /* Determine the window function for the interpolation of second-order transfer function
  in ptr->k1 and ptr->k2 */
  
  /* Window function will have pbi->k_smooth_size elements */
  class_alloc (pwb->k_window, pwb->k_smooth_size*sizeof(double), pbi->error_message);

  for (int index_k=0; index_k < pwb->k_smooth_size; ++index_k) {
  
    double k = pwb->k_smooth_grid[index_k];
    double pk = pbi->pk_pt[index_k];
  
    pwb->k_window[index_k] = 1./pow(k,2);
  }

  /* Inverse window function will have ptr->k_size[ppt->index_md_scalars] elements */
  class_alloc (pwb->k_window_inverse, ptr->k_size[ppt->index_md_scalars]*sizeof(double), pbi->error_message);

  for (int index_k=0; index_k < ptr->k_size[ppt->index_md_scalars]; ++index_k) {

    double k = ptr->k[ppt->index_md_scalars][index_k];
    double pk = pbi->pk[index_k];

    pwb->k_window_inverse[index_k] = pow(k,2);
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
  class_alloc (pwb->T_rescaling_factor, number_of_threads*sizeof(double*), pbi->error_message);
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
      number of k-values for the best sampled (k1,k2) wavemode in the transfer2 module */
    class_calloc_parallel(pwb->k3_grid[thread], ptr2->k3_size_max, sizeof(double), pbi->error_message);
    class_calloc_parallel(pwb->delta_k3[thread], ptr2->k3_size_max, sizeof(double), pbi->error_message);
    class_calloc_parallel(pwb->T_rescaling_factor[thread], ptr2->k3_size_max, sizeof(double), pbi->error_message);
  
  
    /* Allocate memory for the interpolation arrays (used only for the k2 and k3 integrations) */
    class_calloc_parallel (pwb->integral_splines[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);
    class_calloc_parallel (pwb->interpolated_integral[thread], ptr->k_size[ppt->index_md_scalars], sizeof(double), pbi->error_message);
    class_calloc_parallel (pwb->f[thread], pwb->k_smooth_size, sizeof(double), pbi->error_message);
  
  } // end of parallel region
  
  if (abort == _TRUE_) return _FAILURE_;
  
  
  
  /* Allocate array that will contain the unsymmetrised bispectrum. It contains 6 levels: three for
  the fields (e.g. TEE, TTE) and three for the multipoles (l1,l2,l3).  */
  pwb->count_allocated_for_unsymmetrised_bispectrum = 0;
  
  class_alloc (pwb->unsymmetrised_bispectrum, pbi->bf_size*sizeof(double *****), pbi->error_message);
  for (int X=0; X < pbi->bf_size; ++X) {
    class_alloc (pwb->unsymmetrised_bispectrum[X], pbi->bf_size*sizeof(double ****), pbi->error_message);
    for (int Y=0; Y < pbi->bf_size; ++Y) {
      class_alloc (pwb->unsymmetrised_bispectrum[X][Y], pbi->bf_size*sizeof(double ***), pbi->error_message);
      for (int Z=0; Z < pbi->bf_size; ++Z) {
        /* Allocate l1-level */
        class_alloc (pwb->unsymmetrised_bispectrum[X][Y][Z], pbi->l_size*sizeof(double **), pbi->error_message);
        for (int index_l1=0; index_l1<pbi->l_size; ++index_l1) {
          /* Allocate l2-level */
          class_alloc (pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1], pbi->l_size*sizeof(double *), pbi->error_message);
          for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
            /* Allocate l3-level (making sure to use calloc) */
            int l3_size = pbi->l_triangular_size[index_l1][index_l2];
            class_calloc (pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2], l3_size, sizeof(double), pbi->error_message);
            pwb->count_allocated_for_unsymmetrised_bispectrum += l3_size;
          } // end of for(index_l2)
        } // end of for(index_l1)
        
        /* Determine the parity of the considered bispectrum */
        pwb->bispectrum_parity = (pbi->field_parity[X] + pbi->field_parity[Y] + pbi->field_parity[Z])%2;

        /* We should not be computing an odd bispectrum with an all-odd grid, or viceversa, otherwise
        the bispectrum would be indentically zero */
        class_test (
          ((ppr->compute_only_even_ls==_TRUE_) && (pwb->bispectrum_parity == _ODD_)) ||
          ((ppr->compute_only_odd_ls==_TRUE_) && (pwb->bispectrum_parity == _EVEN_)),
          pbi->error_message,
          "computing an odd-parity bispectrum with an even l-grid, or viceversa.");
        
      } // end of for(k)
    } // end of for(j)
  } // end of for(i)
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> allocated ~ %.3g MB (%ld doubles) for the unsymmetrised bispectrum array\n",
      pwb->count_allocated_for_unsymmetrised_bispectrum*sizeof(double)/1e6, pwb->count_allocated_for_unsymmetrised_bispectrum);
    


  // =================================================================================================
  // =                                          Miscellaneous                                        =
  // =================================================================================================

  /* Associate to each field (T,E,...) its transfer function, which was computed in the transfer2.c module */
  for (int X = 0; X < pbi->bf_size; ++X) {
    
    if ((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t)) {
      pwb->index_tt2_of_bf[X] = ptr2->index_tt2_T;
    }
    else if ((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e)) {
      pwb->index_tt2_of_bf[X] = ptr2->index_tt2_E;
    }
    else if ((pbi->has_bispectra_b == _TRUE_) && (X == pbi->index_bf_b)) {
      pwb->index_tt2_of_bf[X] = ptr2->index_tt2_B;
    }
    /* The Rayleigh bispectrum will be computed ignoring its second-order part */
    else if ((pbi->has_bispectra_r == _TRUE_) && (X == pbi->index_bf_r)) {
      pwb->index_tt2_of_bf[X] = -1;
    }
    else {
      class_stop (pbi->error_message, "%s bispectrum not implemented yet", pbi->bf_labels[X]);
    }
  }


  /* Compute the m-dependent coefficient that enters the rescaling of the second-order transfer
  function in the bispectrum integral. */
  for (int index_M3=0; index_M3 < ppr2->m_size; ++index_M3) {

    int M3 = ppr2->m[index_M3];

    double fact_absM3=1, fact_2absM3=1;
    int copy=M3;
    while (copy > 0) fact_absM3 *= copy--;
    copy=2*M3;
    while (copy > 0) fact_2absM3 *= copy--;
    pwb->M3_coefficient[index_M3] = pow(2.,abs(M3)) * fact_absM3 / sqrt(fact_2absM3);

    /* Check that the rescaling coefficient is unity for m=0, sqrt(2) for m=1,
    and 8/sqrt(24) for m=2 */
    class_test (((M3==0) && (fabs(1-pwb->M3_coefficient[index_M3]/1)>_SMALL_))
      || ((M3==1) && (fabs(1-pwb->M3_coefficient[index_M3]/sqrt_2)>_SMALL_))
      || ((M3==2) && (fabs(1-pwb->M3_coefficient[index_M3]/(8./sqrt(24)))>_SMALL_)),
      pbi->error_message,
      "error in the computation of the transfer rescaling");
  }

  return _SUCCESS_;

}
  
  
  
  
  

  
int bispectra2_intrinsic_workspace_free(
    struct perturbs2 * ppt2,
    struct transfers2 * ptr2,
    struct bispectra * pbi,
    struct bispectra_workspace_intrinsic * pwb
    )
{

  free (pwb->r);
  free (pwb->delta_r);

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
    free(pwb->T_rescaling_factor[thread]);
    free(pwb->integral_splines[thread]);
    free(pwb->interpolated_integral[thread]);
    free(pwb->f[thread]);
    
  }  if (abort == _TRUE_) return _FAILURE_;
  
  free(pwb->k3_grid);
  free(pwb->delta_k3);
  free(pwb->T_rescaling_factor);
  free(pwb->integral_splines);
  free(pwb->interpolated_integral);
  free(pwb->f);
 
  /* Free pwb->unsymmetrised bispectrum */
  for (int X=0; X < pbi->bf_size; ++X) {
   for (int Y=0; Y < pbi->bf_size; ++Y) {
     for (int Z=0; Z < pbi->bf_size; ++Z) {
       for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
         for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
           free (pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1][index_l2]);
         } // end of for(index_l2)
         free (pwb->unsymmetrised_bispectrum[X][Y][Z][index_l1]);
       } // end of for(index_l1)
       free (pwb->unsymmetrised_bispectrum[X][Y][Z]);
     } // end of for(k)
     free (pwb->unsymmetrised_bispectrum[X][Y]);
   } // end of for(j)
   free (pwb->unsymmetrised_bispectrum[X]);
  } // end of for(i)
  free (pwb->unsymmetrised_bispectrum);
 
  return _SUCCESS_; 
  
}

  




int bispectra2_intrinsic_integrate_over_k3 (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_tt2_k3,
    int index_M3,
    int offset_L3,
    struct bispectra_workspace_intrinsic * pwb
    )
{

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

      /* Allocate 'k2' level. Make sure to use calloc. */
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
  
  /* We shall keep track of the average size of the integration grid in k3, which is (k1,k2) dependent */
  double average_k3_grid_size = 0;
  
  /* We compute the integral over k3 for all possible l-values */
  for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

    int l3 = pbi->l[index_l3];
    int L3 = abs(l3-pwb->abs_M3) + offset_L3;
    
    /* Skip the (L3,l3,M3) configurations forbidden by the triangular condition */
    if (L3 > (l3+pwb->abs_M3)) {
      pwb->count_memorised_for_integral_over_k3 += 0.5*pwb->k_smooth_size*(pwb->k_smooth_size+1)*pwb->r_size;
      continue;
    }  

    /* Find L3 inside pbs2->l1, the list of l's where we computed the Bessel functions */
    int index_L3 = pbs2->index_l1[L3];

    /* Paranoid android */
    class_test ((ppr2->m_max_2nd_order==0) && (l3!=L3),
      pbi->error_message,
      "inconsistency! for scalar modes, L3 must be equal to l3");

    class_test (pbs2->l1[index_L3]!=L3,
      pbi->error_message,
      "error in the indexing of pbs2->l1. Is the pbs2->extend_l1_using_m parameter true?");

    /* Load the transfer functions from disk */
    if ((ptr2->load_transfers_from_disk == _TRUE_) || (ptr2->store_transfers_to_disk == _TRUE_)) {
      class_call (transfer2_load_transfers_from_disk (
                    ppt2,
                    ptr2,
                    index_tt2_k3 + lm_cls(index_l3, index_M3)),
        ptr2->error_message,
        pbi->error_message);
    }
  
    if (pbi->bispectra_verbose > 2)
      printf("     * computing the k3 integral for l3=%d, index_l3=%d, L3=%d, index_L3=%d\n",
        l3, index_l3, pbs2->l1[index_L3], index_L3);

    abort = _FALSE_;
    #pragma omp parallel shared (abort) private (thread)
    {
  
      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
  
      #pragma omp for schedule (dynamic)
      for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {
  
        double k1 = pwb->k_smooth_grid[index_k1];

        /* We only need to consider those k2's that are equal or smaller than k1,
        as the quadratic sources were symmetrised  in the perturbation2 module */
        for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {

          double k2 = pwb->k_smooth_grid[index_k2]; 

          if (pbi->bispectra_verbose > 4)
            printf ("      \\ computing (index_k1,index_k2)=(%4d,%4d), (k1,k2)=(%g,%g)\n", index_k1, index_k2, k1, k2);
  
          // ===================================================
          // =            Fix the integration domain           =
          // ===================================================

          int dump;

          /* Compute the integration grid in the k3 variable */
          class_call_parallel (transfer2_get_k3_list (
                                 ppr,
                                 ppr2,
                                 ppt2,
                                 pbs,
                                 pbs2,
                                 ptr2,
                                 index_k1,
                                 index_k2,
                                 pwb->k3_grid[thread],  /* output */
                                 &dump   
                                 ),
            ptr2->error_message,
            pbi->error_message);

          /* Get the size of the integration grid. Note that when extrapolation is turned on, the k3-grid will
          also include values that do not satisfty the triangular condition k1 + k2 = k3. */
          int k3_size = ptr2->k_size_k1k2[index_k1][index_k2];
          
          /* Determine the measure for the trapezoidal rule for k3 */  
          pwb->delta_k3[thread][0] = pwb->k3_grid[thread][1] - pwb->k3_grid[thread][0];
            
          for (int index_k3=1; index_k3<(k3_size-1); ++index_k3)
            pwb->delta_k3[thread][index_k3] = pwb->k3_grid[thread][index_k3+1] - pwb->k3_grid[thread][index_k3-1];
            
          pwb->delta_k3[thread][k3_size-1] = pwb->k3_grid[thread][k3_size-1] - pwb->k3_grid[thread][k3_size-2];

          /* Let's be super cautious */
          for (int index_k3=0; index_k3<k3_size; ++index_k3)
            class_test_parallel (pwb->delta_k3[thread][index_k3] < 0,
              pbi->error_message,
              "something went terribly wrong, negative trapezoidal measure for index_k1=%d, index_k2=%d :-/",
              index_k1, index_k2);

          /* Define the pointer to the second-order transfer function as a function of k3. */
          double * transfer = ptr2->transfer[index_tt2_k3 + lm_cls(index_l3,index_M3)]
                              [index_k1]
                              [index_k2];
          
          // ===================================================
          // =                     Integrate                   =
          // ===================================================
  
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
          
            /* It is important to note that the used Bessel here has order L3 rather than l3 */
            class_call_parallel (bessel2_convolution (
                                   ppr,
                                   pbs2,
                                   pwb->k3_grid[thread],
                                   pwb->delta_k3[thread],
                                   k3_size,
                                   transfer,
                                   NULL,
                                   index_L3,
                                   pwb->r[index_r],
                                   &(pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]),
                                   pbi->error_message
                                   ),
              pbi->error_message,
              pbi->error_message);

            /* Check when 'm' is odd and k1=k2 then T(k1,k2,k3) is small with respect to 1 (which is
            the same as assuming that the bispectrum is small with respect to A_s*A_s). */
            /* TODO: It is a good idea to have this check, but here probably it's not the best place.
            Maybe we can move the check in the perturbations module. There, the order of the perturbations
            is O(1) hence we could just use _SMALL_ as an absolute epsilon. */
            if ((index_k1==index_k2) && (pwb->abs_M3%2!=0)) {
              double characteristic_scale = ppm->A_s*ppm->A_s;
              double epsilon = 1e-4*characteristic_scale;
              class_test_permissive (
                fabs (pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]) > epsilon,
                ppt2->error_message,
                "k1=k2, m odd, but bispectrum=%20.12g is larger than the characteristic scale (%20.12g). Problem?",
                pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2], epsilon);
            }

            /* Include the M3-dependent coefficient coming from the rescaling of the transfer function.
            Note that this is always equal to 1 for m=0. */
            if (pwb->abs_M3!=0)
              pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2] *= pwb->M3_coefficient[index_M3];

            /* Multiply the result by the factor 2 coming from the fact that in the bispectrum formula
            the second-order transfer functions appears as
            T(\vec{k1},\vec{k2},\vec{3}) + T(\vec{k2},\vec{k1},\vec{3}) */
            pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2] *= 2;

            /* Update the counters */
            #pragma omp atomic
            average_k3_grid_size += k3_size; 

            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k3;
              
            #pragma omp flush(abort)

            /* Print the integral as a function of r */
            // if ((pwb->abs_M3==1) && (offset_L3==0))
            //   if (l3==200)
            //     /* For bks run they correspond to 0.03728321 and 0.01724184, which are Christian's 13 and 6 */
            //     if ((index_k1==85) && (index_k2==63))
            //     // if ((index_k1==1) && (index_k2==0))
            //       fprintf (stderr, "%17.7g %17.7g\n",
            //         pwb->r[index_r], pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]);

            /* Print the integral as a function of k2 */
            // if (offset_L3==0)
            //  if (index_k1==13) /* To be used when you have the same grid as Christian */
            //   // if (index_k1==160) /* for bks runs it corresponds to 0.1912861 */
            //   if (index_k1==85) /* for bks runs it corresponds to 0.03728321, which is Christian's 13 */
            //   // if (index_k1==118) /* For bks runs, corresponds to 0.1 */
            //   // if (index_k1==95) /* For ref k-sampling, corresponds to 0.1 */
            //     // if (index_r==83) /* 14000 for 250r13to16 */
            //     if (index_r==49) /* 14000 for 99r135to145 */
            //       if ((l3==200) && (pwb->abs_M3==1))
            //         fprintf (stderr, "%17.7g %17.7g\n",
            //           k2, pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]);

          } // end of for(index_r)          
        } // end of for(index_k2)
      } // end of for(index_k1)
    } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */

  
    /* Free the memory associated with the second order transfer function for this (l,m) */
    if ((ptr2->load_transfers_from_disk == _TRUE_) || (ptr2->store_transfers_to_disk == _TRUE_)) {
      class_call (transfer2_free_type_level (
                    ppt2,
                    ptr2,
                    index_tt2_k3 + lm_cls(index_l3, index_M3)),
        ptr2->error_message,
        pbi->error_message);
    }
  
  } // end of for(index_l3);
    
  if (pbi->bispectra_verbose > 2)
    printf("     * memorised ~ %.3g MB (%ld doubles) for the k3-integral array (<k3_size>=%g)\n",
      pwb->count_memorised_for_integral_over_k3*sizeof(double)/1e6,
      pwb->count_memorised_for_integral_over_k3, average_k3_grid_size/pwb->count_memorised_for_integral_over_k3);

  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k3 != pwb->count_allocated_for_integral_over_k3,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k3, pwb->count_memorised_for_integral_over_k3);

  return _SUCCESS_; 
  
} // end of bispectra2_intrinsic_integrate_over_k3
  





int bispectra2_interpolate_over_k2 (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_k1,
    int index_l3,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_intrinsic * pwb
    )
{

  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];

  /* So far, we always assumed that k1>=k2 because the transfer function has the following symmetry:
  T (\vec{k1},\vec{k_2},\vec{k_3}) = (-1)^m * T (\vec{k_2},\vec{k_1},\vec{k_3}). The (-1)^m factor
  comes from the fact that one can exchange \vec{k_1} and \vec{k_2} by performing a rotation around
  \vec{k} of pi radians. This is equivalent to a factor of exp(i*pi*m) = (-1)^m, and also implies
  that T_lm (\vec{k_1},\vec{k_2},\vec{k_3}) = 0 when m is odd and \vec{k_1}=\vec{k_2}.
    
  Now we need to interpolate the k2 dependence of the transfer function, hence we build a temporary
  array where f(index_k1, index_k2) = (-1)^m * integral_over_k3(index_k2, index_k1) when
  index_k1 < index_k2.
  
  We also include a (k1/k2)^m factor that mimicks the rescaling of the transfer function
  under the exchange \vec{k_1}<->\vec{k_2}. This follows from the fact that
  k1*sin(theta_1) = k2*sin(theta_2), which implies that
  T_rescaled(\vec{k2},\vec{k1},\vec{k3}) = T_rescaled(\vec{k1},\vec{k2},\vec{k3}) * (k1/k2)^m */

  for (int index_k2=0; index_k2 < k_pt_size; ++index_k2) {

    if (index_k1 > index_k2)
      f[index_k2] = pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2];
    else 
      f[index_k2] = pwb->integral_over_k3[index_l3][index_r][index_k2][index_k1]
                  * pow (-k_pt[index_k1]/k_pt[index_k2], pwb->abs_M3);
    
    /* Multiply by window function */
    f[index_k2] *= pwb->k_window[index_k2];
                                        
  }
  
  if (ppr->transfers_k2_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
         pbi->error_message,
         pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (int index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    class_test((index_k+1) >= k_pt_size, pbi->error_message,
      "some of the elements in k_tr are larger than the largest k_pt. Stop to avoid seg fault.");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;

    /* Interpolate for each value of l3, r, k1 */
    if (ppr->transfers_k2_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k2_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1] +
        ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Further convolve with the primordial power spectrum */
    interpolated_integral[index_k_tr] *= pbi->pk[index_k_tr];

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
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]/pbi->pk[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  // }


  return _SUCCESS_;

}










int bispectra2_intrinsic_integrate_over_k2 (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_tt_k2,
    struct bispectra_workspace_intrinsic * pwb
    )
{

  /* Integration grid */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  // ======================================================================================================
  // =                                   Allocate memory for INT_l2_l3(r,k1)                              =
  // ======================================================================================================
  
  /* Because we shall recycle the array for more than one fields (T,E...) we make sure we allocate it
  only once. This is achieved by performing the allocation only at the beginning of the loop over
  the field (that is, when pwb->Y==0) */

  if (pwb->Y == 0) {

    /* Initialize counter */
    pwb->count_allocated_for_integral_over_k2 = 0;
  
    /* Allocate l2-level */
    class_alloc (pwb->integral_over_k2, pbi->l_size*sizeof(double ***), pbi->error_message);
  
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
  
      /* Allocate l3-level */
      class_alloc (pwb->integral_over_k2[index_l3], pbi->l_size*sizeof(double **), pbi->error_message);
  
      for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
    
        /* Allocate r-level */
        class_alloc (pwb->integral_over_k2[index_l3][index_l2], pwb->r_size*sizeof(double *), pbi->error_message);
  
        /* Allocate 'k1' level. Make sure to use calloc. */
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {
  
          int k1_size = pwb->k_smooth_size;
          class_calloc (pwb->integral_over_k2[index_l3][index_l2][index_r], k1_size, sizeof(double), pbi->error_message);
  
          /* Increase memory counter */
          pwb->count_allocated_for_integral_over_k2 += k1_size;
  
        } // end of for(index_r)
      } // end of for(index_l2)
    } // end of for(index_l3)
    
    if (pbi->bispectra_verbose > 2)
      printf("     * allocated ~ %.3g MB (%ld doubles) for the k2-integral array (l_size=%d)\n",
        pwb->count_allocated_for_integral_over_k2*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k2, pbi->l_size);

  } // end of if pwb->Y==0
  
  
  // ==============================================================================================================
  // =                                   Compute the INT_l2_l3(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k2 = 0;
  
  /* As for the k3_integral, we parallelize the loop over 'r'.  Here, we set it as the outermost loop
  as we do not need to load the second-order transfer functions from disk. */
  abort = _FALSE_;
  #pragma omp parallel shared (abort) private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k2 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);
  
      for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
          
        for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {
    
          /* Interpolate the integral I_l3(k1,k2,r) that we computed above in the integration grid of k2.
          Note that we pass the integral_splines and interpolated_integral arrays separately rather
          than accessing them from pwb, because they are thread dependent (while pwb isn't). */
          class_call_parallel (bispectra2_interpolate_over_k2 (
                                 ppr,
                                 ppr2,
                                 ppt,
                                 ppt2,
                                 pbs,
                                 pbs2,
                                 ptr,
                                 ptr2,
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
            
          /* Some debug - show the function that was interpolated as a function of r */
          // if ((pbi->l[index_l3]==205) && (pwb->abs_M3==1))
          //   if (pwb->offset_L3==1)
          //     if (index_k1==46) {
          //       index_k2 = 37;
          //       fprintf (stderr, "%17.7g %17.7g\n", pwb->r[index_r], pwb->f[thread][index_k2]);
          //     }
            
          /* Some debug - show the function that was interpolated as a function of k2 */
          // if (pwb->offset_L3==-1)
          //   if (index_k1==46)
          //     if (index_k2==37)
          //     // if (index_r==49) /* 14000 for 99r135to145 */
          //       // if ((pbi->l[index_l3]==205) && (pwb->abs_M3==1))
          //       //   for (index_k2=0; index_k2 < pwb->k_smooth_size; ++index_k2)
          //           fprintf (stderr, "%17.7g %17.7g\n",
          //             pwb->k_smooth_grid[index_k2], pwb->f[thread][index_k2]);

          /* Some debug - show the function after interpolation as a function of k2 */
          // if (pwb->offset_L3==0)
          //   if (index_k1==85) /* for ref runs with l_max=200 and kmax=6 it corresponds to 0.03728321, which is Christian's 13 */
          //   // if (index_k1==118) /* For better k-sampling, corresponds to 0.1 */
          //   // if (index_k1==95) /* For ref k-sampling, corresponds to 0.1 */
          //     // if (index_r==83) /* 14000 for 250r13to16 */
          //     if (index_r==49) /* 14005.1 for 100r135to145 */
          //       if ((pbi->l[index_l3]==200) && (pwb->abs_M3==1))
          //         for (index_k2=0; index_k2 < ptr->k_size[ppt->index_md_scalars]; ++index_k2)
          //           fprintf (stderr, "%17.7g %17.7g\n",
          //             ptr->k[ppt->index_md_scalars][index_k2],
          //             pwb->interpolated_integral[thread][index_k2]/pbi->pk[index_k2]);

          for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {  
  
            /* Define the pointer to the first-order transfer functions as a function of k */
            double * transfer = &(ptr->transfer [ppt->index_md_scalars]
                                                [((ppt->index_ic_ad * ptr->tt_size[ppt->index_md_scalars] + index_tt_k2)
                                                * ptr->l_size[ppt->index_md_scalars] + index_l2) * k_tr_size]);
  
            class_call_parallel (bessel_convolution (
                                   ppr,
                                   pbs,
                                   k_tr,
                                   pbi->delta_k,
                                   k_tr_size,
                                   pwb->interpolated_integral[thread],
                                   transfer,
                                   index_l2,
                                   pwb->r[index_r],
                                   &(pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1]),
                                   pbi->error_message
                                   ),
              pbi->error_message,
              pbi->error_message);

            /* Some debug - Print the integral as a function of r */
            // if ((pwb->abs_M3==1) && (pwb->offset_L3==0))
            //   if ((pbi->l[index_l2]==200) && (pbi->l[index_l3]==200))
            //     /* for l_max=200 and kmax=6 it corresponds to 0.03728321, which is Christian's 13 */
            //     if (index_k1==85)
            //     // if ((index_k1==1) && (index_k2==0))
            //       fprintf (stderr, "%17.7g %17.7g\n",
            //         pwb->r[index_r], pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1]);

            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k2;
  
            #pragma omp flush(abort)
  
          } // end of for(index_l2)
        } // end of for(index_k1)
      } // end of for(index_l3)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  if (pbi->bispectra_verbose > 2)
    printf("     * memorised ~ %.3g MB (%ld doubles) for the k2-integral array (k2_size=%d)\n",
      pwb->count_memorised_for_integral_over_k2*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k2,
      k_tr_size);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k2 != pwb->count_allocated_for_integral_over_k2,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k2, pwb->count_memorised_for_integral_over_k2);
  
  /* Free the memory that was allocated for the I_l3 integral, but only if we have already computed it
  for all the required probes */
  if (pwb->Y == (pbi->bf_size-1)) {
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
      for (int index_r=0; index_r < pwb->r_size; ++index_r) {
        for (int index_k1=0; index_k1<pwb->k_smooth_size; ++index_k1) {
          free (pwb->integral_over_k3[index_l3][index_r][index_k1]);
        } // end of for(index_k1)
        free (pwb->integral_over_k3[index_l3][index_r]);
      } // end of for(index_r)
      free (pwb->integral_over_k3[index_l3]);
    } // end of for(index_l3)
    free (pwb->integral_over_k3);
  } 

  return _SUCCESS_; 

}

        




int bispectra2_interpolate_over_k1 (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_l3,
    int index_l2,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_intrinsic * pwb
    )
{

  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];
  
  /* Define the function to be interpolated, and multiply it by a window function */
  for (int index_k1=0; index_k1 < k_pt_size; ++index_k1) {
    
    f[index_k1] = pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1];
    f[index_k1] *= pwb->k_window[index_k1];
  }
 
  if (ppr->transfers_k1_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
         pbi->error_message,
         pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (int index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    class_test((index_k+1) >= k_pt_size, pbi->error_message,
      "some of the elements in k_tr are larger than the largest k_pt. Stop to avoid seg fault.");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;
      
    /* Interpolate for each value of l2, l3, r */
    if (ppr->transfers_k1_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k1_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1]
        + ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Further convolve with the primordial power spectrum */
    interpolated_integral[index_k_tr] *= pbi->pk[index_k_tr];

    /* Revert the effect of the window function */
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  } // end of for (index_k_tr)


  /* Some debug - print the original array and the interpolation */
  // if ((index_l2==10) && (index_l3==11) && (index_r==49)) {
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
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]/pbi->pk[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  // }

  return _SUCCESS_;

}









    
        
        
int bispectra2_intrinsic_integrate_over_k1 (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_tt_k1,
    int offset_L1,
    struct bispectra_workspace_intrinsic * pwb
    )
{

  /* Integration grid */
  int k_tr_size = ptr->k_size[ppt->index_md_scalars];
  double * k_tr = ptr->k[ppt->index_md_scalars];

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  // ======================================================================================================
  // =                                   Allocate memory for INT_l2_l3(r,k1)                              =
  // ======================================================================================================
  
  /* We shall allocate the array pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min][index_r]
  so that the l1 level is the one satisfying the triangular inequality (|l2-l3| <= l1 <= l2+l3). 
  pwb->integral_over_k1 is recycled by the different iterations in the L1 loop and Z-field loop, hence we allocate
  it only at the first iteration (offset_L1==, pwb->Y==0 and pwb->Z==0) */
  
  if ((offset_L1 == 0) && (pwb->Y==0) && (pwb->Z==0)) {

    /* Initialize counter */
    pwb->count_allocated_for_integral_over_k1 = 0;
  
    /* Allocate l2-level */
    class_alloc (pwb->integral_over_k1, pbi->l_size*sizeof(double ***), pbi->error_message);
  
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
  
      /* Allocate l2-level */
      class_alloc (pwb->integral_over_k1[index_l3], pbi->l_size*sizeof(double **), pbi->error_message);
  
      for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
  
        /* Determine the limits for l1, which come from the triangular inequality |l3-l2| <= l1 <= l3+l2 */
        int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
        int l1_size = pbi->l_triangular_size[index_l3][index_l2];
  
        /* Allocate l1-level */
        class_alloc (pwb->integral_over_k1[index_l3][index_l2], l1_size*sizeof(double *), pbi->error_message);
      
        /* Allocate r-level. Make sure to use calloc. */
        for (int index_l1=index_l1_min; index_l1<(index_l1_min + l1_size); ++index_l1) {
  
          class_calloc (pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min],
                        pwb->r_size,
                        sizeof(double),
                        pbi->error_message);
        
          /* Increase memory counter */
          pwb->count_allocated_for_integral_over_k1 += pwb->r_size;
      
        } // end of for(index_l1)
      } // end of for(index_l2)
    } // end of for(index_l3)
    
    if (pbi->bispectra_verbose > 2)
      printf("     * allocated ~ %.3g MB (%ld doubles) for the k1-integral array (l_size=%d)\n",
        pwb->count_allocated_for_integral_over_k1*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k1, pbi->l_size);
  
  } // end of if(offset_L1==0 and pwb->Z==0)
  
  
  // ==============================================================================================================
  // =                                   Compute the INT_l2_l3(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k1 = 0;

  /* We parallelize the loop over 'r'. */
  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,ppt2,pbs,ptr,ptr2,ppm,pbi,pwb,abort)       \
    private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k1 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);

      /* We compute the integral over k1 for all possible l2-values */
      for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

        int l3 = pbi->l[index_l3];
  
        if (pbi->bispectra_verbose > 3)
          printf("      \\ computing the k1 integral for l3=%d, index_l3=%d, L1=%d, offset_L1=%d\n",
            pbi->l[index_l3], index_l3, pbs2->l1[offset_L1], offset_L1);
  
        /* We compute the integral over k1 for all possible l3-values */
        for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {  

          /* Determine the limits for l1, which come from the triangular inequality |l2-l2| <= l1 <= l2+l2 */
          int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
          int index_l1_max = pbi->index_l_triangular_max[index_l3][index_l2];
          
          /* DISABLED: see comment outside the outside the outermost loop. */
          /* We compute only those l1's whereby l1 >= l2. The reason is that the k1<->k2 symmetry of the
          second-order transfer functions enforces the l2<->l1 symmetry of I_l2_l3_l1(r). Note that
          this is a different symmetry than the full (l2,l3,l1) symmetry of the bispectrum, which we shall
          enforce later. If you comment out the following line together with the next one outside the
          outermost loop, the result should not change. */
          // index_l1_min = MAX (index_l2, index_l1_min);

          /* Interpolate the integral I_l2_l3(k1,r) that we computed above in the integration grid of k1 */
          if (index_l1_max >= index_l1_min) {
            
            class_call_parallel (bispectra2_interpolate_over_k1 (
                          ppr,
                          ppr2,
                          ppt,
                          ppt2,
                          pbs,
                          pbs2,
                          ptr,
                          ptr2,
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
              
          }
          else
            continue;


          /* l1 is the index of the first transfer function that we are going to integrate,
          while L1 is the order of the Bessel function */
          for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {
            
            int l1 = pbi->l[index_l1];
            int L1 = abs(l1-pwb->abs_M3) + offset_L1;

            /* Skip the (L3,l3,M3) configurations forbidden by the triangular condition */
            if (L1 > (l1+pwb->abs_M3)) {
              #pragma omp atomic
              ++pwb->count_memorised_for_integral_over_k1;
              continue;
            }  
      
            /* Find L1 inside pbs2->l1, the list of l's where we have computed the Bessel functions 
            for the extended l-range that is necessary to compute the m!=0 bispectra. */
            int index_L1 = pbs2->index_l1[L1];

            /* Paranoid android */
            class_test_parallel ((ppr2->m_max_2nd_order==0) && (l1!=L1),
              pbi->error_message,
              "inconsistency! for scalar modes, L1 must be equal to l1");
              
            class_test_parallel (pbs2->l1[index_L1]!=L1,
              pbi->error_message,
              "error in the indexing of pbs2->l1. Is the pbs2->extend_l1_using_m parameter true?");

            /* Define the pointer to the first-order transfer functions as a function of k */
            double * transfer = &(ptr->transfer [ppt->index_md_scalars]
                                                [((ppt->index_ic_ad * ptr->tt_size[ppt->index_md_scalars] + index_tt_k1)
                                                * ptr->l_size[ppt->index_md_scalars] + index_l1) * k_tr_size]);
  
            /* It is important to note that the used Bessel here has order L1 rather than l1 */
            class_call_parallel (bessel2_convolution (
                                   ppr,
                                   pbs2,
                                   k_tr,
                                   pbi->delta_k,
                                   k_tr_size,
                                   pwb->interpolated_integral[thread],
                                   transfer,
                                   index_L1,
                                   pwb->r[index_r],
                                   &(pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min][index_r]),
                                   pbi->error_message
                                   ),
              pbi->error_message,
              pbi->error_message);
  
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k1;

            /* Print the integral as a function of r */
            // if ((pwb->abs_M3==1) && (pwb->offset_L3==0) && (offset_L1==0))
            //     if ((l1==200) && (l2==200) && (pbi->l[index_l3]==200))
            //       fprintf (stderr, "%17.7g %17.7g\n",
            //         pwb->r[index_r], pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min][index_r]);

            #pragma omp flush(abort)
  
          } // end of for(index_l1)
        } // end of for(index_l2)
      } // end of for(index_l3)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region

  /* DISABLED:  While the optimization makes sense, the way we implement it below leads to nan's. Investigate. */
  /* Take care of the case when l1 < l2 */
  // for (index_r = 0; index_r < pwb->r_size; ++index_r) {
  //   for (index_l2=0; index_l2<pbi->l_size; ++index_l2) {
  //     for (index_l3=0; index_l3<pbi->l_size; ++index_l3) {
  // 
  //       int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
  //       int index_l1_max = pbi->index_l_triangular_max[index_l3][index_l2];
  // 
  //       /* The following is equivalent to restrict the loop to l1<l2 */
  //       index_l1_max = MIN (index_l2-1, index_l1_max);
  // 
  //       for (index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {
  // 
  //           int index_l2_min = pbi->index_l_triangular_min[index_l3][index_l1];
  //         
  //           pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min][index_r]
  //             = pwb->integral_over_k1[index_l1][index_l3][index_l2-index_l2_min][index_r];
  //         
  //           /* Update the counter */
  //           #pragma omp atomic
  //           ++pwb->count_memorised_for_integral_over_k1;
  //         
  //       }
  //     }
  //   }
  // }

  if (pbi->bispectra_verbose > 2)
    printf("     * memorised ~ %.3g MB (%ld doubles) for the k1-integral array (k1_size=%d)\n",
      pwb->count_memorised_for_integral_over_k1*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k1,
      k_tr_size);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k1 != pwb->count_allocated_for_integral_over_k1,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k1, pwb->count_memorised_for_integral_over_k1);

  /* Free the memory that was allocated for the integral over k2, but do that only when we are at
  the last iteration of the loops on offset_L1 (which is always performed since 2*abs_M3 is always
  even), Z and Y. */
  if ((offset_L1 == (2*pwb->abs_M3)) && (pwb->Y == (pbi->bf_size-1)) && (pwb->Z == (pbi->bf_size-1))) {
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
      for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {      
          free (pwb->integral_over_k2[index_l3][index_l2][index_r]);
        } // end of for(index_r)
        free (pwb->integral_over_k2[index_l3][index_l2]);
      } // end of for(index_l2)
      free (pwb->integral_over_k2[index_l3]);
    } // end of for(index_l3)
    free (pwb->integral_over_k2);
  } // end of if last term in the L1, Y, Z sums
  
  return _SUCCESS_;
}
      




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
    struct bispectra * pbi,
    struct bispectra_workspace_intrinsic * pwb
    )
{
  
  // ==========================================================================================================
  // =                                    Allocate memory for INT_l1_l2_l3(r)                                 =
  // ==========================================================================================================
  
  /* We shall allocate the array pwb->integral_over_r[index_l3][index_l2][index_l1-index_l_triangular_min]
  so that the l1 level is the one satisfying the triangular inequality. 
  pwb->integral_over_r is recycled by the different iterations in the L1 loop and in the Z-field loop, hence
  we allocate it only at the first iteration (offset_L1==0 and pwb->Z==0) */

  if ((pwb->offset_L1 == 0) && (pwb->Y == 0) && (pwb->Z == 0)) {
  
    /* Initialize counter */
    pwb->count_allocated_for_integral_over_r = 0;
  
    /* Allocate l3-level */
    class_alloc (pwb->integral_over_r, pbi->l_size*sizeof(double **), pbi->error_message);
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
      /* Allocate l2-level */
      class_alloc (pwb->integral_over_r[index_l3], pbi->l_size*sizeof(double *), pbi->error_message);
      for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
        /* Allocate l1-level */
        int l1_size = pbi->l_triangular_size[index_l3][index_l2];
        class_alloc (pwb->integral_over_r[index_l3][index_l2], l1_size*sizeof(double), pbi->error_message);
        pwb->count_allocated_for_integral_over_r += l1_size;
      } // end of for(index_l2)
    } // end of for(index_l3)

    if (pbi->bispectra_verbose > 2)
      printf(" -> allocated ~ %.3g MB (%ld doubles) for the r-integral array\n",
        pwb->count_allocated_for_integral_over_r*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_r);

  } // end of if(offset_L1==0 and pwb->Z == 0)
  
  
  // ==========================================================================================================
  // =                                    Compute  INT_l1_l2_l3(r)  integral                                  =
  // ==========================================================================================================
  
  /* We now proceed to the to integrate I(l3,l2,l1,r) over r. The function also multiplies
  the result by the appropriate coefficients.  */
    
  /* We parallelize the outer loop over 'l3'. */
  int abort = _FALSE_;
  #pragma omp parallel shared (ppt,ppt2,pbs,ptr,ptr2,ppm,pbi,pwb,abort)
  {
  
    #pragma omp for schedule (dynamic)
    for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the r-integral for l3=%d, index_l3=%d\n", pbi->l[index_l3], index_l3);
  
      /* We compute the integral over k2 for all possible l3-values */
      for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
  
        /* Determine the limits for l1, which come from the triangular inequality |l2-l3| <= l1 <= l2+l3 */
        int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
        int index_l1_max = pbi->index_l_triangular_max[index_l3][index_l2];
  
        for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {  

          if (pbi->bispectra_verbose > 4)
            printf("     * now considering (index_l3,index_l2,index_l1) = (%d,%d,%d)\n", index_l3, index_l2, index_l1);
        
          /* Shortcut to the main part of the integrand function */
          double * I = pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min];
          double integral = 0;
  
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
            double r = pwb->r[index_r];
            double integrand = r*r * I[index_r];

            if (integrand == 0.)
              continue;
      
            /* Increment the estimate of the integral */
            integral += integrand * pwb->delta_r[index_r];
    
            /* Some debug - output intermediate results on stderr for a custom (l2,l3,l1) configuration */
            // int index_l3 = pbi->l[index_l3];
            // int index_l2 = pbi->l[index_l2];
            // int index_l1 = pbi->l[index_l1];
            // if ( (l3==100) && (l2==200) && (l1==300) ) {
            //   if (index_r==0) {
            //     fprintf(stderr, "##########    l3 = %d, l2 = %d, l1 = %d, n_rows = %d    ##########\n\n",
            //       l3, l2, l1, r, pwb->r_size);
            //     fprintf(stderr, "%12s %17s %17s %17s\n", "r", "r^2*I_l2_l3_l1(r)", "integral","delta_r");
            //   }
            //   else {
            //     fprintf(stderr, "%12.7g %17.7g %17.7g %17.7g\n", r, integrand, integral, pwb->delta_r[index_r]);
            //   }
            // }
   
          } // end of for(index_r)
  
  
          /* Fill the result array and include the factor 1/2 from trapezoidal rule */
          pwb->integral_over_r[index_l3][index_l2][index_l1-index_l1_min] = 0.5 * integral;
  
          /* Some debug - output the integral as a function of r on stderr for a custom (l2,l3,l1) */
          // if ( (l2==l3) && (l3==l1) ) {
          //   fprintf(stderr, "%12d %17.7g\n", l2, pwb->integral_over_r[index_l3][index_l2][index_l1-index_l1_min]);
          // }
            
        } // end of for(index_l1)
      } // end of for(index_l2)
      
      #pragma omp flush(abort)
      
    } // end of for(index_l3)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  /* We can free the memory that was allocated for the integral over k1, as it is no longer needed */
  if ((pwb->offset_L1 == (2*pwb->abs_M3)) && (pwb->Y == (pbi->bf_size-1)) && (pwb->Z == (pbi->bf_size-1))) {
    for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
      for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
        int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
        int index_l1_max = pbi->index_l_triangular_max[index_l3][index_l2];
        for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {  
          free (pwb->integral_over_k1[index_l3][index_l2][index_l1-index_l1_min]);
        } // end of for(index_l1)
        free (pwb->integral_over_k1[index_l3][index_l2]);
      } // end of for(index_l3)
      free (pwb->integral_over_k1[index_l3]);
    } // end of for(index_l2)
    free (pwb->integral_over_k1);
  } // end of last term in the L1 and loop
  
  return _SUCCESS_; 
}




int bispectra2_intrinsic_geometrical_factors (
    struct precision * ppr,
    struct precision2 * ppr2,
    struct perturbs * ppt,
    struct perturbs2 * ppt2,
    struct bessels * pbs,
    struct bessels2 * pbs2,
    struct transfers * ptr,
    struct transfers2 * ptr2,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_M3,
    int offset_L3,
    int offset_L1,
    double *** result, /* out */
    struct bispectra_workspace_intrinsic * pwb
    )
{

  /* Considered azimuthal mode */
  int M3 = pwb->M3;
  int abs_M3 = pwb->abs_M3;

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

  // =======================================================================
  // =                    Allocate geometrical factors                     =
  // =======================================================================

  /* We will need to compute 5 geometrical factors, which correspond to four 3j-symbols
  and one 6j */
  enum geometrical_factors {
    l3_L3_M3,             /* three-j symbol (l3,L3,M3)(M3,0,-M3) */
    l1_L1_M3,             /* three-j symbol (l1,L1,M3)(0,0,0) */
    L1_l2_L3,             /* three-j symbol (L1,l2,L3)(0,0,0) */
    l1_l2_l3,             /* three-j symbol (l1,l2,l3)(0,0,0) */
    l1_l3_l2_L3_L1_M3     /* six-j symbol {l1,l3,l2}{L3,L1,M2} */
  };
  int n_geometrical_factors = 5;

  /* Temporary arrays and values needed to store the results of the 3j and 6j computations */
  int size[number_of_threads][n_geometrical_factors];
  int min[number_of_threads][n_geometrical_factors];
  int max[number_of_threads][n_geometrical_factors];
  double *** value;
  class_alloc (value, number_of_threads*sizeof(double **), pbi->error_message);
  for (thread=0; thread < number_of_threads; ++thread) {
    class_alloc (value[thread], n_geometrical_factors*sizeof(double *), pbi->error_message);
    for (int ii=0; ii < n_geometrical_factors; ++ii)
      class_alloc (value[thread][ii], (2*pbi->l_max+1)*sizeof(double), pbi->error_message);
  }

  // ========================================================================
  // =                             Cycle on l3                            =
  // ========================================================================
  #pragma omp parallel for private (thread)  
  for (int index_l3=0; index_l3 < pbi->l_size; ++index_l3) {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    int l3 = pbi->l[index_l3];
    int L3 = abs(l3-abs_M3) + offset_L3;

    /* Skip the (L3,l3,M3) configurations forbidden by the triangular condition */
    if (L3 > (l3+abs_M3))
      continue;
  
    /* Temporary variables to hold the limits of the 3j's in double precision format */
    double min_D, max_D;
  
    /* Compute the three-j symbol (l3,L3,M3)(M3,0,-M3) for all allowed values of L3.
    Mind the column positions! */
    class_call_parallel (drc3jj (
                           abs_M3, l3, -M3, M3,
                           &min_D, &max_D,
                           value[thread][l3_L3_M3],
                           (2*pbi->l_max+1),
                           pbi->error_message       
                           ),
      pbi->error_message,
      pbi->error_message);
  
    min[thread][l3_L3_M3] = (int)(min_D + _EPS_);
    max[thread][l3_L3_M3] = (int)(max_D + _EPS_);
    size[thread][l3_L3_M3] = max[thread][l3_L3_M3] - min[thread][l3_L3_M3] + 1;
  
    class_test_parallel ((L3 - min[thread][l3_L3_M3]) >= size[thread][l3_L3_M3],
      pbi->error_message,
      "error in the computation of the three-j symbol l3_L3_M3");

    /* Value of l3_L3_M3 in L3 */
    double FACTOR_l3_L3_M3 = value[thread][l3_L3_M3][L3 - min[thread][l3_L3_M3]];

    /* Debug (l3,L3,M3)(M3,0,-M3) */
    // printf ("I(l3,L3,|M3|)(M3,0,-M3) = (%d,%d,%d)(%d,0,%d) = %g\n",
    //     l3, L3, abs_M3, M3, -M3, FACTOR_l3_L3_M3);

    // ============================================================================
    // =                         Sum over -|M3| and +|M3|                         =
    // ============================================================================

    /* For a given value of |M3|, we need to sum over -|M3| and +|M3|. See the long
    comment in bispectra2_intrinsic_init (inside the offset_L3 loop) for details
    on what we do here.  */
    double SUMMED_FACTOR_l3_L3_M3 = FACTOR_l3_L3_M3;
    if (abs_M3 != 0) {
      
      if (pwb->bispectrum_parity == _EVEN_)
        SUMMED_FACTOR_l3_L3_M3 = FACTOR_l3_L3_M3 + alternating_sign(offset_L3)*FACTOR_l3_L3_M3;

      else
        SUMMED_FACTOR_l3_L3_M3 = FACTOR_l3_L3_M3 - alternating_sign(offset_L3)*FACTOR_l3_L3_M3;

    } // end of if(M3!=0)
  
    // =========================================================================
    // =                             Cycle on l2                             =
    // =========================================================================
    for (int index_l2=0; index_l2 < pbi->l_size; ++index_l2) {
      
      int l2 = pbi->l[index_l2];
            
      /* Compute the three-j symbol (l1,l2,l3)(0,0,0) for all allowed values of l1.
      If we are dealing with an odd bispectrum, compute (l1,l2,l3)(2,0,-2) instead */      
      int F = ((pwb->bispectrum_parity == _EVEN_) ? 0:2);

      class_call_parallel (drc3jj (
                             l2, l3, F, -F,
                             &min_D, &max_D,
                             value[thread][l1_l2_l3],
                             (2*pbi->l_max+1),
                             pbi->error_message       
                             ),
        pbi->error_message,
        pbi->error_message);

      min[thread][l1_l2_l3] = (int)(min_D + _EPS_);
      max[thread][l1_l2_l3] = (int)(max_D + _EPS_);
      size[thread][l1_l2_l3] = max[thread][l1_l2_l3] - min[thread][l1_l2_l3] + 1;

      /* Compute the three-j symbol (L1,l2,L3)(0,0,0) for all allowed values of L1 */
      class_call_parallel (drc3jj (
                             l2, L3, 0, 0,
                             &min_D, &max_D,
                             value[thread][L1_l2_L3],
                             (2*pbi->l_max+1),
                             pbi->error_message       
                             ),
        pbi->error_message,
        pbi->error_message);
      
      min[thread][L1_l2_L3] = (int)(min_D + _EPS_);
      max[thread][L1_l2_L3] = (int)(max_D + _EPS_);
      size[thread][L1_l2_L3] = max[thread][L1_l2_L3] - min[thread][L1_l2_L3] + 1;
      
      /* Allowed values of l1 */
      int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
      int index_l1_max = pbi->index_l_triangular_max[index_l3][index_l2];
  
      // =========================================================================
      // =                             Cycle on l1                             =
      // =========================================================================
      for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {
  
        int l1 = pbi->l[index_l1];
        int L1 = abs(l1-abs_M3) + offset_L1;
            
        /* Skip the (L1,l1,M3) configurations forbidden by the triangular condition */
        if (L1 > (l1+abs_M3))
          continue;

        /* Getting paranoid about those offset indices... */
        class_test_parallel ((!is_triangular_int(l1,L1,abs_M3)) || (!is_triangular_int(l3,L3,abs_M3)),
          pbi->error_message,
          "error with either offset_L1 or offset_L3 loop indexing");
      
        /* Enforce the triangular condition given by the (L1,l2,L3)(0,0,0) 3j-symbol. In terms of
        offset_L1 and offset_L3, the constraints excludes the cases when their difference is large
        with respect to l2. Note that there is no need to enforce the four triangular conditions
        arising from the 6j symbol, as they are the same as those for the 3j's. */
        if (!is_triangular_int(L1,l2,L3))
          continue;

        /* For an even/odd, the only non-vanishing contributions come from even/odd l1+l2+l3  */
        short is_even_configuration = ((l1+l2+l3)%2==0);

        /* Value of l1_l2_l3 in L1 */
        class_test_parallel ((l1 - min[thread][l1_l2_l3]) >= size[thread][l1_l2_l3],
          pbi->error_message,
          "error in the computation of the three-j symbol l1_l2_l3");

        double FACTOR_l1_l2_l3 = value[thread][l1_l2_l3][l1 - min[thread][l1_l2_l3]];

        /* The reduced bispectrum is given by the angle averaged bispectrum divided by FACTOR_l1_l2_l3.
        Here we make sure that FACTOR_l1_l2_l3 is not zero. We do not worry if it is zero for even 
        bispectra and odd l1+l2+l3, because in that case it must vanish and we cannot define a
        reduced bispectrum. */
        if (!((pwb->bispectrum_parity == _EVEN_) && (is_even_configuration==_FALSE_)))
          class_test_permissive (fabs(FACTOR_l1_l2_l3) < _MINUSCULE_,
            pbi->error_message,
            "possibility of having nans, caution! (l1,l2,l3)=(%d,%d,%d), F=%d, 3J=%g",
            l1, l2, l3, F, FACTOR_l1_l2_l3);

        /* Debug l1_l2_l3 */
        // printf ("I(l1,l2,l3)(%d,0,%d) = (%d,%d,%d)(%d,0,%d) = %g\n",
        //   F, -F, l1, l2, l3, F, -F, FACTOR_l1_l2_l3);

        /* Value of L1_l2_L3 in L1 */
        class_test_parallel ((L1 - min[thread][L1_l2_L3]) >= size[thread][L1_l2_L3],
          pbi->error_message,
          "error in the computation of the three-j symbol L1_l2_L3 (L1=%d, min=%d, size=%d)",
          L1, min[thread][L1_l2_L3], size[thread][L1_l2_L3]);

        double FACTOR_L1_l2_L3 = value[thread][L1_l2_L3][L1 - min[thread][L1_l2_L3]];
        
        /* Debug (L1,l2,L3)(0,0,0) */
        // if (offset_L1 == 4)
        //   printf ("I(L1,l2,L3)(0,0,0) = (%d,%d,%d)(0,0,0) = %g\n",
        //       L1, l2, L3, FACTOR_L1_l2_L3);

        /* Compute the three-j symbol (l1,L1,M3)(0,0,0) for all allowed values of L1 */
        class_call_parallel (drc3jj (
                               l1, abs_M3, 0, 0,
                               &min_D, &max_D,
                               value[thread][l1_L1_M3],
                               (2*pbi->l_max+1),
                               pbi->error_message       
                               ),
          pbi->error_message,
          pbi->error_message);        
        
        min[thread][l1_L1_M3] = (int)(min_D + _EPS_);
        max[thread][l1_L1_M3] = (int)(max_D + _EPS_);
        size[thread][l1_L1_M3] = max[thread][l1_L1_M3] - min[thread][l1_L1_M3] + 1;
        
        /* Value of l1_L1_M3 in L1 */
        class_test_parallel ((L1 - min[thread][l1_L1_M3]) >= size[thread][l1_L1_M3],
          pbi->error_message,
          "error in the computation of the three-j symbol l1_L1_M3");
          
        double FACTOR_l1_L1_M3 = value[thread][l1_L1_M3][L1 - min[thread][l1_L1_M3]];

        /* Debug (l1,L1,M3)(0,0,0) */
        // printf ("I(l1,L1,M3)(0,0,0) = (%d,%d,%d)(0,0,0) = %g\n",
        //   l1, L1, abs_M3, FACTOR_l1_L1_M3);

        /* Compute the six-j symbol {l1,l3,l2}{L3,L1,|M3|} for all allowed values of L1. In order to fit
        it in the Slatec function, we recast it as {L1,L3,l2}{l3,l1,|M3|} */
        class_call_parallel (drc6j (
                               /*L1,*/ L3, l2, l3, l1, abs_M3,
                               &min_D, &max_D,
                               value[thread][l1_l3_l2_L3_L1_M3],
                               (2*pbi->l_max+1),
                               pbi->error_message       
                               ),
          pbi->error_message,
          pbi->error_message);

        min[thread][l1_l3_l2_L3_L1_M3] = (int)(min_D + _EPS_);
        max[thread][l1_l3_l2_L3_L1_M3] = (int)(max_D + _EPS_);
        size[thread][l1_l3_l2_L3_L1_M3] = max[thread][l1_l3_l2_L3_L1_M3] - min[thread][l1_l3_l2_L3_L1_M3] + 1;

        /* Value of {l1,l3,l2}{L3,L1,|M3|} in L1 */
        class_test_parallel ((L1 - min[thread][l1_l3_l2_L3_L1_M3]) >= size[thread][l1_l3_l2_L3_L1_M3],
          pbi->error_message,
          "error in the computation of the three-j symbol l1_l3_l2_L3_L1_M3");
        
        double FACTOR_l1_l3_l2_L3_L1_M3 =
          value[thread][l1_l3_l2_L3_L1_M3][L1 - min[thread][l1_l3_l2_L3_L1_M3]];

        /* Debug {l1,l3,l2}{L3,L1,|M3|} */
        // if (offset_L3 == 2)
        //   if (offset_L1 == 2)
        //     printf ("I{l1,l3,l2}{L3,L1,|M3|} = {%d,%d,%d}{%d,%d,%d} = %g\n",
        //       l1, l3, l2, L3, L1, abs_M3, FACTOR_l1_l3_l2_L3_L1_M3);

        // ==========================================================================
        // =                        Increment the bispectrum                        =
        // ==========================================================================

        /* Four-dimensional integral */
        double integral = pwb->integral_over_r[index_l3][index_l2][index_l1-index_l1_min];

        /* Prefactor for the reduced bispectrum. This already has a factor
        sqrt((2l1+1.)*(2l2+1.)*(2l3+1.)/(4*_PI_)) taken out, and it
        is multiplied by (2*l1+1)*(2*l2+1)*(2*l3+1) with respect
        to the one in Christian's notes because in Song the transfer functions
        are the Legendre ones. */
        double prefactor;
        
        if (pwb->bispectrum_parity == _EVEN_) {
          prefactor = 8/_PI_CUBE_ * (2*L1+1.) * (2*L3+1.);
        }
        else {
          /* The last line is the P_l -> Y_lm factor, without 2*l3+1 */
          prefactor = - 0.5 / (4*_PI_FOURTH_) * (2*L1+1.) * (2*L3+1.) * (2*l1+1) * (2*l2+1);
        }
        
        /* Product of all the geometrical factors */
        double geometry = SUMMED_FACTOR_l3_L3_M3 * FACTOR_l1_L1_M3 * FACTOR_L1_l2_L3 * FACTOR_l1_l3_l2_L3_L1_M3;

        /* Check parity for T and E-modes */
        if (pwb->bispectrum_parity == _EVEN_) {
          class_test_parallel ((!is_even_configuration) && (fabs(geometry)>_MINUSCULE_),
            pbi->error_message,
            "geometry does not respect parity: G_%d_%d_%d=%g! Possible bug in the sum of 3js and 6js.",
            l1, l2, l3, geometry);
        }
        else {
          /* Check parity for B-modes */
          class_test_parallel ((is_even_configuration) && (fabs(geometry)>_MINUSCULE_),
            pbi->error_message,
            "geometry does not respect parity: G_%d_%d_%d=%g! Possible bug in the sum of 3js and 6js.",
            l1, l2, l3, geometry);
        }

        /* Divide by the 3j symbol (l1,l2,l3)(0,0,0) in order to obtain the reduced bispectrum. Note that
        for odd values of l1+l2+l3, the ratio will give a nan, so in that case we leave the geometry
        factor as it is; this is is not an issue because for odd l1+l2+l3, the geometry factor vanishes.
        For an odd bispectrum, divide by (l1,l2,l3)(-2,0,2), which is different from zero both for 
        odd and even l1+l2+l3. */
        if (pwb->bispectrum_parity == _EVEN_) {
          if (is_even_configuration) 
            geometry /= FACTOR_l1_l2_l3;
        }
        /* DISABLED: we don't need to take the (l1 l2 l3)(2 0 -2) explicitly */
        // else
        //   geometry *= FACTOR_l1_l2_l3;
          



        // ***   Alternating sign   ***

        /* The bispectrum formula has a prefactor of the form (-i)^(l3+L3-l1-L1). The exponent has the same
        parity of offset_L1+offset_L3 (offset_L1=L1-|l1-|M3|| and offset_L3=L3-|l3-|M3||). For intensity,
        both offsets have to be even, which ensures that the prefactor, and thus the whole bispectrum, is
        real.

        For B-modes, offset_L1 is still even but offset_L3 needs to be odd (see the comment in the 
        offset_L3 loop in bispectra2_intrinsic_init). This makes the overall factor always imaginary,
        which is the way it should be (we are computing <conj(B)IE> which is needed to obtain the
        C_l's) */
        int exponent;

        if (pwb->bispectrum_parity == _EVEN_)
          exponent = l3+L3-l1-L1;
        else
          exponent = l3+L3-l1-L1-1;

        class_test_parallel ((exponent%2!=0) && (fabs(geometry)>_MINUSCULE_),
          pbi->error_message,
          "imaginary bispectrum, something is wrong with the 3js");

          
        /* Add the contribution from (l1,l2,l3,M3,L3,L1) to the bispectrum */
        /* TODO: the m_max check here should be removed. We want a unique way to define the
        reduced bispectrum regardless of the m */
        if (pwb->bispectrum_parity == _EVEN_) {
          if (ppr2->m_max_2nd_order != 0) {

            #pragma omp atomic
            result[index_l3][index_l2][index_l1-index_l1_min]
              += alternating_sign(exponent/2) * prefactor * geometry * integral;

          }
          /* We treat m=0 as a special case because in this way we can define the reduced
          bispectrum for m=0 also when l1+l2+l3 is odd. We cannot do this for the m>0 cases
          because there we have to numerically divide the bispectrum formula by the 3j-symbol
          (l1,l2,l3)(0,0,0) rather than extracting it analitically */
          else {
        
            result[index_l3][index_l2][index_l1-index_l1_min]
              = 8/_PI_CUBE_ * integral;
        
          }
          
          /* Check that for scalar modes the integral simplifies to a simpler formula, similar
          to the one in eq. 17 of Fergusson & Shellard (2007) */
          if (abs_M3 == 0) {

            /* Simple formula with the (l1,l2,l3) 3j-symbol analytically extracted */
            double bispectrum_for_m0 = 8/_PI_CUBE_ * integral;

            /* Full formula, with the (l1,l2,l3) 3j-symbol factored out numerically */
            double bispectrum = 0;
            if (is_even_configuration)
              bispectrum = alternating_sign(exponent/2) * prefactor * geometry * integral;          

            if (fabs(bispectrum) > _MINUSCULE_)
              class_test_parallel (fabs(1-bispectrum_for_m0/bispectrum) > ppr->smallest_allowed_variation,
                pbi->error_message,
                "error in the geometrical factor, m=0 not recovered (%g != %g) for b_%d_%d_%d.",
                bispectrum_for_m0, bispectrum, l2, l3, l1);
          }
        } // end of if even parity
        else {

          #pragma omp atomic
          result[index_l3][index_l2][index_l1-index_l1_min]
            += alternating_sign(exponent/2) * prefactor * geometry * integral;
          
        } // end of if odd parity

        /* Comparison with Christian's B-modes */
        // if (pwb->abs_M3==1) {
        //   if ((l1>200) && (l1<310)) {
        // 
        //     double geometry_prefactor = - 0.5 / (4*_PI_FOURTH_) * (2*L1+1.) * (2*L3+1.) / (2.*l3+1);
        //     double integral_prefactor = (2.*l1+1) * (2.*l2+1) * (2.*l3+1);
        //   
        //     #pragma omp critical
        //     fprintf (stderr, "(l1,l2,l3 | L1,L3)=(%3d,%3d,%3d | %3d,%3d): geometry=%10g, integral=%10g\n",
        //       l1,l2,l3,L1,L3,
        //       alternating_sign(exponent/2) * geometry_prefactor * geometry,
        //       integral_prefactor*integral);
        // 
        //   }
        // }
        /* Some debug */
        // printf ("(M3,L3,L1,l1,l2,l3,offset_L3,offset_L1)=(%d,%d,%d,%d,%d,%d,%d,%d):",
        //   ppr2->m[index_M3],L3,L1,l1,l2,l3,offset_L3,offset_L1);
        // printf ("\tprefactor=%g, integral=%g, geometry=%g, product=%g\n",
        //   alternating_sign(exponent/2)*prefactor, integral, geometry,
        //   alternating_sign(exponent/2)*alternating_sign(exponent)*integral*geometry);
        
        /* Mathematica batch debug */
        // double val = alternating_sign(exponent/2)*geometry;
        // if ((abs_M3==2) && (offset_L3==4)) {
        //   if (val==0) {
        //     fprintf (stderr, "Abs[g[%d,%d,%d,%d,%d,%d]] == 0,\n",
        //       abs_M3,L3,L1,l1,l2,l3);
        //   }
        //   else {
        //     fprintf (stderr, "Abs[1-g[%d,%d,%d,%d,%d,%d]/((%g)*10^(%f))],\n",
        //       abs_M3,L3,L1,l1,l2,l3,sign(val),log10(fabs(val)));
        //   }
        // }

        /* Print coefficient */
        // printf ("coeff_%d_%d_%d(M3=%d,L3=%d,L1=%d) = %g\n", l1, l2, l3,
        //   prefactor/((2*l1+1)*(2*l2+1)*(2*l3+1)));

        #pragma omp flush(abort)
  
      } // end of for(index_l1)
    } // end of for(index_l2)
  } // end of for(index_l3) and of parallel region
  if (abort == _TRUE_) return _FAILURE_;


  /* We can free the memory that was allocated for the integral over r, as it is no longer needed */
  if ((offset_L1 == (2*pwb->abs_M3)) && (pwb->Y == (pbi->bf_size-1)) && (pwb->Z == (pbi->bf_size-1))) {
    for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
      for (int index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
        free (pwb->integral_over_r[index_l3][index_l2]);
      } // end of for(index_l2)
      free (pwb->integral_over_r[index_l3]);
    } // end of for(index_l3)
    free (pwb->integral_over_r);
  } // end of last term in the L1 loop

  /* Free 3j values array */
  for (thread=0; thread < number_of_threads; ++thread) {
    for (int ii=0; ii < n_geometrical_factors; ++ii)
      free (value[thread][ii]);
    free (value[thread]);
  }
  free (value);

  return _SUCCESS_;

}