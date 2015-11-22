/** @file spectra2.c Compute spectra of second order pertrubations
 *
 * Christian Fidler, 11.05.2015    
 * 
 *
 */

#include "spectra2.h"

int spectra2_init(
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
     struct lensing * ple,
     struct bispectra * pbi,
     struct spectra * psp
     )
{

  // ====================================================================================
  // =                               Preliminaries checks                                 =
  // ====================================================================================

  /* Check whether we need to compute intrinsic spectra at all */  

  if (!ppt2->has_cls && !ppt2->has_pks) {

    printf_log_if (psp->spectra_verbose, 0,
      "No second-order spectrum requested; spectra2 module skipped.\n");

    return _SUCCESS_;
  }
  else {
    
    printf_log_if (psp->spectra_verbose, 0,
      "Computing second-order spectra.\n");
  }



  // ====================================================================================
  // =                                  Preparations                                    =
  // ====================================================================================

  /* The second-order C_l and P(k) are obtained by solving a convolution integral over
  k1 and k2. We solve this integral by summing over the nodes where the source and
  transfer functions are computed, using the trapezoidal rule. For k1 and k2 the node
  points concide with the smooth sampling in ppt2->k. The C_l need to be integrated
  also along the triangular direction k3; in this direction we use the fine sampling
  of the transfer2 module, ptr2->k. */

  psp->k_size = ppt2->k_size;
  psp->k = ppt2->k;




  // ====================================================================================
  // =                               Angular power spectra                              =
  // ====================================================================================

  if (ppt2->has_cls) {

    class_call (spectra2_cls (
                  ppr, ppr2, pba, pth, ppt, ppt2, pbs,
                  pbs2, ptr, ptr2, ppm, ple, pbi, psp),
      psp->error_message,
      psp->error_message);

  }

  
  
  // ==================================================================================
  // =                               Fourier Power Spectra                            =
  // ==================================================================================


  if (ppt2->has_pks) {

    class_call (spectra2_pks (
                  ppr, ppr2, pba, pth, ppt, ppt2, pbs,
                  pbs2, ptr, ptr2, ppm, ple, pbi, psp),
      psp->error_message,
      psp->error_message);

  }

  return _SUCCESS_;

}



/**
 * Compute the second-order C_l, that is, the angular power spectrum of the
 * CMB at second-order, for all the required fields (TT, TE, BB...).
 *
 * This function computes the statistical average of <a_lm^(2) a_lm^(2)>,
 * where a_lm^(2) is the second-order multipole decomposition of the
 * CMB sky. Note that this function does not compute the <a_lm^(1) a_lm^(3)>
 * contribution to the spectrum, which has the same order of magnitude of
 * <a_lm^(2) a_lm^(2)>, because SONG cannot compute 3rd-order perturbations.
 * 
 * The second-order C_l are obtained by solving a 3D integral in Fourier
 * space (k1,k2,k3) and are stored in the psp->cl array.
 */

int spectra2_cls (
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
     struct lensing * ple,
     struct bispectra * pbi,
     struct spectra * psp
     )
{

  /* SONG supports only scalar perturbations at first order */
  int index_md = ppt->index_md_scalars;

  class_test (!ppt->has_scalars || ppt->md_size > 1,
    psp->error_message,
    "SONG supports only scalar modes at first order");

  class_test (ppt->ic_size[0] > 1,
    psp->error_message,
    "SONG only supports one type of initial condition");


  /* We compute the C_l for all l where we computed the second-order transfer
  functions; these are determined by the user via the parameters l_linstep and
  l_logstep. Note that if lensing is requested, the maximum l in this sampling
  might be smaller than the one for the first-order C_l (psp->l). */
  psp->l_size_song = ptr2->l_size;
  psp->l_song = ptr2->l;


  /* TEMPORARILY DISABLED */  
  // if ((ppt2->k3_sampling == sym_k3_sampling))
  //   printf("Angular power spectra can only be computed in standart k3 sampling. Change sampling strategy\n");
    

  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct) {

    if (psp->cl_type[index_ct] != second_order)
      continue;

    printf_log_if (psp->spectra_verbose, 0,
      " -> computing second-order Cl_%s\n", psp->ct_labels[index_ct]);

    /* Variables used for the output of the parallel computation */
    int abort = _FALSE_;
    

    // ====================================================================================
    // =                                Sum over M, cycle L                               =
    // ====================================================================================

    #pragma omp parallel for schedule (dynamic)
    for (int index_l=0; index_l<psp->l_size_song; ++index_l) {

      int l = psp->l_song[index_l];

      /* Generate grid for k3 integration. We do it inside the l-loop
      to avoid race conditions */
      double * k3_grid = malloc (ptr2->k3_size_max*sizeof(double));

      /* Initialise spectrum */
      double * result = &psp->cl[index_md][index_l * psp->ct_size + index_ct];
      *result = 0;

      /* Array to store the various m contributions */
      double result_m[ppr2->m_size];

      /* Loop over m */
      for (int index_M=0; index_M < ppr2->m_size; ++index_M) {

        int m = ptr2->m[index_M];
        
        result_m[index_M] = 0;
        
        printf_log_if (psp->spectra_verbose, 1,
          "     * computing spectra for l = %d m = %d\n", l, m);


        /* Determine the transfer functions needed to compute the requested Cl */

        int index_tt_1;
        int index_tt_2;
      
        if (psp->has_t2t2 && index_ct==psp->index_ct_t2t2) {
          index_tt_1 = ptr2->index_tt2_T + lm_cls(index_l, index_M);
          index_tt_2 = ptr2->index_tt2_T + lm_cls(index_l, index_M);
          // printf ("l=%3d[%3d], m=%d[%d]: index_tt = %d\n", l, index_l, m, index_M, index_tt_1);
        }
        else if (psp->has_e2e2 && index_ct==psp->index_ct_e2e2) { 
          index_tt_1 = ptr2->index_tt2_E + lm_cls(index_l, index_M);
          index_tt_2 = ptr2->index_tt2_E + lm_cls(index_l, index_M);
          // printf ("l=%3d[%3d], m=%d[%d]: index_tt = %d\n", l, index_l, m, index_M, index_tt_1);
        }
        else if (psp->has_t2e2 && index_ct==psp->index_ct_t2e2) { 
          index_tt_1 = ptr2->index_tt2_T + lm_cls(index_l, index_M);
          index_tt_2 = ptr2->index_tt2_E + lm_cls(index_l, index_M);
          // printf ("l=%3d[%3d], m=%d[%d]: index_tt_1 = %d, index_tt_2 = %d\n", l, index_l, m, index_M, index_tt_1, index_tt_2);
        }
        else if (psp->has_b2b2 && index_ct==psp->index_ct_b2b2) {
          if (m==0) continue; /* scalars do not contribute to b_modes */
          index_tt_1 = ptr2->index_tt2_B + lm_cls(index_l, index_M);
          index_tt_2 = ptr2->index_tt2_B + lm_cls(index_l, index_M);
          // printf ("l=%3d[%3d], m=%d[%d]: index_tt = %d\n", l, index_l, m, index_M, index_tt_1);
        }
        else {
          sprintf (psp->error_message, "could not find transfer type for requested cl");
          abort = _TRUE_;
        }


        /* Load transfer functions if needed */

        if (ppr2->load_transfers_from_disk || ppr2->store_transfers_to_disk) {

          class_call_parallel (transfer2_load_transfers_from_disk (
                                 ppt2,
                                 ptr2,
                                 index_tt_1),
            ptr2->error_message,
            psp->error_message);

          if (index_tt_2 != index_tt_1)
            class_call_parallel (transfer2_load_transfers_from_disk (
                                   ppt2,
                                   ptr2,
                                   index_tt_2),
              ptr2->error_message,
              psp->error_message);
        }



        // =====================================================================================
        // =                                Integrate over k1,k2,k3                            =
        // =====================================================================================

        for (int index_k1 = 0; index_k1 < psp->k_size; ++index_k1) {

          /* Find stepsize */

          double step_k1;
          if (index_k1 == 0) step_k1 = (psp->k[1] - psp->k[0])/2;
          else if (index_k1 == psp->k_size-1) step_k1 = (psp->k[psp->k_size-1] - psp->k[psp->k_size-2])/2;
          else step_k1 = (psp->k[index_k1+1] - psp->k[index_k1-1])/2;

          double k1 = psp->k[index_k1];

          /* Primordial spectrum k1 */

          double spectra_k1;
          class_call_parallel (primordial_spectrum_at_k (
                                 ppm,
                                 index_md,
                                 linear,
                                 k1,
                                 &spectra_k1),
            ppm->error_message,
            psp->error_message);
          spectra_k1 = 2*_PI_*_PI_/(k1*k1*k1) * spectra_k1;

          /* Loop over k2 */

          for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {

            /* Print some info */
            printf_log_if (psp->spectra_verbose, 2,
              "     \\ computing integral for (k1,k2) = (%.3g,%.3g)\n",
              psp->k[index_k1], psp->k[index_k2]);


            /* Find stepsize */

            double step_k2;
            if (index_k1 == 0) step_k2 = 0;
            else if (index_k2 == 0) step_k2 = (psp->k[1] - psp->k[0])/2;
            else if (index_k2 == index_k1) step_k2 = (psp->k[index_k1] - psp->k[index_k1-1])/2;
            else step_k2 = (psp->k[index_k2+1] - psp->k[index_k2-1])/2;

            double k2 = psp->k[index_k2];


            /* Primordial spectrum k2 */

            double spectra_k2;
            class_call_parallel (primordial_spectrum_at_k (
                                   ppm,
                                   index_md,
                                   linear,
                                   k2,
                                   &spectra_k2),
              ppm->error_message,
              psp->error_message);
            spectra_k2 = 2*_PI_*_PI_/(k2*k2*k2) * spectra_k2;


            /* Integration grid in k3 */

            int dump;
            class_call_parallel (transfer2_get_k3_list (
                                   ppr,
                                   ppr2,
                                   ppt2,
                                   pbs,
                                   pbs2,
                                   ptr2,
                                   index_k1,
                                   index_k2,
                                   k3_grid,
                                   &dump
                                   ),
              ptr2->error_message,
              psp->error_message);

            int k3_size = ptr2->k_size_k1k2[index_k1][index_k2];

            class_test_parallel (k3_size < 2,
              psp->error_message,
              "integration grid has less than two elements, cannot use trapezoidal integration");


            /* Define the pointer to the second-order transfer functions as a function of k3.
            Note that this transfer function has already been rescaled according to eq. 6.26
            of http://arxiv.org/abs/1405.2280 in the perturbations.c module.  */

            double * transfer_1 = ptr2->transfer[index_tt_1][index_k1][index_k2];
            double * transfer_2 = ptr2->transfer[index_tt_2][index_k1][index_k2];

            int triangular_first = 0;
            int triangular_last = 0;

            /* Final integration over k3 */

            for (int index_k3 = 0; index_k3 < k3_size; ++index_k3) {

              double k3 = k3_grid[index_k3];

              /* Find stepsize including the triangular condition edge */

              double step_k3;
              if (index_k3 == 0) step_k3 = (k3_grid[1] - k3_grid[0])/2.;
              else if (index_k3 == k3_size -1) step_k3 = (k3_grid[k3_size-1] - k3_grid[k3_size -2])/2.;
              else step_k3 = (k3_grid[index_k3+1] - k3_grid[index_k3 -1])/2.;

              /* Triangular condition */

              if ( k3 < k1-k2 ||  k3 > k1+k2) { // out of triangular region
                step_k3 = 0.;
              }
              if ( k3 >= k1-k2 && triangular_first == 0 && index_k3 < k3_size-1 ) { //first point
                triangular_first = 1;
                step_k3 = (k3_grid[index_k3+1]-k3)/2. + k3 - (k1-k2);
              }
              if ( k3_grid[index_k3+1] > k1+k2 && triangular_last == 0 && index_k3 > 0) {
                //last point, note that this may not be found if the last point is bigger than kmax,
                // however in that case no special treatment for the last point is needed.
                triangular_last = 1;
                step_k3 = (k3 - k3_grid[index_k3-1])/2. + (k1+k2) - k3;
              }
              
              /* If we are in the extrapolated range, skip this k3 */
              if (step_k3 == 0)
                continue;

              /* Counter the rescaling of the sources, if needed; see documentation for
              ppt2->rescale_cmb_sources for more details. */
              double inverse_rescaling = 1;
              class_call_parallel (perturb2_cmb_rescaling (ppt2,k1,k2,k3,m,&inverse_rescaling),
                ppt2->error_message, psp->error_message);
              inverse_rescaling = 1/inverse_rescaling;


              // =====================================================================================
              // =                                Add up contributions                               =
              // =====================================================================================

              double total_factor =
                  // symmetry factor (doing only half plane in k1 k2)
                  2 *
                  // integration weight
                  k1 * k2 * k3 *
                  // integration stepsize
                  step_k1 * step_k2 * step_k3 *
                  // if the sources were rescaled, revert the rescaling
                  inverse_rescaling * 
                  // spectra factors (2l+1) already in transfer def (how does sum over m work? does this need a factor 2 for m neq 0? yes :))
                  2/_PI_ / pow(2*_PI_,3) *
                  // definition of second order perturbation theory
                  // already in transfer definition 1./2./2.*
                  // factor 4 of Delta also already in transfer def
                  // primordial spectra
                  spectra_k1 * spectra_k2;

              double integrand = total_factor * transfer_1[index_k3] * transfer_2[index_k3];

              *result += integrand;
              result_m[index_M] += integrand;

              /* Debug: print the integrand function */
              // if (l==2 && m==0)
              //   printf ("integrand = %g, result =%g\n", integrand, *result);

            } // loop k3
          } // loop k2
        } // loop k1

        if (ppr2->load_transfers_from_disk || ppr2->store_transfers_to_disk) {

          class_call_parallel (transfer2_free_type_level (
                                 ppt2,
                                 ptr2,
                                 index_tt_1),
            ptr2->error_message,
            psp->error_message);

          if (index_tt_2 != index_tt_1)
            class_call_parallel (transfer2_free_type_level (
                                   ppt2,
                                   ptr2,
                                   index_tt_2),
              ptr2->error_message,
              psp->error_message);
        }

        /* Debug: print the M contribution */
        // if (m==0 || m==1)
        //   printf ("%10d %14g %14g\n", l, result_m[index_M], *result);

      } // sum over M
    
      free (k3_grid);
      
      #pragma omp flush(abort)
      
    } // loop over l

    if (abort)
      return _FAILURE_;

  } // loop over ct


  // =====================================================================================
  // =                                    Interpolation                                  =
  // =====================================================================================

  /* Compute Cl derivatives in view of spline interpolation */
  class_call (spectra_cls_spline(
                pba,
                ppt,
                ptr,
                ppm,
                psp,
                index_md),
    psp->error_message,
    psp->error_message);


  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct) {

    if (psp->cl_type[index_ct] != second_order)
      continue;

    /* Make sure that the interpolating the C_l for l>l_max_song yields zer0 */
    psp->l_max_ct[index_md][index_ct] = psp->l_song[psp->l_size_song-1];

    /* Make sure that the C_l and their derivatives vanish for l > l_max */
    for (int index_l=0; index_l<psp->l_size[index_md]; ++index_l) {
      int l = psp->l[index_l];
      if (l > psp->l_song[psp->l_size_song-1]) {
        psp->cl[index_md][index_l*psp->ct_size + index_ct] = 0;
        psp->ddcl[index_md][index_l*psp->ct_size + index_ct] = 0;
        if (psp->compute_cl_derivative == _TRUE_) {
          psp->d_lsq_cl[index_md][index_l*psp->ct_size + index_ct] = 0;
          psp->dd_lsq_cl[index_md][index_l*psp->ct_size + index_ct] = 0;
          psp->spline_d_lsq_cl[index_md][index_l*psp->ct_size + index_ct] = 0;
        }
      }
    }
  }

  /* Debug: test the interpolation */
  // for (int l=2; l <= psp->l_max[index_md]; ++l) {
  //
  //   double * cl = calloc (psp->ct_size, sizeof (double));
  //   double ** cl_md = malloc (psp->md_size * sizeof(double*));
  //   double ** cl_md_ic = malloc (ppt->md_size * sizeof(double*));
  //
  //   for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {
  //     class_alloc(cl_md[index_mode], psp->ct_size*sizeof(double), psp->error_message);
  //     class_alloc(cl_md_ic[index_mode], psp->ic_ic_size[index_mode]*psp->ct_size*sizeof(double), psp->error_message);
  //   }
  //
  //   class_call (spectra_cl_at_l(
  //                 psp,
  //                 (double)l,
  //                 cl,
  //                 cl_md,
  //                 cl_md_ic),
  //     psp->error_message,
  //     psp->error_message);
  //
  //   int index_ct = psp->index_ct_t2t2;
  //
  //   int index_L = -1;
  //   for (int index_l=0; index_l < psp->l_size[index_md]; ++index_l) {
  //     if ((int)(psp->l[index_l]+_EPS_) == l) {
  //       index_L = index_l;
  //       break;
  //     }
  //   }
  //
  //   double cl_node = index_L<0 ? 0:psp->cl[index_md][index_L * psp->ct_size + index_ct];
  //   double ddcl_node = index_L<0 ? 0:psp->ddcl[index_md][index_L * psp->ct_size + index_ct];
  //
  //   fprintf (stderr, "%10d %16g %16g %16g\n", l, cl[psp->index_ct_t2t2], cl_node, ddcl_node);
  //
  // }
  
  return _SUCCESS_;
  
}



/**
 * Compute the second-order P(k), that is, the Fourier power spectrum at second
 * order, for all the required perturbations.
 *
 * For a given perturbation X, this function computes <X(k)^(2) X(k)^(2)>. The
 * same-order contribution <X(k)^(1) X(k)^(3)> is not computed, because it 
 * requires the computation of third-order perturbations.
 * 
 * The second-order P(k) are obtained by solving a 2D convolution integral in
 * Fourier space, and are stored in the psp->ln_pk array.
 */

int spectra2_pks (
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
     struct lensing * ple,
     struct bispectra * pbi,
     struct spectra * psp
     )
{
  
  
  // ====================================================================================
  // =                              Store primordial P(k)                               =
  // ====================================================================================
  
  /* Store the primordial power spectrum inside the spectra structure for faster access */

  class_call (spectra_primordial_power_spectrum(
                pba,
                ppt,
                ptr,
                ppm,
                psp),
    psp->error_message,
    pbi->error_message);



  // ====================================================================================
  // =                                  Compute P(k)                                    =
  // ====================================================================================

  class_call (spectra2_integrate_fourier (
                ppr,
                ppr2,
                ppt,
                ppm,
                ppt2,
                psp),
    psp->error_message,
    psp->error_message);




//    class_call(spectra2_get_k3_size (
//       ppr,
//       ppr2,
//       ppt2,
//       psp
//       ),
//     psp->error_message,
//     psp->error_message);
//
//
//
//   /* Shortcut to the file where we shall print the transfer functions. */
//   FILE * file_sp = psp->spectra_file;
//
//   /* Choose how label & values should be formatted */
//   char format_label[64] = "%13s(%02d) ";
//   char format_value[64] = "%+17e ";
//   char buffer[64];
//
//   // conformal time
//   int index_sp = 1;
//
//   // Print time
//   /*
//   fprintf(file_sp, format_label, "tau",index_sp);
//   index_sp++;
//
//   fprintf(file_sp, format_label, "a", index_sp);
//   index_sp++;
//
//   fprintf(file_sp, format_label, "Y",index_sp);
//   index_sp++;
//
//   for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {
//       sprintf(buffer, "magnet k=%f", psp->k[index_k3] );
//        fprintf(file_sp, format_label, buffer,index_sp );
//        index_sp++;
//      }
//
//   fprintf(file_sp, "\n");
//   */
//
//   //Print k
//
//   fprintf(file_sp, format_label, "k",index_sp);
//   index_sp++;
//
//   for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
//       sprintf(buffer, "magnet tau=%f", ppt2->tau_sampling[index_tau] );
//        fprintf(file_sp, format_label, buffer,index_sp );
//        index_sp++;
//      }
//
//      fprintf(file_sp, "\n");
//
//
//   /* Four loops over k1, k2, transfer type follow */
//
//
//   if ((ppt2->k3_sampling == sym_k3_sampling) ) {
//
//     class_call(
//       spectra2_integrate_fourier_sym(
//           ppr,
//           ppr2,
//           ppt,
//           ppm,
//           ppt2,
//           psp),
//         ptr2->error_message,
//         psp->error_message);
//
//
//   } else {
//     class_call(
//       spectra2_integrate_fourier(
//           ppr,
//           ppr2,
//           ppt,
//           ppm,
//           ppt2,
//           psp),
//         ptr2->error_message,
//         psp->error_message);
//
//   }
//
//
//
//
//
// //  for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {
//  //   fprintf(file_sp, format_value, psp->k[index_k3]);
//  //    fprintf(file_sp, format_value, psp->spectra[ppt2->index_tp2_M + lm(1,1)][index_k3*ppt2->tau_size + 300] );
//  //   fprintf(file_sp, "\n");
// //  }
//   int last_index = 0;
//
//   double *pvecback;
//   class_alloc (pvecback, pba->bg_size*sizeof(double), psp->error_message);
//
//  // print time
//  /*
//   for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
//
//   class_call (background_at_tau(
//                 pba,
//                 ppt2->tau_sampling[index_tau],
//                 pba->long_info,
//                 pba->inter_normal,
//                 &last_index,
//                 pvecback),
//     pba->error_message,
//     psp->error_message);
//
//       double a = pvecback[pba->index_bg_a];
//       double H = pvecback[pba->index_bg_H];
//       double Hc = pvecback[pba->index_bg_H]*a;
//       double Y = a/pba->a_eq;
//
//
//     fprintf(file_sp, format_value, ppt2->tau_sampling[index_tau]);
//
//     fprintf(file_sp, format_value, a);
//
//     fprintf(file_sp, format_value, Y);
//
//     for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {
//
//        fprintf(file_sp, format_value, psp->spectra[ppt2->index_tp2_M + lm(0,0)][index_k3*ppt2->tau_size + index_tau] );
//
//      }
//     fprintf(file_sp, "\n");
//   }
//   */
//   // Print k
//
//   for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {
//
//
//     fprintf(file_sp, format_value, psp->k[index_k3]);
//
//     for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
//
//        fprintf(file_sp, format_value, psp->spectra[ppt2->index_tp2_M + lm(0,0)][index_k3*ppt2->tau_size + index_tau] );
//
//      }
//     fprintf(file_sp, "\n");
//   }
//
//    printf("keq = %f \n",pba->k_eq);
//
//   fclose(psp->spectra_file);


  return _SUCCESS_;

}





int spectra2_integrate_fourier (
			struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct primordial * ppm,
      struct perturbs2 * ppt2,
      struct spectra * psp
      )
{

  /* SONG supports only scalar perturbations at first order */
  int index_md = ppt->index_md_scalars;

  class_test (!ppt->has_scalars || ppt->md_size > 1,
    psp->error_message,
    "SONG supports only scalar modes at first order");

  class_test (ppt->ic_size[0] > 1,
    psp->error_message,
    "SONG only supports one type of initial condition");

  /* Result vector; note that we compute the power spectra in psp->k */
  double (*pk)[psp->k_size*ppt2->tau_size];
  class_calloc (pk, psp->pk_size*psp->k_size*ppt2->tau_size, sizeof (double), psp->error_message);


  // ====================================================================================
  // =                             Find triangular region                               =
  // ====================================================================================

  int (*index_k_triangular_min)[psp->k_size] = calloc (psp->k_size*psp->k_size, sizeof(int));
  int (*index_k_triangular_max)[psp->k_size] = calloc (psp->k_size*psp->k_size, sizeof(int));
  int (*k_triangular_size)[psp->k_size] = calloc (psp->k_size*psp->k_size, sizeof(int));

  for (int index_k1=0; index_k1 < psp->k_size; ++index_k1) {

    double k1 = psp->k[index_k1];

    for (int index_k2=0; index_k2 < psp->k_size; ++index_k2) {

      double k2 = psp->k[index_k2];

      /* Limits imposed on k3 by the triangular condition */
      double k_min = fabs(k1-k2);
      double k_max = k1+k2;

      /* Find the index corresponding to k_min inside psp->k */
      int index_k_min = 0;
      while (psp->k[index_k_min] < k_min && index_k_min < psp->k_size-1)
        ++index_k_min;

      /* Find the index corresponding to k_max inside psp->k */
      int index_k_max = psp->k_size-1;
      while (psp->k[index_k_max] > k_max && index_k_max > 0)
        --index_k_max;    

      index_k_triangular_min[index_k1][index_k2] = index_k_min;
      index_k_triangular_max[index_k1][index_k2] = index_k_max;
      k_triangular_size[index_k1][index_k2] = index_k_max - index_k_min + 1;

      /* Each (k1,k2) pair will have at least one node in its range, given by
      the largest between k1 and k2 */
      class_test ( k_triangular_size[index_k1][index_k2] < 0,
        psp->error_message,
        "k1=%g[%d], k2=%g[%d], found index_k_min=%d, index_k_max=%d\n",
        k1, index_k1, k2, index_k2, index_k_min, index_k_max);

    } // for(index_k2)
  } // for(index_k1)



  // =====================================================================================
  // =                                 Cycle P(k) type                                   =
  // =====================================================================================

  for (int index_pk=0; index_pk<psp->pk_size; ++index_pk) {
    
    if (psp->pk_type[index_pk] != second_order)
      continue;

    if (psp->spectra_verbose > 0)
      printf (" -> computing second-order P(k) %s\n", psp->pk_labels[index_pk]);

    /* Determine the source functions needed to compute the requested Pk */
    int index_tp2;
  
    if (psp->has_pk_delta2_delta2_cdm && index_pk==psp->index_pk_delta2_delta2_cdm) {
      index_tp2 = ppt2->index_tp2_delta_cdm;
    }
    else if (psp->has_pk_magnetic && index_pk==psp->index_pk_magnetic) {
      index_tp2 = ppt2->index_tp2_M;
    }
    else {
      class_stop (psp->error_message, "could not find source type for requested Pk");
    }


    // =====================================================================================
    // =                                   Sum over k1                                     =
    // =====================================================================================

    int abort = _FALSE_;

    #pragma omp parallel for schedule(dynamic)
    for (int index_k1 = 0; index_k1 < psp->k_size; ++index_k1) {

      double k1 = psp->k[index_k1];

      /* Primordial spectrum of the curvature potential phi in k1 */
      double spectra_k1 = psp->pk_pt[index_k1];
      
      if (psp->spectra_verbose > 1)
        printf ("     * computing integral for index_k1=%d of %d, k1=%g\n",
        index_k1, ppt2->k_size, k1);

      if (ppr2->load_sources_from_disk || ppr2->store_sources_to_disk)
        class_call_parallel (perturb2_load_sources_from_disk(ppt2, index_k1),
                               ppt2->error_message,
                               psp->error_message);

      /* Find stepsize */
      double step_k1;
      if (index_k1 == 0) step_k1 = (psp->k[1] - psp->k[0])/2;
      else if (index_k1 == psp->k_size-1) step_k1 = (psp->k[psp->k_size-1] - psp->k[psp->k_size-2])/2;
      else step_k1 = (psp->k[index_k1+1] - psp->k[index_k1-1])/2;



      // =====================================================================================
      // =                                 Cycle over k3                                     =
      // =====================================================================================

      for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {

        double k3 = psp->k[index_k3];


        // --------------------------------------------------------------------------------
        // -                         Determine integration range                          -
        // --------------------------------------------------------------------------------

        /* The lower limit in the k2 integration range is determined by the smallest
        k where we computed the source function and by the smallest value allowed
        by the triangular condition. */
        double k2_min = MAX (psp->k[0], fabs (k1-k3));

        /* Similarly, the upper limit is determined by the triangular condition and by
        the highest k where we computed the source function. */
        double k2_max = MIN (psp->k[psp->k_size-1], k1+k3);
        
        /* We shall compute the integral by means of a weighted sum over the k2 nodes;
        here we determine which nodes ought to be considered, based on the conditions
        imposed above on k2_min and k2_max */
        int index_k2_min = 0;
        while (psp->k[index_k2_min] < k2_min && index_k2_min < psp->k_size-1)
          ++index_k2_min;

        int index_k2_max = psp->k_size-1;
        while (psp->k[index_k2_max] > k2_max && index_k2_max > 0)
          --index_k2_max;

        int k2_size = index_k2_max - index_k2_min + 1;

        /* Uncomment to use the precomputed indices (same result) */
        // int index_k2_min = index_k_triangular_min[index_k1][index_k3];
        // int index_k2_max = MIN (index_k1, index_k_triangular_max[index_k1][index_k3]);
        // int k2_size = index_k2_max - index_k2_min + 1;

        printf_log_if (psp->spectra_verbose, 2,
          "     \\ computing integral for k1=%g[%d], k3=%g[%d]), k2_size=%d\n",
          k1, index_k1, k3, index_k3, k2_size);

        /* Check that the integration range is positive */
        class_test_parallel (k2_max-k2_min<-ppr->smallest_allowed_variation && k2_size>0,
          psp->error_message,
          "Tound negative integration range: k2_max-k2_min=%27.17g",
          k2_max-k2_min);



        // =====================================================================================
        // =                                   Sum over k2                                     =
        // =====================================================================================

        for (int index_k2 = index_k2_min; index_k2 <= index_k2_max; ++index_k2) {

          double k2 = psp->k[index_k2];
          
          /* Primordial spectrum of the curvature potential phi in k2 */
          double spectra_k2 = psp->pk_pt[index_k2];

          /* The indexing in this loop ensures that the triangular condition is respected */
          class_test_parallel (k1+k2<k3 || fabs(k1-k2)>k3,
            psp->error_message,
            "triangular condition should be respected at this point!");


          // --------------------------------------------------------------------------------
          // -                               Triangular step                                -
          // --------------------------------------------------------------------------------

          /* Find the correct stepsize for k2 based on the triangular inequality */

          double step_k2;

          /* If there is only one node in the k2 integration range, we adopt
          the rectangular integration rule between the lower limit (k2_min)
          and the upper limit (k2_max). Note that it is possible for index_k1=0
          that k2_min is larger than k2_max, due to the artificial upper limit
          on k2, hence the MAX. */
          if (k2_size == 1)
            step_k2 = k2_max - k2_min;

          /* If we are close to the left border of the triangular regime:
          We take the regular trapezoidal step (including the 1/2 factor) up
          to the first node, then we extend the step to the left all the way
          to the end of triangular region, up to k2_min = k1-k3, without the
          1/2 factor because this node is the only one covering this extra
          region. */
          else if (index_k2 == index_k2_min)
            step_k2 = (psp->k[index_k2_min+1] - psp->k[index_k2_min])/2 + (psp->k[index_k2_min] - k2_min);

          /* If we are close to the right border of the triangular regime:
          We take the regular trapezoidal step (including the 1/2 factor) up
          to the last node, then we extend the step all the way to the end of
          triangular region, up to k2_max = MIN(k1,k1+k3), without the 1/2
          factor because this node is the only one covering this extra
          region. */
          else if (index_k2 == index_k2_max)
            step_k2 = (psp->k[index_k2_max] - psp->k[index_k2_max-1])/2 + (k2_max - psp->k[index_k2_max]);

          /* If we are not adjacent to the triangular region, we take the
          regular trapezoidal step */
          else
            step_k2 = (psp->k[index_k2+1] - psp->k[index_k2-1])/2;

          
          // =====================================================================================
          // =                                 Cycle over time                                   =
          // =====================================================================================

          for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {


            // ---------------------------------------------------------------------------------
            // -                                Interpolate in k3                              -
            // ---------------------------------------------------------------------------------

            double source;

            class_call_parallel (perturb2_sources_at_k3 (
                                   ppr2,
                                   ppt,
                                   ppt2,
                                   index_tp2,
                                   index_k1,
                                   index_k2,
                                   k3,
                                   -1,
                                   index_tau,
                                   _TRUE_, /* extrapolate in some cases */
                                   _TRUE_, /* give us the vanilla sources */
                                   &source),
              ppt2->error_message,
              psp->error_message);

            
            // ---------------------------------------------------------------------------------
            // -                               Increment integral                              -
            // ---------------------------------------------------------------------------------

            double integrand = 
              // integration weight
               k1*k2/k3 *
              // integration stepsize
               step_k1 * step_k2 *
              // phi angleintegration
               2 * _PI_ *
              // spectra factors
               2 / pow(2*_PI_,3) *
              // definition of second order perturbation theory
               1./2./2.*
              // primordial spectra
               spectra_k1 * spectra_k2 *
              // sources
               source * source
              // 6.652 * 1.e-29 is sigma_t in m^2, 1.878 * 1.e-26 is critical density in kg/m^3
              // 1.602 * 1.e-19 is electron charge ,  1.e-6 is e-10 from CLASS convention
              // (rho_g_CLASS^0 = e+10 h^2/c^2 * Omega_g^0) and e+4 to covert from international units to Gauss
              // (1 G = 10 e-4 kg/ s/ C)
              //*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
              //*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
              ;

            #pragma omp atomic
            pk[index_pk][index_tau*psp->k_size + index_k3] += integrand;

          } // for(index_tau)

        } // for(index_k3)

      } // for(index_k2)

      if (ppr2->load_sources_from_disk || ppr2->store_sources_to_disk)
        class_call_parallel (perturb2_free_k1_level (ppt2, index_k1),
          ppt2->error_message, ppt2->error_message);

    #pragma omp flush(abort)

    } // for(index_k1)
    
    if (abort == _TRUE_)
      return _FAILURE_;

  } // for(index_pk)



  // ====================================================================================
  // =                                   Output P(k)                                    =
  // ====================================================================================

  FILE * output_file = stderr;
  int index_tau_output = ppt2->tau_size - 1;
  

  /* Print labels */
  fprintf (output_file, "%22s ", "k");
  for (int index_pk=0; index_pk < psp->pk_size; ++index_pk) {
    if (psp->pk_type[index_pk] != second_order) continue;
    fprintf (output_file, "%22s ", psp->pk_labels[index_pk]);
  }
  fprintf (output_file, "\n");
    
    
  /* Print data */
  for (int index_k=0; index_k < psp->k_size; ++index_k) {
    fprintf (output_file, "%22.10g ", psp->k[index_k]);
    for (int index_pk=0; index_pk < psp->pk_size; ++index_pk) {
      if (psp->pk_type[index_pk] != second_order) continue;
      fprintf (output_file, "%22.10g ", pk[index_pk][index_tau_output*psp->k_size + index_k]);
    }
    fprintf (output_file, "\n");
  }
    
  free (index_k_triangular_min);
  free (index_k_triangular_max);
  free (k_triangular_size);
  free (pk);
  
  
  
  // ====================================================================================
  // =                                   Interpolate P(k)                               =
  // ====================================================================================

  /* So far we have computed the second-order power spectra in ppt2->k. CLASS however
  outputs P(k) in ppt->k. Here we interpolate SONG power spectrum in CLASS output
  values. */

	return _SUCCESS_;
	
}




int spectra2_integrate_fourier_sym(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct primordial * ppm,
      struct perturbs2 * ppt2,
      struct spectra * psp
      )
{

//   /* Vector that will contain the interpolated sources for a given k1,k2 and source type */
//   double ** interpolated_sources_in_k;
//   class_alloc (interpolated_sources_in_k, ppt2->tp2_size*sizeof(double *), psp->error_message);
//
//   /* Spline coefficients for the source interpolation */
//   double ** sources_k_spline;
//   class_alloc (sources_k_spline, ppt2->tp2_size*sizeof(double *), psp->error_message);
//
//   double step_k1,step_k3;
//   // loop over k1t (here we call k1t the transformed k1, which is indexed by index_k1)
//   for (int index_k1 = 0; index_k1 < psp->k_size; ++index_k1) {
//     if (psp->spectra_verbose > 1)
//       printf (" computing integral for index_k1=%d of %d, k1t=%g\n",
//         index_k1, psp->k_size, psp->k[index_k1]);
//
//     // load sources
//     if ((ppr2->load_sources_from_disk == _TRUE_) || (ppr2->store_sources_to_disk == _TRUE_))
//           class_call(perturb2_load_sources_from_disk(ppt2, index_k1),
//             ppt2->error_message,
//             psp->error_message);
//
//
//
//     double k1t = psp->k[index_k1];
//
//     //loop over k3t
//     for (int index_k3 = 0; index_k3 < psp->k_size; ++index_k3) {
//       if (psp->spectra_verbose > 2)
//         printf(" -> computing integral for (k1t,k3t) = (%.3g,%.3g)\n", psp->k[index_k1], psp->k[index_k3]);
//
//       double step_k3;
//
//       if (index_k3 == 0) step_k3 = (psp->k[1] - psp->k[0])/2.;
//       else if (index_k3 == psp->k_size -1) step_k3 = (psp->k[psp->k_size-1] - psp->k[psp->k_size -2])/2.;
//       else step_k3 = (psp->k[index_k3+1] - psp->k[index_k3 -1])/2.;
//
//       double k3t = psp->k[index_k3];
//
//
//
//       // -----------------------------------------------------------------------
//       // -                      Interpolate sources in k                       -
//       // -----------------------------------------------------------------------
//
//       for (int index_tp=ppt2->index_tp2_M + lm(0,0); index_tp<=ppt2->index_tp2_M + lm(1,1); ++index_tp) {
//
//       // find interpolation for the specific case of k1t,k3t,tp. This is an interpolation in k2t
//
//       class_alloc(
//           interpolated_sources_in_k[index_tp],
//           psp->k_size*ppt2->tau_size*sizeof(double),
//           psp->error_message);
//
//        class_alloc(
//           sources_k_spline[index_tp],
//           ppt2->k_size*ppt2->tau_size*sizeof(double),
//           psp->error_message);
//
//        /* we use the same grid for k that is also used for k1t*/
//
//         class_call (spectra2_interpolate_sources_in_k2(
//                       ppr,
//                       ppr2,
//                       ppt,
//                       ppt2,
//                       psp,
//                       index_k1,
//                       index_k3,
//                       index_tp,
//                       psp->k,                   /*this is actually not used properly in our interpolation, clean up!*/    /* All the k_grid are filled */
//                       sources_k_spline[index_tp],
//                       interpolated_sources_in_k[index_tp]   /* Will be filled with interpolated values */
//                       ),
//           psp->error_message,
//           psp->error_message);
//
//
//         // here we loop over k3, not k2t. This is output.
//         /* We now have a position in k1t,k3t and an interpolation in k2t. This translates into
//         knowing k1,k2 and have an interpolation in k3. So we are ready for integration*/
//         for (int index_k = 0; index_k < psp->k_size; ++index_k) {
//           int print = 0;
//
//           // construct the non transformed variables
//           double k = psp->k[index_k];
//           double k1 =  (-k1t + k3t + 2.*k)/2.;
//           double k2 = (k1t + k3t)/2.;
//
//
//
//           double spectra_k1;
//
//           if (psp->spectra_verbose > 2)
//             printf(" -> analysing k = %f, k1 = %f, k2 = %f \n", k,k1,k2);
//
//
//            // using the transformed variables k1 or k2 might be negative. This removes that contribution.
//           if (k1 > 0) {
//             class_call (primordial_spectrum_at_k(
//                   ppm,
//                   ppt->index_md_scalars,
//                   linear,
//                   k1,
//                   &spectra_k1),
//               ppm->error_message,
//               psp->error_message);
//
//             /* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
//               spectra_k1 = 2*_PI_*_PI_/(k1*k1*k1) * spectra_k1;
//           } else spectra_k1 = 0.;
//
//           double spectra_k2;
//
//           if (k2 >0 ) {
//             class_call (primordial_spectrum_at_k (
//                   ppm,
//                   ppt->index_md_scalars,
//                   linear,
//                   k2,
//                   &spectra_k2),
//               ppm->error_message,
//               psp->error_message);
//
//             /* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
//               spectra_k2 = 2*_PI_*_PI_/(k2*k2*k2) * spectra_k2;
//           } else spectra_k2 = 0.;
//
//           // find the correct stepsize for k1t based on the triangular inequality //
//           double triangular_step_k1;
//           double max;
//           /*the allowed region ranges from 0 to 2k. Here we use the symmetry and restrict k1t from k to 2k. In general it would be
//           better to use the full range but in that case the interpolation has to be reworked. 0..k has the best grid,
//           but the worst k2t interpolation (k..2k) . In any case both ways are asymmetric, either have a good interpoaltion
//           or a good grid for integration. Therefore best go the full way. This does not matter on a linear grid.*/
//
//           if (k1t < k) {triangular_step_k1 = 0.;}
//           else if (k1t >=k && k1t <= 2.*k) {
//             /*keep in mind that we do not want to go beyond kmax. Therfore we check if 2k
//             is larger and then use kmax instead for upper limit*/
//
//             if (2.*k > ppt2->k[ppt2->k_size -1])
//                max = ppt2->k[ppt2->k_size -1];
//             else max = 2.*k;
//             if (( index_k1 == 0 || ppt2->k[index_k1 - 1] < k) && ( index_k1 == ppt2->k_size - 1 || ppt2->k[index_k1 + 1] > 2.*k)) {
//               /*if a single point is both first and last*/
//               triangular_step_k1 = (  max - ppt2->k[index_k1])  +  (ppt2->k[index_k1] - k);
//
//              }
//
//             else if ( index_k1 == 0 || ppt2->k[index_k1 - 1] < k) {
//             /*first legal point*/
//               triangular_step_k1 = (ppt2->k[index_k1] - k) + (ppt2->k[index_k1+1]- ppt2->k[index_k1])/2.;
//             }
//             else if ( index_k1 == ppt2->k_size - 1 || ppt2->k[index_k1 + 1] > 2.*k) {
//             /*last legal point*/
//
//               triangular_step_k1 = ( max - ppt2->k[index_k1]) + (ppt2->k[index_k1]- ppt2->k[index_k1-1])/2.;
//             }
//             else{
//               triangular_step_k1 = (ppt2->k[index_k1+1]- ppt2->k[index_k1-1])/2.;
//              }
//           }
//           else triangular_step_k1 = 0.;
//
//           //loop over time (no integration)
//           for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
//             psp->spectra[index_tp][index_k*ppt2->tau_size + index_tau] +=
//           // (index_k1==53 ? triangular_step_k2: 0.)
//           // symmetry factor (doing only half plane in k1 k3)
//            2.*
//           // integration weight
//            k1*k2/k/2.*
//           // integration stepsize
//            step_k3*triangular_step_k1*
//           // phi integration
//            2 * _PI_*
//            // m = +-
//            // nothing here, we simplify compute the fourier power spectrum off all tp2
//            // then the user has to do the sum over m himself if he so desires.
//           // spectra factors
//            2./2./_PI_/2./_PI_/2./_PI_*
//           // definition of second order perturbation theory
//            1./2./2.*
//           // primordial spectra
//            spectra_k1*spectra_k2*
//           //sources
//            interpolated_sources_in_k[index_tp][index_k*ppt2->tau_size + index_tau]*interpolated_sources_in_k[index_tp][index_k*ppt2->tau_size + index_tau]
//           // 6.652 * 1.e-29 is sigma_t in m^2, 1.878 * 1.e-26 is critical density in kg/m^3
//           // 1.602 * 1.e-19 is electron charge ,  1.e-6 is e-10 from CLASS convention
//           // (rho_g_CLASS^0 = e+10 h^2/c^2 * Omega_g^0) and e+4 to covert from international units to Gauss
//           // (1 G = 10 e-4 kg/ s/ C)
//           //*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
//           //*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
//           ;
//           }
//
//
//         }
//
//       } // end of for (index_tp)
//
//
//
//       /* Free the memory for the interpolated sources */
//       for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
//         free(interpolated_sources_in_k[index_tp]);
//         free(sources_k_spline[index_tp]);
//       }
//
//     } // end of for(index_k2)
//
//
//     class_call (perturb2_free_k1_level (ppt2, index_k1), ppt2->error_message, ppt2->error_message);
//
//
//
//   } // end of for(index_k1)
//
//
//   free (interpolated_sources_in_k);
//   free(sources_k_spline);
  return _SUCCESS_;

}



int spectra2_interpolate_sources_in_k2(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct spectra * psp,
      int index_k1,
      int index_k3,
      int index_tp2,
      double * k_grid,
      double * sources_k_spline,
      double * interpolated_sources_in_k
      )
{
//
//
//
//
//
//   int index_tau;
//
//
//
//   int k_pt_size = ppt2->k_size;
//   double * k_pt = ppt2->k;
//
//   double k1 = ppt2->k[index_k1];
//
//
//
//   /*We Have to rewrite the sources in order to allow for interpolation.*/
//
//    double * rsources;
//
//
//   class_calloc (
//       rsources,
//         ppt2->tau_size*(index_k1+1),
//         sizeof(double),
//         psp->error_message);
//
//   double * derivs;
//
//
//   class_calloc (
//       derivs,
//         ppt2->tau_size*(index_k1+1),
//         sizeof(double),
//         psp->error_message);
//
//
//   double * subk;
//
//   class_calloc (
//       subk,
//         (index_k1+1),
//         sizeof(double),
//         psp->error_message);
//
//   for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2){
//     subk[index_k2] = ppt2->k[index_k2];
//     for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++){
//         rsources[index_tau*(index_k1+1) + index_k2] = sources(index_tau,index_k3);
//     }
//   }
//
//
//   if (ppr2->sources_k3_interpolation == cubic_interpolation && index_k1> 3) {
//
//
//
//     class_call (array_spline_table_columns (
//                   subk,
//                   (index_k1+1),
//                   rsources,
//                   ppt2->tau_size,
//                   derivs,
//                   _SPLINE_EST_DERIV_,
//                   psp->error_message),
//          psp->error_message,
//          psp->error_message);
//   }
//
//
//   free(subk);
//   // =======================================================
//   // =                    Interpolation                    =
//   // =======================================================
//
//
//   /* Interpolate at each needed k2t value corresponding to one k3 and k1t value using the usual spline interpolation algorithm */
//    int index_k = 0;
//   double h = k_pt[index_k+1] - k_pt[index_k];
//
//   int index_k_sp;
//
//    /*this index represents k3*/
//   for (index_k_sp = 0; index_k_sp < ppt2->k_size; ++index_k_sp) {
//
//     double k2t = - ppt2->k[index_k1] + 2.* ppt2->k[index_k_sp];
//     index_k = 0;
//     if (k2t < 0 || k2t > k1 || k2t > ppt2->k[ppt2->k_size-1]) {
//
//       interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 0.; /*no extrapolation out of triangular, might be changed to flat extrapolation for triangular (NOT FOR KMAX !!!)*/
//     } else {
//
//       while ((k_pt[index_k+1] < k2t) && (index_k + 1 < ppt2->k_size) ){
//         index_k++;
//         h = k_pt[index_k+1] - k_pt[index_k];
//       }
//
//       class_test(h==0., psp->error_message, "stop to avoid division by zero");
//
//       double b = (k2t- k_pt[index_k])/h;
//
//       double a = 1.-b;
//
//
//     /* Interpolate for each value of conformal time */
//       if (ppr2->sources_k3_interpolation == linear_interpolation ||  index_k1 <= 3) {
//         for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
//           interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] =
//             a * rsources[index_tau*(index_k1+1) + index_k] + b * rsources[index_tau*(index_k1+1) + index_k+1];
//       }
//       else if (ppr2->sources_k3_interpolation == cubic_interpolation && index_k1 > 3) {
//         for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
//           interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] =
//             a * rsources[index_tau*(index_k1+1) + index_k] + b * rsources[index_tau*(index_k1+1) + index_k+1]
//             + ((a*a*a-a) * derivs[index_tau*(index_k1+1) + index_k]
//             +(b*b*b-b) * derivs[index_tau*(index_k1+1) + index_k+1])*h*h/6.0;
//       }
//
//    }
//   } // end of for (index_k_tr)
//
//
//
//
//
//   free(rsources);
//    free(derivs);

  return _SUCCESS_;

} // spectra_interpolate_sources_in_k2


#undef sources




