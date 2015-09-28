/** @file print_params.c 
 *
 * Test that the new parameters I introduced in CLASS and SONG work correctly.
 * This main needs to be called with a .ini file; the test consists of comparing 
 * the parameters specified in the .ini file against what is printed to screen.
 *
 * Created by Guido W. Pettinari on 15.06.2011
 * Last edited by Guido W. Pettinari on 30.06.2015
 */
 
#include "song.h"

int main(int argc, char **argv) {

  struct precision pr;        /* precision parameters (1st-order) */
  struct precision2 pr2;      /* precision parameters (2nd-order) */
  struct background ba;       /* cosmological background */
  struct thermo th;           /* thermodynamics */
  struct perturbs pt;         /* source functions (1st-order) */
  struct perturbs2 pt2;       /* source functions (2nd-order) */  
  struct transfers tr;        /* transfer functions (1st-order) */
  struct bessels bs;          /* bessel functions (1st-order) */
  struct bessels2 bs2;        /* bessel functions (2nd-order) */
  struct transfers2 tr2;      /* transfer functions (2nd-order) */
  struct primordial pm;       /* primordial spectra */
  struct spectra sp;          /* output spectra (1st-order) */
  struct bispectra bi;        /* bispectra */
  struct fisher fi;           /* fisher matrix */
  struct nonlinear nl;        /* non-linear spectra */
  struct lensing le;          /* lensed spectra */
  struct output op;           /* output files */
  ErrorMsg errmsg;            /* error messages */

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,
    &sp,&nl,&le,&bs,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  
  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&tr,&bs,&bs2,&tr2,&pm,
    &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  printf("~~~~~ Value of my parameters: ~~~~~\n");

  /* Flags */
  // printf("* Flags & switches\n");
  // printf("\tpt2.has_cls = %d\n", pt2.has_cls);
  // printf("\tpt2.has_cmb_temperature = %d\n", pt2.has_cmb_temperature);
  // printf("\tpt2.has_cmb_polarization_e = %d\n", pt2.has_cmb_polarization_e);
  // printf("\tpt2.has_cmb_polarization_b = %d\n", pt2.has_cmb_polarization_b);
  // printf("\tpt2.has_perturbations2 = %d\n", pt2.has_perturbations2);
  // printf("\tpt.has_perturbations2 = %d\n", pt.has_perturbations2);
  // printf("\tpt2.has_polarization2 = %d\n", pt2.has_polarization2);
  // printf("\tpt.has_polarization2 = %d\n", pt.has_polarization2);
  // printf("\tpt2.has_quadratic_sources = %d\n", pt2.has_quadratic_sources);
  // printf("\tpt2.has_quadratic_liouville = %d\n", pt2.has_quadratic_liouville);
  // printf("\tpt2.has_quadratic_collision = %d\n", pt2.has_quadratic_collision);
  // printf("\tpt2.has_perfect_baryons = %d\n", pt2.has_perfect_baryons);
  // printf("\tpt2.has_perfect_cdm = %d\n", pt2.has_perfect_cdm);
  // printf("\tpt.has_cl_cmb_zeta = %d\n", pt.has_cl_cmb_zeta);
  // printf("\tpt.recombination_only_zeta = %d\n", pt.recombination_only_zeta);
  // printf("\tpt2.perturbations2_verbose = %d\n", pt2.perturbations2_verbose);

  
  /* Interpolation */
  // printf("* Interpolation options\n");
  // printf("\tpr.quadsources_time_interpolation = %d\n", pr.quadsources_time_interpolation);
  // printf("\tpr2.sources_time_interpolation = %d\n", pr2.sources_time_interpolation);
  // printf("\tpr2.sources_k3_interpolation = %d\n", pr2.sources_k3_interpolation);
  // printf("\tpr.bessels_interpolation = %d\n", pr.bessels_interpolation);
  // printf("\tpr.transfers_k1_interpolation = %d\n", pr.transfers_k1_interpolation);
  // printf("\tpr.transfers_k2_interpolation = %d\n", pr.transfers_k2_interpolation);
  // printf("\tpr.bispectra_k3_extrapolation = %d\n", pr.bispectra_k3_extrapolation);
  // printf("\tpr.extra_k3_oscillations_right = %g\n", pr.extra_k3_oscillations_right);
  // printf("\tpr.extra_k3_oscillations_left = %g\n", pr.extra_k3_oscillations_left);

  /* Perturbed recombination */
  // printf("* Perturbed recombination\n");
  // printf("\tth.has_perturbed_recombination_stz = %d\n", th.has_perturbed_recombination_stz);
  // printf("\tth.perturbed_recombination_turnx = %g\n", th.perturbed_recombination_turnx);
  // printf("\tth.compute_xe_derivatives = %d\n", th.compute_xe_derivatives);
  // printf("\tpt.has_perturbed_recombination_stz = %d\n", pt.has_perturbed_recombination_stz);
  // printf("\tpt2.has_perturbed_recombination_stz = %d\n", pt2.has_perturbed_recombination_stz);
  // printf("\tpt2.perturbed_recombination_use_approx = %d\n", pt2.perturbed_recombination_use_approx);
  
  /* Equations parameters */
  // printf("\tpt2.phi_prime_eq = %d\n", pt2.phi_prime_eq);
  
  /* LOS effects at first order */
  // printf("* Line-of-sight effects at first order\n");
  // printf("\tpt.has_scattering_in_los = %d\n", pt.has_scattering_in_los);
  // printf("\tpt.has_pure_metric_in_los = %d*\n", pt.has_pure_metric_in_los);
  // printf("\tpt.has_sw = %d\n", pt.has_sw);
  // printf("\tpt.has_isw = %d\n", pt.has_isw);

  /* LOS effects at second order */
  // printf("* Line-of-sight effects at second order\n");
  // printf("\tpt2.has_pure_scattering_in_los = %d\n", pt2.has_pure_scattering_in_los);
  // printf("\tpt2.has_quad_scattering_in_los = %d\n", pt2.has_quad_scattering_in_los);
  // printf("\tpt2.has_pure_metric_in_los = %d*\n", pt2.has_pure_metric_in_los);
  // printf("\tpt2.has_quad_metric_in_los = %d\n", pt2.has_quad_metric_in_los);
  // printf("\tpt2.has_time_delay_in_los = %d\n", pt2.has_time_delay_in_los);
  // printf("\tpt2.has_redshift_in_los = %d\n", pt2.has_redshift_in_los);
  // printf("\tpt2.has_lensing_in_los = %d\n", pt2.has_lensing_in_los);
  // printf("\tpt2.use_delta_tilde_in_los = %d\n", pt2.use_delta_tilde_in_los);
  // printf("\tpt2.has_sw = %d*\n", pt2.has_sw);
  // printf("\tpt2.has_isw = %d*\n", pt2.has_isw);
  // printf("\tpt2.only_early_isw = %d\n", pt2.only_early_isw);
  // printf("\tpt2.has_recombination_only = %d\n", pt2.has_recombination_only);

  /* Time sampling for quadratic sources */
  // printf("* Time sampling for quadratic sources\n");
  // printf("\tpr.perturb_sampling_stepsize_quadsources = %g\n", pr.perturb_sampling_stepsize_quadsources);
  // printf("\tpt.has_custom_timesampling_for_quadsources = %d\n", pt.has_custom_timesampling_for_quadsources);
  // printf("\tpt.custom_tau_ini_quadsources = %g\n", pt.custom_tau_ini_quadsources);
  // printf("\tpt.custom_tau_end_quadsources = %g\n", pt.custom_tau_end_quadsources);
  // printf("\tpt.custom_tau_size_quadsources = %d\n", pt.custom_tau_size_quadsources);
  // printf("\tpt.custom_tau_mode_quadsources = %d\n", pt.custom_tau_mode_quadsources);
  
  /* Time sampling for second-order sources */
  // printf("* Time sampling for second-order sources\n");
  // printf("\tpr2.perturb_sampling_stepsize_song = %g\n", pr2.perturb_sampling_stepsize_song);
  // printf("\tpr2.perturb_sampling_late_time_boost = %g\n", pr2.perturb_sampling_late_time_boost);
  // printf("\tpr2.custom_tau_start_evolution = %g\n", pr2.custom_tau_start_evolution);
  // printf("\tpt2.recombination_max_to_end_ratio = %g\n", pt2.recombination_max_to_end_ratio);
  // printf("\tpt2.has_custom_timesampling = %d\n", pt2.has_custom_timesampling);
  // printf("\tpt2.custom_tau_ini = %g\n", pt2.custom_tau_ini);
  // printf("\tpt2.custom_tau_end = %g\n", pt2.custom_tau_end);
  // printf("\tpt2.custom_tau_size = %d\n", pt2.custom_tau_size);
  // printf("\tpt2.custom_tau_mode = %d\n", pt2.custom_tau_mode);

  /* Wavemode sampling */
  // printf("* Wavemode sampling for second-order sources\n");
  // printf("\tpt2.k_sampling = %d\n", pt2.k_sampling);
  // printf("\tpr2.k_min_tau0 = %g\n", pr2.k_min_tau0);
  // printf("\tpr2.k_max_tau0_over_l_max = %g\n", pr2.k_max_tau0_over_l_max);
  // printf("\tpr2.k_step_sub = %g\n", pr2.k_step_sub);
  // printf("\tpr2.k_step_super = %g\n", pr2.k_step_super);
  // printf("\tpr2.k_logstep_super = %g\n", pr2.k_logstep_super);
  // printf("\tpr2.k_step_transition = %g\n", pr2.k_step_transition);
  // printf("\tpr2.k_min_scalar = %g\n", pr2.k_min_custom);
  // printf("\tpr2.k_max_scalar = %g\n", pr2.k_max_custom);
  // printf("\tpr2.k_size_custom = %d\n", pr2.k_size_custom);
  // printf("\tpr2.q_linstep_song = %g\n", pr2.q_linstep_song);
  // printf("\tpt2.k3_sampling = %d\n", pt2.k3_sampling);
  // printf("\tpr2.k3_size_min = %d\n", pr2.k3_size_min);
  // printf("\tpr2.k3_size = %d\n", pr2.k3_size);
  // printf("\tpr2.tau_linstep_song = %g\n", pr2.tau_linstep_song);
  // printf("\ttr2.k_sampling = %d\n", tr2.k_sampling);
  // printf("\ttr2.tau_sampling = %d\n", tr2.tau_sampling);

  /* Initial conditions */
  // printf("* Initial conditions parameters\n");    
  // printf("\tpt.has_ad_maberty = %d\n", pt.has_ad_maberty);
  // printf("\tpt.has_zero_ic = %d\n", pt.has_zero_ic);  
  // printf("\tpt2.has_ad = %d\n", pt2.has_ad);  
  // printf("\tpt2.has_ad_first_order = %d\n", pt2.has_ad_first_order);  
  // printf("\tpt2.has_zero_ic = %d\n", pt2.has_zero_ic);
  // printf("\tpt2.has_unphysical_ic = %d\n", pt2.has_unphysical_ic);
  // printf("\tpt2.primordial_local_fnl_phi = %g\n", pt2.primordial_local_fnl_phi);

  /* Bispectrum */
  // printf("* Bispectrum parameters\n");
  // printf("\tpt2.has_cmb_bispectra = %d\n", pt2.has_cmb_bispectra);
  // printf("\tbi.has_bispectra = %d\n", bi.has_bispectra);
  // printf("\tbi.has_local_model = %d\n", bi.has_local_model);
  // printf("\tbi.has_equilateral_model = %d\n", bi.has_equilateral_model);
  // printf("\tbi.has_orthogonal_model = %d\n", bi.has_orthogonal_model);
  // printf("\tbi.has_galileon_model = %d\n", bi.has_galileon_model);
  // printf("\tbi.has_intrinsic = %d\n", bi.has_intrinsic);
  // printf("\tbi.has_intrinsic_squeezed = %d\n", bi.has_intrinsic_squeezed);
  // printf("\tbi.has_intrinsic_squeezed_unlensed = %d\n", bi.has_intrinsic_squeezed_unlensed);
  // printf("\tbi.has_quadratic_correction = %d\n", bi.has_quadratic_correction);
  // printf("\tbi.add_quadratic_correction = %d\n", bi.add_quadratic_correction);
  // printf("\tbi.include_lensing_effects = %d\n", bi.include_lensing_effects);
  // printf("\tbi.lensed_intrinsic = %d\n", bi.lensed_intrinsic);
  // printf("\tpr.extend_lensed_cls = %d\n", pr.extend_lensed_cls);
  // printf("\tpr.bispectra_r_sampling = %d\n", pr.bispectra_r_sampling);
  // printf("\tpr.r_left = %g\n", pr.r_left);
  // printf("\tpr.r_right = %g\n", pr.r_right);
  // printf("\tpr.r_min = %g\n", pr.r_min);
  // printf("\tpr.r_max = %g\n", pr.r_max);
  // printf("\tpr.r_size = %d\n", pr.r_size);

  /* Fisher */
  // printf("\tfi.always_interpolate_bispectra = %d\n", fi.always_interpolate_bispectra);
  // printf("\tfi.l_min_estimator = %d\n", fi.l_min_estimator);
  // printf("\tfi.l_max_estimator = %d\n", fi.l_max_estimator);
  // printf("\tfi.bispectra_interpolation = %d\n", fi.bispectra_interpolation);
  // printf("\tfi.f_sky = %g\n", fi.f_sky);
  // for (int index_channel=0; index_channel < fi.n_channels; ++index_channel) {
  //   printf("\tchannel %d: fi.beam = %g arcmins, %g radians\n",
  //   index_channel, fi.beam[index_channel] * 60. / (_PI_/180.), fi.beam[index_channel]);
  //   printf("\t           fi.noise_t = %g uK\n",
  //   index_channel, sqrt(fi.noise_t[index_channel]) * 1e6*ba.T_cmb / fi.beam[index_channel]);
  //   printf("\t           fi.noise_e = %g uK\n",
  //   index_channel, sqrt(fi.noise_e[index_channel]) * 1e6*ba.T_cmb / fi.beam[index_channel]);
  // }
  // printf("\tfi.include_lensing_effects = %d\n", fi.include_lensing_effects);
  // printf("\tfi.compute_lensing_variance_lmax = %d\n", fi.compute_lensing_variance_lmax);
  // printf("\tfi.ignore_t = %d\n", fi.ignore_t);
  // printf("\tfi.ignore_e = %d\n", fi.ignore_e);
  // printf("\tfi.ignore_b = %d\n", fi.ignore_b);
  // printf("\tfi.ignore_r = %d\n", fi.ignore_r);
  // printf("\tfi.squeezed_ratio = %g\n", fi.squeezed_ratio);

  /* Precision parameters - multipoles */
  // printf("* Precision parameters - multipoles\n");
  // printf("\tpr2.m_size = %d, ", pr2.m_size);
  // printf("\tpr2.m_max_2nd_order = %d, ", pr2.m_max_2nd_order);
  // printf("\tpr2.m = ");
  // for (int index_m=0; index_m < pr2.m_size; ++index_m)
  //   printf("%d%s ", pr2.m[index_m], index_m!=(pr2.m_size-1)?",":"\n");
  // printf("\tpr2.l_max_g = %d\n", pr2.l_max_g);
  // printf("\tpr2.l_max_pol_g = %d\n", pr2.l_max_pol_g);
  // printf("\tpr2.l_max_ur = %d\n", pr2.l_max_ur);
  // printf("\tpr2.l_max_g_quadsources = %d\n", pr2.l_max_g_quadsources);
  // printf("\tpr2.l_max_pol_g_quadsources = %d\n", pr2.l_max_pol_g_quadsources);
  // printf("\tpr2.l_max_ur_quadsources = %d\n", pr2.l_max_ur_quadsources);
  // printf("\tpr2.l_max_los_t = %d\n", pr2.l_max_los_t);
  // printf("\tpr2.l_max_los_quadratic_t = %d\n", pr2.l_max_los_quadratic_t);
  // printf("\tpr2.l_max_los_p = %d\n", pr2.l_max_los_p);
  // printf("\tpr2.l_max_los_quadratic_p = %d\n", pr2.l_max_los_quadratic_p);

  /* Precision parameters - integration */
  // printf("* Precision parameters - integration\n");    
  // printf("\tpr.tol_perturb_integration = %g\n", pr.tol_perturb_integration);
  // printf("\tpr2.tol_perturb_integration_song = %g\n", pr2.tol_perturb_integration_song);
  // printf("\tpr2.start_small_k_at_tau_c_over_tau_h_song = %g\n", pr2.start_small_k_at_tau_c_over_tau_h_song);
  // printf("\tpr2.start_large_k_at_tau_h_over_tau_k_song = %g\n", pr2.start_large_k_at_tau_h_over_tau_k_song);

  /* Precision parameters - bessels */
  // printf("* Precision parameters - Bessels\n");
  // printf("\tpr2.bessel_j_cut_song = %g\n", pr2.bessel_j_cut_song);
  // printf("\tpr2.bessel_J_cut_song = %g\n", pr2.bessel_J_cut_song);
  // printf("\tpr2.bessel_x_step_song = %g\n", pr2.bessel_x_step_song);
  // printf("\tpr2.compute_only_even_ls = %d\n", pr2.compute_only_even_ls);
  // printf("\tpr2.compute_only_odd_ls = %d\n", pr2.compute_only_odd_ls);
  // printf("\tbs2.extend_l1_using_m = %d\n", bs2.extend_l1_using_m);

  /* Perturbations output and debug parameters */
  printf("* Debug parameters\n");
  printf("\tpt2.file_verbose = %d\n", pt2.file_verbose);
  printf("\tpt2.tau_out_size = %d\n", pt2.tau_out_size);
  if (pt2.tau_out_size > 0) {
    printf ("\t\t");
    for (int index_tau_out=0; index_tau_out < pt2.tau_out_size; ++index_tau_out)
      printf("%12g ", pt2.tau_out[index_tau_out]);
    printf ("\n");
  }
  printf("\tpt2.z_out_size = %d\n", pt2.z_out_size);
  if (pt2.z_out_size > 0) {
    printf ("\t\t");
    for (int index_z_out=0; index_z_out < pt2.z_out_size; ++index_z_out)
      printf("%12g ", pt2.z_out[index_z_out]);
    printf ("\n");
  }
  printf("\tpt2.k_out_size = %d\n", pt2.k_out_size);
  for (int index_k_out=0; index_k_out < pt2.k_out_size; ++index_k_out)
    printf("\t\t%12g %12g %12g\n", pt2.k1_out[index_k_out], pt2.k2_out[index_k_out], pt2.k3_out[index_k_out]);
  printf("\tpt2.k_index_out_size = %d\n", pt2.k_index_out_size);
  for (int index_k_out=0; index_k_out < pt2.k_index_out_size; ++index_k_out)
    printf("\t\t%12d %12d %12d\n", pt2.k1_index_out[index_k_out], pt2.k2_index_out[index_k_out], pt2.k3_index_out[index_k_out]);
  printf("\tpt2.k_out_mode = %d\n", pt2.k_out_mode);
  printf ("\tpt2.output_class_perturbations = %d\n", pt2.output_class_perturbations);
  printf ("\tpt2.output_quadratic_sources = %d\n", pt2.output_quadratic_sources);

  /* Storage parameters */
  // printf("* Parameters related to the storage of intermediate results\n");
  // printf("\tpr.store_run = %d\n", pr.store_run);
  // printf("\tpr.load_run = %d\n", pr.load_run);
  // printf("\tpr.append_date_to_run = %d\n", pr.append_date_to_run);
  // printf("\tpr.run_dir = %s\n", pr.run_dir);
  // printf("\tpr.data_dir = %s\n", pr.data_dir);
  // printf("\tpr2.store_sources_to_disk = %d\n", pr2.store_sources_to_disk);
  // printf("\tpt2.sources_dir = %s\n", pt2.sources_dir);
  // printf("\tpr2.store_transfers_to_disk = %d\n", pr2.store_transfers_to_disk);
  // printf("\ttr2.transfers_dir = %s\n", tr2.transfers_dir);
  // printf("\tpr.store_bispectra_to_disk = %d\n", pr.store_bispectra_to_disk);
  // printf("\tbi.bispectra_dir = %s\n", bi.bispectra_dir);

  /* Approximations parameters */
  // printf("* Parameters related to the approximations adopted at second-order\n");
  // printf("\tpt2.tight_coupling_approximation = %d\n", pt2.tight_coupling_approximation);
  // printf("\tpt2.tight_coupling_trigger_tau_c_over_tau_h = %g\n", pt2.tight_coupling_trigger_tau_c_over_tau_h);
  // printf("\tpt2.tight_coupling_trigger_tau_c_over_tau_k = %g\n", pt2.tight_coupling_trigger_tau_c_over_tau_k);
  // printf("\tpt2.radiation_streaming_approximation = %d\n", pt2.radiation_streaming_approximation);
  // printf("\tpt2.radiation_streaming_trigger_tau_over_tau_k = %g\n", pt2.radiation_streaming_trigger_tau_over_tau_k);
  // printf("\tpt2.ur_fluid_approximation = %d\n", pt2.ur_fluid_approximation);
  // printf("\tpt2.ur_fluid_trigger_tau_over_tau_k = %g\n", pt2.ur_fluid_trigger_tau_over_tau_k);
  // printf("\tpt2.no_radiation_approximation = %d\n", pt2.no_radiation_approximation);
  // printf("\tpt2.no_radiation_approximation_rho_m_over_rho_r = %g\n", pt2.no_radiation_approximation_rho_m_over_rho_r);


  return _SUCCESS_;

}
