/** @file input2.c
 *
 * Input module for SONG.
 *
 * Created by Guido W. Pettinari on 13.03.2013
 * Based on input.c by the CLASS team (http://class-code.net/)
 */

#include "input2.h" 

/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 */

int input2_init_from_arguments(
            int argc, 
            char **argv,
            struct precision * ppr,
            struct precision2 * ppr2,
            struct background *pba,
            struct thermo *pth,
            struct perturbs *ppt,
            struct perturbs2 *ppt2,
            struct transfers *ptr,
            struct bessels * pbs,
            struct bessels2 * pbs2,
            struct transfers2 *ptr2,
            struct primordial *ppm,
            struct spectra *psp,
            struct nonlinear *pnl,
            struct lensing *ple,
            struct bispectra *pbi,
            struct fisher *pfi,
            struct output *pop,
            ErrorMsg errmsg
            )
{

  /* Define local variables */

  struct file_content fc;
  struct file_content fc_input;
  struct file_content fc_precision;

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];


  /* Initialize the file_content structures */
  fc.size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';


  /* We parsed the names of the parameter files in input.c */

  strcpy (input_file, ppr->ini_filename);
  strcpy (precision_file, ppr->pre_filename);
  

  /* If there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0] != '\0')
    
    class_call(parser_read_file(input_file,&fc_input,errmsg),
         errmsg,
         errmsg);
         

  /* If there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0] != '\0')
    
    class_call(parser_read_file(precision_file,&fc_precision,errmsg),
         errmsg,
         errmsg);


  /* If one or two files were read, merge their contents in a single 'file_content' structure. */

  if ((input_file[0]!='\0') || (precision_file[0]!='\0'))

    class_call(parser_cat(&fc_input,&fc_precision,&fc,errmsg),
         errmsg,
         errmsg);

  class_call(parser_free(&fc_input),errmsg,errmsg);
  class_call(parser_free(&fc_precision),errmsg,errmsg);
  

  /* Now, initialize all parameters given the input 'file_content' structure.  If its size
  is null, all parameters take their default values. */

  class_call (input2_init(
                &fc,
                // ppr->input_file_content, /* TODO: why is this not working? */
                ppr,
                ppr2,
                pba,
                pth,
                ppt,
                ppt2,
                ptr,
                pbs,
                pbs2,
                ptr2,
                ppm,
                psp,
                pnl,
                ple,
                pbi,
                pfi,
                pop,
    errmsg),
    errmsg,
    errmsg);

  class_call(parser_free(&fc),errmsg,errmsg);
    
  return _SUCCESS_;
}




/**
 * Initialize each parameters, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 */

int input2_init (
         struct file_content * pfc,
         struct precision * ppr,
         struct precision2 * ppr2,
         struct background *pba,
         struct thermo *pth,
         struct perturbs *ppt,
         struct perturbs2 *ppt2,         
         struct transfers *ptr,
         struct bessels * pbs,
         struct bessels2 * pbs2,
         struct transfers2 *ptr2,
         struct primordial *ppm,
         struct spectra *psp,
         struct nonlinear * pnl,
         struct lensing *ple,
         struct bispectra *pbi,
         struct fisher *pfi,
         struct output *pop,
         ErrorMsg errmsg
         )
{

  printf("Running SONG version %s\n", _SONG_VERSION_);

  /** Summary: */

  /** - define local variables */

  int flag1,flag2,flag3,flag;
  int int1;
  int found;
  double param1,param2,param3;
  int entries_read;
  double * pointer, * pointer1, * pointer2;
  int * int_pointer, * int_pointer1, * int_pointer2;
  int * pointer_to_int;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  char string3[_ARGUMENT_LENGTH_MAX_];
  char string[_ARGUMENT_LENGTH_MAX_];
  int i;


  /** - set all parameters (input and precision) to default values */

  class_call (input2_default_params(
               pba,
               pth,
               ppt,
               ppt2,
               ptr,
               pbs,
               pbs2,
               ptr2,
               ppm,
               psp,
               pnl,
               ple,
               pbi,
               pfi,
               pop),
    errmsg,
    errmsg);

  class_call (input2_default_precision (ppr2),
    errmsg,
    errmsg);


  
  // ==================================================================================
  // =                                  Parse output                                  =
  // ==================================================================================
 
  class_call (parser_read_string (pfc,"output",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"tBisp") != NULL) || (strstr(string1,"TBISP") != NULL)) {
      ppt2->has_cmb_temperature = _TRUE_;  
    }
    
    if (((strstr(string1,"pBisp") != NULL) || (strstr(string1,"PBISP") != NULL)) /* obsolete */
       ||(strstr(string1,"eBisp") != NULL) || (strstr(string1,"EBISP") != NULL)) {
      ppt2->has_cmb_polarization_e = _TRUE_;
      ppt->has_cl_cmb_temperature = _TRUE_;
    }

    if ((strstr(string1,"bBisp") != NULL) || (strstr(string1,"BBISP") != NULL)) {
      ppt2->has_cmb_polarization_b = _TRUE_;
      ppt->has_cl_cmb_temperature = _TRUE_;
    }

    if ((strstr(string1,"tCl2") != NULL) || (strstr(string1,"TCL2") != NULL)) {
      ppt2->has_perturbations2 = _TRUE_;
      ppt2->has_cmb_spectra = _TRUE_;
      ppt2->has_cmb_temperature = _TRUE_;
    }
    
    if ((strstr(string1,"pCl2") != NULL) || (strstr(string1,"PCL2") != NULL) ||
        (strstr(string1,"eCl2") != NULL) || (strstr(string1,"ECL2") != NULL)) {
      ppt2->has_perturbations2 = _TRUE_;
      ppt2->has_cmb_spectra = _TRUE_;
      ppt2->has_cmb_polarization_e = _TRUE_;
    }

    if ((strstr(string1,"bCl2") != NULL) || (strstr(string1,"BCL2") != NULL)) {
      ppt2->has_perturbations2 = _TRUE_;
      ppt2->has_cmb_spectra = _TRUE_;
      ppt2->has_cmb_polarization_b = _TRUE_;
    }

    if (((strstr(string1,"delta_cdm_pk") != NULL) || (strstr(string1,"pk_delta_cdm") != NULL) || (strstr(string1,"mPk") != NULL))) {
      ppt2->has_pk_delta_cdm = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }

    if (((strstr(string1,"delta_b_pk") != NULL) || (strstr(string1,"pk_delta_b") != NULL) || (strstr(string1,"mPk") != NULL))) {
      ppt2->has_pk_delta_b = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }

    if (((strstr(string1,"magnetic_pk") != NULL) || (strstr(string1,"pk_magnetic") != NULL) || (strstr(string1,"magPk") != NULL))) {
      ppt2->has_pk_magnetic = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }

    if (((strstr(string1,"delta_cdm_bk") != NULL) || (strstr(string1,"bk_delta_cdm") != NULL) || (strstr(string1,"mBisp") != NULL))) {
      ppt2->has_bk_delta_cdm = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }

    if ((strstr(string1,"stop_at_perturbations1") != NULL) || (strstr(string1,"P1") != NULL)) {
      ppt2->stop_at_perturbations1 = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }

    if ((strstr(string1,"stop_at_perturbations2") != NULL) || (strstr(string1,"P2") != NULL)) {
      ppt2->stop_at_perturbations2 = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;      
    }

    if ((strstr(string1,"stop_at_transfers2") != NULL) || (strstr(string1,"T2") != NULL)) {
      ptr2->stop_at_transfers2 = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
      ppt2->has_cmb_spectra = _TRUE_;
    }

    if (strstr(string1,"k_out") != NULL) {
      ppt2->k_out_mode = _TRUE_;
      ppt2->stop_at_perturbations2 = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;
    }
    
  } // end of output parsing

      
  /* Parse the needed types of bispectra (intrinsic, primordial...) */

  class_call (parser_read_string(pfc,"bispectrum_types",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if ((pbi->has_bispectra == _TRUE_) && (flag1 == _TRUE_)) {
 
    /* Intrinsic bispectrum. This is induced by second-order effects in the evolution of the cosmological
    perturbations. */
    if (strstr(string1,"intrinsic") != NULL) {
      ppt2->has_perturbations2 = _TRUE_;
      ppt2->has_cmb_bispectra = _TRUE_;
    }
    
  } // end of bispectrum_types parsing
  

  /* If second-order perturbations are requested, make sure to compute the derivatives of
  the baryon sound-speed. These are needed to obtain the contribution of the effective baryon
  pressure at second-order */
  if (ppt2->has_perturbations2 == _TRUE_)
    pth->compute_cb2_derivatives = _TRUE_;


  // ====================================================================================
  // =                           Perturbations, time sampling                           =
  // ====================================================================================

  class_read_double("tau_start_evolution_2nd_order", ppr2->custom_tau_start_evolution); /* obsolete */
  class_read_double("tau_start_evolution_song", ppr2->custom_tau_start_evolution);

  class_read_double("recombination_max_to_end_ratio", ppt2->recombination_max_to_end_ratio);

  class_read_string_one_of_two(pfc,
    "custom_time_sampling_for_2nd_order_sources",
    "custom_time_sampling_song_sources");   
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_custom_timesampling = _TRUE_;

  class_read_double("custom_tau_ini_2nd_order_sources", ppt2->custom_tau_ini); /* obsolete */
  class_read_double("custom_tau_ini_song_sources", ppt2->custom_tau_ini);

  class_test (ppt2->custom_tau_ini<=0, errmsg, "please choose 'tau_ini' greater than zero.");
  
  class_read_double("custom_tau_end_2nd_order_sources", ppt2->custom_tau_end); /* obsolete */
  class_read_double("custom_tau_end_song_sources", ppt2->custom_tau_end);  
  
  class_read_int("custom_tau_size_2nd_order_sources", ppt2->custom_tau_size); /* obsolete */
  class_read_int("custom_tau_size_song_sources", ppt2->custom_tau_size);

  class_read_string_one_of_two(pfc,
    "custom_tau_mode_2nd_order_sources",
    "custom_tau_mode_song_sources");
    
  if (flag == _TRUE_) {

    if (((strstr(string,"lin") != NULL) || (strstr(string,"LIN") != NULL)))
      ppt2->custom_tau_mode = lin_tau_sampling;

    else if (((strstr(string,"log") != NULL) || (strstr(string,"LOG") != NULL)))
      ppt2->custom_tau_mode = log_tau_sampling;
    
    else
      class_stop(errmsg,         
        "tau_mode_2nd_order_sources=%s not supported. Choose between 'lin' and 'log'",
        string);
  }
   
  /* Precision parameters for the time sampling */
  class_read_double("perturb_sampling_stepsize_for_2nd_order",
    ppr2->perturb_sampling_stepsize_song); /* obsolete */
  class_read_double("start_small_k_at_tau_c_over_tau_h_2nd_order",
    ppr2->start_small_k_at_tau_c_over_tau_h_song); /* obsolete */
  class_read_double("start_large_k_at_tau_h_over_tau_k_2nd_order",
    ppr2->start_large_k_at_tau_h_over_tau_k_song); /* obsolete */

  class_read_double("perturb_sampling_stepsize_song",
    ppr2->perturb_sampling_stepsize_song);
  class_read_double("start_small_k_at_tau_c_over_tau_h_song",
    ppr2->start_small_k_at_tau_c_over_tau_h_song);
  class_read_double("start_large_k_at_tau_h_over_tau_k_song",
    ppr2->start_large_k_at_tau_h_over_tau_k_song);

  class_read_double("perturb_sampling_late_time_boost",
    ppr2->perturb_sampling_late_time_boost);
    

  // ===========================================================================
  // =                        Perturbations, k-sampling                        =
  // ===========================================================================

  /* What type of  sampling for k1 and k2? */

  class_call(parser_read_string(pfc,"sources2_k_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strcmp(string1,"lin") == 0) || (strcmp(string1,"linear") == 0))
      ppt2->k_sampling = lin_k_sampling;

    else if ((strcmp(string1,"log") == 0) || (strcmp(string1,"logarithmic")) == 0)
      ppt2->k_sampling = log_k_sampling;

    else if (strcmp(string1,"smart") == 0)
      ppt2->k_sampling = smart_sources_k_sampling;
      
    else
      class_stop(errmsg,
        "sources2_k_sampling=%s not supported. Choose between 'lin', 'log', 'class' and 'smart'.",
        string1);
  }

  /* Parameters for the smart k1-k2 sampling */

  class_read_double("k_scalar_min_tau0_2nd_order",ppr2->k_min_tau0); /* obsolete */
  class_read_double("k_scalar_max_tau0_over_l_max_2nd_order",ppr2->k_max_tau0_over_l_max); /* obsolete */
  class_read_double("k_scalar_step_sub_2nd_order",ppr2->k_step_sub); /* obsolete */
  class_read_double("k_scalar_linstep_super_2nd_order",ppr2->k_step_super); /* obsolete */
  class_read_double("k_linstep_super_2nd_order",ppr2->k_step_super); /* obsolete */
  class_read_double("k_scalar_logstep_super_2nd_order",ppr2->k_logstep_super); /* obsolete */
  class_read_double("k_scalar_step_transition_2nd_order",ppr2->k_step_transition); /* obsolete */

  class_read_double("k_min_tau0_song",ppr2->k_min_tau0);
  class_read_double("k_max_tau0_over_l_max_song",ppr2->k_max_tau0_over_l_max);
  class_read_double("k_step_sub_song",ppr2->k_step_sub);
  class_read_double("k_step_super_song",ppr2->k_step_super);
  class_read_double("k_logstep_super_song",ppr2->k_logstep_super);
  class_read_double("k_step_transition_song",ppr2->k_step_transition);
  
  class_read_double("k_per_decade_for_pk_song",ppr2->k_per_decade_for_pk);
  class_read_double("k_per_decade_for_bao_song",ppr2->k_per_decade_for_bao);
  class_read_double("k_bao_center",ppr->k_bao_center);
  class_read_double("k_bao_width",ppr->k_bao_width);

  class_test (ppr2->k_logstep_super <= 1,
    ppr2->error_message,
    "logarithmic step must be larger than 1");

  if (ppt2->has_pk_delta_cdm || ppt2->has_bk_delta_cdm || ppt2->has_pk_magnetic) {

    class_call(parser_read_double(pfc,"P_k_max_h/Mpc_song",&param1,&flag1,errmsg),errmsg,errmsg);
    class_call(parser_read_double(pfc,"P_k_max_1/Mpc_song",&param2,&flag2,errmsg),errmsg,errmsg);
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_), errmsg,
      "In input file, you cannot enter both P_k_max_h/Mpc and P_k_max_1/Mpc, choose one");

    if (flag1 == _TRUE_)
      ppt2->k_max_for_pk = param1*pba->h;
    if (flag2 == _TRUE_)
      ppt2->k_max_for_pk = param2;

    /* Shortcut */
    class_read_double("k_max_for_pk",ppt2->k_max_for_pk);    

    // class_call(parser_read_list_of_doubles(
    //              pfc,
    //              "z_pk_song",
    //              &(int1),
    //              &(pointer1),
    //              &flag1,
    //              errmsg),
    //   errmsg,
    //   errmsg);
    //
    // if (flag1 == _TRUE_) {
    //   class_test(int1 > _Z_PK_NUM_MAX_, errmsg, "please increase _Z_PK_NUM_MAX_ output.h to %d", int1);
    //   pop->z_pk_num = int1;
    //   for (i=0; i<int1; i++)
    //     pop2->z_pk[i] = pointer1[i];
    //   free(pointer1);
    // }
    //
    // class_call(parser_read_double(pfc,"z_max_pk",&param1,&flag1,errmsg),errmsg,errmsg);
    //
    // if (flag1==_TRUE_) {
    //   psp->z_max_pk = param1;
    // }
    // else {
    //   psp->z_max_pk = 0.;
    //   for (i=0; i<pop->z_pk_num; i++)
    //     psp->z_max_pk = MAX(psp->z_max_pk,pop->z_pk[i]);
    // }

  } // if LSS

    
  /* Parameters for the custom lin/log sampling */

  class_read_double("k_min_scalars", ppr2->k_min_custom); /* obsolete */
  class_read_double("k_max_scalars", ppr2->k_max_custom); /* obsolete */
  class_read_int("k_size_scalars", ppr2->k_size_custom); /* obsolete */

  class_read_double("k_min_custom_song", ppr2->k_min_custom);
  class_read_double("k_max_custom_song", ppr2->k_max_custom);
  class_read_int("k_size_custom_song", ppr2->k_size_custom);


  /* - k3 sampling of the sources */

  class_call(parser_read_string(pfc,"sources2_k3_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strcmp(string1,"lin") == 0) || (strcmp(string1,"linear") == 0))
      ppt2->k3_sampling = lin_k3_sampling;

    else if ((strcmp(string1,"log") == 0) || (strcmp(string1,"logarithmic")) == 0)
      ppt2->k3_sampling = log_k3_sampling;

    else if ((strcmp(string1,"smart") == 0) || (strcmp(string1,"class") == 0))
      ppt2->k3_sampling = smart_k3_sampling;

    else if ((strcmp(string1,"sym") == 0) || (strcmp(string1,"symmetric") == 0))
      ppt2->k3_sampling = sym_k3_sampling;

    else if (strcmp(string1,"theta_12") == 0)
      ppt2->k3_sampling = theta12_k3_sampling;
      
    else
      class_stop(errmsg,
        "sources2_k3_sampling=%s not supported. Choose between 'lin', 'log', 'smart', 'theta_12'",
        string1);
  }
  
  /* Minimum number of grid points for any (k1,k2) pair, used when k3_sampling is set to smart */
  class_read_int("k3_size_min", ppr2->k3_size_min);
  
  /* Fixed number of grid points for any (k1,k2) pair, used when k3_sampling is set to either lin or log */
  class_read_int("k3_size", ppr2->k3_size);

  if (ppt2->has_cmb_bispectra && ppt2->k3_sampling == sym_k3_sampling) {
    printf ("\nWARNING: symmetric sampling not supported for intrinsic bispectrum; switching to 'smart' sampling\n\n");
    ppt2->k3_sampling = smart_k3_sampling;
  }

  if (ppt2->has_cmb_spectra && ppt2->k3_sampling == sym_k3_sampling) {
    printf ("\nWARNING: symmetric sampling not supported for second-order C_l; switching to 'smart' sampling\n\n");
    ppt2->k3_sampling = smart_k3_sampling;
  }

  /* DISABLED: Uncomment once you sort out the spectra flags */  
  //  if (ppt2->has_cmb_spectra_fourier && ppt2->k3_sampling != sym_k3_sampling) {
  //    printf ("WARNING: you are computing second-order Fourier spectra without a symmetric k3 sampling\
  // (sources2_k3_sampling=sym); expect the accuracy of the resulting P(k) to decrease for small k\n");
  //  }


  // ====================================================================================
  // =                        Perturbations, differential system                        =
  // ====================================================================================

  class_call(parser_read_string(pfc,"quadratic_sources",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_sources = _FALSE_;

  class_call(parser_read_string(pfc,"quadratic_liouville",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_liouville = _FALSE_;    

  class_call(parser_read_string(pfc,"quadratic_collision",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_collision = _FALSE_;    

  if (ppt2->has_quadratic_sources == _FALSE_) {
    ppt2->has_quadratic_liouville = _FALSE_;
    ppt2->has_quadratic_collision = _FALSE_;
  }

  class_call(parser_read_string(pfc,"polarization_second_order",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt2->has_polarization2 = _FALSE_;
  }

  if (ppt2->has_polarization2 == _TRUE_)
    ppt->has_polarization2 = _TRUE_;
  
  class_call(parser_read_string(pfc,"perfect_baryons",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt2->has_perfect_baryons = _FALSE_;
  }

  class_call(parser_read_string(pfc,"perfect_cdm",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt2->has_perfect_cdm = _FALSE_;
  }

  class_read_double("tol_perturb_integration_2nd_order",ppr2->tol_perturb_integration_song); /* obsolete */
  class_read_double("tol_perturb_integration_song",ppr2->tol_perturb_integration_song);


  // ====================================================================================
  // =                      Perturbations, perturbed recombination                      =
  // ====================================================================================

  class_call(parser_read_string(pfc,"perturbed_recombination_song",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    pth->has_perturbed_recombination_stz = _TRUE_;
    ppt->has_perturbed_recombination_stz = _TRUE_;
    ppt2->has_perturbed_recombination_stz = _TRUE_;
    pth->compute_xe_derivatives = _TRUE_; /* Debug purposes */
  }

  /* From which value of x should we start integrating delta_Xe? */
  class_read_double("perturbed_recombination_turnx", pth->perturbed_recombination_turnx);

  /* Should we use the approximation in eq. 3.23 of Senatore, Tassev & Zaldarriaga 2009? */
  class_call(parser_read_string(pfc,"perturbed_recombination_use_approx",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    ppt2->perturbed_recombination_use_approx = _TRUE_;
    
    /* To use the analytical approximation, we only need the derivatives of the background ionization fraction X_e */
    if (pth->has_perturbed_recombination_stz == _TRUE_)
      pth->compute_xe_derivatives = _TRUE_;
      
    pth->has_perturbed_recombination_stz = _FALSE_;
    ppt->has_perturbed_recombination_stz = _FALSE_;
  }
  

  // ===============================================================================
  // =                           Perturbations, LOS sources                        =
  // ===============================================================================

  class_read_string_one_of_two(pfc,
    "include_pure_scattering_in_los_2nd_order",
    "include_pure_scattering_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_pure_scattering_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_quad_scattering_in_los_2nd_order",
    "include_quad_scattering_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_quad_scattering_in_los = _TRUE_;

  class_read_string_one_of_three(pfc,
    "include_metric_in_los_2nd_order",
    "include_metric_song",
    "include_pure_metric_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_pure_metric_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_quad_metric_in_los_2nd_order",
    "include_quad_metric_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_quad_metric_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_time_delay_in_los_2nd_order",
    "include_time_delay_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_time_delay_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_redshift_in_los_2nd_order",
    "include_redshift_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_redshift_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_lensing_in_los_2nd_order",
    "include_lensing_song");
  if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
    ppt2->has_lensing_in_los = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_sachs_wolfe_in_los_2nd_order",
    "include_sachs_wolfe_song");
    if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
      ppt2->has_sw = _TRUE_;

  class_read_string_one_of_two(pfc,
    "include_integrated_sachs_wolfe_in_los_2nd_order",
    "include_integrated_sachs_wolfe_song");
    if ((flag == _TRUE_) && ((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
      ppt2->has_isw = _TRUE_;

  /* Avoid counting twice the same metric effects. SW and ISW override metric. */
  if ((ppt2->has_sw == _TRUE_) || (ppt2->has_isw == _TRUE_))
    ppt2->has_pure_metric_in_los = _FALSE_;

  class_call(parser_read_string(pfc,"use_delta_tilde_in_los",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->use_delta_tilde_in_los = _TRUE_;

  /* Should we define the SW and ISW effects as in Huang and Vernizzi 2013, i.e. using
  the exponential form of the potentials in the metric? This is required to match the
  analytical approximation for the squeezed bispectrum. */
  class_call(parser_read_string(pfc,"use_exponential_potentials",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->use_exponential_potentials = _TRUE_;
  
  /* Should we include only the early part of the ISW effect? */
  if (ppt2->has_isw == _TRUE_) {
    class_call(parser_read_string(pfc,"only_early_isw",&(string1),&(flag1),errmsg),errmsg,errmsg);
    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppt2->only_early_isw = _TRUE_;
  }

  /* If late time effects are not included, stop integrating the second-order system just
  after recombination */
  if ((ppt2->has_pure_metric_in_los == _FALSE_)
   && (ppt2->has_quad_metric_in_los == _FALSE_)
   && ((ppt2->has_isw == _FALSE_) || (ppt2->only_early_isw == _TRUE_))
   && (ppt2->has_time_delay_in_los == _FALSE_)
   && (ppt2->has_redshift_in_los == _FALSE_)
   && (ppt2->has_lensing_in_los == _FALSE_)
   && (ppt2->has_pk_delta_cdm == _FALSE_)
   && (ppt2->has_pk_magnetic == _FALSE_)
   && (ppt2->has_bk_delta_cdm == _FALSE_)
   && (pth->reio_parametrization == reio_none))
    ppt2->has_only_recombination = _TRUE_;
  
  /* Should we stop integrating the second-order system just after recombination, regardless of the
  other flags? */
  class_call(parser_read_string(pfc,"only_recombination",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_only_recombination = _TRUE_;

  /* Should we consider only the sources at reionisation? */
  class_call(parser_read_string(pfc,"only_reionisation",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_only_reionisation = _TRUE_;

  class_call(parser_read_string(pfc,"only_loss_term",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_only_loss_term = _TRUE_;

  class_call(parser_read_string(pfc,"only_gain_term",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_only_gain_term = _TRUE_;

  class_test (ppt2->has_only_loss_term && ppt2->has_only_gain_term,
    errmsg,
    "choose between gain and loss term");

  /* Doesn't make sense not to have polarisation, if you want to compute polarisation */
  class_test ((ppt2->has_polarization2 == _FALSE_) &&
    ((ppt2->has_cmb_polarization_e == _TRUE_) || ((ppt2->has_cmb_polarization_b == _TRUE_))),
    errmsg,
    "please make sure that 'polarization_second_order' is set to 'yes' if you want to compute E-modes or B-modes");

  /* When delta_tilde is activated, the redshift term is automatically included  */
  class_test ((ppt2->use_delta_tilde_in_los==_TRUE_) && (ppt2->has_redshift_in_los==_TRUE_),
    errmsg,
    "the delta_tilde options are not compatible with 'include_redshift_in_los_2nd_order'");

	/* This option should be last. If true, compute only the test source term */
  class_call(parser_read_string(pfc,"use_test_source",&(string1),&(flag1),errmsg),errmsg,errmsg);
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    ppt2->use_test_source = _TRUE_;
    ppt2->has_sw = _FALSE_;
    ppt2->has_isw = _FALSE_;
		ppt2->has_only_recombination = _TRUE_;	
	}


  // ====================================================================================
  // =                        Perturbations, initial conditions                         =
  // ====================================================================================

  /* Only adiabatic or vanishing initial conditions are supported so far */
  class_call(parser_read_string(pfc,"ic_2nd_order",&string1,&flag1,errmsg),
    errmsg,
    errmsg); /* obsolete */
  class_call(parser_read_string(pfc,"ic_song",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {

    /* If no initial conditions are specified, the default is has_ad=_TRUE_; 
    but if they are specified we should reset has_ad to _FALSE_ before reading */
    ppt2->has_ad=_FALSE_;

    if ((strcmp(string1,"ad") == 0) || (strcmp(string1,"AD") == 0))
      ppt2->has_ad=_TRUE_; 

    if ((strstr(string1,"ad_first_order") != NULL) || (strstr(string1,"AD_1ST") != NULL)
      || (strstr(string1,"ad1") != NULL))
      ppt2->has_ad_first_order=_TRUE_; 

    if ((strstr(string1,"zero") != NULL) || (strstr(string1,"ZERO") != NULL))
      ppt2->has_zero_ic=_TRUE_; 

    if (strstr(string1,"unphysical") != NULL)
      ppt2->has_unphysical_ic=_TRUE_; 
      
    class_test(ppt2->has_ad==_FALSE_ && ppt2->has_ad_first_order==_FALSE_
      && ppt2->has_zero_ic==_FALSE_ && ppt2->has_unphysical_ic==_FALSE_,
      errmsg,         
      "ic_2nd_order=%s not supported. Choose between 'ad', 'ad_first_order', 'zero', 'unphysical'",
      string1);
  }

  class_read_double("primordial_local_fnl_phi", ppt2->primordial_local_fnl_phi);


  // ====================================================================================
  // =                           Perturbations, approximations                          =
  // ====================================================================================

  /* Tight coupling.  Note that, contrary to the 1st order case, we read the approximation
  settings directly into the ppt2 structure rather than in the precision one.  We do so in
  order to keep all second-order related quantities in the same structure */
  class_read_int("tight_coupling_approximation_2nd_order",
    ppt2->tight_coupling_approximation); /* obsolete */
  class_read_double("tight_coupling_trigger_tau_c_over_tau_h_2nd_order",  
    ppt2->tight_coupling_trigger_tau_c_over_tau_h); /* obsolete */
  class_read_double("tight_coupling_trigger_tau_c_over_tau_k_2nd_order",
    ppt2->tight_coupling_trigger_tau_c_over_tau_k); /* obsolete */

  class_read_int("tight_coupling_approximation_song",
    ppt2->tight_coupling_approximation);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_h_song",  
    ppt2->tight_coupling_trigger_tau_c_over_tau_h);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_k_song",
    ppt2->tight_coupling_trigger_tau_c_over_tau_k);

  /* Radiation streaming */
  class_read_int("radiation_streaming_approximation_2nd_order",
    ppt2->radiation_streaming_approximation); /* obsolete */
  class_read_double("radiation_streaming_trigger_tau_over_tau_k_2nd_order",
   ppt2->radiation_streaming_trigger_tau_over_tau_k); /* obsolete */

  class_read_int("radiation_streaming_approximation_song",
    ppt2->radiation_streaming_approximation);
  class_read_double("radiation_streaming_trigger_tau_over_tau_k_song",
   ppt2->radiation_streaming_trigger_tau_over_tau_k);
  class_read_double("radiation_streaming_trigger_tau_c_over_tau_song",
   ppt2->radiation_streaming_trigger_tau_c_over_tau);
   
   /* CLASS parameter radiation_streaming_trigger_tau_c_over_tau is used
   not only in the perturbations.c module, but also in the thermodynamics
   module to compute the start of the free-streaming regime. This value
   in the pth structure, in turn, will affect the starting time of the
   second-order RSA. In order to avoid discrepancies, we set CLASS value
   of radiation_streaming_trigger_tau_c_over_tau to match the SONG one. */
   ppr->radiation_streaming_trigger_tau_c_over_tau
     = ppt2->radiation_streaming_trigger_tau_c_over_tau;

  /* Ultra relativistic fluid approximation */
  class_read_int("ur_fluid_approximation_2nd_order",
    ppt2->ur_fluid_approximation); /* obsolete */
  class_read_double("ur_fluid_trigger_tau_over_tau_k_2nd_order",
    ppt2->ur_fluid_trigger_tau_over_tau_k); /* obsolete */

  class_read_int("ur_fluid_approximation_song",
    ppt2->ur_fluid_approximation);
  class_read_double("ur_fluid_trigger_tau_over_tau_k_song",
    ppt2->ur_fluid_trigger_tau_over_tau_k);

  /* No radiation approximation */
  class_read_int("no_radiation_approximation_2nd_order",
  ppt2->no_radiation_approximation); /* obsolete */
  class_read_double("no_radiation_approximation_rho_m_over_rho_r_2nd_order",
    ppt2->no_radiation_approximation_rho_m_over_rho_r); /* obsolete */

  class_read_int("no_radiation_approximation_song",
  ppt2->no_radiation_approximation);
  class_read_double("no_radiation_approximation_rho_m_over_rho_r_song",
    ppt2->no_radiation_approximation_rho_m_over_rho_r);

  if (ppt2->ur_fluid_approximation != ufa2_none) {
    class_test(ppt2->ur_fluid_trigger_tau_over_tau_k
      ==ppt2->radiation_streaming_trigger_tau_over_tau_k, errmsg,
      "please choose different values for precision parameters\
 ur_fluid_trigger_tau_over_tau_k_song and radiation_streaming_trigger_tau_over_tau_k_song\
 , in order to avoid switching two approximation schemes at the same time");
  }


  // ====================================================================================
  // =                        Perturbations, gauges & equations                         =
  // ====================================================================================

  class_read_string_one_of_two(pfc,
    "phi_prime_equation",  /* obsolete */
    "phi_equation_song");

     
  if (flag == _TRUE_) {
  
    if ( (strcmp(string,"poisson") == 0) || (strcmp(string,"POISSON") == 0)||
         (strcmp(string,"Poisson") == 0) || (strcmp(string,"P") == 0) ||
         (strcmp(string,"p") == 0) ) {
      ppt2->phi_eq = poisson;
    }
    else if ( (strcmp(string,"longitudinal") == 0) || (strcmp(string,"LONGITUDINAL") == 0)||
              (strcmp(string,"Longitudinal") == 0) || (strcmp(string,"L") == 0) ||
              (strcmp(string,"l") == 0) ) {
      ppt2->phi_eq = longitudinal;
    }
    else if ( (strcmp(string,"huang") == 0) || (strcmp(string,"HUANG") == 0)||
              (strcmp(string,"Huang") == 0) || (strcmp(string,"H") == 0) ||
              (strcmp(string,"h") == 0) || (strcmp(string,"Z") == 0) ) {
      ppt2->phi_eq = huang;
      ppt->phi_eq = huang;
    }
    else {
      class_stop (errmsg,
      "phi_equation=%s not supported, choose between poisson, longitudinal and huang",
      string);
    }
  }
  
  /* Set the ppt2->lm_extra parameter according to the gauge */
  if (ppt->gauge == newtonian)
    ppt2->lm_extra = 1;

  else if (ppt->gauge == synchronous)
    ppt2->lm_extra = 3;       
  

  // ================================================================================
  // =                            Perturbations, output                             =
  // ================================================================================

  class_read_int("perturbations2_verbose",ppt2->perturbations2_verbose);

  class_read_int("file_verbose", ppt2->file_verbose);


  /* Read values of (k1,k2,k3) for which to write output files */

  class_call(parser_read_list_of_doubles(
               pfc,
              "k1_out",
              &(int1),
              &(pointer1),
              &flag1,
              errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {

    ppt2->k_out_size = int1;
    
    class_calloc (ppt2->k1_out, ppt2->k_out_size, sizeof(double), errmsg);
    class_calloc (ppt2->k2_out, ppt2->k_out_size, sizeof(double), errmsg);
    class_calloc (ppt2->k3_out, ppt2->k_out_size, sizeof(double), errmsg);

    for (i=0; i<int1; i++)
      ppt2->k1_out[i] = pointer1[i];
    
    free (pointer1);

  }

  class_call(parser_read_list_of_doubles(
               pfc,
              "k2_out",
              &(int1),
              &(pointer1),
              &flag1,
              errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {
    
    class_test(int1 != ppt2->k_out_size, errmsg,
      "specify the same number of values in k1_out, k2_out and k3_out",
      int1);
    
    for (i=0; i<int1; i++)
      ppt2->k2_out[i] = pointer1[i];
    
    free (pointer1);

  }

  class_call(parser_read_list_of_doubles(
               pfc,
              "k3_out",
              &(int1),
              &(pointer1),
              &flag1,
              errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {
    
    class_test(int1 != ppt2->k_out_size, errmsg,
      "specify the same number of values in k1_out, k1_out and k3_out",
      int1);
    
    for (i=0; i<int1; i++)
      ppt2->k3_out[i] = pointer1[i];
    
    free (pointer1);

    ppt2->only_k1k2 = _FALSE_;
  }
  else{
    ppt2->only_k1k2 = _TRUE_;
  }


  /* Check that the added k values satisfy the triangular condition */
  if (ppt2->only_k1k2 == _FALSE_){

    for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {

      double k1 = ppt2->k1_out[index_k_out];
      double k2 = ppt2->k2_out[index_k_out];
      double k3 = ppt2->k3_out[index_k_out];

      class_test (((k1+k2)<k3) || (fabs(k1-k2)>k3), errmsg,
                  "the k_out triplet #%d (%g,%g,%g) does not satisfy the triangular condition",
                  index_k_out, k1, k2, k3);
    }
  }


  /* Make sure that the user specified some k if running in k_out_mode */

  class_test ((ppt2->k_out_mode == _TRUE_) && (ppt2->k_out_size <= 0),
    errmsg,
    "you asked to run SONG in k_out_mode (output=k_out) but you did not specify\
 any output k value. Either turn off the k_out_mode or fill the k1_out, k2_out\
 and k3_out parameters");

  /* In the symmetric k-sampling, we apply a transformation to each (k1,k2,k3)
  triplet; we do the same to the output values, in order to get a match between
  what is asked by the user with k_out and the wavemodes being evolved in SONG. */
    
  if (ppt2->k3_sampling == sym_k3_sampling) {
    
    if (!ppt2->k_out_mode) {

      for (int index_k_out=0; index_k_out < ppt2->k_out_size; ++index_k_out) {

        double KT[4] = {0, ppt2->k1_out[index_k_out], ppt2->k2_out[index_k_out], ppt2->k3_out[index_k_out]};
        double K[4];
  
        class_call (symmetric_sampling_inverse (KT, K, errmsg),
          errmsg, errmsg);

        ppt2->k1_out[index_k_out] = K[1];
        ppt2->k2_out[index_k_out] = K[2];
        ppt2->k3_out[index_k_out] = K[3];

        /* Debug: print the output value before and after */
        // printf ("KT[1]=%g, KT[2]=%g, KT[3]=%g\n", KT[1], KT[2], KT[3]);
        // printf ("K [1]=%g, K [2]=%g, K [3]=%g\n", K [1], K [2], K [3]);
        // printf ("\n");

      }
    }
    
    /* The k_out mode makes SONG compute only the k-values specified in the
    k1_out, k2_out and k3_out parameters. This might create problems with
    interpolation if less than 2 output values are requested. */
    
    else {
      
      class_test (ppt2->k_out_size < 2,
        errmsg,
        "to use the sym_sampling and the k_out_mode together, please provide at least\
 two k-output values, lest the k-interpolation fails");
      
      printf ("\nWARNING: using k_out_mode with the symmetric k-sampling yields very unprecise\
 results unless you choose many output k-values, because the symmetric sampling relies on\
 interpolation over k\n\n");
      
    }
  }

  /* Read values of index_k1, index_k2 and index_k3 for which to write output files.
  If SONG is running in k_out_mode, we ignore these settings, lest there is no k
  to compute */
    
  if (ppt2->k_out_mode == _FALSE_) {

    class_call(parser_read_list_of_integers(
                 pfc,
                "k1_index_out",
                &(int1),
                &(int_pointer1),
                &flag1,
                errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {
        
      ppt2->k_index_out_size = int1;

      class_realloc (ppt2->k1_out, ppt2->k1_out, (ppt2->k_out_size+ppt2->k_index_out_size)*sizeof(double), errmsg);
      class_realloc (ppt2->k2_out, ppt2->k2_out, (ppt2->k_out_size+ppt2->k_index_out_size)*sizeof(double), errmsg);
      class_realloc (ppt2->k3_out, ppt2->k3_out, (ppt2->k_out_size+ppt2->k_index_out_size)*sizeof(double), errmsg);

      class_calloc (ppt2->k1_index_out, ppt2->k_index_out_size, sizeof(int), errmsg);
      class_calloc (ppt2->k2_index_out, ppt2->k_index_out_size, sizeof(int), errmsg);
      class_calloc (ppt2->k3_index_out, ppt2->k_index_out_size, sizeof(int), errmsg);

      for (i=0; i<int1; i++)
        ppt2->k1_index_out[i] = int_pointer1[i];
    
      free (int_pointer1);

    }

    class_call(parser_read_list_of_integers(
                 pfc,
                "k2_index_out",
                &(int1),
                &(int_pointer1),
                &flag1,
                errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {
    
      class_test(int1 != ppt2->k_index_out_size, errmsg,
        "specify the same number of values in k1_index_out, k2_index_out and k3_index_out",
        int1);
    
      for (i=0; i<int1; i++)
        ppt2->k2_index_out[i] = int_pointer1[i];
    
      free (int_pointer1);

    }

    class_call(parser_read_list_of_integers(
                 pfc,
                "k3_index_out",
                &(int1),
                &(int_pointer1),
                &flag1,
                errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {
    
      class_test(int1 != ppt2->k_index_out_size, errmsg,
        "specify the same number of values in k1_index_out, k1_index_out and k3_index_out",
        int1);
    
      for (i=0; i<int1; i++)
        ppt2->k3_index_out[i] = int_pointer1[i];
    
      free (int_pointer1);

    }

  }


  /* Should we output the quadratic sources as well? */

  class_call(parser_read_string(pfc,"output_quadratic_sources",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->output_quadratic_sources = _TRUE_;


  /* Create and open output files for the desired k-values */

  int k_out_total_size = ppt2->k_out_size + ppt2->k_index_out_size;

  if (k_out_total_size > 0) {

    class_calloc (ppt2->index_k1_out, k_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->index_k2_out, k_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->index_k3_out, k_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->k_out_data_byte, k_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->k_out_was_swapped, k_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->k_out_files, k_out_total_size, sizeof(FILE *), errmsg);
    class_calloc (ppt2->k_out_files_sources, k_out_total_size, sizeof(FILE *), errmsg);
    class_calloc (ppt2->k_out_paths, k_out_total_size*_FILENAMESIZE_, sizeof(char), errmsg);
    class_calloc (ppt2->k_out_paths_sources, k_out_total_size*_FILENAMESIZE_, sizeof(char), errmsg);

    if (ppt2->output_quadratic_sources) {
      class_calloc (ppt2->k_out_paths_quad, k_out_total_size*_FILENAMESIZE_, sizeof(char), errmsg);
      class_calloc (ppt2->k_out_files_quad, k_out_total_size, sizeof(FILE *), errmsg);
    }

    for (int index_k_out=0; index_k_out < k_out_total_size; ++index_k_out) {

      /* Build ASCII filenames */
      sprintf (ppt2->k_out_paths[index_k_out],
        "%s/perturbations_song_k%03d.txt",
        pop->root,
        index_k_out);

      /* Open ASCII files */
      class_open(ppt2->k_out_files[index_k_out],
        ppt2->k_out_paths[index_k_out],
        "w",
        errmsg);

      /* Build binary filenames, but do not open the files yet */
      sprintf (ppt2->k_out_paths_sources[index_k_out],
        "%s/sources_song_k%03d.dat",
        pop->root,
        index_k_out);
        
      if (ppt2->output_quadratic_sources == _TRUE_) {
    
        /* Build ASCII filenames */
        sprintf (ppt2->k_out_paths_quad[index_k_out],
          "%s/quadsources_song_k%03d.txt",
          pop->root,
          index_k_out);

        /* Open ASCII files */
        class_open(ppt2->k_out_files_quad[index_k_out],
          ppt2->k_out_paths_quad[index_k_out],
          "w",
          errmsg);
      }


      /* Swap k1 and k2 if the user asked for configurations with k1<k2 */

      ppt2->k_out_was_swapped[index_k_out] = _FALSE_;

      sprintf (ppt2->k_out_swap_message,
        "NOTE: You asked for a (k1,k2) pair with k1<k2, but SONG only computes configurations with\
 k1>=k2. This file contains the perturbations with k1 and k2 swapped; multiply them by a factor\
 (-1)^m to obtain the (k1,k2) pair you asked for. For scalar modes, m=0, there is no need for the\
 factor. For details, please refer to Sec. B.2 of http://arxiv.org/abs/1405.2280.");

      if ((index_k_out < ppt2->k_out_size) && (ppt2->k1_out[index_k_out] < ppt2->k2_out[index_k_out])) {

        ppt2->k_out_was_swapped[index_k_out] = _TRUE_;

        double swap = ppt2->k1_out[index_k_out];
        ppt2->k1_out[index_k_out] = ppt2->k2_out[index_k_out];
        ppt2->k2_out[index_k_out] = swap;

      }

      if ((index_k_out >= ppt2->k_out_size) && (ppt2->k1_index_out[index_k_out-ppt2->k_out_size] < ppt2->k2_index_out[index_k_out-ppt2->k_out_size])) {

        ppt2->k_out_was_swapped[index_k_out] = _TRUE_;

        int swap = ppt2->k1_index_out[index_k_out-ppt2->k_out_size];
        ppt2->k1_index_out[index_k_out-ppt2->k_out_size] = ppt2->k2_index_out[index_k_out-ppt2->k_out_size];
        ppt2->k2_index_out[index_k_out-ppt2->k_out_size] = swap;

      }
    } // for k_out
    

    /* Make sure that k3_out is given in ascending order for the output triplets with
    the same k1_out and k2_out, lest we mess up the addition of the k3_out values
    to the k3 grid in perturb2_get_k_lists() */
  
    for (int i=0; i < ppt2->k_out_size; ++i) {

      for (int j=i; j < ppt2->k_out_size; ++j) {

        if ((ppt2->k1_out[i] == ppt2->k1_out[j]) && (ppt2->k2_out[i] == ppt2->k2_out[j])) {
        
          class_test (ppt2->k3_out[i] > ppt2->k3_out[j],
            errmsg, "k3_out[%d] is larger than k3_out[%d]; triplets with equal (k1_out, k2_out) should\
 be given in ascending k3 order, lest the output files do not match the input k_out values", i, j);

          if ((ppt2->k3_out[i] == ppt2->k3_out[j]) && (i!=j))
            printf ("NOTE: output triplets #%d and #%d are equal; perturbation outputs for #%d will be empty\n", i, j, j);

        }
      }
    }

    
    /* Tell CLASS to output perturbations, too, unless otherwise specified with
    the flag output_class_perturbations */

    ppt->store_perturbations = _TRUE_;
    pop->write_perturbations = _TRUE_;

    class_call(parser_read_string(pfc,"output_class_perturbations",&string1,&flag1,errmsg),
  	     errmsg,
  	     errmsg);	

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))) {
      ppt2->output_class_perturbations = _FALSE_;
      ppt->store_perturbations = _FALSE_;
      pop->write_perturbations = _FALSE_;
    }

  } // if k_out


  /* Read values of tau for which to write output files */

  class_call(parser_read_list_of_doubles(
               pfc,
              "tau_out",
              &(int1),
              &(pointer1),
              &flag1,
              errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {
    
    ppt2->tau_out_size = int1;

    class_calloc (ppt2->tau_out, ppt2->tau_out_size, sizeof(double), errmsg);

    for (i=0; i<int1; i++)
      ppt2->tau_out[i] = pointer1[i];
    
    free (pointer1);

  }

  /* Read values of redshift for which to write output files */

  class_call(parser_read_list_of_doubles(
               pfc,
              "z_out",
              &(int1),
              &(pointer1),
              &flag1,
              errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {
        
    ppt2->z_out_size = int1;

    class_realloc (ppt2->tau_out, ppt2->tau_out, (ppt2->tau_out_size + ppt2->z_out_size)*sizeof(double), errmsg);
    class_calloc (ppt2->z_out, ppt2->z_out_size, sizeof(double), errmsg);

    for (i=0; i<int1; i++)
      ppt2->z_out[i] = pointer1[i];
    
    free (pointer1);

  }
  
  
  int tau_out_total_size = ppt2->tau_out_size + ppt2->z_out_size;
  
  if (tau_out_total_size > 0) {
        
    /* Allocate output file arrays */
    
    class_calloc (ppt2->index_tau_out, tau_out_total_size, sizeof(int), errmsg);
    class_calloc (ppt2->tau_out_was_reduced, tau_out_total_size, sizeof(int), errmsg);

    class_alloc (ppt2->tau_out_paths, k_out_total_size*sizeof(char*), errmsg);
    class_alloc (ppt2->tau_out_files, k_out_total_size*sizeof(FILE**), errmsg);
    for (int index_k=0; index_k < k_out_total_size; ++index_k) {
      class_alloc (ppt2->tau_out_paths[index_k], tau_out_total_size*_FILENAMESIZE_*sizeof(char), errmsg);
      class_alloc (ppt2->tau_out_files[index_k], tau_out_total_size*sizeof(FILE*), errmsg);
    }

    class_calloc (ppt2->tau_out_paths_sources, tau_out_total_size*_FILENAMESIZE_, sizeof(char), errmsg);
    class_calloc (ppt2->tau_out_files_sources, tau_out_total_size, sizeof(FILE *), errmsg);


    /* Create and open output files for the desired tau values */

    for (int index_tau_out=0; index_tau_out < ppt2->tau_out_size; ++index_tau_out) {

      for (int index_k_out=0; index_k_out < k_out_total_size; ++index_k_out) {

        /* Build ASCII filenames */
        sprintf (ppt2->tau_out_paths[index_k_out][index_tau_out],
          "%s/perturbations_song_k%03d_tau%03d.txt",
          pop->root,
          index_k_out,
          index_tau_out);

        /* Open ASCII files, but only if some k-values are requested */
        class_open(ppt2->tau_out_files[index_k_out][index_tau_out],
          ppt2->tau_out_paths[index_k_out][index_tau_out],
          "w",
          errmsg);
        
      }

      /* Build binary filenames */
      sprintf (ppt2->tau_out_paths_sources[index_tau_out],
        "%s/sources_song_tau%03d.dat",
        pop->root,
        index_tau_out);

    }

    /* For redshift values we use different file names but same data structures */

    for (int index_z_out=0; index_z_out < ppt2->z_out_size; ++index_z_out) {

      for (int index_k_out=0; index_k_out < k_out_total_size; ++index_k_out) {

        /* Build ASCII filenames */
        sprintf (ppt2->tau_out_paths[index_k_out][ppt2->tau_out_size+index_z_out],
          "%s/perturbations_song_k%03d_z%03d.txt",
          pop->root,
          index_k_out,
          index_z_out);

        /* Open ASCII files, but only if some k-values are requested */
        class_open(ppt2->tau_out_files[index_k_out][ppt2->tau_out_size+index_z_out],
          ppt2->tau_out_paths[index_k_out][ppt2->tau_out_size+index_z_out],
          "w",
          errmsg);

      }

      /* Build binary filenames */
      sprintf (ppt2->tau_out_paths_sources[ppt2->tau_out_size+index_z_out],
        "%s/sources_song_z%03d.dat",
        pop->root,
        index_z_out);
    }

  } // if (tau_out)
  
  

  // =================================================================================
  // =                                Bessel functions                               =
  // =================================================================================

  class_read_int("bessels2_verbose",pbs2->bessels2_verbose);

  /* Make sure that we compute all the needed Bessel functions for the bispectrum */
  if (ppt2->has_cmb_bispectra) 
    pbs2->extend_l1_using_m = _TRUE_;

  /* Minimum x treshold for the spherical Bessels with j_l1(x). These are the
  Bessels that are summed to obtain J_Llm */
  class_read_double("bessel_j_cut_2nd_order", ppr2->bessel_j_cut_song); /* obsolete */
  class_read_double("bessel_j_cut_song", ppr2->bessel_j_cut_song);

  /* Minimum treshold for the functions J_Llm(x).  These are obtained as a
  weighted sum of spherical Bessels j_l1(x) with |L-l| <= l1 <= L+l */
  class_read_double("bessel_J_cut_2nd_order", ppr2->bessel_J_cut_song); /* obsolete */
  class_read_double("bessel_J_cut_song", ppr2->bessel_J_cut_song);

  /* Linear step dx where we are going to sample the j_l1(x) and J_Llm(x) */
  class_read_double("bessel_x_step_2nd_order", ppr2->bessel_x_step_song); /* obsolete */
  class_read_double("bessel_x_step_song", ppr2->bessel_x_step_song);
  

  // =========================================================================================
  // =                                  Transfer functions                                   =
  // =========================================================================================

  class_read_int("transfer2_verbose", ptr2->transfer2_verbose);

  /* - k sampling */
  class_call(parser_read_string(pfc,"transfer2_k_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg); /* obsolete */
  class_call(parser_read_string(pfc,"transfer2_k3_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if (strstr(string1,"bessel") != NULL)
      ptr2->k_sampling = bessel_k3_sampling;

    else if ((strstr(string1,"class") != NULL) || (strstr(string1,"smart") != NULL))
      ptr2->k_sampling = class_transfer2_k3_sampling;
    
    else
      class_stop(errmsg,
        "transfer2_k3_sampling=%s not supported, choose between 'bessel', 'smart' and 'class'.", string1);
  }

  /* If 'transfer2_k3_sampling=class', choose the density of the k3-sampling for
  the transfer functions */
  class_read_double("k_step_trans_scalars_2nd_order", ppr2->q_linstep_song); /* obsolete */
  class_read_double("q_linstep_song",ppr2->q_linstep_song);

  /* Older versions of SONG used the parameter 'k_step_trans_scalars_2nd_order' instead of
  'q_linstep_song' to specify the frequency of the k-sampling for the second-order transfer
  functions. The change is not only in name but in substance, because the two parameters
  have a different definition:
    q_linstep / k_step_trans_scalars = (tau0-tau_rec)/pth->rs_rec ~ 95.
  Here we apply the corrective factor, hoping that the user will update his parameter file
  with q_linstep. For more detail, please refer to the comment in transfer.c on top of the
  definition of q_period.*/

  class_call(parser_read_double(pfc,"k_step_trans_scalars_2nd_order",&param1,&flag1,errmsg),
    errmsg,
    errmsg);

  if (flag1==_TRUE_) {

    ppr2->old_run = _TRUE_;
  
    /* The estimate (tau0-tau_rec)/pth->rs_rec ~ 95 is not exact, so there will be
    a slight difference between the old and the new k-sampling.  For this reason,
    we modify the parameter only if not loading from disk an older run, because
    in that case any difference in the k-sampling will result in errors. */
    if (ppr->load_run == _FALSE_) {
        double ratio = 95;
        ppr2->q_linstep_song = param1 * ratio;
        printf ("\nOBSOLETE PARAMETER: changed k_step_trans_scalars_2nd_order=%g to q_linstep_song=%g (ratio=%g).\n\n",
          param1, ppr2->q_linstep_song, ratio);
    }
  }
  

  /* - time sampling */
  class_call(parser_read_string(pfc,"transfer2_tau_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"sources") != NULL) || (strstr(string1,"source") != NULL))
      ptr2->tau_sampling = sources_tau_sampling;

    else if ((strstr(string1,"custom") != NULL) || (strstr(string1,"smart") != NULL))
      ptr2->tau_sampling = custom_tau_sampling;

    else if (strstr(string1,"bessel") != NULL)
      ptr2->tau_sampling = bessel_tau_sampling;
    
    else
      class_stop(errmsg,
        "transfer2_tau_sampling=%s not supported, choose between 'bessel', 'smart' or 'custom'.", string1);
  }

  /* Specify the density for the time sampling of the transfer function integral. Used only
  if If transfer2_tau_sampling=smart. Older versions of SONG used the parameter tau_step_trans_song,
  which is related to the new one by a 2*pi factor. */

  class_read_double("tau_step_trans_2nd_order", ppr2->tau_linstep_song); /* obsolete */
  if (flag1 == _TRUE_)
    ppr2->tau_linstep_song /= 2*_PI_;
  class_read_double("tau_step_trans_song", ppr2->tau_linstep_song); /* obsolete */
  if (flag1 == _TRUE_)
    ppr2->tau_linstep_song /= 2*_PI_;
  class_read_double("tau_linstep_song", ppr2->tau_linstep_song);
  
  /* If the user sets a negative value, do not add points to the time grid */
  if (ppr2->tau_linstep_song < 0)
    ppr2->tau_linstep_song = _HUGE_;

  

  // =============================================================================================
  // =                                       Bispectra                                           =
  // =============================================================================================

  /* Should we include the quadratic corrections to the intrinsic bispectrum coming from
  the bolometric temperature and from the redshift term? */
  class_call(parser_read_string(pfc,"add_quadratic_correction",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);	

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))) {
    pbi->add_quadratic_correction = _FALSE_;
  }



  // =============================================================================================
  // =                                     Fisher matrices                                       =
  // =============================================================================================



  

  // =============================================================================================
  // =                                      File output                                          =
  // =============================================================================================

  // ----------------------------------------------------------------------------------------
  // -                               Disk storage of sources                                -
  // ----------------------------------------------------------------------------------------

  /* Store to disk the second-order line of sight sources? */
  class_call(parser_read_string(pfc,"store_sources",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);
      
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppr2->store_sources_to_disk = _TRUE_;

  sprintf(ppt2->sources_dir, "%s/sources", ppr->data_dir);

  /* If we are not loading from disk, just create the source directory */
  if ((ppr2->store_sources_to_disk == _TRUE_) && (ppr->load_run == _FALSE_)) {
    
    class_test (mkdir (ppt2->sources_dir, 0777) != 0,
      errmsg,
      "could not create directory '%s', maybe it already exists?", ppt2->sources_dir);
  }
  /* If we are in a run directory, checks if it already contains the source functions */
  else if (ppr->load_run == _TRUE_) {

    struct stat st;
    short sources_dir_exists = (stat(ppt2->sources_dir, &st)==0);

    /* If the sources directory exists, then we shall load the 2nd-order source functions from it */
    if (sources_dir_exists) {
      ppr2->store_sources_to_disk = _FALSE_;
      ppr2->load_sources_from_disk = _TRUE_;
      if (ppt2->perturbations2_verbose > 1)
        printf (" -> found source functions folder in run directory.\n");
    }
    /* Otherwise, create it */
    else if (ppr2->store_sources_to_disk == _TRUE_) {
              
      if (ppt2->perturbations2_verbose > 1)
        printf (" -> source functions folder not found in run directory, will create it.\n");

      class_test (mkdir (ppt2->sources_dir, 0777)!=0,
        errmsg,
        "could not create directory '%s', maybe it already exists?", ppt2->sources_dir);
        
      ppr2->load_sources_from_disk = _FALSE_;
    }
  }

  /* Create/open the status file. The 'a+' mode means that if the file does not exist it will be created,
  but if it exist it won't be erased (append mode) */
  // if (ppr2->store_sources_to_disk == _TRUE_) {
  //   sprintf(ppt2->sources_status_path, "%s/sources_status_file.txt", ppr->data_dir);
  //   class_open(ppt2->sources_status_file, ppt2->sources_status_path, "a+", errmsg);
  // }

  class_test ((ppr2->store_sources_to_disk == _TRUE_) && (ppr2->load_sources_from_disk == _TRUE_),
    errmsg,
    "cannot load and save sources at the same time!");
    

  // ----------------------------------------------------------------------------------------
  // -                               Disk storage of transfers                              -
  // ----------------------------------------------------------------------------------------

  /* Store to disk the second-order transfer functions? */
  class_call(parser_read_string(pfc,"store_transfers",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);
      
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppr2->store_transfers_to_disk = _TRUE_;

  sprintf(ptr2->transfers_dir, "%s/transfers", ppr->data_dir);

  /* If we are not loading from disk, just create the transfer directory */
  if ((ppr2->store_transfers_to_disk == _TRUE_) && (ppr->load_run == _FALSE_)) {
    
    class_test (mkdir (ptr2->transfers_dir, 0777) != 0,
      errmsg,
      "could not create directory '%s', maybe it already exists?", ptr2->transfers_dir);
  }
  /* If we are in a run directory, checks if it already contains the transfer functions */
  else if (ppr->load_run == _TRUE_) {

    struct stat st;
    short transfers_dir_exists = (stat(ptr2->transfers_dir, &st)==0);

    /* If the transfers directory exists, then we shall load the 2nd-order transfer functions from it */
    if (transfers_dir_exists) {
      ppr2->store_transfers_to_disk = _FALSE_;
      ppr2->load_transfers_from_disk = _TRUE_;
      if (ptr2->transfer2_verbose > 1)
        printf (" -> found transfer functions folder in run directory.\n");
    }
    /* Otherwise, create it */
    else if (ppr2->store_transfers_to_disk == _TRUE_) {
              
      if (ptr2->transfer2_verbose > 1)
        printf (" -> transfer functions folder not found in run directory, will create it.\n");

      class_test (mkdir (ptr2->transfers_dir, 0777)!=0,
        errmsg,
        "could not create directory '%s', maybe it already exists?", ptr2->transfers_dir);
        
      ppr2->load_transfers_from_disk = _FALSE_;
    }
  }

  /* Create/open the status file. The 'a+' mode means that if the file does not exist it will be created,
  but if it exist it won't be erased (append mode) */
  // if (ppr2->store_transfers_to_disk == _TRUE_) {
  //   sprintf(ptr2->transfers_status_path, "%s/transfers_status_file.txt", ppr->data_dir);
  //   class_open(ptr2->transfers_status_file, ptr2->transfers_status_path, "a+", errmsg);
  // }

  class_test ((ppr2->store_transfers_to_disk == _TRUE_) && (ppr2->load_transfers_from_disk == _TRUE_),
    errmsg,
    "cannot load and save transfers at the same time!");
    

  // =============================================================================================
  // =                                  Interpolation techniques                                 =
  // =============================================================================================
  
  class_call(parser_read_string(pfc,"sources_time_interpolation",&string1,&flag1,errmsg),
       errmsg,
       errmsg); 

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr2->sources_time_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL)
    || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr2->sources_time_interpolation = cubic_interpolation;
    
    else
      class_stop (errmsg,         
        "sources_time_interpolation=%s not supported, choose between 'linear', 'cubic' and 'spline'",
        string1);

  }

  class_call(parser_read_string(pfc,"sources_k3_interpolation",&string1,&flag1,errmsg),
       errmsg,
       errmsg); 

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr2->sources_k3_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL)
    || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr2->sources_k3_interpolation = cubic_interpolation;
    
    else
      class_stop(errmsg,         
        "sources_k3_interpolation=%s not supported, choose between 'linear', 'cubic' and 'spline'",
        string1);
  }


  // =========================================================================================
  // =                                   Angular scales                                      =
  // =========================================================================================

  /* Read l_max for the Boltzmann hierarchies */
  class_read_int("l_max_g_2nd_order", ppr2->l_max_g); /* obsolete */
  class_read_int("l_max_pol_g_2nd_order", ppr2->l_max_pol_g); /* obsolete */
  class_read_int("l_max_ur_2nd_order", ppr2->l_max_ur);   /* obsolete */
  
  class_read_int("l_max_g_song", ppr2->l_max_g);
  class_read_int("l_max_pol_g_song", ppr2->l_max_pol_g);    
  class_read_int("l_max_ur_song", ppr2->l_max_ur);        

  /* Compute the maximum l_max for the Boltzmann hierarchy */
  ppr2->l_max_boltzmann = 0;  
  ppr2->l_max_boltzmann = MAX (ppr2->l_max_boltzmann, ppr2->l_max_g);
  if (ppt2->has_polarization2 == _TRUE_)
    ppr2->l_max_boltzmann = MAX (ppr2->l_max_boltzmann, ppr2->l_max_pol_g);
  if (pba->Omega0_ur!=0)
    ppr2->l_max_boltzmann = MAX (ppr2->l_max_boltzmann, ppr2->l_max_ur);


  /* Read l_max for the quadratic sources in the Boltzmann hierarchies. If the user specified
  a negative value for one of them, set it to the corresponding l_max_XXX_song (see above).
  Also make sure that each of them is not larger than the corresponding l_max_XXX_song. The
  following lines must go below the definitions of l_max_g_song, etc. */
  class_read_int("l_max_g_quadsources", ppr2->l_max_g_quadsources);
  if ((ppr2->l_max_g_quadsources<0) || (ppr2->l_max_g_quadsources>ppr2->l_max_g))
    ppr2->l_max_g_quadsources = ppr2->l_max_g;

  class_read_int("l_max_pol_g_quadsources", ppr2->l_max_pol_g_quadsources);
  if ((ppr2->l_max_pol_g_quadsources<0) || (ppr2->l_max_pol_g_quadsources>ppr2->l_max_pol_g))
    ppr2->l_max_pol_g_quadsources = ppr2->l_max_pol_g;

  class_read_int("l_max_ur_quadsources", ppr2->l_max_ur_quadsources);
  if ((ppr2->l_max_ur_quadsources<0) || (ppr2->l_max_ur_quadsources>ppr2->l_max_ur))
    ppr2->l_max_ur_quadsources = ppr2->l_max_ur;

  /* Read l_max for the line of sight integration */
  class_read_int("l_max_T_los", ppr2->l_max_los_t); /* obsolete */
  class_read_int("l_max_E_los", ppr2->l_max_los_p); /* obsolete */
  class_read_int("l_max_B_los", ppr2->l_max_los_p); /* obsolete */
  class_read_int("l_max_los_t", ppr2->l_max_los_t);
  class_read_int("l_max_los_p", ppr2->l_max_los_p);

  /* Same as above, but for the quadratic terms */
  class_read_int("l_max_T_quadratic_los", ppr2->l_max_los_quadratic_t); /* obsolete */
  class_read_int("l_max_E_quadratic_los", ppr2->l_max_los_quadratic_p); /* obsolete */
  class_read_int("l_max_B_quadratic_los", ppr2->l_max_los_quadratic_p); /* obsolete */
  class_read_int("l_max_los_quadratic_t", ppr2->l_max_los_quadratic_t);
  class_read_int("l_max_los_quadratic_p", ppr2->l_max_los_quadratic_p);

  /* Compute the maximum l_max for the line of sight integration */
  ppr2->l_max_los = 0;
  ppr2->l_max_los_quadratic = 0;
  
  if (ppt2->has_cmb_temperature == _TRUE_) {
    ppr2->l_max_los = MAX (ppr2->l_max_los, ppr2->l_max_los_t);
    ppr2->l_max_los_quadratic = MAX (ppr2->l_max_los_quadratic, ppr2->l_max_los_quadratic_t);
  }

  if ((ppt2->has_cmb_polarization_e == _TRUE_) || (ppt2->has_cmb_polarization_b == _TRUE_)) {
    ppr2->l_max_los = MAX (ppr2->l_max_los, ppr2->l_max_los_p);
    ppr2->l_max_los_quadratic = MAX (ppr2->l_max_los_quadratic, ppr2->l_max_los_quadratic_p);
  }


  // =======================================================================================
  // =                                 Azimuthal modes                                     =
  // =======================================================================================

  /* Read the list of requested azimuthal 'm' values. They must be given in ascending order. */

  class_call (parser_read_list_of_integers (
    pfc,"modes_2nd_order",&(ppr2->m_size),&(int_pointer1),&flag1, errmsg),
    errmsg,
    errmsg); /* obsolete */

  class_call (parser_read_list_of_integers (
    pfc,"modes_song",&(ppr2->m_size),&(int_pointer2),&flag2, errmsg),
    errmsg,
    errmsg);

  /* string2 wins over string1 */
  flag = _TRUE_;
  if (flag2 == _TRUE_)
    int_pointer = int_pointer2;
  else if (flag1 == _TRUE_)
    int_pointer = int_pointer1;
  else
    flag = _FALSE_;

  if (flag == _TRUE_) {
    for (i=0; i < ppr2->m_size; ++i)
      ppr2->m[i] = int_pointer[i];
  }

  if (flag1 == _TRUE_)
    free (int_pointer1);
  else if (flag2 == _TRUE_)
    free (int_pointer2);

  /* Check that the m-list is strictly ascending */
  for (i=0; i < (ppr2->m_size-1); ++i)
    class_test (ppr2->m[i] >= ppr2->m[i+1],
      errmsg,
      "the m-list provided  in modes_song is not strictly ascending");

  /* Maximum 'm' that will be computed */
  ppr2->m_max_song = ppr2->m[ppr2->m_size-1];

  if (ppt2->has_cmb_polarization_b && ppr2->m_max_song == 0)
    printf ("\nWARNING: will ignore any requested B-mode output. Second-order B-modes vanish\
 for scalar perturbations; to obtain a non-vanishing B-mode result, please include in the\
 modes_song parameter at least one value larger than zero\n\n");
  
  /* Check that the m's are positive */
  class_test (ppr2->m[0] < 0,
    errmsg,
    "the m-list provided in modes_song has negative numbers in it");

  /* Check that m_max is smaller than limit */  
  class_test (ppr2->m_max_song > (_MAX_NUM_AZIMUTHAL_-1),
    errmsg,
    "the maximum value of m cannot exceed %d, please choose modes_song\
 accordingly or increase the macro _MAX_NUM_AZIMUTHAL_",
    _MAX_NUM_AZIMUTHAL_);

  /* Find out the index in ppr2->m corresponding to a given m. */
	for(int m=0; m<=ppr2->m_max_song; ++m) {

		ppr2->index_m[m] = -1;

		for (int index_m=0; index_m<ppr2->m_size; ++index_m)
			if (m==ppr2->m[index_m]) ppr2->index_m[m] = index_m;
	}

	/* Build a logical array to check whether a given m is in m_vec.
  First, intialise it to _FALSE_. */
	for(int m=0; m<_MAX_NUM_AZIMUTHAL_; ++m)
    ppr2->compute_m[m] = _FALSE_;

	for(int m=0; m<=ppr2->m_max_song; ++m)
		for (int index_m=0; index_m<ppr2->m_size; ++index_m)
			if (m==ppr2->m[index_m]) ppr2->compute_m[m] = _TRUE_;

  /* The largest l that will be needed by SONG. The first term (pbs->l_max) is the maximum
  multipole we shall compute the spectra and bispectra in. The second term represent
  the additional l's where we need to the Bessel functions in order to compute the projection
  functions (ppr2->l_max_los) and the bispectrum integral (ppr2->m_max_song). This line
  must be after setting ppr2->l_max_los and ppr2->l_max_boltzmann. */
  int l_max = pbs->l_max + MAX (ppr2->l_max_los, ppr2->m_max_song);
  
  /* In some debug scenarios, we will need to evolve a lot of Boltzmann moments, more than
  the maximum l for the C_l */
  l_max = MAX (l_max, ppr2->l_max_boltzmann);

  /* Determine for each l (starting from 0) the position of its maximum m in ppr2->m. The resulting
  array ppr2->index_m_max will be used throughout the code to cycle through the m's allowed for a given l. */
  class_alloc (ppr2->index_m_max, (l_max+1)*sizeof(int), ppr2->error_message);
  for (int l=0; l<=l_max; ++l) {
   
    /* Ignore the l's that are smaller than the smallest m. It is important that the value here is -1
    because (i) a cycle starting from 0 won't start and (ii) a size is computed as index_max+1, hence
    -1 means 0 size */
    if (l < ppr2->m[0]) {
      ppr2->index_m_max[l] = -1;
      continue;
    }
    
    /* Cycle till you find the largest possible m for this l */
    int index_m = ppr2->m_size-1;
    while (ppr2->m[index_m] > l) index_m--;
    ppr2->index_m_max[l] = index_m;
    
    /* Some debug */
    // printf ("for l=%d, the maximum allowed m in ppr2->m is m=%d\n",
    //   l, ppr2->m[ppr2->index_m_max[l]]);
  }

  /* For certain bispectrum types, the 3j-symbol (l1,l2,l3)(0,0,0) does not appear explicitly
  in the bispectrum formula and therefore it cannot be pulled out analytically. This is the
  case for the intrinsic bispectrum when m>0, or for the CMB-lensing and quadratic bispectra
  in presence of polarisation. To circumvent this issue, we choose to have 
  an l-grid where all the l's are even. This is not completely satisfactory because half of
  the configurations (those with even l1+l2+l3 but two odd components, like 2,3,3 or 2,3,7)
  will be always skipped. The alternative, however, is worse: if the 1D l-grid has both
  even and odd l's, it can happen that all of the configurations are skipped! Think of
  having a grid with step 2 starting from an odd value: you will always get odd l's, which
  means that l1+l2+l3 is always odd too. (A potential solution is to have a step of delta_l=1
  up to some even l and then carry on with an even linear step).
  An exception to this rule is when ppr->l_linstep=1, that is, when we take all l's (even
  and odd) in our 1D l-list. In this case, nothing can be skipped and the Fisher estimator
  will give an exact result. */
  if ((ppr->l_linstep!=1) && (pbi->has_intrinsic==_TRUE_) && (ppr2->m_max_song>0)) {
       
    printf ("\n");
    printf ("   *@^#?!?! FORCING THE COMPUTATION OF A GRID OF EVEN L'S\n");
    printf ("\n");      
    ppr->compute_only_even_ls = _TRUE_;
  
    // printf ("\n");
    // printf ("   *@^#?!?! FORCING THE COMPUTATION OF A GRID OF ODD L'S\n");
    // printf ("\n");
    // ppr->compute_only_odd_ls = _TRUE_;
  
    class_test_permissive (
      ((ppt2->has_cmb_polarization_b==_TRUE_) && (ppr->compute_only_even_ls==_TRUE_)) ||
      ((ppt2->has_cmb_polarization_e==_TRUE_) && (ppr->compute_only_odd_ls==_TRUE_)) ||
      ((ppt2->has_cmb_temperature==_TRUE_) && (ppr->compute_only_odd_ls==_TRUE_)),
      errmsg,
      "careful, your choice of parity is wrong!");
  }
  
  
  // ====================================================================================
  // =                                 Quadratic sources                               =
  // ====================================================================================

  /* Do we want to rescale all CMB multipoles with a factor 1/sin(theta_1)^m? The
  rescaling does not affect the m>0 transfer functions. It is needed to compute the
  intrinsic bispectrum, so by default it is active. See the header file perturbations2.h
  for details on the rescaling. */
  ppt2->rescale_cmb_sources = _TRUE_;
  
  /* Uncomment if you want the output functions to output non-rescaled functions */
  if ((ppt2->stop_at_perturbations2 == _TRUE_) || (ptr2->stop_at_transfers2 == _TRUE_))
    ppt2->rescale_cmb_sources = _FALSE_;

  /* For m=0, the rescaling is equal to unity, so we do not apply it. This is an
  optimisation that you can comment out in case you want to test the rescaling. */
  if (ppr2->m_max_song == 0)
    ppt2->rescale_cmb_sources = _FALSE_;


  // ====================================================================================
  // =                              Compute pbs2->xx_max                                =
  // ====================================================================================

  /* Determine pbs2->xx_max, the upper limit of the x-domain of the projection functions
  J_Llm(x). These appear in the second-order line-of-sight integral with x = k*(tau0-tau),
  therefore we set pbs2->xx_max = k_max*tau0, where tau0 is the conformal age of the
  Universe */
  pbs2->xx_max = ppr2->k_max_tau0_over_l_max * ppt->l_scalar_max;

  /* Accomodate for the case where the k-sampling is set manually (usually for
  debug purposes) */
  if ((ppt2->k_sampling == lin_k_sampling) || (ppt2->k_sampling == log_k_sampling)) {
    double tau0_guess = 18000;
    pbs2->xx_max = ppr2->k_max_custom * tau0_guess;
  }

  /* Accomodate for the case where the user asked to output custom k-values */    
  if (ppt2->k_out_size > 0) {
    double k_max = 0;
    for (int i=0; i < ppt2->k_out_size; ++i)
      k_max = MAX (MAX (MAX (ppt2->k1_out[i],ppt2->k2_out[i]), ppt2->k3_out[i]), k_max);
    double tau0_guess = 18000;
    pbs2->xx_max = MAX (pbs2->xx_max, k_max * tau0_guess);
  }

  /* Copy the step size in xx to the bessel2 structure */
  pbs2->xx_step = ppr2->bessel_x_step_song;

  /* Extend pbs2->xx_max to avoid potential out-of-bounds errors in the interpolation
  of J_Llm(x) */
  pbs2->xx_max += pbs2->xx_step;
  pbs2->xx_max *= 1.05;

  /* Extend the domain of the Bessel functions j_l(x) computed in bessel.c (pbs->x_max) to
  match the domain of the projection functions J_Llm(x) computed in bessel2.c (pbs2->xx_max).
  The reason is that the computation of J_Llm(x) requires j_l(x). */
  pbs->x_max = MAX (pbs->x_max, pbs2->xx_max);

  /* Determine the maximum value of L for which we should compute the projection functions
  J_Llm(x) in the second-order Bessel module. To explain the inclusion of ppr2->m_max_song,
  refer to the long comment in bessel2_get_l1_list(). */
  pbs2->L_max = MAX (ppr2->l_max_los, ppr2->m_max_song);

  return _SUCCESS_;

}


/** 
 * All default parameter values (for input parameters)
 */

int input2_default_params (
       struct background *pba,
       struct thermo *pth,
       struct perturbs *ppt,
       struct perturbs2 *ppt2,       
       struct transfers *ptr,
       struct bessels * pbs,
       struct bessels2 * pbs2,
       struct transfers2 *ptr2,     
       struct primordial *ppm,
       struct spectra *psp,
       struct nonlinear * pnl,
       struct lensing *ple,
       struct bispectra *pbi,
       struct fisher *pfi,
       struct output *pop
       ) {


  // ==============================================================
  // =                   perturbed recombination                  =
  // ==============================================================

  pth->has_perturbed_recombination_stz = _FALSE_;
  pth->perturbed_recombination_turnx = 36;
  ppt->has_perturbed_recombination_stz = _FALSE_;
  ppt2->has_perturbed_recombination_stz = _FALSE_;
  ppt2->perturbed_recombination_use_approx = _FALSE_;
  

  // ===========================================================
  // =                     perturb1 structure                  =
  // ===========================================================

  ppt->gauge = newtonian;

  /* Make sure CLASS does not add points to the k grid  */
  ppt->k_output_values_num = -1;
  ppt->store_perturbations = _FALSE_;
  pop->write_perturbations = _FALSE_;


  // ===========================================================
  // =                     perturb2 structure                  =
  // ===========================================================
  
  /* - Flags */
  
  ppt2->has_cmb_temperature = _FALSE_;
  ppt2->has_cmb_polarization_e = _FALSE_; 
  ppt2->has_cmb_polarization_b = _FALSE_; 
  ppt2->has_pk_delta_cdm = _FALSE_;
  ppt2->has_pk_delta_b = _FALSE_;
  ppt2->has_pk_magnetic = _FALSE_;
  ppt2->has_bk_delta_cdm = _FALSE_;

  ppt2->has_perturbations2 = _FALSE_;
  ppt2->has_polarization2  = _TRUE_;
  ppt2->has_quadratic_sources = _TRUE_;
  ppt2->has_quadratic_liouville = _TRUE_;
  ppt2->has_quadratic_collision = _TRUE_;
  ppt2->has_perfect_baryons = _TRUE_;
  ppt2->has_perfect_cdm = _TRUE_;
  ppt2->rescale_cmb_sources = _TRUE_;
  ppt2->compute_quadsources_derivatives = _FALSE_;

  ppt2->rescale_cmb_sources = _FALSE_;

  ppt2->has_pure_scattering_in_los = _FALSE_;
  ppt2->has_quad_scattering_in_los = _FALSE_;
  ppt2->has_pure_metric_in_los = _FALSE_;
  ppt2->has_quad_metric_in_los = _FALSE_;

  ppt2->has_time_delay_in_los = _FALSE_;
  ppt2->has_redshift_in_los = _FALSE_;
  ppt2->has_lensing_in_los = _FALSE_;

  ppt2->use_delta_tilde_in_los = _FALSE_;

  ppt2->has_sw = _FALSE_;
  ppt2->use_exponential_potentials = _FALSE_;
  ppt2->has_isw = _FALSE_;
  ppt2->only_early_isw = _FALSE_;

	ppt2->use_test_source = _FALSE_;

  ppt2->has_only_recombination = _FALSE_;
  ppt2->has_only_reionisation = _FALSE_;
  ppt2->has_only_loss_term = _FALSE_;
  ppt2->has_only_gain_term = _FALSE_;
  
  ppt2->has_cmb_spectra = _FALSE_;
  ppt2->has_cmb_bispectra = _FALSE_;
  
  ppt2->stop_at_perturbations2 = _FALSE_;
  ppt2->stop_at_perturbations1 = _FALSE_;
  

  /* - Initial conditions */

  ppt2->has_ad = _TRUE_;      
  ppt2->has_zero_ic = _FALSE_;
  ppt2->has_ad_first_order = _FALSE_;
  ppt2->has_unphysical_ic = _FALSE_;
  ppt2->primordial_local_fnl_phi = 0;
  

  /* - Approximations at second order */

  ppt2->tight_coupling_approximation = tca2_first_order_pitrou;
  ppt2->tight_coupling_trigger_tau_c_over_tau_h = 0.01;
  ppt2->tight_coupling_trigger_tau_c_over_tau_k = 0.007;
  
  ppt2->radiation_streaming_approximation = rsa2_MD;
  ppt2->radiation_streaming_trigger_tau_over_tau_k = 90;
  ppt2->radiation_streaming_trigger_tau_c_over_tau = 10;
  
  ppt2->ur_fluid_approximation = ufa2_none;
  ppt2->ur_fluid_trigger_tau_over_tau_k = 15;
  
  ppt2->no_radiation_approximation = nra2_none;
  ppt2->no_radiation_approximation_rho_m_over_rho_r = 100;


  /* - Choose equation for the time potential */

  ppt2->phi_eq = huang;


  /* - l sampling */

  ppt2->lm_extra = 1;
  

  /* - Time sampling */

  ppt2->recombination_max_to_end_ratio = 1000;

  ppt2->has_custom_timesampling = _FALSE_;
  ppt2->custom_tau_ini  = 1;
  ppt2->custom_tau_end  = 0;      // 0 means: integrate all the way to today
  ppt2->custom_tau_size = 2000;
  ppt2->custom_tau_mode = lin_tau_sampling;

  
  /* - k sampling */

  ppt2->k_sampling = smart_sources_k_sampling;
  ppt2->k3_sampling = smart_k3_sampling;
  ppt2->k_max_for_pk = 0.1;

  /* - Technical parameters */

  ppt2->perturbations2_verbose = 0;
  
  ppt2->file_verbose = 1;

  ppt2->tau_out_size = 0;
  ppt2->z_out_size = 0;

  ppt2->k_out_size = 0;
  ppt2->k_index_out_size = 0;
  ppt2->output_class_perturbations = _TRUE_;
  ppt2->output_quadratic_sources = _FALSE_;
  ppt2->k_out_mode = _FALSE_;
          

  // =============================================================
  // =                     bessel2 structure                     =
  // =============================================================

  pbs2->bessels2_verbose = 0;
  pbs2->extend_l1_using_m = _FALSE_;

  
  // ============================================================
  // =                     transfer2 structure                  =
  // ============================================================

  ptr2->transfer2_verbose = 0;
  ptr2->k_sampling = class_transfer2_k3_sampling;
  ptr2->tau_sampling = sources_tau_sampling;
  ptr2->stop_at_transfers2 = _FALSE_;


  // ============================================================
  // =                     bispectra structure                  =
  // ============================================================

  pbi->add_quadratic_correction = _TRUE_;
  
    

  // ========================================================
  // =                    Fisher structure                  =
  // ========================================================



  return _SUCCESS_;

}


/** 
 * Initialize the precision parameter structure. 
 * 
 * All precision parameters used in the other moduels are listed here
 * and assigned here a default value.
 *
 * @param ppr Input/Ouput: a precision_params structure pointer  
 * @return the error status
 *
 */

int input2_default_precision ( struct precision2 * ppr2 ) {


  // ===================================================================
  // =                              Tolerance                          =
  // ===================================================================

  ppr2->tol_perturb_integration_song=1e-4;


  // ==================================================================
  // =                          Multipole limits                      =
  // ==================================================================

  ppr2->m_max_song=0;

  ppr2->l_max_g=8;
  ppr2->l_max_pol_g=8;
  ppr2->l_max_ur=8; 

  ppr2->l_max_g_quadsources=-1;
  ppr2->l_max_pol_g_quadsources=-1;
  ppr2->l_max_ur_quadsources=-1;

  ppr2->l_max_los_t=ppr2->l_max_los_quadratic_t=2;
  ppr2->l_max_los_p=ppr2->l_max_los_quadratic_p=2;

  /* By default, compute only the scalar (m=0) modes */
  ppr2->m_size = 1;
  ppr2->m[0] = 0;


  // ======================================================================
  // =                            Time samplings                          =
  // ======================================================================

  ppr2->custom_tau_start_evolution = 0;
  ppr2->perturb_sampling_stepsize_song = 0.4;
  ppr2->perturb_sampling_late_time_boost = 1;
  ppr2->start_small_k_at_tau_c_over_tau_h_song = 0.0015;  /* decrease to start earlier in time */
  ppr2->start_large_k_at_tau_h_over_tau_k_song = 0.07;    /* decrease to start earlier in time */



  // =====================================================================
  // =                             k samplings                           =
  // =====================================================================

  /* Smart sampling parameters */
  ppr2->k_min_tau0 = 0.1;
  ppr2->k_max_tau0_over_l_max = 2;
  ppr2->k_step_sub = 0.1;
  ppr2->k_step_super = 0.025;
  ppr2->k_logstep_super = 1.8;
  ppr2->k_step_transition = 0.2;
  ppr2->k_per_decade_for_pk = 10;
  ppr2->k_per_decade_for_bao = 70;
  ppr2->k_bao_center = 3;
  ppr2->k_bao_width = 4;

  /* Transfer function k-sampling (used only if ptr2->k_sampling == class_transfer2_k3_sampling) */
  ppr2->q_linstep_song = 0.45;

  /* Transfer function tau-sampling (used only if ptr2->tau_sampling == custom_tau_sampling) */
  ppr2->tau_linstep_song = 0.6;
  
  /* Scalars */
  ppr2->k_min_custom = 1e-4;
  ppr2->k_max_custom = 0.1;  
  ppr2->k_size_custom = 10;

  /* k-triangular */
  ppr2->k3_size_min = 5;
  ppr2->k3_size = 100;
  

  // ====================================================================
  // =                           Bessel functions                       =
  // ====================================================================

  ppr2->bessel_j_cut_song = 1e-12;
  ppr2->bessel_J_cut_song = 1e-6;
  ppr2->bessel_x_step_song = 0.2;



  // ======================================================================
  // =                    Interpolation and integration                   =
  // ======================================================================
  
  ppr2->sources_time_interpolation = linear_interpolation;
  ppr2->sources_k3_interpolation = linear_interpolation;

  
  // ===============================================================================
  // =                                  Technical stuff                            =
  // ===============================================================================

  ppr2->old_run = _FALSE_;
  ppr2->store_sources_to_disk = _FALSE_;
  ppr2->load_sources_from_disk = _FALSE_;
  ppr2->store_transfers_to_disk = _FALSE_;
  ppr2->load_transfers_from_disk = _FALSE_;

  return _SUCCESS_;

}


int input2_free (struct precision2 * ppr2)
{
  
  free(ppr2->index_m_max);
    
  return _SUCCESS_;
  
}




int song_version(
      char * version
      )
{
  
  sprintf(version,"%s",_SONG_VERSION_);
  return _SUCCESS_;
}


