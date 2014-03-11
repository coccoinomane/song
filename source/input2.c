/** @file input2.c Documented input module for SONG.
 *
 * Guido W. Pettinari, 13.03.2013
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
            struct bessels * pbs,
            struct bessels2 * pbs2,
            struct transfers *ptr,
            struct transfers2 *ptr2,
            struct primordial *ppm,
            struct spectra *psp,
            struct bispectra *pbi,
            struct fisher *pfi,
            struct nonlinear *pnl,
            struct lensing *ple,
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

  class_call (input2_init(&fc,
                ppr,
                ppr2,
                pba,
                pth,
                ppt,
                ppt2,
                pbs,
                pbs2,
                ptr,
                ptr2,
                ppm,
                psp,
                pbi,
                pfi,
                pnl,
                ple,
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
         struct bessels * pbs,
         struct bessels2 * pbs2,
         struct transfers *ptr,
         struct transfers2 *ptr2,
         struct primordial *ppm,
         struct spectra *psp,
         struct bispectra *pbi,
         struct fisher *pfi,
         struct nonlinear * pnl,
         struct lensing *ple,
         struct output *pop,
         ErrorMsg errmsg
         )
{
          
  /** Summary: */

  /** - define local variables */

  int flag1,flag2,flag3;
  int int1;
  double param1,param2,param3;
  int entries_read;
  int * int_pointer;
  int * pointer_to_int;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];

  int i;


  /** - set all parameters (input and precision) to default values */

  class_call (input2_default_params(
               pba,
               pth,
               ppt,
               ppt2,
               pbs,
               pbs2,
               ptr,
               ptr2,
               ppm,
               psp,
               pbi,
               pfi,
               pnl,
               ple,
               pop),
    errmsg,
    errmsg);

  class_call (input2_default_precision (ppr2),
    errmsg,
    errmsg);

  
  // ========================================================
  // =                      Parse output                    =
  // ========================================================
 
  class_call (parser_read_string (pfc,"output",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"tBisp") != NULL) || (strstr(string1,"tBispectrum") != NULL) || (strstr(string1,"tB") != NULL)) {
      ppt2->has_cmb_temperature = _TRUE_;  
    }
    
    if (((strstr(string1,"pBisp") != NULL) || (strstr(string1,"pBispectrum") != NULL) || (strstr(string1,"PBISP") != NULL))
       ||(strstr(string1,"eBisp") != NULL) || (strstr(string1,"eBispectrum") != NULL) || (strstr(string1,"EBISP") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_; /* The intensity C_l's are needed for the delta_tilde transformation */
      ppt2->has_cmb_polarization_e = _TRUE_;  
    }

    if ((strstr(string1,"bBisp") != NULL) || (strstr(string1,"bBispectrum") != NULL) || (strstr(string1,"bBisp") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt2->has_cmb_polarization_b = _TRUE_;
    }

    if ((strstr(string1,"rBisp") != NULL) || (strstr(string1,"rBispectrum") != NULL) || (strstr(string1,"rB") != NULL)) {
      /* Second order Rayleigh is not supported (and not needed because we only compute <T^(2)R^(1)R^(1)>) */
    }
    
    /* Compute only first-order transfer functions */
    if (strstr(string1,"early_transfers1") != NULL) {
      ppt2->has_early_transfers1_only = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;    
    }    

    /* Compute only second-order early transfer functions */
    if (strstr(string1,"early_transfers2") != NULL) {
      ppt2->has_early_transfers2_only = _TRUE_;
      ppt2->has_perturbations2 = _TRUE_;      
      ppt2->has_cls = _FALSE_;
    }
    else {
      /* Compute only second-order transfer functions today */
      if (strstr(string1,"transfers2") != NULL) {
        ptr2->has_transfers2_only = _TRUE_;
        ppt2->has_perturbations2 = _TRUE_;
        ppt2->has_cls = _TRUE_;
      }
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
      ppt2->has_bispectra = _TRUE_;
    }
    
  } // end of bispectrum_types parsing
  


  // ====================================================================================
  // =                           Perturbations, time sampling                           =
  // ====================================================================================


  class_read_double("tau_start_evolution_2nd_order", ppt2->tau_start_evolution);

  class_read_double("recombination_max_to_end_ratio", ppt2->recombination_max_to_end_ratio);


  /* Do we need a custom time-sampling? */
  class_call(parser_read_string(pfc,"custom_time_sampling_for_2nd_order_sources",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);
   
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_custom_timesampling = _TRUE_;


  class_read_double("custom_tau_ini_2nd_order_sources", ppt2->custom_tau_ini);

  class_test (ppt2->custom_tau_ini<=0, errmsg, "please choose 'tau_ini' greater than zero.");
  
  class_read_double("custom_tau_end_2nd_order_sources", ppt2->custom_tau_end);  
  
  class_read_int("custom_tau_size_2nd_order_sources", ppt2->custom_tau_size);

  class_call(parser_read_string(pfc,"custom_tau_mode_2nd_order_sources",&string1,&flag1,errmsg),
       errmsg,
       errmsg); 

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"lin") != NULL) || (strstr(string1,"LIN") != NULL)))
      ppt2->custom_tau_mode = lin_tau_sampling;

    else if (((strstr(string1,"log") != NULL) || (strstr(string1,"LOG") != NULL)))
      ppt2->custom_tau_mode = log_tau_sampling;

    else if (((strstr(string1,"class") != NULL) || (strstr(string1,"CLASS") != NULL)))
      ppt2->custom_tau_mode = class_tau_sampling;
    
    else
      class_stop(errmsg,         
        "You wrote: tau_mode_2nd_order_sources=%s. Could not identify any of the supported time samplings ('lin', 'log', 'class') in such input",
        string1);

  }
   
  /* Precision parameters for the time sampling */
  class_read_double("perturb_sampling_stepsize_for_2nd_order", ppr2->perturb_sampling_stepsize_2nd_order);
  class_read_double("start_small_k_at_tau_c_over_tau_h_2nd_order", ppr2->start_small_k_at_tau_c_over_tau_h_2nd_order);
  class_read_double("start_large_k_at_tau_h_over_tau_k_2nd_order", ppr2->start_large_k_at_tau_h_over_tau_k_2nd_order);
  


  // ===========================================================================
  // =                        Perturbations, k-sampling                        =
  // ===========================================================================


  // ******        Sources k1-k2 sampling       *******


  /* Parse the method to use */
  class_call(parser_read_string(pfc,"sources2_k_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strcmp(string1,"lin") == 0) || (strcmp(string1,"linear") == 0))
      ppt2->k_sampling = lin_k_sampling;

    else if ((strcmp(string1,"log") == 0) || (strcmp(string1,"logarithmic")) == 0)
      ppt2->k_sampling = log_k_sampling;

    else if (strcmp(string1,"class") ==0)
      ppt2->k_sampling = class_sources_k_sampling;

    else if (strcmp(string1,"smart") == 0)
      ppt2->k_sampling = smart_sources_k_sampling;
      
    else
      class_stop(errmsg,
        "Could not recognize the value given for the 'sources2_k_sampling' option. Choose between 'lin', 'log' or 'smart'.", "");

  }

  /* Parameters for the smart k1-k2 sampling */
  class_read_double("k_scalar_min_tau0_2nd_order",ppr2->k_scalar_min_tau0_2nd_order);
  class_read_double("k_scalar_max_tau0_over_l_max_2nd_order",ppr2->k_scalar_max_tau0_over_l_max_2nd_order);
  class_read_double("k_scalar_step_sub_2nd_order",ppr2->k_scalar_step_sub_2nd_order);
  class_read_double("k_scalar_linstep_super_2nd_order",ppr2->k_scalar_linstep_super_2nd_order);
  class_read_double("k_scalar_logstep_super_2nd_order",ppr2->k_scalar_logstep_super_2nd_order);
  class_read_double("k_scalar_step_transition_2nd_order",ppr2->k_scalar_step_transition_2nd_order);

  class_test (ppr2->k_scalar_logstep_super_2nd_order <= 1,
    ppr2->error_message,
    "a logarithmic step must be larger than 1");
    

  /* Parameters for the custom lin/log sampling */
  class_read_double("k_min_scalars", ppr2->k_min_scalars);
  class_read_double("k_max_scalars", ppr2->k_max_scalars);
  class_read_int("k_size_scalars", ppr2->k_size_scalars);



  // ******        Sources k3 sampling       *******

  /* Parse the method for the k3 sampling */
  class_call(parser_read_string(pfc,"sources2_k3_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strcmp(string1,"lin") == 0) || (strcmp(string1,"linear") == 0))
      ppt2->k3_sampling = lin_k3_sampling;

    else if ((strcmp(string1,"log") == 0) || (strcmp(string1,"logarithmic")) == 0)
      ppt2->k3_sampling = log_k3_sampling;

    else if ((strcmp(string1,"smart") == 0) || (strcmp(string1,"class") ==0))
      ppt2->k3_sampling = smart_k3_sampling;

    else if (strcmp(string1,"theta_12") == 0)
      ppt2->k3_sampling = theta12_k3_sampling;

    else if (strcmp(string1,"theta_13") == 0)
      ppt2->k3_sampling = theta13_k3_sampling;
      
    else
      class_stop(errmsg,
        "Could not recognize the value given for the 'sources2_k3_sampling' option. Choose between 'lin', 'log', 'smart', 'theta_12', 'theta_13'.", "");

  }
  
  /* Minimum number of grid points for any (k1,k2) pair, used when 'k3_sampling' is set to smart */
  class_read_int("k3_size_min", ppr2->k3_size_min);
  
  /* Fixed number of grid points for any (k1,k2) pair, used when 'k3_sampling' is set to either lin or log */
  class_read_int("k3_size", ppr2->k3_size);






  // ====================================================================================
  // =                        Perturbations, differential system                        =
  // ====================================================================================

  // *** Set ppt2->has_quadratic_sources
  class_call(parser_read_string(pfc,"quadratic_sources",&(string1),&(flag1),errmsg),errmsg,errmsg);
            
  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_sources = _FALSE_;


  // *** Set ppt2->has_quadratic_liouville
  class_call(parser_read_string(pfc,"quadratic_liouville",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_liouville = _FALSE_;    


  // *** Set ppt2->has_quadratic_collision
  class_call(parser_read_string(pfc,"quadratic_collision",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quadratic_collision = _FALSE_;    

  /* Set quadliouville and quadcollision to _FALSE_ if the user asked for no quadratic sources */
  if (ppt2->has_quadratic_sources == _FALSE_) {
    ppt2->has_quadratic_liouville = _FALSE_;
    ppt2->has_quadratic_collision = _FALSE_;
  }

  // *** Set ppt2->has_time_delay_in_liouville
  class_call(parser_read_string(pfc,"include_time_delay_in_liouville",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_time_delay_in_liouville = _FALSE_;

  // *** Set ppt2->has_redshift_in_liouville
  class_call(parser_read_string(pfc,"include_redshift_in_liouville",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_redshift_in_liouville = _FALSE_;

  // *** Set ppt2->has_lensing_in_liouville
  class_call(parser_read_string(pfc,"include_lensing_in_liouville",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_lensing_in_liouville = _FALSE_;


  // *** Set ppt->polarization and ppt2->has_polarization2
  class_call(parser_read_string(pfc,"polarization_second_order",&(string1),&(flag1),errmsg),errmsg,errmsg);
     
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    ppt2->has_polarization2 = _TRUE_;
  }

  /* This fires up if the used did not specify a value for 'polarization_second_order', which by
  default is true. */
  if (ppt2->has_polarization2 == _TRUE_)
    ppt->has_polarization2 = _TRUE_;
  

  // *** Set ppt2->has_perfect_baryons
  class_call(parser_read_string(pfc,"perfect_baryons",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt2->has_perfect_baryons = _FALSE_;
  }


  // *** Set ppt2->has_perfect_cdm
  class_call(parser_read_string(pfc,"perfect_cdm",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt2->has_perfect_cdm = _FALSE_;
  }

  // ====================================================================================
  // =                      Perturbations, perturbed recombination                      =
  // ====================================================================================


  // *** Set ppt->has_perturbed_recombination
  class_call(parser_read_string(pfc,"perturbed_recombination",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    pth->has_perturbed_recombination = _TRUE_;
    ppt->has_perturbed_recombination = _TRUE_;
    ppt2->has_perturbed_recombination = _TRUE_;
    pth->compute_xe_derivatives = _TRUE_; /* Debug purposes */
  }

  /* From which value of x should we start integrating delta_Xe? */
  class_read_double("perturbed_recombination_turnx", pth->perturbed_recombination_turnx);

  /* Should we use the approximation in eq. 3.23 of Senatore, Tassev & Zaldarriaga 2009? */
  class_call(parser_read_string(pfc,"perturbed_recombination_use_approx",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    ppt2->perturbed_recombination_use_approx = _TRUE_;
    
    /* To use the analytical approximation, we only need the derivatives of the background ionization fraction X_e */
    if (pth->has_perturbed_recombination == _TRUE_)
      pth->compute_xe_derivatives = _TRUE_;
      
    pth->has_perturbed_recombination = _FALSE_;
    ppt->has_perturbed_recombination = _FALSE_;
  }
  

  // ===============================================================================
  // =                           Perturbations, LOS sources                        =
  // ===============================================================================


  // *** Set ppt2->has_pure_scattering_in_los
  class_call(parser_read_string(pfc,"include_pure_scattering_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_pure_scattering_in_los = _FALSE_;

  // *** Set ppt2->has_photon_monopole_in_los
  class_call(parser_read_string(pfc,"include_photon_monopole_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_photon_monopole_in_los = _FALSE_;
  
  // *** Set ppt2->has_quad_scattering_in_los
  class_call(parser_read_string(pfc,"include_quad_scattering_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
    ppt2->has_quad_scattering_in_los = _FALSE_;

  /* Metric sources exist only for temperature */
  if (ppt2->has_cmb_temperature == _TRUE_) {

    // *** Set ppt2->has_metric_in_los
    class_call(parser_read_string(pfc,"include_metric_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
      ppt2->has_metric_in_los = _FALSE_;

    // *** Set ppt2->has_quad_metric_in_los
    class_call(parser_read_string(pfc,"include_quad_metric_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))
      ppt2->has_quad_metric_in_los = _FALSE_;
  }

  // *** Set ppt2->has_time_delay_in_los
  class_call(parser_read_string(pfc,"include_time_delay_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_time_delay_in_los = _TRUE_;

  // *** Set ppt2->has_redshift_in_los
  class_call(parser_read_string(pfc,"include_redshift_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_redshift_in_los = _TRUE_;

  // *** Set ppt2->has_lensing_in_los
  class_call(parser_read_string(pfc,"include_lensing_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_lensing_in_los = _TRUE_;

  // *** Set ppt2->use_delta_tilde_in_los
  class_call(parser_read_string(pfc,"use_delta_tilde_in_los",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->use_delta_tilde_in_los = _TRUE_;

  /* Metric sources exist only for temperature */
  if (ppt2->has_cmb_temperature == _TRUE_) {

    // *** Set ppt2->has_sachs_wolfe_in_los
    class_call(parser_read_string(pfc,"include_sachs_wolfe_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppt2->has_sachs_wolfe_in_los = _TRUE_;

    // *** Set ppt2->has_integrated_sachs_wolfe_in_los
    class_call(parser_read_string(pfc,"include_integrated_sachs_wolfe_in_los_2nd_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppt2->has_integrated_sachs_wolfe_in_los = _TRUE_;


    /* Avoid counting twice the same metric effect */
    if ((ppt2->has_sachs_wolfe_in_los == _TRUE_) || (ppt2->has_integrated_sachs_wolfe_in_los == _TRUE_))
      ppt2->has_metric_in_los = _FALSE_;
    
    /* Set the integration_by_parts flag on if we need to obtain non-standard LOS source terms */
    if ((ppt2->has_sachs_wolfe_in_los == _TRUE_) || (ppt2->has_integrated_sachs_wolfe_in_los == _TRUE_))
      ppt2->has_integration_by_parts_of_los = _TRUE_;    
  }

  /* If effects that are not peaked at recombination are included, we need to extend the integration range up to today */
  if ((ppt2->has_metric_in_los == _FALSE_) && (ppt2->has_integrated_sachs_wolfe_in_los == _FALSE_)
      && (ppt2->has_quad_metric_in_los == _FALSE_) && (ppt2->has_time_delay_in_los == _FALSE_)
      && (ppt2->has_redshift_in_los == _FALSE_) && (ppt2->has_lensing_in_los == _FALSE_))
    ppt2->has_recombination_only = _TRUE_;
  

  /* Should we define the SW effect as in Huang and Vernizzi 2013, i.e. including the -psi*psi contribution? */
  if (ppt2->has_sachs_wolfe_in_los == _TRUE_) {
    
    // *** Set ppt2->has_sachs_wolfe_in_los
    class_call(parser_read_string(pfc,"use_zhiqi_sw",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppt2->use_zhiqi_sw = _TRUE_;
    
    class_test_permissive (ppt2->has_quad_metric_in_los == _TRUE_,
      errmsg,
      "flag 'use_zhiqi_sw' and 'quad_metric' are incompatible: you are double counting the psi*psi terms!");
    
  }

  /* Should we include only the early ISW effect? */
  if (ppt2->has_integrated_sachs_wolfe_in_los == _TRUE_) {
    
    // *** Set ppt2->has_sachs_wolfe_in_los
    class_call(parser_read_string(pfc,"only_early_isw",&(string1),&(flag1),errmsg),errmsg,errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

      ppt2->only_early_isw = _TRUE_;

      /* Unless other late time effects are included, stop integrating just after recombination */
      if ((ppt2->has_metric_in_los == _FALSE_) && (ppt2->has_quad_metric_in_los == _FALSE_)
          && (ppt2->has_time_delay_in_los == _FALSE_) && (ppt2->has_redshift_in_los == _FALSE_) && (ppt2->has_lensing_in_los == _FALSE_))
        ppt2->has_recombination_only = _TRUE_;
      else
        class_test_permissive (1==1,
          errmsg,
          "flag 'only_early_isw' not implemented yet for 'quad_metric=yes'. Will just include the full ISW.");

    }

  }



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
    ppt2->has_sachs_wolfe_in_los = _FALSE_;
    ppt2->has_integrated_sachs_wolfe_in_los = _FALSE_;
		ppt2->has_recombination_only = _TRUE_;	
		ppt2->has_integration_by_parts_of_los = _FALSE_;
		
	}









  // ====================================================================================
  // =                        Perturbations, initial conditions                         =
  // ====================================================================================


  /* Only adiabatic or vanishing initial conditions are supported so far */

  class_call(parser_read_string(pfc,"ic_2nd_order",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_) {

    /* if no initial conditions are specified, the default is has_ad=_TRUE_; 
       but if they are specified we should reset has_ad to _FALSE_ before reading */
    ppt2->has_ad=_FALSE_;

    if ((strcmp(string1,"ad") == 0) || (strcmp(string1,"AD") == 0))
      ppt2->has_ad=_TRUE_; 

    if ((strstr(string1,"ad_first_order") != NULL) || (strstr(string1,"AD_1ST") != NULL) || (strstr(string1,"ad1") != NULL))
      ppt2->has_ad_first_order=_TRUE_; 

    if ((strstr(string1,"zero") != NULL) || (strstr(string1,"AD_ZERO") != NULL))
      ppt2->has_zero_ic=_TRUE_; 

    if (strstr(string1,"unphysical") != NULL)
      ppt2->has_unphysical_ic=_TRUE_; 
      
    class_test(ppt2->has_ad==_FALSE_ && ppt2->has_ad_first_order==_FALSE_ && ppt2->has_zero_ic==_FALSE_ && ppt2->has_unphysical_ic==_FALSE_,
      errmsg,         
      "You wrote: ic_2nd_order=%s. Could not identify any of the supported initial conditions ('ad', 'ad_first_order', 'zero', 'unphysical') in such input",string1);

  }

  /* What is the value of the primordial non-Gaussianity of the local type? */
  class_read_double("primordial_local_fnl_phi", ppt2->primordial_local_fnl_phi);






  // *** Set ppt2->match_final_time_los
  class_call(parser_read_string(pfc,"match_final_time_los",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->match_final_time_los = _TRUE_;



  // ====================================================================================
  // =                           Perturbations, approximations                          =
  // ====================================================================================
  

  /* Tight coupling.  Note that, as opposite to the 1st order case, we read the approximation settings
  directly into the ppt2 structure rather than in the precision one.  We do so in order to keep all
  2nd-order related quantities in the same structure */
  class_read_int("tight_coupling_approximation_2nd_order", ppt2->tight_coupling_approximation);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_h_2nd_order", ppt2->tight_coupling_trigger_tau_c_over_tau_h);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_k_2nd_order", ppt2->tight_coupling_trigger_tau_c_over_tau_k);


  /* Radiation streaming */
  class_read_int("radiation_streaming_approximation_2nd_order", ppt2->radiation_streaming_approximation);
  class_read_double("radiation_streaming_trigger_tau_over_tau_k_2nd_order", ppt2->radiation_streaming_trigger_tau_over_tau_k);

  class_read_int("ur_fluid_approximation_2nd_order", ppt2->ur_fluid_approximation);
  class_read_double("ur_fluid_trigger_tau_over_tau_k_2nd_order", ppt2->ur_fluid_trigger_tau_over_tau_k);

  /* No radiation approximation */
  class_read_int("no_radiation_approximation_2nd_order", ppt2->no_radiation_approximation);
  class_read_double("no_radiation_approximation_rho_m_over_rho_r_2nd_order", ppt2->no_radiation_approximation_rho_m_over_rho_r);

  class_test(ppt2->ur_fluid_trigger_tau_over_tau_k == ppt2->radiation_streaming_trigger_tau_over_tau_k,
    errmsg,
    "please choose different values for precision parameters ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");




  // ====================================================================================
  // =                        Perturbations, gauges & equations                         =
  // ====================================================================================

  class_call( parser_read_string(pfc,
                "phi_prime_equation",
                &(string1),
                &(flag1),
                errmsg),
    errmsg,
    errmsg);
     
  if (flag1 == _TRUE_) {
  
    if ( (strcmp(string1,"poisson") == 0) || (strcmp(string1,"POISSON") == 0)||
         (strcmp(string1,"Poisson") == 0) || (strcmp(string1,"P") == 0) || (strcmp(string1,"p") == 0) ) {
      ppt2->phi_prime_eq = poisson;
    }
    else if ( (strcmp(string1,"longitudinal") == 0) || (strcmp(string1,"LONGITUDINAL") == 0)||
         (strcmp(string1,"Longitudinal") == 0) || (strcmp(string1,"L") == 0) || (strcmp(string1,"l") == 0) ) {
      ppt2->phi_prime_eq = longitudinal;
    }
    else {
      class_stop (errmsg, "unrecognized parameter '%s' for the field 'phi_prime_equation'", string1);
    }
  }
  
  /*  Set the ppt2->lm_extra parameter according to the gauge */
  if (ppt->gauge == 0)
    ppt2->lm_extra = 1;

  else if (ppt->gauge == 1)
    ppt2->lm_extra = 3;       
  
  
  

  // ================================================================================
  // =                            Perturbations, debug                              =
  // ================================================================================

  class_read_int("perturbations2_verbose",ppt2->perturbations2_verbose);

  // *** Set ppt2->has_debug_files  (should be before the filenames are read)
  class_call(parser_read_string(pfc,
       "dump_debug_files",&(string1),&(flag1),errmsg),errmsg,errmsg);
            
  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt2->has_debug_files = _TRUE_;


  if (ppt2->has_debug_files == _TRUE_) {

    /* Where to store the evolved transfer functions (optional, used only if in debug mode) */
    class_call(parser_read_string(pfc,
         "transfers_filename",&(string1),&(flag1),errmsg),errmsg,errmsg);  
    
    if ((flag1 == _TRUE_) && (string1 != NULL) && (ppt2->has_debug_files==_TRUE_)) {
      strcpy(ppt2->transfers_filename, string1);
    }
    class_open(ppt2->transfers_file,ppt2->transfers_filename,"w",errmsg);

    /* Where to store the quadratic sources (optional, used only if in debug mode) */
    class_call(parser_read_string(pfc,"quadsources_filename",&(string1),&(flag1),errmsg),errmsg,errmsg);  
    
    if ((flag1 == _TRUE_) && (string1 != NULL) && (ppt2->has_debug_files==_TRUE_)) {
      strcpy(ppt2->quadsources_filename, string1);
    }
    class_open(ppt2->quadsources_file,ppt2->quadsources_filename,"w",errmsg);

    /* Where to store the quadratic part of the liouville operator (optional, used only if in debug mode) */
    class_call(parser_read_string(pfc,"quadliouville_filename",&(string1),&(flag1),errmsg),errmsg,errmsg);  
    
    if ((flag1 == _TRUE_) && (string1 != NULL) && (ppt2->has_debug_files==_TRUE_)) {
      strcpy(ppt2->quadliouville_filename, string1);
    }
    class_open(ppt2->quadliouville_file,ppt2->quadliouville_filename,"w",errmsg);


    /* Where to store the quadratic part of the collision term (optional, used only if in debug mode) */
    class_call(parser_read_string(pfc,"quadcollision_filename",&(string1),&(flag1),errmsg),errmsg,errmsg);  
    
    if ((flag1 == _TRUE_) && (string1 != NULL) && (ppt2->has_debug_files==_TRUE_)) {
      strcpy(ppt2->quadcollision_filename, string1);
    }
    class_open(ppt2->quadcollision_file,ppt2->quadcollision_filename,"w",errmsg);

    /* Set which wavemode to print to file (optional, used only if in debug mode) */
    class_read_int("index_k1_debug",
      ppt2->index_k1_debug);

    class_read_int("index_k2_debug",
      ppt2->index_k2_debug);

    class_read_int("index_k3_debug",
      ppt2->index_k3_debug);

    /* Set how many angular modes to print to file (optional, used only if in debug mode) */
    class_read_int("l_max_debug",
      ppt2->l_max_debug);
  
  } // end of if has_debug_files




  // =================================================================================
  // =                                Bessel functions                               =
  // =================================================================================

  class_read_int("bessels2_verbose",pbs2->bessels2_verbose);

  /* Make sure that we compute all the needed Bessel functions for the bispectrum */
  if (ppt2->has_bispectra) 
    pbs2->extend_l1_using_m = _TRUE_;

  /* Minimum x treshold for the spherical Bessels with j_l1(x).  These are the
    Bessels that are summed to obtain J_Llm */
  class_read_double("bessel_j_cut_2nd_order", ppr2->bessel_j_cut_2nd_order);

  /* Minimum treshold for the functions J_Llm(x).  These are obtained as a
    weighted sum of spherical Bessels j_l1(x) with |L-l| <= l1 <= L+l */
  class_read_double("bessel_J_cut_2nd_order", ppr2->bessel_J_cut_2nd_order);

  /* Linear step dx where we are going to sample the j_l1(x) and J_Llm(x) */
  class_read_double("bessel_x_step_2nd_order", ppr2->bessel_x_step_2nd_order);

  /* Tolerance for the integration of the 2nd-order system.  This parameter goes
    directly into the evolver as the parameter 'rtol' */
  class_read_double("tol_perturb_integration_2nd_order",ppr2->tol_perturb_integration_2nd_order);










  // =========================================================================================
  // =                                  Transfer functions                                   =
  // =========================================================================================

  class_read_int("transfer2_verbose", ptr2->transfer2_verbose);

  // ******        Tranfer function k-sampling       *******
  class_call(parser_read_string(pfc,"transfer2_k_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if (strstr(string1,"bessel") != NULL)
      ptr2->k_sampling = bessel_k_sampling;

    else if ((strstr(string1,"class") != NULL) || (strstr(string1,"smart") != NULL))
      ptr2->k_sampling = class_transfer2_k_sampling;
    
    else
      class_stop(errmsg,
        "Could not recognize the value given for the 'transfer2_k_sampling' option. Choose between 'bessel' or 'class'.", "");
  }

  /* If "class", choose the density of the k3-sampling for the transfer functions */
  class_read_double("k_step_trans_scalars_2nd_order", ppr2->k_step_trans_scalars_2nd_order);


  // ******        Tranfer function tau-sampling       *******
  class_call(parser_read_string(pfc,"transfer2_tau_sampling",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if (strstr(string1,"bessel") != NULL)
      ptr2->tau_sampling = bessel_tau_sampling;

    else if ((strstr(string1,"custom") != NULL) || (strstr(string1,"smart") != NULL))
      ptr2->tau_sampling = custom_transfer2_tau_sampling;
    
    else
      class_stop(errmsg,
        "Could not recognize the value given for the 'transfer2_tau_sampling' option. Choose between 'bessel' or 'smart'.", "");

  }

  /* If "custom", choose the density of the tau-sampling for the transfer functions */
  class_read_double("tau_step_trans_2nd_order", ppr2->tau_step_trans_2nd_order);



  

  // =============================================================================================
  // =                                       Bispectra                                           =
  // =============================================================================================





  // =============================================================================================
  // =                                   Fisher matrices                                         =
  // =============================================================================================



  

  // =============================================================================================
  // =                           Storage of intermediate results                                 =
  // =============================================================================================


  /* Which intermediate results should we store? */
  if ( (ppr->store_run == _TRUE_) || (ppr->load_run == _TRUE_) ) {
    

    // *** Store line of sight sources?
    class_call(parser_read_string(pfc,"store_sources",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);
     
    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppt2->store_sources_to_disk = _TRUE_;



    // *** Store transfer functions?
    class_call(parser_read_string(pfc,"store_transfers",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);
     
    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ptr2->store_transfers_to_disk = _TRUE_;


  } // end of if (store_run)


  /* The flag ppt2->store_sources_to_disk has two different meanings according to whether ppr->load_run is false
    or true.  If ppr->load_run == _TRUE_, then store_sources_to_disk == _TRUE_ indicates that the run that we are
    loading contains the sources.  If ppr->load_run == _FALSE_, then store_sources_to_disk == _TRUE_ means that we shall
    save the sources to disk.  When ppt2->store_sources_to_disk is _FALSE_, we assume that there is no storing
    nor loading of the sources.  As a consequence, when ppt2->store_sources_to_disk == _FALSE_ the usage of RAM
    is maximized, as all the sources need to be stored in ppt2->sources at the same time.  This is why we define
    the negation of ppt2->store_sources_to_disk as ppt2->keep_sources_in_memory.
    
    The same reasoning applies to ptr2->store_transfers_to_disk. */
  ppt2->load_sources_from_disk = _FALSE_;
  if ((ppr->load_run==_TRUE_) && (ppt2->store_sources_to_disk==_TRUE_))
    ppt2->load_sources_from_disk = _TRUE_;


  ptr2->load_transfers_from_disk = _FALSE_;
  if ((ppr->load_run==_TRUE_) && (ptr2->store_transfers_to_disk==_TRUE_))
    ptr2->load_transfers_from_disk = _TRUE_;







  // =============================================================================================
  // =                                  Interpolation techniques                                 =
  // =============================================================================================


  // ****   Interpolation of 2nd-order sources   ****
  
  class_call(parser_read_string(pfc,"sources_time_interpolation",&string1,&flag1,errmsg),
       errmsg,
       errmsg); 

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr2->sources_time_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr2->sources_time_interpolation = cubic_interpolation;
    
    else
      class_stop (errmsg,         
        "You wrote: sources_time_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",
        string1);

  }


  class_call(parser_read_string(pfc,"sources_k3_interpolation",&string1,&flag1,errmsg),
       errmsg,
       errmsg); 

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr2->sources_k3_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr2->sources_k3_interpolation = cubic_interpolation;
    
    else
      class_stop(errmsg,         
        "You wrote: sources_k3_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",
        string1);

  }





  // =========================================================================================
  // =                                   Angular scales                                      =
  // =========================================================================================


  // *** Read l_max for the Boltzmann hierarchies
  class_read_int("l_max_g_2nd_order", ppr2->l_max_g_2nd_order);
  class_read_int("l_max_pol_g_2nd_order", ppr2->l_max_pol_g_2nd_order);    
  class_read_int("l_max_ur_2nd_order", ppr2->l_max_ur_2nd_order);        
  class_read_int("l_max_g_ten_2nd_order", ppr2->l_max_g_ten_2nd_order);       
  class_read_int("l_max_pol_g_ten_2nd_order", ppr2->l_max_pol_g_ten_2nd_order);            
   
              

  // *** Read l_max for the quadratic sources

  /* If the user specified a negative value for one of them, set it to the corresponding
  l_max_2nd_order (see above). Also make sure that each of them is not larger than the
  corresponding l_max_2nd_order. The following lines MUST go below the definitions of
  l_max_g_2nd_order, etc. */

  class_read_int("l_max_g_quadsources", ppr2->l_max_g_quadsources);
  if ((ppr2->l_max_g_quadsources<0) || (ppr2->l_max_g_quadsources>ppr2->l_max_g_2nd_order))
    ppr2->l_max_g_quadsources = ppr2->l_max_g_2nd_order;

  class_read_int("l_max_pol_g_quadsources", ppr2->l_max_pol_g_quadsources);
  if ((ppr2->l_max_pol_g_quadsources<0) || (ppr2->l_max_pol_g_quadsources>ppr2->l_max_pol_g_2nd_order))
    ppr2->l_max_pol_g_quadsources = ppr2->l_max_pol_g_2nd_order;

  class_read_int("l_max_ur_quadsources", ppr2->l_max_ur_quadsources);
  if ((ppr2->l_max_ur_quadsources<0) || (ppr2->l_max_ur_quadsources>ppr2->l_max_ur_2nd_order))
    ppr2->l_max_ur_quadsources = ppr2->l_max_ur_2nd_order;

  class_read_int("l_max_g_ten_quadsources", ppr2->l_max_g_ten_quadsources);
  if ((ppr2->l_max_g_ten_quadsources<0) || (ppr2->l_max_g_ten_quadsources>ppr2->l_max_g_ten_2nd_order))
    ppr2->l_max_g_ten_quadsources = ppr2->l_max_g_ten_2nd_order;

  class_read_int("l_max_pol_g_ten_quadsources", ppr2->l_max_pol_g_ten_quadsources);
  if ((ppr2->l_max_pol_g_ten_quadsources<0) || (ppr2->l_max_pol_g_ten_quadsources>ppr2->l_max_pol_g_ten_2nd_order))
    ppr2->l_max_pol_g_ten_quadsources = ppr2->l_max_pol_g_ten_2nd_order;

  // *** Read l_max for the line of sight integration
  class_read_int("l_max_T_los", ppr2->l_max_los_t); /* obsolete, use l_max_los_t */
  class_read_int("l_max_E_los", ppr2->l_max_los_p); /* obsolete, use l_max_los_p */
  class_read_int("l_max_B_los", ppr2->l_max_los_p); /* obsolete, use l_max_los_p */
  class_read_int("l_max_los_t", ppr2->l_max_los_t);
  class_read_int("l_max_los_p", ppr2->l_max_los_p);


  // *** Same as above, but for the quadratic terms
  class_read_int("l_max_T_quadratic_los", ppr2->l_max_los_quadratic_t); /* obsolete, use l_max_los_quadratic_t */
  class_read_int("l_max_E_quadratic_los", ppr2->l_max_los_quadratic_p); /* obsolete, use l_max_los_quadratic_p */
  class_read_int("l_max_B_quadratic_los", ppr2->l_max_los_quadratic_p); /* obsolete, use l_max_los_quadratic_p */
  class_read_int("l_max_los_quadratic_t", ppr2->l_max_los_quadratic_t);
  class_read_int("l_max_los_quadratic_p", ppr2->l_max_los_quadratic_p);


  // *** Compute the maximum l_max for the line of sight integration  
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
                pfc,
                "modes_2nd_order",
                &(ppr2->m_size),
                &(int_pointer),
                &flag1,
                errmsg),
    errmsg,
    errmsg);
    
  if (flag1 == _TRUE_) {
    for (i=0; i < ppr2->m_size; ++i)
      ppr2->m[i] = int_pointer[i];
    free (int_pointer);
  }

  /* Maximum 'm' that will be computed */
  ppr2->m_max_2nd_order = ppr2->m[ppr2->m_size-1];
  
  /* Check that the m-list is strictly ascending */
  for (i=0; i < (ppr2->m_size-1); ++i)
    class_test (ppr2->m[i] >= ppr2->m[i+1],
      errmsg,
      "the m-list provided  in 'modes_2nd_order' is not strictly ascending");

  /* Check that the m's positive */
  class_test (ppr2->m[0] < 0,
    errmsg,
    "the m-list provided in 'modes_2nd_order' has negative numbers in it");

  /* Check that m_max is smaller than limit */  
  class_test (ppr2->m_max_2nd_order > (_MAX_NUM_AZIMUTHAL_-1),
    errmsg,
    "the maximum value of the azimuthal number 'm' cannot exceed %d, please choose 'modes_2nd_order' accordingly",
    _MAX_NUM_AZIMUTHAL_);

  /* For m>0 we cannot compute the reduced bispectrum when l1+l2+l3 is odd */
  // if ((pbi->has_intrinsic==_TRUE_) && (
  //   (ppr2->m_max_2nd_order>0) ||
  //   (ppt2->has_cmb_polarization_e) ||
  //   (ppt2->has_cmb_polarization_b))) {
  // // if ((pbi->has_intrinsic==_TRUE_) && (
  // //   (ppr2->m_max_2nd_order>0) ||
  // // if (pbi->has_intrinsic==_TRUE_) {
  //      
  //   printf ("\n");
  //   printf ("   *@^#?!?! FORCING THE COMPUTATION OF A GRID OF EVEN L'S\n");
  //   printf ("\n");      
  //   ppr->compute_only_even_ls = _TRUE_;
  // 
  //   // printf ("\n");
  //   // printf ("   *@^#?!?! FORCING THE COMPUTATION OF A GRID OF ODD L'S\n");
  //   // printf ("\n");
  //   // ppr->compute_only_odd_ls = _TRUE_;
  // 
  //   class_test_permissive (
  //     ((ppt2->has_cmb_polarization_b==_TRUE_) && (ppr->compute_only_even_ls==_TRUE_)) ||
  //     ((ppt2->has_cmb_polarization_e==_TRUE_) && (ppr->compute_only_odd_ls==_TRUE_)) ||
  //     ((ppt2->has_cmb_temperature==_TRUE_) && (ppr->compute_only_odd_ls==_TRUE_)),
  //     errmsg,
  //     "careful, your choice of parity is wrong!");
  // }

  /* Find out the index in ppr2->m corresponding to a given m. */
  class_alloc (ppr2->index_m, (ppr2->m_max_2nd_order+1)*sizeof(int), errmsg);

	for(int m=0; m<=ppr2->m_max_2nd_order; ++m) {

		ppr2->index_m[m] = -1;

		for (int index_m=0; index_m<ppr2->m_size; ++index_m)
			if (m==ppr2->m[index_m]) ppr2->index_m[m] = index_m;
	}

	/* Build a logical array to check whether a given m is in m_vec.
  First, intialise it to _FALSE_. */
	for(int m=0; m<_MAX_NUM_AZIMUTHAL_; ++m)
    ppr2->compute_m[m] = _FALSE_;

	for(int m=0; m<=ppr2->m_max_2nd_order; ++m)
		for (int index_m=0; index_m<ppr2->m_size; ++index_m)
			if (m==ppr2->m[index_m]) ppr2->compute_m[m] = _TRUE_;

  /* Set the scalar, vector, tensor flags */
  ppt2->has_scalars = ppr2->compute_m[0];
  ppt2->has_vectors = ppr2->compute_m[1];
  ppt2->has_tensors = ppr2->compute_m[2];

  // *** Determine what m's are allowed for a given l 
  
  /* The largest l's we will ever use in the code. The first term (pbs->l_max) is the maximum multipole
  we shall compute the bispectra in. The second term represent the additional l's where we need to
  the Bessel functions in order to compute the projection functions (ppr2->l_max_los) and the bispectrum
  integral (ppr2->m_max_2nd_order). 
  Must be after setting ppr2->l_max_los. */
  class_test (pbs->l_max==0., errmsg, "pbs->l_max=0. Did you forget to set ppt2->has_cls==_TRUE_ in input.c?");
  int l_max = pbs->l_max + MAX (ppr2->l_max_los, ppr2->m_max_2nd_order);
  class_alloc (ppr2->index_m_max, (l_max+1)*sizeof(int), ppr2->error_message);

  /* Determine for each l (starting from 0) the position of its maximum m in ppr2->m. The resulting
  array ppr2->index_m_max will be used throughout the code to cycle through the m's allowed for a given l. */
  int l;
  for (l=0; l<=l_max; ++l) {
   
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


  // ==============================================================================
  // =                          Quadratic sources rescaling                       =
  // ==============================================================================

  /* Do we want to rescale all multipoles with a factor 1/sin(theta_1)^m? The rescaling does not
  affect the m>0 transfer functions. It is needed to compute the intrinsic bispectrum, so by default
  it is active. See the header file perturbations2.h for details on the rescaling. */
  ppt2->rescale_quadsources = _TRUE_;

  /* Uncomment if you want the m=0 sources to be computed without the rescaling, when they are
  the only requested sources. */
  // if (ppr2->m_max_2nd_order == 0)
  //   ppt2->rescale_quadsources = _TRUE_;
  
  /* Uncomment if you want the output functions to output non-rescaled functions */
  if ((ppt2->has_early_transfers2_only == _TRUE_) || (ptr2->has_transfers2_only == _TRUE_))
    ppt2->rescale_quadsources = _FALSE_;


  
  // =============================================================================================
  // =                               Stuff that needs to be at the end                           =
  // =============================================================================================

  // -----------------------------------------------
  // -             Update pbs->x_max               -
  // -----------------------------------------------

  /* Determine maximum k-value needed in the line of sight integration at second-order. This is
  needed to determine the maximum sampling of the Bessel functions. */
  double k_max;
  double tau0_inf = 8000.;
  double tau0_sup = 15000.;
  
  if ((ppt2->k_sampling == class_sources_k_sampling) || (ppt2->k_sampling == smart_sources_k_sampling))
    /* Automatic k-sampling */
    k_max = ppr2->k_scalar_max_tau0_over_l_max_2nd_order * ppt->l_scalar_max / tau0_inf;
  else
    /* Manual k-sampling */
    k_max = ppr2->k_max_scalars;

  /* Maximum argument for the Bessel functions in the line of sight integration.  For the time being, we do
  set x_max  here, even if we shouldn't as tau0 cannot be accessed by this module.  This is why we define
  a superior limit for tau0. */
  double x_max = MAX (pbs->x_max, k_max * MAX(tau0_sup, pbi->r_max));
  

  /* We could have done this way, but then it is not obvious how to include the r_max in the
  class_sources_k_sampling case */
  // /* Determine maximum argument for the Bessel functions in the line of sight integration. */
  // double x_max;
  // 
  // if (ppt2->k_sampling == class_sources_k_sampling)
  //   /* At second order, we include a 2 factor in x_max, because ppr2->k_scalar_max_tau0_over_l_max_2nd_order is used
  //     to build the grid in k1 and k2, while the maximum k-value needed will be k = k1+k2 (when cosk1k2 = 1)  */
  //   x_max = 2 * ppr2->k_scalar_max_tau0_over_l_max_2nd_order * ppt->l_scalar_max;
  // else
  //   /* When we specify the k-range manually, also x_max should be set by hand to be k_max * tau0.  For the time
  //     being, we do set x_max here, even if we shouldn't as tau0 cannot be accessed by this module.  This is why
  //     we define a superior limit for tau0.  Note that we also include a 2 factor in x_max, because k_max_scalars
  //     refers to k1 and k2, while the maximum k-value needed will be k = k1+k2 (when cosk1k2 = 1). */
  //   x_max = 2 * ppr2->k_max_scalars * MAX(tau0_sup, pbi->r_max);


  /* Set the actual limit for x_max */
  printf("# Temporary message: Setting pbs->x_max from %g to %g\n",
    ((int)(pbs->x_max * 1.1 / pbs->x_step)+1)*pbs->x_step,
    ((int)(x_max * 1.1 / pbs->x_step)+1)*pbs->x_step);
  
  pbs->x_max = ((int)(x_max * 1.1 / pbs->x_step)+1)*pbs->x_step;


  /* Find the maximum x for which we shall compute J_Llm(x) and j_l1(x),
  needed to compute the LOS integral and the spectra at second order. */

  /* In CLASS, the maximum values are determined in the following way:
   - l_max is determined by the params.ini parameters l_scalar_max, l_tensor_max, delta_l_max
   - k_max_cl is determined by k_scalar_max_tau0_over_l_max, given that tau0 is fixed for a given cosmology
   - x_max is determined by l_max and k_scalar_max_tau0_over_l_max, since the latter is just x_max*l_max.
   
   In SONG, for the time being, we choose not to change neither l_max or k_max.  However this can only
   work if you change your get_k_list in perturbations2.c to be equal to the k_list in perturbations.c     */

  /* Note that pbs->x_max should be larger/equal than pbs->xx_max, as we shall interpolate the Bessels
  with bessels_at_x in the bispectra module, with x that can be as large as conformal_age * k_max */

  pbs2->xx_max = pbs->x_max;
  pbs2->xx_step = ppr2->bessel_x_step_2nd_order;
  
  /* pbs2->L_max is used to determine the maximum value of L for which we should compute the
  projection functions in the second-order Bessel module. To explain the inclusion of
  ppr2->m_max_2nd_order, refer to the first long comment in bessel2_get_l1_list. */
  pbs2->L_max = MAX (ppr2->l_max_los, ppr2->m_max_2nd_order);

  return _SUCCESS_;

}


/** 
 * All default parameter values (for input parameters)
 *
 * @param pba Input : pointer to background structure 
 * @param pth Input : pointer to thermodynamics structure 
 * @param ppt Input : pointer to perturbation structure
 * @param pbs Input : pointer to bessels structure
 * @param ptr Input : pointer to transfer structure 
 * @param ppm Input : pointer to primordial structure 
 * @param psp Input : pointer to spectra structure
 * @param pop Input : pointer to output structure
 * @return the error status
 */

int input2_default_params (
       struct background *pba,
       struct thermo *pth,
       struct perturbs *ppt,
       struct perturbs2 *ppt2,       
       struct bessels * pbs,
       struct bessels2 * pbs2,
       struct transfers *ptr,
       struct transfers2 *ptr2,     
       struct primordial *ppm,
       struct spectra *psp,
       struct bispectra *pbi,
       struct fisher *pfi,
       struct nonlinear * pnl,
       struct lensing *ple,
       struct output *pop
       ) {


  // ==============================================================
  // =                   perturbed recombination                  =
  // ==============================================================

  pth->has_perturbed_recombination = _FALSE_;
  pth->perturbed_recombination_turnx = 36;
  ppt->has_perturbed_recombination = _FALSE_;
  ppt2->has_perturbed_recombination = _FALSE_;
  ppt2->perturbed_recombination_use_approx = _FALSE_;
  

  // ===========================================================
  // =                     perturb1 structure                  =
  // ===========================================================


  // ===========================================================
  // =                     perturb2 structure                  =
  // ===========================================================
  
  // *** Flags
  ppt2->perturbations2_verbose = 0;
  ppt2->has_perturbations2 = _FALSE_;
  ppt2->has_polarization2  = _TRUE_;
  ppt2->has_quadratic_sources = _TRUE_;
  ppt2->has_quadratic_liouville = _TRUE_;
  ppt2->has_quadratic_collision = _TRUE_;
  ppt2->has_perfect_baryons = _TRUE_;
  ppt2->has_perfect_cdm = _TRUE_;
  ptr2->has_transfers2_only = _FALSE_;
  ppt2->store_sources_to_disk = _FALSE_;
  ppt2->rescale_quadsources = _TRUE_;

  ppt2->rescale_quadsources = _FALSE_;

  ppt2->has_time_delay_in_liouville = _TRUE_;
  ppt2->has_redshift_in_liouville = _TRUE_;
  ppt2->has_lensing_in_liouville = _TRUE_;
  
  ppt2->has_pure_scattering_in_los = _TRUE_;
  ppt2->has_photon_monopole_in_los = _TRUE_;
  ppt2->has_quad_scattering_in_los = _TRUE_;
  ppt2->has_metric_in_los = _TRUE_;
  ppt2->has_quad_metric_in_los = _TRUE_;

  ppt2->has_time_delay_in_los = _FALSE_;
  ppt2->has_redshift_in_los = _FALSE_;
  ppt2->has_lensing_in_los = _FALSE_;

  ppt2->use_delta_tilde_in_los = _FALSE_;
  
  ppt2->has_integration_by_parts_of_los = _FALSE_;
  ppt2->has_sachs_wolfe_in_los = _FALSE_;
  ppt2->use_zhiqi_sw = _FALSE_;
  ppt2->has_integrated_sachs_wolfe_in_los = _FALSE_;
  ppt2->only_early_isw = _FALSE_;

	ppt2->use_test_source = _FALSE_;

  ppt2->has_recombination_only = _FALSE_;

  /* Possible outputs at 2nd order */
  ppt2->has_cmb_temperature = _FALSE_;
  ppt2->has_cmb_polarization_e = _FALSE_; 
  ppt2->has_cmb_polarization_b = _FALSE_; 
  ppt2->has_pk_matter = _FALSE_;
  
  ppt2->has_cls = _FALSE_;
  ppt2->has_bispectra = _FALSE_;
  
  ppt2->has_early_transfers2_only = _FALSE_;
  ppt2->has_early_transfers1_only = _FALSE_;
  
  // *** Initial conditions
  ppt2->has_ad = _TRUE_;      
  ppt2->has_zero_ic = _FALSE_;
  ppt2->has_ad_first_order = _FALSE_;
  ppt2->has_unphysical_ic = _FALSE_;
  ppt2->primordial_local_fnl_phi = 0;




  // *** Approximations at second order

  ppt2->tight_coupling_approximation = tca2_none;
  ppt2->tight_coupling_trigger_tau_c_over_tau_h = 0.015;
  ppt2->tight_coupling_trigger_tau_c_over_tau_k = 0.01;
  
  ppt2->radiation_streaming_approximation = rsa2_none;
  ppt2->radiation_streaming_trigger_tau_over_tau_k = 45;
  
  ppt2->ur_fluid_approximation = ufa2_none;
  ppt2->ur_fluid_trigger_tau_over_tau_k = 15;
  
  ppt2->no_radiation_approximation = nra2_none;
  ppt2->no_radiation_approximation_rho_m_over_rho_r = 5;


  // *** Choose equations
  ppt2->phi_prime_eq = poisson;


  // *** Time sampling
  ppt2->tau_start_evolution = 5;

  ppt2->recombination_max_to_end_ratio = 1000;

  ppt2->has_custom_timesampling = _FALSE_;
  ppt2->custom_tau_ini  = 1;
  ppt2->custom_tau_end  = 0;      // 0 means: integrate all the way to today
  ppt2->custom_tau_size = 2000;
  ppt2->custom_tau_mode = lin_tau_sampling;

  ppt2->match_final_time_los = _FALSE_;
  
  // *** K sampling
  ppt2->k_sampling = lin_k_sampling;
  ppt2->k3_sampling = lin_k3_sampling;


  // *** Technical parameters

  ppt2->has_debug_files = _FALSE_;
  strcpy(ppt2->transfers_filename,"output/transfers2.txt");
  strcpy(ppt2->quadsources_filename,"output/quadsources.txt");
  strcpy(ppt2->quadliouville_filename,"output/quadliouville.txt");
  strcpy(ppt2->quadcollision_filename,"output/quadcollision.txt");

  ppt2->index_k1_debug = 0;
  ppt2->index_k2_debug = 0;
  ppt2->index_k3_debug = 0;
  ppt2->l_max_debug = 5;
  

  
  // =============================================================
  // =                     bessel2 structure                     =
  // =============================================================
  pbs2->bessels2_verbose = 0;
  pbs2->extend_l1_using_m = _FALSE_;

  
  // ============================================================
  // =                     transfer2 structure                  =
  // ============================================================
  ptr2->transfer2_verbose = 0;
  ptr2->k_sampling = class_transfer2_k_sampling;
  ptr2->tau_sampling = custom_transfer2_tau_sampling;
  ptr2->store_transfers_to_disk = _FALSE_;



  // ============================================================
  // =                     bispectra structure                  =
  // ============================================================


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


  // ******      Tolerance       *******

  ppr2->tol_perturb_integration_2nd_order=1.e-5;




  // ******      Multipole limits       *******

  ppr2->m_max_2nd_order=0;

  ppr2->l_max_g_2nd_order=6;
  ppr2->l_max_g_quadsources=6;
  ppr2->l_max_pol_g_2nd_order=6;
  ppr2->l_max_pol_g_quadsources=6;
  ppr2->l_max_ur_2nd_order=6; 
  ppr2->l_max_ur_quadsources=6;   
  ppr2->l_max_g_ten_2nd_order=6;
  ppr2->l_max_g_ten_quadsources=6;
  ppr2->l_max_pol_g_ten_2nd_order=6;
  ppr2->l_max_pol_g_ten_quadsources=6;

  ppr2->l_max_los_t=4;
  ppr2->l_max_los_quadratic_t=4;
  ppr2->l_max_los_p=4;
  ppr2->l_max_los_quadratic_p=4;

  
  /* By default, compute only the scalar (m=0) modes */
  ppr2->m_size = 1;
  ppr2->m[0] = 0;


  // ******      Time samplings       *******

  ppr2->perturb_sampling_stepsize_2nd_order = 0.08;
  ppr2->start_small_k_at_tau_c_over_tau_h_2nd_order = 0.0015;  /* decrease to start earlier in time */
  ppr2->start_large_k_at_tau_h_over_tau_k_2nd_order = 0.07;    /* decrease to start earlier in time */



  // ******        k-sampling      *******
  
  /* CLASS smart sampling */
  ppr2->k_scalar_min_tau0_2nd_order = 1.;
  ppr2->k_scalar_max_tau0_over_l_max_2nd_order = 2.;
  ppr2->k_scalar_step_sub_2nd_order = 0.05;
  ppr2->k_scalar_linstep_super_2nd_order = 0.0025;
  ppr2->k_scalar_logstep_super_2nd_order = 1.2;
  ppr2->k_scalar_step_transition_2nd_order = 0.2;

  /* Transfer function k-sampling (used only if ptr2->k_sampling == class_transfer2_k_sampling) */
  ppr2->k_step_trans_scalars_2nd_order = 0.04;

  /* Transfer function tau-sampling (used only if ptr2->tau_sampling == custom_transfer2_tau_sampling) */
  ppr2->tau_step_trans_2nd_order = 0.25;

  
  /* Scalars */
  ppr2->k_min_scalars = 1e-4;
  ppr2->k_max_scalars = 0.1;  
  ppr2->k_size_scalars = 10;


  /* k-triangular */
  ppr2->k3_size_min = 5;
  ppr2->k3_size = 100;
  



  // *******      Bessel functions       ********
  ppr2->bessel_j_cut_2nd_order = 1e-12;
  ppr2->bessel_J_cut_2nd_order = 1e-6;
  ppr2->bessel_x_step_2nd_order = 0.3;



  // *******     Interpolation and integration      *********
  ppr2->sources_time_interpolation = linear_interpolation;
  ppr2->sources_k3_interpolation = linear_interpolation;


  return _SUCCESS_;

}


int precision2_free (struct precision2 * ppr2)
{
  
  free(ppr2->m);
  free(ppr2->compute_m);
  free(ppr2->index_m);
  free(ppr2->index_m_max);
    
  return _SUCCESS_;
  
}




int song_version(
      char * version
      ) {
  
  sprintf(version,"%s",_SONG_VERSION_);
  return _SUCCESS_;
}


