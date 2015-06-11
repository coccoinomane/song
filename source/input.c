/** @file input.c Documented input module.
 *
 * Julien Lesgourgues, 27.08.2010
 */

#include "input.h"

/**
 * Use this routine to extract initial parameters from files 'xxx.ini'
 * and/or 'xxx.pre'. They can be the arguments of the main() routine.
 *
 * If class is embedded into another code, you will probably prefer to
 * call directly input_init() in order to pass input parameters
 * through a 'file_content' structure.
 */

int input_init_from_arguments(
			      int argc,
			      char **argv,
			      struct precision * ppr,
			      struct background *pba,
			      struct thermo *pth,
			      struct perturbs *ppt,
			      struct bessels * pbs,
			      struct transfers *ptr,
			      struct primordial *ppm,
			      struct spectra *psp,
			      struct bispectra *pbi,
			      struct fisher *pfi,
			      struct nonlinear *pnl,
            struct lensing *ple,
			      struct output *pop,
			      ErrorMsg errmsg
			      ) {

  /** Summary: */

  /** - define local variables */

  // *** MY MODIFICATIONS ****
  class_alloc (ppr->input_file_content, sizeof(struct file_content), errmsg);
  struct file_content * pfc = ppr->input_file_content;
  // *** ORIGINAL CLASS
  // struct file_content fc;
  // *** END OF MY MODIFICATIONS ****
  struct file_content fc_input;
  struct file_content fc_precision;

  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  int i;
  char extension[5];

  /** - Initialize the two file_content structures (for input
      parameters and precision parameters) to some null content. If no
      arguments are passed, they will remain null and inform
      init_params() that all parameters take default values. */

  pfc->size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';


  // *** MY MODIFICATIONS ***

  /* Check that the first argument exists as a file or as a directory */
  struct stat st;
  class_test ((argc>1) && (stat (argv[1], &st) != 0),
		errmsg,
		"the specified argument does not exist ('%s')", argv[1]);

  /* Check whether the first argument is a directory, and if this is the case use the parameter files
  that are stored inside  */
  ppr->load_run = (S_ISDIR (st.st_mode) != 0);

  if ((argc == 2) && (ppr->load_run == _TRUE_)) {

    /* From now on, we shall assume that the specified directory contains a previously stored run */
    strcpy(ppr->run_dir, argv[1]);

    /* Store the paths of the parameter files into the strings 'input_file' and 'precision_file'. Later on
    we shall also stored them into the structure fields ppr->ini_filename and ppr->pre_filename. */
    sprintf (input_file, "%s/run_params.ini", ppr->run_dir);
    sprintf (precision_file, "%s/run_params.pre", ppr->run_dir);

    /* It is mandatory that the run directory contains a params.ini and a params.pre file.  We
    now check that they exist */
    class_test (stat (input_file, &st) != 0,
      errmsg,
      "the run directory does not contain the parameter file '%s'", input_file);

    class_test (stat (precision_file, &st) != 0,
      errmsg,
      "the run directory does not contain the parameter file '%s'", input_file);

    printf("# We shall load the run contained in the folder '%s'.\n", ppr->run_dir);

  }

  /** If some arguments are passed, identify eventually some 'xxx.ini'
      and 'xxx.pre' files, and store their name. */

  if (ppr->load_run == _FALSE_) {

    if (argc > 1) {
      for (i=1; i<argc; i++) {
        strncpy(extension,(argv[i]+strlen(argv[i])-4),4);
        extension[4]='\0';
        if (strcmp(extension,".ini") == 0) {
  	class_test(input_file[0] != '\0',
  		   errmsg,
  		   "You have passed more than one input file with extension '.ini', choose one.");
  	strcpy(input_file,argv[i]);
        }
        if (strcmp(extension,".pre") == 0) {
  	class_test(precision_file[0] != '\0',
  		   errmsg,
  		   "You have passed more than one precision with extension '.pre', choose one.");
  	strcpy(precision_file,argv[i]);
        }
      }
    }

  } // end of if(ppr->load_run == _FALSE_)



  /** - if there is an 'xxx.ini' file, read it and store its content. */

  if (input_file[0] != '\0') {

    class_call(parser_read_file(input_file,&fc_input,errmsg),
         errmsg,
         errmsg);

    /* Save the parameter filename into the ppr structure */
    strcpy (ppr->ini_filename, input_file);

  }

  /** - if there is an 'xxx.pre' file, read it and store its content. */

  if (precision_file[0] != '\0') {

    class_call(parser_read_file(precision_file,&fc_precision,errmsg),
         errmsg,
         errmsg);

    /* Save the parameter filename into the ppr structure */
    strcpy (ppr->pre_filename, precision_file);

  }


  // *** END OF MY MODIFICATIONS ***







  /** - if one or two files were read, merge their contents in a
      single 'file_content' structure. */

  if ((input_file[0]!='\0') || (precision_file[0]!='\0'))

    class_call(parser_cat(&fc_input,&fc_precision,pfc,errmsg),
	       errmsg,
	       errmsg);

  class_call(parser_free(&fc_input),errmsg,errmsg);
  class_call(parser_free(&fc_precision),errmsg,errmsg);

  /** - now, initialize all parameters given the input 'file_content'
      structure.  If its size is null, all parameters take their
      default values. */

  class_call(input_init(
      pfc,
			ppr,
			pba,
			pth,
			ppt,
			pbs,
			ptr,
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


  // *** MY MODIFICATIONS ***

  /* Commented the following line so that fc stays in ppr->input_file_content */
  // class_call(parser_free(&fc),errmsg,errmsg);
  // *** END OF MY MODIFICATIONS ***

  return _SUCCESS_;
}




/**
 * Initialize each parameters, first to its default values, and then
 * from what can be interpreted from the values passed in the input
 * 'file_content' structure. If its size is null, all parameters keep
 * their default values.
 */

int input_init(
	       struct file_content * pfc,
	       struct precision * ppr,
	       struct background *pba,
	       struct thermo *pth,
	       struct perturbs *ppt,
	       struct bessels * pbs,
	       struct transfers *ptr,
	       struct primordial *ppm,
	       struct spectra *psp,
	       struct bispectra *pbi,
         struct fisher *pfi,
	       struct nonlinear * pnl,
         struct lensing *ple,
	       struct output *pop,
	       ErrorMsg errmsg
	       ) {

  /** Summary: */

  /** - define local variables */

  int flag1,flag2,flag3;
  double param1,param2,param3;
  int N_ncdm=0,n,entries_read;
  int int1,fileentries;
  double fnu_factor;
  double * pointer1;
  int * pointer_to_int;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];
  char string[_ARGUMENT_LENGTH_MAX_];
  char buffer[1024];

  double Omega_tot;

  int i;

  FILE * param_output;
  FILE * param_unused;
  char param_output_name[_LINE_LENGTH_MAX_];
  char param_unused_name[_LINE_LENGTH_MAX_];

  double sigma_B; /**< Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

  double rho_ncdm;

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);

  /** - set all parameters (input and precision) to default values */

  class_call(input_default_params(pba,
				  pth,
				  ppt,
				  pbs,
				  ptr,
				  ppm,
				  psp,
				  pbi,
          pfi,
				  pnl,
				  ple,
				  pop),
	     errmsg,
	     errmsg);

  class_call(input_default_precision(ppr),
	     errmsg,
	     errmsg);

  /** - if entries passed in file_content structure, carefully read
      and interpret each of them, and tune accordingly the relevant
      input parameters */




  /** (a) background parameters */

  /* h (dimensionless) and [H0/c] in Mpc^{-1} = h / 2999.7 */
  class_call(parser_read_double(pfc,"H0",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"h",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
	     errmsg,
	     "In input file, you cannot enter both h and H0, choose one");
  if (flag1 == _TRUE_) {
    pba->H0 = param1 * 1.e3 / _c_;
    pba->h = param1 / 100.;
  }
  if (flag2 == _TRUE_) {
    pba->H0 = param2 *  1.e5 / _c_;
    pba->h = param2;
  }

  /* Omega_0_g (photons) and T_cmb */
  class_call(parser_read_double(pfc,"T_cmb",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"Omega_g",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"omega_g",&param3,&flag3,errmsg),
	     errmsg,
	     errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
	     errmsg,
	     "In input file, you can only enter one of T_cmb, Omega_g or omega_g, choose one");

  if (class_none_of_three(flag1,flag2,flag3)) {
    pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  }
  else {

    if (flag1 == _TRUE_) {
      /* Omega0_g = rho_g / rho_c0, each of them expressed in Kg/m/s^2 */
      /* rho_g = (4 sigma_B / c) T^4 */
      /* rho_c0 = 3 c^2 H0^2 / (8 pi G) */
      pba->Omega0_g = (4.*sigma_B/_c_*pow(param1,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
      pba->T_cmb=param1;
    }

    if (flag2 == _TRUE_) {
      pba->Omega0_g = param2;
      pba->T_cmb=pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*sigma_B/_c_),0.25);
    }

    if (flag3 == _TRUE_) {
      pba->Omega0_g = param3/pba->h/pba->h;
      pba->T_cmb = pow(pba->Omega0_g * (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_) / (4.*sigma_B/_c_),0.25);
    }
  }

  Omega_tot = pba->Omega0_g;

  /* Omega_0_b (baryons) */
  class_call(parser_read_double(pfc,"Omega_b",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"omega_b",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_b or omega_b, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_b = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_b = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_b;

  /* Omega_0_ur (ultra-relativistic species / massless neutrino) */
  class_call(parser_read_double(pfc,"N_eff",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"Omega_ur",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"omega_ur",&param3,&flag3,errmsg),
	     errmsg,
	     errmsg);
  class_test(class_at_least_two_of_three(flag1,flag2,flag3),
	     errmsg,
	     "In input file, you can only enter one of N_eff, Omega_ur or omega_ur, choose one");

  if (class_none_of_three(flag1,flag2,flag3)) {
    pba->Omega0_ur = 3.04*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  }
  else {

    if (flag1 == _TRUE_) {
      pba->Omega0_ur = param1*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
    }
    if (flag2 == _TRUE_) {
      pba->Omega0_ur = param2;
    }
    if (flag3 == _TRUE_) {
      pba->Omega0_ur = param3/pba->h/pba->h;
    }
  }

  Omega_tot += pba->Omega0_ur;

  /* Omega_0_cdm (CDM) */
  class_call(parser_read_double(pfc,"Omega_cdm",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"omega_cdm",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
	     errmsg,
	     "In input file, you can only enter one of Omega_cdm or omega_cdm, choose one");
  if (flag1 == _TRUE_)
    pba->Omega0_cdm = param1;
  if (flag2 == _TRUE_)
    pba->Omega0_cdm = param2/pba->h/pba->h;

  Omega_tot += pba->Omega0_cdm;

  /* non-cold relics (ncdm) */
  class_read_int("N_ncdm",N_ncdm);
  if ((flag1 == _TRUE_) && (N_ncdm > 0)){
    pba->N_ncdm = N_ncdm;
    /* Precision parameters for ncdm has to be read now since they are used here:*/
    class_read_double("tol_M_ncdm",ppr->tol_M_ncdm);
    class_read_double("tol_ncdm",ppr->tol_ncdm);
    class_read_double("tol_ncdm_bg",ppr->tol_ncdm_bg);

    /* Read temperatures: */
    class_read_list_of_doubles_or_default("T_ncdm",pba->T_ncdm,pba->T_ncdm_default,N_ncdm);

    /* Read chemical potentials: */
    class_read_list_of_doubles_or_default("ksi_ncdm",pba->ksi_ncdm,pba->ksi_ncdm_default,N_ncdm);

    /* Read degeneracy of each ncdm species: */
    class_read_list_of_doubles_or_default("deg_ncdm",pba->deg_ncdm,pba->deg_ncdm_default,N_ncdm);

    /* Read mass of each ncdm species: */
    class_read_list_of_doubles_or_default("m_ncdm",pba->m_ncdm_in_eV,0.0,N_ncdm);

    /* Read Omega of each ncdm species: */
    class_read_list_of_doubles_or_default("Omega_ncdm",pba->Omega0_ncdm,0.0,N_ncdm);

    /* Read omega of each ncdm species: (Use pba->M_ncdm temporarily)*/
    class_read_list_of_doubles_or_default("omega_ncdm",pba->M_ncdm,0.0,N_ncdm);

    /* Check for duplicate Omega/omega entries, missing mass definition and
       update pba->Omega0_ncdm:*/
    for(n=0; n<N_ncdm; n++){
      /* pba->M_ncdm holds value of omega */
      if (pba->M_ncdm[n]!=0.0){
	class_test(pba->Omega0_ncdm[n]!=0,errmsg,
		   "Nonzero values for both Omega and omega for ncdm species %d are specified!",n);
	pba->Omega0_ncdm[n] = pba->M_ncdm[n]/pba->h/pba->h;
      }
      if ((pba->Omega0_ncdm[n]==0.0) && (pba->m_ncdm_in_eV[n]==0.0)) {
	/* this is the right place for passing the default value of
	   the mass (all parameters must have a default value; most of
	   them are defined in input_default_params{}, but the ncdm mass
	   is a bit special and there is no better place for setting its
	   default value). We put an aribitrary value m << 10^-3 eV,
	   i.e. the ultra-relativistic limit.*/
	pba->m_ncdm_in_eV[n]=1.e-5;
      }
    }

    /* Check if filenames for interpolation tables are given: */
    class_read_list_of_integers_or_default("use_ncdm_psd_files",pba->got_files,_FALSE_,N_ncdm);

    if (flag1==_TRUE_){
      for(n=0,fileentries=0; n<N_ncdm; n++){
	if (pba->got_files[n] == _TRUE_) fileentries++;
      }

      if (fileentries > 0) {

	/* Okay, read filenames.. */
	class_call(parser_read_list_of_strings(pfc,"ncdm_psd_filenames",
					       &entries_read,&(pba->ncdm_psd_files),&flag2,errmsg),
		   errmsg,
		   errmsg);
	class_test(flag2 == _FALSE_,errmsg,
		   "Input use_ncdm_files is found, but no filenames found!");
	class_test(entries_read != fileentries,errmsg,
		   "Numer of filenames found, %d, does not match number of _TRUE_ values in use_ncdm_files, %d",
		   entries_read,fileentries);
      }
    }
/* Read (optional) p.s.d.-parameters:*/
    parser_read_list_of_doubles(pfc,
				"ncdm_psd_parameters",
				&entries_read,
				&(pba->ncdm_psd_parameters),
				&flag2,
				errmsg);

    class_call(background_ncdm_init(ppr,pba),
	       pba->error_message,
	       errmsg);

    /* We must calculate M from omega or vice versa if one of them is missing.
       If both are present, we must update the degeneracy parameter to
       reflect the implicit normalisation of the distribution function.*/
    for (n=0; n < N_ncdm; n++){
      if (pba->m_ncdm_in_eV[n] != 0.0){
	/* Case of only mass or mass and Omega/omega: */
	pba->M_ncdm[n] = pba->m_ncdm_in_eV[n]/_k_B_*_eV_/pba->T_ncdm[n]/pba->T_cmb;
	class_call(background_ncdm_momenta(pba->q_ncdm_bg[n],
					   pba->w_ncdm_bg[n],
					   pba->q_size_ncdm_bg[n],
					   pba->M_ncdm[n],
					   pba->factor_ncdm[n],
					   0.,
					   NULL,
					   &rho_ncdm,
					   NULL,
					   NULL,
					   NULL),
		   pba->error_message,
		   errmsg);
	if (pba->Omega0_ncdm[n] == 0.0){
	  pba->Omega0_ncdm[n] = rho_ncdm/pba->H0/pba->H0;
	}
	else{
	  fnu_factor = (pba->H0*pba->H0*pba->Omega0_ncdm[n]/rho_ncdm);
	  pba->factor_ncdm[n] *= fnu_factor;
	  /* dlnf0dlnq is already computed, but it is
	     independent of any normalisation of f0.
	     We don't need the factor anymore, but we
	     store it nevertheless:*/
	  pba->deg_ncdm[n] *=fnu_factor;
	}
      }
      else{
	/* Case of only Omega/omega: */
	class_call(background_ncdm_M_from_Omega(ppr,pba,n),
		   pba->error_message,
		   errmsg);
	printf("M_ncdm:%g\n",pba->M_ncdm[n]);
	pba->m_ncdm_in_eV[n] = _k_B_/_eV_*pba->T_ncdm[n]*pba->M_ncdm[n]*pba->T_cmb;
      }
      pba->Omega0_ncdm_tot += pba->Omega0_ncdm[n];
      //printf("Adding %g to total Omega..\n",pba->Omega0_ncdm[n]);
    }
  }
  Omega_tot += pba->Omega0_ncdm_tot;

  /* Omega_0_k (curvature) */
  class_read_double("Omega_k",pba->Omega0_k);

  /* Omega_0_lambda (cosmological constant), Omega0_fld (dark energy fluid) */
  class_call(parser_read_double(pfc,"Omega_Lambda",&param1,&flag1,errmsg),
	     errmsg,
	     errmsg);
  class_call(parser_read_double(pfc,"Omega_fld",&param2,&flag2,errmsg),
	     errmsg,
	     errmsg);
  class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
	     errmsg,
	     "In input file, you can enter only two out of Omega_Lambda, Omega_de, Omega_k, the third one is inferred");

  if ((flag1 == _FALSE_) && (flag2 == _FALSE_)) {
    pba->Omega0_lambda = 1.+pba->Omega0_k-pba->Omega0_g-pba->Omega0_ur-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_ncdm_tot;
  }
  else {
    if (flag1 == _TRUE_) {
      pba->Omega0_lambda= param1;
      pba->Omega0_fld = 1. + pba->Omega0_k - param1 - Omega_tot;
    }
    if (flag2 == _TRUE_) {
      pba->Omega0_lambda= 1. + pba->Omega0_k - param2 - Omega_tot;
      pba->Omega0_fld = param2;
    }
  }

  if (pba->Omega0_fld != 0.) {
    class_read_double("w0_fld",pba->w0_fld);
    class_read_double("wa_fld",pba->wa_fld);
    class_read_double("cs2_fld",pba->cs2_fld);

    class_test(pba->w0_fld<=-1.,
	       errmsg,
	       "Your choice w_fld=%g is not valid, it will lead to instabilities or division by zero\n",
	       pba->w0_fld);

    class_test(pba->w0_fld+pba->w0_fld>=1./3.,
	       errmsg,
	       "Your choice for w0_fld+wa_fld=%g is suspicious, ther would not be radiation domination at early times\n",
	       pba->w0_fld+pba->wa_fld);


  }


  class_test(pba->Omega0_k != 0.,
	     errmsg,
	     "Open/close case not written yet");

  /* scale factor today (arbitrary) */
  class_read_double("a_today",pba->a_today);

  /** (b) assign values to thermodynamics cosmological parameters */

  /* primordial helium fraction */
  class_call(parser_read_string(pfc,"YHe",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"BBN") != NULL) || (strstr(string1,"bbn") != NULL)) {
      pth->YHe = _BBN_;
    }
    else {
      class_read_double("YHe",pth->YHe);
    }

  }

  /* recombination parameters */
  class_call(parser_read_string(pfc,"recombination",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"HYREC") != NULL) || (strstr(string1,"hyrec") != NULL) || (strstr(string1,"HyRec") != NULL)) {
      pth->recombination = hyrec;
    }

  }

  // *** MY MODIFICATIONS ***
  /* Read the frequency for the Cl's of the Rayleigh scattering, and compute the ratio between
  the Rayleigh and Thomson interaction rates at such frequency. The formula is given in arXiv:1307.8148. */
  class_read_double("rayleigh_frequency",pth->rayleigh_frequency);
  // *** END OF MY MODIFICATIONS ***

  /* reionization parametrization */
  class_call(parser_read_string(pfc,"reio_parametrization",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {
    flag2=_FALSE_;
    if (strcmp(string1,"reio_none") == 0) {
      pth->reio_parametrization=reio_none;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_camb") == 0) {
      pth->reio_parametrization=reio_camb;
      flag2=_TRUE_;
    }
    if (strcmp(string1,"reio_bins_tanh") == 0) {
      pth->reio_parametrization=reio_bins_tanh;
      flag2=_TRUE_;
    }

    class_test(flag2==_FALSE_,
	       errmsg,
	       "could not identify reionization_parametrization value, check that it is one of 'reio_none', 'reio_camb', 'reio_bins_tanh', ...");
  }

  // *** MY_MODIFICATIONS ***
  if (pth->reio_parametrization != reio_none)
    ppr->has_reionization = _TRUE_;
  // *** END OF MY_MODIFICATIONS ***

  /* reionization parameters if reio_parametrization=reio_camb */
  if (pth->reio_parametrization == reio_camb) {
    class_call(parser_read_double(pfc,"z_reio",&param1,&flag1,errmsg),
	       errmsg,
	       errmsg);
    class_call(parser_read_double(pfc,"tau_reio",&param2,&flag2,errmsg),
	       errmsg,
	       errmsg);
    class_test(((flag1 == _TRUE_) && (flag2 == _TRUE_)),
	       errmsg,
	       "In input file, you can only enter one of z_reio or tau_reio, choose one");
    if (flag1 == _TRUE_) {
      pth->z_reio=param1;
      pth->reio_z_or_tau=reio_z;
    }
    if (flag2 == _TRUE_) {
      pth->tau_reio=param2;
      pth->reio_z_or_tau=reio_tau;
    }

    class_read_double("reionization_exponent",pth->reionization_exponent);
    class_read_double("reionization_width",pth->reionization_width);
    class_read_double("helium_fullreio_redshift",pth->helium_fullreio_redshift);
    class_read_double("helium_fullreio_width",pth->helium_fullreio_width);

  }

  /* reionization parameters if reio_parametrization=reio_bins_tanh */
  if (pth->reio_parametrization == reio_bins_tanh) {
    class_read_int("binned_reio_num",pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_z",pth->binned_reio_z,pth->binned_reio_num);
    class_read_list_of_doubles("binned_reio_xe",pth->binned_reio_xe,pth->binned_reio_num);
    class_read_double("binned_reio_step_sharpness",pth->binned_reio_step_sharpness);
  }


  /** (c) define which perturbations and sources should be computed, and down to which scale */

  ppt->has_perturbations = _FALSE_;
  ppt->has_cls = _FALSE_;

  class_call(parser_read_string(pfc,"output",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {


    // ==================================================================================
    // =                                  Class outputs                                 =
    // ==================================================================================

    if ((strstr(string1,"tCl") != NULL) || (strstr(string1,"TCl") != NULL) || (strstr(string1,"TCL") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"pCl") != NULL) || (strstr(string1,"PCl") != NULL) || (strstr(string1,"PCL") != NULL)) {
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    // *** MY MODIFICATIONS ***
    if ((strstr(string1,"rCl") != NULL) || (strstr(string1,"RCl") != NULL) || (strstr(string1,"RCL") != NULL)) {
      pth->has_rayleigh_scattering = _TRUE_;
      ppt->has_cl_cmb_rayleigh = _TRUE_;
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;

    }

    if ((strstr(string1,"zCl") != NULL) || (strstr(string1,"ZCl") != NULL) || (strstr(string1,"ZCL") != NULL)) {
      ppt->has_cl_cmb_zeta = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }
    // *** END OF MY MODIFICATIONS ***

    if ((strstr(string1,"lCl") != NULL) || (strstr(string1,"LCl") != NULL) || (strstr(string1,"LCL") != NULL)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"dCl") != NULL) || (strstr(string1,"DCl") != NULL) || (strstr(string1,"DCL") != NULL)) {
      ppt->has_cl_density=_TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if ((strstr(string1,"mPk") != NULL) || (strstr(string1,"MPk") != NULL) || (strstr(string1,"MPK") != NULL)) {
      ppt->has_pk_matter=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }

    if ((strstr(string1,"mTk") != NULL) || (strstr(string1,"MTk") != NULL) || (strstr(string1,"MTK") != NULL)) {
      ppt->has_matter_transfers=_TRUE_;
      ppt->has_perturbations = _TRUE_;
    }


    // *** MY MODIFICATIONS ***

    // ========================================================================================
    // =                                 Primordial bispectra                                 =
    // ========================================================================================

    if ((strstr(string1,"mBisp") != NULL) || (strstr(string1,"mBispectrum") != NULL) || (strstr(string1,"mB") != NULL)) {
      ppt->has_perturbations = _TRUE_;
      ppt->has_bispectra = _TRUE_;
      pbi->has_bispectra = _TRUE_;
    }

    if ((strstr(string1,"tBisp") != NULL) || (strstr(string1,"tBispectrum") != NULL) || (strstr(string1,"tB") != NULL)) {
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_bi_cmb_temperature = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
      ppt->has_bispectra = _TRUE_;
      pbi->has_bispectra = _TRUE_;
    }

    if ((strstr(string1,"pBisp") != NULL) || (strstr(string1,"pBispectrum") != NULL) || (strstr(string1,"PBISP") != NULL)
      ||(strstr(string1,"eBisp") != NULL) || (strstr(string1,"eBispectrum") != NULL) || (strstr(string1,"EBISP") != NULL)) {
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_bi_cmb_polarization = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
      ppt->has_bispectra = _TRUE_;
      pbi->has_bispectra = _TRUE_;
    }

    if ((strstr(string1,"rBisp") != NULL) || (strstr(string1,"rBispectrum") != NULL) || (strstr(string1,"rB") != NULL)) {
      pth->has_rayleigh_scattering = _TRUE_;
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_cl_cmb_rayleigh = _TRUE_;
      ppt->has_bi_cmb_rayleigh = _TRUE_;
      ppt->has_perturbations = _TRUE_;
      ppt->has_cls = _TRUE_;
      ppt->has_bispectra = _TRUE_;
      pbi->has_bispectra = _TRUE_;
    }

    if (strstr(string1,"fisher") != NULL) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_; /* The covariance matrix needs lensed C_l's */
      pfi->has_fisher = _TRUE_;
    }

    // ================================================================================
    // =                                 Song outputs                                 =
    // ================================================================================

    if (strstr(string1,"early_transfers1") != NULL) {
      ppt->has_perturbations = _TRUE_;
      ppt->has_perturbations2 = _TRUE_;
    }

    if (strstr(string1,"early_transfers2") != NULL) {
      ppt->has_perturbations = _TRUE_;
      ppt->has_perturbations2 = _TRUE_;
      ppt->has_cls = _TRUE_; /* otherwise CLASS sets pbs->l_max=0 */
    }

    if (strstr(string1,"transfers2") != NULL) {
      ppt->has_perturbations = _TRUE_;
      ppt->has_perturbations2 = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    if (strstr(string1,"tCl2") != NULL) {
      ppt->has_perturbations = _TRUE_;
      ppt->has_perturbations2 = _TRUE_;
      ppt->has_cl_cmb_temperature = _TRUE_;
      ppt->has_cl_cmb_polarization = _TRUE_;
      ppt->has_cls = _TRUE_;
    }

    // *** END OF MY MODIFICATIONS ***

  } // end of output parsing


  // *** MY MODIFICATIONS ***

  /* Parse the types of bispectra (intrinsic, primordial...) to compute */

  class_call (parser_read_string(pfc,"bispectrum_types",&string1,&flag1,errmsg),
    errmsg,
    errmsg);

  if ((pbi->has_bispectra == _TRUE_) && (flag1 == _TRUE_)) {

    /* Compute primordial bispectrum with a local shape function */
    if (strstr(string1,"local") != NULL) {
      pbi->has_local_model = _TRUE_;
    }

    /* Compute primordial bispectrum with an equilateral shape function */
    if (strstr(string1,"equilateral") != NULL) {
      pbi->has_equilateral_model = _TRUE_;
    }

    /* Compute primordial bispectrum with an orthogonal shape function */
    if (strstr(string1,"orthogonal") != NULL) {
      pbi->has_orthogonal_model = _TRUE_;
    }

    /* Compute primordial bispectrum with the two Galileon shapes in arXiv:1303.2125 */
    if (strstr(string1,"galileon") != NULL) {
      pbi->has_galileon_model = _TRUE_;
    }

    /* Compute the approximation of the intrinsic bispectrum in the squeezed limit, the lensed version */
    if ((strstr(string1,"intrinsic_squeezed") != NULL)
    || (strstr(string1,"i_squeezed") != NULL) || (strstr(string1,"i-squeezed") != NULL)) {
      pbi->has_intrinsic_squeezed = _TRUE_;
      ppt->has_cl_cmb_zeta = _TRUE_;
      psp->compute_cl_derivative = _TRUE_;
    }

    /* Compute the approximation of the intrinsic bispectrum in the squeezed limit, the unlensed version */
    if ((strstr(string1,"unlensed_intrinsic_squeezed") != NULL)
    || (strstr(string1,"i_u_squeezed") != NULL) || (strstr(string1,"i-u-squeezed") != NULL)) {
      pbi->has_intrinsic_squeezed_unlensed = _TRUE_;
      ppt->has_cl_cmb_zeta = _TRUE_;
    }

    /* Compute the approximation of the local bispectrum in the squeezed limit */
    if ((strstr(string1,"local_squeezed") != NULL)
    || (strstr(string1,"l_squeezed") != NULL) || (strstr(string1,"l-squeezed") != NULL)) {
      pbi->has_local_squeezed = _TRUE_;
      ppt->has_cl_cmb_zeta = _TRUE_;
      psp->compute_cl_derivative = _TRUE_;
    }

    /* Compute a test oscillating bispectrum */
    if ((strstr(string1,"cosine_shape") != NULL) || (strstr(string1,"cosine") != NULL)) {
      pbi->has_cosine_shape = _TRUE_;
    }

    /* Compute the CMB-lensing bispectrum */
    if ((strstr(string1,"cmb_lensing") != NULL) || (strstr(string1,"cmb-lensing") != NULL)
    || (strstr(string1,"CMB_lensing") != NULL) || (strstr(string1,"CMB-lensing") != NULL)
    || (strstr(string1,"isw_lensing") != NULL) || (strstr(string1,"isw-lensing") != NULL)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      pbi->has_cmb_lensing = _TRUE_;
    }

    /* Compute the approximation of the CMB-lensing bispectrum in the squeezed limit */
    if ((strstr(string1,"c_squeezed") != NULL) || (strstr(string1,"c-squeezed") != NULL)) {
      ppt->has_cl_cmb_lensing_potential = _TRUE_;
      pbi->has_cmb_lensing_squeezed = _TRUE_;
    }

    /* Compute the quadratic correction bispectrum (effectively the CMB four-point function) */
    if (strstr(string1,"quadratic") != NULL) {
      pbi->has_quadratic_correction = _TRUE_;
    }

    /* Intrinsic bispectrum. This is induced by second-order effects in the evolution of the cosmological
    perturbations. */
    if (strstr(string1,"intrinsic") != NULL) {
      pbi->has_intrinsic = _TRUE_;
      pbi->has_quadratic_correction = _TRUE_;
      ppt->has_perturbations2 = _TRUE_;
      ppt->has_cls = _TRUE_;
      /* We also require to be able to compute the analytical approximation in the squeezed limit,
      which is used as a window function for the interpolation of the intrinsic bispectrum */
      ppt->has_cl_cmb_zeta = _TRUE_;
      psp->compute_cl_derivative = _TRUE_;
    }

  } // end of bispectrum_types parsing

  // *** END OF MY MODIFICATIONS ***





  /* Set the modes (has_scalars...) and initial conditions */
  if (ppt->has_perturbations == _TRUE_) {

    class_call(parser_read_string(pfc,"modes",&string1,&flag1,errmsg),
	       errmsg,
	       errmsg);

    if (flag1 == _TRUE_) {

      /* if no modes are specified, the default is has_scalars=_TRUE_;
	 but if they are specified we should reset has_scalars to _FALSE_ before reading */
      ppt->has_scalars=_FALSE_;

      if ((strstr(string1,"s") != NULL) || (strstr(string1,"S") != NULL))
	ppt->has_scalars=_TRUE_;

      if ((strstr(string1,"v") != NULL) || (strstr(string1,"V") != NULL))
	ppt->has_vectors=_TRUE_;

      if ((strstr(string1,"t") != NULL) || (strstr(string1,"T") != NULL))
	ppt->has_tensors=_TRUE_;

      class_test(class_none_of_three(ppt->has_scalars,ppt->has_vectors,ppt->has_tensors),
		 errmsg,
		 "You wrote: modes=%s. Could not identify any of the modes ('s', 'v', 't') in such input",string1);
    }

    if (ppt->has_scalars == _TRUE_) {

      class_call(parser_read_string(pfc,"ic",&string1,&flag1,errmsg),
		 errmsg,
		 errmsg);

      if (flag1 == _TRUE_) {

	/* if no initial conditions are specified, the default is has_ad=_TRUE_;
	   but if they are specified we should reset has_ad to _FALSE_ before reading */
	ppt->has_ad=_FALSE_;
	ppt->has_ad_maberty=_FALSE_;
	ppt->has_zero_ic=_FALSE_;


	if ((strcmp(string1,"ad") == 0) || (strcmp(string1,"AD") == 0))
	  ppt->has_ad=_TRUE_;

  // *** MY MODIFICATIONS ***
	if ((strcmp(string1,"ad_maberty") == 0) || (strcmp(string1,"AD_MABERTY") == 0))
	  ppt->has_ad_maberty=_TRUE_;

	if ((strcmp(string1,"zero") == 0) || (strcmp(string1,"AD_ZERO") == 0))
	  ppt->has_zero_ic=_TRUE_;
  // *** END OF MY MODIFICATIONS ***

	if ((strcmp(string1,"bi") == 0) || (strcmp(string1,"BI") == 0))
	  ppt->has_bi=_TRUE_;

	if ((strcmp(string1,"cdi") == 0) || (strcmp(string1,"CDI") == 0))
	  ppt->has_cdi=_TRUE_;

	if ((strcmp(string1,"nid") == 0) || (strcmp(string1,"NID") == 0))
	  ppt->has_nid=_TRUE_;

	if ((strcmp(string1,"niv") == 0) || (strcmp(string1,"NIV") == 0))
	  ppt->has_niv=_TRUE_;

  // *** MY MODIFICATIONS ***
	class_test(ppt->has_ad==_FALSE_ && ppt->has_ad_maberty==_FALSE_ && ppt->has_zero_ic==_FALSE_ && ppt->has_bi ==_FALSE_ && ppt->has_cdi ==_FALSE_ && ppt->has_nid ==_FALSE_ && ppt->has_niv ==_FALSE_,
  // *** END OF MY MODIFICATIONS ***
		   errmsg,
		   "You wrote: ic=%s. Could not identify any of the initial conditions ('ad', 'ad_maberty', 'bi', 'cdi', 'nid', 'niv') in such input",string1);

      }
    }

    else {

      class_test(ppt->has_cl_cmb_lensing_potential == _TRUE_,
		 errmsg,
		 "Inconsistency: you want C_l's for cmb lensing potential, but no scalar modes\n");

      class_test(ppt->has_pk_matter == _TRUE_,
		 errmsg,
		 "Inconsistency: you want P(k) of matter, but no scalar modes\n");

    }

  }


  // *** Set the gauge
  class_call(parser_read_string(pfc,"gauge",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"newtonian") != NULL) || (strstr(string1,"Newtonian") != NULL) || (strstr(string1,"new") != NULL)) {
      ppt->gauge = 0;
    }

    if ((strstr(string1,"synchronous") != NULL) || (strstr(string1,"sync") != NULL) || (strstr(string1,"Synchronous") != NULL)) {
      ppt->gauge = 1;
    }
  }





  /** (d) define the primordial spectrum */

  class_call(parser_read_string(pfc,"P_k_ini type",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {
    flag2=_FALSE_;
    if (strcmp(string1,"analytic_Pk") == 0) {
      ppm->primordial_spec_type = analytic_Pk;
      flag2=_TRUE_;
    }
    class_test(flag2==_FALSE_,
	       errmsg,
	       "could not identify primordial spectrum type, check that it is one of 'analytic_pk', ...");
  }

  if (ppm->primordial_spec_type == analytic_Pk) {

    class_read_double("k_pivot",ppm->k_pivot);

    if (ppt->has_scalars == _TRUE_) {

      class_read_double("A_s",ppm->A_s);

      if (ppt->has_ad == _TRUE_) {

	class_read_double("n_s",ppm->n_s);
	class_read_double("alpha_s",ppm->alpha_s);

      }

      if (ppt->has_bi == _TRUE_) {

	class_read_double("f_bi",ppm->f_bi);
	class_read_double("n_bi",ppm->n_bi);
	class_read_double("alpha_bi",ppm->alpha_bi);

      }

      if (ppt->has_cdi == _TRUE_) {

	class_read_double("f_cdi",ppm->f_cdi);
	class_read_double("n_cdi",ppm->n_cdi);
	class_read_double("alpha_cdi",ppm->alpha_cdi);

      }

      if (ppt->has_nid == _TRUE_) {

	class_read_double("f_nid",ppm->f_nid);
	class_read_double("n_nid",ppm->n_nid);
	class_read_double("alpha_nid",ppm->alpha_nid);

      }

      if (ppt->has_niv == _TRUE_) {

	class_read_double("f_niv",ppm->f_niv);
	class_read_double("n_niv",ppm->n_niv);
	class_read_double("alpha_niv",ppm->alpha_niv);

      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_bi == _TRUE_)) {
	class_read_double_one_of_two("c_ad_bi","c_bi_ad",ppm->c_ad_bi);
	class_read_double_one_of_two("n_ad_bi","n_bi_ad",ppm->n_ad_bi);
	class_read_double_one_of_two("alpha_ad_bi","alpha_bi_ad",ppm->alpha_ad_bi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
	class_read_double_one_of_two("c_ad_cdi","c_cdi_ad",ppm->c_ad_cdi);
	class_read_double_one_of_two("n_ad_cdi","n_cdi_ad",ppm->n_ad_cdi);
	class_read_double_one_of_two("alpha_ad_cdi","alpha_cdi_ad",ppm->alpha_ad_cdi);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_nid == _TRUE_)) {
	class_read_double_one_of_two("c_ad_nid","c_nid_ad",ppm->c_ad_nid);
	class_read_double_one_of_two("n_ad_nid","n_nid_ad",ppm->n_ad_nid);
	class_read_double_one_of_two("alpha_ad_nid","alpha_nid_ad",ppm->alpha_ad_nid);
      }

      if ((ppt->has_ad == _TRUE_) && (ppt->has_niv == _TRUE_)) {
	class_read_double_one_of_two("c_ad_niv","c_niv_ad",ppm->c_ad_niv);
	class_read_double_one_of_two("n_ad_niv","n_niv_ad",ppm->n_ad_niv);
	class_read_double_one_of_two("alpha_ad_niv","alpha_niv_ad",ppm->alpha_ad_niv);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_cdi == _TRUE_)) {
	class_read_double_one_of_two("c_bi_cdi","c_cdi_bi",ppm->c_bi_cdi);
	class_read_double_one_of_two("n_bi_cdi","n_cdi_bi",ppm->n_bi_cdi);
	class_read_double_one_of_two("alpha_bi_cdi","alpha_cdi_bi",ppm->alpha_bi_cdi);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
	class_read_double_one_of_two("c_bi_nid","c_nid_bi",ppm->c_bi_nid);
	class_read_double_one_of_two("n_bi_nid","n_nid_bi",ppm->n_bi_nid);
	class_read_double_one_of_two("alpha_bi_nid","alpha_nid_bi",ppm->alpha_bi_nid);
      }

      if ((ppt->has_bi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
	class_read_double_one_of_two("c_bi_niv","c_niv_bi",ppm->c_bi_niv);
	class_read_double_one_of_two("n_bi_niv","n_niv_bi",ppm->n_bi_niv);
	class_read_double_one_of_two("alpha_bi_niv","alpha_niv_bi",ppm->alpha_bi_niv);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_nid == _TRUE_)) {
	class_read_double_one_of_two("c_cdi_nid","c_nid_cdi",ppm->c_cdi_nid);
	class_read_double_one_of_two("n_cdi_nid","n_nid_cdi",ppm->n_cdi_nid);
	class_read_double_one_of_two("alpha_cdi_nid","alpha_nid_cdi",ppm->alpha_cdi_nid);
      }

      if ((ppt->has_cdi == _TRUE_) && (ppt->has_niv == _TRUE_)) {
	class_read_double_one_of_two("c_cdi_niv","c_niv_cdi",ppm->c_cdi_niv);
	class_read_double_one_of_two("n_cdi_niv","n_niv_cdi",ppm->n_cdi_niv);
	class_read_double_one_of_two("alpha_cdi_niv","alpha_niv_cdi",ppm->alpha_cdi_niv);
      }

      if ((ppt->has_nid == _TRUE_) && (ppt->has_niv == _TRUE_)) {
	class_read_double_one_of_two("c_nid_niv","c_niv_nid",ppm->c_nid_niv);
	class_read_double_one_of_two("n_nid_niv","n_niv_nid",ppm->n_nid_niv);
	class_read_double_one_of_two("alpha_nid_niv","alpha_niv_nid",ppm->alpha_nid_niv);
      }

    }

    if (ppt->has_tensors == _TRUE_) {

      class_read_double("r",ppm->r);

      class_call(parser_read_string(pfc,"n_t",&string1,&flag1,errmsg),
		 errmsg,
		 errmsg);

      if (flag1 == _TRUE_) {

	if ((strstr(string1,"SCC") != NULL) || (strstr(string1,"scc") != NULL)) {
	  ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);
	}
	else {
	  class_read_double("n_t",ppm->n_t);
	}

      }

      class_call(parser_read_string(pfc,"alpha_t",&string1,&flag1,errmsg),
		 errmsg,
		 errmsg);

      if (flag1 == _TRUE_) {

	if ((strstr(string1,"SCC") != NULL) || (strstr(string1,"scc") != NULL)) {
	  ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);
	}
	else {
	  class_read_double("alpha_t",ppm->alpha_t);
	}

      }
    }
  }

  /** (e) parameters for final spectra */

  // *** MY MODIFICATIONS ***

  /* Read l_scalar_max also if ppt->has_cls is false, so that we can have a meaningful k-sampling when debugging the
    second-order early transfer functions using the has_early_transfers2_only flag. */

  if (ppt->has_scalars == _TRUE_) {
    class_read_double("l_max_scalars",ppt->l_scalar_max);
  }

  if (ppt->has_tensors == _TRUE_) {
    class_read_double("l_max_tensors",ppt->l_tensor_max);
  }
  // *** ORIGINAL CLASS

  // if (ppt->has_cls == _TRUE_) {
  //
  //   if (ppt->has_scalars == _TRUE_) {
  //     class_read_double("l_max_scalars",ppt->l_scalar_max);
  //   }
  //
  //   if (ppt->has_tensors == _TRUE_) {
  //     class_read_double("l_max_tensors",ppt->l_tensor_max);
  //   }
  // }

  // *** END OF MY MODIFICATIONS ***


  if ((ppt->has_scalars == _TRUE_) &&
      (ppt->has_cls == _TRUE_) &&
      (ppt->has_cl_cmb_lensing_potential == _TRUE_)) {

    class_call(parser_read_string(pfc,
				  "lensing",
				  &(string1),
				  &(flag1),
				  errmsg),
	       errmsg,
	       errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
      ple->has_lensed_cls = _TRUE_;
    }
  }

  if (ppt->has_pk_matter == _TRUE_) {

    class_call(parser_read_double(pfc,"P_k_max_h/Mpc",&param1,&flag1,errmsg),
	       errmsg,
	       errmsg);
    class_call(parser_read_double(pfc,"P_k_max_1/Mpc",&param2,&flag2,errmsg),
	       errmsg,
	       errmsg);
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),
	       errmsg,
	       "In input file, you cannot enter both P_k_max_h/Mpc and P_k_max_1/Mpc, choose one");
    if (flag1 == _TRUE_) {
      ppt->k_scalar_kmax_for_pk=param1*pba->h;
    }
    if (flag2 == _TRUE_) {
      ppt->k_scalar_kmax_for_pk=param2;
    }

    class_call(parser_read_list_of_doubles(pfc,
					   "z_pk",
					   &(int1),
					   &(pointer1),
					   &flag1,
					   errmsg),
	       errmsg,
	       errmsg);

    if (flag1 == _TRUE_) {
      class_test(int1 > _Z_PK_NUM_MAX_,
		 errmsg,
		 "you want to write some output for %d different values of z, hence you should increase _Z_PK_NUM_MAX_ in include/output.h to at least this number",
		 int1);
      pop->z_pk_num = int1;
      for (i=0; i<int1; i++) {
	pop->z_pk[i] = pointer1[i];
      }
      free(pointer1);
    }

    class_call(parser_read_double(pfc,"z_max_pk",&param1,&flag1,errmsg),
	       errmsg,
	       errmsg);

    if (flag1==_TRUE_) {
      psp->z_max_pk = param1;
    }
    else {
      psp->z_max_pk = 0.;
      for (i=0; i<pop->z_pk_num; i++)
	psp->z_max_pk = MAX(psp->z_max_pk,pop->z_pk[i]);
    }
  }

  /* deal with selection functions */
  if(ppt->has_cl_density == _TRUE_) {

    class_call(parser_read_string(pfc,
				  "selection",
				  &(string1),
				  &(flag1),
				  errmsg),
	       errmsg,
	       errmsg);

    if (flag1 == _TRUE_) {
      if (strstr(string1,"gaussian") != NULL) {
	ppt->selection=gaussian;
      }
      else if (strstr(string1,"tophat") != NULL) {
	ppt->selection=tophat;
      }
      else {
	class_stop("In selection function input: type is unclear","");
      }
    }

    class_call(parser_read_list_of_doubles(pfc,
					   "selection_mean",
					   &(int1),
					   &(pointer1),
					   &flag1,
					   errmsg),
	       errmsg,
	       errmsg);

    if ((flag1 == _TRUE_) && (int1>0)) {

      class_test(int1 > _SELECTION_NUM_MAX_,
		 errmsg,
		 "you want to compute density Cl's for %d different bins, hence you should increase _SELECTION_NUM_MAX_ in include/transfer.h to at least this number",
		 int1);

      ppt->selection_num = int1;
      for (i=0; i<int1; i++) {
	class_test((pointer1[i] < 0.) || (pointer1[i] > 1000.),
		   errmsg,
		   "input of selection functions: you asked for a mean redshift equal to %e, sounds odd",
		   pointer1[i]);
	ppt->selection_mean[i] = pointer1[i];
      }
      free(pointer1);
      /* first set all widths to default; correct eventually later */
      for (i=1; i<int1; i++) {
	class_test(ppt->selection_mean[i]<=ppt->selection_mean[i-1],
		   errmsg,
		   "input of selection functions: the list of mean redshifts must be passed in growing order; you entered %e before %e",ppt->selection_mean[i-1],ppt->selection_mean[i]);
	ppt->selection_width[i] = ppt->selection_width[0];
      }

      class_call(parser_read_list_of_doubles(pfc,
					     "selection_width",
					     &(int1),
					     &(pointer1),
					     &flag1,
					     errmsg),
	       errmsg,
	       errmsg);

      if ((flag1 == _TRUE_) && (int1>0)) {

	if (int1==1) {
	  for (i=0; i<ppt->selection_num; i++) {
	    ppt->selection_width[i] = pointer1[0];
	  }
	}
	else if (int1==ppt->selection_num) {
	  for (i=0; i<int1; i++) {
	    ppt->selection_width[i] = pointer1[i];
	  }
	}
	else {
	  class_stop(ptr->error_message,
		     "In input for selection function, you asked for %d bin centers and %d bin widths; number of bins unclear; you should pass either one bin width (common to all bins) or %d bin witdths",
		     ppt->selection_num,int1,ppt->selection_num);
	}
      }
    }
  }

  class_read_string("root",pop->root);

  class_call(parser_read_string(pfc,
				"headers",
				&(string1),
				&(flag1),
				errmsg),
	     errmsg,
	     errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL))) {
    pop->write_header = _FALSE_;
  }

  class_call(parser_read_string(pfc,"format",&string1,&flag1,errmsg),
	       errmsg,
	       errmsg);

  if (flag1 == _TRUE_) {

      if ((strstr(string1,"class") != NULL) || (strstr(string1,"CLASS") != NULL))
	pop->output_format = class_format;
      else {
	if ((strstr(string1,"camb") != NULL) || (strstr(string1,"CAMB") != NULL))
	  pop->output_format = camb_format;
	else
	  class_stop(errmsg,
		     "You wrote: format=%s. Could not identify any of the possible formats ('class', 'CLASS', 'camb', 'CAMB')",string1);
      }
  }

  class_call(parser_read_string(pfc,
				"bessel file",
				&(string1),
				&(flag1),
				errmsg),
	     errmsg,
	     errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {
    pbs->bessel_always_recompute = _FALSE_;
  }

  /** (f) parameter related to the non-linear spectra computation */

  class_call(parser_read_string(pfc,
				"non linear",
				&(string1),
				&(flag1),
				errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"halofit") != NULL) || (strstr(string1,"Halofit") != NULL) || (strstr(string1,"HALOFIT") != NULL)) {
      pnl->method=nl_halofit;
    }
    if ((strstr(string1,"trg") != NULL) || (strstr(string1,"TRG") != NULL)) {
      pnl->method=nl_trg;
    }
    if ((strstr(string1,"one-loop") != NULL) || (strstr(string1,"oneloop") != NULL) || (strstr(string1,"one loop") != NULL)) {
      pnl->method=nl_trg_one_loop;
    }
    if ((strstr(string1,"test linear") != NULL) || (strstr(string1,"test-linear") != NULL)) {
      pnl->method=nl_trg_linear;
    }
  }

  class_test((pnl->method>nl_none) && (ppt->has_pk_matter==_FALSE_),
	     errmsg,
	     "it is not consistent to ask for non-linear power spectrum but not for linear one: you should include mPk in the 'output' entry");

  if (pnl->method==nl_trg) {

    class_call(parser_read_string(pfc,
				  "non linear ic",
				  &(string1),
				  &(flag1),
				  errmsg),
	       errmsg,
	       errmsg);

    if ((strstr(string1,"linear") != NULL) || (strstr(string1,"lin") != NULL)) {
      pnl->ic=nl_lin;
    }
  }

  /** (g) amount of information sent to standard output (none if all set to zero) */

  class_read_int("background_verbose",
		 pba->background_verbose);

  class_read_int("thermodynamics_verbose",
		 pth->thermodynamics_verbose);

  class_read_int("perturbations_verbose",
		 ppt->perturbations_verbose);

  class_read_int("bessels_verbose",
		 pbs->bessels_verbose);

  class_read_int("transfer_verbose",
		 ptr->transfer_verbose);

  class_read_int("primordial_verbose",
		 ppm->primordial_verbose);

  class_read_int("spectra_verbose",
		 psp->spectra_verbose);

  class_read_int("nonlinear_verbose",
		 pnl->nonlinear_verbose);

  class_read_int("lensing_verbose",
		 ple->lensing_verbose);

  class_read_int("output_verbose",
		 pop->output_verbose);




  /** (h) all precision parameters */

  /** h.1. parameters related to the background */

  class_read_double("a_ini_over_a_today_default",ppr->a_ini_over_a_today_default);
  class_read_double("back_integration_stepsize",ppr->back_integration_stepsize);
  class_read_double("tol_background_integration",ppr->tol_background_integration);
  class_read_double("tol_initial_Omega_r",ppr->tol_initial_Omega_r);
  class_read_double("tol_ncdm_initial_w",ppr->tol_ncdm_initial_w);

  /** h.2. parameters related to the thermodynamics */

  class_read_string("sBBN file",ppr->sBBN_file);

  class_read_double("recfast_z_initial",ppr->recfast_z_initial);

  class_read_int("recfast_Nz0",ppr->recfast_Nz0);
  class_read_double("tol_thermo_integration",ppr->tol_thermo_integration);

  class_read_int("recfast_Heswitch",ppr->recfast_Heswitch);
  class_read_double("recfast_fudge_He",ppr->recfast_fudge_He);

  class_read_int("recfast_Hswitch",ppr->recfast_Hswitch);
  class_read_double("recfast_fudge_H",ppr->recfast_fudge_H);
  if (ppr->recfast_Hswitch == _TRUE_) {
    class_read_double("recfast_delta_fudge_H",ppr->recfast_delta_fudge_H);
    class_read_double("recfast_AGauss1",ppr->recfast_AGauss1);
    class_read_double("recfast_AGauss2",ppr->recfast_AGauss2);
    class_read_double("recfast_zGauss1",ppr->recfast_zGauss1);
    class_read_double("recfast_zGauss2",ppr->recfast_zGauss2);
    class_read_double("recfast_wGauss1",ppr->recfast_wGauss1);
    class_read_double("recfast_wGauss2",ppr->recfast_wGauss2);
  }

  class_read_double("recfast_z_He_1",ppr->recfast_z_He_1);
  class_read_double("recfast_delta_z_He_1",ppr->recfast_delta_z_He_1);
  class_read_double("recfast_z_He_2",ppr->recfast_z_He_2);
  class_read_double("recfast_delta_z_He_2",ppr->recfast_delta_z_He_2);
  class_read_double("recfast_z_He_3",ppr->recfast_z_He_3);
  class_read_double("recfast_delta_z_He_3",ppr->recfast_delta_z_He_3);
  class_read_double("recfast_x_He0_trigger",ppr->recfast_x_He0_trigger);
  class_read_double("recfast_x_He0_trigger2",ppr->recfast_x_He0_trigger2);
  class_read_double("recfast_x_He0_trigger_delta",ppr->recfast_x_He0_trigger_delta);
  class_read_double("recfast_x_H0_trigger",ppr->recfast_x_H0_trigger);
  class_read_double("recfast_x_H0_trigger2",ppr->recfast_x_H0_trigger2);
  class_read_double("recfast_x_H0_trigger_delta",ppr->recfast_x_H0_trigger_delta);
  class_read_double("recfast_H_frac",ppr->recfast_H_frac);

  class_read_string("Alpha_inf hyrec file",ppr->hyrec_Alpha_inf_file);
  class_read_string("R_inf hyrec file",ppr->hyrec_R_inf_file);
  class_read_string("two_photon_tables hyrec file",ppr->hyrec_two_photon_tables_file);

  class_read_double("reionization_z_start_max",ppr->reionization_z_start_max);
  class_read_double("reionization_sampling",ppr->reionization_sampling);
  class_read_double("reionization_optical_depth_tol",ppr->reionization_optical_depth_tol);
  class_read_double("reionization_start_factor",ppr->reionization_start_factor);

  class_read_int("thermo_rate_smoothing_radius",ppr->thermo_rate_smoothing_radius);

  /** h.3. parameters related to the perturbations */

  class_read_int("evolver",ppr->evolver);
  class_read_int("pk_definition",ppr->pk_definition);
  class_read_double("k_scalar_min_tau0",ppr->k_scalar_min_tau0);
  class_read_double("k_scalar_max_tau0_over_l_max",ppr->k_scalar_max_tau0_over_l_max);
  class_read_double("k_scalar_step_sub",ppr->k_scalar_step_sub);
  class_read_double("k_scalar_step_super",ppr->k_scalar_step_super);
  class_read_double("k_scalar_step_transition",ppr->k_scalar_step_transition);
  class_read_double("k_scalar_k_per_decade_for_pk",ppr->k_scalar_k_per_decade_for_pk);
  class_read_double("k_scalar_k_per_decade_for_bao",ppr->k_scalar_k_per_decade_for_bao);
  class_read_double("k_scalar_bao_center",ppr->k_scalar_bao_center);
  class_read_double("k_scalar_bao_width",ppr->k_scalar_bao_width);
  class_read_double("k_tensor_min_tau0",ppr->k_tensor_min_tau0);
  class_read_double("k_tensor_max_tau0_over_l_max",ppr->k_tensor_max_tau0_over_l_max);
  class_read_double("k_tensor_step_sub",ppr->k_tensor_step_sub);
  class_read_double("k_tensor_step_super",ppr->k_tensor_step_super);
  class_read_double("k_tensor_step_transition",ppr->k_tensor_step_transition);
  class_read_double("start_small_k_at_tau_c_over_tau_h",ppr->start_small_k_at_tau_c_over_tau_h);
  class_read_double("start_large_k_at_tau_h_over_tau_k",ppr->start_large_k_at_tau_h_over_tau_k);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_h",ppr->tight_coupling_trigger_tau_c_over_tau_h);
  class_read_double("tight_coupling_trigger_tau_c_over_tau_k",ppr->tight_coupling_trigger_tau_c_over_tau_k);
  class_read_double("start_sources_at_tau_c_over_tau_h",ppr->start_sources_at_tau_c_over_tau_h);

  class_read_int("tight_coupling_approximation",ppr->tight_coupling_approximation);

  /** derivatives of baryon sound speed only computed if some non-minimal tight-coupling schemes is requested */
  if ((ppr->tight_coupling_approximation == (int)first_order_CLASS) || (ppr->tight_coupling_approximation == (int)second_order_CLASS)) {
    pth->compute_cb2_derivatives = _TRUE_;
  }

  class_read_int("l_max_g",ppr->l_max_g);
  class_read_int("l_max_pol_g",ppr->l_max_pol_g);
  class_read_int("l_max_ur",ppr->l_max_ur);
  if (pba->N_ncdm>0)
    class_read_int("l_max_ncdm",ppr->l_max_ncdm);
  class_read_int("l_max_g_ten",ppr->l_max_g_ten);
  class_read_int("l_max_pol_g_ten",ppr->l_max_pol_g_ten);

  class_read_double("curvature_ini",ppr->curvature_ini);

  class_read_double("entropy_ini",ppr->entropy_ini);
  class_read_double("gw_ini",ppr->gw_ini);
  class_read_double("perturb_integration_stepsize",ppr->perturb_integration_stepsize);
  class_read_double("tol_tau_approx",ppr->tol_tau_approx);
  class_read_double("tol_perturb_integration",ppr->tol_perturb_integration);
  class_read_double("perturb_sampling_stepsize",ppr->perturb_sampling_stepsize);

  class_read_double("selection_cut_at_sigma",ppr->selection_cut_at_sigma);
  class_read_double("l_switch_limber_for_cl_density_over_z",ppr->l_switch_limber_for_cl_density_over_z);

  class_read_int("radiation_streaming_approximation",ppr->radiation_streaming_approximation);
  class_read_double("radiation_streaming_trigger_tau_over_tau_k",ppr->radiation_streaming_trigger_tau_over_tau_k);
  class_read_double("radiation_streaming_trigger_tau_c_over_tau",ppr->radiation_streaming_trigger_tau_c_over_tau);

  class_read_int("ur_fluid_approximation",ppr->ur_fluid_approximation);
  class_read_int("ncdm_fluid_approximation",ppr->ncdm_fluid_approximation);
  class_read_double("ur_fluid_trigger_tau_over_tau_k",ppr->ur_fluid_trigger_tau_over_tau_k);
  class_read_double("ncdm_fluid_trigger_tau_over_tau_k",ppr->ncdm_fluid_trigger_tau_over_tau_k);

  class_test(ppr->ur_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
	     errmsg,
	     "please choose different values for precision parameters ur_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

  if (pba->N_ncdm>0) {

    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->radiation_streaming_trigger_tau_over_tau_k,
	       errmsg,
	       "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

    class_test(ppr->ncdm_fluid_trigger_tau_over_tau_k==ppr->ur_fluid_trigger_tau_over_tau_k,
	       errmsg,
	       "please choose different values for precision parameters ncdm_fluid_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to avoid switching two approximation schemes at the same time");

  }






  /** h.4. parameter related to the Bessel functions */

  class_read_double("l_logstep",ppr->l_logstep);
  class_read_int("l_linstep",ppr->l_linstep);
  class_read_double("bessel_x_step",ppr->bessel_x_step);
  class_read_double("bessel_j_cut",ppr->bessel_j_cut);
  class_read_double("bessel_tol_x_min",ppr->bessel_tol_x_min);
  class_read_string("bessel_file_name",ppr->bessel_file_name);

  /** h.5. parameter related to the primordial spectra */

  class_read_double("k_per_decade_primordial",ppr->k_per_decade_primordial);

  /** h.6. parameter related to the transfer functions */

  class_read_double("k_step_trans_scalars",ppr->k_step_trans_scalars);
  class_read_double("k_step_trans_tensors",ppr->k_step_trans_tensors);
  class_read_int("transfer_cut",ppr->transfer_cut);
  class_read_double("transfer_cut_threshold_osc",ppr->transfer_cut_threshold_osc);
  class_read_double("transfer_cut_threshold_cl",ppr->transfer_cut_threshold_cl);
  class_read_double("l_switch_limber",ppr->l_switch_limber);

  /** h.7. parameters related to nonlinear calculations */

  class_read_double("halofit_dz",ppr->halofit_dz);
  class_read_double("halofit_min_k_nonlinear",ppr->halofit_min_k_nonlinear);
  class_read_double("halofit_sigma_precision",ppr->halofit_sigma_precision);

  class_read_int("double escape",ppr->double_escape);
  class_read_double("z_ini",ppr->z_ini);
  class_read_int("eta_size",ppr->eta_size);
  class_read_double("k_L",ppr->k_L);
  class_read_double("k_min",ppr->k_min);
  class_read_double("logstepx_min",ppr->logstepx_min);
  class_read_double("logstepk1",ppr->logstepk1);
  class_read_double("logstepk2",ppr->logstepk2);
  class_read_double("logstepk3",ppr->logstepk3);
  class_read_double("logstepk4",ppr->logstepk4);
  class_read_double("logstepk5",ppr->logstepk5);
  class_read_double("logstepk6",ppr->logstepk6);
  class_read_double("logstepk7",ppr->logstepk7);
  class_read_double("logstepk8",ppr->logstepk8);
  class_read_double("k_growth_factor",ppr->k_growth_factor);
  class_read_double("k_scalar_max_for_pk_nl",ppr->k_scalar_max_for_pk_nl);

  if ((pnl->method==nl_trg_one_loop) ||
      (pnl->method==nl_trg)) {

    /* when using the trg module, the following parameters need to
       be changed */

    ppt->k_scalar_kmax_for_pk
      = MAX(
	    ppt->k_scalar_kmax_for_pk,
	    ppr->k_scalar_max_for_pk_nl*pba->h);

    psp->z_max_pk = ppr->z_ini+1.;

  }

  /** h.8. parameter related to lensing */

  class_read_int("accurate_lensing",ppr->accurate_lensing);
  class_read_int("delta_l_max",ppr->delta_l_max);
  if (ppr->accurate_lensing == _TRUE_) {
    class_read_int("num_mu_minus_lmax",ppr->num_mu_minus_lmax);
    class_read_int("tol_gauss_legendre",ppr->tol_gauss_legendre);
  }



  /* check various l_max */

  pbs->l_max=0;
  pbs->x_max=0;

  // *** MY MODIFICATIONS
  if ((ppt->has_cls == _TRUE_) || (pbi->has_bispectra==_TRUE_)) {
  // *** PREVIOUSLY
  // if (ppt->has_cls == _TRUE_) {
  // *** END OF MY MODIFICATIONS

    if (ppt->has_scalars == _TRUE_) {

      /**if (ppt->has_cl_cmb_temperature == _TRUE_ ||
          ppt->has_cl_cmb_polarization == _TRUE_||
          ppt->has_cl_cmb_lensing_potential == _TRUE_)*/
      pbs->l_max=MAX(ppt->l_scalar_max,pbs->l_max);

      pbs->x_max=MAX(pbs->l_max*ppr->k_scalar_max_tau0_over_l_max,pbs->x_max);

    }

    if (ppt->has_tensors == _TRUE_) {
      pbs->l_max=MAX(ppt->l_tensor_max,pbs->l_max);

      pbs->x_max=MAX(pbs->l_max*ppr->k_tensor_max_tau0_over_l_max,pbs->x_max);
    }
  }

  pbs->x_step = ppr->bessel_x_step;

  pbs->x_max = ((int)(pbs->x_max * 1.01 / pbs->x_step)+1)*pbs->x_step;



  /** (i) eventually write all the read parameters in a file */

  class_call(parser_read_string(pfc,"write parameters",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))) {

    sprintf(param_output_name,"%s%s",pop->root,"parameters.ini");
    sprintf(param_unused_name,"%s%s",pop->root,"unused_parameters");

    class_open(param_output,param_output_name,"w",errmsg);
    class_open(param_unused,param_unused_name,"w",errmsg);

    fprintf(param_output,"# List of input/precision parameters actually read\n");
    fprintf(param_output,"# (all other parameters set to default values)\n");
    fprintf(param_output,"# Obtained with CLASS %s (for developpers: svn version %s)\n",_VERSION_,_SVN_VERSION_);
    fprintf(param_output,"#\n");
    fprintf(param_output,"# This file can be used as the input file of another run\n");
    fprintf(param_output,"#\n");

    fprintf(param_unused,"# List of input/precision parameters passed\n");
    fprintf(param_unused,"# but not used (just for info)\n");
    fprintf(param_unused,"#\n");

    for (i=0; i<pfc->size; i++) {
      if (pfc->read[i] == _TRUE_)
	fprintf(param_output,"%s = %s\n",pfc->name[i],pfc->value[i]);
      else
	fprintf(param_unused,"%s = %s\n",pfc->name[i],pfc->value[i]);
    }
    fprintf(param_output,"#\n");

    fclose(param_output);
    fclose(param_unused);
  }






  // *** MY MODIFICATIONS ***


  // ==========================================================
  // =                     Perturbations                      =
  // ==========================================================

  /* Extra k-sampling parameter to better sample to largest scales */
  class_read_double("k_scalar_logstep_super",ppr->k_scalar_logstep_super);


  // ******             First-order LOS effects             ******

  // *** Set ppt->has_scattering_in_los
  class_call(parser_read_string(pfc,"include_scattering_in_los_1st_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt->has_scattering_in_los = _FALSE_;
  }

  // *** Set ppt->has_photon_monopole_in_los
  class_call(parser_read_string(pfc,"include_photon_monopole_in_los_1st_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt->has_photon_monopole_in_los = _FALSE_;
  }

  // *** Set ppt->has_metric_in_los
  class_call(parser_read_string(pfc,"include_metric_in_los_1st_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") == NULL) && (strstr(string1,"Y") == NULL)) {
    ppt->has_metric_in_los = _FALSE_;
  }



  // *** Set ppt->has_sw
  class_call(parser_read_string(pfc,"include_sachs_wolfe_in_los_1st_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)) {
    ppt->has_sw = _TRUE_;
  }

  // *** Set ppt->has_isw
  class_call(parser_read_string(pfc,"include_integrated_sachs_wolfe_in_los_1st_order",&(string1),&(flag1),errmsg),errmsg,errmsg);

  if ((flag1 == _TRUE_) && (strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)) {
    ppt->has_isw = _TRUE_;
  }

  /* Avoid counting twice the same metric effect */
  if ((ppt->has_sw == _TRUE_) || (ppt->has_isw == _TRUE_))
    ppt->has_metric_in_los = _FALSE_;

  /* If effects that are not peaked at recombination are included, we need to extend the integration range up to today */
  // if ((ppt->has_metric_in_los == _FALSE_) && (ppt->has_isw == _FALSE_))
  //   ppt->has_recombination_only = _TRUE_;


  /* The zeta-T correlation is not implemented yet in synchronous gauge */
  if (ppt->gauge == synchronous)
    ppt->has_cl_cmb_zeta = _FALSE_;

  /* Should the curvature perturbation zeta only include contributions from recombination?
  I.e. should we ignore reionisation when computing zeta? */
  if ((ppt->has_cl_cmb_zeta == _TRUE_) && (pth->reio_parametrization!=reio_none)) {

    class_call(parser_read_string(pfc,"recombination_only_zeta",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"n") != NULL) || (strstr(string1,"N") != NULL)))
      ppt->recombination_only_zeta = _FALSE_;
  }


  // *****      Tau sampling for quadratic sources      ******

  /* Frequency of sampling. This parameter is overridden if the user specifies a custom timesampling. */
  class_read_double("perturb_sampling_stepsize_for_quadsources", ppr->perturb_sampling_stepsize_quadsources);



  /* Do we need a custom time-sampling? */
  class_call(parser_read_string(pfc,"custom_time_sampling_for_quadsources",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppt->has_custom_timesampling_for_quadsources = _TRUE_;


  class_read_double("custom_tau_ini_quadsources", ppt->custom_tau_ini_quadsources);

  class_test (ppt->custom_tau_ini_quadsources<=0, errmsg, "please choose 'custom_tau_ini_quadsources' greater than zero.");

  class_read_double("custom_tau_end_quadsources", ppt->custom_tau_end_quadsources);

  class_read_int("custom_tau_size_quadsources", ppt->custom_tau_size_quadsources);

  class_call(parser_read_string(pfc,"custom_tau_mode_quadsources",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"lin") != NULL) || (strstr(string1,"LIN") != NULL)))
      ppt->custom_tau_mode_quadsources = lin_tau_sampling;

    else if (((strstr(string1,"log") != NULL) || (strstr(string1,"LOG") != NULL)))
      ppt->custom_tau_mode_quadsources = log_tau_sampling;

    else if (((strstr(string1,"class") != NULL) || (strstr(string1,"CLASS") != NULL)))
      ppt->custom_tau_mode_quadsources = class_tau_sampling;

    else
      class_test(1==1,
    		   errmsg,
    		   "You wrote: custom_tau_mode_quadsources=%s. Could not identify any of the supported time samplings ('lin', 'log', 'class') in such input",string1);

  }

  /* This is used only for the interpolation of ppt->quadsources when the custom timesampling is chosen
    to be logarithmic */
  ppt->custom_log_tau_ini_quadsources = log(ppt->custom_tau_ini_quadsources);

  /* Define time step, either linear or logarithmic (used only for interpolation purposes) */
  if (ppt->custom_tau_mode_quadsources == lin_tau_sampling)
    ppt->custom_tau_step_quadsources = (ppt->custom_tau_end_quadsources - ppt->custom_tau_ini_quadsources)/(ppt->custom_tau_size_quadsources-1);

  else if (ppt->custom_tau_mode_quadsources == log_tau_sampling)
    ppt->custom_tau_step_quadsources = log(ppt->custom_tau_end_quadsources/ppt->custom_tau_ini_quadsources)/(ppt->custom_tau_size_quadsources-1);



  /* Interpolation of ppt->quadsources */

  class_call(parser_read_string(pfc,"quadsources_time_interpolation",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr->quadsources_time_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr->quadsources_time_interpolation = cubic_interpolation;

    else
      class_test(1==1,
    		   errmsg,
    		   "You wrote: quadsources_time_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",string1);

  }





  // =================================================
  // =                     Bessels                   =
  // =================================================

  class_call(parser_read_string(pfc,"bessels_interpolation",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr->bessels_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr->bessels_interpolation = cubic_interpolation;

    else
      class_test(1==1,
    		   errmsg,
    		   "You wrote: bessels_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",string1);

  }




  // ===================================================
  // =                     Bispectra                   =
  // ===================================================

  class_read_int("bispectra_verbose",
    pbi->bispectra_verbose);




  /* Which kind of technique should we use for the k3-integration of the bispectrum? */
  class_call(parser_read_string(pfc,"bispectra_k3_extrapolation",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    if (strstr(string1,"no_extrapolation") != NULL)
      ppr->bispectra_k3_extrapolation = no_k3_extrapolation;

    else if (strstr(string1,"flat_extrapolation") != NULL)
      ppr->bispectra_k3_extrapolation = flat_k3_extrapolation;

    else if (strstr(string1,"linear_extrapolation") != NULL)
      ppr->bispectra_k3_extrapolation = linear_k3_extrapolation;

    else
      class_stop(errmsg,
        "Could not recognize the value given for the 'bispectra_k3_extrapolation' option. Choose between 'no_extrapolation', 'flat_extrapolation'\
or 'linear_extrapolation'.", "");

  }

  /* How much to extend the k3 range in view of the bispectrum integration? */
  class_read_double("extra_k3_oscillations_right", ppr->extra_k3_oscillations_right);
  class_read_double("extra_k3_oscillations_left", ppr->extra_k3_oscillations_left);



  // ****   l-interpolation of bispectra   ****

  class_call(parser_read_string(pfc,"bispectra_interpolation",&string1,&flag1,errmsg),
         errmsg,
         errmsg);

  if (flag1 == _TRUE_) {

    /* For trilinear interpolation, we need an all-even grid */
    if ((strstr(string1,"trilinear") != NULL) || (strstr(string1,"tri") != NULL)) {
      pfi->bispectra_interpolation = trilinear_interpolation;
      ppr->compute_only_even_ls = _TRUE_;
    }

    /* For trilinear interpolation, we need an all-even grid */
    else if ((strstr(string1,"smart") != NULL) || (strstr(string1,"SMART") != NULL)) {
      pfi->bispectra_interpolation = smart_interpolation;
      ppr->compute_only_even_ls = _TRUE_;
    }

    else if (strstr(string1,"sum") != NULL) {
      pfi->bispectra_interpolation = sum_over_all_multipoles;
    }
    else if ((strcmp(string1,"mesh_2d") == 0) || (strcmp(string1,"mesh_2D") == 0)
    || (strcmp(string1,"mesh") == 0) || (strcmp(string1,"mesh_3d") == 0) || (strcmp(string1,"mesh_3D") == 0)) {
      pfi->bispectra_interpolation = mesh_interpolation_2D;
    }
    else {
      class_test(1==1,
      errmsg,
      "You wrote: bispectra_interpolation=%s. Could not identify any of the supported interpolation techniques ('trilinear', 'mesh_3d', 'mesh_2d', 'sum') in such input",string1);
    }

  }

  // ****   Interpolation of transfer functions   ****

  class_call(parser_read_string(pfc,"transfers_k1_interpolation",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr->transfers_k1_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr->transfers_k1_interpolation = cubic_interpolation;

    else
      class_test(1==1,
    		   errmsg,
    		   "You wrote: transfers_k1_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",string1);

  }


  class_call(parser_read_string(pfc,"transfers_k2_interpolation",&string1,&flag1,errmsg),
	     errmsg,
	     errmsg);

  if (flag1 == _TRUE_) {

    if (((strstr(string1,"linear") != NULL) || (strstr(string1,"LINEAR") != NULL)))
      ppr->transfers_k2_interpolation = linear_interpolation;

    else if (((strstr(string1,"cubic") != NULL) || (strstr(string1,"CUBIC") != NULL) || (strstr(string1,"spline") != NULL) || (strstr(string1,"SPLINE") != NULL)))
      ppr->transfers_k2_interpolation = cubic_interpolation;

    else
      class_test(1==1,
    		   errmsg,
    		   "You wrote: transfers_k2_interpolation=%s. Could not identify any of the supported interpolation techniques ('linear', 'cubic') in such input",string1);

  }




  // ***  Update x_max  ***

  if (pbi->has_bispectra == _TRUE_) {

    /* Take all multipoles if requested */
    if (pfi->bispectra_interpolation == sum_over_all_multipoles)
      ppr->l_linstep = 1;

    /* Check that l_max is even for trilinear interpolation */
    if (pfi->bispectra_interpolation == trilinear_interpolation)
      class_test ((pbs->l_max%2)!=0,
        errmsg,
        "For trilinear interpolation of the bispectra, ensure that both 'l_max_scalars' and 'l_max_tensors' are even");

    /* Parameters related to the r-integration in the bispectrum integral  */
    pbi->r_min = 13000;
    class_read_double("r_min", pbi->r_min);

    double tau0_sup = 15000.;
    pbi->r_max = tau0_sup;
    class_read_double("r_max", pbi->r_max);

    pbi->r_size = 100;
    class_read_int("r_size", pbi->r_size);


    /* TODO: update x_max to take into account the k3 extrapolation (use following piece of code)*/
    // if (ppr->bispectra_k3_extrapolation != no_k3_extrapolation) {
    //
    //   /* The physical range is the one dictated by the triangular condition: k1 + k2 = k3 */
    //   double physical_k3_range = k_max_pt - k_min_pt;
    //
    //   /* The extended range allows for the Bessel function j(k3*(tau0 - tau_rec)) to develop at least
    //     ppr->extra_k3_oscillations oscillations in k3 in order to stabilize the bispectrum integral. */
    //   double extended_k3_range = 2.*_PI_*ppr->extra_k3_oscillations/(ptr2->tau0 - ptr2->tau_rec);
    //
    //   if (physical_k3_range < extended_k3_range)
    //     k_max_tr += (extended_k3_range - physical_k3_range);
    //
    //   // printf("physical_k3_range=%g, extended_k3_range=%g\n", physical_k3_range, extended_k3_range);
    //   // printf("PRE:  k_max_tr = %g\n", k_max_pt);
    //   // printf("POST: k_max_tr = %g\n", k_max_tr);
    // }


    /* Determine maximum k-value needed in the line of sight integration at first order. This is
      needed to determine the maximum sampling of the Bessel functions */
    double tau0_inf = 8000.;
    double k_max = ppr->k_scalar_max_tau0_over_l_max * ppt->l_scalar_max / tau0_inf;

    /* Maximum argument for the Bessel functions in the line of sight integration.  For the time being, we do
      set x_max  here, even if we shouldn't as tau0 cannot be accessed by this module.  This is why we define
      a superior limit for tau0. */
    double x_max = MAX (pbs->x_max, k_max * MAX(tau0_sup, pbi->r_max));

    // printf("# Temporary message: Setting pbs->x_max from %g to %g\n",
    //   ((int)(pbs->x_max * 1.1 / pbs->x_step)+1)*pbs->x_step,
    //   ((int)(x_max * 1.1 / pbs->x_step)+1)*pbs->x_step);

    pbs->x_max = ((int)(x_max * 1.1 / pbs->x_step)+1)*pbs->x_step;


  } // end of if(has_bispectra)


  // =====================================================================
  // =                              C_l's                                =
  // =====================================================================

  /* Should we include the lensing effects on the bispectrum and on the Fisher matrix
  estimator? If so, extend the lensed C_l's all the way to l_max by setting
  the flag 'ppr->extend_lensed_cls'. The default behaviour is to compute the lensed C_l's
  only up to l_max - ppr->delta_l_max. */
  if (ple->has_lensed_cls == _TRUE_) { /* equivalent to setting lensing=yes in parameter file */

    if (pfi->has_fisher == _TRUE_) {
      pfi->include_lensing_effects = _TRUE_;
      pbi->has_cmb_lensing_kernel = _TRUE_;  /* needed to compute the effect of lensing variance */
      ppr->extend_lensed_cls = _TRUE_;       /* needed to include the lensing noise in the covariance matrix */
    }

    if ((pbi->has_bispectra == _TRUE_) &&
    ((pbi->has_cmb_lensing == _TRUE_)              /* analytic bispectrum that needs the lensed C_l's */
    || (pbi->has_cmb_lensing_squeezed == _TRUE_)   /* analytic bispectrum that needs the lensed C_l's */
    || (pbi->has_cmb_lensing_kernel == _TRUE_)     /* analytic bispectrum that needs the lensed C_l's */
    || (pbi->has_intrinsic_squeezed == _TRUE_))) { /* analytic bispectrum that needs the lensed C_l's */
      pbi->include_lensing_effects = _TRUE_;
      pbi->lensed_intrinsic = _TRUE_;
      ppr->extend_lensed_cls = _TRUE_;
    }
  }

  /* If requested, do not include the lensed C_l's in the squeezed approximation for the intrinsic
  bispectrum */
  if (pbi->include_lensing_effects == _TRUE_) {

    class_call(parser_read_string(pfc,"lensed_intrinsic",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"n") != NULL) || (strstr(string1,"N") != NULL)))
      pbi->lensed_intrinsic = _FALSE_;
  }

  if (pfi->include_lensing_effects == _TRUE_) {

    class_call(parser_read_string(pfc,"fisher_lensvar_lmax",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);
    class_call(parser_read_string(pfc,"compute_lensing_variance_lmax",&(string2),&(flag2),errmsg), /* obsolete */
        errmsg,
        errmsg);

    /* string1 wins over string2 */
    if (flag1 == _TRUE_)
      strcpy (string, string1);
    else if (flag2 == _TRUE_)
      strcpy (string, string2);

    if (((flag1 == _TRUE_)||(flag2 == _TRUE_))
    &&((strstr(string,"y") != NULL) || (strstr(string,"Y") != NULL)))
      pfi->compute_lensing_variance_lmax = _TRUE_;

  }


  // =====================================================================
  // =                           Fisher matrix                           =
  // =====================================================================

  class_read_int("fisher_verbose",
    pfi->fisher_verbose);

  class_call(parser_read_string(pfc,"always_interpolate_bispectra",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    pfi->always_interpolate_bispectra = _TRUE_;

  /* Minimum and maximum multipole to consider in the Fisher sum */
  class_read_int("fisher_l_min",pfi->l_min_estimator);
  class_read_int("fisher_lmin",pfi->l_min_estimator);
  class_read_int("fisher_l_max",pfi->l_max_estimator);
  class_read_int("fisher_lmax",pfi->l_max_estimator);

  /* Read the experiment sky coverage. */
  class_read_double("experiment_f_sky", pfi->f_sky);

  /* Read the experiment beam at FWHM in arcminutes and convert it to radians. */
  class_call (parser_read_list_of_doubles (pfc,
           "experiment_beam_fwhm",
           &(int1),
           &(pointer1),
           &flag1,
           errmsg),
       errmsg,
       errmsg);

  if (flag1 == _TRUE_) {

    class_test(int1 > _N_FREQUENCY_CHANNELS_MAX_,
      errmsg,
      "you specified too many frequency bands for the experiment, increase _N_FREQUENCY_CHANNELS_MAX_ in include/fisher.h or specify \
less than %d values for 'experiment_beam_fwhm'", _N_FREQUENCY_CHANNELS_MAX_);

    pfi->n_channels = int1;

    for (i=0; i<int1; i++)
      pfi->beam[i] = pointer1[i]/60. * _PI_/180.;

    free(pointer1);
  }


  /* Read the TEMPERATURE noise for the experiment in uK*arcminutes, and convert it to the actual noise */
  if ((pfi->has_fisher == _TRUE_) && (ppt->has_bi_cmb_temperature == _TRUE_)) {

    class_call (parser_read_list_of_doubles (pfc,
                 "experiment_noise_t",
                 &(int1),
                 &(pointer1),
                 &flag1,
                 errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {

      class_test(int1 != pfi->n_channels,
        errmsg,
        "the number of entries for 'experiment_beam_fwhm' and 'experiment_noise_t' must be equal");

      /* Compute the actual noise, as in astro-ph/0506396v2. A negative value implies infinite noise,
      which is equivalent to ignoring the frequency channel for this field */
      for (i=0; i<int1; i++) {
        if (pointer1[i]>=0)
          pfi->noise_t[i] = pow(pointer1[i]/1e6*pfi->beam[i]/pba->T_cmb,2);
        else
          pfi->noise_t[i] = _HUGE_;
      }

      free(pointer1);
    }
  } // end of T noise

  /* Read the POLARIZATION noise for the experiment in uK, and convert it to Kelvins^2 */
  if ((pfi->has_fisher == _TRUE_) && (ppt->has_cl_cmb_polarization == _TRUE_)) {

    class_call (parser_read_list_of_doubles (pfc,
                 "experiment_noise_e",
                 &(int1),
                 &(pointer1),
                 &flag1,
                 errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {

      class_test(int1 != pfi->n_channels,
        errmsg,
        "the number of entries for 'experiment_beam_fwhm' and 'experiment_noise_e' must be equal");

      /* Compute the actual noise, as in astro-ph/0506396v2. A negative value implies infinite noise,
      which is equivalent to ignoring the frequency channel for this field */
      for (i=0; i<int1; i++) {
        if (pointer1[i]>=0)
          pfi->noise_e[i] = pow(pointer1[i]/1e6*pfi->beam[i]/pba->T_cmb,2);
        else
          pfi->noise_e[i] = _HUGE_;
      }

      free(pointer1);
    }
  } // end of E noise

  /* Read the RAYLEIGH noise for the experiment in uK*arcminutes, and convert it to the actual noise */
  if ((pfi->has_fisher == _TRUE_) && (ppt->has_cl_cmb_rayleigh == _TRUE_)) {

    class_call (parser_read_list_of_doubles (pfc,
                 "experiment_noise_r",
                 &(int1),
                 &(pointer1),
                 &flag1,
                 errmsg),
      errmsg,
      errmsg);

    if (flag1 == _TRUE_) {

      class_test(int1 != pfi->n_channels,
        errmsg,
        "the number of entries for 'experiment_beam_fwhm' and 'experiment_noise_r' must be equal");

      /* Compute the actual noise, as in astro-ph/0506396v2. A negative value implies infinite noise,
      which is equivalent to ignoring the frequency channel for this field */
      for (i=0; i<int1; i++) {
        if (pointer1[i]>=0)
          pfi->noise_r[i] = pow(pointer1[i]/1e6*pfi->beam[i]/pba->T_cmb,2);
        else
          pfi->noise_r[i] = _HUGE_;
      }

      free(pointer1);
    }
  } // end of T noise


  class_call(parser_read_string(pfc,"fisher_ignore",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);
  class_call(parser_read_string(pfc,"ignored_fields_in_fisher",&(string2),&(flag2),errmsg), /* osbolete */
      errmsg,
      errmsg);

  /* string1 wins over string2 */
  if (flag1 == _TRUE_)
    strcpy (string, string1);
  else if (flag2 == _TRUE_)
    strcpy (string, string2);

  if ((flag1 == _TRUE_) || (flag2 == _TRUE_)) {

    if ((strstr(string,"t") != NULL) || (strstr(string,"T") != NULL))
      pfi->ignore_t = _TRUE_;

    if ((strstr(string,"e") != NULL) || (strstr(string,"E") != NULL))
      pfi->ignore_e = _TRUE_;

    if ((strstr(string,"b") != NULL) || (strstr(string,"B") != NULL))
      pfi->ignore_b = _TRUE_;

    if ((strstr(string,"r") != NULL) || (strstr(string,"R") != NULL))
      pfi->ignore_r = _TRUE_;

  }

  class_read_double("fisher_squeezed_ratio",pfi->squeezed_ratio);

  // ==========================================================================================
  // =                                 Create run directory                                   =
  // ==========================================================================================

  /* The run directory 'ppr->run_dir' is the directory where the parameter files
  (run_params.ini and run_params.pre), data folders (sources, transfers, bispectra) and
  result files (cl.dat, fisher.dat, etc.) will be stored. */

  class_call(parser_read_string(pfc,"store_run",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppr->store_run = _TRUE_;

  /* Create a run directory if asked (ppr->store_run == _TRUE_) and if it does not already exist
  (ppr->load_run == _FALSE_). The latter case is always true unless we are running CLASS with a single
  argument and that argument is the run directory. */
  if ((ppr->store_run == _TRUE_) && (ppr->load_run == _FALSE_)) {

    /* Read the name of the run directory from the parameter file */
    class_call(parser_read_string(pfc,
         "run_directory",
         &(string1),
         &(flag1),
         errmsg),
         errmsg,
         errmsg);

    if ((flag1 == _TRUE_) && (string1 != NULL)) {

      /* Expand shell symbols (such as ~ and ..) and environment variables in the path */
      wordexp_t exp_result;
      class_test (wordexp(string1, &exp_result, 0)!=0, errmsg, "error in word expansion");
      strcpy(string1, exp_result.we_wordv[0]);
      wordfree(&exp_result);

      strcpy(ppr->run_dir, string1);

      /* By default, we set the data_dir to be the same as the run_dir, that is, we
      store/read the data to/from the same directory where we store/read the
      parameter files and the result files. */
      strcpy(ppr->data_dir, string1);
    }

    /* Check that the directory does not exist */
    struct stat st;
    stat (ppr->run_dir, &st);
    class_test (S_ISDIR (st.st_mode) != 0,
      errmsg,
      "target directory '%s' already exists, choose another one", ppr->run_dir);

    /* Should the date be appended to the run directory? */
    class_call(parser_read_string(pfc,"append_date",&(string1),&(flag1),errmsg),
        errmsg,
        errmsg);

    if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
      ppr->append_date_to_run = _TRUE_;


    /* Determine the actual run directory according to the append_date variable */
    if (ppr->append_date_to_run == _TRUE_) {

      time_t now;
      struct tm *d;
      char suffix[64];

      time(&now);
      d = localtime(&now);

      strftime(suffix, 64, "%Y-%m-%d_%H-%M-%S", d);
      sprintf(ppr->run_dir, "%s_%s", ppr->run_dir, suffix);
    }

    /* Create the directory for the run */
    class_test (mkdir (ppr->run_dir, 0777)!=0,
      errmsg,
      "could not create run directory '%s'. Parent directory doesn't exist? Actual directory already exists?",
      ppr->run_dir);

    /* Print some information to screen */
    printf("# We shall store the current run to the folder %s.\n", ppr->run_dir);

    /* Copy the parameter files to the run directory */
    char new_ini_filepath[_FILENAMESIZE_], new_pre_filepath[_FILENAMESIZE_], command[3*_FILENAMESIZE_];

    sprintf(new_ini_filepath, "%s/run_params.ini", ppr->run_dir);
    sprintf(command, "cp %s %s", ppr->ini_filename, new_ini_filepath);
    if (strcmp(ppr->ini_filename, new_ini_filepath) != 0) system(command);

    sprintf(new_pre_filepath, "%s/run_params.pre", ppr->run_dir);
    sprintf(command, "cp %s %s", ppr->pre_filename, new_pre_filepath);
    if (strcmp(ppr->pre_filename, new_pre_filepath) != 0) system(command);

  } // end of if not load_run

  /* In any case, store or load run, make the root coincide with the run directory, so that the output
  files (cl.dat, fisher.dat, etc.) will be dumped there */
  if ((ppr->store_run == _TRUE_) || (ppr->load_run == _TRUE_))
    sprintf (pop->root, "%s/", ppr->run_dir);

  /* Set an environment variable for the run directory, so that it can be used in the
  parameter file by the user, for example to set the location of the data directory */
  setenv ("SONG_RUN", ppr->run_dir, 1);
  setenv ("SONG_RUN_PATH", ppr->run_dir, 1);
  setenv ("SONG_RUN_DIR", ppr->run_dir, 1);
  setenv ("SONG_RUN_FOLDER", ppr->run_dir, 1);


  // ============================================================================================
  // =                                   Read data directory                                    =
  // ============================================================================================

  /* The data directory 'ppr->data_dir' is the directory where the data folders (sources, transfers,
  bispectra) will be read from. By default it coincides with the run directory 'ppr->data_dir' */

  class_call(parser_read_string(pfc,"data_directory",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);

  if ((flag1 == _TRUE_) && (string1 != NULL)) {

    /* Expand shell symbols (such as ~ and ..) and environment variables in the path */
    wordexp_t exp_result;
    class_test (wordexp(string1, &exp_result, 0)!=0, errmsg, "error in word expansion");
    strcpy(string1, exp_result.we_wordv[0]);
    wordfree(&exp_result);

    /* Check that the data directory exists, but only if it is going to be needed */
    struct stat st;
    stat (string1, &st);
    class_test (S_ISDIR (st.st_mode) == 0,
  		errmsg,
  		"the data directory does not exist, change the parameter 'data_dir' ('%s')", string1);

    strcpy (ppr->data_dir, string1);
  }



  // =============================================================================================
  // =                               Disk storage of bispectra                                   =
  // =============================================================================================

  /* Store to disk the bispectra? */
  class_call(parser_read_string(pfc,"store_bispectra",&(string1),&(flag1),errmsg),
      errmsg,
      errmsg);

  if ((flag1 == _TRUE_) && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)))
    ppr->store_bispectra_to_disk = _TRUE_;

  sprintf(pbi->bispectra_dir, "%s/bispectra", ppr->data_dir);

  /* If we are not loading from disk, just create the bispectra directory */
  if ((ppr->store_bispectra_to_disk == _TRUE_) && (ppr->load_run == _FALSE_)) {

    class_test (mkdir (pbi->bispectra_dir, 0777) != 0,
      errmsg,
      "could not create directory '%s', maybe it already exists?", pbi->bispectra_dir);
  }
  /* If we are in a run directory, checks if it already contains the bispectra */
  else if (ppr->load_run == _TRUE_) {

    struct stat st;
    short bispectra_dir_exists = (stat(pbi->bispectra_dir, &st)==0);

    /* If the bispectra directory exists, then we shall load the 2nd-order bispectra from it */
    if (bispectra_dir_exists) {
      ppr->store_bispectra_to_disk = _FALSE_;
      ppr->load_bispectra_from_disk = _TRUE_;
      if (pbi->bispectra_verbose > 1)
        printf (" -> found bispectra folder in run directory.\n");
    }
    /* Otherwise, create it */
    else if (ppr->store_bispectra_to_disk == _TRUE_) {

      if (pbi->bispectra_verbose > 1)
        printf (" -> bispectra folder not found in run directory, will create it.\n");

      class_test (mkdir (pbi->bispectra_dir, 0777)!=0,
        errmsg,
        "could not create directory '%s', maybe it already exists?", pbi->bispectra_dir);

      ppr->load_bispectra_from_disk = _FALSE_;
    }
  }

  /* Create/open the status file. The 'a+' mode means that if the file does not exist it will be created,
  but if it exist it won't be erased (append mode) */
  if (ppr->store_bispectra_to_disk == _TRUE_) {
    // sprintf(pbi->bispectra_status_path, "%s/bispectra_status_file.txt", ppr->data_dir);
    // class_open(pbi->bispectra_status_file, pbi->bispectra_status_path, "a+", errmsg);
  }

  class_test ((ppr->store_bispectra_to_disk == _TRUE_) && (ppr->load_bispectra_from_disk == _TRUE_),
    errmsg,
    "cannot load and save bispectra at the same time!");


  // =======================================================================================
  // =                               Even or odd l-grid?                                   =
  // =======================================================================================

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
  will give an exact result.

  If you are interested only in the Fisher matrix (with mesh interpolation) and only for
  analytical bispectra, it is not necessary to set 'ppr->compute_only_even_ls = _TRUE_'.
  In fact, in this case the Fisher module does not need the precomputed bispectrum in
  pbi->bispectra, as it obtains the bispectrum from scratch using the closed formulas for
  the analytical bispectrum, which are quick to evaluate. However, the computed bispectrum
  in pbi->bispectra will still be subject to the pitfall described above (i.e. it
  might have most zero entries if the 1D l-grid is chosen such that l1+l2+l3 is
  mostly odd). Therefore, you might consider setting 'ppr->compute_only_even_ls = _TRUE_' by
  hand in the functions that rely on pbi->bispectra, such as 'print_bispectra' */
  // if ((ppr->l_linstep!=1)
  // && (pfi->bispectra_interpolation != mesh_interpolation_2D)
  // && (pfi->bispectra_interpolation != mesh_interpolation_3D)
  // && (ppt->has_bi_cmb_polarization==_TRUE_)
  // && ((pbi->has_quadratic_correction==_TRUE_)
  // || (pbi->has_cmb_lensing ==_TRUE_)
  // || (pbi->has_cmb_lensing_squeezed ==_TRUE_)
  // || (pbi->has_cmb_lensing_kernel ==_TRUE_))) {
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
  // }


  // *** END OF MY MODIFICATIONS ***

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

int input_default_params(
			 struct background *pba,
			 struct thermo *pth,
			 struct perturbs *ppt,
			 struct bessels * pbs,
			 struct transfers *ptr,
			 struct primordial *ppm,
			 struct spectra *psp,
			 struct bispectra *pbi,
       struct fisher *pfi,
			 struct nonlinear * pnl,
			 struct lensing *ple,
			 struct output *pop
			 ) {

  double sigma_B; /**< Stefan-Boltzmann constant in W/m^2/K^4 = Kg/K^4/s^3 */

  sigma_B = 2. * pow(_PI_,5) * pow(_k_B_,4) / 15. / pow(_h_P_,3) / pow(_c_,2);




  /** - background structure */

  pba->h = 0.704;
  pba->H0 = pba->h * 1.e5 / _c_;
  pba->T_cmb = 2.726;
  pba->Omega0_g = (4.*sigma_B/_c_*pow(pba->T_cmb,4.)) / (3.*_c_*_c_*1.e10*pba->h*pba->h/_Mpc_over_m_/_Mpc_over_m_/8./_PI_/_G_);
  pba->Omega0_ur = 3.04*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  pba->Omega0_b = 0.02253/0.704/0.704;
  pba->Omega0_cdm = 0.1122/0.704/0.704;
  pba->N_ncdm = 0;
  pba->Omega0_ncdm_tot = 0.;
  pba->ksi_ncdm_default = 0.;
  pba->ksi_ncdm = NULL;
  pba->T_ncdm_default = pow(4.0/11.0,1.0/3.0);
  pba->T_ncdm = NULL;
  pba->deg_ncdm_default = 1.;
  pba->deg_ncdm = NULL;
  pba->ncdm_psd_parameters = NULL;
  pba->ncdm_psd_files = NULL;

  pba->Omega0_k = 0.;
  pba->Omega0_lambda = 1.+pba->Omega0_k-pba->Omega0_g-pba->Omega0_ur-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_ncdm_tot;
  pba->Omega0_fld = 0.;
  pba->a_today = 1.;
  pba->w0_fld=-1.;
  pba->wa_fld=0.;
  pba->cs2_fld=1.;

  /** - thermodynamics structure */

  pth->YHe=_BBN_;
  pth->recombination=recfast;
  pth->reio_parametrization=reio_camb;
  pth->reio_z_or_tau=reio_z;
  pth->z_reio=10.3;
  pth->tau_reio=0.085;
  pth->reionization_exponent=1.5;
  pth->reionization_width=1.5;
  pth->helium_fullreio_redshift=3.5;
  pth->helium_fullreio_width=0.5;

  pth->binned_reio_num=0;
  pth->binned_reio_z=NULL;
  pth->binned_reio_xe=NULL;
  pth->binned_reio_step_sharpness = 0.3;

  pth->compute_cb2_derivatives=_FALSE_;

  // *** MY MODIFICATIONS ***
  pth->has_rayleigh_scattering = _FALSE_;
  // *** END OF MY MODIFICATIONS ***

  /** - perturbation structure */
  ppt->has_cl_cmb_temperature = _FALSE_;
  ppt->has_cl_cmb_polarization = _FALSE_;
  // *** MY MODIFICATIONS ***
  ppt->has_cl_cmb_rayleigh = _FALSE_;
  ppt->has_bi_cmb_temperature = _FALSE_;
  ppt->has_bi_cmb_polarization = _FALSE_;
  ppt->has_bi_cmb_rayleigh = _FALSE_;
  // *** END OF MY MODIFICATIONS ***
  ppt->has_cl_cmb_lensing_potential = _FALSE_;
  ppt->has_cl_density = _FALSE_;
  ppt->has_pk_matter = _FALSE_;
  ppt->has_matter_transfers = _FALSE_;

  ppt->has_ad=_TRUE_;
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;

  ppt->has_scalars=_TRUE_;
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;

  ppt->l_scalar_max=2500;
  ppt->l_tensor_max=500;
  ppt->k_scalar_kmax_for_pk=0.1;

  // *** MY MODIFICATIONS ***
  pth->rayleigh_frequency = 143;
  // *** END OF MY MODIFICATIONS ***

  ppt->gauge=0;



  /** - bessels structure */

  pbs->l_max = MAX(ppt->l_scalar_max,ppt->l_tensor_max);
  pbs->bessel_always_recompute = _TRUE_;

  /** - primordial structure */

  ppm->primordial_spec_type = analytic_Pk;
  ppm->k_pivot = 0.002;
  ppm->A_s = 2.42e-9;
  ppm->n_s = 0.967;
  ppm->alpha_s = 0.;
  ppm->f_bi = 1.;
  ppm->n_bi = 1.;
  ppm->alpha_bi = 0.;
  ppm->f_cdi = 1.;
  ppm->n_cdi = 1.;
  ppm->alpha_cdi = 0.;
  ppm->f_nid = 1.;
  ppm->n_nid = 1.;
  ppm->alpha_nid = 0.;
  ppm->f_niv = 1.;
  ppm->n_niv = 1.;
  ppm->alpha_niv = 0.;
  ppm->c_ad_bi = 0.;
  ppm->n_ad_bi = 0.;
  ppm->alpha_ad_bi = 0.;
  ppm->c_ad_cdi = 0.;
  ppm->n_ad_cdi = 0.;
  ppm->alpha_ad_cdi = 0.;
  ppm->c_ad_nid = 0.;
  ppm->n_ad_nid = 0.;
  ppm->alpha_ad_nid = 0.;
  ppm->c_ad_niv = 0.;
  ppm->n_ad_niv = 0.;
  ppm->alpha_ad_niv = 0.;
  ppm->c_bi_cdi = 0.;
  ppm->n_bi_cdi = 0.;
  ppm->alpha_bi_cdi = 0.;
  ppm->c_bi_nid = 0.;
  ppm->n_bi_nid = 0.;
  ppm->alpha_bi_nid = 0.;
  ppm->c_bi_niv = 0.;
  ppm->n_bi_niv = 0.;
  ppm->alpha_bi_niv = 0.;
  ppm->c_cdi_nid = 0.;
  ppm->n_cdi_nid = 0.;
  ppm->alpha_cdi_nid = 0.;
  ppm->c_cdi_niv = 0.;
  ppm->n_cdi_niv = 0.;
  ppm->alpha_cdi_niv = 0.;
  ppm->c_nid_niv = 0.;
  ppm->n_nid_niv = 0.;
  ppm->alpha_nid_niv = 0.;
  ppm->r = 1.;
  ppm->n_t = -ppm->r/8.*(2.-ppm->r/8.-ppm->n_s);
  ppm->alpha_t = ppm->r/8.*(ppm->r/8.+ppm->n_s-1.);

  /** - transfer structure */

  ppt->selection_num=1;
  ppt->selection=gaussian;
  ppt->selection_mean[0]=1.;
  ppt->selection_width[0]=0.1;

  /** - output structure */

  pop->z_pk_num = 1;
  pop->z_pk[0] = 0.;
  sprintf(pop->root,"output/");
  pop->write_header = _TRUE_;
  pop->output_format = class_format;

  /** - spectra structure */

  psp->z_max_pk = pop->z_pk[0];

  /** - nonlinear structure */

  /** - lensing structure */

  ple->has_lensed_cls = _FALSE_;

  /** - nonlinear structure */

  pnl->method = nl_none;
  pnl->ic = nl_pt;

  /** - all verbose parameters */

  pba->background_verbose = 0;
  pth->thermodynamics_verbose = 0;
  ppt->perturbations_verbose = 0;
  pbs->bessels_verbose = 0;
  ptr->transfer_verbose = 0;
  ppm->primordial_verbose = 0;
  psp->spectra_verbose = 0;
  pnl->nonlinear_verbose = 0;
  ple->lensing_verbose = 0;
  pop->output_verbose = 0;









  // *** MY MODIFICATIONS ***


  // ===========================================================
  // =                     Perturbs structure                  =
  // ===========================================================

  /* What is in the line-of-sight integral? */
  ppt->has_scattering_in_los = _TRUE_;
  ppt->has_photon_monopole_in_los = _TRUE_;
  ppt->has_metric_in_los = _TRUE_;
  ppt->has_sw = _FALSE_;
  ppt->has_isw = _FALSE_;

  /* Alternative initial conditions at 1st order */
  ppt->has_ad_maberty=_FALSE_;
  ppt->has_zero_ic=_FALSE_;

  ppt->has_cl_cmb_zeta = _FALSE_; /* Compute zeta variable? */
  ppt->recombination_only_zeta = _TRUE_; /* Is zeta evaluated exclusively at recombination? */
  ppt->has_bispectra = _FALSE_; /* Do we need bispectra? */

  // ==========================================================
  // =                 Perturbed recombination                =
  // ==========================================================

  ppt->has_perturbed_recombination = _FALSE_;
  pth->has_perturbed_recombination = _FALSE_;
  pth->compute_xe_derivatives = _FALSE_;
  pth->perturbed_recombination_turnx = 36;


  // ===============================================
  // =                Bessel structure             =
  // ===============================================



  // ================================================
  // =                Spectra structure             =
  // ================================================

  psp->compute_cl_derivative = _FALSE_;


  // =============================================================
  // =                     Bispectra structure                   =
  // =============================================================

  pbi->bispectra_verbose = 0;
  pbi->has_bispectra = _FALSE_;
  pbi->has_local_model = _FALSE_;
  pbi->has_equilateral_model = _FALSE_;
  pbi->has_orthogonal_model = _FALSE_;
  pbi->has_galileon_model = _FALSE_;
  pbi->has_intrinsic_squeezed = _FALSE_;
  pbi->has_intrinsic_squeezed_unlensed = _FALSE_;
  pbi->has_local_squeezed = _FALSE_;
  pbi->has_cosine_shape = _FALSE_;
  pbi->has_cmb_lensing = _FALSE_;
  pbi->has_cmb_lensing_squeezed = _FALSE_;
  pbi->has_cmb_lensing_kernel = _FALSE_;
  pbi->has_quadratic_correction = _FALSE_;
  pbi->include_lensing_effects = _FALSE_;
  pbi->lensed_intrinsic = _FALSE_;

  // ========================================================
  // =                    Fisher structure                  =
  // ========================================================

  pfi->fisher_verbose = 0;
  pfi->always_interpolate_bispectra = _FALSE_;
  pfi->has_fisher = _FALSE_;
  pfi->l_min_estimator = 2;
  pfi->l_max_estimator = 10000000;
  pfi->bispectra_interpolation = mesh_interpolation_2D;
  pfi->f_sky = 1;
  pfi->n_channels = 1;
  pfi->beam[0] = 0;
  pfi->noise_t[0] = 0;
  pfi->noise_e[0] = 0;
  pfi->noise_r[0] = 0;
  pfi->ignore_t = _FALSE_;
  pfi->ignore_e = _FALSE_;
  pfi->ignore_b = _FALSE_;
  pfi->ignore_r = _FALSE_;
  pfi->include_lensing_effects = _FALSE_;
  pfi->squeezed_ratio = 0;


  // ==========================================================
  // =               Specific to second order                 =
  // ==========================================================

  ppt->has_perturbations2 = _FALSE_;
  ppt->has_polarization2  = _FALSE_;

  ppt->has_custom_timesampling_for_quadsources = _FALSE_;
  ppt->custom_tau_ini_quadsources  = 0.1;
  ppt->custom_tau_end_quadsources  = 0;
  ppt->custom_tau_size_quadsources = 2000;
  ppt->custom_tau_mode_quadsources = lin_tau_sampling;

  pbi->has_intrinsic = _FALSE_;


  // *** END OF MY MODIFICATIONS


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

int input_default_precision ( struct precision * ppr ) {

  /** Summary: */

  /**
   * - parameters related to the background
   */

  ppr->a_ini_over_a_today_default = 1.e-14;
  ppr->back_integration_stepsize = 7.e-3;
  ppr->tol_background_integration = 1.e-2;

  ppr->tol_initial_Omega_r = 1.e-4;
  ppr->tol_M_ncdm = 1.e-7;
  ppr->tol_ncdm = 1.e-3;
  ppr->tol_ncdm_bg = 1.e-5;
  ppr->tol_ncdm_initial_w=1.e-3;

  /**
   * - parameters related to the thermodynamics
   */

  /* for bbn */
  sprintf(ppr->sBBN_file,"bbn/sBBN.dat");

  /* for recombination */

  ppr->recfast_z_initial=1.e4;

  ppr->recfast_Nz0=20000;
  ppr->tol_thermo_integration=1.e-2;

  ppr->recfast_Heswitch=6;                 /* from recfast 1.4 */
  ppr->recfast_fudge_He=0.86;              /* from recfast 1.4 */

  ppr->recfast_Hswitch = _TRUE_;           /* from recfast 1.5 */
  ppr->recfast_fudge_H = 1.14;             /* from recfast 1.4 */
  ppr->recfast_delta_fudge_H = -0.035;     /* from recfast 1.5 */
  ppr->recfast_AGauss1 = -0.14;            /* from recfast 1.5 */
  ppr->recfast_AGauss2 =  0.05;            /* from recfast 1.5 */
  ppr->recfast_zGauss1 =  7.28;            /* from recfast 1.5 */
  ppr->recfast_zGauss2 =  6.75;            /* from recfast 1.5 */
  ppr->recfast_wGauss1 =  0.18;            /* from recfast 1.5 */
  ppr->recfast_wGauss2 =  0.33;            /* from recfast 1.5 */

  ppr->recfast_z_He_1 = 8000.;             /* from recfast 1.4 */
  ppr->recfast_delta_z_He_1 = 50.;         /* found to be OK on 3.09.10 */
  ppr->recfast_z_He_2 = 5000.;             /* from recfast 1.4 */
  ppr->recfast_delta_z_He_2 = 100.;        /* found to be OK on 3.09.10 */
  ppr->recfast_z_He_3 = 3500.;             /* from recfast 1.4 */
  ppr->recfast_delta_z_He_3 = 50.;         /* found to be OK on 3.09.10 */
  ppr->recfast_x_He0_trigger = 0.995;      /* raised from 0.99 to 0.995 for smoother Helium */
  ppr->recfast_x_He0_trigger2 = 0.995;     /* raised from 0.985 to same as previous one for smoother Helium */
  ppr->recfast_x_He0_trigger_delta = 0.05; /* found to be OK on 3.09.10 */
  ppr->recfast_x_H0_trigger = 0.995;       /* raised from 0.99 to 0.995 for smoother Hydrogen */
  ppr->recfast_x_H0_trigger2 = 0.995;      /* raised from 0.98 to same as previous one for smoother Hydrogen */
  ppr->recfast_x_H0_trigger_delta = 0.05;  /* found to be OK on 3.09.10 */

  ppr->recfast_H_frac=1.e-3;               /* from recfast 1.4 */

  sprintf(ppr->hyrec_Alpha_inf_file,"hyrec/Alpha_inf.dat");
  sprintf(ppr->hyrec_R_inf_file,"hyrec/R_inf.dat");
  sprintf(ppr->hyrec_two_photon_tables_file,"hyrec/two_photon_tables.dat");

  /* for reionization */

  ppr->reionization_z_start_max = 50.;
  ppr->reionization_sampling=1.e-2;
  ppr->reionization_optical_depth_tol=1.e-2;
  ppr->reionization_start_factor=8.;

  /* general */

  ppr->thermo_rate_smoothing_radius=50;

  /**
   * - parameters related to the perturbations
   */

  ppr->evolver = ndf15;
  ppr->pk_definition = delta_m_squared;

  ppr->k_scalar_min_tau0=1.;
  ppr->k_scalar_max_tau0_over_l_max=2.;
  ppr->k_scalar_step_sub=0.05;
  ppr->k_scalar_step_super=0.0025;
  ppr->k_scalar_step_transition=0.2;

  ppr->k_scalar_k_per_decade_for_pk=10.;
  ppr->k_scalar_k_per_decade_for_bao=70.;
  ppr->k_scalar_bao_center=3.;
  ppr->k_scalar_bao_width=4.;

  ppr->k_tensor_min_tau0=1.4;
  ppr->k_tensor_max_tau0_over_l_max = 2.;
  ppr->k_tensor_step_sub=0.1;
  ppr->k_tensor_step_super=0.0025;
  ppr->k_tensor_step_transition=0.2;

  ppr->start_small_k_at_tau_c_over_tau_h = 0.0015;  /* decrease to start earlier in time */
  ppr->start_large_k_at_tau_h_over_tau_k = 0.07;  /* decrease to start earlier in time */
  ppr->tight_coupling_trigger_tau_c_over_tau_h=0.015; /* decrease to switch off earlier in time */
  ppr->tight_coupling_trigger_tau_c_over_tau_k=0.01; /* decrease to switch off earlier in time */
  ppr->start_sources_at_tau_c_over_tau_h = 0.008; /* decrease to start earlier in time */
  ppr->tight_coupling_approximation=(int)compromise_CLASS;

  ppr->l_max_g=10;
  ppr->l_max_pol_g=8;
  ppr->l_max_ur=12;
  ppr->l_max_ncdm=12;
  ppr->l_max_g_ten=5;
  ppr->l_max_pol_g_ten=5;

  ppr->curvature_ini=1.; /* initial curvature; used to fix adiabatic initial conditions; must remain fixed to one as long as the primordial adiabatic spectrum stands for the curvature power spectrum */

  ppr->entropy_ini=1.;   /* initial entropy; used to fix isocurvature initial conditions; must remain fixed to one as long as the primordial isocurvature spectrum stands for an entropy power spectrum */
  ppr->gw_ini=0.25; /* to match normalization convention for GW in most of literature and ensure standard definition of r */

  ppr->perturb_integration_stepsize=0.5;

  ppr->tol_tau_approx=1.e-5;
  ppr->tol_perturb_integration=1.e-4;
  ppr->perturb_sampling_stepsize=0.08;

  ppr->selection_cut_at_sigma=5.;
  ppr->l_switch_limber_for_cl_density_over_z=40.;

  ppr->radiation_streaming_approximation = rsa_MD_with_reio;
  ppr->radiation_streaming_trigger_tau_over_tau_k = 45.;
  ppr->radiation_streaming_trigger_tau_c_over_tau = 5.;

  ppr->ur_fluid_approximation = ufa_CLASS;
  ppr->ur_fluid_trigger_tau_over_tau_k = 15.;

  ppr->ncdm_fluid_approximation = ncdmfa_CLASS;
  ppr->ncdm_fluid_trigger_tau_over_tau_k = 16.;




  /**
   * - parameter related to the Bessel functions
   */

  ppr->l_logstep=1.15;
  ppr->l_linstep=40;

  ppr->bessel_x_step=0.5;
  ppr->bessel_j_cut=1.e-5;
  ppr->bessel_tol_x_min =1.e-4;
  sprintf(ppr->bessel_file_name,"bessels.dat");

  /**
   * - parameter related to the primordial spectra
   */

  ppr->k_per_decade_primordial = 10.;

  /**
   * - parameter related to the transfer functions
   */

  ppr->k_step_trans_scalars=0.004;
  ppr->k_step_trans_tensors=0.004;
  ppr->transfer_cut=tc_osc;
  ppr->transfer_cut_threshold_osc=0.007; /* 03.12.10 for chi2plT0.01 */
  ppr->transfer_cut_threshold_cl=1.e-8; /* 14.12.10 for chi2plT0.01 */

  ppr->l_switch_limber=10.;


  /**
   * - parameters related to trg module
   */

  ppr->halofit_dz=0.1;
  ppr->halofit_min_k_nonlinear=0.0035;
  ppr->halofit_sigma_precision=0.05;
  ppr->double_escape=2;
  ppr->z_ini = 35.;
  ppr->eta_size = 101;
  ppr->k_L = 1.e-3;
  ppr->k_min = 1.e-4;
  ppr->logstepx_min = 1.04;
  ppr->logstepk1 = 1.11;
  ppr->logstepk2 = 0.09;
  ppr->logstepk3 = 300.;
  ppr->logstepk4 = 0.01;
  ppr->logstepk5 = 1.02;
  ppr->logstepk6 = 0.;
  ppr->logstepk7 = 0.;
  ppr->logstepk8 = 0.;
  ppr->k_growth_factor = 0.1;
  ppr->k_scalar_max_for_pk_nl = 1000.;

  /**
   * - parameter related to lensing
   */

  ppr->accurate_lensing=_FALSE_;
  ppr->num_mu_minus_lmax=70;
  ppr->delta_l_max=250;

  /**
   * - automatic estimate of machine precision
   */

  get_machine_precision(&(ppr->smallest_allowed_variation));

  class_test(ppr->smallest_allowed_variation < 0,
	     ppr->error_message,
	     "smallest_allowed_variation = %e < 0",ppr->smallest_allowed_variation);

  ppr->tol_gauss_legendre = ppr->smallest_allowed_variation;



  // *** MY MODIFICATIONS ***

  /* k-sampling */
  ppr->k_scalar_logstep_super=1.2;

  /* l-sampling */
  ppr->compute_only_even_ls = _FALSE_;
  ppr->compute_only_odd_ls = _FALSE_;
  ppr->extend_lensed_cls = _FALSE_;

  /* Quadsources time-sampling */
  ppr->perturb_sampling_stepsize_quadsources = 0.03;

  /* Interpolation and integration */
  ppr->quadsources_time_interpolation = cubic_interpolation;
  ppr->bessels_interpolation = linear_interpolation;
  ppr->transfers_k1_interpolation = linear_interpolation;
  ppr->transfers_k2_interpolation = linear_interpolation;
  ppr->bispectra_k3_extrapolation = flat_k3_extrapolation;
  ppr->extra_k3_oscillations_left = 50;
  ppr->extra_k3_oscillations_right = 50;

  /* Storage of intermediate results */
  ppr->store_run == _FALSE_;
  ppr->append_date_to_run = _FALSE_;
  ppr->store_bispectra_to_disk = _FALSE_;
  ppr->load_bispectra_from_disk = _FALSE_;

  /* Do not specify default values for ppr->run_dir and for ppr->load_run, as
  these variables are not set in input_init but in the parent function,
  input_init_from_arguments */

  /* By default, assume that the data folders (sources, transfers, bispectra) are in the
  run directory */
  strcpy (ppr->data_dir, ppr->run_dir);

  /* Reionisation flag */
  ppr->has_reionization = _FALSE_;

  // *** END OF MY MODIFICATIONS





  return _SUCCESS_;

}







int class_version(
		  char * version
		  ) {

  sprintf(version,"%s",_VERSION_);
  return _SUCCESS_;
}








/**
 * Computes automatically the machine precision.
 *
 * @param smallest_allowed_variation a pointer to the smallest allowed variation
 *
 * Returns the smallest
 * allowed variation (minimum epsilon * _TOLVAR_)
 */

int get_machine_precision(double * smallest_allowed_variation) {
  double one, meps, sum;

  one = 1.0;
  meps = 1.0;
  do {
    meps /= 2.0;
    sum = one + meps;
  } while (sum != one);
  meps *= 2.0;

  *smallest_allowed_variation = meps * _TOLVAR_;

  return _SUCCESS_;

}
