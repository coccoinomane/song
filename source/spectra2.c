/** @file spectra2.c Compute spectra of second order pertrubations
 *
 * Christian Fidler, 11.05.2015    
 * 
 *
 */

#include "spectra2.h"

#ifdef sources
#undef sources
#endif

#define sources(INDEX_TAU,INDEX_K_TRIANGULAR) \
  ppt2->sources[index_tp2]\
               [index_k1]\
               [index_k2]\
               [(INDEX_TAU)*k_pt_size + (INDEX_K_TRIANGULAR)]




int spectra2_init(
			struct primordial * ppm,
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
      struct spectra2 * psp2
      )
{

  // =================================================================================
  // =                              Preliminary checks                               =
  // =================================================================================

	if ((ppt2->has_perturbations2 == _FALSE_)
     || ((ppt2->has_cls == _FALSE_) && (ppt2->has_pk_matter == _FALSE_)) ){

    if (psp2->spectra2_verbose > 0) 
      printf("Second-order spectra module skipped.\n");
      
    return _SUCCESS_;
  }

	
	
	// ==================================================================================
  // =                               k-sampling                                       =
  // ==================================================================================

	// We use the same sampling as used in ppt2, for k3 we use the sampling of tr2 for angular power spectra
	
	psp2->k_size = ppt2->k_size;
  psp2->k = ppt2->k;
  
  class_calloc(psp2->k3_grid, ptr2->k3_size_max, sizeof(double), psp2->error_message);
 
  // ==================================================================================
  // =                               Allocate workspaces                              =
  // ==================================================================================


 	
  class_alloc (psp2 -> spectra, ppt2->tp2_size*sizeof(double *),psp2->error_message);
  
  for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
         
       class_alloc(
          psp2->spectra[index_tp],
          psp2->k_size*ppt2->tau_size*sizeof(double),
          psp2->error_message);
          for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) {
						for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
							psp2->spectra[index_tp][index_k3*ppt2->tau_size + index_tau] = 0.;
						} 
					}
  }
  
 
  if (psp2->spectra2_verbose > 1) 
  	printf("allocated workspace \n");
  
  
  // ==================================================================================
  // =                               Angular power spectra                            =
  // ==================================================================================

		double step_k1, step_k2, step_k3;
		double ClTT,ClTE,ClEE,ClBB;
		int index_tt2_T,index_tt2_E,index_tt2_B;
		
		
		if (ppt2->has_cls == _TRUE_) {
		
			if ((ppt2->k3_sampling == sym_k3_sampling)) {
				printf("Angular power spectra can only be computed in standart k3 sampling. Change sampling strategy\n");
			}
  
   	
			
   
		
      // -------------------------------------------------------------------------------------
      // -                              Sum over M, cycle L                              -
      // -------------------------------------------------------------------------------------
    for (int index_l=0; index_l<ptr2->l_size; ++index_l) {
  		
  		int l = ptr2->l[index_l];
  		
  		//if parallelisation os wanted it probably should be done here
  		
  		if (ppt2->has_cmb_temperature == _TRUE_) {
  			ClTT = 0.;
  		}
  		if (ppt2-> has_cmb_polarization_e == _TRUE_){
  			ClEE = 0.;
  		}
  		if (ppt2-> has_cmb_polarization_b == _TRUE_){
  			ClBB = 0.;
  		}
  		if ((ppt2-> has_cmb_temperature == _TRUE_) && (ppt2->has_cmb_polarization_e == _TRUE_)) {
  			ClTE = 0.;
  		}
  		
  		//loop over M
  		
      for (int index_M=0; index_M < ppr2->m_size; ++index_M) {
				if (psp2->spectra2_verbose > 2)
					printf("computing angular spectra for l = %d m = %d\n",ptr2->l[index_l],ptr2->m[index_M]);
				int m = ptr2->m[index_M];
				
				//load transfer functions
				
				if (ppt2-> has_cmb_temperature == _TRUE_) {
  				index_tt2_T = ptr2->index_tt2_T + lm_cls(index_l, index_M);
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
     				class_call (transfer2_load_transfers_from_disk (
            	      ppt2,
                    ptr2,
                    index_tt2_T),
        			ptr2->error_message,
        			psp2->error_message);
    			}
    		}
    		
    		if (ppt2->has_cmb_polarization_e == _TRUE_) {
  				index_tt2_E = ptr2->index_tt2_E + lm_cls(index_l, index_M);
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
     				class_call (transfer2_load_transfers_from_disk (
            	      ppt2,
                    ptr2,
                    index_tt2_E),
        			ptr2->error_message,
        			psp2->error_message);
    			}
    		}
    		
    		if (ppt2->has_cmb_polarization_b == _TRUE_) {
  				index_tt2_B = ptr2->index_tt2_B + lm_cls(index_l, index_M);
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
     				class_call (transfer2_load_transfers_from_disk (
            	      ppt2,
                    ptr2,
                    index_tt2_B),
        			ptr2->error_message,
        			psp2->error_message);
    			}
    		}
  			
				
  // ==================================================================================
  // =                               integrations over k1,k2,k3                       =
  // ==================================================================================
				
	
				for (int index_k1 = 0; index_k1 < psp2->k_size; ++index_k1) {
					
					// find stepsize
				
					if (index_k1 == 0) step_k1 = (psp2->k[1] - psp2->k[0])/2.;
        	else if (index_k1 == psp2->k_size -1) step_k1 = (psp2->k[psp2->k_size-1] - psp2->k[psp2->k_size -2])/2.;
        	else step_k1 = (psp2->k[index_k1+1] - psp2->k[index_k1 -1])/2.;

					double k1 = psp2->k[index_k1];
					
					//primordial spectrum k1
					
					double spectra_k1;
    			class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k1,
                  &spectra_k1),
      			ppm->error_message,
      			psp2->error_message);
    			spectra_k1 = 2*_PI_*_PI_/(k1*k1*k1) * spectra_k1;
    			
    			// Loop over k2
  
					for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
					
						/* Print some info */
      			if (psp2->spectra2_verbose > 3)
        			printf(" -> computing integral for (k1,k2) = (%.3g,%.3g)\n", psp2->k[index_k1], psp2->k[index_k2]);
        			
        		// find stepsize 
        			
						if (index_k1 == 0) step_k2 = 0.;
						else if (index_k2 == 0) step_k2 = (psp2->k[1] - psp2->k[0])/2.;
        		else if (index_k2 == index_k1) step_k2 = (psp2->k[index_k1] - psp2->k[index_k1-1])/2.;
        		else step_k2 = (psp2->k[index_k2+1] - psp2->k[index_k2 -1])/2.;
        
        		double k2 = psp2->k[index_k2];
        		
        		//primordial spectrum k2
        		
						double spectra_k2;
    				class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k2,
                  &spectra_k2),
      				ppm->error_message,
      				psp2->error_message);
	    			spectra_k2 = 2*_PI_*_PI_/(k2*k2*k2) * spectra_k2;
  
						int dump;
						class_call(transfer2_get_k3_list (
                                 ppr,
                                 ppr2,
                                 ppt2,
                                 pbs,
                                 pbs2,
                                 ptr2,
                                 index_k1,
                                 index_k2,
                                 psp2->k3_grid,  
                                 &dump   
                                 ),
        	    ptr2->error_message,
          	  psp2->error_message);
          	  
          	  
          	int k3_size = ptr2->k_size_k1k2[index_k1][index_k2];
          
          	class_test(k3_size < 2,
            		psp2->error_message,
            		"integration grid has less than two elements, cannot use trapezoidal integration");
        
						
         		/* Define the pointer to the second-order transfer function as a function of k3.
          	Note that this transfer function has already been rescaled according to eq. 6.26
          	of http://arxiv.org/abs/1405.2280 in the perturbations.c module.  */
          	
          	double * transferT = ptr2->transfer[index_tt2_T][index_k1][index_k2];
          	double * transferE = ptr2->transfer[index_tt2_E][index_k1][index_k2];
          	double * transferB = ptr2->transfer[index_tt2_B][index_k1][index_k2];
          	
          	int triangular_first = 0;
          	int triangular_last = 0;	
          	
          	//final integration over k3
          	
          	for (int index_k3 = 0; index_k3 < k3_size; ++index_k3) {
     					double k3 = psp2->k3_grid[index_k3];
     					
     					//find stepsize including the triangular condition edge
     					
     					if (index_k3 == 0) step_k3 = (psp2->k3_grid[1] - psp2->k3_grid[0])/2.;
        			else if (index_k3 == k3_size -1) step_k3 = (psp2->k3_grid[k3_size-1] - psp2->k3_grid[k3_size -2])/2.;
        			else step_k3 = (psp2->k3_grid[index_k3+1] - psp2->k3_grid[index_k3 -1])/2.;
     					
     					// triangular condition
     					if ( k3 < k1-k2 ||  k3 > k1+k2) { // out of triangular region
     						step_k3 = 0.;
     					}
     					if ( k3 >= k1-k2 && triangular_first == 0 && index_k3 < k3_size-1 ) { //first point
     						triangular_first = 1;
     						step_k3 = (psp2->k3_grid[index_k3+1]-k3)/2. + k3 - (k1-k2);
     						
     					}
     					if ( psp2->k3_grid[index_k3+1] > k1+k2 && triangular_last == 0 && index_k3 > 0) {
     						//last point, note that this may not be found if the last point is bigger than kmax, 
     						// however in that case no special treatment for the last point is needed. 
     						triangular_last = 1;
     						step_k3 = (k3 - psp2->k3_grid[index_k3-1])/2. + (k1+k2) - k3;
     						
     					}
     					
     					
     					// This counters the rescaling performed for the bispectrum 
     					double scale = pow(sqrt( 1. - pow((k3*k3 + k1*k1 - k2*k2)/2./k3/k1,2)),m);
     					
     					
     					// add up contributions
     					
     					double total_factor = 
     							// symmetry factor (doing only half plane in k1 k2)
					 				2.*
									// integration weight	
					 				k1*k2*k3*
									// integration stepsize
					 				step_k1*step_k2*step_k3*
					 				// if the sources were rescaled we have to invert this for the cl spectrum rescaling
					 				//the step_k3 condition makes sure that this is in triangular as the computation of angles fails out of trinagular region
					 		 		(ppt2->rescale_quadsources == _TRUE_ && step_k3 != 0. ? pow(scale,2) : 1. )*
									// spectra factors (2l+1) already in transfer def (how does sum over m work? does this need a factor 2 for m neq 0? yes :))
					 				2. /_PI_ /2./_PI_/2./_PI_/2./_PI_ *
					 				// definition of second order perturbation theory
					 				// already in transfer defenition 1./2./2.*
					 				// factor 4 of Delta also already in transfer def
									// primordial spectra
					 	  		spectra_k1*spectra_k2;
     					
     					if (ppt2-> has_cmb_temperature == _TRUE_) 
  							ClTT += total_factor*transferT[index_k3]*transferT[index_k3] ;
  						if (ppt2-> has_cmb_polarization_e == _TRUE_) 
  							ClEE += total_factor*transferE[index_k3]*transferE[index_k3] ;
  						if (ppt2-> has_cmb_polarization_b == _TRUE_) 
  							ClBB += total_factor*transferB[index_k3]*transferB[index_k3] ;
  						if ((ppt2-> has_cmb_temperature == _TRUE_) && (ppt2->has_cmb_polarization_e == _TRUE_)) 
  							ClTE += total_factor*transferT[index_k3]*transferE[index_k3] ;
     					     					
     					
     					
						} // loop k3
						
					
					} // loop k2 
					
					
				} // loop k1
				
				if (ppt2-> has_cmb_temperature == _TRUE_) {
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
  	    		class_call (transfer2_free_type_level (
                    ppt2,
                    ptr2,
                    index_tt2_T),
        		ptr2->error_message,
        		psp2->error_message);	
    	 	 }
				}
				
				
				if (ppt2-> has_cmb_polarization_e == _TRUE_) {
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
  	    		class_call (transfer2_free_type_level (
                    ppt2,
                    ptr2,
                    index_tt2_E),
        		ptr2->error_message,
        		psp2->error_message);	
    	 	 }
				}
				
				
				if (ppt2-> has_cmb_polarization_b == _TRUE_) {
					if ((ppr2->load_transfers_from_disk == _TRUE_) || (ppr2->store_transfers_to_disk == _TRUE_)) {
  	    		class_call (transfer2_free_type_level (
                    ppt2,
                    ptr2,
                    index_tt2_B),
        		ptr2->error_message,
        		psp2->error_message);	
    	 	 }
				}
	    } // Sum over M
	    printf("%d %g %g %g %g \n",l,
	    	(ppt2-> has_cmb_temperature == _TRUE_?ClTT:0.),
	    	((ppt2-> has_cmb_polarization_e == _TRUE_) && (ppt2-> has_cmb_temperature == _TRUE_) ?ClTE:0.),
	    	(ppt2-> has_cmb_polarization_e == _TRUE_?ClEE:0.),
	   	 	(ppt2-> has_cmb_polarization_b == _TRUE_?ClBB:0.));
  	} // Loop over l
		free(psp2->k3_grid);
		
	}	
		
  // ==================================================================================
  // =                               Fourier Power Spectra		                        =
  // ==================================================================================

  if ((ppt2->has_pk_matter == _TRUE_)){
  
 	class_call(spectra2_get_k3_size (
      ppr,
      ppr2,
      ppt2,
      psp2
      ),
    psp2->error_message,
    psp2->error_message);
  


	/* Shortcut to the file where we shall print the transfer functions. */
  FILE * file_sp = psp2->spectra_file;
  
  /* Choose how label & values should be formatted */
  char format_label[64] = "%13s(%02d) ";
  char format_value[64] = "%+17e ";      
  char buffer[64];
  
  // conformal time
  int index_sp = 1;
	
	// Print time 
	/*
	fprintf(file_sp, format_label, "tau",index_sp);
	index_sp++;
  
  fprintf(file_sp, format_label, "a", index_sp);
  index_sp++;
  
  fprintf(file_sp, format_label, "Y",index_sp);		
  index_sp++;
  	
  for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) { 
  		sprintf(buffer, "magnet k=%f", psp2->k[index_k3] );
 			fprintf(file_sp, format_label, buffer,index_sp );
 			index_sp++;		
 		}
  
  fprintf(file_sp, "\n");
  */
  
  //Print k 
  
  fprintf(file_sp, format_label, "k",index_sp);
	index_sp++;

	for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) { 
  		sprintf(buffer, "magnet tau=%f", ppt2->tau_sampling[index_tau] );
 			fprintf(file_sp, format_label, buffer,index_sp );
 			index_sp++;		
 		} 
 				
 		fprintf(file_sp, "\n");
	
	
  /* Four loops over k1, k2, transfer type follow */
  
 
  if ((ppt2->k3_sampling == sym_k3_sampling) ) {
  	
  	class_call(
  		spectra2_integrate_fourier_sym(
  				ppr,
      		ppr2,
      		ppt,
      		ppm,
      		ppt2,
      		psp2),
      	ptr2->error_message,
        psp2->error_message);
          	  
      	
	} else {
		class_call(
  		spectra2_integrate_fourier(
  				ppr,
      		ppr2,
      		ppt,
      		ppm,
      		ppt2,
      		psp2),
      	ptr2->error_message,
        psp2->error_message);
	
	}

	
  
 

//	for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) {     			
 // 	fprintf(file_sp, format_value, psp2->k[index_k3]);
 //		fprintf(file_sp, format_value, psp2->spectra[ppt2->index_tp2_M + lm(1,1)][index_k3*ppt2->tau_size + 300] );	
 // 	fprintf(file_sp, "\n");
//	}
	int last_index = 0;

	double *pvecback;
  class_alloc (pvecback, pba->bg_size*sizeof(double), psp2->error_message);

 // print time
 /*
	for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) { 
	
	class_call (background_at_tau(
                pba,
                ppt2->tau_sampling[index_tau],
                pba->long_info, 
                pba->inter_normal, 
                &last_index,
                pvecback),
    pba->error_message,
    psp2->error_message);
  
 	 	double a = pvecback[pba->index_bg_a];
 	 	double H = pvecback[pba->index_bg_H];
 	 	double Hc = pvecback[pba->index_bg_H]*a;
 	 	double Y = a/pba->a_eq;  
	
	    			
  	fprintf(file_sp, format_value, ppt2->tau_sampling[index_tau]);
  	
  	fprintf(file_sp, format_value, a);
  	
  	fprintf(file_sp, format_value, Y);
  	
  	for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) { 
  		
 			fprintf(file_sp, format_value, psp2->spectra[ppt2->index_tp2_M + lm(0,0)][index_k3*ppt2->tau_size + index_tau] );	
 			
 		}
  	fprintf(file_sp, "\n");
	}
  */
  // Print k 
  
  for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) { 

	    			
  	fprintf(file_sp, format_value, psp2->k[index_k3]);
  	
  	for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) { 
  		
 			fprintf(file_sp, format_value, psp2->spectra[ppt2->index_tp2_M + lm(0,0)][index_k3*ppt2->tau_size + index_tau] );	
 			
 		}
  	fprintf(file_sp, "\n");
	}
  
  }
  
  
 	printf("keq = %f \n",pba->k_eq);

  fclose(psp2->spectra_file);
  return _SUCCESS_;
}



/**
 * This function interpolates sources S(k1, k2, k, tau) at the needed
 * values of k and tau.  Note that the time sampling for the integration
 * grid could be computed once and for all, as the k grid is fixed (the
 * nodes do not depend on k1 and k2).
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ppt2                  Input : pointer to 2nd-order perturbation structure
 * @param pbs                   Input : pointer to Bessel structure
 * @param psp2                  Input : pointer to 2nd-order spectra structure
 * @return the error status
 */

int spectra2_interpolate_sources_in_k(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2,
      int index_k1,
      int index_k2,
      int index_tp2,
      double * k_grid,
      double * sources_k_spline,
      double * interpolated_sources_in_k
      )
{


  /* Shortcuts */
  int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  double * k_pt = ppt2->k3[index_k1][index_k2];
  int k_sp_size = psp2->k_size;
  double * k_sp = k_grid;
  
   
  /* Cycle index */
  int index_tau;
  
  
       
	
   if (ppr2->sources_k3_interpolation == cubic_interpolation) {

    class_call (array_spline_table_columns (
                  ppt2->k3[index_k1][index_k2],
                  ppt2->k3_size[index_k1][index_k2],
                  ppt2->sources[index_tp2][index_k1][index_k2],
                  ppt2->tau_size,
                  sources_k_spline,
                  _SPLINE_EST_DERIV_,
                  psp2->error_message),
         psp2->error_message,
         psp2->error_message);
  }
 
  // =======================================================
  // =                    Interpolation                    =
  // =======================================================

  /* Limits where for which we shall interpolate the sources */
  int physical_size = psp2->k_physical_size_k1k2[index_k1][index_k2];
  int first_physical_index = psp2->k_physical_start_k1k2[index_k1][index_k2];
  int last_physical_index = first_physical_index + physical_size - 1;

  /* Interpolate at each k value using the usual spline interpolation algorithm */
  int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
  
  int index_k_sp;

    
  for (index_k_sp = first_physical_index; index_k_sp <= last_physical_index; ++index_k_sp) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_sp[index_k_sp])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., psp2->error_message, "stop to avoid division by zero");
    
    double b = (k_sp[index_k_sp] - k_pt[index_k])/h;
    double a = 1.-b;
      
    /* Interpolate for each value of conformal time */
    if (ppr2->sources_k3_interpolation == linear_interpolation) {
      for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1);
    }
    else if (ppr2->sources_k3_interpolation == cubic_interpolation) {
      for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
        interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
          a * sources(index_tau,index_k) + b * sources(index_tau,index_k+1)
          + ((a*a*a-a) * sources_k_spline[index_tau*k_pt_size + index_k]
          +(b*b*b-b) * sources_k_spline[index_tau*k_pt_size + index_k+1])*h*h/6.0;
    } 
  } // end of for (index_k_tr)

//Extrapolation either set to zero or flat

	for (index_k_sp = 0; index_k_sp < first_physical_index; ++index_k_sp) {
    for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++) {
    	interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
    		interpolated_sources_in_k[first_physical_index*ppt2->tau_size + index_tau];
    }
  } // end of for (index_k_tr) 

	for (index_k_sp = last_physical_index+1; index_k_sp < k_sp_size; ++index_k_sp) {
    for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++) {
    	interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
    		interpolated_sources_in_k[last_physical_index*ppt2->tau_size + index_tau];
    }
  } // end of for (index_k_tr) 
 
  return _SUCCESS_;
  
} // end of spectra_interpolate_sources_in_k


int spectra2_get_k3_size (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2
      )
{

	 /* Allocate k1 level */
  int k1_size = psp2->k_size;
 
  class_alloc(psp2->k_true_physical_start_k1k2, k1_size*sizeof(int *), psp2->error_message);
  class_alloc(psp2->k_physical_start_k1k2, k1_size*sizeof(int *), psp2->error_message);
  class_alloc(psp2->k_physical_size_k1k2, k1_size*sizeof(int *), psp2->error_message);

	int index_k1, index_k2;

  for(index_k1=0; index_k1<k1_size; ++index_k1) {
  
  
  
    /* Allocate k2 level */
    int k2_size = index_k1 + 1;

		class_alloc(psp2->k_true_physical_start_k1k2[index_k1], k2_size*sizeof(int), psp2->error_message);
	  class_alloc(psp2->k_physical_start_k1k2[index_k1], k2_size*sizeof(int), psp2->error_message);
    class_alloc(psp2->k_physical_size_k1k2[index_k1], k2_size*sizeof(int), psp2->error_message);
   
    /* Fill k_size_k1k2, k_min and k_max */
    for(index_k2=0; index_k2<=index_k1; ++index_k2) {
 
      
      int k_pt_size = ppt2->k3_size[index_k1][index_k2];
  		double k_min_pt = ppt2->k3[index_k1][index_k2][0];
 			double k_max_pt = ppt2->k3[index_k1][index_k2][k_pt_size-1];

  		int index_k_sp;

 			 // *** Count the number of necessary values

  			/* First point */
  
  			index_k_sp = 0;
  
  			double k = psp2->k[index_k_sp];
   
  			

 				while (k < k_min_pt && index_k_sp<psp2->k_size -1) {
  				index_k_sp++;
				  k = psp2->k[index_k_sp];
			  }
	
 		 	/* The regime where the triangular condition is satisfied starts here */
  		/*This needs allocation!!!!*/
  		psp2->k_physical_start_k1k2[index_k1][index_k2] = index_k_sp;
	
			index_k_sp = psp2->k_size - 1;
			k = psp2->k[index_k_sp];
			while (k > k_max_pt && index_k_sp > 1) {
  			index_k_sp--;
	  		k = psp2->k[index_k_sp];  
  		}

  		psp2->k_physical_size_k1k2[index_k1][index_k2] = index_k_sp- psp2->k_physical_start_k1k2[index_k1][index_k2]+1;
  		if (psp2->k_physical_size_k1k2[index_k1][index_k2] < 1) printf("Alert empty k3 range for k1 = %f and k2 = %f \n", psp2->k[index_k1],psp2->k[index_k2]);
  		
  		
  	// Here we compute the physical size based on the excat triangular condition
  	
  	 	
  			/* First point */
  
  			index_k_sp = 0; /*?????*/
  
  			k = psp2->k[index_k_sp];
   
  		
 				while (k < (psp2->k[index_k1] - psp2->k[index_k2]) && index_k_sp<psp2->k_size -1) {
  				index_k_sp++;
				  k = psp2->k[index_k_sp];
			  }
	
 		 	
  		psp2->k_true_physical_start_k1k2[index_k1][index_k2] = index_k_sp;
	
		
  		
  		
  	} // end of for(index_k2)

  } // end of for(index_k1)

  
  

  /*add some debug*/   
  return _SUCCESS_;
  
  
} // end of transfer2_get_k3_size


int spectra2_integrate_fourier(
			struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct primordial * ppm,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2

)
{

		/* Vector that will contain the interpolated sources for a given k1,k2 and source type */
  double ** interpolated_sources_in_k;
  class_alloc (interpolated_sources_in_k, ppt2->tp2_size*sizeof(double *), psp2->error_message);
  
  /* Spline coefficients for the source interpolation */
	double ** sources_k_spline;
  class_alloc (sources_k_spline, ppt2->tp2_size*sizeof(double *), psp2->error_message);
  
  double step_k1,step_k2;
  // loop over k1
  for (int index_k1 = 0; index_k1 < psp2->k_size; ++index_k1) {
		if (psp2->spectra2_verbose > 1)
      printf (" computing integral for index_k1=%d of %d, k1=%g\n",
        index_k1, psp2->k_size, psp2->k[index_k1]);
        
        
        if ((ppr2->load_sources_from_disk == _TRUE_) || (ppr2->store_sources_to_disk == _TRUE_))
      		class_call(perturb2_load_sources_from_disk(ppt2, index_k1),
          	ppt2->error_message,
          	psp2->error_message);
        
				// compute stepsize in k1
        if (index_k1 == 0) step_k1 = (psp2->k[1] - psp2->k[0])/2.;
        else if (index_k1 == psp2->k_size -1) step_k1 = (psp2->k[psp2->k_size-1] - psp2->k[psp2->k_size -2])/2.;
        else step_k1 = (psp2->k[index_k1+1] - psp2->k[index_k1 -1])/2.;

				double k1 = psp2->k[index_k1];

				//compute primordial spektrum
				double spectra_k1;
    		class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k1,
                  &spectra_k1),
      			ppm->error_message,
      			psp2->error_message);

    				/* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    				spectra_k1 = 2*_PI_*_PI_/(k1*k1*k1) * spectra_k1;
  
  	//loop over k2
    for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
			/* Print some info */
      if (psp2->spectra2_verbose > 2)
        printf(" -> computing integral for (k1,k2) = (%.3g,%.3g)\n", psp2->k[index_k1], psp2->k[index_k2]);
				if (index_k1 == 0) step_k2 = 0.;
				else if (index_k2 == 0) step_k2 = (psp2->k[1] - psp2->k[0])/2.;
        else if (index_k2 == index_k1) step_k2 = (psp2->k[index_k1] - psp2->k[index_k1-1])/2.;
        else step_k2 = (psp2->k[index_k2+1] - psp2->k[index_k2 -1])/2.;
        
        
        
        double k2 = psp2->k[index_k2];
				double spectra_k2;
    		class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k2,
                  &spectra_k2),
      			ppm->error_message,
      			psp2->error_message);

    				/* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    				spectra_k2 = 2*_PI_*_PI_/(k2*k2*k2) * spectra_k2;
  
 				
 			if (psp2->spectra2_verbose > 2)
        printf(" -> stepsizes (step_k1,step_k2) = (%.3g,%.3g)\n", step_k1,step_k2);
				
 			
      
      // -----------------------------------------------------------------------
      // -                      Interpolate sources in k                       -
      // -----------------------------------------------------------------------

      for (int index_tp=ppt2->index_tp2_M + lm(0,0); index_tp<=ppt2->index_tp2_M + lm(1,1); ++index_tp) {
      	
      
       class_alloc(
          interpolated_sources_in_k[index_tp],
          psp2->k_size*ppt2->tau_size*sizeof(double),
          psp2->error_message);
       
       class_alloc(
          sources_k_spline[index_tp],
          ppt2->k3_size[index_k1][index_k2]*ppt2->tau_size*sizeof(double),
          psp2->error_message);
       
       /* we use the same grid for k that is also used for k1*/
       	// we interpolate the sources in k3
        class_call (spectra2_interpolate_sources_in_k(
                      ppr,
                      ppr2,
                      ppt,
                      ppt2,
                      psp2,
                      index_k1,
                      index_k2,
                      index_tp,
                      psp2->k,                       /* All the k_grid are filled */ 	
                      sources_k_spline[index_tp],					    
                      interpolated_sources_in_k[index_tp]   /* Will be filled with interpolated values */
                      ),
          psp2->error_message,
          psp2->error_message);
      
  
      	// We now have a point in k1,k2 and interpolation in k3, we are ready to perform the integrtion
      	for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) {
					
					
      		
      		// find the correct stepsize for k2 based on the triangular inequality //
      		double triangular_step_k2;
      		
      		
      		// first region, excluded by triangular equality
      		if (psp2->k[index_k1] < psp2->k[index_k3] /2.) {triangular_step_k2 = 0.;
      		}
      		
      		// second region with growing width
      		else if (psp2->k[index_k1] < psp2->k[index_k3]) {
      		// This region k2 starts from k - k1
      			// is k_2 next to the boundary?
      			if (index_k2 < psp2->k_true_physical_start_k1k2[index_k3][index_k1] ) {
      				triangular_step_k2 = 0.;
      				
      			} 
      			else if (index_k2 == psp2->k_true_physical_start_k1k2[index_k3][index_k1]) {
      				if (index_k2 == index_k1) {
      					triangular_step_k2 = (2.*psp2->k[index_k1] - psp2->k[index_k3]) ; 
      					
      				}
      				else {

      					triangular_step_k2 = (psp2->k[index_k2+1] - psp2->k[index_k2])/2.
      					 + (psp2->k[index_k2] - (psp2->k[index_k3]- psp2->k[index_k1]));
      					
      					
      					}
      				}
      			else {triangular_step_k2 = step_k2 ;
      				
      			}
      		}
      		
      		//final region with constant width
      		else {
      		// This region starts from k1 - k 
      			if (index_k2 < psp2->k_true_physical_start_k1k2[index_k1][index_k3] ) {
      				triangular_step_k2 = 0.;

      			} 
      		  else if (index_k2 == psp2->k_true_physical_start_k1k2[index_k1][index_k3]) {
      				if (index_k2 == index_k1) {
      				      			
      					triangular_step_k2 = psp2->k[index_k3] ;
      				
      				}
      				else {
      				      			
      					triangular_step_k2 = (psp2->k[index_k2+1]-psp2->k[index_k2])/2.
      					 + (psp2->k[index_k2] -psp2->k[index_k1] + psp2->k[index_k3]) ;
      					
      				}
      			}
      			else {triangular_step_k2 = step_k2;
      			      			

      			}
      		}					
					
					//loop over time, no integration
					for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
						psp2->spectra[index_tp][index_k3*ppt2->tau_size + index_tau] += 
					// (index_k1==53 ? triangular_step_k2: 0.)
					// symmetry factor (doing only half plane in k1 k2)
					 2.*
					// integration weight	
					 k1*k2/psp2->k[index_k3]*
					// integration stepsize
					 step_k1*triangular_step_k2*
					// phi integration
					 2 * _PI_*
					// spectra factors
					 2./2./_PI_/2./_PI_/2./_PI_* 
					// definition of second order perturbation theory
					 1./2./2.*
					 // m = +- not applicable here
					// psi od phi
					// primordial spectra
					 spectra_k1*spectra_k2*
					//sources
					 interpolated_sources_in_k[index_tp][index_k3*ppt2->tau_size + index_tau]*interpolated_sources_in_k[index_tp][index_k3*ppt2->tau_size + index_tau]
					// 6.652 * 1.e-29 is sigma_t in m^2, 1.878 * 1.e-26 is critical density in kg/m^3
					// 1.602 * 1.e-19 is electron charge ,  1.e-6 is e-10 from CLASS convention 
					// (rho_g_CLASS^0 = e+10 h^2/c^2 * Omega_g^0) and e+4 to covert from international units to Gauss
					// (1 G = 10 e-4 kg/ s/ C)
					//*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
					//*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
					;
					} 
					
		      
      	}
      
      } // end of for (index_tp)

                    
     
      /* Free the memory for the interpolated sources */
      for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
        free(interpolated_sources_in_k[index_tp]);
        free(sources_k_spline[index_tp]);
      }

    } // end of for(index_k2)

   
    class_call (perturb2_free_k1_level (ppt2, index_k1), ppt2->error_message, ppt2->error_message);
      	

    
  } // end of for(index_k1)

	free (interpolated_sources_in_k);
  free(sources_k_spline);
  
	return _SUCCESS_;
	
}





int spectra2_integrate_fourier_sym(
			struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct primordial * ppm,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2
)
{

	/* Vector that will contain the interpolated sources for a given k1,k2 and source type */
  double ** interpolated_sources_in_k;
  class_alloc (interpolated_sources_in_k, ppt2->tp2_size*sizeof(double *), psp2->error_message);
  
  /* Spline coefficients for the source interpolation */
	double ** sources_k_spline;
  class_alloc (sources_k_spline, ppt2->tp2_size*sizeof(double *), psp2->error_message);
  
	double step_k1,step_k3;
	// loop over k1t (here we call k1t the transformed k1, which is indexed by index_k1)
	for (int index_k1 = 0; index_k1 < psp2->k_size; ++index_k1) {
		if (psp2->spectra2_verbose > 1)
      printf (" computing integral for index_k1=%d of %d, k1t=%g\n",
        index_k1, psp2->k_size, psp2->k[index_k1]);
        
    // load sources    
    if ((ppr2->load_sources_from_disk == _TRUE_) || (ppr2->store_sources_to_disk == _TRUE_))
      		class_call(perturb2_load_sources_from_disk(ppt2, index_k1),
          	ppt2->error_message,
          	psp2->error_message);
        

        
		double k1t = psp2->k[index_k1];
		
		//loop over k3t 		
    for (int index_k3 = 0; index_k3 < psp2->k_size; ++index_k3) {
			if (psp2->spectra2_verbose > 2)
        printf(" -> computing integral for (k1t,k3t) = (%.3g,%.3g)\n", psp2->k[index_k1], psp2->k[index_k3]);
        
      double step_k3;
        
			if (index_k3 == 0) step_k3 = (psp2->k[1] - psp2->k[0])/2.;
      else if (index_k3 == psp2->k_size -1) step_k3 = (psp2->k[psp2->k_size-1] - psp2->k[psp2->k_size -2])/2.;
      else step_k3 = (psp2->k[index_k3+1] - psp2->k[index_k3 -1])/2.;

    	double k3t = psp2->k[index_k3];
        
 			
      
      // -----------------------------------------------------------------------
      // -                      Interpolate sources in k                       -
      // -----------------------------------------------------------------------

      for (int index_tp=ppt2->index_tp2_M + lm(0,0); index_tp<=ppt2->index_tp2_M + lm(1,1); ++index_tp) {
      
      // find interpolation for the specific case of k1t,k3t,tp. This is an interpolation in k2t	
      
      class_alloc(
          interpolated_sources_in_k[index_tp],
          psp2->k_size*ppt2->tau_size*sizeof(double),
          psp2->error_message);
       
       class_alloc(
          sources_k_spline[index_tp],
          ppt2->k_size*ppt2->tau_size*sizeof(double),
          psp2->error_message);
       
       /* we use the same grid for k that is also used for k1t*/
       
        class_call (spectra2_interpolate_sources_in_k2(
                      ppr,
                      ppr2,
                      ppt,
                      ppt2,
                      psp2,
                      index_k1,
                      index_k3,
                      index_tp,
                      psp2->k,                   /*this is actually not used properly in our interpolation, clean up!*/    /* All the k_grid are filled */ 	
                      sources_k_spline[index_tp],					    
                      interpolated_sources_in_k[index_tp]   /* Will be filled with interpolated values */
                      ),
          psp2->error_message,
          psp2->error_message);
      
  
      	// here we loop over k3, not k2t. This is output. 
      	/* We now have a position in k1t,k3t and an interpolation in k2t. This translates into
      	knowing k1,k2 and have an interpolation in k3. So we are ready for integration*/
      	for (int index_k = 0; index_k < psp2->k_size; ++index_k) {
					int print = 0;
					
					// construct the non transformed variables
					double k = psp2->k[index_k];
					double k1 =  (-k1t + k3t + 2.*k)/2.;
					double k2 = (k1t + k3t)/2.;
					 
					
					
        	double spectra_k1;
        	
        	if (psp2->spectra2_verbose > 2)
        		printf(" -> analysing k = %f, k1 = %f, k2 = %f \n", k,k1,k2);
       
       			
       		// using the transformed variables k1 or k2 might be negative. This removes that contribution.
        	if (k1 > 0) {
    				class_call (primordial_spectrum_at_k(
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k1,
                  &spectra_k1),
      				ppm->error_message,
      				psp2->error_message);

    				/* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    					spectra_k1 = 2*_PI_*_PI_/(k1*k1*k1) * spectra_k1;
  				} else spectra_k1 = 0.;
        
					double spectra_k2;
    			
    			if (k2 >0 ) {
    				class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k2,
                  &spectra_k2),
      				ppm->error_message,
      				psp2->error_message);

    				/* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    					spectra_k2 = 2*_PI_*_PI_/(k2*k2*k2) * spectra_k2;
  				} else spectra_k2 = 0.;
      		
      		// find the correct stepsize for k1t based on the triangular inequality //
      		double triangular_step_k1;
      		double max;
      		/*the allowed region ranges from 0 to 2k. Here we use the symmetry and restrict k1t from k to 2k. In general it would be 
      		better to use the full range but in that case the interpolation has to be reworked. 0..k has the best grid, 
      		but the worst k2t interpolation (k..2k) . In any case both ways are asymmetric, either have a good interpoaltion 
      		or a good grid for integration. Therefore best go the full way. This does not matter on a linear grid.*/
      		
      		if (k1t < k) {triangular_step_k1 = 0.;}
      		else if (k1t >=k && k1t <= 2.*k) {
      			/*keep in mind that we do not want to go beyond kmax. Therfore we check if 2k 
      			is larger and then use kmax instead for upper limit*/
      			
      			if (2.*k > ppt2->k[ppt2->k_size -1]) 
      			 	max = ppt2->k[ppt2->k_size -1]; 
      			else max = 2.*k; 
      			if (( index_k1 == 0 || ppt2->k[index_k1 - 1] < k) && ( index_k1 == ppt2->k_size - 1 || ppt2->k[index_k1 + 1] > 2.*k)) {
      				/*if a single point is both first and last*/
      				triangular_step_k1 = (  max - ppt2->k[index_k1])  +  (ppt2->k[index_k1] - k); 
      					 
      		 	}
      		
      			else if ( index_k1 == 0 || ppt2->k[index_k1 - 1] < k) {
      			/*first legal point*/
      				triangular_step_k1 = (ppt2->k[index_k1] - k) + (ppt2->k[index_k1+1]- ppt2->k[index_k1])/2.;
      			}
      			else if ( index_k1 == ppt2->k_size - 1 || ppt2->k[index_k1 + 1] > 2.*k) {
      			/*last legal point*/
      				
      				triangular_step_k1 = ( max - ppt2->k[index_k1]) + (ppt2->k[index_k1]- ppt2->k[index_k1-1])/2.;
      			}
      			else{
      				triangular_step_k1 = (ppt2->k[index_k1+1]- ppt2->k[index_k1-1])/2.;
      		 	}
      		}
      		else triangular_step_k1 = 0.;
      		
      		//loop over time (no integration)
      		for (int index_tau = 0; index_tau < ppt2->tau_size; ++index_tau) {
						psp2->spectra[index_tp][index_k*ppt2->tau_size + index_tau] += 
					// (index_k1==53 ? triangular_step_k2: 0.)
					// symmetry factor (doing only half plane in k1 k3)
					 2.*
					// integration weight	
					 k1*k2/k/2.* 
					// integration stepsize
					 step_k3*triangular_step_k1*
					// phi integration
					 2 * _PI_*
					 // m = +-
					 // nothing here, we simplify compute the fourier power spectrum off all tp2
					 // then the user has to do the sum over m himself if he so desires. 
					// spectra factors
					 2./2./_PI_/2./_PI_/2./_PI_* 
					// definition of second order perturbation theory
					 1./2./2.*
					// primordial spectra
					 spectra_k1*spectra_k2*
					//sources
					 interpolated_sources_in_k[index_tp][index_k*ppt2->tau_size + index_tau]*interpolated_sources_in_k[index_tp][index_k*ppt2->tau_size + index_tau]
					// 6.652 * 1.e-29 is sigma_t in m^2, 1.878 * 1.e-26 is critical density in kg/m^3
					// 1.602 * 1.e-19 is electron charge ,  1.e-6 is e-10 from CLASS convention 
					// (rho_g_CLASS^0 = e+10 h^2/c^2 * Omega_g^0) and e+4 to covert from international units to Gauss
					// (1 G = 10 e-4 kg/ s/ C)
					//*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
					//*_c_*_c_*_c_* 6.652 * 1.e-29 * 1.e-6 * 1.878 * 1.e-26 / 1.602 / 1.e-19
					;
					} 
					
		      
      	}
      
      } // end of for (index_tp)

                    
     
      /* Free the memory for the interpolated sources */
      for (int index_tp=0; index_tp<ppt2->tp2_size; ++index_tp) {
        free(interpolated_sources_in_k[index_tp]);
        free(sources_k_spline[index_tp]);
      }

    } // end of for(index_k2)

   
    class_call (perturb2_free_k1_level (ppt2, index_k1), ppt2->error_message, ppt2->error_message);
      	

    
  } // end of for(index_k1)
  
  
  free (interpolated_sources_in_k);
  free(sources_k_spline);
  return _SUCCESS_;
	
}

int spectra2_interpolate_sources_in_k2(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2,
      int index_k1,
      int index_k3,
      int index_tp2,
      double * k_grid,
      double * sources_k_spline,
      double * interpolated_sources_in_k
      )
{


 
  
  
  int index_tau;
  
 
  
  int k_pt_size = ppt2->k_size;
  double * k_pt = ppt2->k;

	double k1 = ppt2->k[index_k1];
  
  
  
  /*We Have to rewrite the sources in order to allow for interpolation.*/
 
 	double * rsources;
 

	class_calloc (
      rsources,
        ppt2->tau_size*(index_k1+1),
        sizeof(double),
        psp2->error_message);
  
  double * derivs;
 

	class_calloc (
      derivs,
        ppt2->tau_size*(index_k1+1),
        sizeof(double),
        psp2->error_message);
    
   
  double * subk;
  
  class_calloc (
      subk,
        (index_k1+1),
        sizeof(double),
        psp2->error_message);      	
				
	for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2){
		subk[index_k2] = ppt2->k[index_k2];
		for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++){
				rsources[index_tau*(index_k1+1) + index_k2] = sources(index_tau,index_k3);
		}
	}
	
	
  if (ppr2->sources_k3_interpolation == cubic_interpolation && index_k1> 3) {

	
        
    class_call (array_spline_table_columns (
                  subk,
                  (index_k1+1),
                  rsources,
                  ppt2->tau_size,
                  derivs,
                  _SPLINE_EST_DERIV_,
                  psp2->error_message),
         psp2->error_message,
         psp2->error_message);
  }
 

  free(subk);
  // =======================================================
  // =                    Interpolation                    =
  // =======================================================

  
  /* Interpolate at each needed k2t value corresponding to one k3 and k1t value using the usual spline interpolation algorithm */
 	int index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
  
  int index_k_sp;

	 /*this index represents k3*/ 
  for (index_k_sp = 0; index_k_sp < ppt2->k_size; ++index_k_sp) {
    
    double k2t = - ppt2->k[index_k1] + 2.* ppt2->k[index_k_sp];
    index_k = 0;
    if (k2t < 0 || k2t > k1 || k2t > ppt2->k[ppt2->k_size-1]) {
    	
    	interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 0.; /*no extrapolation out of triangular, might be changed to flat extrapolation for trinagular (NOT FOR KMAX !!!)*/
    } else {
   
    	while ((k_pt[index_k+1] < k2t) && (index_k + 1 < ppt2->k_size) ){
     	 index_k++;
     	 h = k_pt[index_k+1] - k_pt[index_k];
    	}
    
    	class_test(h==0., psp2->error_message, "stop to avoid division by zero");
    
    	double b = (k2t- k_pt[index_k])/h; 
 
    	double a = 1.-b;
    	
      
    /* Interpolate for each value of conformal time */
    	if (ppr2->sources_k3_interpolation == linear_interpolation ||  index_k1 <= 3) {
      	for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
       	 interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
       	   a * rsources[index_tau*(index_k1+1) + index_k] + b * rsources[index_tau*(index_k1+1) + index_k+1];
    	}
    	else if (ppr2->sources_k3_interpolation == cubic_interpolation && index_k1 > 3) {
     	 for (index_tau = 0; index_tau < ppt2->tau_size; index_tau++)
     	   interpolated_sources_in_k[index_k_sp*ppt2->tau_size + index_tau] = 
     	     a * rsources[index_tau*(index_k1+1) + index_k] + b * rsources[index_tau*(index_k1+1) + index_k+1]
      	    + ((a*a*a-a) * derivs[index_tau*(index_k1+1) + index_k]
      	    +(b*b*b-b) * derivs[index_tau*(index_k1+1) + index_k+1])*h*h/6.0;
    	}
    
   } 
  } // end of for (index_k_tr)



		

	free(rsources);
 	free(derivs);
 
  return _SUCCESS_;
  
} // end of spectra_interpolate_sources_in_k


#undef sources




