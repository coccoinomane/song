
#ifndef __SPECTRA2__
#define __SPECTRA2__

#include "primordial.h"
#include "perturbations.h"
#include "perturbations2.h"
#include "transfer2.h"
#include "transfer.h"
#include "bessel2.h"

struct spectra2 {




  
  double ** spectra; 
  
  
 
  /* For a given (k1,k2), index of the first value of k3 that satisfies the triangular condition. All
  entries must be equal zero when no extrapolation is used */
  int ** k_physical_start_k1k2;
	int ** k_true_physical_start_k1k2;
	
	// physical start is based on ppt2 k sampling, true physical start is based on the triangular condition and usually contains a bit more based on the k smapling. 

  /* For a given (k1,k2), number of k3 values that satisfy the triangular condition */
  int ** k_physical_size_k1k2;

 	int k_size;
 	double * k; // k grid for fourier spectra (not numerical parameter)
 	
 	double * k3_grid; // k grid for angular power spectra 
 	


  // =================================================================================
  // =                        Storage of intermediate results                        =
  // =================================================================================
  

 char spectra_filename[_FILENAMESIZE_];  
 FILE * spectra_file;

 

  // ====================================================================================
  // =                                 Debug parameters                                 =
  // ====================================================================================

  short spectra2_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */



};


int spectra2_init(
 			struct primordial * ppm,
      struct precision * ppr,
      struct precision2 * ppr2,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct bessels *pbs,
      struct bessels2 *pbs2,
      struct transfers * ptr,
      struct transfers2 * ptr2,
      struct spectra2 * psp2
      );

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
      );
      

int spectra2_interpolate_sources_in_k2(
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs * ppt,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2,
      int index_kt1,
      int index_kt3,
      int index_tp2,
      double * k_grid,
      double * sources_k_spline,
      double * interpolated_sources_in_k
      );

int spectra2_get_k3_size (
      struct precision * ppr,
      struct precision2 * ppr2,
      struct perturbs2 * ppt2,
      struct spectra2 * psp2
      );
      
#endif
