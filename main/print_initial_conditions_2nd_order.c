/** @file print_initial_conditions_2nd_order.c 
 * Created by Guido Walter Pettinari, 7.08.2011
 * Last modified by Guido Walter Pettinari, 7.08.2011
 *
 * Print the initial conditions of the second-order Boltzmann-Einstein differential system as a function of time.
 * 
 * This script accepts the following arguments:
 * ./print_initial_conditions_2nd_order <ini file> [pre file] <index_k1> <index_k2> <index_cosk1k2>
 *
 */
 
#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions (1st-order) */
  struct perturbs2 pt2;       /* for source functions (2nd-order) */  
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&pt2,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  // Set verbosity to zero (we want to be able to send the output to a file and plot it)
  ba.background_verbose = 0;
  th.thermodynamics_verbose = 0;
  pt.perturbations_verbose = 0;
  pt2.perturbations2_verbose = 0;  

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (perturb2_init(&pr,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }


  // ===================
  // = Parse arguments =
  // ===================
  // Print arguments (debug)
  // int jj=0;
  // for (jj=0; jj<argc; ++jj)
  //   printf("argv[%d] = %s\n", jj, argv[jj]);
  
  int index_mode = 0;
  int index_ic = 0;
  int index_k1, index_k2, index_cosk1k2;
  
  // Allow for optional .pre file
  if (argc==5) {
    index_k1 = atoi(argv[2]);
    index_k2 = atoi(argv[3]);
    index_cosk1k2 = atoi(argv[4]);
  }
  else if (argc==6) {
    index_k1 = atoi(argv[3]);
    index_k2 = atoi(argv[4]);
    index_cosk1k2 = atoi(argv[5]);
  }
  else {
    printf("usage:   %s <ini file> [pre file] <index_k1> <index_k2> <index_cosk1k2>\n", argv[0]);
    return _FAILURE_;
  }
  
  // Sizes associated to the non-running indices in the ***sources table
  int k_size = pt2.k_size[index_mode];
  int cosk1k2_size = pt2.cosk1k2_size[index_mode];    
  
  // Account for overshooting of k
  if(index_k1 > k_size-1) index_k1 = k_size-1;
  if(index_k1 < 0) index_k1 = 0;
  if(index_k2 > k_size-1) index_k2 = k_size-1;
  if(index_k2 < 0) index_k2 = 0;
  if(index_cosk1k2 > cosk1k2_size-1) index_cosk1k2 = cosk1k2_size-1;
  if(index_cosk1k2 < 0) index_cosk1k2 = 0;
  
  // Find k-values
  double k1 = pt2.k[index_mode][index_k1];
  double k2 = pt2.k[index_mode][index_k2];
  double cosk1k2 = pt2.cosk1k2[index_mode][index_cosk1k2];    
  
  
  // =================
  // = Print columns =
  // =================
  // Hard-coded option to choose whether to have the complete list of perturbations,
  // or just the ones set by hand in perturbs2_initial_conditions
  int print_all_pt = _FALSE_;
  // These are the printed perturbs if print_all_pt = _FALSE_ :
  // 0 - delta_g
  // 1 - theta_g
  // 2 - delta_b
  // 3 - theta_b
  // 4 - delta_cdm
  // 5 - delta_neutrinos
  // 6 - theta_neutrinos
  // 7 - shear_neutrinos
  // 8 - metric variable eta
  int reduced_pt2_size = 9;
  
  // Some debug info
  printf("# k1 = %g, k2 = %g, cosk1k2 = %g\n", k1, k2, cosk1k2);
  printf("# index_mode = %d\n", index_mode);
  printf("# index_ic = %d\n", index_ic);
  
  // We need to create a workspace structure in order to extract the initial conditions,
  // since the function perturb2_initial_conditions writes on the perturb_vector structure
  // contained inside a workspace structure, via the function perturb2_vector_init
  struct perturb2_workspace * ppw2;
  class_alloc(ppw2, sizeof(struct perturb2_workspace), pt2.error_message);
  class_call(perturb2_workspace_init(&pr,
        &ba,
        &th,
        &pt,
        &pt2,
        index_mode,
        ppw2
        ),
    pt2.error_message,
    pt2.error_message);    
  
  // The perturbs2_vector_init function is called inside the main time loop in order to fill ppw2->pv with
  // the initial conditions at a specific time.  We call it here for a first time in order to get ppw2->pv->pt2_size
  // and consequently fill the array reduced_pt2_indices with the indices of the variables to print.
  class_call(perturb2_vector_init(&pr,
          &ba,
          &th,
          &pt,
          &pt2,
          index_mode,
          index_ic,
          k1,
          k2,
          cosk1k2,
          ba.tau_table[0],
          ppw2,
          NULL),
      pt2.error_message,
      pt2.error_message);
  
  // Indices of the variables to print.  The indices index the vector ppw2->pv->y.
  int reduced_pt2_indices[] = {
      ppw2->pv->index_pt2_delta_g,
      ppw2->pv->index_pt2_theta_g,
      ppw2->pv->index_pt2_delta_b,
      ppw2->pv->index_pt2_theta_b,
      ppw2->pv->index_pt2_delta_cdm,
      ppw2->pv->index_pt2_delta_ur,
      ppw2->pv->index_pt2_theta_ur,
      ppw2->pv->index_pt2_shear_ur,
      ppw2->pv->index_pt2_eta
  };  
  
  // More debug information
  int size = print_all_pt ? ppw2->pv->pt2_size : reduced_pt2_size;
  printf("# considering %d out of %d pt2 variables\n", size, ppw2->pv->pt2_size);
  
  // Print a first row with the name of the variables
  char reduced_pt2_labels[][64] = 
    {"delta_g", "theta_g", "delta_b", "theta_b", "delta_cdm", "delta_ur", "theta_ur", "shear_ur", "eta"};
  int ii;
  // First row is time
  if (!print_all_pt) {
    printf("%13s ", "tau");
    for(ii=0; ii<reduced_pt2_size; ++ii)
      printf("%13s ", reduced_pt2_labels[ii]);
    printf("\n");
  }
  
  
  // ================
  // = Loop on time =
  // ================
  // We are going to use the same time sampling used to determine background quantities
  int index_tau;
  double tau;
  for(index_tau=0; index_tau<ba.bt_size; ++index_tau) {
    
    tau = ba.tau_table[index_tau];
    
    class_call(perturb2_vector_init(&pr,
          &ba,
          &th,
          &pt,
          &pt2,
          index_mode,
          index_ic,
          k1,
          k2,
          cosk1k2,
          tau,
          ppw2,
          NULL),
      pt.error_message,
      pt2.error_message);
      
    // First column is time.  The plus indicates that even plus signs will be printed.
    printf ("%+6e ", tau);
   
    // Columns from 2 to pt2_size+1 are the perturbed quantities
    int index_pt;
  
    if (print_all_pt) {
      for (index_pt=0; index_pt < ppw2->pv->pt2_size; ++index_pt) {
        double temp = ppw2->pv->y[index_pt];
        printf ("%+6e ", temp);
      }
    }
    else {
      for (index_pt=0; index_pt < reduced_pt2_size; ++index_pt) {
        double temp = ppw2->pv->y[reduced_pt2_indices[index_pt]];
        printf ("%+6e ", temp);
      }
      
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_delta_g]);      // 0 - delta_g
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_theta_g]);      // 1 - theta_g
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_delta_b]);      // 2 - delta_b
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_theta_b]);      // 3 - theta_b
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_delta_cdm]);    // 4 - delta_cdm
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_delta_ur]);     // 5 - delta_neutrinos
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_theta_ur]);     // 6 - theta_neutrinos
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_shear_ur]);     // 7 - shear_neutrinos
      // printf("%+6e ", ppw2->pv->y[ppw2->pv->index_pt2_eta]);          // 8 - metric variable eta
    }
    
    // Insert a carriage return
    printf("\n");
    
  }
  
  class_call(perturb2_workspace_free(&pt2,index_mode,ppw2),
      pt2.error_message,
      pt2.error_message);
  
  
  /****** all calculations done, now free the structures ******/
  
  if (perturb2_free(&pr2,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }
  
  if (perturb_free(&pr,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }
  
  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }
  
  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  
  return _SUCCESS_;

}
