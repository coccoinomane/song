/** @file print_k.c 
 *
 * Print to screen the k-sampling used to solve the differential system in the
 * perturbations.c module.
 *
 * Arguments:
 * -# The .ini file
 * -# The .pre file (optional)
 * -# The index of the mode (scalar, vector, tensor)
 * 
 * Created by Guido Walter Pettinari on 24/06/2011
 * Last edited by Guido Walter Pettinari on 30/06/2015
 *
 * IMPORTANT: this file won't compile; it has to be updated it
 * for the latest version of SONG.
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

  // Print arguments (debug)
  // int jj=0;
  // printf ("argc = %d\n", argc);
  // for (jj=0; jj<argc; ++jj)
  //   printf("argv[%d] = %s\n", jj, argv[jj]);

  // Parse arguments.  First argument is always the ini file for CLASS, while last argument is always index_mode.
  int index_mode = 0;
  if (argc == 2)
    index_mode = 0;
  else if (argc == 3)
    index_mode = atoi(argv[2]);
  else if (argc == 4)
    index_mode = atoi(argv[3]);
  else {
    printf ("usage:     print_k <ini file> [<pre file>] <index mode>\n");
    return 1;
  }
  printf ("# index mode = %d\n", index_mode);

  // Compute perturbations
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

  // Check that the index_mode is correct
  if ((index_mode >= pt.md_size) || (index_mode < 0)) {
    printf ("ERROR: index_mode=%d is out of range. It should be between 0 and %d.\n", index_mode, pt.md_size-1);
    return 1;
  }
      
  // Print k sampling
  int ii;

  for (ii=0; ii<pt.k_size[index_mode]; ++ii)
    printf ("%g\n", pt.k[index_mode][ii]);



  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

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
