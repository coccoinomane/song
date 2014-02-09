/** @file print_k.c 
 * Created by Guido Walter Pettinari on 01/08/2011
 * Last modified by Guido Walter Pettinari on 10/08/2011
 *
 * Print to screen the k-sampling used to solve the differential system in the perturbations2 module.
 * The first argument is the .ini file, the second (optional) argument is the .pre file, the third argument
 * is the index of the mode (S, V, T) for which you want the k-sampling to be printed, while the third
 * argument is either 0 or 1 and selects respectively between k and cosk1k2
 *
 * usage:     %s <ini file> [<pre file>] <index mode S,V,T> <index wavemode k,cosk1k2>
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


  // ===================
  // = Parse arguments =
  // ===================
  // Print arguments (debug)
  // int jj=0;
  // printf ("argc = %d\n", argc);
  // for (jj=0; jj<argc; ++jj)
  //   printf("argv[%d] = %s\n", jj, argv[jj]);

  // Parse arguments
  int index_mode, index_wavemode;
  if (argc == 4) {
    index_mode = atoi(argv[2]);
    index_wavemode = atoi(argv[3]); }
  else if (argc == 5) {
    index_mode = atoi(argv[3]);
    index_wavemode = atoi(argv[4]); }
  else {
    printf ("usage:     %s <ini file> [<pre file>] <index mode S,V,T> <index wavemode k,cosk1k2>\n", argv[0]);
    return _FAILURE_;
  }



  // =========================
  // = Compute perturbations =
  // =========================
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

  if (perturb2_init(&pr,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  // Check that the index_mode is correct
  if ((index_mode >= pt2.md2_size) || (index_mode < 0)) {
    printf ("ERROR: index_mode=%d is out of range. It should be between 0 and %d.\n", index_mode, pt2.md2_size-1);
    return _FAILURE_;
  }

  // Print some useful info
  printf ("# %s modes,  ", pt2.md2_labels[index_mode]);
      
      
  // ====================
  // = Print k-sampling =
  // ====================
  // Print k-sampling according to the chosen wavemode (k or cosk1k2)
  int ii;
  if (index_wavemode==0) {
    printf("k1\n");
    for (ii=0; ii<pt2.k_size[index_mode]; ++ii)
      printf ("%g\n", pt2.k[index_mode][ii]); }
  else if (index_wavemode==1) {
    printf("cosk1k2\n");
    for (ii=0; ii<pt2.cosk1k2_size[index_mode]; ++ii)
      printf ("%g\n", pt2.cosk1k2[index_mode][ii]); }
  else {
    printf ("ERROR: index_wavemode=%d is out of range. It should be between 0 and %d.\n", index_wavemode, 1);
    return _FAILURE_;
  }



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
