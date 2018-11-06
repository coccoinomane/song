/** @file print_k.c 
 * Guido Walter Pettinari 04/08/2011
 *
 * Print to screen the tau-sampling of the sources in the perturbations2 module.
 * The only argument is the .ini file,
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
  struct nonlinear nl;        /* non-linear spectra */
  struct lensing le;          /* lensed spectra */
  struct bispectra bi;        /* bispectra */
  struct fisher fi;           /* fisher matrix */
  struct output op;           /* output files */
  ErrorMsg errmsg;            /* error messages */

  if (input_init_from_arguments(argc,argv,&pr,&ba,&th,&pt,&pt2,&bs,&tr,&pm,&sp,&nl,&le,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (input2_init_from_arguments(argc,argv,&pr,&pr2,&ba,&th,&pt,&pt2,&tr,&bs,&bs2,&tr2,&pm,
    &sp,&nl,&le,&bi,&fi,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  // Set verbosity to zero (we want to be able to send the output to a file and plot it)
  ba.background_verbose = 0;
  th.thermodynamics_verbose = 0;
  pt.perturbations_verbose = 0;
  pt2.perturbations2_verbose = 0;

  // Print arguments (debug)
  // int jj=0;
  // printf ("argc = %d\n", argc);
  // for (jj=0; jj<argc; ++jj)
  //   printf("argv[%d] = %s\n", jj, argv[jj]);

  // Make sure to stop before computing second-order sources
  pt2.stop_at_perturbations1 = _TRUE_;

  // Compute perturbations
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb2_init(&pr,&pr2,&ba,&th,&pt,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_init \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  // Print tau-sampling according to the chosen wavemode (k1, k2, or cosk1k2)
  int ii;
  printf("pt2.tau_size = %d\n", pt2.tau_size);
  for (ii=0; ii<pt2.tau_size; ++ii) {
    printf ("%g\n", pt2.tau_sampling[ii]);
  }



  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

  if (perturb2_free(&pr2,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
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
