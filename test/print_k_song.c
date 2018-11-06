/** @file print_k_song.c 
 *
 * Print to screen the k-sampling used to solve the second-order differential
 * system in the perturbations2.c module.
 *
 * Arguments:
 * -# The .ini file
 * -# The .pre file (optional)
 * -# The index of the mode (scalar, vector, tensor)
 * -# Choose between 0 and 1 to print either k or cosk1k2, respectively.
 *
 * usage:     %s <ini file> [<pre file>] <index mode S,V,T> <index wavemode k,cosk1k2>
 *
 * Created by Guido Walter Pettinari on 01/08/2011
 * Last modified by Guido Walter Pettinari on 30/06/2015
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


  // ===================
  // = Parse arguments =
  // ===================

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

  // // Check that the index_mode is correct
  // if ((index_mode >= pt2.md2_size) || (index_mode < 0)) {
  //   printf ("ERROR: index_mode=%d is out of range. It should be between 0 and %d.\n", index_mode, pt2.md2_size-1);
  //   return _FAILURE_;
  // }

  // // Print some useful info
  // printf ("# %s modes,  ", pt2.md2_labels[index_mode]);
      
      
  // // ====================
  // // = Print k-sampling =
  // // ====================
  // // Print k-sampling according to the chosen wavemode (k or cosk1k2)
  // int ii;
  // if (index_wavemode==0) {
  //   printf("k1\n");
  //   for (ii=0; ii<pt2.k_size[index_mode]; ++ii)
  //     printf ("%g\n", pt2.k[index_mode][ii]); }
  // else if (index_wavemode==1) {
  //   printf("cosk1k2\n");
  //   for (ii=0; ii<pt2.cosk1k2_size[index_mode]; ++ii)
  //     printf ("%g\n", pt2.cosk1k2[index_mode][ii]); }
  // else {
  //   printf ("ERROR: index_wavemode=%d is out of range. It should be between 0 and %d.\n", index_wavemode, 1);
  //   return _FAILURE_;
  // }


  // =================================================================================
  // =                                  Free memory                                  =
  // =================================================================================

  if (perturb2_free(&pr2,&pt2) == _FAILURE_) {
    printf("\n\nError in perturb2_free \n=>%s\n",pt2.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
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
