/** @file input2.h Documented includes for input2 module */

#ifndef __INPUT2__
#define __INPUT2__

#include "input.h"
#include "common2.h"
#include "bessel2.h"
#include "perturbations2.h"
#include "transfer2.h"

/* macros for reading parameter values with routines from the parser */
#define class_read_string_one_of_two(pfc,name1,name2)                   \
  do {                                                                  \
    class_call(parser_read_string(pfc,name1,&(string1),&(flag1),errmsg),\
         errmsg,                                                        \
         errmsg);                                                       \
    class_call(parser_read_string(pfc,name2,&(string2),&(flag2),errmsg),\
         errmsg,                                                        \
         errmsg);                                                       \
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),                  \
         errmsg,                                                        \
         "In input file, you can only enter one of %s, %s, choose one", \
         name1,name2);                                                  \
    /* string2 wins over string 1 */                                    \
    flag = _TRUE_;                                                      \
    if (flag2 == _TRUE_)                                                \
      strcpy (string, string2);                                         \
    else if (flag1 == _TRUE_)                                           \
      strcpy (string, string1);                                         \
    else                                                                \
      flag = _FALSE_;                                                   \
  } while(0);


/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

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
		 struct bessels *pbs,
		 struct bessels2 *pbs2,
     struct transfers2 *ptr2,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct bispectra *pbi,
     struct fisher *pfi,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input2_init(
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
		 struct nonlinear *pnl,
		 struct lensing *ple,
		 struct bispectra *pbi,
     struct fisher *pfi,
		 struct output *pop,
		 ErrorMsg errmsg
		 );

  int input2_default_params(
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
			   struct output *pop
			   );
  
  int input2_default_precision(
			      struct precision2 * ppr2
			      );

  int input2_free(
            struct precision2 * ppr2
            );

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
