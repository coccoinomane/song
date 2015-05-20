#ifndef __SLATEC_3J__
#define __SLATEC_3J__

#include "common.h"

int drc3jm (double L1, double L2, double L3, double M1, /* IN, parameters */
            double *M2MIN, double *M2MAX,               /* OUT, limits of m2 */
            double * THRCOF,                            /* OUT, result */
            int NDIM,                                   /* IN, length of result */
            // int * IER,                                  /* OUT, error code, 0 ok */
            ErrorMsg errmsg);


int drc3jj (double L2, double L3, double M2, double M3, /* IN, parameters */
            double *L1MIN, double *L1MAX,               /* OUT, limits of l1 */
            double * THRCOF,                            /* OUT, result */
            int NDIM,                                   /* IN, length of result */
            // int * IER,                                  /* OUT, error code, 0 ok */
            ErrorMsg errmsg);

int drc6j (double L2, double L3, double L4, double L5, double L6, /* IN, parameters */
           double *L1MIN, double *L1MAX,               /* OUT, limits of l1 */
           double * SIXCOF,                            /* OUT, result */
           int NDIM,                                   /* IN, length of result */
           // int * IER,                                  /* OUT, error code, 0 ok */
           ErrorMsg errmsg);
           
#endif