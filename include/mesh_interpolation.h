#ifndef __MESH__
#define __MESH__

#include "common.h"

/* Header and definitions needed to compile */
// #include "math.h"
// #include "stdio.h"
// #define MIN(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
// #define MAX(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */
// #define _SUCCESS_ 0
// #define _FAILURE_ 0


struct mesh_interpolation_workspace {

  long int n_points;
  long int n_boxes;
  double l_max;
  double link_length;
  double soft_coeff;
  double group_length;
  
  int *** grid;
  double ***** mesh;
  
  /* Should we compute the grid? If not, the user should overwrite by hand the grid field
    with a precomputed grid. */
  short compute_grid;
  
  /* Counters to keep track of the memory usage */
  long int n_allocated_in_mesh;
  long int n_allocated_in_grid;

  /* Counters to keep track of which meshes are used */
  long int count_interpolations;
  long int count_range_extensions;

  
};




/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif


  int mesh_sort (
      struct mesh_interpolation_workspace * pw,
      double ** values
      );

  int mesh_int (
      struct mesh_interpolation_workspace * pw,
      double x,
      double y,
      double z,
      double * interpolated_value
      );

  int mesh_free (
      struct mesh_interpolation_workspace * pw
      );

  double distance (
      double * vec1,
      double * vec2
  );

#ifdef __cplusplus
}
#endif


#endif