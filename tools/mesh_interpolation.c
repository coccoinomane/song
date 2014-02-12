#include "mesh_interpolation.h"






int mesh_sort (
    struct mesh_interpolation_workspace * pw,
    double ** vals
    )
{

  /*This routine presorts an unsorted mesh. 
  In: ** vals: vals[n][0] is the value of the nth point, vals[n][1] is the x coordinate, [2] is the y, [3] is z.  
  num_points: the number of points
  link_length: the local region influencing the values, on an in homogenous grid this should correspond roughly to the largest distance of two neighbouring points
  soft_coeff: softens the link_length. The default value on most meshs should be 0.5
  group_length: down weights clusters of many close points. In an in homogenous grid this should correspond to the shortest distance of two points
  max_l: the mesh assumes the arguments in each direction should be between 0 and max_l
  Out: grid: this specifies how many points fall in each bin of the new sorted mesh. It will be allocated by this routine automatically
  Return value: the new sorted mesh 
  */

  /* First we create a grid, that is we count the number of particles in each box of side
  equal to the linking length. This is stored in grid[ix][iy][iz]. Then we create the mesh,
  an array that contains the information about the particules contained in each box. It is
  accessed as mesh[ix][iy][ik][n] where 'n' is the ID of the particle in the box, which
  goes from 0 to grid[ix][iy][iz]-1.  The information stored in the mesh has 5 elements.
  The first four are x,y,z and f(x,y,z). The last is the density of points in the box
  around the node n; for each other node in the box, the density gets a contribution of
  e^(-distance squared/grouping_length) so that a group of particles clustered
  on a scale smaller than the grouping length count as one particle. The density is
  needed only to downweight those points that are clustered so that they count as one
  (see mesh_int).
  */
  

  /* Read values from structure */
  long int num_points = pw->n_points;
  double max_l = pw->l_max;
  double link_length = pw->link_length;
  double group_length = pw->group_length;
  double soft_coeff = pw->soft_coeff;

  /* Initialize counters */
  pw->n_allocated_in_grid = 0;
  pw->n_allocated_in_mesh = 0;

  long int i,j,k,m,n,ix,iy,iz,status;
  int *** counter;
  double safe = 1.; /*this can be increased to improve the accuracy, but 1 is already very good if link_length is chosen properly*/
  pw->n_boxes = ceil(max_l/(safe*link_length*(1.+soft_coeff)));
  long int n_boxes = pw->n_boxes;


  /* Bin the support points in a grid based on the linking length */
  if (pw->compute_grid==_TRUE_) {

    pw->grid = (int***) malloc(n_boxes*sizeof(int**));
    for(i=0; i<n_boxes; i++) {
      pw->grid[i] = (int**) malloc(n_boxes*sizeof(int*));
      for(j=0; j<n_boxes; j++) {
        pw->grid[i][j] = (int *) calloc(n_boxes,sizeof(int));
        #pragma omp atomic
        pw->n_allocated_in_grid += n_boxes;
      }
    }
    
    #pragma omp parallel for private (i,ix,iy,iz)
    for (i=0; i<num_points; i++){

      ix = floor(vals[i][1]/(safe*link_length*(1.+soft_coeff)));
      iy = floor(vals[i][2]/(safe*link_length*(1.+soft_coeff)));
      iz = floor(vals[i][3]/(safe*link_length*(1.+soft_coeff)));
    
      /* Skip the data that is larger than max_l */
      if ((ix >= n_boxes) || (iy >= n_boxes) || (iz >= n_boxes))
        continue;
    
      #pragma omp atomic
      pw->grid[ix][iy][iz]++;
    }
    
  } // end of if compute grid


  // *** ALLOCATE COUNTER
  counter = (int***) malloc(n_boxes*sizeof(int**));
  for(i=0;i<n_boxes;i++){
    counter[i] = (int**) malloc(n_boxes*sizeof(int*));
    for(j=0;j<n_boxes;j++){
      counter[i][j] = (int *) calloc(n_boxes,sizeof(int));
    }
  }
  
  

  // *** ALLOCATE MESH
  pw->mesh = (double*****) malloc(n_boxes*sizeof(double****));
  for(i=0;i<n_boxes;i++){
    pw->mesh[i] = (double****) malloc(n_boxes*sizeof(double***));
    for(j=0;j<n_boxes;j++){
      pw->mesh[i][j] = (double ***) malloc(n_boxes*sizeof(double**));
      for(k=0;k<n_boxes;k++){
        pw->mesh[i][j][k] = (double **) malloc((pw->grid[i][j][k])*sizeof(double*));
        for(m=0; m < pw->grid[i][j][k]; m++){
          pw->mesh[i][j][k][m] = (double *) calloc(5, sizeof(double)); 
          #pragma omp atomic
          pw->n_allocated_in_mesh += 5;
        }
      }
    }
  }
  
  
  // *** STORE VALUES IN THE MESH
  #pragma omp parallel for private (i,ix,iy,iz)
  for(i=0; i<num_points; i++){

    ix = floor(vals[i][1]/(safe*link_length*(1.+soft_coeff)));
    iy = floor(vals[i][2]/(safe*link_length*(1.+soft_coeff)));
    iz = floor(vals[i][3]/(safe*link_length*(1.+soft_coeff)));

    /* Skip the data that is larger than max_l */
    if ((ix >= n_boxes) || (iy >= n_boxes) || (iz >= n_boxes))
      continue;

    int Q; for (Q=0; Q < 4; ++Q) {
      pw->mesh[ix][iy][iz][counter[ix][iy][iz]][Q] = vals[i][Q];
      // if ((ix==1) && (iy==1) && (iz==1))
      //   printf("pw->mesh[ix=%d][iy=%d][iz=%d][counter=%d][Q=%d] = %g\n",
      //     ix, iy, iz, counter[ix][iy][iz], Q, pw->mesh[ix][iy][iz][counter[ix][iy][iz]][Q]);
    }

    #pragma omp atomic
    counter[ix][iy][iz]++;
  }

  /* Counter not needed anymore */
  for(i = 0; i<n_boxes;i++) {
    for(j = 0; j<n_boxes; j++)
        free (counter[i][j]);
    free (counter[i]);
  }
  free (counter);



  // *** COMPUTE LOCAL DENSITY OF POINTS
  
  /* Compute the local density of points and store it in the [4]-entry of each point.
  This is a loop over the support points, which means that it can be quick if there
  are not too many. */

  #pragma omp parallel for private (i,j,k,m,n)
  for(i = 0; i<n_boxes;i++) {
    for(j = 0; j<n_boxes; j++) {
      for(k = 0; k<n_boxes; k++) {
        for(m = 0; m<pw->grid[i][j][k]; m++){
          for(n = 0; n<pw->grid[i][j][k];n++){ 

            /* The 4th level of mesh is just the n-th particle in the ijk box */
            double dist = distance(pw->mesh[i][j][k][m],pw->mesh[i][j][k][n]);
            double density = exp(-dist*dist/pow(group_length,2));
                  
            /* The 5th level of mesh is the local density around the m-th particle of the ijk box */
            #pragma omp atomic
            pw->mesh[i][j][k][m][4] += density;
                        
          }
        }
      }
    }
  }
  
  /* Some debug */
  // for(i=0;i<n_boxes;i++){
  //     for(j=0;j<n_boxes;j++){
  //     for(k=0;k<n_boxes;k++){
  //       for(m=0;m<(pw->grid)[i][j][k];m++){
  // 
  //           int Q;
  //           for (Q=0; Q < 5; ++Q) {
  //             printf("%g ", pw->mesh[i][j][k][m][Q]);
  //           }
  //           printf("\n");
  //         }
  //       }
  //     }
  //   }
  
  // printf("~*~*~*~*~ Executing line %d of function %s\n", __LINE__, __func__); fflush(stdout);
  // printf("pw->mesh[0][0][0][0][0] = %g\n", pw->mesh[0][0][0][0][0]);
  // printf("pw->mesh[0][0][0][0][1] = %g\n", pw->mesh[0][0][0][0][1]);
  
  return _SUCCESS_;
  
}






int mesh_int (
    struct mesh_interpolation_workspace * pw,
    double x,
    double y,
    double z,
    double * interpolated_value
    )
{
  
  /*this routine interpolates a point on the sorted mesh
  In: ***** vals: sorted mesh as returned by mesh sort
  link_length: the local region influencing the values, on an in homogenous grid this should correspond roughly to the largest distance of two neighbouring points
  soft_coeff: softens the link_length. The default value on most meshs should be 0.5
  group_length: down weights clusters of many close points. In an in homogenous grid this should correspond to the shortest distance of two points
  *** grid: number of points in each bin, as returned by mesh_sort
  l1,l2,l3: the x,y,z coordinates of the point where the interpolation is needed
  max_l: the mesh assumes the arguments in each direction should be between 0 and max_l
  Out: grid: this specifies how many points fall in each bin of the new sorted mesh. It will be allocated by this routine automatically
  Return value: value of the interpolation
  
  IMPORTANT: All parameters like max_l, link_length etc need to be identical to the values given to mesh_sort. Otherwise you will get either random values or the code will crash 

  */
  
  /* Read values from structure */
  double max_l = pw->l_max;
  double link_length = pw->link_length;
  double group_length = pw->group_length;
  double soft_coeff = pw->soft_coeff;
  double ***** mesh = pw->mesh;
  int *** grid = pw->grid;
  
  
  /* Check bounds */
  // double min_l = mesh[0][0][0][0][1];
  // if ((l1 > max_l) || (l2 > max_l) || (l3 > max_l) || (l1 < min_l) || (l2 < min_l) || (l3 < min_l))
  //   printf ("###### WARNING: interpolation out of bounds: l1,l2,l3 should all be within [%g,%g] (you gave %g,%g,%g)\n", min_l, max_l, l1, l2, l3);
  // class_test ((l1 > max_l) || (l2 > max_l) || (l3 > max_l) || (l1 < min_l) || (l2 < min_l) || (l3 < min_l),
  //             pbi->error_message,
  //             "interpolation out of bounds: l1,l2,l3 should all be within [%g,%g] (you gave %g,%g,%g)", min_l, max_l, l1, l2, l3);
  
  int ix,iy,iz,i;
  double safe = 1.;
  long int n_boxes = pw->n_boxes;
  
  /* Locate the box where (x,y,z) is */
  ix = floor(x/(safe*link_length*(1.+soft_coeff)));
  iy = floor(y/(safe*link_length*(1.+soft_coeff)));
  iz = floor(z/(safe*link_length*(1.+soft_coeff)));

  /* Start by considering only the boxes adjacent to the one which contains (x,y,z) */
  int ixmin = (int) (MAX(0,ix-1));
  int iymin = (int) (MAX(0,iy-1));
  int izmin = (int) (MAX(0,iz-1));
  int ixmax = (int) (MIN(n_boxes-1,ix+1));
  int iymax = (int) (MIN(n_boxes-1,iy+1));
  int izmax = (int) (MIN(n_boxes-1,iz+1));
  

  /* Some debug */
  // printf("%d:%d:%d %d:%d:%d %d:%d:%d \n",ixmin,ix,ixmax,iymin,iy,iymax,izmin,iz,izmax);
  
  double norm = 0.;
  double result = 0.;
  double density = 0.;
  double weight;
  long int num_points = 0;
  int status;

  double r[4];
  r[1] = x;
  r[2] = y;
  r[3] = z;


loop:
  result = 0;
  density = 0;
  num_points = 0;

  /* Consider only the boxes close to  */
  for (ix = ixmin; ix<=ixmax; ix++) {
    for(iy = iymin; iy<=iymax; iy++) {
      for(iz = izmin; iz<=izmax; iz++) {
        for(i = 0; i < grid[ix][iy][iz]; i++) {

          /* This part will compute the true density around the point. It will give better accuracy, but slow the code down a lot.*/
          // double ix2, iy2, iz2, i2;
          // for(ix2 = ixmin; ix2<=ixmax;ix2++){for(iy2 = iymin; iy2<=iymax;iy2++){for(iz2 = izmin; iz2<=izmax;iz2++){for(i2 = 0; i2< grid[ix2][iy2][iz2];i2++){ 
          //   status = gsl_sf_exp_e(-pow(distance(vals[ix][iy][iz][i],vals[ix2][iy2][iz2][i2]),2)/pow(group_length,2),&jresult);
          //   if (status == GSL_ERANGE || jresult.val != jresult.val)
          //     density = 0.0; 
          //   else
          //     density = jresult.val;
          //   
          //   vals[ix][iy][iz][i][4] += density;
          // }}}}


          num_points++;
          
          double value = mesh[ix][iy][iz][i][0];
          double dist = distance(mesh[ix][iy][iz][i], r);
          double density = mesh[ix][iy][iz][i][4];
       
          /* We weight the contribution from each point with a 1/distance law, so that the closer the point
          the stronger the influence on (x,y,z). We also include a 1/density factor so that the contribution
          from a point in a high-density regions is penalised. The objective is to make these clusters of points
          to count as one. */
          weight =
            (1./(dist+link_length/10000000.))/density
            * (1.00000001 - erf((dist-link_length)/(link_length*soft_coeff)));
    
          result += weight * value;
          norm += weight;

          /* Some debug */
          // printf("norm:%f, density:%f, val:%f, dist:%f, val@dist:%g\n",
          //   weight, mesh[ix][iy][iz][i][4], result/norm, dist, value);
  
        }
      }
    }
  }

  /* If you have a irregular mesh, some points might be very isolated. In this case the region of
  interpolation needs to be increased */
  if (num_points < 1) {
    // printf("Increasing Range!!! \n");
    #pragma omp atomic
    pw->count_range_extensions++;
    ixmin = (int) MAX(0,ixmin-1);
    iymin = (int) MAX(0,iymin-1);
    izmin = (int) MAX(0,izmin-1);
    ixmax = (int) MIN(n_boxes-1,ixmax+1);
    iymax = (int) MIN(n_boxes-1,iymax+1);
    izmax = (int) MIN(n_boxes-1,izmax+1);
    goto loop; 
  }

  /* Update counter */
  #pragma omp atomic
  pw->count_interpolations++;

  /* Return value */
  *interpolated_value = result/norm;
  
  return _SUCCESS_;
    
}


int mesh_free (
    struct mesh_interpolation_workspace * pw
    )
{
  
  int i, j, k, m;
  
  if ((pw->n_boxes > 0) && (pw->n_points > 0)) {

    /* Free mesh */
    for(i = 0; i<pw->n_boxes;i++) {
      for(j = 0; j<pw->n_boxes; j++) {
        for(k = 0; k<pw->n_boxes; k++) {
          for(m = 0; m<pw->grid[i][j][k]; m++)
            free (pw->mesh[i][j][k][m]);
  
          free (pw->mesh[i][j][k]);
        }
        free (pw->mesh[i][j]);
      }
      free (pw->mesh[i]);
    }
    free (pw->mesh);

  
    /* Free grid */
    if (pw->compute_grid==_TRUE_) {
      for(i = 0; i<pw->n_boxes; i++) {
        for(j = 0; j<pw->n_boxes; j++)
          free (pw->grid[i][j]);
        free (pw->grid[i]);
      }
      free (pw->grid);
    }

  }

  free (pw);

  return _SUCCESS_;
  
}





double distance (double * vec1, double * vec2){
  
  return  sqrt ( (vec1[1]-vec2[1])*(vec1[1]-vec2[1])
                +(vec1[2]-vec2[2])*(vec1[2]-vec2[2])
                +(vec1[3]-vec2[3])*(vec1[3]-vec2[3]) );
}





