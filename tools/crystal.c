#include "crystal.h"


double distance (double * vec1, double * vec2){
  
  return  sqrt ( (vec1[1]-vec2[1])*(vec1[1]-vec2[1])
                +(vec1[2]-vec2[2])*(vec1[2]-vec2[2])
                +(vec1[3]-vec2[3])*(vec1[3]-vec2[3]) );
}



int mesh_sort (
    int type,
    struct bispectra * pbi,
    int index_bt,
    double ** vals
    )
{

  /*This routine presorts an unsorted mesh. 
  In: ** vals: vals[n][0] is the value of the nth point, vals[n][1] is the x coordinate, [2] is the y, [3] is z. vals[n][4] is just a placeholder and can be initialized to anything  
  num_points: the number of points
  link_length: the local region influencing the values, on an in homogenous grid this should correspond roughly to the largest distance of two neighbouring points
  soft_coeff: softens the link_length. The default value on most meshs should be 0.5
  group_length: down weights clusters of many close points. In an in homogenous grid this should correspond to the shortest distance of two points
  max_l: the mesh assumes the arguments in each direction should be between 0 and max_l
  Out: grid: this specifies how many points fall in each bin of the new sorted mesh. It will be allocated by this routine automatically
  Return value: the new sorted mesh 
   
  */

  /* Some debug */
  // int index_conf;
  // for (index_conf=0; index_conf < num_points; ++index_conf) {
  //   int i;
  //   for (i=0; i < 5; ++i) {
  //     printf("%g ", vals[index_conf][i]);
  //   }
  //   printf("\n");
  // }

  /* Read values from structure */
  int num_points = pbi->n_total_configurations;
  double max_l = (double)(pbi->l[pbi->l_size-1]);
  double link_length;
  double soft_coeff;
  double group_length;
  double ***** mesh;
  int *** grid;

  /* Fine mesh */
  if (type == 1) {
    link_length = pbi->mesh_1_link_length;
    soft_coeff = pbi->mesh_1_soft_coeff;
    group_length = pbi->mesh_1_group_length;
    mesh = pbi->bispectrum_mesh_1[index_bt];
    grid = pbi->bispectrum_grid_1[index_bt];
  }
  else if (type == 2) {
    link_length = pbi->mesh_2_link_length;
    soft_coeff = pbi->mesh_2_soft_coeff;
    group_length = pbi->mesh_2_group_length;
    mesh = pbi->bispectrum_mesh_2[index_bt];
    grid = pbi->bispectrum_grid_2[index_bt];
  }
  else
    printf ("ERROR\n");

  int i,j,k,m,n,ix,iy,iz,status;
  int *** counter;
  double safe = 1.; /*this can be increased to improve the accuracy, but 1 is already very good if link_length is chosen properly*/
  double density;
  int ndata = ceil(max_l/ (safe*link_length*(1.+soft_coeff)));

  grid = (int***) malloc(ndata*sizeof(int**));

  counter = (int***) malloc(ndata*sizeof(int**));
  for(i=0;i<ndata;i++){
    (grid)[i] = (int**) malloc(ndata*sizeof(int*));
    counter[i] = (int**) malloc(ndata*sizeof(int*));
    for(j=0;j<ndata;j++){
      (grid)[i][j] = (int *) malloc(ndata*sizeof(int));
      counter[i][j] = (int *) malloc(ndata*sizeof(int));
    }
  }
  
  for(i=0;i<ndata;i++){
    for(j=0;j<ndata;j++){
      for(k=0;k<ndata;k++){
        (grid)[i][j][k] = 0;
        counter[i][j][k] = 0;
      }
    }
  }
  for (i =0;i<num_points;i++){
    ix = floor(vals[i][1]/(safe*link_length*(1.+soft_coeff)));
    iy = floor(vals[i][2]/(safe*link_length*(1.+soft_coeff)));
    iz = floor(vals[i][3]/(safe*link_length*(1.+soft_coeff)));
    (grid)[ix][iy][iz]++;
  }
  // double ***** sorted_vals;
  mesh = (double*****) malloc(ndata*sizeof(double****));
  for(i=0;i<ndata;i++){
    mesh[i] = (double****) malloc(ndata*sizeof(double***));
    for(j=0;j<ndata;j++){
      mesh[i][j] = (double ***) malloc(ndata*sizeof(double**));
      for(k=0;k<ndata;k++){
        /*printf("I %d J %d K %d SIZE %d \n",i,j,k,(grid)[i][j][k]);*/
        mesh[i][j][k] = (double **) malloc(((grid)[i][j][k])*sizeof(double*));
        
        for(m=0;m<(grid)[i][j][k];m++){
          mesh[i][j][k][m] = (double *) malloc(5*sizeof(double)); 
        }
      }
    }
  }
  for(i = 0;i<num_points;i++){
    ix = floor(vals[i][1]/(safe*link_length*(1.+soft_coeff)));
    iy = floor(vals[i][2]/(safe*link_length*(1.+soft_coeff)));
    iz = floor(vals[i][3]/(safe*link_length*(1.+soft_coeff)));

    int Q;
    for (Q=0; Q < 5; ++Q)
      mesh[ix][iy][iz][counter[ix][iy][iz]][Q] = vals[i][Q];

    counter[ix][iy][iz]++;
  }
  
  gsl_sf_result jresult;
  gsl_error_handler_t *old_handler;
  old_handler = gsl_set_error_handler_off();

  /* The following part of the code computes the local density of points and stores them in the [4]-entry of each point */
  
  for(i = 0; i<ndata;i++) {
    for(j = 0; j<ndata; j++) {
      for(k = 0; k<ndata; k++) {
        for(m = 0; m< (grid)[i][j][k]; m++){
          for(n = 0; n< (grid)[i][j][k];n++){ 

            status = gsl_sf_exp_e(-pow(distance(mesh[i][j][k][m],mesh[i][j][k][n]),2)/pow(group_length,2),&jresult);
            if (status == GSL_ERANGE || jresult.val != jresult.val)
              density = 0.0;
            else
              density = jresult.val;
      
            mesh[i][j][k][m][4] += density;
            
          }
        }
      }
    }
  }

  /* Some debug */
  // for(i=0;i<ndata;i++){
  //     for(j=0;j<ndata;j++){
  //     for(k=0;k<ndata;k++){
  //       for(m=0;m<(grid)[i][j][k];m++){
  // 
  //           int Q;
  //           for (Q=0; Q < 5; ++Q) {
  //             printf("%g ", mesh[i][j][k][m][Q]);
  //           }
  //           printf("\n");
  //         }
  //       }
  //     }
  //   }


  return _SUCCESS_;
  
}



int mesh_int (
    int type,
    struct bispectra * pbi,
    int index_bt,
    double l1,
    double l2,
    double l3,
    double *interpolated_value
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
  double max_l = (double)(pbi->l[pbi->l_size-1]);
  double link_length;
  double soft_coeff;
  double group_length;
  double ***** mesh;
  int *** grid;
  
  /* Fine mesh */
  if (type == 1) {
    link_length = pbi->mesh_1_link_length;
    soft_coeff = pbi->mesh_1_soft_coeff;
    group_length = pbi->mesh_1_group_length;
    mesh = pbi->bispectrum_mesh_1[index_bt];
    grid = pbi->bispectrum_grid_1[index_bt];
  }
  else if (type == 2) {
    link_length = pbi->mesh_2_link_length;
    soft_coeff = pbi->mesh_2_soft_coeff;
    group_length = pbi->mesh_2_group_length;
    mesh = pbi->bispectrum_mesh_2[index_bt];
    grid = pbi->bispectrum_grid_2[index_bt];
  }
  else
    printf ("ERROR\n");

  
  /* Check bounds */
  double min_l = mesh[0][0][0][0][1];
  class_test ((l1 > max_l) || (l2 > max_l) || (l3 > max_l) || (l1 < min_l) || (l2 < min_l) || (l3 < min_l),
              pbi->error_message,
              "interpolation out of bounds: l1,l2,l3 should all be within [%g,%g] (you gave %g,%g,%g)", min_l, max_l, l1, l2, l3);
  
  int ix,iy,iz,i,ix2,iy2,iz2,i2;
  double safe = 1.;
  int ndata = ceil(max_l/ (safe*link_length*(1.+soft_coeff)));
  
  ix = floor(l1/(safe*link_length*(1.+soft_coeff)));
  iy = floor(l2/(safe*link_length*(1.+soft_coeff)));
  iz = floor(l3/(safe*link_length*(1.+soft_coeff)));

  int ixmin = (int) fmax(0,ix-1);
  int iymin = (int) fmax(0,iy-1);
  int izmin = (int) fmax(0,iz-1);
  int ixmax = (int) fmin(ndata-1,ix+1);
  int iymax = (int) fmin(ndata-1,iy+1);
  int izmax = (int) fmin(ndata-1,iz+1);
  
  double x[5];
  x[1] = l1;
  x[2] = l2;
  x[3] = l3;

  /* Some debug */
  // printf("%d:%d:%d %d:%d:%d %d:%d:%d \n",ixmin,ix,ixmax,iymin,iy,iymax,izmin,iz,izmax);
  
  double norm = 0.;
  double result = 0.;
  double density = 0.;
  double lnorm;
  int num_points = 0;
  int status;
  gsl_sf_result jresult;
  gsl_error_handler_t *old_handler;
  old_handler = gsl_set_error_handler_off();


loop:
  result = 0;
  density = 0;
  num_points = 0;
  for (ix = ixmin; ix<=ixmax; ix++) {
    for(iy = iymin; iy<=iymax; iy++) {
      for(iz = izmin; iz<=izmax; iz++) {
        for(i = 0; i< grid[ix][iy][iz]; i++) {

          /* This part will compute the true density around the point. It will give better accuracy, but slow the code down a lot.*/
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
          double dist = distance(mesh[ix][iy][iz][i], x);
          double value = mesh[ix][iy][iz][i][0];
       
          lnorm =   (1./(dist+link_length/10000000))/mesh[ix][iy][iz][i][4]
              * (1.00000001 - gsl_sf_erf( 
                  (dist-link_length)/(link_length*soft_coeff) 
              ));
    
          result += lnorm * value;
          norm += lnorm;

          /* Some debug */
          // printf("norm:%f, density:%f, val:%f, dist:%f, val@dist:%g\n",
          //   lnorm, mesh[ix][iy][iz][i][4], result/norm, dist, value);
  
        }
      }
    }
  }

  /* If you have a irregular mesh, some points might be very isolated. In this case the region of interpolation needs
    to be increased */
  if (num_points < 1) {
    printf("Increasing Range!!! \n");
    ixmin = (int) fmax(0,ixmin-1);
    iymin = (int) fmax(0,iymin-1);
    izmin = (int) fmax(0,izmin-1);
    ixmax = (int) fmin(ndata-1,ixmax+1);
    iymax = (int) fmin(ndata-1,iymax+1);
    izmax = (int) fmin(ndata-1,izmax+1);
    goto loop; 
  }

  *interpolated_value = result/norm;
  
  return _SUCCESS_;
    
}
