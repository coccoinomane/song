/** @file print_couplings.c
 *
 * Print the coupling factors used in SONG to build the Boltzmann
 * hierarchies.
 *
 * The coupling factors are defined in detail in section A.4.1 
 * of http://arxiv.org/abs/1405.2280.
 *
 * Created by Guido W. Pettinari on 15.06.2011
 * Last edited by Guido W. Pettinari on 30.06.2015
 */

#include "song.h"

int main (int argc, char const *argv[]) {

  int l, m1, m;
  int l_size = 5;

  for(l=0; l<=l_size; ++l) {
    for(m=-l; m<=l; ++m) {
      

      m1 = m;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", coupling_c_minus(l,m1,m));
      printf("c_plus  = %4.4e", coupling_c_plus(l,m1,m));      

      printf("\n");

      m1 = m+1;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", coupling_c_minus(l,m1,m));
      printf("c_plus  = %4.4e", coupling_c_plus(l,m1,m));      

      printf("\n");

      m1 = m-1;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", coupling_c_minus(l,m1,m));
      printf("c_plus  = %4.4e", coupling_c_plus(l,m1,m));      

      printf("\n\n");

    }
  }

  return 0;
}