#include "common.h"

int main (int argc, char const *argv[]) {

  int l, m1, m;
  int l_size = 5;

  for(l=0; l<=l_size; ++l) {
    for(m=-l; m<=l; ++m) {
      

      m1 = m;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", c_minus(l,m1,m));
      printf("c_plus  = %4.4e", c_plus(l,m1,m));      

      printf("\n");

      m1 = m+1;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", c_minus(l,m1,m));
      printf("c_plus  = %4.4e", c_plus(l,m1,m));      

      printf("\n");

      m1 = m-1;
      printf("(%d,%d,%d) : ", l, m1, m);
      printf("c_minus = %4.4e,   ", c_minus(l,m1,m));
      printf("c_plus  = %4.4e", c_plus(l,m1,m));      

      printf("\n\n");

    }
  }



}