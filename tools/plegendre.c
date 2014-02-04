#include "stdio.h"
#include "math.h"
#define _PI_ 3.1415926535897932384626433832795e0

// Compute Legendre polynomials.
// Taken from Numerical Recipes, Third edition, pag.294, Press et al. 2002.
double plegendre_lm(int l, int m, double x);


/// ========================
/// = Legendre polynomials =
/// ========================
// From Numerical Recipes, Third edition, pag.294, Press et al. 2002.
// The function returns the associated Legendre polynomial using a recurrence relation algorithm.
// The returned polynomials include the same normalisation constant found in the definition
// of the spherical harmonics, that is sqrt( (2l+1)/(4pi) * (l-m)!/(l+m)! )
double plegendre_lm(int l, int m, double x) {
  
  int i,ll;
  double fact,oldfact,pll,pmm,pmmp1,omx2;

  if (m < 0 || m > l || abs(x) > 1.0) {
    printf("ERROR: Bad arguments in routine plegendre_lm");
 }

  pmm=1.0;
  // Compute P_m^m
  if (m > 0) {
    omx2=(1.0-x)*(1.0+x);
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= omx2*fact/(fact+1.0);
      fact += 2.0;
    }
  }
  pmm=sqrt((2*m+1)*pmm/(4.0*_PI_));
  if (m & 1)
    pmm=-pmm;
  if (l == m)
    return pmm;
  else {
    // Compute P_m+1^m
    pmmp1=x*sqrt(2.0*m+3.0)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      // Compute P_l^m, l>m+1
      oldfact=sqrt(2.0*m+3.0);
      for (ll=m+2;ll<=l;ll++) {
        fact=sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
        pll=(x*pmmp1-pmm/oldfact)*fact;
        oldfact=fact;
        pmm=pmmp1;
        pmmp1=pll; 
      }
      return pll; 
    }
  }
}
