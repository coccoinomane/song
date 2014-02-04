#include "song.h"

int main (int argc, char const *argv[])
{
  
  if (argc<3) {
    printf ("%s needs three arguments: l,m,x\n", argv[0]);
    return _FAILURE_;
  }
  
  int l = atoi (argv[1]);
  int m = atoi (argv[2]);
  double x = atof (argv[3]);
    
  double plm;
  double plm_rescaled_analytically;
  double plm_rescaled;
  double plm_ratio;
  
  if (m >= 0) {  
    plm = plegendre_lm (l,m,x);
    plm_rescaled_analytically = plegendre_lm_rescaled_analytically (l,m,x);
    plm_rescaled = plegendre_lm_rescaled (l,m,x);
    plm_ratio = plm/pow(1.0-x*x, m/2.0);
  }
  else {
    plm = alternating_sign(m) * plegendre_lm (l,abs(m),x);
    plm_rescaled_analytically = plm * pow(1.0-x*x, abs(m)/2.0);
    plm_rescaled = plegendre_lm_rescaled (l,m,x);
    plm_ratio = plm * pow(1.0-x*x, abs(m)/2.0);
  }
  
  printf ("l=%d,m=%d,x=%.17g\n", l, m, x);
  printf ("plm                       = %.17g\n", plm);
  printf ("plm_rescaled_analytically = %.17g\n", plm_rescaled_analytically);
  printf ("plm_rescaled              = %.17g\n", plm_rescaled);
  printf ("plm_ratio                 = %.17g\n", plm_ratio);
  
  return 0;
}