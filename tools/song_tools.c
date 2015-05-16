/** @file song_tools.c 
 * Created by Guido Walter Pettinari on 12.08.2011
 * Last edited by Guido Walter Pettinari on 1.09.2011
 *
 * Support function for the perturbations2 module
 */

#include "song_tools.h"
#include "slatec_3j_C.h"
#include "math.h"


// ==============================================================================
// =                                3J and 6J symbols                           =
// ==============================================================================

/** 
 * Compute the 3j symbol
 * (    l1     l2     l3   )
 * ( -m2-m3    m2     m3   )
 * 
 * The result shall be stored into 'threej'.
 *
 */
int threej_single(
       int l1, int l2, int l3, int m2, int m3, // In
       double *threej,                         // Out
       ErrorMsg errmsg       
       )
{

  /* Limits from the triangular condition */
  int l1_min = abs(l2-l3);
  int l1_max = l2+l3;

	/* Check input */
	class_test (l1>l1_max || l1<l1_min,
    errmsg,
		"ERROR, %s: 'l1' should be between %d and %d.\n", __func__, l1_min, l1_max);

  /* The Slatec function DRC3JJ computes the 3j symbol for any allowed value
  of l1.  Here we compute the number of such values and allocate memory
  for the result array. */
  int out_size = l1_max - l1_min + 1;
  double *result = (double*)calloc(out_size,sizeof(double));

  /* Do the actual computation */
  double l1_min_D;
  double l1_max_D;
  
  class_call (drc3jj (
                l2, l3, m2, m3,
                &l1_min_D, &l1_max_D,
                result,
                out_size,
                errmsg       
                ),
    errmsg,
    errmsg);
    
  l1_min = (int)(l1_min_D+_EPS_);

  /* Find the index corresponding to the requested l1 */
  *threej = result[l1-l1_min];
  
  free (result);
  
  /* Print the result */
  // printf("3J(%d,%g,%g)(%g,%g,%g) = %g\n", l1,l2,l3,m1,m2,m3,result[l1-l1_min]);
  
  return _SUCCESS_;
  
}


/** 
 *  Compute the ratio between the 3j symbols
 *  (    l1+2*N1     l2+2*N2     l3+2*N3   )
 *  (          0           0           0   )
 * and
 *  (    l1     l2     l3   )
 *  (     0      0      0   )  
 * using the recursive relation in Schulten & Gordon, 1961,
 * http://scitation.aip.org/content/aip/journal/jmp/16/10/10.1063/1.522426?ver=pdfcov.
 *
 * The function can be easily generalised to arbitrary azimuthal numbers m1,m2,m3
 * by making use of the threej_B function from the same paper and coded below.
 * 
 * The relation is valid only for even values of l1+l2+l3, but the function will return
 * also when this is not true.
 * 
 * TO DO:
 *
 * 1) Allow odd increments. I want to be able to compute stuff like 3J(l1+3 l2 l2+3)/3J(l1 l2 l2).
 * 
 * 2) Allow for l1+l2+l3 to be odd.
 *
 */
int threej_ratio_L (
      int L1, int L2, int L3,            // In
      int N1, int N2, int N3,            // In
      double *result,                    // Out, should be allocated with M+1 elements
      ErrorMsg errmsg
      )
{
  
  /* Debug - print arguments */
  // printf ("{L1, L2, L3, N1, N2, N3} = {%d,%d,%d,%d,%d,%d};\n",
  //   L1, L2, L3, N1, N2, N3);

  class_test (!is_triangular_int (L1,L2,L3),
    errmsg,
    "arguments violate triangular condition, (L1,L2,L3)=(%d,%d,%d)\n",
    L1, L2, L3);
    
  class_test (((L1+2*N1)<0) || ((L2+2*N2)<0) || ((L3+2*N3)<0),
    errmsg,
    "arguments give rise to negative multipoles");

  double A_factor;
    
  /* First element in the recursion is 3J[0,0,0]/3J[0,0,0]=1 */
  *result = 1;

  /* Initialise arrays of N and L values */
  int L[] = {0, L1, L2, L3};
  int n[] = {0, N1, N2, N3};
  int l[] = {0, L1+2*N1, L2+2*N2, L3+2*N3};
    
  /* We obtain 3J[L1+2*n1,L2+2*n2,L3+2*n3] by using the recurrence relation in Schulten & Gordon, 1961.
  The starting point is 3J[L1,L2,L3], which is multiplied by products of the A function to obtain 
  3J[L1+2*n1,L2+2*n2,L3+2*n3]. The number of products is equal to to n1+n2+n3. The order of the
  recurrences is not unique; for example, 3J[L1+2,L2+2,L3] can be obtained either as
  A(L1,L2+2,L3)*A(L2,L3,L1)*3J[L1,L2,L3] or A(L1,L2,L3)*A(L2,L3,L1+2)*3J[L1,L2,L3]. We pick the
  ordering such that the triangular condition on the arguments of A is always satisfied. */

  /* Cycle until all counters reach the zero, corresponding to 3J[L1,L2,L3] */
  while ((n[1]!=0) || (n[2]!=0) || (n[3]!=0)) {
 
    /* Debug - print the current N1,N2,N3 */
    // printf ("{n1,n2,n3}={%d,%d,%d}\n", n[1], n[2], n[3]);
    
    /* Initialise temporary vectors */
    int n_nonzero = 0;
    int non_zero_positions[] = {0, 0, 0, 0};
    int l_next[] = {0, l[1], l[2], l[3]};

    /* Find which columns have to be incremented/decremented. We must find at
    least one (n_nonzero>=1) because of the 'while' condition. */
    for (int i=1; i <= 3; ++i) {
      if (n[i] != 0) {
        n_nonzero++;
        non_zero_positions[n_nonzero] = i;
      }
    }
    
    /* By default, we try to increment/decrement the first non-zero column first */
    int i = 1;
    int pos = non_zero_positions[i];
    l_next[pos] -= 2*SIGN(n[pos]);

    /* If incrementing/decrementing the chosen column would violate the triangular
    condition, try with the next non-zero column. */
    while (!is_triangular_int (l_next[1], l_next[2], l_next[3])) {
      l_next[pos] += 2*SIGN(n[pos]); /* Undo previous attempt */
      pos = non_zero_positions[++i]; /* Increment position */
      class_test (i>n_nonzero, errmsg, "could not find valid configurations, bug?");
      l_next[pos] -= 2*SIGN(n[pos]); /* Try with the next non-zero column */
    }

    /* Compute the two positions complementary to 'pos' using the modulo function */
    int pos2 = (pos+1)%3;
    if (pos2==0) pos2=3;
    int pos3 = (pos+2)%3;
    if (pos3==0) pos3=3;

    /* Compute the ratio between 3J[l] and 3J[l_next] */
    class_call (threej_A_factor (L[pos], l[pos2], l[pos3], n[pos], &A_factor, errmsg), errmsg, errmsg);
    *result *= A_factor;
    n[pos] = SIGN(n[pos])*(abs(n[pos])-1);
    l[pos] = l_next[pos];
    
  }; // end of while
  
  class_test (isnan (*result), 
    errmsg,
    "encountered nan for {L1, L2, L3, N1, N2, N3} = {%d,%d,%d,%d,%d,%d};\n",
    L1, L2, L3, N1, N2, N3);

  /* Debug - print result */
  // printf ("*result = %g\n", *result);
    
  /* Debug - Mathematica batch output */
  // fprintf (stderr, "Abs[1-directRatioL[%d.,%d.,%d.,%d.,%d.,%d.]/((%d)*10^(%f))],\n",
  //   L1, L2, L3, N1, N2, N3, SIGN(*result), log10(fabs(*result)));

  return _SUCCESS_;

}

/** 
 *  Compute the ratio between the 3j symbols
 *  (    l1+2*N1    l2     l3   )
 *  (     0          0      0   )
 * and
 *  (    l1     l2     l3   )
 *  (     0      0      0   )  
 * for all values of 'n' in [0,abs(N)], using the recursive relation in  Schulten & Gordon, 1961
 * http://scitation.aip.org/content/aip/journal/jmp/16/10/10.1063/1.522426?ver=pdfcov.
 *
 * The function can be easily generalised to arbitrary azimuthal numbers m1,m2,m3
 * by making use of the threej_B function (ibidem).
 * 
 * The ratio only exists for even values of l1+l2+l3, but the function will return
 * also when this is not true. The returned value can be thought as continuously
 * interpolating the ratio for odd values of l1+l2+l3.
 *
 * Caution: this function was never tested!
 *  
 */
int threej_ratio_L_recursive (
      int l1, int l2, int l3, int N1,     // In
      double *result,                     // Out, should be allocated with M+1 elements
      ErrorMsg errmsg
      )
{

  class_test (!is_triangular_int (l1,l2,l3),
    errmsg,
    "arguments violate triangular condition, (l1,l2,l3)=(%d,%d,%d)\n",
    l1, l2, l3);

    
  class_test ((l1+2*N1)<0,
    errmsg,
    "arguments give rise to negative multipoles");

  /* First element in the recursion is 3J[0,0,0]/3J[0,0,0]=1 */
  result[0] = 1.;

  /* Recursive relation for the element 'n' */
  for (int n=1; n <= abs(N1); ++n) {
    double A_factor;
    class_call (threej_A_factor (l1, l2, l3, SIGN(N1)*n, &A_factor, errmsg), errmsg, errmsg);
    result[n] = result[n-1] * A_factor;
  }

  /* Debug - show results */
  // for (int n=0; n < N1+1; ++n) {
  //   printf ("result[%d] = %g\n", n, result[n]);
  // }

  return _SUCCESS_;

}

/**
 *  Compute the ratio between the 3j symbols
 *  (    l1+2*N1     l2     l3   )
 *  (     0           0      0   )
 * and
 *  (    l1     l2     l3   )
 *  (     0      0      0   )  .
 *
 * This function calls 'threej_ratio_L_recursive' and only outputs the last element of the
 * result array. See comments in 'threej_ratio_L_recursive' for more detail.
 *
 */
int threej_ratio_L1 (
      int l1, int l2, int l3, int N1,     // In
      double *result,                    // Out
      ErrorMsg errmsg
      )
{

  double ratio[N1+1];

  class_call (threej_ratio_L_recursive (l1,l2,l3,N1,&(ratio[0]),errmsg),
    errmsg, errmsg);

  *result = ratio[N1];

  return _SUCCESS_;

}


/** 
 *  Compute the ratio between the 3j symbols
 *  (    l1     l2     l3   )
 *  (     0      m     -m   )
 * and
 *  (    l1     l2     l3   )
 *  (     0      0      0   )  
 * for the values of 'm' in [0,M], using the recursive relation in  Schulten & Gordon, 1961
 * http://scitation.aip.org/content/aip/journal/jmp/16/10/10.1063/1.522426?ver=pdfcov.
 *
 * The ratio only exists for even values of l1+l2+l3, but the function will return
 * also when this is not true. The returned value can be thought as continuously
 * interpolating the ratio for odd values of l1+l2+l3.
 *  
 */
int threej_ratio_M_recursive (
      int l1, int l2, int l3, int M,     // In
      double *result,                    // Out, should be allocated with M+1 elements
      ErrorMsg errmsg
      )
{

  class_test ((M>l2) || (M>l3) || (l1>(l2+l3)) || (l1<abs(l2-l3)),
    errmsg,
    "the arguments violate one or more 3J conditions");

  /* First element in the recursion is 3J[0,0,0]/3J[0,0,0]=1 */
  result[0] = 1;

  if (M > 0) {
  
    /* Second element in the recurision is 3J[0,1,-1]/3J[0,0,0] */
    double C_0 = sqrt(l2*l3*(l2+1.)*(l3+1.));
    double D_0 = l2*(l2+1.) + l3*(l3+1.) - l1*(l1+1.);
    double C_1 = C_0;
    result[1] = -D_0/(C_0+C_1);

    /* Recursive relation for the element 'm' */
    for (int m=2; m <= M; ++m) {
      
      double C, C_minus1, D_minus1;
      
      class_call (threej_C(l1,l2,l3,m,-m,&C,errmsg),errmsg,errmsg);
      class_call (threej_C(l1,l2,l3,m-1,-m+1,&C_minus1,errmsg),errmsg,errmsg);
      class_call (threej_D(l1,l2,l3,m-1,-m+1,&D_minus1,errmsg),errmsg,errmsg);
      
      result[m] = - (result[m-2]*C_minus1 + result[m-1]*D_minus1) / C;
    }
    
  } // end of if(M>0)

  /* Debug - show results */
  // for (int m=0; m < M+1; ++m) {
  //   printf ("result[%d] = %g\n", m, result[m]);
  // }

  return _SUCCESS_;

}

/** 
 *  Compute the ratio between the 3j symbols
 *  (    l1     l2     l3   )
 *  (     0      M     -M   )
 * and
 *  (    l1     l2     l3   )
 *  (     0      0      0   ) .
 *
 * This function calls 'threej_ratio_M_recursive' and only outputs the last element of
 * the result array. See comments in 'threej_ratio_M_recursive' for more detail.
 *
 */
int threej_ratio_M (
      int l1, int l2, int l3, int M,     // In
      double *result,                    // Out
      ErrorMsg errmsg
      )
{

  double ratio[M+1];
  
  class_call (threej_ratio_M_recursive (l1,l2,l3,M,&(ratio[0]),errmsg),
    errmsg, errmsg);
    
  *result = ratio[M];

  return _SUCCESS_;

}



/** 
 * Support function for 'threej_ratio_L_recursive'. The first argument
 * (l1) is the multipole with a zero in the lower row.
 */
int threej_A (
      int l1, int l2, int l3,
      int m1,
      double *result,
      ErrorMsg errmsg
      )
{

  class_test (!is_triangular_int (l1,l2,l3),
    errmsg,
    "arguments violate triangular condition, (l1,l2,l3)=(%d,%d,%d)\n",
    l1, l2, l3);
    
  class_test (m1>l1,
    errmsg,
    "m1>l1 is not a valid 3J combination, (l1,m1,l2,l3)=(%d,%d,%d,%d)\n",
    l1, m1, l2, l3);

  double l1_squared = l1*l1;
  double l2_minus_l3 = l2-l3;
  double l2_plus_l3_plus_1 = l2+l3+1;
  
  *result = sqrt((l1_squared-l2_minus_l3*l2_minus_l3)
                *(l2_plus_l3_plus_1*l2_plus_l3_plus_1-l1_squared)
                *(l1_squared-m1*m1));
  
  return _SUCCESS_;

}

/** 
 * Support function for 'threej_ratio_L_recursive'. The first argument
 * (l1) is the multipole with a zero in the lower row.
 */
int threej_B (
      int l1, int l2, int l3,
      int m1, int m2, int m3,
      double *result,
      ErrorMsg errmsg
      )
{

  class_test ((m1>l1) || (m2>l2) || (m3>l3) || (l1>(l2+l3)) || (l1<abs(l2-l3)),
    errmsg,
    "the arguments violate one or more 3J conditions");
   
  *result = -(2*l1+1.)*(l2*(l2+1.)*m1 - l3*(l3+1.)*m1 - l1*(l1+1.)*(m3-m2));
  
  return _SUCCESS_;

}

/** 
 * Support function for 'threej_ratio_M_recursive'. The first argument
 * (l1) is the one for which we are computing the recursion.
 */
int threej_C (
      int l1, int l2, int l3,
      int m2, int m3,
      double *result,
      ErrorMsg errmsg
      )
{
   
  *result = sqrt((l2-m2+1.)*(l2+m2)*(l3+m3+1.)*(l3-m3));

  return _SUCCESS_;

}

/** 
 * Support function for 'threej_ratio_M_recursive'. The first argument
 * (l1) is the one for which we are computing the recursion.
 */
int threej_D (
      int l1, int l2, int l3,
      int m2, int m3,
      double *result,
      ErrorMsg errmsg
      )
{
   
  *result = l2*(l2+1.) + l3*(l3+1.) - l1*(l1+1.) + 2.*m2*m3;

  return _SUCCESS_;
  
}

/** 
 * Support function for 'threej_ratio_L_recursive'. The first argument
 * (l1) is the one for which we are computing the recursion.
 * The values of l1,l2,l3 do not need to satisfy a triangular condition.
 */
int threej_A_factor (
      int l1, int l2, int l3,
      int n1,
      double *result,
      ErrorMsg errmsg
      )
{

  // printf ("%s: {l1,l2,l3,n1}={%d,%d,%d,%d}\n",
  //   __func__, l1, l2, l3, n1);

  double A;

  /* Forward recursion */
  if (n1 > 0) {

    int l1_plus_2n = l1 + 2*n1;
    int l1_plus_2n_minus_1 = l1 + 2*n1 - 1;

    class_call (threej_A (l1_plus_2n_minus_1, l2, l3, 0, &A, errmsg), errmsg, errmsg);
    double num = -l1_plus_2n * A;
    
    class_call (threej_A (l1_plus_2n, l2, l3, 0, &A, errmsg), errmsg, errmsg);
    double den = l1_plus_2n_minus_1 * A;

    class_test (fabs(den) < _MINUSCULE_,
      errmsg,
      "num=%g, den=%g, num/den=%g, risk of nans for (l1,l2,l3)=(%d,%d,%d)\n",
      num, den, num/den, l1, l2, l3);
    
    *result = num/den;
   
  }
  /* Backward recursion */
  else if (n1 < 0) {
    
    int l1_plus_2n_plus_1 = l1 + 2*n1 + 1;
    int l1_plus_2n_plus_2 = l1 + 2*n1 + 2;

    class_call (threej_A (l1_plus_2n_plus_2, l2, l3, 0, &A, errmsg), errmsg, errmsg);
    double num = -l1_plus_2n_plus_1 * A;
    
    class_call (threej_A (l1_plus_2n_plus_1, l2, l3, 0, &A, errmsg), errmsg, errmsg);
    double den = l1_plus_2n_plus_2 * A;

    class_test (fabs(den) < _MINUSCULE_,
      errmsg,
      "num=%g, den=%g, num/den=%g, risk of nans for (l1,l2,l3)=(%d,%d,%d)\n",
      num, den, num/den, l1, l2, l3);
    
    *result = num/den;
    
  }
  /* Enforce 3J[0,0,0]=3J[0,0,0] */
  else {
    *result = 1;
  }
  
  return _SUCCESS_;
  
}


// ======================================================================
// =                           Bessel functions                         =
// ======================================================================

/**
 * Compute spherical Bessel function j_l(x) for a given l and x.
 *
 * Inspired from Numerical Recipies. This is the same as the function
 * bessel_j in the Bessel structure, but without requiring pbs as an
 * argument, so that it can be called from anywhere.
 */
double spherical_bessel_j(
       int l,
       double x
       )
{

    double nu,nu2,beta,beta2;
    double x2,sx,sx2,cx;
    double cotb,cot3b,cot6b,secb,sec2b;
    double trigarg,expterm,fl;
    double l3,cosb;
    double jl;
    
    if (l < 0) {
      printf ("ERROR, %s: argument 'l' smaller than zero\n", __func__);
      return -1;
    };

    if (x < 0) {
      printf ("ERROR, %s: argument 'x' smaller than zero\n", __func__);
      return -1;
    };
  
    fl = (double)l;
  
    x2 = x*x;
  
    /************* Use closed form for l<7 **********/
  
    if (l < 7) {
  
      sx=sin(x);
      cx=cos(x);
  
      if(l == 0) {
        if (x > 0.1) jl=sx/x;
        else jl=1.-x2/6.*(1.-x2/20.);
        return jl;
      }
  
      if(l == 1) {
        if (x > 0.2) jl=(sx/x -cx)/x;
        else jl=x/3.*(1.-x2/10.*(1.-x2/28.));
        return jl;
      }
  
      if (l == 2) {
        if (x > 0.3) jl=(-3.*cx/x-sx*(1.-3./x2))/x;
        else jl=x2/15.*(1.-x2/14.*(1.-x2/36.));
        return jl;
      }
  
      if (l == 3) {
        if (x > 0.4) jl=(cx*(1.-15./x2)-sx*(6.-15./x2)/x)/x;
        else jl=x*x2/105.*(1.-x2/18.*(1.-x2/44.));
        return jl;
      }
  
      if (l == 4) {
        if (x > 0.6) jl=(sx*(1.-45./x2+105./x2/x2) +cx*(10.-105./x2)/x)/x;
        else jl=x2*x2/945.*(1.-x2/22.*(1.-x2/52.));
        return jl;
      }
      
      if (l == 5) {
        if (x > 1) jl=(sx*(15.-420./x2+945./x2/x2)/x -cx*(1.0-105./x2+945./x2/x2))/x;
        else jl=x2*x2*x/10395.*(1.-x2/26.*(1.-x2/60.));
        return jl;
      }
  
      if (l == 6) {
        if (x > 1) jl=(sx*(-1.+(210.-(4725.-10395./x2)/x2)/x2)+
         cx*(-21.+(1260.-10395./x2)/x2)/x)/x;
        else jl=x2*x2*x2/135135.*(1.-x2/30.*(1.-x2/68.));
        return jl;
      }
  
    }
  
    else {
  
      if (x <= 1.e-40) {
        jl=0.0;
        return jl;
      }
  
      nu= fl + 0.5;
      nu2=nu*nu;
  
      if ((x2/fl) < 0.5) {
        jl=exp(fl*log(x/nu/2.)+nu*(1-log(2.))-(1.-(1.-3.5/nu2)/nu2/30.)/12./nu)
  /nu*(1.-x2/(4.*nu+4.)*(1.-x2/(8.*nu+16.)*(1.-x2/(12.*nu+36.))));
        return jl;
      }
  
      if ((fl*fl/x) < 0.5) {
  
        beta = x - _PI_/2.*(fl+1.);
        jl = (cos(beta)*(1.-(nu2-0.25)*(nu2-2.25)/8./x2*(1.-(nu2-6.25)*(nu2-12.25)/48./x2))
       -sin(beta)*(nu2-0.25)/2./x* (1.-(nu2-2.25)*(nu2-6.25)/24./x2*(1.-(nu2-12.25)*(nu2-20.25)/80./x2)) )/x;
        
        return jl;
  
      }
  
      l3 = pow(nu,0.325);
  
      if (x < nu-1.31*l3) {
  
        cosb=nu/x;
        sx=sqrt(nu2-x2);
        cotb=nu/sx;
        secb=x/nu;
        beta=log(cosb+sx/x);
        cot3b=cotb*cotb*cotb;
        cot6b=cot3b*cot3b;
        sec2b=secb*secb;
        expterm=((2.+3.*sec2b)*cot3b/24.
         - ((4.+sec2b)*sec2b*cot6b/16.
     + ((16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b/5760.
        + (32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128./nu)*cot6b/nu)/nu)/nu;
        jl=sqrt(cotb*cosb)/(2.*nu)*exp(-nu*beta+nu/cotb-expterm);
  
        return jl;
  
      }
  
      if (x > nu+1.48*l3) {
  
        cosb=nu/x;
        sx=sqrt(x2-nu2);
        cotb=nu/sx;
        secb=x/nu;
        beta=acos(cosb);
        cot3b=cotb*cotb*cotb;
        cot6b=cot3b*cot3b;
        sec2b=secb*secb;
        trigarg=nu/cotb-nu*beta-_PI_/4.
  -((2.+3.*sec2b)*cot3b/24.
    +(16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b*cot6b/5760./nu2)/nu;
        expterm=((4.+sec2b)*sec2b*cot6b/16.
         -(32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b*cot6b/128./nu2)/nu2;
        jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);
  
        return jl;
      }
      
      /* last possible case */
  
      beta=x-nu;
      beta2=beta*beta;
      sx=6./x;
      sx2=sx*sx;
      secb=pow(sx,1./3.);
      sec2b=secb*secb;
      jl=(_GAMMA1_*secb + beta*_GAMMA2_*sec2b
   -(beta2/18.-1./45.)*beta*sx*secb*_GAMMA1_
   -((beta2-1.)*beta2/36.+1./420.)*sx*sec2b*_GAMMA2_
   +(((beta2/1620.-7./3240.)*beta2+1./648.)*beta2-1./8100.)*sx2*secb*_GAMMA1_
   +(((beta2/4536.-1./810.)*beta2+19./11340.)*beta2-13./28350.)*beta*sx2*sec2b*_GAMMA2_
   -((((beta2/349920.-1./29160.)*beta2+71./583200.)*beta2-121./874800.)*
     beta2+7939./224532000.)*beta*sx2*sx*secb*_GAMMA1_)*sqrt(sx)/12./sqrt(_PI_);
  
      // *** MY MODIFICATIONS
      // printf("j_l(%d,%20.15g) = %15g\n", l, x, jl);
  
      return jl;
  
    }
  
    printf ("ERROR, %s: value of l=%d or x=%e out of bounds\n", __func__,l,x);
    return -1;

}

// ====================================================================================
// =                                Coupling factors                                  =
// ====================================================================================

/**
 * 'C' and 'D' coupling coefficients that appear in the in the Boltzmann hierarchy for the photon
 * temperature (see eqs. A.11 and 2.18 of Beneke & Fidler 2011). They are basically
 * compact forms of 3j symbols that naturally arise when peforming the multipole expansion
 * of Boltzmann equation.
 */
double coupling_c_minus (int l, int m1, int m) {

  if ((abs(m1-m)>1) && (abs(m)>l))
    printf("WARNING: 'c_minus' called with wrong arguments: l=%d, m1=%d, m=%d.\n", l, m1, m);

  double result;

  if (m1 == m+1)
    if (l==m)
      result = 0;
    else
      result = sqrt((l-1.-m)*(l-m)) / (sqrt_2*(2*l-1.));

  if (m1 == m-1)
    if (l==-m)
      result = 0;
    else
      result = sqrt((l-1.+m)*(l+m))/(sqrt_2*(2*l-1.));

  if (m1 == m)
    if (m == 0)
      result = l/(2*l-1.);
    else
      result = sqrt((double)l*l-m*m)/(2*l-1.);
  
  // *** Some debug
  // printf("c_minus(%d,%d,%d)=%g\n", l, m1, m, result);
  
  return result;
      
}


double coupling_c_plus (int l, int m1, int m) {
  
  if ((abs(m1-m)>1) && (abs(m)>l))
    printf("WARNING: 'c_plus' called with wrong arguments: l=%d, m1=%d, m=%d.\n", l, m1, m);

  double result;

  if (m1 == m+1)
    result = -sqrt((l+1.+m)*(l+2+m))/(sqrt_2*(2*l+3.));

  if (m1 == m-1)
    result = -sqrt((l+1.-m)*(l+2-m))/(sqrt_2*(2*l+3.));

  if (m1 == m)
    if (m == 0)
      result = (l+1)/(2*l+3.);
    else
      result = sqrt((l+1.)*(l+1)-m*m)/(2*l+3.);
 
  return result;
      
}


double coupling_d_zero (int l, int m1, int m) {

  if ((abs(m1-m)>1) && (abs(m)>l))
    printf("WARNING: 'd_zero' called with wrong arguments: l=%d, m1=%d, m=%d.\n", l, m1, m);

  double result;

  if (m1 == m+1)
    if (l==m)
      result = 0;
    else
      result = -sqrt(2*(l+1.+m)*(l-m)) / (l*(l+1.));

  if (m1 == m-1)
    if (l==-m)
      result = 0;
    else
      result = sqrt(2*(l+1.-m)*(l+m)) / (l*(l+1.));

  if (m1 == m)
    if (m == 0)
      result = 0;
    else
      result = -2*m/(l*(l+1.));
        
  // *** Some debug
  // printf("d_zero(%d,%d,%d)=%g\n", l, m1, m, result);
  
  return result;
      
}

double coupling_d_minus (int l, int m1, int m) {
  
  return sqrt(l*l-4.)/l * coupling_c_minus(l,m1,m);
  
}

double coupling_d_plus (int l, int m1, int m) {
 
  return sqrt( (l-1.)*(l+3.) )/(l+1.) * coupling_c_plus(l,m1,m);
  
}



/**
 * Compute the Gaunt-like coupling coefficient:
 * 
 *    (-1)^m3 * (2*l3+1) * ( l1  l2  l3 ) * (  l1   l2      l3  )
 *                         ( 0    F  -F )   (  m1   m2  -m1-m2  )
 *   
 *     = ( l1  l2 | l3 ) * (  l1   l2  |    l3  )
 *       ( 0    F |  F )   (  m1   m2  | m1+m2  )
 *
 * for all allowed values of l1 and m2, with m3=-m1-m2, where (   ) denotes
 * a 3j symbol and (   | ) a Clebsch-Gordan coefficient. The result is stored in
 * 'result' as result[l1-l1_min][m2-m2_min], which should be preallocated, where
 * l1_min and m2_min are outputs of this function.
 *
 */
int coupling_general (
  int l2, int l3, int m1, int F,
  double * three_j_000, /* should be preallocated with at least l2_max doubles */
  int three_j_000_size,
  double * three_j_mmm, /* should be preallocated with at least m1_max doubles */
  int three_j_mmm_size,
  int * l1_min, int * l1_max,
  int * m2_min, int * m2_max,
  double ** result,     /* should be preallocated with at least l2_max*m1_max doubles */
  ErrorMsg errmsg 
  )
{
  
  /* Test input */
  class_test (abs(m1)>l2+l3, errmsg, "m1 is out of bounds: abs(%d)>%d+%d", m1, l2, l3);
  class_test ((l2<F) || (l3<F), errmsg, "l2 and l3 out of bounds, should be smaller than F=%d", F);
 
  /* Temporary values needed for the computation of the 3j symbol */
  double l1_min_D, l1_max_D;
  double m2_min_D, m2_max_D;

  /* Compute 
  * (  l1  l2  l3  )
  * (  0   F   -F   )
  * for all allowed values of l2 */
  class_call (drc3jj (
                l2, l3, F, -F,
                &l1_min_D, &l1_max_D,
                three_j_000,
                three_j_000_size,
                errmsg       
                ),
    errmsg,
    errmsg);
    
  *l1_min = (int)(l1_min_D+_EPS_);
  *l1_max = (int)(l1_max_D+_EPS_);
    
  /* Adjust the range of l1, so to exclude from the output those values of l1 that are
  smaller than the requested m1. Note that even after increasing l1_min, it will always be
  smaller than l1_max because, above, we have already ensured that abs(m1)<=l2+l3=l_max. */
  double * three_j_000_correct = three_j_000;
  if (abs(m1) > *l1_min) {
    three_j_000_correct += abs(m1) - *l1_min;
    *l1_min = abs(m1);
  }
  
  /* Paranoid check */
  class_test (*l1_min>*l1_max, errmsg, "paranoid android");

  /* LOOP ON L1 - made in such a way that it always respects the triangular condition */
  for (int l1=*l1_min; l1 <= *l1_max; ++l1) {

    class_test ((F==0) && ((l1+l2+l3)%2!=0) && (three_j_000_correct[l1-*l1_min]!=0),
      errmsg,
      "error in computation of 3j, odd is not zero! three_j(%d,%d,%d)=%g", l1, l2, l3, three_j_000_correct[l1-*l1_min]);

    /* Compute
    * (  l1   l2   l3       )
    * (  m1   m2   -m1-m2  )
    * for all allowed values of m2 */
    class_call (drc3jm (
                  l1, l2, l3, m1,
                  &m2_min_D, &m2_max_D,
                  three_j_mmm,
                  three_j_mmm_size,
                  errmsg       
                  ),
      errmsg,
      errmsg);

    *m2_min = (int)(m2_min_D + SIGN(m2_min_D)*_EPS_);
    *m2_max = (int)(m2_max_D + SIGN(m2_max_D)*_EPS_);

    /* LOOP ON M2 */
    for (int m2 = *m2_min; m2 <= *m2_max; ++m2) {

      int m3 = -m1-m2;

      /* Coupling coefficient for this (l1,l2,l3,m1,m2) */
      result[l1-*l1_min][m2-*m2_min] = 
        ALTERNATING_SIGN(m3) * (2*l3+1.) * three_j_000_correct[l1-*l1_min] * three_j_mmm[m2-*m2_min];
  
      /* Debug - print result */
      // if ((l1==2)&&(l2==2)&&(l3==1)&&(m1==1)&&(m3==0))
      //   printf ("*** C(l1=%d,l2=%d,l3=%d,m1=%d,m2=%d,m3=%d)=%g\n",
      //     l1, l2, l3, m1, m2, m3, result[l1-*l1_min][m2-*m2_min]);
  
    } // end of for(m)            
  } // end of for(l)
  
  return _SUCCESS_;
}


// ============================================================================================
// =                                 Legendre polynomials                                     =
// ============================================================================================


/**
 * Associate Legendre Polynomials from Numerical Recipes, Third edition, pag.294, Press et al. 2002.
 * The function returns the associated Legendre polynomial using a recurrence relation algorithm.
 * The returned polynomials include the same normalisation constant found in the definition
 * of the spherical harmonics, that is sqrt( (2l+1)/(4pi) * (l-m)!/(l+m)! )
 */
double plegendre_lm (int l, int m, double x) {
  
  int i,ll;
  double fact,oldfact,pll,pmm,pmmp1,omx2;

  if (m < 0 || m > l || fabs(x) > 1.0) {
    printf("ERROR: Bad arguments in routine plegendre_lm\n");
 }

  pmm=1.0;
  /* Compute P_m^m */
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
    /* Compute P_m+1^m */
    pmmp1=x*sqrt(2.0*m+3.0)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      /* Compute P_l^m, l>m+1 */
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


/**
 * Return the value of
 * P_lm (x) / (1-x^2)^(m/2)
 * where P_lm(x) is the associated Legendre polynomial, computed using a recurrence relation
 * algorithm. The returned polynomials include the same normalisation constant found in the
 * definition of the spherical harmonics, that is sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!).
 * The function accepts only positive values of m. For negative values, just use the normal
 * plegendre_lm and multiply it by (1-x^2)^(m/2).
 *
 * The credits for this function go to Press et al. 2002. ("Numerical Recipes, Third edition",
 * pag.294) for the computation of the Legendre polynomial, and to Wolfgang Ehrhardt
 * (https://groups.google.com/forum/#!topic/sci.math.num-analysis/13ZmcOtpWzg) for the
 * trick to rescale the function by (1-x^2)^(m/2) analitically.
 */
double plegendre_lm_rescaled_analytically (int l, int m, double x) {
  
  int i,ll;
  double fact,oldfact,pll,pmm,pmmp1,omx2;

  if (m < 0 || m > l || fabs(x) > 1.0) {
    printf("ERROR: Bad arguments in routine plegendre_lm\n");
 }

  pmm=1.0;
  /* Compute P_m^m/(1-x*x)^(m/2) */
  if (m > 0) {
    omx2=1;
    /* Uncomment the following line to compute the usual associate Legendre polynomials */
    // omx2=(1.0-x)*(1.0+x);
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
    /* Compute P_m+1^m */
    pmmp1=x*sqrt(2.0*m+3.0)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      /* Compute P_l^m, l>m+1 */
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





double plegendre_lm_rescaled (int l, int m, double x) {
  
  if (m == 0) {
    return plegendre_lm (l, m, x);
  }
  else if (m > 0) {
    return plegendre_lm_rescaled_analytically (l, m, x);
  }
  /* When m is negative there is a factor (-1)^m coming from the transformation rule of the
  associated Legendre polynomials. Note also that the sin(theta)^m factor flips from 
  denominator to the numerator.
  The negative-m rescaled polynomials are proportional to (1-x*x)^m, hence they go
  very fast to zero for x->1. We include that case by hand for optimisation purposes. */
  else if (m < 0) {
    if (fabs(x)!=1)
      return ALTERNATING_SIGN(m) * plegendre_lm (l, abs(m), x) * pow(1.0-x*x, 0.5*abs(m));
    else
      return 0;
  }

}



/**
 * Legendre Polynomials from alglib-3.6.0 (http://www.alglib.net/specialfunctions/polynomials/legendre.php)
 */
double plegendre (int n, double x)
{
    double a;
    double b;
    int i;
    double result;


    result = 1;
    a = 1;
    b = x;
    if( n==0 )
    {
        result = a;
        return result;
    }
    if( n==1 )
    {
        result = b;
        return result;
    }
    for(i=2; i<=n; i++)
    {
        result = ((2*i-1)*x*b-(i-1)*a)/i;
        a = b;
        b = result;
    }
    return result;
}








// =============================================================================
// =                        Multipole related functions                        =
// =============================================================================



/** 
 * Return the index corresponding to an l,m pair, considering m in [0,l].
 * If m_max is specified, then it is assumed that m is in [0,min(l,m_max)].
 * Example:
 *   l,m    multipole2offset_l_m(l, m, m_max=2)
 *   0,0 -> 0 
 *   1,0 -> 1 
 *   1,1 -> 2 
 *   2,0 -> 3 
 *   2,1 -> 4 
 *   2,2 -> 5 
 *   3,0 -> 6 
 *   3,1 -> 7 
 *   3,2 -> 8 
 *   4,0 -> 9 
 *   4,1 -> 10
 *   4,2 -> 11
 *   5,0 -> 12
 *   5,1 -> 13
 *   5,2 -> 14
 *   6,0 -> 15
 *   6,1 -> 16
 *   6,2 -> 17
 */
int multipole2offset_l_m(int l, int m, int m_max) {
  
  /* Check input (comment for more speed) */
  if (m > m_max) {
    printf("ERROR, %s: l = %3d, m = %3d, m_max = %3d  ---->  'm' out of range\n", __func__, l, m, m_max);
    return -1;
  }
  if (m > l) {
    printf("ERROR, %s: l = %3d, m = %3d, m_max = %3d  ---->  'm' cannot be larger than 'l'\n", __func__, l, m, m_max);
    return -1;
  }
  if (l < 0) {
    printf("ERROR, %s: l = %3d  ---->  'l' cannot be negative.\n", __func__, l);
    return -1;
  }
  
  int offset;
  
  /* If l<=m_max, then the index of the pair l,m is simply l + (l-1) + ... + 1 + m. */
  if (l <= m_max)
    offset = ((l+1)*l)/2 + m;
  /* If l>m_max, then any extra l>m_max contributes with a fixed offset of (l-m_max-1)*(m_max+1) */
  else
    offset = ((m_max+2)*(m_max+1))/2         /* Number of multipoles up to l=m_max+1 */
           + (l-m_max-1)*(m_max+1) + m;      /* Contribution coming from those l's that only have m_max+1 multipoles */

  /* Some debug */
  // printf("l = %3d, m = %3d, m_max = %3d  ---->  offset = %d\n", l, m, m_max, offset);

  return offset;

}



/**
 * Return the number of elements in the l,m hierarchy with l in [0,l_max] and m in [0,min(l,m_max)]
 * Example:
 *   l_max    m_max              lm_number_of_elements(l_max, m_max)
 *   0        0         ->       1 
 *   1        0         ->       2
 *   1        1         ->       3 
 *   1        2         ->       3 
 *   2        0         ->       3
 *   2        1         ->       5
 *   2        2         ->       6 
 *   3        0         ->       4
 *   3        1         ->       7
 *   3        2         ->       9
 *   3        3         ->       10
 */
int size_l_m(int l_max, int m_max) {
  
  return multipole2offset_l_m(l_max, MIN(l_max,m_max), m_max) + 1;

}





/**
 * Return the index corresponding to an L,M pair, considering that the M-index is constrained to the
 * values contained in m_vec.  Note that if m_vec = (0,1,2, ..., m_max), this function should give the
 * same result as multipole2offset_l_m(L, M, m_max).
 *
 * Examples:
 *	
 *    L,M    multipole2offset_l_indexm(L, M, [1,2])
 *    1,1 -> 0
 *    2,1 -> 1
 *    2,2 -> 2 
 *    3,1 -> 3
 *    3,2 -> 4
 *    4,1 -> 5
 *    4,2 -> 6 
 *		5,1 -> 7
 *		5,2 -> 8
 *
 *    L,M    multipole2offset_l_indexm(L, M, [0,2,4])
 *    0,0 -> 0
 *    1,0 -> 1
 *    2,0 -> 2 
 *    2,2 -> 3 
 *    3,0 -> 4
 *    3,2 -> 5
 *    4,0 -> 6
 *    4,2 -> 7 
 *		4,4 -> 8
 *		5,0 -> 9
 *		5,2 -> 10
 *		5,4 -> 11
 */
int multipole2offset_l_indexm (int L, int M, int * m_vec, int m_size) {
  

	/* Declaration of variables */
	int m=0, l=0, index_m=0;

	/* Build a logical array to check whether a given m is in m_vec */
	int m_max = m_vec[m_size-1];
	int m_in_vec[m_max+1];

	for(m=0; m<=m_max; ++m) {

		m_in_vec[m] = 0;

		for (index_m=0; index_m<m_size; ++index_m)
			if (m==m_vec[index_m]) ++m_in_vec[m];
			
		// *** Some debug
		// printf("m_in_vec[%d] = %d\n", m, m_in_vec[m]);
	}




	// *** Test input

	/* Check that M is contained in m_vec */
  if ((m_in_vec[M] == 0) || (M > m_max)) {
    printf("ERROR, %s: l = %3d, m = %3d, m_max = %3d  ---->  'm' not contained in m_vec.\n", __func__, L, M, m_max);
    return -1;
  }
	/* Check that there are no duplicate entries in m_vec */
	for(m=0; m<=m_max; ++m) {
	  if (m_in_vec[m] > 1) {
	    printf("ERROR, %s: l = %3d, m = %3d, m_max = %3d  ---->  m_vec has duplicate entries.\n", __func__, L, M, m_max);
	    return -1;
		}
  }
	if (M > L) {
    printf("ERROR, %s: l = %3d, m = %3d, m_max = %3d  ---->  'm' cannot be larger than 'l'\n", __func__, L, M, m_max);
    return -1;
  }
  if (L < 0) {
    printf("ERROR, %s: l = %3d  ---->  'l' cannot be negative.\n", __func__, L);
    return -1;
  }




  // *** Compute offset
  int offset = 0;

	for(l=0; l<=L; ++l) {
		for(m=0; m<=MIN(l,m_max); ++m) {

			/* During the last cycle on the l's, do not consider m's larger than M */
			if ( (l==L) && (m>M) ) continue;
			
			/* Count only those m's that are included in the list m_vec */
			if (m_in_vec[m] == 0) continue;

			/* What survives has to be counted */
			offset++;
			
			// *** Some debug
			// printf("l=%d, m=%d, offset=%d\n", l, m, offset-1);
		}
	}
	
  return offset - 1;

}







/**
 * Inverse function of multipole2offset_l_indexm. Given an offset, l_max and a list of m's,
 * determine the unique (l,m) pair associated with that offset. The output is written into
 * L, index_M.
 */
int offset2multipole_l_indexm (int offset, int l_max, int * m_vec, int m_size,
                               int * L, int * index_M)
{
    
  // *** Test input
	if (offset < 0) {
    printf("ERROR in %s: offset = %3d --> 'offset' should be a positive integer\n", __func__, offset);
    *L = *index_M = -1;
    return _FAILURE_;
  }

  int max_offset = size_l_indexm (l_max, m_vec, m_size) - 1;
	if (offset > max_offset) {
    printf("ERROR in %s: offset = %3d, max_offset = %3d --> 'offset' cannot be larger than 'max_offset'\n",
      __func__, offset, max_offset);
    *L = *index_M = -1;
    return _FAILURE_;
  }

  // *** Find (l,m) by calling multipole2offset_l_indexm
  int current_offset = -1;

  for (int l=0; l<=l_max; ++l) {
    
    for (int index_m=0; index_m<m_size; ++index_m) {
      
      int m = m_vec[index_m];

      /* Skip bad values of l and m */
      if (l < m) continue;
      current_offset = multipole2offset_l_indexm (l, m, m_vec, m_size);
      
      /* Debug multipole2offset */
      // printf("multipole2offset (%d,%d) = %d\n", l, m, current_offset);
      
      /* Check whether we found the correct (l,m) */
      if (current_offset == offset) {
        *index_M = index_m;
        break;
      }

    } // end of for(index_l)
    
    if (current_offset == offset) {
      *L = l;
      break;
    }
    
  } // end of for(index_m)

  /* Test that the output makes sense */
  if ((*L > l_max) || (*L<0) || (*index_M >= m_size) || (index_M<0)) {

    printf("ERROR in %s: result (L,index_M)=(%d,%d) is out of bounds L=[%d,%d], index_M=[%d,%d]\n",
      __func__, *L, *index_M, 0, l_max, 0, m_size-1);

    return _FAILURE_;
  }

  return _SUCCESS_;  
  
}





/**
 * Return the number of elements in the l,m hierarchy considering that the m-index is constrained to the
 * values contained in m_vec 
 */
int size_l_indexm(int l_max, int * m_vec, int m_size) {
    
  /* Return 0 and a warning if the smallest m is larger than l_max */
  if (m_vec[0] > l_max) {
    // printf ("WARNING in lm_number_of_elements_mlist: the smallest m=%d is larger than l_max=%d. Returning 0\n",
    //   m_vec[0], l_max);
    return 0;
  }
    
  /* Determine the maximum allowed value of m as the largest element in m that is also smaller than l_max */
  int index_m_max = m_size-1;
  while (m_vec[index_m_max] > l_max) --index_m_max;
  int m_max = m_vec[index_m_max];
	return multipole2offset_l_indexm (l_max, m_max, m_vec, m_size) + 1;

}






/**
 * Return the index corresponding to an L,M pair, considering that the L-index and M-index are
 * constrained to the values contained in l_vec and m_vec, respectively. Note that if
 * l_vec = (0,1,2,...,l_max) and m_vec = (0,1,2, ..., m_max), this function should give the
 * same result as multipole2offset_l_m (L, M, m_max).
 *
 * Examples:
 *	
 *    L,M    multipole2offset_indexl_indexm(L, M, [0,2,4], [1,2])
 *    2,1 -> 0
 *    2,2 -> 1 
 *    4,1 -> 2
 *    4,2 -> 3
 *
 *    L,M    multipole2offset_indexl_indexm(L, M, [1,2], [0,2,4])
 *    1,0 -> 0
 *    2,0 -> 1 
 *    2,2 -> 2
 *
 *    L,M    multipole2offset_indexl_indexm(L, M, [1,2,3,5,7,8], [0,2,4])
 *    1,0 -> 0
 *    2,0 -> 1
 *    2,2 -> 2
 *    3,0 -> 3
 *    3,2 -> 4 
 *    5,0 -> 5
 *    5,2 -> 6
 *    5,4 -> 7 
 *    7,0 -> 8
 *    7,2 -> 9
 *    7,4 -> 10 
 *    8,0 -> 11
 *    8,2 -> 12
 *    8,4 -> 13
 *
 */
int multipole2offset_indexl_indexm(int L, int M, int * l_vec, int l_size, int * m_vec, int m_size) {
  
	int index_l=0, index_m=0;


  /* Find the position of L inside l_vec */
	int index_L=-1;
	for(index_l=0; index_l<l_size; ++index_l)
		if (l_vec[index_l] == L) index_L = index_l;

  /* Find the position of M inside m_vec */
	int index_M=-1;
	for(index_m=0; index_m<m_size; ++index_m)
		if (m_vec[index_m] == M) index_M = index_m;
  
	// *** Some debug
  // printf("L=%d,index_L=%d,M=%d,index_M=%d\n",L,index_L,M,index_M);
  
  // ************* Check the input ***************
  
  /* Check that L is contained in l_vec */
  if (index_L==-1) {
    printf("ERROR in %s: l=%d not contained in l_vec.\n", __func__, L);
    return -1;
  }
  
  /* Check that M is contained in m_vec */
  if (index_M==-1) {
    printf("ERROR in %s: m=%d not contained in m_vec.\n", __func__, M);
    return -1;
  }
  
  /* Check that l_vec is sorted in strictly ascending order */    
  for(index_l=0; index_l<(l_size-1); ++index_l) {
    if (l_vec[index_l] >= l_vec[index_l+1]) {
      printf("ERROR in %s: l_vec is not sorted in strictly ascending order.\n", __func__);
      return -1;
    }
  }
  
  /* Check that m_vec is sorted in strictly ascending order */    
  for(index_m=0; index_m<(m_size-1); ++index_m) {
    if (m_vec[index_m] >= m_vec[index_m+1]) {
      printf("ERROR in %s: m_vec is not sorted in strictly ascending order.\n", __func__);
      return -1;
    }
  }
  
  /* Check that M is smaller than or equal to L */
  if (M > L) {
    printf("ERROR in %s: l = %3d, m = %3d  ---->  'm' cannot be larger than 'l'\n", __func__, L, M);
    return -1;
  }
  
  /* Check that L is positive */
  if (L < 0) {
    printf("ERROR in %s: l = %3d  ---->  'l' cannot be negative.\n", __func__, L);
    return -1;
  }


  // ***********  Compute offset *************
  int offset = 0;

	for(index_l=0; index_l<=index_L; ++index_l) {

		for(index_m=0; index_m<m_size; ++index_m) {

			/* Always ensure that m <= l */
			if (l_vec[index_l] < m_vec[index_m]) continue;

			/* During the last cycle on the l's, do not consider m's larger than M */
			if ( (index_l==index_L) && (index_m>index_M) ) continue;

			/* What survives has to be counted */
			offset++;
			
			// *** Some debug
      // printf("l=%d, m=%d, offset=%d\n", l_vec[index_l], m_vec[index_m], offset-1);
		}
	}

	// *** Some debug
  // printf("L=%d,index_L=%d,M=%d,index_M=%d,offset=%d\n",L,index_L,M,index_M,offset-1);

  return offset - 1;

}




/**
 * Return the number of elements in the l,m hierarchy considering that the l-index
 * and m-index are constrained to the values contained in l_vec and m_vec, respectively.
 */
int size_indexl_indexm(int * l_vec, int l_size, int * m_vec, int m_size) {

  int l_max = l_vec[l_size-1];
  int m_max = m_vec[m_size-1];
	return multipole2offset_indexl_indexm(l_max, m_max, l_vec, l_size, m_vec, m_size) + 1;
	
}



/**
 * Inverse function of multipole2offset_indexl_indexm. Given an offset and two lists of possible l's and m's,
 * determine the unique (l,m) pair associated with that offset. The output is written into
 * index_L, index_M.
 */
int offset2multipole_indexl_indexm(int offset, int * l_vec, int l_size, int * m_vec, int m_size,
                                   int * index_L, int * index_M)
{
    
  // *** Test input
	if (offset < 0) {
    printf("ERROR in %s: offset = %3d --> 'offset' should be a positive integer\n", __func__, offset);
    *index_L = *index_M = -1;
    return _FAILURE_;
  }

  int max_offset = size_indexl_indexm (l_vec, l_size, m_vec, m_size) - 1;
	if (offset > max_offset) {
    printf("ERROR in %s: offset = %3d, max_offset = %3d --> 'offset' cannot be larger than 'max_offset'\n",
      __func__, offset, max_offset);
    *index_L = *index_M = -1;
    return _FAILURE_;
  }
  

  // *** Find (l,m) by calling multipole2offset
  int current_offset = -1;

  for (int index_l=0; index_l<l_size; ++index_l) {

    int l = l_vec[index_l];
    
    for (int index_m=0; index_m<m_size; ++index_m) {
      
      int m = m_vec[index_m];

      /* Skip bad values of l and m */
      if (l < m) continue;
      current_offset = multipole2offset_indexl_indexm (l, m, l_vec, l_size, m_vec, m_size);
      
      /* Debug multipole2offset */
      // printf("multipole2offset (%d,%d) = %d\n", l, m, current_offset);
      
      /* Check whether we found the correct (l,m) */
      if (current_offset == offset) {
        *index_M = index_m;
        break;
      }

    } // end of for(index_l)
    
    if (current_offset == offset) {
      *index_L = index_l;
      break;
    }
    
  } // end of for(index_m)

  /* Test that the output makes sense */
  if ((*index_L >= l_size) || (*index_L<0) || (*index_M >= m_size) || (index_M<0)) {

    printf("ERROR in %s: result (index_L,index_M)=(%d,%d) is out of bounds index_L=[%d,%d], index_M=[%d,%d]\n",
      __func__, *index_L, *index_M, 0, l_size-1, 0, m_size-1);

    return _FAILURE_;
  }

  return _SUCCESS_;  
  
}



/* Index the massive hierarchy with the tree indices n,l,m.  'n' can be any positive integer,
 * l should be a positive integer smaller than 'n', and 'm' should be a positive integer smaller
 * than 'l'.
 * Example with l_max=2, m_max=1:
 *   n,l,m               multipole2offset_unconstrained_n_l_m(n,l,m,2,1)
 *   0,0,0      ->       0 
 *   1,0,0      ->       1
 *   1,1,0      ->       2
 *   1,1,1      ->       3
 *   2,0,0      ->       4 
 *   2,1,0      ->       5    
 *   2,1,1      ->       6 
 *   2,2,0      ->       7
 *   2,2,1      ->       8
 *   3,0,0      ->       9
 *   3,1,0      ->       10
 *   3,1,1      ->       11
 *   3,2,0      ->       12
 *   3,2,1      ->       13
 *   4,0,0      ->       14
 *   4,1,0      ->       15 
 *   4,1,1      ->       16
 *   4,2,0      ->       17
 *   4,2,1      ->       18
 */
int multipole2offset_unconstrained_n_l_m(int n, int l, int m, int l_max, int m_max) {
  
  // Check input
  if (l > n) {
    printf("ERROR, %s: n = %3d, l = %3d  ---->  'l' cannot be larger than 'n'\n", __func__, n, l);
    return _FAILURE_;
  }

  int n_i=0;
  int offset=0;

  for(n_i=0; n_i<n; ++n_i) {
    /* Each value 'n_i' of 'n' contributes with an lm hierarchy where l_max=n_i.
      However, we impose for l to be upper-limited by l_max, hence the actual
      hierarchy stops at MIN(l_max,n_i) */
    offset += size_l_m (MIN(n_i,l_max), m_max);
  }
  
  /* Note that if n=0, then we do not enter in the above loop and we recover the lm_offset */
  return offset + multipole2offset_l_m (l,m,m_max);
  
}




/* This function is the same as 'multipole2offset_unconstrained_n_l_m', but it only allows
 * l to be equal to 0 or n. This is a dirty trick to be able to address the following list
 * of (n,l,m) combinations:
 *  (0,0,0)
 *  (1,1,0)
 *  (2,0,0)
 *  (2,2,0)
 *  (1,1,1)
 *  (2,2,1)
 *  (2,2,2)
 * which are the closed system of beta-moments needed to evolve the baryon and CDM hierarchies.
 * We truncate the hierarchy at n=2 because the quadratic sources for n>2 involve the first-order
 * higher moments (like the shear) which vanish for a perfect fluid.
 *
 * Example with l_max=3, m_max=2:
 *   n,l,m               multipole2offset_n_l_m(n,l,m,2,1)
 *   0,0,0      ->       0 
 *   1,1,0      ->       1
 *   1,1,1      ->       2
 *   2,0,0      ->       3 
 *   2,2,0      ->       4
 *   2,2,1      ->       5
 *   2,2,2      ->       6
 *   3,3,0      ->       7
 *   3,3,1      ->       8
 *   3,3,2      ->       9
 *   4,0,0      ->       10
 *   5,0,0      ->       11
 *   6,0,0      ->       12
 */
int multipole2offset_n_l_m(int n, int l, int m, int l_max, int m_max) {
  
  /* Check input */
  if ( (l!=0) && (l!=n) ) {
    printf("ERROR, %s: l = %3d, n = %3d  ---->  'l' can be only equal to 'n' or '0'.\n", __func__, l, n);
    return _FAILURE_;
  }

  if ( (l==0) && (n%2!=0) ) {
    printf("ERROR, %s: l = 0, n = %3d  ---->  the input l=0 is allowed only when 'n' is even.\n", __func__, n);
    return _FAILURE_;    
  }

  if ( m>MIN(l,m_max) ) {
    printf("ERROR, %s: m = %3d, l = %3d, m_max = %3d ---->  'm' must be smaller than MIN(l,m_max).\n",
      __func__, m, l, m_max);
    return _FAILURE_;
  }

  int n_i=0;
  int offset=0;

  /* If n=0 (and hence l=m=0), we do not enter the below loop and return offset=0 */
  if (n==0)
    return offset;

  /* Each value 'n_i' of 'n' contributes with its l=0 and l=n_i multipoles.
  We skip the n_i=0 iteration because we want nlm(0,0,0) = 0. */
  for(n_i=1; n_i<n; ++n_i) {

    /* l=0 contribution, but only if 'n_i' is even.  We do not want to consider
    the 1,0,0 multipole because it does not matter for the cold matter hierarchy. */
    if (n_i%2==0)
      offset += 1;
    
    /* l=n_i contribution */
    if(l_max>=n_i)
      offset += MIN(m_max,n_i) + 1;
    
  }
  
  /* n_i=n contribution */
  if (n%2==0)  
    offset += 1;
  if (l==n)
    offset += MIN(m,m_max) + 1;
  
  return offset;
  
}


/* Return the number of elements in a massive hierarchy with n_max, l_max, m_max,
 * where the constraints are: l=0 or l=n, if 'n' even there is no l=n moment.
 */   
int size_n_l_m (int n_max, int l_max, int m_max) {
 
  return multipole2offset_n_l_m (n_max, MIN(l_max,n_max), MIN(l_max,m_max), l_max, m_max) + 1;
  
}
  
  
  
/* This function is the same as 'multipole2offset_n_l_m', but it only allows m to be in m_vec.
 * This is a dirty trick to be able to address the following list of (n,l,m) combinations:
 *  (0,0,0)
 *  (1,1,0)
 *  (2,0,0)
 *  (2,2,0)
 *  (1,1,1)
 *  (2,2,1)
 *  (2,2,2)
 * which are the closed system of beta-moments needed to evolve the baryon and CDM hierarchies.
 * We truncate the hierarchy at n=2 because the quadratic sources for n>2 involve the first-order
 * higher moments (like the shear) which vanish for a perfect fluid.
 *
 */
int multipole2offset_n_l_indexm(int N, int L, int M, int l_max, int * m_vec, int m_size) {


	/* Declaration of variables */
	int n=0, l=0, m=0, index_m=0;

  /* Find the position of M inside m_vec */
	int index_M=-1;
	for(index_m=0; index_m<m_size; ++index_m)
		if (m_vec[index_m] == M) index_M = index_m;


  /* Check that M is contained in m_vec */
  if (index_M==-1) {
    printf("ERROR in %s: m=%d not contained in m_vec.\n", __func__, M);
    return -1;
  }
  
  /* Check that m_vec is sorted in strictly ascending order */    
  for(index_m=0; index_m<(m_size-1); ++index_m) {
    if (m_vec[index_m] >= m_vec[index_m+1]) {
      printf("ERROR in %s: m_vec is not sorted in strictly ascending order.\n", __func__);
      return -1;
    }
  }
  
  /* Check that M is smaller than or equal to L */
  if (M > L) {
    printf("ERROR in %s: l = %3d, m = %3d  ---->  'm' cannot be larger than 'l'\n", __func__, L, M);
    return -1;
  }
  
  /* Check that L is positive */
  if (L < 0) {
    printf("ERROR in %s: l = %3d  ---->  'l' cannot be negative.\n", __func__, L);
    return -1;
  }

  /* L can only be equal to 0 or N */
  if ( (L!=0) && (L!=N) ) {
    printf("ERROR in %s: l = %3d, n = %3d  ---->  'l' can be only equal to 'n' or '0'.\n", __func__, L, N);
    return _FAILURE_;
  }

  /* L can be zero only when N is even */
  if ( (L==0) && (N%2!=0) ) {
    printf("ERROR in %s: l = 0, n = %3d  ---->  the input l=0 is allowed only when 'n' is even.\n", __func__, N);
    return _FAILURE_;    
  }

  // *** Compute offset
  int offset = 0;

  for (n=0; n<=N; ++n) {

  	for (l=0; l<=MIN(n,l_max); ++l) {
      
			/* During the last cycle on the n's, do not consider l's larger than L */
			if ((n==N) && (l>L)) continue;
      
      /* l can only be equal to n or to zero */
      if ((l!=0) && (l!=n)) continue;
      
      /* L can be zero only when N is even */
      if ((l==0) && (n%2!=0) ) continue;
      
  		for (index_m=0; index_m<m_size; ++index_m) {

  			/* Always ensure that m <= l */
  			if (l < m_vec[index_m]) continue;

        /* During the last cycle on the l's, do not consider m's larger than M */
  			if ((l==L) && (index_m>index_M)) continue;

  			/* What survives has to be counted */
  			offset++;
			
  			// *** Some debug
        // printf("n=%d, l=%d, m=%d, index_m=%d, offset=%d\n", n, l, m_vec[index_m], index_m, offset-1);

  		} // end of for(m)
  	} // end of for(l)
  } // end of for(n)
	
  return offset - 1;

}

  
/* Return the number of elements in a massive hierarchy with n_max, l_max, m_max,
 * where the constraints are: l=0 or l=n, if 'n' even there is no l=n moment.
 */   
int size_n_l_indexm (int n_max, int l_max, int * m_vec, int m_size) {
 
  /* Return zero if the smallest m is larger than l_max */
  if (m_vec[0] > l_max)
    return 0;
    
  /* Determine the maximum allowed value of m as the largest element in m that is also smaller than l_max */
  int index_m_max = m_size-1;
  while (m_vec[index_m_max] > l_max) --index_m_max;
  int m_max = m_vec[index_m_max];

  return multipole2offset_n_l_indexm(n_max, MIN(n_max,l_max), MIN(l_max,m_max), l_max, m_vec, m_size) + 1;
  
}
  
  
  
  
  
  
  
  
  
  
// =====================================================================================
// =                                     Interpolation                                 =
// =====================================================================================
  
  
  
/**
 * Function to calculate second derivatives of ppt->sources at the nodes.
 * This function reflects the very specific indexing pattern of ppt->sources, and therefore
 * is not easily recycleable outside CLASS.
 */
int spline_sources_derivs_two_levels(
			     double * x, /* vector of size tau_size */
			     int tau_size,
			     double ** y_array,
			     int tp_size,   
			     double ** ddy_array,
			     short spline_mode,
			     ErrorMsg errmsg
			     )
{

  double * p;
  double * qn;
  double * un; 
  double * u;
  double sig;
  int index_tau;
  int index_tp;
  double dy_first;
  double dy_last;

  u = malloc((tau_size-1) * tp_size * sizeof(double));
  p = malloc(tp_size * sizeof(double));
  qn = malloc(tp_size * sizeof(double));
  un = malloc(tp_size * sizeof(double));

  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }
  if (p == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate p",__func__,__LINE__);
    return _FAILURE_;
  }
  if (qn == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate qn",__func__,__LINE__);
    return _FAILURE_;
  }
  if (un == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate un",__func__,__LINE__);
    return _FAILURE_;
  }

  // *** Print some debug info
  // printf("tp_size = %d\n", tp_size);
  // printf("Allocated vectors u, p, qn, un\n");
  
  // Equivalent to index_tau, while index_tp is equivalent to index_tp
  index_tau=0;

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {
      // ddy_array[index_tau*tp_size+index_tp] = u[index_tau*tp_size+index_tp] = 0.0;
      ddy_array[index_tp][index_tau] = u[index_tau*tp_size+index_tp] = 0.0;
    }
      printf("DONE ***");    
  }
  else if (spline_mode == _SPLINE_EST_DERIV_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {

    	dy_first = 
    	  ((x[2]-x[0])*(x[2]-x[0])*
         // (y_array[1*tp_size+index_tp]-y_array[0*tp_size+index_tp])-
    	   (y_array[index_tp][1]
    	    - y_array[index_tp][0]) -
    	   (x[1]-x[0])*(x[1]-x[0])*
         // (y_array[2*tp_size+index_tp]-y_array[0*tp_size+index_tp]))/
         (y_array[index_tp][2] - 
           y_array[index_tp][0]))/           
    	  ((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));

      // ddy_array[index_tau*tp_size+index_tp] = -0.5;
    	ddy_array[index_tp][index_tau] = -0.5;        

    	u[index_tau*tp_size+index_tp] = 
    	  (3./(x[1] -  x[0]))*
        // ((y_array[1*tp_size+index_tp]-y_array[0*tp_size+index_tp])/
    	  ((y_array[index_tp][1] - 
    	    y_array[index_tp][0])/          
    	   (x[1] - x[0])-dy_first);

    }
  }
  else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
  }


  for (index_tau=1; index_tau < tau_size-1; index_tau++) {

    sig = (x[index_tau] - x[index_tau-1])/(x[index_tau+1] - x[index_tau-1]);

    for (index_tp=0; index_tp < tp_size; index_tp++) {

      // p[index_tp] = sig * ddy_array[(index_tau-1)*tp_size+index_tp] + 2.0;
      p[index_tp] = sig * ddy_array[index_tp][(index_tau-1)] + 2.0;        

      // ddy_array[index_tau*tp_size+index_tp] = (sig-1.0)/p[index_tp];
      ddy_array[index_tp][index_tau] = (sig-1.0)/p[index_tp];        

      u[index_tau*tp_size+index_tp] = 
        // (y_array[(index_tau+1)*tp_size+index_tp] - y_array[index_tau*tp_size+index_tp])
        (y_array[index_tp][(index_tau+1)] -
         y_array[index_tp][index_tau])
        / (x[index_tau+1] - x[index_tau])
        // - (y_array[index_tau*tp_size+index_tp] - y_array[(index_tau-1)*tp_size+index_tp])
        - (y_array[index_tp][index_tau] -
           y_array[index_tp][(index_tau-1)])
        / (x[index_tau] - x[index_tau-1]);

      u[index_tau*tp_size+index_tp] = (6.0 * u[index_tau*tp_size+index_tp] /
        (x[index_tau+1] - x[index_tau-1]) 
        - sig * u[(index_tau-1)*tp_size+index_tp]) / p[index_tp];
        
      /* Some debug */
      // printf("y_array[index_tp][index_tau] = %g\n", y_array[index_tp][index_tau]);
        
    }

  }

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {
      qn[index_tp]=un[index_tp]=0.0;
    }

  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_tp=0; index_tp < tp_size; index_tp++) {

        dy_last = 
          ((x[tau_size-3]-x[tau_size-1])*(x[tau_size-3]-x[tau_size-1])*
           // (y_array[(tau_size-2)*tp_size+index_tp]-y_array[(tau_size-1)*tp_size+index_tp])-
           (y_array[index_tp][(tau_size-2)]
           - y_array[index_tp][(tau_size-1)])-     
           (x[tau_size-2]-x[tau_size-1])*(x[tau_size-2]-x[tau_size-1])*
           // (y_array[(tau_size-3)*tp_size+index_tp]-y_array[(tau_size-1)*tp_size+index_tp]))/
           (y_array[index_tp][(tau_size-3)]
           - y_array[index_tp][(tau_size-1)]))/
          ((x[tau_size-3]-x[tau_size-1])*(x[tau_size-2]-x[tau_size-1])*(x[tau_size-3]-x[tau_size-2]));

        qn[index_tp]=0.5;

        un[index_tp]=
          (3./(x[tau_size-1] - x[tau_size-2]))*
          // (dy_last-(y_array[(tau_size-1)*tp_size+index_tp] - y_array[(tau_size-2)*tp_size+index_tp])/
          (dy_last-(y_array[index_tp][(tau_size-1)]
                    - y_array[index_tp][(tau_size-2)])/
           (x[tau_size-1] - x[tau_size-2])); 

      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
  
  index_tau=tau_size-1;

  for (index_tp=0; index_tp < tp_size; index_tp++) {
    // ddy_array[index_tau*tp_size+index_tp] = 
    ddy_array[index_tp][index_tau] =
      (un[index_tp] - qn[index_tp] * u[(index_tau-1)*tp_size+index_tp]) /
      // (qn[index_tp] * ddy_array[(index_tau-1)*tp_size+index_tp] + 1.0);
      (qn[index_tp] * ddy_array[index_tp][(index_tau-1)] + 1.0);        
  }

  for (index_tau=tau_size-2; index_tau >= 0; index_tau--) {
    for (index_tp=0; index_tp < tp_size; index_tp++) {

      // ddy_array[index_tau*tp_size+index_tp] = ddy_array[index_tau*tp_size+index_tp]
      // * ddy_array[(index_tau+1)*tp_size+index_tp] + u[index_tau*tp_size+index_tp];
      ddy_array[index_tp][index_tau] =
        ddy_array[index_tp][index_tau]
        * ddy_array[index_tp][(index_tau+1)]
        + u[index_tau*tp_size+index_tp];

      /* Some debug */
      // printf("ddy_array[index_tp][index_tau] = %g\n", ddy_array[index_tp][index_tau]);

    }
  }

  free(qn);
  free(un);
  free(p);
  free(u);

  return _SUCCESS_;
}
  
  
  
  
int spline_sources_interpolate_two_levels(
			     double * x_array,
			     int tau_size,
			     double ** y_array,
			     double ** ddy_array,
			     int tp_size,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to tp_size */
			     ErrorMsg errmsg
			     )
{
 
  int inf,sup,mid;
  double h,a,b;
  
  inf=0;
  sup=tau_size-1;
  

  // =================
  // = Table look-up =
  // =================
  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    /* Table lookup with bisection */
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}
    }
  }
  else {

    if (x < x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }
    if (x > x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    /* Table lookup with bisection */
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }
  }


  // =================
  // = Interpolation =
  // =================

  /* Store the index for use in the next call of the closeby function */
  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  int index_tp;
  for (index_tp=0; index_tp<result_size; index_tp++) {

    /* Cubic interpolation, see eq. 3.3.3 of Numerical Recipes in C.
      To have plain linear interpolation, comment the last two lines. */
    result[index_tp] = 
      a * y_array[index_tp][inf] +
      b * y_array[index_tp][sup] +
      ((a*a*a-a)* ddy_array[index_tp][inf] + 
       (b*b*b-b)* ddy_array[index_tp][sup])*h*h/6.;


  }
  
  
  
  return _SUCCESS_;
}

  
  
  
  
int spline_sources_interpolate_two_levels_growing_closeby(
			     double * x_array,
			     int tau_size,
			     double ** y_array,
			     double ** ddy_array,
			     int tp_size,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to tp_size */
			     ErrorMsg errmsg
			     )
{

  int inf,sup;
  double h,a,b;

  inf = *last_index;
  class_test(inf<0 || inf>(tau_size-1),
      errmsg,
      "*lastindex=%d out of range [0:%d]\n",inf,tau_size-1);
      
  // *** Look at the left of last_index
  while (x < x_array[inf]) {
    inf--;
    if (inf < 0) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
       x,x_array[0]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  // *** Look at the right of last_index
  while (x > x_array[sup]) {
    sup++;
    if (sup > (tau_size-1)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
       x,x_array[tau_size-1]);
      return _FAILURE_;
    }
  }
  inf = sup-1;
  
  // *** Store the position for later use
  *last_index = inf;

  // *** Actual interpolation
  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  int index_tp;
  for (index_tp=0; index_tp<result_size; index_tp++) {
    
    /* Cubic interpolation, see eq. 3.3.3 of Numerical Recipes in C.
      To have plain linear interpolation, comment the last two lines. */
    result[index_tp] = 
      a * y_array[index_tp][inf] +
      b * y_array[index_tp][sup] +
      ((a*a*a-a)* ddy_array[index_tp][inf] + 
       (b*b*b-b)* ddy_array[index_tp][sup])*h*h/6.;

  }

  return _SUCCESS_;
}
  
  
  
  


/**
 * Compute and store the second-derivatives at the node points, needed for later cubic spline interpolation.
 * This function is a readaptation of CLASS1 "array_spline_table_lines" that takes into account the 
 * unusual indexing used for the ppt->sources array, that is:
 *  sources[index_mode]
 *    [index_ic * ppt->tp_size[index_mode] + index_type]
 *    [index_tau * ppt->k_size[index_mode] + index_k]
 * 
 * Note that ppt->dd_sources has the same kind of indexing.
 */

int spline_sources_derivs(
			     double * x, /* vector of size tau_size */
			     int tau_size,
			     double *** y_array, /* array of size tau_size*tp_size with elements 
						  y_array[index_tau*tp_size+index_tp] */
			     int tp_size,   
			     double *** ddy_array, /* array of size tau_size*tp_size */
			     short spline_mode,
           int index_mode,
           int index_ic,           
           int index_k,
           int k_size,
			     ErrorMsg errmsg
			     )
{

  double * p;
  double * qn;
  double * un; 
  double * u;
  double sig;
  int index_tau;
  int index_tp;
  double dy_first;
  double dy_last;

  u = malloc((tau_size-1) * tp_size * sizeof(double));
  p = malloc(tp_size * sizeof(double));
  qn = malloc(tp_size * sizeof(double));
  un = malloc(tp_size * sizeof(double));

  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }
  if (p == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate p",__func__,__LINE__);
    return _FAILURE_;
  }
  if (qn == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate qn",__func__,__LINE__);
    return _FAILURE_;
  }
  if (un == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate un",__func__,__LINE__);
    return _FAILURE_;
  }

  // *** Print some debug info
  // printf("tp_size = %d, index_mode = %d, index_ic = %d, index_k = %d\n", tp_size, index_mode, index_ic, index_k);
  // printf("Allocated vectors u, p, qn, un\n");
  
  // Equivalent to index_tau, while index_tp is equivalent to index_tp
  index_tau=0;

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {
      // ddy_array[index_tau*tp_size+index_tp] = u[index_tau*tp_size+index_tp] = 0.0;
      ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] = u[index_tau*tp_size+index_tp] = 0.0;
    }
      printf("DONE ***");    
  }
  else if (spline_mode == _SPLINE_EST_DERIV_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {

    	dy_first = 
    	  ((x[2]-x[0])*(x[2]-x[0])*
         // (y_array[1*tp_size+index_tp]-y_array[0*tp_size+index_tp])-
    	   (y_array[index_mode][index_ic*tp_size + index_tp][1*k_size+index_k]
    	    - y_array[index_mode][index_ic*tp_size + index_tp][0*k_size+index_k]) -
    	   (x[1]-x[0])*(x[1]-x[0])*
         // (y_array[2*tp_size+index_tp]-y_array[0*tp_size+index_tp]))/
         (y_array[index_mode][index_ic*tp_size + index_tp][2*k_size+index_k] - 
           y_array[index_mode][index_ic*tp_size + index_tp][0*k_size+index_k]))/           
    	  ((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));

      // ddy_array[index_tau*tp_size+index_tp] = -0.5;
    	ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] = -0.5;        

    	u[index_tau*tp_size+index_tp] = 
    	  (3./(x[1] -  x[0]))*
        // ((y_array[1*tp_size+index_tp]-y_array[0*tp_size+index_tp])/
    	  ((y_array[index_mode][index_ic*tp_size + index_tp][1*k_size+index_k] - 
    	    y_array[index_mode][index_ic*tp_size + index_tp][0*k_size+index_k])/          
    	   (x[1] - x[0])-dy_first);

    }
  }
  else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
  }


  for (index_tau=1; index_tau < tau_size-1; index_tau++) {

    sig = (x[index_tau] - x[index_tau-1])/(x[index_tau+1] - x[index_tau-1]);

    for (index_tp=0; index_tp < tp_size; index_tp++) {

      // p[index_tp] = sig * ddy_array[(index_tau-1)*tp_size+index_tp] + 2.0;
      p[index_tp] = sig * ddy_array[index_mode][index_ic*tp_size + index_tp][(index_tau-1)*k_size+index_k] + 2.0;        

      // ddy_array[index_tau*tp_size+index_tp] = (sig-1.0)/p[index_tp];
      ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] = (sig-1.0)/p[index_tp];        

      u[index_tau*tp_size+index_tp] = 
        // (y_array[(index_tau+1)*tp_size+index_tp] - y_array[index_tau*tp_size+index_tp])
        (y_array[index_mode][index_ic*tp_size + index_tp][(index_tau+1)*k_size+index_k] -
         y_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k])
        / (x[index_tau+1] - x[index_tau])
        // - (y_array[index_tau*tp_size+index_tp] - y_array[(index_tau-1)*tp_size+index_tp])
        - (y_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] -
           y_array[index_mode][index_ic*tp_size + index_tp][(index_tau-1)*k_size+index_k])
        / (x[index_tau] - x[index_tau-1]);

      u[index_tau*tp_size+index_tp] = (6.0 * u[index_tau*tp_size+index_tp] /
        (x[index_tau+1] - x[index_tau-1]) 
        - sig * u[(index_tau-1)*tp_size+index_tp]) / p[index_tp];
    }

  }

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_tp=0; index_tp < tp_size; index_tp++) {
      qn[index_tp]=un[index_tp]=0.0;
    }

  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_tp=0; index_tp < tp_size; index_tp++) {

        dy_last = 
          ((x[tau_size-3]-x[tau_size-1])*(x[tau_size-3]-x[tau_size-1])*
           // (y_array[(tau_size-2)*tp_size+index_tp]-y_array[(tau_size-1)*tp_size+index_tp])-
           (y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-2)*k_size+index_k]
           - y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-1)*k_size+index_k])-     
           (x[tau_size-2]-x[tau_size-1])*(x[tau_size-2]-x[tau_size-1])*
           // (y_array[(tau_size-3)*tp_size+index_tp]-y_array[(tau_size-1)*tp_size+index_tp]))/
           (y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-3)*k_size+index_k]
           - y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-1)*k_size+index_k]))/
          ((x[tau_size-3]-x[tau_size-1])*(x[tau_size-2]-x[tau_size-1])*(x[tau_size-3]-x[tau_size-2]));

        qn[index_tp]=0.5;

        un[index_tp]=
          (3./(x[tau_size-1] - x[tau_size-2]))*
          // (dy_last-(y_array[(tau_size-1)*tp_size+index_tp] - y_array[(tau_size-2)*tp_size+index_tp])/
          (dy_last-(y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-1)*k_size+index_k]
                    - y_array[index_mode][index_ic*tp_size + index_tp][(tau_size-2)*k_size+index_k])/
           (x[tau_size-1] - x[tau_size-2])); 

      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
  
  index_tau=tau_size-1;

  for (index_tp=0; index_tp < tp_size; index_tp++) {
    // ddy_array[index_tau*tp_size+index_tp] = 
    ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] =
      (un[index_tp] - qn[index_tp] * u[(index_tau-1)*tp_size+index_tp]) /
      // (qn[index_tp] * ddy_array[(index_tau-1)*tp_size+index_tp] + 1.0);
      (qn[index_tp] * ddy_array[index_mode][index_ic*tp_size + index_tp][(index_tau-1)*k_size+index_k] + 1.0);        
  }

  for (index_tau=tau_size-2; index_tau >= 0; index_tau--) {
    for (index_tp=0; index_tp < tp_size; index_tp++) {

      // ddy_array[index_tau*tp_size+index_tp] = ddy_array[index_tau*tp_size+index_tp]
      // * ddy_array[(index_tau+1)*tp_size+index_tp] + u[index_tau*tp_size+index_tp];
      ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k] =
        ddy_array[index_mode][index_ic*tp_size + index_tp][index_tau*k_size+index_k]
        * ddy_array[index_mode][index_ic*tp_size + index_tp][(index_tau+1)*k_size+index_k]
        + u[index_tau*tp_size+index_tp];


    }
  }

  free(qn);
  free(un);
  free(p);
  free(u);

  return _SUCCESS_;
}
 
 
/**
 * Find the interpolated value of all the sources at a certain time.  The array with the second derivatives
 * at the nodes (double *** dd_array) must have been already computed by using the function
 * "spline_sources_derivs". This function is a readaptation of CLASS1 "array_interpolate_spline" that takes
 * into account the unusual indexing used for the ppt->sources array, that is:
 *   sources[index_mode]
 *     [index_ic * ppt->tp_size[index_mode] + index_type]
 *     [index_tau * ppt->k_size[index_mode] + index_k]
 */
int spline_sources_interpolate(
			     double * x_array,
			     int tau_size,
			     double *** y_array,
			     double *** ddy_array,
			     int tp_size,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to tp_size */
           int index_mode,
           int index_ic,           
           int index_k,
           int k_size,
			     ErrorMsg errmsg
			     )
{
 
  int inf,sup,mid;
  double h,a,b;
  
  inf=0;
  sup=tau_size-1;
  

  // =================
  // = Table look-up =
  // =================
  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    // Table lookup with bisection
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}
    }
  }
  else {

    if (x < x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }
    if (x > x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    // Table lookup with bisection
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }
  }


  // =================
  // = Interpolation =
  // =================

  /* Store the index for use in the next call of the closeby function */
  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  int index_tp;
  for (index_tp=0; index_tp<result_size; index_tp++) {

    /* Cubic interpolation, see eq. 3.3.3 of Numerical Recipes in C.
      To have plain linear interpolation, comment the last two lines. */
    result[index_tp] = 
      a * y_array[index_mode][index_ic*tp_size + index_tp][inf*k_size+index_k] +
      b * y_array[index_mode][index_ic*tp_size + index_tp][sup*k_size+index_k] +
      ((a*a*a-a)* ddy_array[index_mode][index_ic*tp_size + index_tp][inf*k_size+index_k] + 
       (b*b*b-b)* ddy_array[index_mode][index_ic*tp_size + index_tp][sup*k_size+index_k])*h*h/6.;


  }
  
  return _SUCCESS_;
}
 
 
 
 
/* Same as spline_sources_interpolate, but faster if x is arranged in growing order, and the point x is
presumably very close to the previous point x from the last call of this function. */

int spline_sources_interpolate_growing_closeby(
			     double * x_array,
			     int tau_size,
			     double *** y_array,
			     double *** ddy_array,
			     int tp_size,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to tp_size */
           int index_mode,
           int index_ic,           
           int index_k,
           int k_size,
			     ErrorMsg errmsg
			     ) {

  int inf,sup;
  double h,a,b;

  inf = *last_index;
  class_test(inf<0 || inf>(tau_size-1),
      errmsg,
      "*lastindex=%d out of range [0:%d]\n",inf,tau_size-1);
      
  // *** Look at the left of last_index
  while (x < x_array[inf]) {
    inf--;
    if (inf < 0) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
       x,x_array[0]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  // *** Look at the right of last_index
  while (x > x_array[sup]) {
    sup++;
    if (sup > (tau_size-1)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
       x,x_array[tau_size-1]);
      return _FAILURE_;
    }
  }
  inf = sup-1;
  
  // *** Store the position for later use
  *last_index = inf;

  // *** Actual interpolation
  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  int index_tp;
  for (index_tp=0; index_tp<result_size; index_tp++) {
    
    // Cubic interpolation, see eq. 3.3.3 of Numerical Recipes in C.
    // To have plain linear interpolation, comment the last two lines.
    result[index_tp] = 
      a * y_array[index_mode][index_ic*tp_size + index_tp][inf*k_size+index_k] +
      b * y_array[index_mode][index_ic*tp_size + index_tp][sup*k_size+index_k] +
      ((a*a*a-a)* ddy_array[index_mode][index_ic*tp_size + index_tp][inf*k_size+index_k] + 
       (b*b*b-b)* ddy_array[index_mode][index_ic*tp_size + index_tp][sup*k_size+index_k])*h*h/6.;

  }

  return _SUCCESS_;
}




 /**
  * Takes the same input as 'array_interpolate_spline' but does linear interpolation, by
  * ignoring the spline matrix. Useful for debugging the splines thoroughout CLASS.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_fake(
			     double * x_array,
			     int n_lines,
			     double * array,
			     double * array_splined,
			     int n_columns,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to n_columns */
			     ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double h,a,b;
  
  inf=0;
  sup=n_lines-1;
  
  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    if (x > x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) = 
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i);

  return _SUCCESS_;
}



 
// ============================================================================
// =                                 Sampling related                         =
// ============================================================================ 
 

/**
 * Fill an array with logarithmically spaced points. The array should already be
 * allocated with 'n_points' doubles
 */
int log_space (double * x, double x_min, double x_max, int n_points) {

	int j;

	x[0] = x_min;
  x[n_points-1]=x_max; 
	double step = pow (x_max/x_min, 1.0/(n_points-1));	

  for (j=1; j<n_points-1; ++j)  
    x[j] = x[j-1]*step;    

	return _SUCCESS_;
	
}



/**
 * Fill an array with linearly spaced points. The array should already be
 * allocated with 'n_points' doubles
 */
int lin_space (double * x, double x_min, double x_max, int n_points) {

	int j;

	x[0]=x_min;
  x[n_points-1]=x_max;
	double step = (x_max-x_min)/(n_points-1);

	for (j=1; j<n_points-1; ++j)
	    x[j] = x_min + j*step;
	    
	return _SUCCESS_;
		
}




// ===================================================================================
// =                                Matrix operations                                =
// ===================================================================================

/* Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html) */


/**
 * Recursive definition of determinate using expansion by minors.
 * Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
 */
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   }
   else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   }
   else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   }
   else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         // det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         det += ALTERNATING_SIGN(j1+2) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}


/*
   Find the cofactor matrix of a square matrix
   Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         // b[i][j] = pow(-1.0,i+j+2.0) * det;
         b[i][j] = ALTERNATING_SIGN(i+j+2) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}


/*
   Transpose of a square matrix, do it in place
   Credits to Christopher M. Brown (http://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html)
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}


/** 
 * Inverse of a square matrix. The inversion can be done in-place, that is, the
 * the user can select in=out.
 */
void InverseMatrix(double **in,int n,double **out)
{

  int i,j;

  /* Obvious case where n=1 */
  if (n==1) {
    out[0][0] = 1/in[0][0];
    return;
  }
  
  /* Compute determinant */
  double det = Determinant(in, n);

  /* Compute cofactor matrix */
  double ** cofactor = malloc(n*sizeof(double *));
  for (i=0;i<n;i++)
    cofactor[i] = malloc(n*sizeof(double));
  CoFactor(in, n, cofactor);

  /* Compute inverse matrix (note that we take the transpose of cofactor) */
  for (i=0;i<n;i++)
     for (j=0;j<n;j++)
        out[i][j] = cofactor[j][i]/det;
  
  for (i=0;i<n;i++)
     free(cofactor[i]);
  free(cofactor);

}


/** 
 * Print a matrix to standard output.
 */
void PrintMatrix(double **in,int n)
{

  int i,j;

  for (int i = 0; i < n; ++i) {
    printf ("( ");
    for (int j = 0; j < n; ++j)
      printf ("%18.10e ", in[i][j]);
    printf (")\n");
  }

  return;

}






// ======================================================================
// =                          Assert functions                          =
// ======================================================================

/** 
 * Check whether a triad of numbers (l1,l2,l3) satisfies the triangular condition
 * |l1-l2| <= l3 <= l1+l2 (integer input)
 */
int is_triangular_int (int l1, int l2, int l3) {
  
  if ((l3>=abs(l1-l2)) && (l3<=(l1+l2)))
    return _TRUE_;
  else
    return _FALSE_;  
  
}


/** 
 * Check whether a triad of numbers (l1,l2,l3) satisfies the triangular condition
 * |l1-l2| <= l3 <= l1+l2 (floating point input)
 */
int is_triangular_double (double l1, double l2, double l3) {

  if ((l3>=fabs(l1-l2)) && (l3<=(l1+l2)))
    return _TRUE_;
  else
    return _FALSE_;  

}



// ========================================================================
// =                                 Misc                                 =
// ========================================================================

/** 
 * Identity function.
 */  
double identity_double (double x) {
  return x;
}

/**
 * Return +1 if positive, -1 if negative, 0 if zero.
 */
int sign_int (int x) {

  if (x>0) {
    return +1;
  }
  else if (x<0) {
    return -1;
  }
  else {
    return 0;
  }
}

/**
 * Given three numbers, return their list that orders them in
 * ascending order. For example, 
 * ordering(3,1,2) = (2,3,1).
 * The input array (n) and the output array (ordering) must
 * be pre-allocated with 4 values, as we start counting from 1.
 */
int ordering_int (
      int * n,            /* In */
      int * ordering,     /* Out */
      ErrorMsg errmsg
      )
{
  
  if (n[1] >= n[2]) {
    if (n[2] >= n[3]) {
      ordering[1] = 3;
      ordering[2] = 2;
      ordering[3] = 1;
    }
    else /* if n[3] > n[2] */ {
      if (n[1] >= n[3]) {
        ordering[1] = 2;
        ordering[2] = 3;
        ordering[3] = 1;
      }
      else /* if n[3] > n[1] */ {
        ordering[1] = 2;
        ordering[2] = 1;
        ordering[3] = 3;
      }
    }
  }
  else /* if n[2] > n[1] */ {
    if (n[1] >= n[3]) {
      ordering[1] = 3;
      ordering[2] = 1;
      ordering[3] = 2;
    }
    else /* if n[3] > n[1] */ {
      if (n[2] >= n[3]) {
        ordering[1] = 1;
        ordering[2] = 3;
        ordering[3] = 2;
      }
      else /* if n[3] > n[2] */ {
        ordering[1] = 1;
        ordering[2] = 2;
        ordering[3] = 3;
      }
    }
  }
  
  return _SUCCESS_;
  
}


/**
 * Reorder a vector of three elements according to an ordering vector.
 * The latter can be the output of 'ordering_int'.
 * For example, 
 * reorder( (4,7,3), (2, 1, 3) ) = (2,3,1).
 * The modification is made in-place. Both arrays must be pre-allocated
 * with 4 values, as we start counting from 1.
 */
int reorder_int (
      int * n,             /* In/Out */
      int * ordering,      /* In */
      ErrorMsg errmsg
      )
{
  
  /* Check that the ordering vector is valid */
  class_test ((ordering[1]*ordering[2]*ordering[3]!=6)
  || (ordering[1]+ordering[2]+ordering[3]!=6)
  || (MAX(ordering[1],MAX(ordering[2],ordering[3]))!=3)
  || (MIN(ordering[1],MIN(ordering[2],ordering[3]))!=1),
   errmsg, 
   "ordering vector not valid")
  
  /* Temporary array to store the values of n */
  double n_[] = {0, n[1], n[2], n[3]};
   
  /* Reorder 'n' */
  n[1] = n_[ordering[1]];
  n[2] = n_[ordering[2]];
  n[3] = n_[ordering[3]];
  
  return _SUCCESS_;
  
}
  
  
  
  
  
  
  
  
  
  
  


