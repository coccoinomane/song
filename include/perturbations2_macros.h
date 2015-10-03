/** @file perturbations2.h
 * 
 * Shortcuts for the perturbations2 module.
 *
 * The second-order equations evolved in the perturbations2.c module are
 * lengthy. In this file we define preprocessor macros that allow
 * us to write them down in SONG in a human readable way.
 *
 * For example, using the macros defined below we can write the polarisation
 * equations for m=0,1,2 as
 *
 *   dE(2,m) = -k * (d_plus(2,m,m)*E(3,m) + d_zero(2,m,m)*B(2,m))
 *             - kappa_dot*(E(2,m) + sqrt_6*Pi) ,
 *
 * which is a form very similar to the "published" version in eq. 4.146 of
 * http://arxiv.org/abs/1405.2280.
 * 
 */

#ifndef __PERTURBATIONS2_MACROS__
#define __PERTURBATIONS2_MACROS__


// ------------------------------------------------------------------------------------
// -                              Multipoles shorthands                               -
// ------------------------------------------------------------------------------------

/**
 * Shorthand to access the ppt2->sources array.
 *
 * It can be used only when ppw2 is defined and after perturb2_geometrical_corner()
 * and perturb2_get_k_lists() have been called.
 */
#define sources(index_type) ppt2->sources[(index_type)]\
                            [ppw2->index_k1]\
                            [ppw2->index_k2]\
                            [index_tau*ppt2->k3_size[ppw2->index_k1][ppw2->index_k2] + ppw2->index_k3]


/**
 * Shorthand used to index the first level of the following arrays for massless species:
 * ppt2->sources, pv->y and pv->dy.
 *
 * For example, the (l,m) multipole of the photon hierarchy is found at:
 *
 * ppw2->pv->y[ppw2->pv->index_pt2_monopole_g + lm(l,m)]
 * ppw2->pv->dy[ppw2->pv->index_pt2_monopole_g + lm(l,m)]
 *
 * while the line of sight source for the (l,m) multipole of the E-modes is in
 *
 * ppt2->sources[ppt2>index_tp2_E + lm(l,m)].
 */
#define lm(l,m) ppt2->lm_array[l][ppr2->index_m[m]]

/**
 * Shorthand used to index the first level of the pv->y and pv->dy arrays for the massive
 * species (baryons and cold dark matter).
 *
 * For example, the (n,l,m) moment of the baryon hierarchy is found at:
 *
 * ppw2->pv->y[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)]
 * ppw2->pv->dy[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)]
 *
 * NB: n and l cannot exceed 2 because we only consider the n=0,1,2 moments for
 * the massive species.
 */
#define nlm(n,l,m) ppt2->nlm_array[n][l][ppr2->index_m[m]]

/**
 * Shorthand used to index ppw2->rotation_1 and ppw2->rotation_2, the arrays required to
 * to compute the quadratic sources in the desired values of \vec{k1} and \vec{k2}.
 *
 * We cannot use lm(l,m) for this purpose because the range in (l,m) is larger for
 * the quadratic sources than for the second-order perturbations; for more details,
 * see perturb2_get_lm_lists().
 */
#define lm_quad(l,m) ppt2->lm_array_quad[l][m]


/**
 * Shorthands to access y, the vector of evolved perturbations.
 *
 * Using these macros allows to lighten the notation of the second-order equations,
 * especially in perturb2_derivs(). Make sure to use the macros only in those functions
 * that have access to dy.
 * 
 * The macros will return zero if the requested multipole is not being evolved. For
 * example, E(l,m) and B(l,m) will always return zero if polarisation is turned off,
 * and all photon moments with l>1 will return zero if the tight coupling approximation
 * is turned on. This behaviour allows more flexibility in writing the equations,
 * while avoiding segmentation faults. For example, it allows to write the equations for
 * the photon intensity hierarchy in the same way regardless of whether polarization is
 * on or off.
 * 
 * If an approximation is active, some of the macros will not be available. For
 * example, in the tight coupling regime, the photon quadrupole is stored in
 * ppw2->I_2m[m] rather than in I(2,m), because it is not evolved. Make sure to use
 * the following arrays, which are updated in perturb2_workspace_at_tau(), rather
 * than the corresponding macros:
 * 
 * - ppw2->I_00         photon intensity monopole (depends on RSA)
 * - ppw2->N_00         neutrino dipole (depends on RSA)
 * - ppw2->N_1m[m]      neutrino dipole (depends on RSA)
 * - ppw2->I_1m[m]      photon intensity dipole (depends on TCA and RSA)
 * - ppw2->I_2m[m]      photon  intensity quadrupole (depends on TCA)
 * - ppw2->E_2m[m]      E-polarisation quadrupole (depends on TCA)
 * - ppw2->B_2m[m]      B-polarisation quadrupole (depends on TCA)
 * - ppw2->C_1m[m]      Dipole collision term: kappa_dot * (4/3*b_11m-I_1m) (depends on TCA and RSA)
 * - ppw2->b_200        baryon "pressure" (depends on perfect fluid approximation)
 * - ppw2->b_22m[m]     baryon quadrupole (depends on perfect fluid approximation)
 * - ppw2->cdm_200      CDM "pressure" (depends on perfect fluid approximation)
 * - ppw2->cdm_22m[m]   CDM quadrupole (depends on perfect fluid approximation)
 * 
 * NB: These macros cannot be used to write on y, but only to access it.
 */
//@{
#define I(l,m) ( (((l)>ppw2->pv->l_max_g) || (abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_g + lm(l,m)] )
#define E(l,m) ( ((ppt2->has_polarization2 == _FALSE_)||((l)>ppw2->pv->l_max_pol_g)||(abs(m)>l)||(l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_E + lm(l,m)] )
#define B(l,m) ( ((ppt2->has_polarization2 == _FALSE_)||((l)>ppw2->pv->l_max_pol_g)||(abs(m)>l)||(l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_B + lm(l,m)] )
#define N(l,m) ( (((l)>ppw2->pv->l_max_ur) || (abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_ur + lm(l,m)] )                              
#define b(n,l,m) ( (((n)>ppw2->pv->n_max_b) || ((l)>ppw2->pv->l_max_b) || (abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)] )
#define cdm(n,l,m) ( (((n)>ppw2->pv->n_max_cdm) || ((l)>ppw2->pv->l_max_cdm) || (abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_cdm + nlm(n,l,m)] )                              
//@}
                                                         
/**
 * Shorthands to access and write dy, the vector with the time derivatives of the evolved
 * perturbations.
 *
 * Using these macros allows to lighten the notation of the second-order equations,
 * especially in perturb2_derivs(). Make sure to use the macros only in those functions
 * that have access to dy.
 */
//@{
#define dI(l,m) dy[ppw2->pv->index_pt2_monopole_g + lm(l,m)]
#define dE(l,m) dy[ppw2->pv->index_pt2_monopole_E + lm(l,m)]
#define dB(l,m) dy[ppw2->pv->index_pt2_monopole_B + lm(l,m)]
#define dN(l,m) dy[ppw2->pv->index_pt2_monopole_ur + lm(l,m)]
#define db(n,l,m) dy[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)]
#define dcdm(n,l,m) dy[ppw2->pv->index_pt2_monopole_cdm + nlm(n,l,m)]
//@}               



// ------------------------------------------------------------------------------------
// -                           Quadratic sources shorthands                           -
// ------------------------------------------------------------------------------------


/**
 * Shorthands to access and write the quadratic sources of the Boltzmann equation.
 */
//@{
#define dI_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_g + lm(l,m)]
#define dE_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_E + lm(l,m)]
#define dB_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_B + lm(l,m)]
#define db_qs2(n,l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_b + nlm(n,l,m)]
#define dcdm_qs2(n,l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_cdm + nlm(n,l,m)]
#define dN_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_ur + lm(l,m)]
//@}


/**
 * Shorthands to access and write the collisional quadratic sources of the Boltzmann
 * equation.
 */
//@{
#define dI_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_g + lm(l,m)]
#define dE_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_E + lm(l,m)]
#define dB_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_B + lm(l,m)]
#define db_qc2(n,l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_b + nlm(n,l,m)]
//@}


/**
 * Shorthand to access the rotation coefficients.
 *
 * The rotation coefficients are used to obtain the first-order perturbations in an
 * arbitrary k configuration, starting from the perturbations computed in the
 * perturbations.c module, which have k aligned with the z-axis, according to the
 * rotation formula
 *
 *   X(\vec{k1})_lm = rot_1(l,m) X_raw(k1)_l0
 *   X(\vec{k2})_lm = rot_2(l,m) X_raw(k2)_l0,
 *
 * where X_raw(k1)_l0 is the multipole as computed by the perturbations.c module.
 *
 * The rotation coefficients rot_1(l,m) and rot_2(l,m) are basically the associated
 * Legendre polynomials P_lm(cos\theta_1) and P_lm(cos\theta_2). For more details on
 * them, refer to Sec. B.1 of http://arxiv.org/abs/1405.2280.
 *
 * The macros will return zero if accessed with l<0 and abs(m)>l; this will give us
 * the freedom to write the equations inside a loop over (l,m).
 *
 * NB: These macros can be used only after perturb2_geometrical_corner() has been
 * called.
 */
//@{
#define rot_1(l,m) ( (((l)<0)||(abs(m)>(l))) ? 0 : ( (m)<0 ? ppw2->rotation_1_minus[lm_quad(l,abs(m))] : ppw2->rotation_1[lm_quad(l,m)]))
#define rot_2(l,m) ( (((l)<0)||(abs(m)>(l))) ? 0 : ( (m)<0 ? ppw2->rotation_2_minus[lm_quad(l,abs(m))] : ppw2->rotation_2[lm_quad(l,m)]))
//@}


/**
 * Shorthands to access the first-order massless species, as outputted by the
 * first-order system.
 * 
 * The suffix _raw highlights that these perturbations were computed with the k-vector
 * aligned with the symmetry axis of the spherical harmonics.
 *
 * The perturbations in an arbitrary k-direction is obtained by multiplying the raw
 * perturbation with the rotation coefficients. For example,
 * 
 * I_1(l,m) = rot_1(l,m) * I_1_raw(l).
 *
 * NB: These macros can be used only if the vectors with the first order perturbations,
 * ppw2->pvec_sources1 and ppw2->pvec_sources2, have been filled with a previous call
 * to perturb_song_sources_at_tau().
 */
//@{
#define I_1_raw(l) ( l < 0 ? 0 : pvec_sources1[ppt->index_qs_monopole_g + l] )
#define I_2_raw(l) ( l < 0 ? 0 : pvec_sources2[ppt->index_qs_monopole_g + l] )
#define E_1_raw(l) ( ((l < 0)||(ppt2->has_polarization2==_FALSE_)) ? 0 : pvec_sources1[ppt->index_qs_monopole_E + l] )
#define E_2_raw(l) ( ((l < 0)||(ppt2->has_polarization2==_FALSE_)) ? 0 : pvec_sources2[ppt->index_qs_monopole_E + l] )
#define N_1_raw(l) ( l < 0 ? 0 : pvec_sources1[ppt->index_qs_monopole_ur + l] )
#define N_2_raw(l) ( l < 0 ? 0 : pvec_sources2[ppt->index_qs_monopole_ur + l] )
//@}


/**
 * Shorthands to access the first-order massless species.
 *
 * These are obtained from the "raw" perturbations by applying the rotation coefficients
 * rot_1(l,m) and rot_2(l,m). See the description of rot1 and rot2 for details.
 *
 * NB: These macros can be used only if the vectors with the first order perturbations,
 * ppw2->pvec_sources1 and ppw2->pvec_sources2, have been filled with a previous call
 * to perturb_song_sources_at_tau().
 */
//@{
#define I_1(l,m) (rot_1(l,m))*(I_1_raw(l))
#define I_2(l,m) (rot_2(l,m))*(I_2_raw(l))
#define E_1(l,m) (rot_1(l,m))*(E_1_raw(l))
#define E_2(l,m) (rot_2(l,m))*(E_2_raw(l))
#define N_1(l,m) (rot_1(l,m))*(N_1_raw(l))
#define N_2(l,m) (rot_2(l,m))*(N_2_raw(l))
//@}


/* NOT IMPLEMENTED YET!
Set this coefficient to 2 if you expanded the perturbations as X ~ X^(1) + 1/2 * X^(2),
or to unity if you use X ~ X^(1) + X^(2) instead. This feature is not fully implented
yet, hence for the time being keep it equal to 2. */
#define quad_coefficient 2


// ------------------------------------------------------------------------------------
// -                              Coupling coefficients                               -
// ------------------------------------------------------------------------------------

/**
 * Shorthands for the coupling coefficients C and D.
 * 
 * See documentation in perturbations2.h.
 */
//@{
#define c_minus(l,m1,m) ppt2->c_minus[lm(l,m)][m-m1+1]
#define c_plus(l,m1,m) ppt2->c_plus[lm(l,m)][m-m1+1]
#define d_minus(l,m1,m) ppt2->d_minus[lm(l,m)][m-m1+1]
#define d_plus(l,m1,m) ppt2->d_plus[lm(l,m)][m-m1+1]
#define d_zero(l,m1,m) ppt2->d_zero[lm(l,m)][m-m1+1]
//@}


/**
 * Shorthands to access and write the summed coupling coefficients.
 *
 * See documentation in perturbations2.h.                              
 */
//@{
#define c_minus_12(l,m) ppw2->c_minus_product_12[lm(l,m)]
#define c_minus_21(l,m) ppw2->c_minus_product_21[lm(l,m)]
#define c_plus_12(l,m)  ppw2->c_plus_product_12[lm(l,m)]  
#define c_plus_21(l,m)  ppw2->c_plus_product_21[lm(l,m)]  
#define c_minus_11(l,m) ppw2->c_minus_product_11[lm(l,m)]
#define c_minus_22(l,m) ppw2->c_minus_product_22[lm(l,m)]
#define c_plus_11(l,m)  ppw2->c_plus_product_11[lm(l,m)]  
#define c_plus_22(l,m)  ppw2->c_plus_product_22[lm(l,m)]  

#define r_minus_12(l,m) ppw2->r_minus_product_12[lm(l,m)]
#define r_minus_21(l,m) ppw2->r_minus_product_21[lm(l,m)]
#define r_plus_12(l,m)  ppw2->r_plus_product_12[lm(l,m)]  
#define r_plus_21(l,m)  ppw2->r_plus_product_21[lm(l,m)] 
 
#define d_minus_12(l,m) ppw2->d_minus_product_12[lm(l,m)]
#define d_minus_21(l,m) ppw2->d_minus_product_21[lm(l,m)]
#define d_plus_12(l,m)  ppw2->d_plus_product_12[lm(l,m)]  
#define d_plus_21(l,m)  ppw2->d_plus_product_21[lm(l,m)]
#define d_minus_11(l,m) ppw2->d_minus_product_11[lm(l,m)]
#define d_minus_22(l,m) ppw2->d_minus_product_22[lm(l,m)]
#define d_plus_11(l,m)  ppw2->d_plus_product_11[lm(l,m)]  
#define d_plus_22(l,m)  ppw2->d_plus_product_22[lm(l,m)]  

#define d_zero_12(l,m)  ppw2->d_zero_product_12[lm(l,m)]  
#define d_zero_21(l,m)  ppw2->d_zero_product_21[lm(l,m)]  
#define d_zero_11(l,m)  ppw2->d_zero_product_11[lm(l,m)]  
#define d_zero_22(l,m)  ppw2->d_zero_product_22[lm(l,m)]  

#define k_minus_12(l,m) ppw2->k_minus_product_12[lm(l,m)]
#define k_minus_21(l,m) ppw2->k_minus_product_21[lm(l,m)]
#define k_plus_12(l,m)  ppw2->k_plus_product_12[lm(l,m)]  
#define k_plus_21(l,m)  ppw2->k_plus_product_21[lm(l,m)]  
#define k_minus_11(l,m) ppw2->k_minus_product_11[lm(l,m)]
#define k_minus_22(l,m) ppw2->k_minus_product_22[lm(l,m)]
#define k_plus_11(l,m)  ppw2->k_plus_product_11[lm(l,m)]  
#define k_plus_22(l,m)  ppw2->k_plus_product_22[lm(l,m)]  

#define k_zero_12(l,m)  ppw2->k_zero_product_12[lm(l,m)]  
#define k_zero_21(l,m)  ppw2->k_zero_product_21[lm(l,m)]  
#define k_zero_11(l,m)  ppw2->k_zero_product_11[lm(l,m)]  
#define k_zero_22(l,m)  ppw2->k_zero_product_22[lm(l,m)]
//@}


// ------------------------------------------------------------------------------------
// -                                 Debug shortcuts                                  -
// ------------------------------------------------------------------------------------
                              
#define fprintf_k_debug(args...) {                                  \
  if((ppt2->k_out_size > 0) &&                                      \
     (ppw2->index_k1==ppt2->index_k1_out[0]) &&                     \
     (ppw2->index_k2==ppt2->index_k2_out[0]) &&                     \
     (ppw2->index_k3==ppt2->index_k3_out[0])) {                     \
    fprintf (args);                                                 \
  }                                                                 \
}

#define printf_k_debug(args...) {                                   \
  fprintf_k_debug (stdout, args)                                    \
}

#define p1(x) {                                             \
  fprintf_k_debug (stderr, "%25.17g\n", x);                            \
}

#define p2(x, y) {                                         \
  fprintf_k_debug (stderr, "%25.17g %25.17g\n", x, y);                    \
}

#define p3(x, y, z) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g\n", x, y, z);            \
}

#define p4(a,b,c,d) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g\n", a, b, c, d);            \
}

#define p5(a,b,c,d,e) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e);            \
}

#define p6(a,b,c,d,e,f) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f);            \
}

#define p7(a,b,c,d,e,f,g) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g);            \
}

#define p8(a,b,c,d,e,f,g,h) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g, h);            \
}

#define p9(a,b,c,d,e,f,g,h,i) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g, h, i);            \
}

#define p10(a,b,c,d,e,f,g,h,i,l) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g, h, i, l);            \
}

#define p11(a,b,c,d,e,f,g,h,i,l,m) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g, h, i, l, m);            \
}

#define p12(a,b,c,d,e,f,g,h,i,l,m,n) {                                     \
  fprintf_k_debug (stderr, "%25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g %25.17g\n", a, b, c, d, e, f, g, h, i, l, m, n);            \
}



#endif