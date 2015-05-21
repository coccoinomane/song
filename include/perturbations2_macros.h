// ------------------------------------------------------------------------------------
// -                          Shortcuts for 2nd-order moments                         -
// ------------------------------------------------------------------------------------

/* Shortcut to access the ppt2->sources array. It must be used only if the workspace ppw2 is
defined, and after the functions perturb2_geometrical_corner and perturb2_get_k_lists have
been called. */
#define sources(index_type) ppt2->sources[(index_type)]\
                            [ppw2->index_k1]\
                            [ppw2->index_k2]\
                            [index_tau*ppt2->k3_size[ppw2->index_k1][ppw2->index_k2] + ppw2->index_k3]


/* Define the photon temperature multipoles. Note that we set them to zero by default when
the no-radiation approximation is switched on. */
#define I(l,m) ( ((abs(m)>l) || (l<0) || ((l)>ppw2->pv->l_max_g)) ? 0 : y[ppw2->pv->index_pt2_monopole_g + lm(l,m)] )
#define dI(l,m) dy[ppw2->pv->index_pt2_monopole_g + lm(l,m)]

/* Define the photon polarization multipoles.  We include a check on the E and B modes, so
that they are non-zero only if polarization is requested.  This will allow us to write
the photon temperature hierarchy in the same way regardless of whether polarization is on
or off. */
#define E(l,m) ( ((ppt2->has_polarization2 == _FALSE_)||(abs(m)>l)||(l<0)||((l)>ppw2->pv->l_max_pol_g)) ? 0 : y[ppw2->pv->index_pt2_monopole_E + lm(l,m)] )
#define dE(l,m) dy[ppw2->pv->index_pt2_monopole_E + lm(l,m)]

#define B(l,m) ( ((ppt2->has_polarization2 == _FALSE_)||(abs(m)>l)||(l<0)||((l)>ppw2->pv->l_max_pol_g)) ? 0 : y[ppw2->pv->index_pt2_monopole_B + lm(l,m)] )
#define dB(l,m) dy[ppw2->pv->index_pt2_monopole_B + lm(l,m)]

/* Define the neutrinos multipoles */
#define N(l,m) ( ((abs(m) > l) || ((l)>ppw2->pv->l_max_ur)) ? 0 : y[ppw2->pv->index_pt2_monopole_ur + lm(l,m)] )
#define dN(l,m) dy[ppw2->pv->index_pt2_monopole_ur + lm(l,m)]

/* Define the baryon beta-moments */
#define b(n,l,m) ( ((abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)] )
#define db(n,l,m) dy[ppw2->pv->index_pt2_monopole_b + nlm(n,l,m)]

/* Define the CDM beta-moments */
#define cdm(n,l,m) ( ((abs(m)>l) || (l<0)) ? 0 : y[ppw2->pv->index_pt2_monopole_cdm + nlm(n,l,m)] )
#define dcdm(n,l,m) dy[ppw2->pv->index_pt2_monopole_cdm + nlm(n,l,m)]

// ------------------------------------------------------------------------------------
// -                          Shortcuts for quadratic sources                         -
// ------------------------------------------------------------------------------------

/* Set this coefficient to 2 if you expanded the perturbations as X ~ X^(1) + 1/2 * X^(2),
or to unity if you use X ~ X^(1) + X^(2) instead. This feature is not fully implented
yet, hence for the time being keep it equal to 2 */
#define quad_coefficient 2

/* Shortcuts to access the full quadratic sources, as computed by perturb2_quadratic_sources or
interpolated by perturb2_quadratic_sources_at_tau. */
#define dI_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_g + lm(l,m)]
#define dE_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_E + lm(l,m)]
#define dB_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_B + lm(l,m)]
#define db_qs2(n,l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_b + nlm(n,l,m)]
#define dcdm_qs2(n,l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_cdm + nlm(n,l,m)]
#define dN_qs2(l,m) ppw2->pvec_quadsources[ppw2->index_qs2_monopole_ur + lm(l,m)]

/* Same as above, but only for the collision term. These sources are never interpolated, and
the split is useful to include the quadratic part of the collision term in the line of
sight sources. */
#define dI_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_g + lm(l,m)]
#define dE_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_E + lm(l,m)]
#define dB_qc2(l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_B + lm(l,m)]
#define db_qc2(n,l,m) ppw2->pvec_quadcollision[ppw2->index_qs2_monopole_b + nlm(n,l,m)]

// ------------------------------------------------------------------------------------
// -                          Shortcuts for 1st-order moments                         -
// ------------------------------------------------------------------------------------

/* Rotation coefficients for the first-order quantities. They are defined as
sqrt(4pi/(2l+1)) Y_lm(theta,phi)
and they appear in the the rotation formula
Delta_lm(\vec{k1}) = sqrt(4pi/(2l+1)) Y_lm(theta,phi) Delta_l(k1)
where Delta_l is the first-order multipole computed with k1 aligned with the zenith.
                              
Note that We set the multipoles to zero when they are accessed with l<0 and abs(m)>l.
This will give us some freedom in writing the equations inside loops over (l,m).

IMPORTANT: These macros can be used inside perturb2_solve only and after perturb2_geometrical_corner
has been called. */
#define rot_1(l,m) ( (((l) < 0) || (abs(m)>(l))) ? 0 : ( (m) < 0 ? ppw2->rotation_1_minus[lm_quad(l,abs(m))] : ppw2->rotation_1[lm_quad(l,m)]))
#define rot_2(l,m) ( (((l) < 0) || (abs(m)>(l))) ? 0 : ( (m) < 0 ? ppw2->rotation_2_minus[lm_quad(l,abs(m))] : ppw2->rotation_2[lm_quad(l,m)]))

/* First order multipoles as computed by the first-order system. The 'tilde' suffix denotes
any 1st-order quantity computed with the symmetry axis aligned with the wavemode in its argument.
For example, I_1_tilde(l) is the (l,0) 1st-order moment in k1, computed with k1 aligned
with the z-axis. The E-modes are set to be different from zero only if polarization is
requested.  We have no need for the B-modes because they vanish at 1st order.

IMPORTANT: These macros can be used whenever a pvec_sources1 vector exists. This happens
either after a call to with perturb_song_sources_at_tau (as in perturb_source_term
or perturb2_save_early_transfers) or when pvec_sources1 and pvec_sources2 directly access the
quadratic sources stored in ppt->quadsources (as in perturb2_quadratic_sources). */
#define I_1_tilde(l) ( l < 0 ? 0 : pvec_sources1[ppt->index_qs_monopole_g + l] )
#define I_2_tilde(l) ( l < 0 ? 0 : pvec_sources2[ppt->index_qs_monopole_g + l] )

#define E_1_tilde(l) ( ((l < 0)||(ppt2->has_polarization2==_FALSE_)) ? 0 : pvec_sources1[ppt->index_qs_monopole_E + l] )
#define E_2_tilde(l) ( ((l < 0)||(ppt2->has_polarization2==_FALSE_)) ? 0 : pvec_sources2[ppt->index_qs_monopole_E + l] )

#define N_1_tilde(l) ( l < 0 ? 0 : pvec_sources1[ppt->index_qs_monopole_ur + l] )
#define N_2_tilde(l) ( l < 0 ? 0 : pvec_sources2[ppt->index_qs_monopole_ur + l] )

/* Rotated first-order multipoles.

IMPORTANT: These macros can be used whenever I_1_tilde can (see above) */
#define I_1(l,m) (rot_1(l,m))*(I_1_tilde(l))
#define I_2(l,m) (rot_2(l,m))*(I_2_tilde(l))
#define E_1(l,m) (rot_1(l,m))*(E_1_tilde(l))
#define E_2(l,m) (rot_2(l,m))*(E_2_tilde(l))
#define N_1(l,m) (rot_1(l,m))*(N_1_tilde(l))
#define N_2(l,m) (rot_2(l,m))*(N_2_tilde(l))

/* Define shorthands for the coupling coefficients. This is the only place where the arrays like
ppt2->c_minus are used */
#define c_minus(l,m1,m) ppt2->c_minus[lm(l,m)][m-m1+1]
#define c_plus(l,m1,m) ppt2->c_plus[lm(l,m)][m-m1+1]
#define d_minus(l,m1,m) ppt2->d_minus[lm(l,m)][m-m1+1]
#define d_plus(l,m1,m) ppt2->d_plus[lm(l,m)][m-m1+1]
#define d_zero(l,m1,m) ppt2->d_zero[lm(l,m)][m-m1+1]

/* Define shorthands for the summed coupling coefficients */
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
  
  






