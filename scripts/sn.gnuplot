set terminal pdf size 7in,4in

# =========================
# = SIGNAL TO NOISE IDEAL =
# =========================

root_t = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/wmap7_better_t/"
root_e = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/wmap7_better_p/"
root_te = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/L_convergence/wmap7_better_t+p_302L2000/"

do for [file in "lmax lmin"] {
  do for [bt in "intrinsic local equilateral"] {
    do for [label in "SN SN_cum"] {
      set title "S/N of ".bt." bispectrum as a function of ".file offset screen 0.25,0 font "Times-Roman,12"
      set output label."_".bt."_".file.".pdf"
      set multiplot layout 1,2
        unset logscale; set logscale x;
        plot [2:2000] [1e-4:] root_t."fisher_".file.".dat" u 1:(column(label."_".bt)*$1) w li lw 5 t "TTT",\
                              root_e."fisher_".file.".dat" u 1:(column(label."_".bt)*$1) w li lw 5 t "EEE",\
                              root_te."fisher_".file.".dat" u 1:(column(label."_".bt)*$1) w li lw 5 t "T+E, 8 bispectra"
        set logscale
        set title " "
        replot
      unset multiplot
    }
  }
}

# =======
# = CLS =
# =======

set output "cls.pdf"  
set multiplot layout 1,2
  unset title
  set logscale
  plot [2:2000] root_te."cl.dat" u 1:2 w li lw 5
  plot [2:2000] root_te."cl.dat" u 1:3 w li lw 5
unset multiplot

# ===============================
# = SIGNAL TO NOISE EXPERIMENTS =
# ===============================

root_planck = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/wmap7_ref_t+p_PlanckNoise2500/"
root_core = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/wmap7_ref_t+p_Core3000/"
root_prism = "/Users/coccoinomane/data/song/runs_from_honda/r+z/pre_polarised_redshift_term/r+z/wmap7_ref_t+p_PrismNoise3000/"

do for [file in "lmax"] {
  set title "S/N as a function of ".file offset 0,0 font "Times-Roman,20"
  set output "sn_".file."_experiments.pdf"
  set multiplot layout 1,2
    unset logscale
    plot [2:3000] [1e-4:] root_planck."fisher_".file.".dat" u 1:(column("SN_cum_intrinsic")*$1) w li lw 5 t "Planck T+E",\
                          root_core."fisher_".file.".dat" u 1:(column("SN_cum_intrinsic")*$1) w li lw 5 t "Core T+E",\
                          root_prism."fisher_".file.".dat" u 1:(column("SN_cum_intrinsic")*$1) w li lw 5 t "Prism T+E"
    set logscale
    replot
  unset multiplot
}

set terminal wxt