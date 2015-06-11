# Usage:
#
#   gnuplot analytic_approx_script.gnuplot -e "root='/<run_dir>/<filename>'"
#
# The filename is the common root of the files to be printed, without extension.
#
# For example, 
#
#   gnuplot -e "root='test_run/squeezed_small_6/squeezed_small_6_intrinsic'" analytic_approx_script.gnuplot 
#
# will plot the following files:
#
#   squeezed_small_6_intrinsic_ttt.txt
#   squeezed_small_6_intrinsic_tte.txt
#   squeezed_small_6_intrinsic_tet.txt
#   squeezed_small_6_intrinsic_ett.txt
#   squeezed_small_6_intrinsic_eet.txt
#   squeezed_small_6_intrinsic_ete.txt
#   squeezed_small_6_intrinsic_tee.txt
#   squeezed_small_6_intrinsic_eee.txt

set macros

# Determine parameters
types = "ttt tte tet ett eet ete tee eee"
l_min = 6
l_max = 2000
frac_min = 1e-3
frac_max = 10

# Bispectrum column
bcol = "($9)"

# Uncomment if $9 is the column with the bolometric temperature
# bispectrum and $14 and $15 are the columns with the temperature
# and redshift quadratic corrections, respectively
# bcol = "($9-$14-$15)"

# Analytic approximation column
acol = "($10)"


# ==========================================
# = PLOT INTRINSIC VS ANALYTICAL BISPECTRA =
# ==========================================

unset multiplot

# Uncomment to save to file
set terminal pdf size 9in,6in
set output root.".pdf"

# Plot
unset logscale
set multiplot layout 4,2
do for [bf in types]  {
  file = root."_".bf.".txt"
  plot [l_min:l_max] file u 1:@bcol w li lw 4 title bf, "" u 1:@acol w li lw 4 title "approx"
}
unset multiplot
set output

# ============================================
# = PLOT FRACTIONAL DIFFERENCE IN ONE FIGURE =
# ============================================

# Uncomment to save to file
set terminal pdf size 9in,6in
set output root."_frac.pdf"

# Exclude redundant bispectra
types = "ttt tte ett eet tee eee"

# Plot
set logscale
set multiplot layout 3,2
do for [bf in types]  {
  file = root."_".bf.".txt"
  plot [l_min:l_max] [frac_min:frac_max] file u 1:(abs(1-(@bcol/@acol))) w li lw 4 title bf
}
unset multiplot
set output

# ==========================================
# = PLOT FRACTIONAL DIFFERENCE IN ONE PLOT =
# ==========================================

# Uncomment to save to file
set terminal pdf size 9in,6in
set output root."_frac_oneplot.pdf"

# Exclude redundant bispectra
types = "ttt tte ett eet tee eee"

# Plot
set logscale
plot [l_min:l_max] [frac_min:frac_max] for [bf in types] root."_".bf.".txt" u 1:(abs(1-(@bcol/@acol))) w li lw 4 title bf
set output
