# Arguments
# root: should be in the format path/file_name without extension.
#       example: gnuplot analytic_approx_script.gnuplot -e "root='<run dir>/squeezed_small_6_intrinsic'"

# Determine parameters
types = "ttt tte tet ett eet ete tee eee"

# ==========================================
# = PLOT INTRINSIC VS ANALYTICAL BISPECTRA =
# ==========================================

# Uncomment to save to file
set terminal pdf size 9in,6in
set output root.".pdf"

# Plot
unset logscale
set multiplot layout 4,2
do for [bf in types]  {
  file = root."_".bf.".txt"
  plot [:] file u 1:9 w li lw 4 title bf, "" u 1:10 w li lw 4 title "approx"
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
  plot [2:] [1e-5:100] file u 1:(abs(1-($9/$10))) w li lw 4 title bf
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
plot [2:] [1e-5:100] for [bf in types] root."_".bf.".txt" u 1:(abs(1-($9/$10))) w li lw 4 title bf
set output
