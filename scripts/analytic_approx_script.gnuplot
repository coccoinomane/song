# Arguments
# root: should be in the format path/file_name without extension.
#       example: -e "root='/Users/coccoinomane/data/song/runs_from_honda/r+z/wmap7_better_t+p/squeezed_small_6_intrinsic'"

# Determine parameters
types = "ttt tte tet ett eet ete tee eee"

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
