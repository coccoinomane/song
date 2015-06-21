# Print to file squeezed configurations of the bispectra, fixing
# the smallest multipole, for all 8 temperature and polarisation
# combinations.
#
# usage:    ./bisp_squeezed.sh <l> <run_directory_list>
#
# The bispectra will be read from the run directory, and the output
# files will be stored in the folder run_directory/squeezed_small_<l>.
#
# Each output file will contain the bispectrum configurations with
# L1 fixed to <l> and with varying values of L2=L3.
#
# For example,
# 
#   ./bisp_squeezed.sh 6 test_run
#
# will create the following 8 files in the folder test_run/squeezed_small_6:
#
#   squeezed_small_6_ttt.txt
#   squeezed_small_6_tte.txt
#   squeezed_small_6_tet.txt
#   squeezed_small_6_ett.txt
#   squeezed_small_6_eet.txt
#   squeezed_small_6_ete.txt
#   squeezed_small_6_tee.txt
#   squeezed_small_6_eee.txt
#
# Modify the 'probes' variable below to choose a different set of bispectra to
# print.
#
# If you give a list of folders in <run_directory>, the script will
# process each of them.

#! /bin/sh

# Which bispectra to print?
probes = "ttt tte tet ett eet ete tee eee"

# Parse arguments
l_long=$1
shift
runs=$@

# Compile SONG
make print_bispectra -j4 > /dev/null;

# Print bispectra to file
for run in $runs; do
  out_dir=$run/squeezed_small_${l_long}
  mkdir $out_dir
  for bf in $probes; do
    for bt in "intrinsic"; do # add bispectra types here 
      type_bispectrum="${bt}_${bf}"
      cmd="time ./print_bispectra $type_bispectrum squeezed_small_scale $l_long $run 2> $out_dir/squeezed_small_${l_long}_${type_bispectrum}.txt"
      echo $cmd
      eval $cmd
    done
  done
  echo "Generating PDF file with squeezed approximation..."
  gnuplot -e "root='$out_dir/squeezed_small_${l_long}_intrinsic'" scripts/analytic_approx_script.gnuplot
done
