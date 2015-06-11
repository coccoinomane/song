# Print to file squeezed configurations of the bispectra, fixing
# the smallest multipole, for all 8 temperature and polarisation
# combinations.
#
# usage:    ./script.sh <l> <run_directory>
#
# The bispectra will be read from the run directory, and the output
# files will be stored in the folder run_directory/squeezed_small_<l>.
#
# Each output file will contain the bispectrum configurations with
# L1 fixed to <l> and with varying values of L2=L3.
#
# For example,
# 
#   ./script.sh 6 test_run
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

#! /bin/sh

l_long=$1
shift
runs=$@

make print_bispectra -j4 > /dev/null;

for run in $runs; do
  out_dir=$run/squeezed_small_${l_long}
  mkdir $out_dir
  for bf in ttt tte tet ett eet ete tee eee; do
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
