#! /bin/sh

l_long=$1
shift
folders=$@

make print_bispectra -j4;

for folder in $folders; do
  for bf in ttt tte tet ett eet ete tee eee; do
    for bt in "intrinsic"; do # add bispectra types here 
      type_bispectrum="${bt}_${bf}"
      cmd="time ./print_bispectra $type_bispectrum squeezed_small_scale $l_long $folder 2> $folder/squeezed_small_${l_long}_${type_bispectrum}.txt"
      echo $cmd
      eval $cmd
    done
  done
  echo "Generating PDF file with squeezed approximation..."
  gnuplot -e "root='$folder/squeezed_small_${l_long}_intrinsic'" scripts/analytic_approx_script.gnuplot
done