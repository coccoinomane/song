#! /bin/sh
l1=$1
shift
l2=$1
shift
folders=$@

make print_bispectra -j4;

for folder in $folders; do
  for bf in ttt tte tet ett eet ete tee eee; do
    for bt in "local" "intrinsic"; do
      type_bispectrum="${bt}_${bf}"
      cmd="time ./print_bispectra $type_bispectrum triangular $l1 $l2 $folder 2> $folder/triangular_${l1}_${l2}_${type_bispectrum}.txt"
      echo $cmd
      eval $cmd
    done
  done
done
