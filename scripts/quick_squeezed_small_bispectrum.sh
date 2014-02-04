#! /bin/sh

type_bispectrum=$1
shift
l_long=$1
shift
folders=$@

make print_bispectra -j4;

for folder in $folders; do
  cmd="time ./print_bispectra $type_bispectrum squeezed_small_scale $l_long $folder 2> $folder/squeezed_small_${l_long}_${type_bispectrum}.txt"
  echo $cmd
  eval $cmd
done
