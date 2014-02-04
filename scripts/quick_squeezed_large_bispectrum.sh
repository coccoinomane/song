#! /bin/sh

n_bispectrum=$1
shift
l_short=$1
shift
folders=$@

if [ $n_bispectrum != "0" ];
then
  suffix="_$n_bispectrum"
else
  suffix=""
fi

make -j2;

for folder in $folders; do
  cmd="time ./print_bispectra $n_bispectrum squeezed_large_scale $l_short $folder 2> $folder/squeezed_large_$l_short$suffix.txt"
  echo $cmd
  eval $cmd
done