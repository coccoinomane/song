#! /bin/sh

n_bispectrum=$1
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
  cmd="time ./print_bispectra $n_bispectrum equilateral $folder 2> $folder/equilateral$suffix.txt"
  echo $cmd
  eval $cmd
done
