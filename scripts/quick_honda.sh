#! /bin/sh

threads=$1
shift
folders=$@

root="~/guest/data/song/runs"
make clean;
make -f makefile_honda -j8 honda > /dev/null;
export OMP_NUM_THREADS=$threads;

for folder in $folders; do
  cmd="time ./honda $root/$folder 2> $root/$folder/fisher.txt"
  echo $cmd
  eval $cmd
done
