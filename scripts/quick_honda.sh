#! /bin/sh

threads=$1
shift
folders=$@

root="~/guest/data/song/runs"
make clean;
make -f makefile_honda -j8 honda > /dev/null;
export OMP_NUM_THREADS=$threads;

for folder in $folders; do
  cmd="unbuffer time ./honda $root/$folder 2>&1 | tee $root/$folder/log.txt"
  echo "$cmd"
  eval "$cmd"
done
