#! /bin/sh

folders=$@

make -j2 class;

for folder in $folders; do
  cmd="time ./class $folder"
  echo $cmd
  eval $cmd
done
