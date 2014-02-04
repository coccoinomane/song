#! /bin/bash

# if [[ $# < 2 ]]; then
#   echo "ERROR: provide at least 2 arguments"
#   exit 1
# fi
# 
# name=$1
# shift
folders=$@

# if [ $name != ""];
#   then
#   name="fnl.txt"
# fi

make -j2;

for folder in $folders; do
  cmd="time ./song $folder 2> $folder/fnl.txt"
  echo $cmd
  eval $cmd
done