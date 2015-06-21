## Execute SONG sequentially on all folders give as arguments.
##
## The first argument is the number of threads to use for the parallel
## computation; all remaining arguments are interpreted as a list
## of run folders on which SONG will be run.
##
## The script will dump all messages from SONG in a log.txt file
## inside each folder.

#! /bin/sh

threads=$1
shift
folders=$@

make clean;
make -f makefile -j8 song > /dev/null;
export OMP_NUM_THREADS=$threads;

for folder in $folders; do
  cmd="unbuffer time ./song $folder 2>&1 | tee $folder/log.txt"
  echo "$cmd"
  eval "$cmd"
done
