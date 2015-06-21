# Print to file a triangular slice of the bispectrum, fixing
# the two smallest multipoles.
#
# usage:    ./bisp_triangular.sh <bispectrum_type> <l1> <l2> <run_directory>
#
# The bispectra will be read from the run directory, and the output
# files will be stored in the following file:
#
#   run_directory/triangular_<l1>_<l2>_<bispectrum_type>.txt
#
# The output file will contain the bispectrum configurations B_L1_L2_L3
# with L1=<l1>, L2=<l2> and L3 varying according to the triangular
# condition and the l-list computed in the run.
#
# For example,
# 
#   ./bisp_triangular.sh local_tte 1000 850 test_run
#
# will extract from the SONG run in the folder 'test_run' the bispectrum
# configurations with L1=1000, L2=850 for the TTE local primordial bispectrum,
# and print them to the file test_run/triangular_1000_850_local_tte.txt
# for all computed values of L3.
#
# If you give a list of folders in <run_directory>, the script will
# process each of them.

#! /bin/sh

# Parse arguments
type_bispectrum=$1
shift
l1=$1
shift
l2=$1
shift
folders=$@

# Compile SONG
make print_bispectra -j4 > /dev/null;

# Print bispectra to file
for folder in $folders; do
  cmd="time ./print_bispectra $type_bispectrum triangular $l1 $l2 $folder 2> $folder/triangular_${l1}_${l2}_${type_bispectrum}.txt"
  echo $cmd
  eval $cmd
done
