# Modify a SONG run so that sources & transfer functions are 
# not computed the next time it is used.
#
# Useful only when we want to reuse a run that contains bispectra to
# compute Fisher matrices, without having to recompute the source &
# transfer functions.

#! /bin/bash

folders=$@

for folder in $folders
do
  # Create sources and transfer folders
  mkdir -p $folder/sources
  mkdir -p $folder/transfers
  # Update parameter file
  sed -i '.bak' 's/^store_sources.*/store_sources = yes/' $folder/run_params.ini
  sed -i '.bak' 's/^store_transfers.*/store_transfers = yes/' $folder/run_params.ini
  # Remove date, if any
  mv $folder ${folder/_2013*/}    
done
