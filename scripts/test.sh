#!/bin/bash

bla=false
if [ ! false ]; then
  echo ciao
fi

if [ $# -eq 2 ]; then
  run="param_run"
  [ ! -f "$1" ] && echo "ini file not found" && exit 1
  paramsini="$PWD/$1"
  [ ! -f "$2" ] && echo "pre file not found" && exit 1
  paramspre="$PWD/$2"
elif [ $# -eq 1 ]; then
  run="load_run"
  [ ! -d "$1" ] && echo "directory not found" && exit 1
  rundir="$1"
else
  echo "specify at least one parameter"
  exit 1
fi

