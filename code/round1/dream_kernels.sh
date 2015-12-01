#!/usr/bin/env bash
# ------------------------------------------------------------------
shopt -s extglob

abspath_script="$(readlink -f -e "$0")"
script_absdir="$(dirname "$abspath_script")"
script_name="$(basename "$0" .sh)"

#if [ $# -eq 0 ] 
#    then
#        cat "$script_absdir/${script_name}_help.txt"
#        exit 1
#fi

TEMP=$(getopt -o hk: -l help,whatever: -n "$script_name.sh" -- "$@")

if [ $? -ne 0 ] 
then
  echo "Terminating..." >&2
  exit -1
fi

eval set -- "$TEMP"


# Defaults
whatever=2

while true
do
  case "$1" in
    -h|--help)    
      cat "$script_absdir/${script_name}_help.txt"
      exit
      ;;  
    -k|--whatever)  
      whatever="$2"
      shift 2
      ;;  
   --) 
      shift
      break
      ;;  
    *)  
      echo "$script_name.sh:Internal error!"
      exit 1
      ;;  
  esac
done


# Run all R scripts to obtain kernels in parallel
kernel_genex.R &
kernel_cnv.R &
kernel_mutation.R &
#kernel_methylation.R

wait %1 %2 %3 || exit $?
