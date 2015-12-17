#!/usr/bin/env bash
# ------------------------------------------------------------------
shopt -s extglob

abspath_script="$(readlink -f -e "$0")"
script_absdir="$(dirname "$abspath_script")"
script_name="$(basename "$0" .sh)"

if [ $# -eq 0 ] 
    then
        cat "${script_absdir}/${script_name}_help.txt"
        exit 1
fi

TEMP=$(getopt -o ht:r: -l help,threads:,round: -n "$script_name.sh" -- "$@")

if [ $? -ne 0 ] 
then
  echo "Terminating..." >&2
  exit -1
fi

eval set -- "$TEMP"


# Defaults
threads=2
round=2

while true
do
  case "$1" in
    -h|--help)    
      cat "${script_absdir}/${script_name}_help.txt"
      exit
      ;;  
    -t|--threads)  
      threads="$2"
      shift 2
      ;;  
    -r|--round)  
      round="$2"
      shift 2
      ;;  
   --) 
      shift
      break
      ;;  
    *)  
      echo "${script_name}.sh:Internal error!"
      exit 1
      ;;  
  esac
done

# Input
suffixes="$1"

# Run cross-validation in parallel
export round

cat "${suffixes}" | xargs -i -n 1 --max-proc "${threads}" bash -c \
	'./R/whatever.R ../../data/round${round}/kernels_train_test/kernel_train_{} \
	../../data/round${round}/kernels_train_test/kernel_test_{}' \
	1>> ../../data/round2/cross_validation.txt  2>>"${script_name}.log"
#	'echo ../../data/round${round}/kernels_train_test/kernel_train_{}' \
#	1>> ../../data/round2/cross_validation.txt  2>>"${script_name}.log"
    

