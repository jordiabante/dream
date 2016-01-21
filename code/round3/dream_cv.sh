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

TEMP=$(getopt -o ht:o:a:r: -l help,threads:,outdir:,algorithm:,round: -n "$script_name.sh" -- "$@")

if [ $? -ne 0 ] 
then
  echo "Terminating..." >&2
  exit -1
fi

eval set -- "$TEMP"


# Defaults
threads=2
round=3
algorithm="ridge"

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
    -o|--outdir)  
      outdir="$2"
      shift 2
      ;;  
    -a|--algorithm)  
      algorithm="$2"
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

# Output
outdir="${outdir}/${algorithm}"

mkdir -p "$outdir"

# Run cross-validation in parallel
export round
export algorithm

cat "${suffixes}" | xargs -i -n 1 --max-proc "${threads}" bash -c \
	'./R/${algorithm}_cross_validation_parameter.R \
	../../data/round${round}/kernel_train_test_product/kernel_train_{} \
	../../data/round${round}/kernel_train_test_product/kernel_test_{} \
	1> ${outdir}/{}_${algorithm}_cv.txt'
