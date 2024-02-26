#!/bin/bash

ml biology
ml bowtie2

# Default parameters
split_files=0 # assume fastq is paired-end and interleaved
usage_message="Usage: ./blast_reads.sh [-i | -u | -1 -2] PATH_TO_READS -r PATH_TO_REFERENCE -o OUTPUT_FILE"

# Read input parameters
while getopts "i:r:u:1:2:o:" OPTION
do
  case $OPTION in
    i)
      fastq_file=${OPTARG}
      ;;
    r)
      ref_file=${OPTARG}
      ;;
    u)
      unpaired_file=${OPTARG}
      split_files=1
      ;;
    1)
      forward_file=${OPTARG}
      split_files=2
      ;;
    2)
      reverse_file=${OPTARG}
      ;;
    o)
      output_file=${OPTARG}
      ;;
    \?)
      echo $usage_message 1>&2
      exit
      ;;
  esac
done

if [ $OPTIND -eq 1 ]; then
  echo $usage_message 1>&2
  exit
elif [ -z "$output_file" ]; then
  echo "Error: output file required!" 1>&2
  echo $usage_message 1>&2
  exit
fi

blastn -db ../data/reference_genomes/Synechococcus_OS -task blastn -out blast_default.out -query ../results/tests/UncmicMusRedA1C3_FD.fasta
