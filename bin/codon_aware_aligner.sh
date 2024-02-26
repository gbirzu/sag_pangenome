#!/bin/bash -l
set -e
set -u

usage_message="Usage: codon_aware_aligner.sh -i INPUT_FILE -o OUTPUT_FILE [-t NUM_THREADS]"
num_threads=1

# Read input parameters
while getopts "i:o:t:" OPTION
do
  case $OPTION in
    i)
      input_file=${OPTARG}
      ;;
    o)
      output_file=${OPTARG}
      ;;
    t)
      num_threads=${OPTARG}
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      echo $usage_message 1>&2
      exit
      ;;
    ?)
      echo $usage_message 1>&2
      exit
      ;;
  esac
done

echo ${num_threads}

out_aa=$(echo ${output_file} | sed 's/\.fna/.faa/g')
python3 pg_process_alignment_files.py -i ${input_file} -o ${out_aa} --translate_sequences

out_aa_aln=$(echo ${out_aa} | sed 's/\.faa/_aln.faa/g')
mafft --thread ${num_threads} --quiet --reorder --auto ${out_aa} > ${out_aa_aln}

python3 pg_process_alignment_files.py -i ${out_aa_aln} -o ${output_file} --back_translate

if [ $? -eq 0 ]
then
    # Clean up
    rm -f ${out_aa}
    rm -f ${out_aa_aln}
fi

