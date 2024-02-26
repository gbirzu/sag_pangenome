#!/bin/bash -l

usage_message="Usage: align_file_by_index.sh -I INPUT_DIR -O OUTPUT_DIR -i INDEX"

# Read input parameters
while getopts "I:O:i:" OPTION
do
  case $OPTION in
    I)
      input_dir=${OPTARG}
      ;;
    O)
      output_dir=${OPTARG}
      ;;
    i)
      index=${OPTARG}
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      echo $usage_message 1>&2
      exit
      ;;
    *)
      echo $usage_message 1>&2
      exit
      ;;
  esac
done

if [ $OPTIND -ne 7 ]
then
    echo "Error: Missing argument(s)."
    echo $usage_message 1>&2
    exit
fi

echo $input_dir, $output_dir, $index

