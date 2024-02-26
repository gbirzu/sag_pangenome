import argparse
import os
import sys
import glob
import time
import codon_aware_align as codon_align
import seq_processing_utils as seq_utils
import pangenome_utils as pg_utils
import numpy as np
import alignment_tools as align_utils
from pangenome_utils import PangenomeMap
from Bio import AlignIO


def translate_sequences(f_in, f_out=None, in_ext='fna', out_ext='faa'):
    input_seqs = codon_align.read_input_seqs(f_in)
    if f_out is not None:
        codon_align.write_aa_file(input_seqs, f_out)
    else:
        f_out = f_in.replace(f'.{in_ext}', f'.{out_ext}')
        codon_align.write_aa_file(input_seqs, f_out)

def back_translate_alignment(f_in, nucl_seqs_dir, f_out=None, verbose=False):
    og_id = f_in.split('/')[-1].replace('_aln.faa', '')
    f_nucl = f'{nucl_seqs_dir}{og_id}.fna'
    if os.path.exists(f_nucl):
        if verbose:
            print(f'Back-translating {f_in} alignment...')

        aa_aln = AlignIO.read(f_in, 'fasta')
        nucl_seqs = codon_align.read_input_seqs(f_nucl)
        nucl_aln = seq_utils.back_translate(aa_aln, nucl_seqs)

        # Write sequences to file
        if f_out is None:
            f_out = f_in.replace('.faa', '.fna')
        AlignIO.write(nucl_aln, f_out, 'fasta')

    else:
        print(f'ERROR: No sequences file found for {og_id} in {nucl_seqs_dir}!')
        sys.exit(-1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--input_dir', default=None, help='Input dir with FASTA file for each protein cluster.')
    parser.add_argument('-O', '--output_dir', help='Output directory.')
    parser.add_argument('-e', '--ext', default='fna', help='FASTA file extension.')
    parser.add_argument('-i', '--input_file', default=None, help='FASTA file with cluster sequences.')
    parser.add_argument('-f', '--file_list', default=None, help='Text file with file paths to process.')
    parser.add_argument('-o', '--output_file', default=None, help='FASTA file with output seqeunces.')
    parser.add_argument('-s', '--seq_type', default='nucl', help='Sequence types ["prot" | "nucl"].')
    parser.add_argument('-v', '--verbose', action='store_true', help='Run in verbose mode.')
    parser.add_argument('--translate_sequences', action='store_true', help='Translate input sequences and save results to file(s).')
    parser.add_argument('--back_translate', action='store_true', help='Back translate AA alignments to nucleotides.')
    parser.add_argument('--split_mixed_orthogroup_alignments', action='store_true', help='Split OG alignments for mixed clusters with fine-scale subclusters.')
    args = parser.parse_args()

    if args.translate_sequences == True:
        if args.input_file:
            translate_sequences(args.input_file, f_out=args.output_file, in_ext=args.ext)

        elif args.input_dir:
            if args.verbose:
                print(f'Translating FASTA files in {args.input_dir}')

            input_files = glob.glob(f'{args.input_dir}*.{args.ext}')
            for f_in in input_files:
                translate_sequences(f_in, in_ext=args.ext)

        elif args.file_list:
            with open(args.file_list, 'r') as in_handle:
                input_files = [line.strip() for line in in_handle.readlines()]
                for f_in in input_files:
                    translate_sequences(f_in, in_ext=args.ext)

        else:
            print('No input given. Use -i, -f, or -I (see help).')
            sys.exit(-1)

    if args.back_translate == True:
        if (args.file_list is not None) and (args.input_dir is not None):
                with open(args.file_list, 'r') as in_handle:
                    input_files = [line.strip() for line in in_handle.readlines()]
                    for f_in in input_files:
                        back_translate_alignment(f_in, args.input_dir)

        elif (args.input_file is not None) and (args.input_dir is not None):
            back_translate_alignment(args.input_file, args.input_dir, f_out=args.output_file)

        else:
            print('ERROR: Need both list of input files and input directory for original sequences to back translate alignments.')
            sys.exit(-1)


