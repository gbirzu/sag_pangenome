import argparse
import os
import sys
import codon_aware_align as codon_align
import seq_processing_utils as seq_utils
import pangenome_utils as pg_utils
import numpy as np


def translate_sequences(f_in, f_out=None, in_ext='fna', out_ext='faa'):
    input_seqs = codon_align.read_input_seqs(f_in)
    if f_out is not None:
        codon_align.write_aa_file(input_seqs, f_out)
    else:
        f_out = f_in.replace(f'.{in_ext}', f'.{out_ext}')
        codon_align.write_aa_file(input_seqs, f_out)

def back_translate_alignment(f_in, f_nucl, f_out=None):
    if os.path.exists(f_nucl):
        aa_aln = seq_utils.read_alignment(f_in)
        nucl_seqs = codon_align.read_input_seqs(f_nucl)
        nucl_aln = seq_utils.back_translate(aa_aln, nucl_seqs)

        # Write sequences to file
        if f_out is None:
            f_out = f_in.replace('.faa', '.fna')
        seq_utils.write_alignment(nucl_aln, f_out)

    else:
        print(f'ERROR: Sequences file {f_nucl} not found!')
        sys.exit(-1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--direction', default='forward', help='Direction of translation, with "nucl" -> "prot" being "forward" ["forward", "backward"].')
    parser.add_argument('-e', '--ext', default='fna', help='FASTA file extension.')
    parser.add_argument('-i', '--input_file', default=None, help='FASTA file with cluster sequences.')
    parser.add_argument('-o', '--output_file', default=None, help='FASTA file with output sequences.')
    parser.add_argument('-n', '--nucleotide_file', default=None, help='Needed for "backward" translation only.')
    args = parser.parse_args()

    if args.direction == 'forward':
        translate_sequences(args.input_file, f_out=args.output_file, in_ext=args.ext)
    elif args.direction == 'backward':
        back_translate_alignment(args.input_file, args.nucleotide_file, f_out=args.output_file)
    else:
        print(f'ERROR: Invalid translation direction {args.direction}. Choose between "forward" and "backward".')
        sys.exit(-1)




