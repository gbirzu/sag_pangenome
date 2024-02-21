import argparse
import pandas as pd
import pickle
import glob
import utils 
import seq_processing_utils as seq_utils

def extract_seq_records(annotations_df, contig_records, type='nucl', format='genes'):
    seq_records = {}
    if format == 'genes':
        gene_idx = annotations_df.loc[annotations_df['type'] == 'CDS', :].index # include CDS only
        for gene_id in gene_idx: 
            contig_id = annotations_df.loc[gene_id, 'contig']
            contig_seq = contig_records[contig_id]
            x_start, x_end = annotations_df.loc[gene_id, ['start_coord', 'end_coord']].astype(int)
            gene_seq = contig_seq.seq[x_start - 1:x_end]
            if annotations_df.loc[gene_id, 'strand'] == '-':
                gene_seq = gene_seq.reverse_complement()
            if type == 'prot':
                gene_seq = gene_seq.translate(table=11)
            seq_records[gene_id] = seq_utils.make_SeqRecord(gene_seq, id=gene_id, description='')

    elif format == 'contigs':
        for contig_id in contig_records:
            seq_records[contig_id] = contig_records[contig_id]

    return seq_records


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--input_dir', default=None, help='Directory with input GFF3-FASTA files.')
    parser.add_argument('-O', '--output_dir', default='', help='Directory in which gene seq FASTA files are saved.')
    parser.add_argument('-f', '--output_format', default='genes', help='Parsing for output FASTA file ["genes", "contigs"].')
    parser.add_argument('-i', '--input_file', default=None, help='Path of input GFF3-FASTA file.')
    parser.add_argument('-o', '--output_file', default=None, help='Path of output gene seq FASTA file.')
    parser.add_argument('-s', '--seq_type', default='nucl', help='Sequence type ["nucl", "prot"].')
    parser.add_argument('-v', '--verbose', action='store_true', help='Run in verbose mode.')
    args = parser.parse_args()

    if args.input_file is not None:
        annotations_df, contig_records = utils.read_annotation_file(args.input_file, format='gff-fasta')
        seq_records = extract_seq_records(annotations_df, contig_records, type=args.seq_type, format=args.output_format)

        if args.output_file is not None:
            seq_utils.write_seqs_dict(seq_records, args.output_file)
        else:
            print('Warning: --output_file not set.')
            print(annotations_df)
        
    elif args.input_dir is not None:
        annotation_files = sorted(glob.glob(f'{args.input_dir}*.gff'))
        for f_in in annotation_files:
            annotations_df, contig_records = utils.read_annotation_file(f_in, format='gff-fasta')
            seq_records = extract_seq_records(annotations_df, contig_records, type=args.seq_type)

            sample_id = f_in.split('/')[-1].replace('.gff', '')
            f_out = f'{args.output_dir}{sample_id}_filtered_genes.fna'
            seq_utils.write_seqs_dict(seq_records, f_out)

        


