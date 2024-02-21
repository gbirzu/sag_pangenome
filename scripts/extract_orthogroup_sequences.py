import argparse
import numpy as np
import pangenome_utils as pg_utils
import seq_processing_utils as seq_utils


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-A', '--annotations_dir', help='Directory with SAG annotations in .gff format.')
    parser.add_argument('-O', '--output_dir', help='Output directory.')
    parser.add_argument('-g', '--orthogroup_table', help='Input orthogroup table.')
    parser.add_argument('-s', '--seq_type', default='nucl', help='Output sequence type ["prot" | "nucl"].')
    args = parser.parse_args()

    # Read FASTA format
    if args.seq_type == 'nucl':
        ext = 'fna'
        seq_type = args.seq_type
    elif args.seq_type == 'prot':
        ext = 'faa'
        seq_type = args.seq_type
    else:
        print(f'Unknown sequence type: {args.seq_type}! Defaulting to "nucl".')
        seq_type = 'nucl'
        ext = 'fna'

    pangenome_map = pg_utils.PangenomeMap(args.annotations_dir, f_orthogroup_table=args.orthogroup_table)
    og_table = pangenome_map.og_table
    sag_ids = pangenome_map.get_sag_ids()

    for cluster_id in og_table.index:
        cluster_genes = np.concatenate([s.split(';') for s in og_table.loc[cluster_id, sag_ids].dropna()])

        seq_records = []
        for seq_id in cluster_genes:
            record = pangenome_map.extract_gene_record(seq_id, type=seq_type)
            seq_records.append(record)

        # Verify all sequences were added
        assert len(cluster_genes) == og_table.loc[cluster_id, 'num_seqs']
        assert len(seq_records) == len(cluster_genes)

        # Write sequences
        seq_utils.write_seqs(seq_records, f'{args.output_dir}{cluster_id}.{ext}', 'fasta')

