import argparse
import pangenome_utils as pg_utils

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='File with merged BLASTP results.')
    parser.add_argument('-o', '--output_file', help='Output file (for single query file only).')
    args = parser.parse_args()

    pg_utils.write_protein_graph(args.input_file, args.output_file)


