import argparse
import subprocess
import seq_processing_utils as seq_utils
from Bio import SeqIO
from Bio import AlignIO


def read_input_seqs(in_fasta, return_seq_order=False):
    seq_records = SeqIO.parse(in_fasta, 'fasta')
    nucl_seqs = {}
    for rec in seq_records:
        nucl_seqs[rec.id] = rec

    if return_seq_order == True:
        rec_ids = []
        with open(in_fasta, 'r') as in_handle:
            for line in in_handle.readlines():
                line = line.rstrip()
                if line.startswith('>'):
                    rec_id = line.split(' ')[0].strip('>')
                    rec_ids.append(rec_id)
        return nucl_seqs, rec_ids
    else:
        return nucl_seqs

def write_aa_file(nucl_seqs, output_file, rec_ids=None):
    if rec_ids is None:
        records = [rec.translate(table='Bacterial', id=rec.id, description=rec.description) for rec in nucl_seqs.values()]
    else:
        nucl_seq_records = [nucl_seqs[rec_id] for rec_id in rec_ids]
        records = [rec.translate(table='Bacterial', id=rec.id, description=rec.description) for rec in nucl_seq_records]
    SeqIO.write(records, output_file, 'fasta')

def run_tests(args):
    print('Testing codon-aware aligner...')

    test_dir = '../results/tests/locus_alignments/'
    f_seqs = f'{test_dir}rplE.fna'
    f_out = f'{test_dir}rplE_codon_aln_test.fna'
    f_aa_seqs = f_seqs.replace(f'.fna', '.faa')

    nucl_seqs = read_input_seqs(f_seqs)
    write_aa_file(nucl_seqs, f_aa_seqs)

    # Perform alignment
    f_aa_aln = f_out.replace(f'.fna', '.faa')
    if args.aligner == 'mafft':
        with open(f_aa_aln, 'w') as out_handle:
            aln_out = subprocess.run([args.aligner_bin, '--quiet', '--reorder', '--thread', f'{args.num_threads}', '--auto', f_aa_seqs], stdout=out_handle, stderr=subprocess.STDOUT)
        aa_aln = seq_utils.read_mafft_alignment(f_aa_aln)
    else:
        aln_out = subprocess.run([args.aligner_bin, '-in', f_aa_seqs, '-out', f_aa_aln], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        aa_aln = AlignIO.read(f_aa_aln, 'fasta')
    if args.verbose:
        print(aln_out)

    nucl_aln = seq_utils.back_translate(aa_aln, nucl_seqs)
    AlignIO.write(nucl_aln, f_out, 'fasta')

    if args.cleanup == True:
        # Clean up
        subprocess.call(' '.join(['rm', '-f', f_aa_seqs]), shell=True)
        subprocess.call(' '.join(['rm', '-f', f_aa_aln]), shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_fasta', type=str, help='Path to input FASTA file (must have .fasta extension).')
    parser.add_argument('-a', '--aligner', type=str, default='muscle', help='Aligner to use (MUSCLE by default).')
    parser.add_argument('-b', '--aligner_bin', type=str, default='muscle', help='Path to aligner binary program.')
    parser.add_argument('-e', '--extension', default='fasta', help='FASTA extension.')
    parser.add_argument('-o', '--output_file', type=str, default=None, help='Output file for nucleotide alignments.')
    parser.add_argument('-c', '--cleanup', action='store_true', help='Remove intermediate files at the end.')
    parser.add_argument('-n', '--num_threads', default=1, type=int, help='Number of threads.')
    parser.add_argument('-t', '--test', action='store_true', help='Run tests.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print MUSCLE output.')
    args = parser.parse_args()

    if args.test:
        run_tests(args)
        exit()

    if args.output_file is None:
        f_out = args.in_fasta.replace(f'.{args.extension}', f'_aln.{args.extension}')
    else:
        f_out = args.output_file

    #nucl_seqs, rec_ids = read_input_seqs(args.in_fasta, return_seq_order=True)
    nucl_seqs = read_input_seqs(args.in_fasta)
    f_aa_seqs = args.in_fasta.replace(f'.{args.extension}', '.faa')
    #write_aa_file(nucl_seqs, f_aa_seqs, rec_ids=rec_ids)
    write_aa_file(nucl_seqs, f_aa_seqs)

    # Perform alignment
    f_aa_aln = f_out.replace(f'.{args.extension}', '.faa')
    if args.aligner == 'mafft':
        with open(f_aa_aln, 'w') as out_handle:
            #aln_out = subprocess.run([args.aligner_bin, '--quiet', '--thread', f'{args.num_threads}', '--auto', f_aa_seqs], stdout=out_handle, stderr=subprocess.STDOUT)
            aln_out = subprocess.run([args.aligner_bin, '--quiet', '--reorder', '--thread', f'{args.num_threads}', '--auto', f_aa_seqs], stdout=out_handle, stderr=subprocess.STDOUT)
        aa_aln = seq_utils.read_mafft_alignment(f_aa_aln)
    else:
        aln_out = subprocess.run([args.aligner_bin, '-in', f_aa_seqs, '-out', f_aa_aln], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        aa_aln = AlignIO.read(f_aa_aln, 'fasta')
    if args.verbose:
        print(aln_out)

    nucl_aln = seq_utils.back_translate(aa_aln, nucl_seqs)
    AlignIO.write(nucl_aln, f_out, 'fasta')

    if args.cleanup == True:
        # Clean up
        subprocess.call(' '.join(['rm', '-f', f_aa_seqs]), shell=True)
        subprocess.call(' '.join(['rm', '-f', f_aa_aln]), shell=True)
