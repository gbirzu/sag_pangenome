from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def write_seqs_dict(seq_dict, output_file, format='fasta'):
    records = list(seq_dict.values())
    SeqIO.write(records, output_file, 'fasta')

def write_seqs(seq_records, output_file, format='fasta'):
    SeqIO.write(seq_records, output_file, 'fasta')

def make_Seq(seq_str, **kwargs):
    '''
    Wrapper for Biopython Seq constructor
    '''
    return Seq(seq_str, **kwargs)

def make_SeqRecord(seq, **kwargs):
    '''
    Wrapper for Biopython SeqRecord constructor
    '''
    return SeqRecord(seq, **kwargs)

if __name__ == '__main__':
    print('No tests to run...')

