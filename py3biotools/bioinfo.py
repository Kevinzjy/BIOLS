import os
import sys
import pandas as pd
try:
    from commands import  getstatusoutput
except ImportError:
    from subprocess import getstatusoutput


class StupidError(Exception): 
    # Constructor or Initializer 
    def __init__(self, value): 
        self.value = value 
  
    # __str__ is to print() the value 
    def __str__(self): 
        return(repr(self.value)) 


class GTFParser(object):
    """
    Class for parsing annotation gtf
    """

    def __init__(self, content):
        self.contig = content[0]
        self.source = content[1]
        self.type = content[2]
        self.start, self.end = int(content[3]), int(content[4])
        self.strand = content[6]
        self.attr_string = content[8]

    @property
    def attr(self):
        """
        Parsing attribute column in gtf file
        """
        import re
        field = {}
        for attr_values in [re.split(r'\s+', i.strip()) for i in self.attr_string.split(';')[:-1]]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field

    def __repr__(self):
        return '{} {}:{}-{}:{}'.format(self.type, self.contig, self.start, self.end, self.strand)

    def seq(self, fasta):
        tmp_seq = fasta.seq(self.contig, self.start-1, self.end)
        if self.strand == '+':
            return tmp_seq
        else:
            return revcomp(tmp_seq)


class GFFParser(GTFParser):
    @property
    def attr(self):
        """
        Parsing attribute column in gtf file
        """
        import re
        field = {}
        for attr_values in [re.split(r'=', i.strip()) for i in self.attr_string.split(';')]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field


class Faidx(object):
    """
    API for unify pysam::Fastafile and mappy2::Aligner
    """
    def __init__(self, infile):
        """Input reference genome fasta"""
        import pysam
        if not os.path.exists(infile + '.fai'):
            sys.stderr.write('Index of reference genome not found, generating ...')
        try:
            faidx = pysam.FastaFile(infile)
        except ValueError:
            raise StupidError('Cannot generate index of reference genome, index file is missing')
        except IOError:
            raise StupidError('Cannot generate index of reference genome, file could not be opened')

        self.faidx = faidx
        self.contig_len = {contig: faidx.get_reference_length(contig) for contig in faidx.references}

    def seq(self, contig, start, end):
        """
        Return sequence of given coordinate
        """
        return self.faidx.fetch(contig, start, end)

    def close(self):
        self.faidx.close()


class Fasta(object):
    """
    load fasta file into memory
    """
    def __init__(self, infile):
        faidx = Faidx(infile)
        self.contig_len = faidx.contig_len
        self.genome = {ctg: faidx.seq(ctg, 0, size) for ctg, size in self.contig_len.items()}
        faidx.close()

    def seq(self, contig, start, end):
        if contig not in self.genome:
            return None
        return self.genome[contig][start:end]


def load_fasta(fname, is_gz=False):
    import gzip
    from .utils import to_str
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        if to_str(line).startswith('>'):
            if seq_id is not None:
                sequences[seq_id] = seq
            seq_id = to_str(line).rstrip().lstrip('>')
            seq = ''
        else:
            seq += to_str(line).rstrip()
    sequences[seq_id] = seq
    f.close()
    return sequences


def load_fastq(fname, is_gz=False):
    import gzip
    from py3biotools.utils import to_str
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        seq_id = to_str(line).rstrip().lstrip('@').split(' ')[0]
        seq = to_str(f.readline()).rstrip()
        sep = to_str(f.readline()).rstrip()
        qual = to_str(f.readline()).rstrip()
        sequences[seq_id] = (seq, qual)
    f.close()
    return sequences


def index_annotation(gtf):
    """
    Generate binned index for element in gtf
    """
    from .utils import tree

    gtf_index = defaultdict(dict)
    intron_index = defaultdict(dict)
    splice_site_index = tree()

    last_exon = None
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            # only include gene and exon feature for now
            if content[2] not in ['gene', 'exon']:
                continue

            parser = GTFParser(content)

            # Extract splice site
            if content[2] == 'exon':
                splice_site_index[parser.contig][parser.start][parser.strand]['start'] = 1
                splice_site_index[parser.contig][parser.end][parser.strand]['end'] = 1

                # Load intron
                if last_exon is not None and last_exon.attr['transcript_id'] == parser.attr['transcript_id']:
                    intron_start = last_exon.end if last_exon.strand == '+' else last_exon.start
                    intron_end = parser.start if parser.strand == '+' else parser.end
                    intron_strand = parser.strand

                    intron_start, intron_end = min(intron_start, intron_end), max(intron_start, intron_end)
                    start_div, end_div = intron_start // 500, intron_end // 500
                    for i in range(start_div, end_div + 1):
                        intron_index[parser.contig].setdefault(i, []).append((intron_start, intron_end, intron_strand))

                last_exon = parser

            # Binned index
            start_div, end_div = parser.start // 500, parser.end // 500
            for i in range(start_div, end_div + 1):
                gtf_index[parser.contig].setdefault(i, []).append(parser)

    return gtf_index, intron_index, splice_site_index


def revcomp(seq):
    """
    Convert sequence to reverse complementary
    """
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def get_bsj(seq, bsj):
    """Return transformed sequence of given BSJ"""
    return seq[bsj:] + seq[:bsj]
    

def get_n50(sequence_lengths):
    """
    Get n50 of sequence lengths
    """
    sequence_lengths = sorted(sequence_lengths, reverse=True)
    total_bases = sum(sequence_lengths)
    target_bases = total_bases * 0.5
    bases_so_far = 0
    for sequence_length in sequence_lengths:
        bases_so_far += sequence_length
        if bases_so_far >= target_bases:
            return sequence_length
    return 0
