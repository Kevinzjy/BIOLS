#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
import sys
import click
import logging
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parents[1]))
from py3biotools import *

formatter = logging.Formatter("%(asctime)-15s [%(levelname)-5s] %(message)s", "[%a %Y-%m-%d %H:%M:%S]")
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setFormatter(formatter)

LOGGER = logging.getLogger(__file__)
LOGGER.setLevel(logging.DEBUG)
LOGGER.addHandler(console_handler)


def load_fasta(infile):
    from py3biotools.bioinfo import Fasta
    fasta = Fasta(str(infile))
    return fasta


def load_annotation(infile):
    from collections import defaultdict
    from py3biotools.bioinfo import GTFParser
    gene_index = defaultdict(dict)
    tscp_index = defaultdict(list)
    exon_index = defaultdict(list)

    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            # only include gene and exon feature for now
            if content[2] not in ['gene', 'transcript', 'exon']:
                continue
            parser = GTFParser(content)

            # Extract splice site
            if content[2] == 'gene':
                gene_index[parser.attr['gene_id']] = parser
            if content[2] == 'transcript':
                tscp_index[parser.attr['gene_id']].append(parser)
            if content[2] == 'exon':
                exon_index[parser.attr['transcript_id']].append(parser)
    return gene_index, tscp_index, exon_index


def ensure_intron(x, y):
    if x.start < y.start:
        return (x.end, y.start)
    else:
        return (y.end, x.start)


GENE_INDEX = None
TSCP_INDEX = None
EXON_INDEX = None


def initializer(gene_index, tscp_index, exon_index):
    global GENE_INDEX, TSCP_INDEX, EXON_INDEX
    GENE_INDEX = gene_index
    TSCP_INDEX = tscp_index
    EXON_INDEX = exon_index


def get_introns(chunk):
    import pyranges as pr
    from operator import itemgetter
    gene_intron = {}

    for gene_id in chunk:
        introns = []
        for tscp in TSCP_INDEX[gene_id]:
            for i in np.arange(1, len(EXON_INDEX[tscp.attr['transcript_id']])):
                introns.append(ensure_intron(EXON_INDEX[tscp.attr['transcript_id']][i-1], EXON_INDEX[tscp.attr['transcript_id']][i]))
        intron_df = pr.PyRanges(pd.DataFrame({
            'Chromosome': GENE_INDEX[gene_id].contig,
            'Start': list(map(itemgetter(0), introns)),
            'End': list(map(itemgetter(1), introns)),
            'Strand': GENE_INDEX[gene_id].strand,
        }))
        gene_intron[gene_id] = intron_df
    return gene_intron


def format_bed(x):
    tmp_line = [x['Chromosome'], x['Start']+1, x['End']-1, '{}:{}-{}'.format(x['Chromosome'], x['Start']+1, x['End']-1), '.', x['Strand']]
    return '\t'.join([str(i) for i in tmp_line])


def get_intron_sequence(fasta, x):
    from py3biotools.bioinfo import revcomp
    seq = fasta.seq(x['Chromosome'], x['Start']+1-1, x['End']-1)
    if x['Strand'] == '-':
        seq = revcomp(seq)
    return seq


@click.command()
@click.option('--genome', '-g', default=None, help='Genome FASTA')
@click.option('--anno', '-a', default=None, help='Annotation GTF')
@click.option('--bed', '-b', default=None, help='Output bed file')
@click.option('--fasta', '-f', default=None, help='Output FASTA')
@click.option('--cpus', '-t', default=6, help='Threads')
def main(genome, anno, bed, cpus, fasta):
    from py3biotools.logger import ProgressBar

    LOGGER.info('Running {} to extract all intron position and sequence'.format(Path(__file__)))
    if genome is None or anno is None:
        LOGGER.error('Wrong argument specified, please run with --help to see further instructions.')
    if bed is None and fasta is None:
        LOGGER.error('No output specified, please run with --help to see further instructions.')

    # Load genome
    LOGGER.info('Loading genome fasta: {}'.format(Path(genome).absolute()))
    faidx = load_fasta(Path(genome).absolute())

    # Load reference
    LOGGER.info('Loading genome annotation: {}'.format(Path(anno).absolute()))
    gene_index, tscp_index, exon_index = load_annotation(Path(anno).absolute())

    # Get introns with multiprocessing
    pool = Pool(int(cpus), initializer, (gene_index, tscp_index, exon_index, ))
    jobs = []
    for tmp in grouper(list(gene_index), 1000):
        chunk = [i for i in tmp if i is not None]
        jobs.append(pool.apply_async(get_introns, (chunk, )))
    pool.close()

    gene_intron = {}
    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    for job in jobs:
        tmp_intron = job.get()
        gene_intron.update(tmp_intron)
        cnt += 1
        prog.update(100*cnt/len(jobs))
    pool.join()
    prog.update(100)

    # Merge introns of protein coding genes
    LOGGER.info('Merge intron of all protein coding genes')
    all_annotated_introns = []
    for gene_id in gene_intron:
        if gene_intron[gene_id].df.shape[0] == 0:
            continue
        if 'gene_type' in gene_index[gene_id].attr and gene_index[gene_id].attr['gene_type'] != 'protein_coding':
            continue
        elif 'gene_biotype' in gene_index[gene_id].attr and gene_index[gene_id].attr['gene_biotype'] != 'protein_coding':
            continue
        else:
            pass
        all_annotated_introns.append(gene_intron[gene_id].df)
    all_annotated_introns = pd.concat(all_annotated_introns)
    all_annotated_introns = all_annotated_introns.drop_duplicates()
    all_annotated_introns.index = all_annotated_introns.apply(lambda x: '{}:{}-{}'.format(x['Chromosome'], x['Start']+1, x['End']-1), axis=1)

    # Output
    if bed is not None:
        LOGGER.info('Output bed: {}'.format(Path(bed).absolute()))
    if bed is not None:        
        with open(bed, 'w') as out:
            for i in all_annotated_introns.apply(format_bed, axis=1).values:
                out.write(i + '\n')

    if fasta is not None:
        LOGGER.info('Output FASTA: {}'.format(Path(fasta).absolute()))
        if fasta is not None:
            with open(fasta, 'w') as out:
                for idx, row in all_annotated_introns.iterrows():
                    intron_id = '{}:{}-{}'.format(row['Chromosome'], row['Start']+1, row['End']-1)
                    seq = get_intron_sequence(faidx, row)
                    out.write('>{}\n{}\n'.format(intron_id, seq))


if __name__ == '__main__':
    main()
