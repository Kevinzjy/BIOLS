from collections import defaultdict


def collect_kmers(seq, k=8, use_hpc=False, is_circular=False):
    """
    split sequence into kmers
    :param seq: sequence
    :param k: k-mer size
    :param use_hpc: True if use homopolymer-compressed kmers
    :return: dict of kmers and occurence position
    """
    kmers = {}
    kmer_occ = defaultdict(list)

    split_func = split_hpc_kmers if use_hpc else split_kmers
    if is_circular:
        tmp_kmers, tmp_kmer_occ = collect_kmers(seq * 2, k, use_hpc, False)
        for i in tmp_kmers:
            if i < len(seq):
                continue
            kmers[i - len(seq)] = tmp_kmers[i]
            kmer_occ[tmp_kmers[i]].append(i - len(seq))
    else:
        for i, kmer in split_func(seq, k):
            kmers[i] = kmer
            kmer_occ[kmer].append(i)

    return kmers, kmer_occ


def split_kmers(seq, k=8):
    for x in range(len(seq)):
        if x >= k - 1:
            yield x, seq[x-k+1:x+1]


def split_hpc_kmers(seq, k=8):
    hpc = [seq[0], ]
    for x, (i, j) in enumerate(zip(seq[:-1], seq[1:])):
        if i != j:
            hpc.append(j)
            if len(hpc) >= k:
                yield x + 1, ''.join(hpc[-k:])


def compress_seq(seq):
    hpc = [seq[0], ]
    for x, (i, j) in enumerate(zip(seq[:-1], seq[1:])):
        if i != j:
            hpc.append(j)
    return ''.join(hpc)