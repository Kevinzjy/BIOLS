import os
import sys
import re
import itertools
import time


def ranking(ranks, names, order=1):
    """Ranking a list using another list as key

    Paramters
    ---------
    ranks : list,
        list of keys for ranking

    names : list,
        list of labels corresponding to each key in ranks

    order : 1 or -1
        if order equals -1, the order is reversed

    """
    import numpy as np
    from sklearn.preprocessing import MinMaxScaler
    minmax = MinMaxScaler()
    ranks = minmax.fit_transform(order * np.array([ranks]).T).T[0]
    ranks = map(lambda x: round(x, 2), ranks)
    return dict(zip(names, ranks))


def check_file(file_name):
    """Check if a file exists

    Parameters
    ----------
    file_name : str
        Name of file

    Returns
    -------
    str
        Absolute path of file
    """
    if os.path.exists(file_name) and os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        sys.exit('File: {}, not found'.format(file_name))


def check_dir(dir_name):
    """Check if directory exists

    Parameters
    ----------
    dir_name : str
        Name of directory, create file if dir not exists

    Returns
    -------
    str
        Absolute path of directory
    """
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            # Output directory already exists
            pass
        else:
            sys.exit('Directory: {}, clashed with existed files'.format(dir_name))
    else:
        os.mkdir(dir_name)
    return os.path.abspath(dir_name)


def to_str(bytes_or_str):
    """
    Return Instance of str
    """
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    """
    Return Instance of bytes
    """
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def empty_iter(iterable):
    """
    Return None if the iterable object is empty
    """
    import itertools
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return itertools.chain([first], iterable)


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def flatten(x):
    """
    Flatten list of lists
    """
    import itertools

    flatted_list = list(itertools.chain(*x))
    return flatted_list


def pairwise(iterable):
    """"
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    from itertools import tee
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def tree():
    """
    Tree structure
    """
    from collections import defaultdict
    return defaultdict(tree)


def generate_random_key(length):
    """
    Generate key of specific length
    """
    import random
    import string
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(length))


def hash_str(f):
    """
    Return sha256 of input string
    """
    import hashlib
    return hashlib.sha256(str(f).encode()).hexdigest()


def sorted_iters(iters, key, reverse=False):
    """
    Return iters with minimum values of given keys

    Parameters
    ----------
    iters : iterable object
        Iterable object
    key : key
        key to sort
    reverse : boolean
        if True, return maximum values

    Returns
    -------
    str
        Absolute path of directory
    """
    from operator import itemgetter
    x = sorted(iters, key=itemgetter(key), reverse=reverse)
    return [i for i in x if i[key] == x[0][key]]


def download(url, outfile):
    import requests
    from contextlib import closing
    from py3biotools.logger import ProgressBar

    with closing(requests.get(url, stream=True)) as response, open(outfile, 'wb') as out:
        chunk_size = 1024
        try:
            content_size = int(response.headers['content-length'])
            sys.stderr.write('Downloading file to {}, total size: {}\n'.format(outfile, content_size))
            prog = ProgressBar()
            cnt = 0
            for data in response.iter_content(chunk_size=chunk_size):
                out.write(data)
                cnt += len(data)
                prog.update(100 * cnt / content_size)
        except KeyError:
            out.write(response.content)
    return 1


def merge_intervals(blocks):
    from operator import itemgetter
    tmp = sorted(blocks, key=itemgetter(0, 1))
    merged = []
    last_st, last_en = tmp[0][0], tmp[0][1]
    for x in tmp[1:]:
        st, en = x[0], x[1]
        if st <= last_en + 1:
            last_en = max(en, last_en)
            last_st = min(st, last_st)
        else:
            merged.append([last_st, last_en])
            last_st, last_en = st, en
    merged.append([last_st, last_en])
    return merged


def cluster_peaks(x, gap=5):
    if len(x) == 0:
        return []
    x_sorted = sorted(x)
    clustered = [[x_sorted[0], x_sorted[0]], ]
    for i in x_sorted[1:]:
        if i-clustered[-1][1] > gap:
            clustered.append([i,i])
        else:
            clustered[-1][1] = i
    return clustered