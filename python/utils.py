import os
import sys
import re
import itertools
import time


class ProgressBar(object):
    """A simple progress bar with date stamp"""
    def __init__(self, width=50):
        """Init the ProgressBar object

        Paramters
        ---------
        width : integer, optional (default=50)
            The width of progress bar
        """
        self.last_x = -1
        self.width = width

    def update(self, x):
        """Update progress bar
        
        Paramters
        ---------
        x : integer, (0 <= x <= 100)
            the percentage of progress in [0, 100], if x equals input 
            from the last time, the progress bar will not be updated
        """

        assert 0 <= x <= 100
        if self.last_x == int(x):
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()
        if x == 100:
            sys.stderr.write('\n')


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
