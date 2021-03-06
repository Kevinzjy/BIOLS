import os
import sys
import re
import pickle
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
from multiprocessing import Pool
from .utils import *

# Plotting parameters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Arial"

sys.stderr.write("""Imported modules for python bioinformatics analysis
Author: Jinyang Zhang
Email: zhangjinyang@biols.ac.cn
Version: 1.0\n""")
sys.stderr.write(sys.version + '\n')
