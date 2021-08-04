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
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['lines.linewidth'] = 1.5

plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.linewidth'] = 1.0
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14

plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1

plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.title_fontsize'] = 14

plt.rcParams['figure.titlesize'] = 14
plt.rcParams['figure.figsize'] = 4, 3
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.bbox'] = "tight"

sys.stderr.write("""Imported modules for python bioinformatics analysis
Author: Jinyang Zhang
Email: zhangjinyang@biols.ac.cn
Version: 1.0\n""")
sys.stderr.write(sys.version + '\n')
