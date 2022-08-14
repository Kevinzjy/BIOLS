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

plt.rcParams['lines.linewidth'] = .75

plt.rcParams['axes.linewidth'] = .75
plt.rcParams['grid.linewidth'] = .75
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['axes.titlesize'] = 7

plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['xtick.major.width'] = .75
plt.rcParams['ytick.major.width'] = .75
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['xtick.minor.width'] = .75
plt.rcParams['ytick.minor.width'] = .75

plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['legend.title_fontsize'] = 7

plt.rcParams['figure.titlesize'] = 7
plt.rcParams['figure.figsize'] = 1.25, 1.25
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.bbox'] = "tight"

sys.stderr.write("""Imported modules for python bioinformatics analysis
Author: Jinyang Zhang
Email: zhangjinyang@biols.ac.cn
Version: 1.0\n""")
sys.stderr.write(sys.version + '\n')
