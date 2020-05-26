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
from utils import *

# Plotting parameters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Arial"

print('Imported modules for python bioinformatics analysis')
print('Author: JinyangZhang')
print('Version:')
print(sys.version)


