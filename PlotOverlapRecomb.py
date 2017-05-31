# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:51:25 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes in regions of high recombination

# usage PlotOverlapRecomb.py

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
import json
import random
import copy
import sys
import os
import math
import numpy as np
from scipy import stats
from HsaNestedGenes import *










chr     start   end     name
chr1    1696659 1821625 CS_chr1_1
chr1    2491060 2709164 CS_chr1_2
chr1    3709412 3781096 CS_chr1_3
chr1    5940650 6048126 CS_chr1_4
chr1    6358604 6481842 CS_chr1_5
chr1    6705356 6783967 CS_chr1_6
chr1    6832235 6956195 CS_chr1_7
chr1    7256943 7310351 CS_chr1_8
chr1    7986759 8168564 CS_chr1_9
chr1    8429598 8495519 CS_chr1_10
chr1    8522553 8878809 CS_chr1_11
chr1    9594670 9651279 CS_chr1_12
chr1    10006172        10229158        CS_chr1_13
chr1    10282166        10479678        CS_chr1_14
chr1    10491209        10561282        CS_chr1_15
chr1    10570721        10636459        CS_chr1_16
chr1    11128654        11319587        CS_chr1_17
chr1    12302072        12564127        CS_chr1_18
chr1    14005069        14109114        CS_chr1_19
chr1    15814186        15907908        CS_chr1_20
chr1    15922371        15988900        CS_chr1_21
