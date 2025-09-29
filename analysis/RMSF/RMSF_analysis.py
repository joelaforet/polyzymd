'''Functions used to compute the RMSF of various portions of the protein in a simulation'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'


import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms, align

from biopandas.pdb import PandasPdb

import json
from collections import defaultdict
from tqdm import tqdm
from matplotlib import cm
from matplotlib.cm import get_cmap

# Suppress some MDAnalysis warnings about writing PDB files
import warnings
warnings.filterwarnings('ignore')

from typing import Union


