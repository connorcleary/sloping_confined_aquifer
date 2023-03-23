import sys
import os
from tempfile import TemporaryDirectory

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from flopy.utils.gridgen import Gridgen

# temporary directory
temp_dir = TemporaryDirectory()
workspace = temp_dir.name

