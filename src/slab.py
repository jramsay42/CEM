"""
    Modelling a wave incident upon a 1 foot thick dielectric slab.
    mu_r = 2, eps_r = 6 surrounded by air.
"""

from math import sqrt

import numpy as np
import matplotlib.pyplot as plot
from scipy.constants import mu_0
from scipy.constants import epsilon_0
from scipy.constants import c

from Util import *

print(calculate_wavelength_step(1e9, 2, 6, 20))
