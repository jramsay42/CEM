"""
    Modelling a wave incident upon a 1 foot thick dielectric slab.
    mu_r = 2, eps_r = 6 surrounded by air.
"""

import math

import numpy as np
import matplotlib.pyplot as plot
from scipy.constants import mu_0
from scipy.constants import epsilon_0
from scipy.constants import c

from Util import calculate_wavelength_step

""" CONSTANTS """
SLAB_THICKNESS = 0.3048 #m
MU_R = 2
EPS_R = 6
MAX_FREQUENCY = 1e9
NUMBER_WAVELENGTHS_DESIRED = 20
SPACING = 10
N_MAX = math.sqrt(MU_R * EPS_R)

""" PARAMETER CALCULATION """
NUMBER_STEPS = math.ceil(SLAB_THICKNESS / calculate_wavelength_step(MAX_FREQUENCY, MU_R, EPS_R, NUMBER_WAVELENGTHS_DESIRED))
Z_STEP_SIZE = SLAB_THICKNESS / NUMBER_STEPS

""" DETERMINING NUMBER OF CELLS """
SIZE_OF_GRID = NUMBER_STEPS + 2 * SPACING + 3 #Some extra cells for source and recording values
MU_Z = np.arange(0, SIZE_OF_GRID * Z_STEP_SIZE, Z_STEP_SIZE) #mu position vector
EPS_Z = np.arange(0, SIZE_OF_GRID * Z_STEP_SIZE, Z_STEP_SIZE) #eps position vector

""" COMPUTE POSITION OF THE SLAB """
SLAB_START = 2 + SPACING + 1 #Index the slab starts at
SLAB_END = SLAB_START + NUMBER_STEPS - 1 #Index the slab ends at

""" MATERIAL GRID """
MU_Z[SLAB_START:SLAB_END] = MU_R #Put in the slab
EPS_Z[SLAB_START: SLAB_END] = EPS_R  #Put in the slab
N_Z = MU_Z.size #How many Z cells there are for the slab

""" SETUP FDTD """
T_STEP_SIZE = Z_STEP_SIZE/ (2 * c)
TAU = 1 / (2 * MAX_FREQUENCY)
T_0 = 6 * TAU #Rule of thumb
PROPAGATION_TIME = N_MAX * N_Z * Z_STEP_SIZE / c
T = 12 * TAU + 5 * PROPAGATION_TIME #Rule of thumb
N_T = math.ceil(T / T_STEP_SIZE)

T_VECTOR = np.arange(0, T_STEP_SIZE * N_T, T_STEP_SIZE)

""" COMPUTE SOURCE FUNCTIONS"""
TIME_STAGGERING = 1.0 * Z_STEP_SIZE / (2 * c) + T_STEP_SIZE / 2
E_Y_SRC = np.exp(-np.square(((T_VECTOR - T_0) / TAU)))
H_X_SRC = -np.sqrt(1) * np.exp(-np.square((T_VECTOR - T_0 + TIME_STAGGERING) / TAU))

""" INITIALIZE THE FOURIER TRANSFORMS """
N_FREQ = 100
F_VECTOR = np.arange(0, MAX_FREQUENCY, MAX_FREQUENCY / N_FREQ)
K = np.exp(-1j * 2 * np.pi * T_STEP_SIZE * F_VECTOR)
REF = np.zeros(N_FREQ)
TRN = np.zeros(N_FREQ)
SRC = np.zeros(N_FREQ)

""" RUN FDTD """


""" VISUALIZE SIMULATION """


