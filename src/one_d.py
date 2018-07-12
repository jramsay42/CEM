"""
    1-D FD-TD

    E and H fields are on a staggered Yee-grid.
    E and H are staggered in time by DELTA_T/2
    H is normalised to sqrt(mu/epsilon)H

    Thus,
        Curl (E) = -mu_r/C_0 dH/dt
        Curl (H) = e_r/C_0 dE/dt

    Expanding this (with linear, diagonal tensors),

        dEz/dy - dEy/dz = -mu_xx/C_0 dHx/dt
        dEx/dz - dEz/dx = -mu_yy/C_0 dHy/dt
        dEy/dx - dEx/dy = -mu_zz/C_0 dHz/dt

        dHz/dy - dHy/dz = eps_xx/C_0 dEx/dt
        dHx/dz - dHz/dx = eps_yy/C_0 dEy/dt
        dHy/dx - dHx/dy = eps_zz/C_0 dEz/dt
 """

import numpy as np
import matplotlib.pyplot as plot
from scipy.constants import mu_0
from scipy.constants import epsilon_0
from scipy.constants import c

"""
    For one dimension, let's assume d/dy and d/dx terms are zero. I.e.
    variation is in the z-direction.

    Thus,

        dEy/dx = mu_xx/C_0 dHx/dt
        dEx/dz = -mu_yy/C_0 dHy/dt
        0      = mu_zz/C_0 dHz/dt

        -dHy/dz = eps_xx/C_0 dEx/dt
        dHx/dz  = eps_yy/C_0 dEy/dt
        0       = eps_zz/C_0 dEz/dt

    Our FD update equations are then,

        E_y(k, t + DELTA_T) = E_y(k, t) + (C_0 * DELTA_T) / eps_yy(k) * (H_x(k, t + DELTA_T / 2) - H_x(k-1, t + DELTA_T / 2)) / DELTA_Z

        H_x(k, t + DELTA_T / 2) = H_x(k, t - DELTA_T / 2) + (C_0 * DELTA_T) / mu_xx(k) * (E_y(k+1, t) - E_y(k, t)) / DELTA_Z

        E_x(k, t + DELTA_T) = E_x(k, t) + (C_0 * DELTA_T) / eps_xx(k) * -(H_y(k, t + DELTA_T / 2) - H_y(k-1, t + DELTA_T / 2)) / DELTA_Z

        H_y(k, t + DELTA_T / 2) = H_y(k, t - DELTA_T / 2) + (C_0 * DELTA_T) / mu_yy(k) * -(E_y(k+1, t) - E_y(k, t)) / DELTA_Z

"""

#Init simulation
#Init parameters
#Compute grid resolution
#Build device on grid
#Compute time step
#Compute source
#Init FTs
#Compute update coefficients
#Init fields to zero
#Init boundary terms to zero

DELTA_T = 0.001
DELTA_Z = 0.001

T = np.arange(0, 1, DELTA_T) #time vector
Z = np.arange(0, 1, DELTA_Z) #position vector

EPS_YY = np.squeeze(np.full((1, Z.size), 1)) #relative permitivity
MU_XX = np.squeeze(np.full((1, Z.size), 1)) #relative permeability

e_x = np.zeros(Z.size)
e_y = np.zeros(Z.size)
h_x = np.zeros(Z.size)
h_y = np.zeros(Z.size)

M_EY = c * DELTA_T / EPS_YY #Ey update vector
M_HX = c * DELTA_T / MU_XX #Hx update vector


#Main FD-TD loop
for t_idx in range(T.size):
#Record H at boundary

#Update H (Dirchlet Boundary Conditions)
    for z_idx in range(Z.size - 1):
        h_x[z_idx] = h_x[z_idx] + M_HX[z_idx] * (e_y[z_idx + 1] - e_y[z_idx]) / DELTA_Z 
    h_x[Z.size - 1] = h_x[Z.size - 1] + M_HX[Z.size - 1] * (0 - e_y[Z.size - 1])

#Handle H source

#Record E at boundary

#Update E (Dirchlet Boundary Conditions)
    e_y[1] = e_y[1] + M_EY[1] * (h_x[1] - 0)
    for z_idx in range(1, Z.size - 1):
        e_y[z_idx] = e_y[z_idx] + M_EY[z_idx] * (h_x[z_idx + 1] - h_x[z_idx]) / DELTA_Z

#Handle E source

#Update FTs

#Visualize simulation