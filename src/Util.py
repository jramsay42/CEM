""" CEM Utility routines. """

from math import sqrt

from scipy.constants import c

def calculate_grid_resolution(mu_r, eps_r):
    """ 
        Calculates the required grid resolution. 
    """
    return sqrt(mu_r * eps_r)

def calculate_minimum_wavelength(max_frequency, mu_r, eps_r):
    """
        Calculates the smallest possible wavelength.
    """
    return c / (max_frequency * calculate_grid_resolution(mu_r, eps_r))

def calculate_wavelength_step(max_frequency, mu_r, eps_r, num_lambda):
    """ 
        Calculates the wavelength step required
    """
    return calculate_minimum_wavelength(max_frequency, mu_r, eps_r) / num_lambda