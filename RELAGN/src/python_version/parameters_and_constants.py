# Import global modules. #
import astropy.units as units, numpy as np

# Declare constants. #
m_sun_cgs = 1.99e33 * units.g  # Solar mass in cgs.
c_cgs = 3e10 * units.cm * units.s ** (-1)  # Speed of light in cgs.
G_cgs = 6.67e-8 * units.cm ** 3 * units.g ** (-1) * units.s ** (-2)  # Gravitational constant in cgs.

# Declare tolerance parameters. #
rtol, atol = 1e-10, 1e-10
max_iteration = 20  # Number of maximum iterations allowed.

# Declare accretion disc parameters. #
num_bins = 5000  # Number of radial bins the intra-ISCO and extra-ISCO regions are divided into.
num_steps = 10  # Factor to adjust the accretion disc's extent when calculating Toomre Q parameter.