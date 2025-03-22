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
viscosity = 0.1  # Dimensionless viscosity parameter.
rotation = 'counter-rotation'  # Accretion disc orientation relative to the black hole spin vector.
num_bins = 5000  # Number of radial bins the intra-ISCO and extra-ISCO regions are divided into.
num_steps = 10  # Factor to adjust the accretion disc's extent when calculating Toomre Q parameter.

# Declare black hole parameters. #
bh_mass = 1e7  # Dimensionless black hole mass.
spin = 0.5  # Dimensionless spin parameter.
accretion_rate = 1e5  # Dimensionless black hole accretion rate.

# List of dimensionless spins. #
step = 0.001
spins = np.arange(start=0.000, stop=0.998 + step, step=step)
spins = spins[spins <= 0.998]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# RELAGN SED example with values corresponding to a snapshot of an isolated galaxy formation simulation
# spin = 5.60032948001510155756e-01  # Dimensionless
# BH_Mass = 8.93032240454057365218e-04  # In 1e10 Msun / h.
# BH_Mdot = 3.12010333482085414070e-04  # In 1e10 Msun / 9.8e8 yr.
# Edd_Mdot = 3.12010333482085448764e-02  # In 1e10 Msun / 9.8e8 yr.
#
# hubble = 0.71
# bh_mass_astro = BH_Mass * 1e10 / hubble  # In Msun.
# bh_mass_cgs = bh_mass_astro * m_sun_cgs  # In g.
# bh_mass = bh_mass_astro / 3  # Dimensionless black hole mass.
#
# bh_mdot_astro = BH_Mdot * 1e10 / 9.8e8  # In Msun / yr.
# bh_mdot_cgs = bh_mdot_astro * m_sun_cgs / (units.yr.to('s') * units.s)  # In g /  s
# accretion_rate = bh_mdot_cgs / (1e17 * units.g / units.s)  # Dimensionless black hole accretion rate.
#
# bh_mdot_edd = (Edd_Mdot * 1e10 / 9.8e8) * m_sun_cgs / (units.yr.to('s') * units.s)
