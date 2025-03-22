# Import global modules. #
import time
# Import local modules. #
import combinations, profiles, radii
from parameters_and_constants import viscosity, bh_mass, accretion_rate, spin, spins, rotation

start_time = time.time()  # Start the time to calculate the global runtime.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Build an accretion disc and store the ranges of validity of each regime (required for all profiles below!). #
radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

# Plot radial profiles of accretion disc properties. #
profiles.surface_density(viscosity, bh_mass, accretion_rate, spin, rotation)
profiles.temperature(viscosity, bh_mass, accretion_rate, spin, rotation)
profiles.stability(viscosity, bh_mass, accretion_rate, spin, rotation)
profiles.regional_solutions(viscosity, bh_mass, accretion_rate, spins, rotation)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

combinations.surface_densities_viscosities([0.01, 0.1, 1], bh_mass, accretion_rate, spin, rotation)
combinations.surface_densities_bh_masses(viscosity, [1e6, 1e7, 1e8], accretion_rate, spin, rotation)
combinations.surface_densities_accretion_rates(viscosity, bh_mass, [1e4, 1e5, 1e6], spin, rotation)
combinations.surface_densities_spins(viscosity, bh_mass, accretion_rate, [0, 0.5, 0.9], rotation)

combinations.stabilities_viscosities([0.01, 0.1, 1], bh_mass, accretion_rate, spin, rotation)
combinations.stabilities_bh_masses(viscosity, [1e6, 1e7, 1e8], accretion_rate, spin, rotation)
combinations.stabilities_accretion_rates(viscosity, bh_mass, [1e4, 1e5, 1e6], spin, rotation)
combinations.stabilities_spins(viscosity, bh_mass, accretion_rate, [0, 0.5, 0.9], rotation)

combinations.correction_functions(viscosity, bh_mass, accretion_rate, spins)
combinations.parameter_space([0.01, 0.1, 1], [1e6, 1e7, 1e8], [1e4, 1e5, 1e6], [0, 0.5, 0.9])

print('Finished "main.py" in %.4s s' % (time.time() - start_time))  # Print the global runtime.
