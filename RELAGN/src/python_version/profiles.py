# Import global modules. #
import matplotlib.pyplot as plt, numpy as np, os, time, inspect
# Import local modules. #
import properties, radii, utilities
from parameters_and_constants import m_sun_cgs, num_bins, atol, rtol

# Extract the global path where the scripts are stored and create global paths to save data and plots. #
global_scripts_path = os.path.dirname(os.path.realpath(__file__))
global_data_path = global_scripts_path + '/data/'


def ranges_of_validity(viscosity, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Calculate the ranges of validity for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: None
    """
    local_time = time.time()  # Start the time to calculate the local runtime.

    # Extract this function's name and create a local path to save the data. #
    function_name = inspect.currentframe().f_code.co_name
    local_data_path = global_data_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the data exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_data_path): os.makedirs(local_data_path)

    # Declare a dictionary to store the ranges of validity. #
    ranges_of_validity = {'intra-ISCO': [], 'Gas-ES': [], 'Rad-ES': [], 'Gas-FF': []}

    # Calculate evenly log-spaced 'radial_bins' from 'r_photon' to 'r_isco'. #
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    r_isco = radii.r_isco(bh_mass_cgs, spin, rotation).value
    intra_radial_bins = np.logspace(start=np.log10(radii.r_photon(bh_mass_cgs, spin, rotation).value),
                                    stop=np.log10(r_isco), num=num_bins)
    ranges_of_validity['intra-ISCO'].append(intra_radial_bins[1:-1])  # 'r_photon' and 'r_isco' are excluded.

    # Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost'. #
    r_outermost = radii.r_outermost(bh_mass_cgs).value
    radial_bins = np.logspace(start=np.log10(r_isco), stop=np.log10(r_outermost), num=num_bins)
    x_bins = (radial_bins / radii.r_grav(bh_mass_cgs).value) ** (1 / 2)  # Dimensionless radial bins.

    # Check for every 'regime' if its 'ratios' (i.e. 'pressure_ratio' and 'opacity_ratio') are both less than one. If
    # yes, then store the corresponding 'radial_bins'. #
    mask_GES, = np.where(
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[0] < 1)
        &
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-ES')[1] < 1))
    mask_RES, = np.where(
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[0] < 1)
        &
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Rad-ES')[1] < 1))
    mask_GFF, = np.where(
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[0] < 1)
        &
        (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, 'Gas-FF')[1] < 1))

    # If two or more regimes are valid in the same radial bin, split them based on their 'pressure_ratio'. #
    mask_GES, mask_RES, mask_GFF = utilities.resolve_overlapping_regimes(viscosity, bh_mass, accretion_rate, spin,
                                                                         rotation, x_bins, mask_GES, mask_RES, mask_GFF)

    # Before 'resolve_lone_regimes', 'resolve_empty_regimes' (i.e. assign empty bins to adjacent regimes). #
    mask_GES, mask_RES, mask_GFF = utilities.resolve_empty_regimes(viscosity, bh_mass, accretion_rate, spin, rotation,
                                                                   x_bins, mask_GES, mask_RES, mask_GFF)

    # If a regime is only valid in one radial bin, replace it with an adjacent valid regime. #
    mask_GES, mask_RES, mask_GFF = utilities.resolve_lone_regimes(viscosity, bh_mass, accretion_rate, spin, rotation,
                                                                  x_bins, mask_GES, mask_RES, mask_GFF)

    # Final sanity check, more info in their print statement. #
    if len(mask_GES) + len(mask_RES) + len(mask_GFF) != num_bins:
        print('Terminating "ranges_of_validity"! Radial bin issues still exist after all "resolve_"', sep=''), exit()

    # Loop over all regimes and record their 'ranges_of_validity'. #
    for regime, mask_regime in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], [mask_GES, mask_RES, mask_GFF]):
        if len(mask_regime) > 0:
            # Split 'mask_regime' into sub-arrays of sequential numbers when another 'regime' interferes. #
            split_points = np.flatnonzero(np.abs(np.diff(mask_regime)) != 1) + 1
            split_arrays = np.split(mask_regime, split_points)
            for i in range(len(split_arrays)):
                ranges_of_validity[regime].append(radial_bins[split_arrays[i]])

            # Check for every 'regime' if the first radial bin and 'r_isco', and the last radial bin and 'r_outermost'
            # are equal within a tolerance. If yes, assign 'r_isco' as the starting and 'r_outermost' as the stopping
            # radii. #
            if np.isclose(radial_bins[mask_regime[0]], r_isco, rtol=rtol, atol=atol):
                ranges_of_validity[regime][0][0] = r_isco
                np.save(local_data_path + 'isco_regime', regime)
            if np.isclose(radial_bins[mask_regime[-1]], r_outermost, rtol=rtol, atol=atol):
                ranges_of_validity[regime][-1][-1] = r_outermost
                np.save(local_data_path + 'outermost_regime', regime)

    # Save the data and exit the script. #
    np.save(local_data_path + 'ranges_of_validity', ranges_of_validity)
    print('Exiting "ranges_of_validity" from "profiles.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
