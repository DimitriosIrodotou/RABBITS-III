# Import global modules. #
import astropy.units as units, numpy as np, scipy.integrate as integrate, time, os, inspect
# Import local modules. #
import profiles, properties, utilities
from parameters_and_constants import c_cgs, G_cgs, m_sun_cgs, max_iteration, num_steps

# Extract the global path where the scripts are stored and create a global path to save data. #
global_scripts_path = os.path.dirname(os.path.realpath(__file__))
global_data_path = global_scripts_path + '/data/'


def r_grav(mass, verbose=False):
    """
    Calculate the gravitational radius based on Eq. .
    :param mass: Black hole mass in gram.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Gravitational radius in centimeter.
    """
    local_time = time.time()
    print('Entering "r_grav" from "radii.py" \n' if verbose else '', end='')

    if mass.unit == 'g':  # Check if 'mass' is in gram.
        radius = G_cgs * mass * c_cgs ** (-2)
    else:
        print('Terminating "r_grav"! "mass" is in units of "', mass.unit, '" instead of "g"', sep=''), exit()

    print('Exiting "r_grav" from "radii.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
    return radius if radius.unit == 'cm' else (print('Terminating "r_grav"! "radius" is in units of "', radius.unit,
                                                     '" instead of "cm"', sep=''), exit())


def r_photon(mass, spin, rotation, vectorized=False, verbose=False):
    """
    Calculate the radius of the photon orbit on Eq. 2.18 from 1972ApJ...178..347B.
    :param mass: Black hole mass in gram.
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Radius in centimeter.
    """
    local_time = time.time()
    print('Entering "r_photon" from "regions.py" \n' if verbose else '', end='')

    # 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. #
    if rotation == 'co-rotation':
        alignment = +1
    elif rotation == 'counter-rotation':
        alignment = -1
    else:
        print('Terminating "r_isco"! "rotation" is "', rotation, '" instead of "co-rotation" or "counter-rotation"',
              sep=''), exit()

    # If 'spin' is an array of spins, use vectorization. #
    if vectorized:
        condition_spin = ((0 <= spin) & (spin <= 0.998)).all()
    else:
        condition_spin = 0 <= spin <= 0.998

    if condition_spin:  # Check if 'spin' is between 0 and 0.998.
        radius = 2 * r_grav(mass) * (1 + np.cos(2 / 3 * np.arccos(-1 * np.sign(alignment) * spin)))
    else:
        print('Terminating "r_photon"! "spin" is "', spin,
              '" instead of greater than or equal to 0 and less than or equal to 0.998.', sep=''), exit()

    print('Exiting "r_photon" from "regions.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
    return radius if radius.unit == 'cm' else (print('Terminating "r_photon"! "radius" is in units of "',
                                                     radius.unit, '" instead of "cm"', sep=''), exit())


def r_isco(mass, spin, rotation, vectorized=False, verbose=False):
    """
    Calculate the radius of the innermost stable circular orbit (ISCO) based on Eq. .
    :param mass: Black hole mass in gram.
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: ISCO radius in centimeter.
    """
    local_time = time.time()
    print('Entering "r_isco" from "radii.py" \n' if verbose else '', end='')

    # 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. #
    if rotation == 'co-rotation':
        alignment = +1
    elif rotation == 'counter-rotation':
        alignment = -1
    else:
        print('Terminating "r_isco"! "rotation" is "', rotation, '" instead of "co-rotation" or "counter-rotation"',
              sep=''), exit()

    # If 'spin' is an array of spins, use vectorization. #
    if vectorized:
        condition_spin = ((0 <= spin) & (spin <= 0.999)).all()
    else:
        condition_spin = 0 <= spin <= 0.998

    if condition_spin:  # Check if 'spin' is between 0 and 0.998.
        z_1 = 1 + (1 - spin ** 2) ** (1 / 3) * ((1 + spin) ** (1 / 3) + (1 - spin) ** (1 / 3))
        z_2 = (3 * spin ** 2 + z_1 ** 2) ** (1 / 2)
        z = 3 + z_2 - np.sign(alignment) * ((3 - z_1) * (3 + z_1 + 2 * z_2)) ** (1 / 2)
        radius = r_grav(mass) * z
    else:
        print('Terminating "r_isco"! "spin" is "', spin,
              '" instead of greater than or equal to 0 and less than or equal to 0.998.', sep=''), exit()

    print('Exiting "r_isco" from "radii.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
    return radius if radius.unit == 'cm' else (print('Terminating "r_isco"! "radius" is in units of "', radius.unit,
                                                     '" instead of "cm"', sep=''), exit())


def r_outermost(mass, verbose=False):
    """
    Calculate the radius of the outermost edge of the accretion disc based on Eq. .
    :param mass: Black hole mass in gram.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Outermost radius in centimeter.
    """
    local_time = time.time()
    print('Entering "r_outermost" from "radii.py" \n' if verbose else '', end='')

    if mass.unit == 'g':  # Check if 'mass' is in gram.
        radius = 10 ** 15.78 * (mass / (1e9 * m_sun_cgs)) ** 0.8 * units.cm
    else:
        print('Terminating "r_outermost"! "mass" is in units of "', mass.unit, '" instead of "g"', sep=''), exit()

    print('Exiting "r_outermost" from "radii.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')

    # Extract this function's name and load the 'extent_factor'. #
    function_name = inspect.currentframe().f_code.co_name
    local_data_path = global_data_path + function_name + '/'
    extent_factor = np.load(local_data_path + 'extent_factor.npy')

    return extent_factor * radius if radius.unit == 'cm' else (
        print('Terminating "r_outermost"! "radius" is in units of "', radius.unit, '" instead of "cm"', sep=''), exit())


def accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Calculate the radial extent of the accretion disc for given black hole and accretion disc properties in
    parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: None
    """
    extent_factor, iterations, flag = 1, 0, 'repeat'  # Declare local variables.

    # Check if the path to save the data exists. If not, then create the corresponding directory and save the initial
    # 'extent_factor'. #
    extent_data_path = global_data_path + 'r_outermost/'
    if not os.path.exists(extent_data_path): os.makedirs(extent_data_path)
    np.save(extent_data_path + 'extent_factor', extent_factor)

    while flag == 'repeat' or extent_factor == 1:
        # Calculate, save, and load the 'ranges_of_validity' of different regimes. #
        profiles.ranges_of_validity(viscosity, bh_mass, accretion_rate, spin, rotation)
        ranges_data_path = global_data_path + 'ranges_of_validity/'
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        outermost_regime = np.load(ranges_data_path + 'outermost_regime.npy').item()

        # Calculate the 'toomre_parameter' inside the 'ranges_of_validity' of the 'outermost_regime'. Select only the
        # outermost radial bins (i.e. entry -1 in 'ranges_of_validity[outermost_regime]') which are expected to be
        # gravitationally unstable. #
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
        x_bins = (ranges_of_validity[outermost_regime][-1] / r_grav(bh_mass_cgs).value) \
                 ** (1 / 2)  # Dimensionless bins.
        toomre_parameter = properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin,
                                                       rotation, outermost_regime)

        # Remove any NaN entries that will prevent the accretion disc from reaching a self-gravitational stability. #
        toomre_parameter = toomre_parameter[~np.isnan(toomre_parameter)]

        # Check if the first and last 'toomre_parameter' is below one. If yes, that means the 'extent' has
        # increased more than enough, so set it to so that the next 'r_outermost' is 90% of the current, hence the
        # current 'outermost_regime' does not exist when calculating the new 'ranges_of_validity' in the next do-while
        # iteration. #
        if toomre_parameter[0] < 1 and toomre_parameter[-1] < 1:
            extent_factor, flag = 1, 'repeat'
            np.save(extent_data_path + 'extent_factor', extent_factor)
            extent_factor = 0.90 * max(ranges_of_validity[outermost_regime][-1]) / r_outermost(bh_mass_cgs)
            np.save(extent_data_path + 'extent_factor', extent_factor.value)

        # Check if the first and last 'toomre_parameter' are between one. If yes, that means the  instability happened
        # between radial bins so 'interpolate' the solution. #
        if toomre_parameter[-1] < 1 < toomre_parameter[0]:
            index = np.argmin(np.abs(toomre_parameter - 1))
            # Check if the 'toomre_parameter' that is closest to one (i.e. stability threshold) is above. If yes, then
            # use the value after 'index' (should be below one) to interpolate and find which 'r' results in
            # 'toomre_parameter' of one. If not, use the value before (should be above one).
            if toomre_parameter[index].value > 1:
                r = utilities.interpolate(ranges_of_validity[outermost_regime][-1][index + 1],
                                          toomre_parameter[index + 1].value,
                                          ranges_of_validity[outermost_regime][-1][index],
                                          toomre_parameter[index].value)
            elif toomre_parameter[index].value < 1:
                r = utilities.interpolate(ranges_of_validity[outermost_regime][-1][index - 1],
                                          toomre_parameter[index - 1].value,
                                          ranges_of_validity[outermost_regime][-1][index],
                                          toomre_parameter[index].value)
            else:
                r = ranges_of_validity[outermost_regime][-1][index]
            extent_factor = 1
            np.save(extent_data_path + 'extent_factor', extent_factor)
            extent_factor, flag = r / r_outermost(bh_mass_cgs).value, "done"
            np.save(extent_data_path + 'extent_factor', extent_factor)

        # Check if the last 'toomre_parameter' is above one. If yes, that means the 'extent' needs to be increased
        # by a factor of 'num_steps'. #
        if toomre_parameter[-1] > 1:
            extent_factor, flag = extent_factor * num_steps, 'repeat'
            np.save(extent_data_path + 'extent_factor', extent_factor)

        iterations += 1
        if iterations == max_iteration:
            print('Terminating "accretion_disc_extent"! "max_iteration" has been reached.', sep=''), exit()

    # Calculate, save the latest 'ranges_of_validity'. #
    profiles.ranges_of_validity(viscosity, bh_mass, accretion_rate, spin, rotation)
