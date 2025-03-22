# Import global modules. #
import astropy.units as units, numpy as np, time, os
# Import local modules. #
import radii
from parameters_and_constants import c_cgs, G_cgs, m_sun_cgs

# Extract the global path where the scripts are stored and create a global path to save data. #
global_scripts_path = os.path.dirname(os.path.realpath(__file__))
global_data_path = global_scripts_path[0:-7] + 'data/'


def correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the relativistic correction functions for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Correction functions.
    """
    local_time = time.time()
    print('Entering "correction_functions" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "correction_functions"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. #
        if rotation == 'co-rotation':
            alignment = +1
        elif rotation == 'counter-rotation':
            alignment = -1
        else:
            print('Terminating "correction_functions"! "rotation" is "', rotation,
                  '" instead of "co-rotation" or "counter-rotation"', sep=''), exit()

        # If 'spin' is an array of spins use vectorization. #
        if vectorized:
            condition_spin = ((0 <= spin) & (spin <= 0.999)).all()
        else:
            condition_spin = 0 <= spin <= 0.998

        if condition_spin:  # Check if 'spin' is between 0 and 0.998.
            # Dimensionless radius of the ISCO. #
            x_isco = (radii.r_isco(3 * mass * m_sun_cgs, spin, rotation, vectorized) /
                      radii.r_grav(3 * mass * m_sun_cgs)) ** (1 / 2)

            # Calculate the correction functions and their value at the ISCO (_isco). #
            cal_c = 1 - 3 * x ** (-2) + np.sign(alignment) * 2 * spin * x ** (-3)
            cal_c_isco = 1 - 3 * x_isco ** (-2) + np.sign(alignment) * 2 * spin * x_isco ** (-3)
            cal_d = 1 - 2 * x ** (-2) + spin ** 2 * x ** (-4)
            cal_d_isco = 1 - 2 * x_isco ** (-2) + spin ** 2 * x_isco ** (-4)
            cal_f = np.sign(alignment) * (1 - np.sign(alignment) * 2 * spin * x ** (-3) + spin ** 2 * x ** (-4))
            cal_f_isco = np.sign(alignment) * (
                    1 - np.sign(alignment) * 2 * spin * x_isco ** (-3) + spin ** 2 * x_isco ** (-4))
            cal_g = 1 - 2 * x ** (-2) + np.sign(alignment) * spin * x ** (-3)
            cal_g_isco = 1 - 2 * x_isco ** (-2) + np.sign(alignment) * spin * x_isco ** (- 3)
            cal_r = cal_c ** (-1) * cal_f ** 2 - spin ** 2 * x ** (-2) * (cal_c ** (-1 / 2) * cal_g - 1)
            cal_r_isco = cal_c_isco ** (-1) * cal_f_isco ** 2 - spin ** 2 * x_isco ** (-2) * (
                    cal_c_isco ** (-1 / 2) * cal_g_isco - 1)

            x_1 = 2 * np.cos(np.arccos(np.sign(alignment) * spin) / 3 - np.pi / 3)
            x_2 = 2 * np.cos(np.arccos(np.sign(alignment) * spin) / 3 + np.pi / 3)
            x_3 = -2 * np.cos(np.arccos(np.sign(alignment) * spin) / 3)

            if regime == 'Gas-ES':
                h_0 = 0.00184101952815808 * \
                      viscosity ** (1 / 8) * mass ** (-3 / 8) * accretion_rate ** (1 / 4) * x_isco ** (1 / 8) * \
                      cal_c_isco ** (-1 / 8) * cal_r_isco ** (-1 / 2)
            elif regime == 'Rad-ES':
                h_0 = 0.002
            elif regime == 'Gas-FF':
                h_0 = 0.0013166295323149043 * \
                      viscosity ** (1 / 17) * mass ** (-5 / 17) * accretion_rate ** (3 / 17) * x_isco ** (5 / 17) * \
                      cal_c_isco ** (-1 / 17) * cal_d_isco ** (-1 / 34) * cal_r_isco ** (-8 / 17)

            elif regime == 'intra-ISCO':
                # Local path to load the 'outermost_regime' from. #
                ranges_data_path = global_data_path + 'ranges_of_validity/'
                isco_regime = np.load(ranges_data_path + 'isco_regime.npy')

                if isco_regime == 'Gas-ES':
                    h_0 = 0.00184101952815808 * \
                          viscosity ** (1 / 8) * mass ** (-3 / 8) * accretion_rate ** (1 / 4) * x_isco ** (1 / 8) * \
                          cal_c_isco ** (-1 / 8) * cal_r_isco ** (-1 / 2)
                elif isco_regime == 'Rad-ES':
                    h_0 = 0.002
                elif isco_regime == 'Gas-FF':
                    h_0 = 0.0013166295323149043 * \
                          viscosity ** (1 / 17) * mass ** (-5 / 17) * accretion_rate ** (3 / 17) * x_isco ** (5 / 17) * \
                          cal_c_isco ** (-1 / 17) * cal_d_isco ** (-1 / 34) * cal_r_isco ** (-8 / 17)
                else:
                    print('Terminating "correction_functions"! "isco_regime" is "', regime,
                          '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

            else:
                print('Terminating "correction_functions"! "regime" is "', regime,
                      '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

            cal_p_isco = 2 ** (-1 / 2) * viscosity * x_isco * h_0 * cal_d_isco ** (1 / 2) * cal_r_isco ** (1 / 2)

            cal_p = cal_p_isco + x - x_isco - np.sign(alignment) * (3 * spin / 2) * np.log(x / x_isco) \
                    - (3 * (x_1 - np.sign(alignment) * spin) ** 2) / (x_1 * (x_1 - x_2) * (x_1 - x_3)) \
                    * np.log((x - x_1) / (x_isco - x_1)) \
                    - (3 * (x_2 - np.sign(alignment) * spin) ** 2) / (x_2 * (x_2 - x_3) * (x_2 - x_1)) \
                    * np.log((x - x_2) / (x_isco - x_2)) \
                    - (3 * (x_3 - np.sign(alignment) * spin) ** 2) / (x_3 * (x_3 - x_1) * (x_3 - x_2)) \
                    * np.log((x - x_3) / (x_isco - x_3))
        else:
            print('Terminating "correction_functions"! "spin" is "', spin,
                  '" instead of greater than or equal to 0 and less than or equal to 0.998.', sep=''), exit()

        print('Exiting "correction_functions" from "properties.py" after %.4s '
              % (time.time() - local_time) + 's \n' if verbose else '', end='')
        return cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco


def ratios(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the pressure and opacity ratios for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Pressure and opacity ratios.
    """
    local_time = time.time()
    print('Entering "ratios" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "ratios"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            pressure_ratio = 68.98449619493478 * \
                             viscosity ** (1 / 10) * mass ** (-7 / 10) * accretion_rate ** (4 / 5) * x ** (-29 / 10) * \
                             cal_c ** (-9 / 10) * cal_d ** (1 / 10) * cal_p ** (4 / 5) * cal_r ** (-1 / 2)
            opacity_ratio = 4.359450454931594e-6 * mass * accretion_rate ** (-1) * x ** 4 \
                            * cal_c * cal_p ** (-1) * cal_r ** (1 / 2)
        elif regime == 'Rad-ES':
            pressure_ratio = 0.000025300042906895175 * \
                             viscosity ** (-1 / 4) * mass ** (7 / 4) * accretion_rate ** (-2) * x ** (29 / 4) * \
                             cal_c ** (9 / 4) * cal_d ** (-1 / 4) * cal_p ** (-2) * cal_r ** (5 / 4)
            opacity_ratio = 2.1927664368169804e-8 * \
                            viscosity ** (-1 / 8) * mass ** (15 / 8) * accretion_rate ** (-2) * x ** (61 / 8) * \
                            cal_c ** (17 / 8) * cal_d ** (-1 / 8) * cal_p ** (-2) * cal_r ** (9 / 8)
        elif regime == 'Gas-FF':
            pressure_ratio = 0.26699048406463 * \
                             viscosity ** (1 / 10) * mass ** (-1 / 4) * accretion_rate ** (7 / 20) * x ** (-11 / 10) * \
                             cal_c ** (-9 / 20) * cal_d ** (1 / 10) * cal_p ** (7 / 20) * cal_r ** (-11 / 40)
            opacity_ratio = 478.9433271561746 * mass ** (-1 / 2) * accretion_rate ** (1 / 2) * x ** (-2) * \
                            cal_c ** (-1 / 2) * cal_p ** (1 / 2) * cal_r ** (-1 / 4)
        else:
            print('Terminating "ratios"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        print('Exiting "ratios" from "properties.py" after %.4s ' % (
                time.time() - local_time) + 's \n' if verbose else '', end='')
        return pressure_ratio, opacity_ratio


def densities(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the density profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Density in gram centimeter^-3.
    """
    local_time = time.time()
    print('Entering "densities" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "densities"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            density = 8.1 * units.g * units.cm ** (-3) * \
                      viscosity ** (-7 / 10) * mass ** (-11 / 10) * accretion_rate ** (2 / 5) * x ** (-37 / 10) * \
                      cal_c ** (3 / 10) * cal_d ** (-7 / 10) * cal_p ** (2 / 5) * cal_r ** (1 / 2)
        elif regime == 'Rad-ES':
            density = 2.5e-5 * units.g * units.cm ** (-3) * viscosity ** (-1) * mass * accretion_rate ** (-2) * x ** 5 * \
                      cal_c ** 3 * cal_d ** (-1) * cal_p ** (-2) * cal_r ** 2
        elif regime == 'Gas-FF':
            density = 51 * units.g * units.cm ** (-3) * \
                      viscosity ** (-7 / 10) * mass ** (-5 / 4) * accretion_rate ** (11 / 20) * x ** (-43 / 10) * \
                      cal_c ** (3 / 20) * cal_d ** (-7 / 10) * cal_p ** (11 / 20) * cal_r ** (17 / 40)
        else:
            print('Terminating "densities"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        print('Exiting "densities" from "properties.py" after %.4s ' % (
                time.time() - local_time) + 's \n' if verbose else '', end='')
        return density if density.unit == 'g /cm3' else \
            (print('Terminating "densities"! "density" is in units of "', density.unit,
                   '" instead of "g / cm3"', sep=''), exit())


def surface_densities(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the surface density profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Surface density in gram centimeter^-2.
    """
    local_time = time.time()
    print('Entering "surface_densities" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "surface_densities"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            surface_density = 49610.3 * units.g * units.cm ** (-2) * \
                              viscosity ** (-4 / 5) * mass ** (-2 / 5) * accretion_rate ** (3 / 5) * x ** (-9 / 5) * \
                              cal_c ** (1 / 5) * cal_d ** (-4 / 5) * cal_p ** (3 / 5)

        elif regime == 'Rad-ES':
            surface_density = 10.4248 * units.g * units.cm ** (-2) * \
                              viscosity ** (-1) * mass * accretion_rate ** (-1) * x ** 4 * \
                              cal_c ** 2 * cal_d ** (-1) * cal_p ** (-1) * cal_r

        elif regime == 'Gas-FF':
            surface_density = 170462 * units.g * units.cm ** (-2) * \
                              viscosity ** (-4 / 5) * mass ** (-1 / 2) * accretion_rate ** (7 / 10) * x ** (-11 / 5) * \
                              cal_c ** (1 / 10) * cal_d ** (-4 / 5) * cal_p ** (7 / 10) * cal_r ** (-1 / 20)

        elif regime == 'intra-ISCO':

            # Convert the dimensionless radius to Boyer-Lindquist. #
            r_grav = radii.r_grav(3 * mass * m_sun_cgs).value
            r_isco = radii.r_isco(3 * mass * m_sun_cgs, spin, rotation, vectorized).value
            r = x ** 2 * r_grav

            # Dimensionless radius of the ISCO. #
            x_isco = (r_isco / r_grav) ** (1 / 2)

            # Local path to load the 'outermost_regime' from. #
            ranges_data_path = global_data_path + 'ranges_of_validity/'
            isco_regime = np.load(ranges_data_path + 'isco_regime.npy')

            if isco_regime == 'Gas-ES':
                surface_density_isco = 49610.3 * units.g * units.cm ** (-2) * viscosity ** (-4 / 5) * mass ** (-2 / 5) \
                                       * accretion_rate ** (3 / 5) * x_isco ** (-9 / 5) * cal_c_isco ** (1 / 5) \
                                       * cal_d_isco ** (-4 / 5) * cal_p_isco ** (3 / 5)
                radial_velocity_isco = np.abs(-725088) * units.cm / units.s * viscosity ** (4 / 5) * mass ** (-3 / 5) \
                                       * accretion_rate ** (2 / 5) * x_isco ** (-1 / 5) * cal_c_isco ** (-1 / 5) \
                                       * cal_d_isco ** (4 / 5) * cal_p_isco ** (-3 / 5)

            elif isco_regime == 'Rad-ES':
                surface_density_isco = 10.4248 * units.g * units.cm ** (-2) * viscosity ** (-1) * mass \
                                       * accretion_rate ** (-1) * x_isco ** 4 * cal_c_isco ** 2 \
                                       * cal_d_isco ** (-1) * cal_p_isco ** (-1) * cal_r_isco
                radial_velocity_isco = np.abs(-3.45059e9) * units.cm / units.s * viscosity * mass ** (-2) \
                                       * accretion_rate ** 2 * x_isco ** (-6) * cal_c_isco ** (-2) \
                                       * cal_d_isco * cal_p_isco * cal_r_isco ** (-1)

            elif isco_regime == 'Gas-FF':
                surface_density_isco = 170462 * units.g * units.cm ** (-2) * viscosity ** (-4 / 5) \
                                       * mass ** (-1 / 2) * accretion_rate ** (7 / 10) * x_isco ** (-11 / 5) \
                                       * cal_c_isco ** (1 / 10) * cal_d_isco ** (-4 / 5) * cal_p_isco ** (7 / 10) \
                                       * cal_r_isco ** (-1 / 20)
                radial_velocity_isco = np.abs(-211025) * units.cm / units.s * viscosity ** (4 / 5) * mass ** (-1 / 2) \
                                       * accretion_rate ** (3 / 10) * x_isco ** (1 / 5) * cal_c_isco ** (-1 / 10) \
                                       * cal_d_isco ** (4 / 5) * cal_p_isco ** (-7 / 10) * cal_r_isco ** (1 / 20)
            else:
                print('Terminating "correction_functions"! "isco_regime" is "', regime,
                      '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

            epsilon = radial_velocity_isco / c_cgs * np.sqrt(3 * r_isco / (2 * r_grav))
            surface_density = surface_density_isco * r_isco / r \
                              * (epsilon ** (-1) * (r_isco / r - 1) ** (3 / 2) + 1) ** (-1)

        else:
            print('Terminating "surface_densities"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES", "Gas-FF", or "intra-ISCO"', sep=''), exit()

        print('Exiting "surface_densities" from "properties.py" after %.4s ' % (
                time.time() - local_time) + 's \n' if verbose else '', end='')
        return surface_density if surface_density.unit == 'g / cm2' else (
            print('Terminating "surface_densities"! "surface_density" is in units of "', surface_density.unit,
                  '" instead of "g / cm2"', sep=''), exit())


def opening_angles(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the opening angle profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: Opening angle = height / radial coordinate.
    """
    local_time = time.time()
    print('Entering "opening_angles" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "opening_angles"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            opening_angle = 7e-3 \
                            * viscosity ** (-1 / 10) * mass ** (-3 / 10) * accretion_rate ** (1 / 5) * x ** (-1 / 10) * \
                            cal_c ** (-1 / 10) * cal_d ** (-1 / 10) * cal_p ** (1 / 5) * cal_r ** (-1 / 2)
        elif regime == 'Rad-ES':
            opening_angle = 0.5 * mass ** (-1) * accretion_rate * x ** (-3) * cal_c ** (-1) * cal_p * cal_r ** (-1)
        elif regime == 'Gas-FF':
            opening_angle = 3.8e-3 \
                            * viscosity ** (-1 / 10) * mass ** (-1 / 4) * accretion_rate ** (3 / 20) * x ** (1 / 10) * \
                            cal_c ** (-1 / 20) * cal_d ** (-1 / 10) * cal_p ** (3 / 20) * cal_r ** (-19 / 40)
        else:
            print('Terminating "opening_angles"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        print('Exiting "opening_angles" from "properties.py" after %.4s ' % (time.time() - local_time) + 's \n'
              if verbose else '', end='')
        return opening_angle if opening_angle.unit == '' else \
            (print('Terminating "opening_angles"! "opening_angle" is in units of "', opening_angle.unit,
                   '" instead of dimensionless', sep=''), exit())


def pressures(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the pressure profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the duration of the function.
    :return: Pressure in dyn centimeter^-2.
    """
    local_time = time.time()
    print('Entering "pressures" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "pressures"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            pressure = 1.8e17 * units.dyn * units.cm ** (-2) * \
                       viscosity ** (-9 / 10) * mass ** (-17 / 10) * accretion_rate ** (4 / 5) * x ** (-59 / 10) * \
                       cal_c ** (1 / 10) * cal_d ** (-9 / 10) * cal_p ** (4 / 5) * cal_r ** (1 / 2)
        elif regime == 'Rad-ES':
            pressure = 2.6e15 * units.dyn * units.cm ** (-2) * viscosity ** (-1) * mass ** (-1) * x ** (-3) * \
                       cal_c * cal_d ** (-1) * cal_r
        elif regime == 'Gas-FF':
            pressure = 3.3e17 * units.dyn * units.cm ** (-2) * \
                       viscosity ** (-9 / 10) * mass ** (-7 / 4) * accretion_rate ** (17 / 20) * x ** (-61 / 10) * \
                       cal_c ** (1 / 20) * cal_d ** (-9 / 10) * cal_p ** (17 / 20) * cal_r ** (19 / 40)
        else:
            print('Terminating "pressures"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        print('Exiting "pressures" from "properties.py" after %.4s ' % (
                time.time() - local_time) + 's \n' if verbose else '', end='')
        return pressure if pressure.unit == 'dyn / cm2' else \
            (print('Terminating "pressures"! "pressure" is in units of "', pressure.unit,
                   '" instead of "dyn / cm2"', sep=''), exit())


def temperatures(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the temperature profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the duration of the function.
    :return: Temperature in K.
    """
    local_time = time.time()
    print('Entering "temperatures" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "temperatures"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        if regime == 'Gas-ES':
            temperature = 2.6e8 * units.K * \
                          viscosity ** (-1 / 5) * mass ** (-3 / 5) * accretion_rate ** (2 / 5) * x ** (-11 / 5) * \
                          cal_c ** (-1 / 5) * cal_d ** (-1 / 5) * cal_p ** (2 / 5)
        elif regime == 'Rad-ES':
            temperature = 3.2e7 * units.K * viscosity ** (-1 / 4) * mass ** (-1 / 4) * x ** (-3 / 4) * \
                          cal_c ** (1 / 4) * cal_d ** (-1 / 4) * cal_r ** (1 / 4)
        elif regime == 'Gas-FF':
            temperature = 7.7e7 * units.K * \
                          viscosity ** (-1 / 5) * mass ** (-1 / 2) * accretion_rate ** (3 / 10) * x ** (-9 / 5) * \
                          cal_c ** (-1 / 10) * cal_d ** (-1 / 5) * cal_p ** (3 / 10) * cal_r ** (1 / 20)
        else:
            print('Terminating "temperatures"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        print('Exiting "temperatures" from "temperature_properties.py" after %.4s ' % (
                time.time() - local_time) + 's \n' if verbose else '', end='')
        return temperature if temperature.unit == 'K' else \
            (print('Terminating "temperatures"! "temperature" is in units of "', temperature.unit,
                   '" instead of "K"', sep=''), exit())


def feedback_efficiencies(viscosity, mass, accretion_rate, mass_cgs, accretion_rate_cgs, spin, rotation, regime,
                          vectorized=False, verbose=False):
    """
    Calculate the radiative feedback efficiencies based on Eq. .
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param mass_cgs: Black hole mass in gram.
    :param accretion_rate_cgs: Black hole accretion rate in gram second^-1.
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the duration of the function.
    :return: Eddington accretion rate in gram second^-1.
    """
    local_time = time.time()
    print('Entering "feedback_efficiencies" from "properties.py" \n' if verbose else '', end='')

    if mass_cgs.unit == 'g':  # Check if 'mass_cgs' is in gram.
        # 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. #
        if rotation == 'co-rotation':
            alignment = +1
        elif rotation == 'counter-rotation':
            alignment = -1
        else:
            print('Terminating "feedback_efficiencies"! "rotation" is "', rotation,
                  '" instead of "co-rotation" or "counter-rotation"', sep=''), exit()

        # Spin-dependent radiative efficiency based on Eq. . #
        e_isco = (1 - 2 * radii.r_grav(mass_cgs) / (
                3 * radii.r_isco(mass_cgs, spin, rotation, vectorized))) ** (1 / 2)
        epsilon_alpha_bullet = 1 - e_isco

        # Calculate the correction functions and their value at the ISCO (_isco). #
        x_isco = (radii.r_isco(mass_cgs, spin, rotation, vectorized) / radii.r_grav(mass_cgs)) ** (1 / 2)
        cal_c_isco = 1 - 3 * x_isco ** (-2) + np.sign(alignment) * 2 * spin * x_isco ** (-3)
        cal_d_isco = 1 - 2 * x_isco ** (-2) + spin ** 2 * x_isco ** (-4)
        cal_f_isco = np.sign(alignment) * (
                1 - np.sign(alignment) * 2 * spin * x_isco ** (-3) + spin ** 2 * x_isco ** (-4))
        cal_g_isco = 1 - 2 * x_isco ** (-2) + np.sign(alignment) * spin * x_isco ** (- 3)
        cal_r_isco = cal_c_isco ** (-1) * cal_f_isco ** 2 - spin ** 2 * x_isco ** (-2) * (
                cal_c_isco ** (-1 / 2) * cal_g_isco - 1)

        if regime == 'Gas-ES':
            h_0 = 0.00184101952815808 * \
                  viscosity ** (1 / 8) * mass ** (-3 / 8) * accretion_rate ** (1 / 4) * x_isco ** (1 / 8) * \
                  cal_c_isco ** (-1 / 8) * cal_r_isco ** (-1 / 2)
        elif regime == 'Rad-ES':
            h_0 = 0.002
        elif regime == 'Gas-FF':
            h_0 = 0.0013166295323149043 * \
                  viscosity ** (1 / 17) * mass ** (-5 / 17) * accretion_rate ** (3 / 17) * x_isco ** (5 / 17) * \
                  cal_c_isco ** (-1 / 17) * cal_d_isco ** (-1 / 34) * cal_r_isco ** (-8 / 17)
        else:
            print('Terminating "feedback_efficiencies"! "regime" is "', regime,
                  '" instead of "Gas-ES","Rad-ES" or "Gas-FF"', sep=''), exit()

        cal_p_isco = 2 ** (-1 / 2) * viscosity * x_isco * h_0 * cal_d_isco ** (1 / 2) * cal_r_isco ** (1 / 2)

        # Relativistic angular velocity at the ISCO based on Eq. . #
        omega_isco = np.sign(alignment) * (G_cgs * mass_cgs) ** (1 / 2) / \
                     (radii.r_isco(mass_cgs, spin, rotation, vectorized) ** (3 / 2) +
                      np.sign(alignment) * spin * radii.r_grav(mass_cgs) ** (3 / 2))

        # Relativistic specific angular momentum at the ISCO based on Eq. . #
        z_1 = 1 + (1 - spin ** 2) ** (1 / 3) * ((1 + spin) ** (1 / 3) + (1 - spin) ** (1 / 3))
        z_2 = (3 * spin ** 2 + z_1 ** 2) ** (1 / 2)
        z = 3 + z_2 - np.sign(alignment) * ((3 - z_1) * (3 + z_1 + 2 * z_2)) ** (1 / 2)
        j_isco = np.sign(alignment) * G_cgs * mass_cgs / c_cgs * \
                 (z ** 2 - np.sign(alignment) * 2 * spin * z ** (1 / 2) + spin ** 2) / \
                 (z * (z - 3 + np.sign(alignment) * 2 * spin * z ** (-1 / 2)) ** (1 / 2))

        # Additional torque at the ISCO based on Eq. . #
        g_isco = mass_cgs * accretion_rate_cgs * cal_p_isco * (c_cgs ** 2 * e_isco - omega_isco * j_isco) ** (-1)

        # Total radiative efficiency based on Eq. . #
        epsilon_r = epsilon_alpha_bullet + G_cgs * g_isco * omega_isco / (c_cgs * accretion_rate_cgs)

    else:
        print('Terminating "feedback_efficiencies"! "mass_cgs" is in units of "', mass_cgs.unit, '" instead of "g"',
              sep=''), exit()

    print('Exiting "feedback_efficiencies" from "radii.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
    return epsilon_r, epsilon_alpha_bullet


def fluxes(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized=False, verbose=False):
    """
    Calculate the flux profiles for different regimes based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate based on Eq. .
    :param spin: Dimensionless spin parameter based on Eq. .
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param vectorized: Enable vectorization.
    :param verbose: Print when entering/exiting and the duration of the function.
    :return: Flux in erg centimeter^-2 second^-1.
    """
    local_time = time.time()
    print('Entering "fluxes" from "properties.py" \n' if verbose else '', end='')

    try:  # Check if 'mass' and 'accretion_rate' are dimensionless.
        print('Terminating "fluxes"! "mass" is in units of "', mass.unit,
              '" and "accretion_rate" is in units of "', accretion_rate.unit,
              '" instead of both being dimensionless', sep=''), exit()
    except AttributeError:
        # Get the correction functions. #
        cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, cal_p, \
            cal_p_isco = correction_functions(x, viscosity, mass, accretion_rate, spin, rotation, regime, vectorized)

        fluxes = 5.48793e25 * units.erg * units.cm ** (-2) * units.s ** (-1) * mass ** (-2) * accretion_rate * \
                 x ** (-7) * cal_c ** (-1) * cal_p

        print('Exiting "fluxes" from "fluxes.py" after %.4s ' % (time.time() - local_time) + 's \n'
              if verbose else '', end='')
        return fluxes if fluxes.unit == 'erg / (cm2 s)' else \
            (print('Terminating "fluxes"! "fluxes" is in units of "', fluxes.unit,
                   '" instead of "erg / (cm2 s)"', sep=''), exit())


def toomre_parameter(x, viscosity, bh_mass, accretion_rate, spin, rotation, regime, verbose=False):
    """
    Calculate the local relativistic Toomre Q parameter based on Eq. .
    :param x: Dimensionless radial coordinate.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param regime: Select between 'Gas-ES','Rad-ES' or 'Gas-FF'
    :param verbose: Print when entering/exiting and the local runtime.
    :return: None
    """
    local_time = time.time()
    print('Entering "toomre_parameter" from "properties.py" \n' if verbose else '', end='')

    # 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. #
    if rotation == 'co-rotation':
        alignment = +1
    elif rotation == 'counter-rotation':
        alignment = -1
    else:
        print('Terminating "r_isco_vect"! "rotation" is "', rotation,
              '" instead of "co-rotation" or "counter-rotation"', sep=''), exit()
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    omega = c_cgs ** 3 / ((G_cgs * bh_mass_cgs) * (x ** 3 + alignment * spin))
    kappa = omega * np.sqrt(1 - 6 * x ** (-2) + alignment * 8 * spin * x ** (-3) - 3 * spin ** 2 * x ** (-4))
    sound_speed = np.sqrt(pressures(x, viscosity, bh_mass, accretion_rate, spin, rotation, regime)
                          / densities(x, viscosity, bh_mass, accretion_rate, spin, rotation, regime))

    print('Exiting "toomre_parameter" from "properties.py" after %.4s ' % (time.time() - local_time) + 's \n'
          if verbose else '', end='')
    return kappa * sound_speed / (np.pi * G_cgs
                                  * surface_densities(x, viscosity, bh_mass, accretion_rate, spin, rotation, regime))
