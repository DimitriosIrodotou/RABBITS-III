# Import global modules. #
from matplotlib.lines import Line2D as Line2D
import matplotlib.pyplot as plt, numpy as np, os, time, inspect
# Import local modules. #
import plot_utilities, properties, radii, utilities, astropy.units as units
from parameters_and_constants import m_sun_cgs, num_bins, atol, rtol

# Extract the global path where the scripts are stored and create global paths to save data and plots. #
global_scripts_path = os.path.dirname(os.path.realpath(__file__))
global_data_path = global_scripts_path[0:-7] + 'data/'
global_plot_path = global_scripts_path[0:-7] + 'plots/'


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


def surface_density(viscosity, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    """
    local_time = time.time()  # Start the time to calculate the local runtime.

    # Extract this function's name and create local paths to load the data and save the plot. #
    function_name = inspect.currentframe().f_code.co_name
    ranges_data_path, local_plot_path = global_data_path + 'ranges_of_validity/', global_plot_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the plot exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_plot_path): os.makedirs(local_plot_path)

    # Generate a figure and set the axis parameters. #
    figure, axis = plt.subplots(1, figsize=(10, 10))
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    plot_utilities.set_axis(axis=axis, x_lim=[1e12, 1e19], y_lim=[1e-3, 1e6], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\Sigma_\circledcirc / \mathrm{(g\; cm^{-2})}$',
                            which='major', size=30)

    # Add text. #
    figure.text(x=0.74, y=0.7, s=r'$\alpha=%s$' % str(viscosity) + '\n' r'$\alpha_\bullet= %s$' % str(spin)
                                 + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                                 + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
    lines, labels = [], []
    ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
    for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF', 'intra-ISCO'],
                              ['tab:orange', 'tab:purple', 'tab:brown', 'tab:red']):

        for i in range(len(ranges_of_validity[regime])):
            x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) ** (1 / 2)  # Dimensionless bins.

            axis.plot(ranges_of_validity[regime][i],
                      properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, regime),
                      color=colour, lw=5)

        # If the 'regime' exists, create a label and a line for it for the legend. #
        if len(ranges_of_validity[regime]) > 0:
            labels.append('$\mathrm{%s}$' % regime)
            lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='lower right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def temperature(viscosity, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    """
    local_time = time.time()  # Start the time to calculate the local runtime.

    # Extract this function's name and create local paths to load the data and save the plot. #
    function_name = inspect.currentframe().f_code.co_name
    ranges_data_path, local_plot_path = global_data_path + 'ranges_of_validity/', global_plot_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the plot exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_plot_path): os.makedirs(local_plot_path)

    # Generate a figure and set the axis parameters. #
    figure, axis = plt.subplots(1, figsize=(10, 10))
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    plot_utilities.set_axis(axis=axis, x_lim=[1e12, 1e19], y_lim=[1e2, 1e5], x_scale='log', y_scale='log',
                            x_label='$r/\mathrm{cm}$', y_label='$T/\mathrm{K}$', which='major', size=30)

    # Add text. #
    figure.text(x=0.74, y=0.7, s=r'$\alpha=%s$' % str(viscosity) + '\n' r'$\alpha_\bullet= %.5s$' % str(spin)
                                 + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                                 + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
    lines, labels = [], []
    ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
    for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'],
                              ['tab:orange', 'tab:purple', 'tab:brown']):

        for i in range(len(ranges_of_validity[regime])):
            x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) ** (1 / 2)  # Dimensionless bins.

            sigma = properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, regime)
            temperature = properties.temperatures(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, regime)

            if regime == 'Gas-ES' or regime == 'Rad-ES':
                kappa = 0.40 * units.cm ** 2 * units.g ** (-1)
            elif regime == 'Gas-FF':
                density = properties.densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, regime)
                kappa = 0.64e23 * density / (units.g / units.cm ** 3) * (temperature / units.K) ** (-7 / 2) \
                        * units.cm ** 2 * units.g ** (-1)

            # Convert the mid-plane temperature to the effective temperature. #
            t_eff = (8 * temperature ** 4 / (3 * kappa * sigma)) ** (1 / 4)

            axis.scatter(ranges_of_validity[regime][i], t_eff, color=colour)

        # If the 'regime' exists, create a label and a line for it for the legend. #
        if len(ranges_of_validity[regime]) > 0:
            labels.append('$\mathrm{%s}$' % regime)
            lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # RELAGN SED example with values corresponding to a snapshot of an isolated galaxy formation simulation.
    # data = np.loadtxt('data_disc.txt', dtype=float)
    # axis.plot(data[:, 0], data[:, 1], lw=4, ls='dashed', color='teal')
    # data = np.loadtxt('data_warm.txt', dtype=float)
    # axis.plot(data[:, 0], data[:, 1], lw=4, ls='dashed', color='plum')
    # data = np.loadtxt('data_hot.txt', dtype=float)
    # axis.plot(data[:, 0], data[:, 1], lw=4, ls='dashed', color='dodgerblue')

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='lower right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.pdf',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def stability(viscosity, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the Toomre Q profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spin: Dimensionless spin parameter.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    """
    local_time = time.time()  # Start the time to calculate the local runtime.

    # Extract this function's name and create local paths to load the data and save the plot. #
    function_name = inspect.currentframe().f_code.co_name
    ranges_data_path, local_plot_path = global_data_path + 'ranges_of_validity/', global_plot_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the plot exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_plot_path): os.makedirs(local_plot_path)

    # Generate a figure and set the axis parameters. #
    figure, axis = plt.subplots(1, figsize=(10, 10))
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-2, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$',
                            y_label='$c_\mathrm{s}\; \kappa\; / (\pi\ G\ \Sigma_\circledcirc)$', which='major', size=30)

    # Add text. #
    figure.text(x=0.74, y=0.7, s=r'$\alpha=%s$' % str(viscosity) + '\n' r'$\alpha_\bullet= %s$' % str(spin)
                                 + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                                 + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
    lines, labels = [], []
    ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
    for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):

        for i in range(len(ranges_of_validity[regime])):
            x_bins = (ranges_of_validity[regime][i][1:] / radii.r_grav(bh_mass_cgs).value) ** (
                    1 / 2)  # Dimensionless bins.

            axis.plot(ranges_of_validity[regime][i][1:],
                      properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation, regime),
                      color=colour, lw=5)

        # If the 'regime' exists, create a label and a line for it for the legend. #
        if len(ranges_of_validity[regime]) > 0:
            labels.append('$\mathrm{%s}$' % regime)
            lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='lower right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def regional_solutions(viscosity, bh_mass, accretion_rate, spins, rotation, verbose=False):
    """
    Plot the regional solutions profile for a given viscosity, mass and accretion rate, and all dimensionless spins. #
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spins: List of dimensionless spin parameters.
    :param rotation: Select between 'co-rotation' or 'counter-rotation'.
    :param verbose: Print when entering/exiting and the local runtime.
    :return: None
    """
    local_time = time.time()  # Start the time to calculate the local runtime.

    # Extract this function's name and create local paths to load the data and save the plot. #
    function_name = inspect.currentframe().f_code.co_name
    ranges_data_path, local_plot_path = global_data_path + 'ranges_of_validity/', global_plot_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the plot exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_plot_path): os.makedirs(local_plot_path)

    # Generate a figure and set the axis parameters. #
    figure, axis = plt.subplots(1, figsize=(10, 10))
    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19],
                            y_lim=[0.0, 1.07], x_label='$r / \mathrm{cm}$', y_label=r'$\alpha_\bullet$',
                            x_scale='log', semilog_x=True, which='major', size=30)

    # Add text. #
    figure.text(x=1e12, y=1.01, s=r'$\alpha=%s$' % str(viscosity), fontsize=30, transform=axis.transData)
    figure.text(x=1e14, y=1.01, s=r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass))), fontsize=30,
                transform=axis.transData)
    figure.text(x=1e16, y=1.01, s=r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transData)

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])
    start, stop = np.ceil(axis.get_ylim()[0]), np.floor(axis.get_ylim()[1])
    axis.set_yticks(ticks=np.linspace(start=start, stop=stop, num=5 * int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$%.3s$' % name for name in axis.get_yticks()])

    # Create the legends. #
    labels = ['$\mathrm{Gas-ES}$', '$\mathrm{Rad-ES}$', '$\mathrm{Gas-FF}$', '$\mathrm{intra-ISCO}$']
    lines = [Line2D(xdata=[0], ydata=[0], color='tab:orange', lw=0),
             Line2D(xdata=[0], ydata=[0], color='tab:purple', lw=0),
             Line2D(xdata=[0], ydata=[0], color='tab:brown', lw=0), Line2D(xdata=[0], ydata=[0], color='tab:red', lw=0)]
    legend = axis.legend(handles=lines, labels=labels, loc='lower left', bbox_to_anchor=(-0.15, 0),
                         labelcolor='linecolor', fontsize=30, frameon=False)

    # Calculate the dimensionless radius from the ISCO to the outermost radius. #
    rs = np.logspace(start=np.log10(radii.r_isco(bh_mass_cgs, spins, rotation, vectorized=True).value),
                     stop=np.log10(radii.r_outermost(bh_mass_cgs).value), num=len(spins))
    x_bins = (rs / radii.r_grav(bh_mass_cgs).value) ** (1 / 2)
    x = x_bins.flatten() ** 2 * radii.r_grav(bh_mass_cgs)
    y = (np.ones(shape=[len(spins), len(spins)]) * spins).flatten()

    # Loop over all regimes. #
    for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):
        # Calculate the range of validity of a 'regime'. #
        mask_regime, = np.where(
            (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spins, rotation, regime, vectorized=True)[
                 0].flatten() < 1)
            &
            (properties.ratios(x_bins, viscosity, bh_mass, accretion_rate, spins, rotation, regime, vectorized=True)[
                 1].flatten() < 1))

        # Check if the 'regime' is anywhere valid. If yes, then plot the corresponding spin-radius plane. #
        if len(mask_regime) > 0:
            axis.scatter(x=x[mask_regime], y=y[mask_regime], color=colour, s=2, edgecolor='none', rasterized=True)

    # Calculate the dimensionless radius for the intra-ISCO region but exclude 'r_photon' and 'r_isco'. #
    intra_radial_bins = np.logspace(start=np.log10(radii.r_photon(bh_mass_cgs, spins, rotation, vectorized=True).value),
                                    stop=np.log10(radii.r_isco(bh_mass_cgs, spins, rotation, vectorized=True).value),
                                    num=len(spins) + 2)
    intra_y = (np.ones(shape=[len(spins), len(spins)]) * spins).flatten()
    axis.scatter(x=intra_radial_bins[1:-1], y=intra_y, color='tab:red', s=2, edgecolor='none', rasterized=True)

    # Plot the dimensionless radius of the photon orbit. #
    # axis.plot(radii.r_photon(bh_mass_cgs, spins, rotation, vectorized=True), spins, linestyle=(0, (5, 10)), lw=5,
    #           color='k', label='$R_\mathrm{ph}$')

    # Plot the dimensionless radius of the ISCO. #
    # axis.plot(radii.r_isco(bh_mass_cgs, spins, rotation, vectorized=True), spins, linestyle=(0, (15, 10)), lw=5,
    #           color='k', label='$R_\mathrm{isco}$')

    # Plot the outermost radius. #
    # axis.vlines(x=radii.r_outermost(bh_mass_cgs).value, ymin=min(spins), ymax=max(spins), ls='-', lw=5, color='k',
    #             label='$R_\mathrm{outer}$')

    # Create the legend, save and close the figure, exit the script. #
    # axis.add_artist(legend)
    # axis.legend(loc='center right', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')
