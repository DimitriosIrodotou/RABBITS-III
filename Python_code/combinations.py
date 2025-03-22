# Import global modules. #
from matplotlib.lines import Line2D as Line2D
import matplotlib.pyplot as plt, numpy as np, os, time, inspect
# Import local modules. #
import plot_utilities, properties, radii
from parameters_and_constants import m_sun_cgs

# Extract the global path where the scripts are stored and create global paths to save data and plots. #
global_scripts_path = os.path.dirname(os.path.realpath(__file__))
global_data_path = global_scripts_path[0:-7] + 'data/'
global_plot_path = global_scripts_path[0:-7] + 'plots/'


def surface_densities_viscosities(viscosities, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosities: Dimensionless viscosity parameters.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-3, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\Sigma_\circledcirc / \mathrm{(g\; cm^{-2})}$',
                            which='major', size=30)

    # Add text. #
    figure.text(x=0.5, y=0.1,
                s=r'$\alpha = \{%s, %s, %s\}$' % (str(viscosities[0]), str(viscosities[1]), str(viscosities[2]))
                  + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                  + '\n' r'$m_\bullet = 10^{%s}$' % str(int(np.log10(bh_mass)))
                  + '\n' r'$\dot{m}_\bullet = 10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$0.01$', xy=(1e15, 1e2), xytext=(1e17, 3e3), fontsize=25)
    axis.annotate(r'$0.1$', xy=(1e15, 1e2), xytext=(4e17, 2e2), fontsize=25)
    axis.annotate(r'$1$', xy=(1e15, 1e2), xytext=(2e18, 1e1), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for viscosity in viscosities:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF', 'intra-ISCO'],
                                  ['tab:orange', 'tab:purple', 'tab:brown', 'tab:red']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i],
                          properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                       regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    print(local_plot_path + function_name)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def surface_densities_bh_masses(viscosity, bh_masses, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_masses: Dimensionless black hole masses.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-3, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\Sigma_\circledcirc / \mathrm{(g\; cm^{-2})}$',
                            which='major', size=30)

    # Add text. #
    figure.text(x=0.5, y=0.1, s=r'$\alpha = %s$' % str(viscosity) + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                                + '\n' r'$m_\bullet = \{10^{%s}, 10^{%s}, 10^{%s}\}$'
                                % (str(int(np.log10(bh_masses[0]))), str(int(np.log10(bh_masses[1]))),
                                   str(int(np.log10(bh_masses[2]))))
                                + '\n' r'$\dot{m}_\bullet = 10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$10^6$', xy=(1e15, 1e2), xytext=(1e17, 5e1), fontsize=25)
    axis.annotate(r'$10^7$', xy=(1e15, 1e2), xytext=(4e17, 7e1), fontsize=25)
    axis.annotate(r'$10^8$', xy=(1e15, 1e2), xytext=(1e18, 2e2), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend. #

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for bh_mass in bh_masses:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF', 'intra-ISCO'],
                                  ['tab:orange', 'tab:purple', 'tab:brown', 'tab:red']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i],
                          properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                       regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def surface_densities_accretion_rates(viscosity, bh_mass, accretion_rates, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rates: Dimensionless black hole accretion rates.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-3, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\Sigma_\circledcirc / \mathrm{(g\; cm^{-2})}$',
                            which='major', size=30)

    # Add text. #
    figure.text(x=0.5, y=0.1, s=r'$\alpha = %s$' % str(viscosity) + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                                + '\n' r'$m_\bullet = 10^{%s}$' % str(int(np.log10(bh_mass)))
                                + '\n' r'$\dot{m}_\bullet = \{10^{%s}, 10^{%s}, 10^{%s}\}$'
                                % (str(int(np.log10(accretion_rates[0]))), str(int(np.log10(accretion_rates[1]))),
                                   str(int(np.log10(accretion_rates[2])))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$10^4$', xy=(1e15, 1e2), xytext=(2e18, 1e1), fontsize=25)
    axis.annotate(r'$10^4$', xy=(1e15, 1e2), xytext=(3e12, 1e-3), fontsize=25)
    axis.annotate(r'$10^5$', xy=(1e15, 1e2), xytext=(5e17, 1e2), fontsize=25)
    axis.annotate(r'$10^5$', xy=(1e15, 1e2), xytext=(3e12, 1e-2), fontsize=25)
    axis.annotate(r'$10^6$', xy=(1e15, 1e2), xytext=(2e17, 1e3), fontsize=25)
    axis.annotate(r'$10^6$', xy=(1e15, 1e2), xytext=(3e12, 1e-1), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for accretion_rate in accretion_rates:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF', 'intra-ISCO'],
                                  ['tab:orange', 'tab:purple', 'tab:brown', 'tab:red']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i],
                          properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                       regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def surface_densities_spins(viscosity, bh_mass, accretion_rate, spins, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spins: Dimensionless spin parameters.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-3, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\Sigma_\circledcirc / \mathrm{(g\; cm^{-2})}$',
                            which='major', size=30)

    # Add text. #
    figure.text(x=0.5, y=0.1, s=r'$\alpha=%s$' % str(viscosity)
                                + '\n' r'$\alpha_\bullet = \{ %s, %s, %s \}$'
                                % (str(spins[0]), str(spins[1]), str(spins[2]))
                                + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                                + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$0.0$', xy=(1e15, 1e2), xytext=(2e13, 3e-3), fontsize=25)
    axis.annotate(r'$0.5$', xy=(1e15, 1e2), xytext=(5e12, 4e-3), fontsize=25)
    axis.annotate(r'$0.9$', xy=(1e15, 1e2), xytext=(2e12, 1e-2), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend. #

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for spin in spins:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF', 'intra-ISCO'],
                                  ['tab:orange', 'tab:purple', 'tab:brown', 'tab:red']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i],
                          properties.surface_densities(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                       regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def stabilities_viscosities(viscosities, bh_mass, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the stability profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosities: Dimensionless viscosity parameters.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-1, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$',
                            y_label='$c_\mathrm{s}\; \kappa\; / (\pi\ G\ \Sigma_\circledcirc)$', which='major', size=30)

    # Add text. #
    figure.text(x=0.1, y=0.05,
                s=r'$\alpha = \{%s, %s, %s\}$' % (str(viscosities[0]), str(viscosities[1]), str(viscosities[2]))
                  + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                  + '\n' r'$m_\bullet = 10^{%s}$' % str(int(np.log10(bh_mass)))
                  + '\n' r'$\dot{m}_\bullet = 10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$0.01$', xy=(1e15, 1e2), xytext=(4e16, 4e-1), fontsize=25)
    axis.annotate(r'$0.1$', xy=(1e15, 1e2), xytext=(2e17, 4e-1), fontsize=25)
    axis.annotate(r'$1$', xy=(1e15, 1e2), xytext=(1e18, 4e-1), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for viscosity in viscosities:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i][1:] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i][1:],
                          properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                      regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def stabilities_bh_masses(viscosity, bh_masses, accretion_rate, spin, rotation, verbose=False):
    """
    Plot the stability profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_masses: Dimensionless black hole masses.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-1, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$',
                            y_label='$c_\mathrm{s}\; \kappa\; / (\pi\ G\ \Sigma_\circledcirc)$', which='major', size=30)

    # Add text. #
    figure.text(x=0.1, y=0.05, s=r'$\alpha = %s$' % str(viscosity) + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                                 + '\n' r'$m_\bullet = \{10^{%s}, 10^{%s}, 10^{%s}\}$'
                                 % (str(int(np.log10(bh_masses[0]))), str(int(np.log10(bh_masses[1]))),
                                    str(int(np.log10(bh_masses[2]))))
                                 + '\n' r'$\dot{m}_\bullet = 10^{%s}$' % str(int(np.log10(accretion_rate))),
                fontsize=30, transform=axis.transAxes)

    axis.annotate(r'$10^6$', xy=(1e15, 1e2), xytext=(2e12, 3e6), fontsize=25)
    axis.annotate(r'$10^7$', xy=(1e15, 1e2), xytext=(2e13, 6e5), fontsize=25)
    axis.annotate(r'$10^8$', xy=(1e15, 1e2), xytext=(2e14, 1.5e5), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for bh_mass in bh_masses:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i][1:] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i][1:],
                          properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                      regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def stabilities_accretion_rates(viscosity, bh_mass, accretion_rates, spin, rotation, verbose=False):
    """
    Plot the surface density profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rates: Dimensionless black hole accretion rates.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-1, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$',
                            y_label='$c_\mathrm{s}\; \kappa\; / (\pi\ G\ \Sigma_\circledcirc)$', which='major', size=30)

    # Add text. #
    figure.text(x=0.1, y=0.05, s=r'$\alpha = %s$' % str(viscosity) + '\n' r'$\alpha_\bullet = %s$' % str(spin)
                                 + '\n' r'$m_\bullet = 10^{%s}$' % str(int(np.log10(bh_mass)))
                                 + '\n' r'$\dot{m}_\bullet = \{10^{%s}, 10^{%s}, 10^{%s}\}$'
                                 % (str(int(np.log10(accretion_rates[0]))), str(int(np.log10(accretion_rates[1]))),
                                    str(int(np.log10(accretion_rates[2])))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$10^4$', xy=(1e15, 1e2), xytext=(1e18, 3e-1), fontsize=25)
    axis.annotate(r'$10^5$', xy=(1e15, 1e2), xytext=(2.5e17, 3e-1), fontsize=25)
    axis.annotate(r'$10^6$', xy=(1e15, 1e2), xytext=(7e16, 3e-1), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for accretion_rate in accretion_rates:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i][1:] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i][1:],
                          properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                      regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def stabilities_spins(viscosity, bh_mass, accretion_rate, spins, rotation, verbose=False):
    """
    Plot the stability profile for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spins: Dimensionless spin parameters.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[1e-1, 1e7], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$',
                            y_label='$c_\mathrm{s}\; \kappa\; / (\pi\ G\ \Sigma_\circledcirc)$', which='major', size=30)

    # Add text. #
    figure.text(x=0.1, y=0.05, s=r'$\alpha=%s$' % str(viscosity)
                                 + '\n' r'$\alpha_\bullet = \{ %s, %s, %s \}$'
                                 % (str(spins[0]), str(spins[1]), str(spins[2]))
                                 + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                                 + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))), fontsize=30,
                transform=axis.transAxes)

    axis.annotate(r'$0.0$', xy=(1e15, 1e2), xytext=(3e13, 1e6), fontsize=25)
    axis.annotate(r'$0.5$', xy=(1e15, 1e2), xytext=(1e13, 2e6), fontsize=25)
    axis.annotate(r'$0.9$', xy=(1e15, 1e2), xytext=(3e12, 1e6), fontsize=25)
    lines, labels, regimes = [], [], []  # Declare empty lists to store information for the plot's legend.

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    for spin in spins:
        bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

        # Build an accretion disc and store the ranges of validity of each regime. #
        radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

        # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
        ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()
        for regime, colour in zip(['Gas-ES', 'Rad-ES', 'Gas-FF'], ['tab:orange', 'tab:purple', 'tab:brown']):

            for i in range(len(ranges_of_validity[regime])):
                x_bins = (ranges_of_validity[regime][i][1:] / radii.r_grav(bh_mass_cgs).value) \
                         ** (1 / 2)  # Dimensionless bins.

                axis.plot(ranges_of_validity[regime][i][1:],
                          properties.toomre_parameter(x_bins, viscosity, bh_mass, accretion_rate, spin, rotation,
                                                      regime), color=colour, lw=3)

            # If the 'regime' exists, create a label and a line for it for the legend. #
            if len(ranges_of_validity[regime]) > 0 and regime not in regimes:
                regimes.append(regime)
                labels.append('$\mathrm{%s}$' % regime)
                lines.append(Line2D(xdata=[0], ydata=[0], color=colour, lw=0))

    # Create the legend, save and close the figure, exit the script. #
    axis.legend(handles=lines, labels=labels, loc='upper right', labelcolor='linecolor', fontsize=25, frameon=False)
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def correction_functions(viscosity, bh_mass, accretion_rate, spins, verbose=False):
    """
    Plot the correction function for given black hole and accretion disc properties in parameters_and_constants.py.
    :param viscosity: Dimensionless viscosity parameter.
    :param bh_mass: Dimensionless black hole mass.
    :param accretion_rate: Dimensionless black hole accretion rate.
    :param spins: Dimensionless spin parameters.
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
    plot_utilities.set_axis(axis=axis, x_lim=[1e11, 1e19], y_lim=[2e-2, 2e0], x_scale='log', y_scale='log',
                            x_label='$r / \mathrm{cm}$', y_label='$\mathcal{C}$', which='major', size=30)

    # Add text. #
    figure.text(x=2e15, y=3e-2, fontsize=30, transform=axis.transData,
                s=r'$\alpha=%s$' % str(viscosity) + '\n' r'$m_\bullet= 10^{%s}$' % str(int(np.log10(bh_mass)))
                  + '\n' r'$\dot{m}_\bullet=10^{%s}$' % str(int(np.log10(accretion_rate))))

    # Customise tick labels. #
    start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
    axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
    start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
    axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
    axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

    # Define colour maps for 'co-rotation' and 'counter-rotation'. #
    cmap_co_rotation, cmap_counter_rotation = plt.cm.YlGnBu, plt.cm.YlOrRd
    colours_co_rotation = cmap_co_rotation(np.linspace(start=0, stop=1, num=len(spins)))
    colours_counter_rotation = cmap_counter_rotation(np.linspace(start=0, stop=1, num=len(spins)))

    # Customise the colour bar. #
    cbaxes = axis.inset_axes(bounds=(1e14, 1e-1, 1e16, 0.02), transform=axis.transData)
    plot_utilities.create_colorbar(axis=cbaxes, plot=plt.cm.ScalarMappable(cmap=cmap_counter_rotation.reversed()),
                                   label=r'$-\alpha_\bullet$', size=30, orientation='horizontal', ticks=[0, 0.5, 0.998],
                                   ticklabels=['$0.998$', '$0.5$', ''])
    cbaxes = axis.inset_axes(bounds=(1e16, 1e-1, 01e18, 0.02), transform=axis.transData)
    plot_utilities.create_colorbar(axis=cbaxes, plot=plt.cm.ScalarMappable(cmap=cmap_co_rotation),
                                   label=r'$\alpha_\bullet$', size=30, orientation='horizontal', ticks=[0, 0.5, 0.998],
                                   ticklabels=['$0.0$', '$0.5$', '$0.998$'])

    bh_mass_cgs = 3 * bh_mass * m_sun_cgs  # Black hole mass in gram.

    for colour_index, spin in enumerate(spins):
        print(spin)
        # Loop over the two orientation
        for rotation, coloumap in zip(['co-rotation', 'counter-rotation'],
                                      [colours_co_rotation, colours_counter_rotation]):
            # Build an accretion disc and store the ranges of validity of each regime. #
            radii.accretion_disc_extent(viscosity, bh_mass, accretion_rate, spin, rotation)

            # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
            ranges_of_validity = np.load(ranges_data_path + 'ranges_of_validity.npy', allow_pickle=True).item()

            # For each 'regime' plot the corresponding surface density profile inside its 'ranges_of_validity'. #
            for regime in ['Gas-ES', 'Rad-ES', 'Gas-FF']:
                for i in range(len(ranges_of_validity[regime])):
                    # Dimensionless bins. #
                    x_bins = (ranges_of_validity[regime][i] / radii.r_grav(bh_mass_cgs).value) ** (1 / 2)

                    cal_c, cal_c_isco, cal_d, cal_d_isco, cal_f, cal_f_isco, cal_g, cal_g_isco, cal_r, cal_r_isco, \
                        cal_p, cal_p_isco = properties.correction_functions(x_bins, viscosity, bh_mass, accretion_rate,
                                                                            spin, rotation, regime)

                    axis.plot(ranges_of_validity[regime][i], cal_c, color=coloumap[colour_index], lw=3)

    # Create the legend, save and close the figure, exit the script. #
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')


def parameter_space(viscosities, bh_masses, accretion_rates, spins, verbose=False):
    """
    Plot the parameter space in a radar chart.
    :param viscosities: Dimensionless viscosity parameters.
    :param bh_masses: Dimensionless black hole masses.
    :param accretion_rates: Dimensionless black hole accretion rates.
    :param spins: Dimensionless spin parameters.
    :param verbose: Print when entering/exiting and the local runtime.
    :return:
    """
    local_time = time.time()

    # Extract this function's name and create local paths to load the data and save the plot. #
    function_name = inspect.currentframe().f_code.co_name
    ranges_data_path, local_plot_path = global_data_path + 'ranges_of_validity/', global_plot_path + function_name + '/'
    print('Entering "' + function_name + '()" from "' + os.path.basename(__file__) + '" \n' if verbose else '', end='')

    # Check if the path to save the plot exists. If not, then create the corresponding directory. #
    if not os.path.exists(local_plot_path): os.makedirs(local_plot_path)

    # Generate a figure and set the axis parameters. #
    figure = plt.figure(figsize=(10, 10))

    values = [['$10^{%.2s}$' % bh_mass for bh_mass in np.log10(bh_masses).astype(int)],
              ['$%.3s$' % spin for spin in spins],
              ['$%.4s$' % viscosity for viscosity in viscosities],
              ['$10^{%.2s}$' % accretion_rate for accretion_rate in np.log10(accretion_rates).astype(int)]]
    axis_labels = [r'$m_\bullet$', r'$\alpha_\bullet$', r'$\alpha$', r'$\dot{m}_\bullet$']

    axes = [figure.add_axes([0.05, 0.05, 0.95, 0.95], projection="polar", label="axes%d" % i) for i in
            range(len(axis_labels))]
    angles = np.arange(0, 360, 360.0 / len(axis_labels))
    axes[0].set_thetagrids(angles, labels=axis_labels, fontsize=30)

    # Set the polar and orthogonal axes line widths. #
    x_gridlines, y_gridlines = axes[0].xaxis.get_gridlines(), axes[0].yaxis.get_gridlines()
    for x_gridline, y_gridline in zip(x_gridlines, y_gridlines):
        x_gridline.set_linewidth(2)
        y_gridline.set_linewidth(2)
        x_gridline.set_linestyle('solid')
        y_gridline.set_linestyle('dashed')

    # Make the rest of the axes invisible. #
    for axis in axes[1:]:
        axis.patch.set_visible(False)
        axis.xaxis.set_visible(False)

    # Make the rest of the axes invisible. #
    for axis, angle, value in zip(axes, angles, values):
        axis.set_ylim(0, len(values))
        axis.spines["polar"].set_visible(False)
        axis.set_theta_offset(np.deg2rad(180 / len(axis_labels)))
        axis.set_rgrids(range(1, len(values)), angle=angle, labels=value, fontsize=40)

    # Set the colour of labels and values of a given axis. #
    preferred_axis = 2
    polar_labels, polar_values = axes[0].get_xmajorticklabels(), axes[preferred_axis].get_ymajorticklabels()

    # Axis label colour #
    polar_labels[preferred_axis].set_color('red')

    # Axis values colour. #
    for axis in axes:
        polar_value = axis.get_ymajorticklabels()
        # polar_value[0].set_color('red')
        polar_value[1].set_color('blue')
        # polar_value[2].set_color('red')

    for polar_value in polar_values:
        polar_value.set_color('red')

    # Move the labels radially out. #
    for tick in axes[0].xaxis.get_major_ticks(): tick.set_pad(15)

    # Add text. #
    figure.text(x=0.85, y=0.55, s=r'$Top$', color='red', fontsize=30)
    figure.text(x=0.72, y=0.55, s=r'$Mid.$', color='red', fontsize=30)
    figure.text(x=0.6, y=0.55, s=r'$Bot.$', color='red', fontsize=30)

    # Move the labels radially out. #
    for tick in axes[0].xaxis.get_major_ticks(): tick.set_pad(15)

    # Create the legends. #
    labels = ['$\mathrm{Gas\!-\!ES}$', '$\mathrm{Rad\!-\!ES}$', '$\mathrm{Gas\!-\!FF}$', '$\mathrm{intra\!-\!ISCO}$']
    lines = [Line2D(xdata=[0], ydata=[0], color='tab:orange', lw=0),
             Line2D(xdata=[0], ydata=[0], color='tab:purple', lw=0),
             Line2D(xdata=[0], ydata=[0], color='tab:brown', lw=0),
             Line2D(xdata=[0], ydata=[0], color='tab:red', lw=0)]
    axes[0].legend(handles=lines, labels=labels, loc='upper center', ncols=4, labelcolor='linecolor', fontsize=25,
                   columnspacing=-0.5, frameon=False)

    # Save and close the figure, exit the script. #
    plt.savefig(fname=local_plot_path + function_name + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.png',
                bbox_inches='tight')
    plt.close()
    print('Exiting after %.4s ' % (time.time() - local_time) + 's \n' if verbose else '', end='')
