# %%
# Importing packages
import numpy as np
import matplotlib.pyplot as plt

# adding relagn to pythonpath
import os
import sys

mdir = os.path.abspath('.')
mdir = mdir.replace('/examples', '/src/python_version')
sys.path.append(mdir)

from relagn import relagn
import radii, plot_utilities, time, astropy.units as units

# Simulation-specific properties, update as required. #
hubble = 0.7
viscosity = 0.1
m_sun_cgs = 1.99e33 * units.g  # Solar mass in cgs, see Table 1

# Get properties from the relativistic_accretion_logs files. For example, for a random accretion event we have #
rotation = 'co-rotation'
spin_parameter = 5.60032948001510155756e-01  # Dimensionless
BH_Mass = 8.93032240454057365218e-04  # In 1e10 Msun / h.
BH_Mdot = 3.12010333482085414070e-04  # In 1e10 Msun / 9.8e8 yr.
Edd_Mdot = 3.12010333482085448764e-02  # In 1e10 Msun / 9.8e8 yr.

# Convert the above properties from GADGET units to astrophysical and to dimensionless units. #
bh_mass_astro = BH_Mass * 1e10 / hubble  # In Msun.
bh_mass_cgs = bh_mass_astro * m_sun_cgs  # In g.
bh_mass_dmsnless = bh_mass_astro / 3  # Dimensionless black hole mass.

bh_mdot_astro = BH_Mdot * 1e10 / 9.8e8  # In Msun / yr.
bh_mdot_cgs = bh_mdot_astro * m_sun_cgs / (units.yr.to('s') * units.s)  # In g /  s
bh_mdot_dmsnless = bh_mdot_cgs / (1e17 * units.g / units.s)  # Dimensionless black hole accretion rate.

# Create an accretion disc for the BH-AD properties defined above. #
radii.accretion_disc_extent(viscosity, bh_mass_dmsnless, bh_mdot_dmsnless, spin_parameter, rotation)
r_out = radii.r_outermost(bh_mass_cgs) / radii.r_grav(bh_mass_cgs)  # Outermost disc radius in units of r_grav.

# Generate the figure. #
figure, axis = plt.subplots(1, figsize=(10, 10))
plot_utilities.set_axis(axis=axis, x_lim=[1e13, 1e21], y_lim=[3e40, 1e44], x_scale='log', y_scale='log',
                        x_label=r'$\mathrm{Frequency / Hz}$', y_label=r'$\nu F_{\nu} / \mathrm{(erg/s)}$',
                        which='major', size=30)

# Calculate the SED for the updated and the original (method='NT') version of RELAGN.
sagn = relagn(viscosity, rotation, bh_mass_dmsnless, bh_mdot_dmsnless, M=bh_mass_astro, log_rout=np.log10(r_out.value),
              log_mdot=np.log10(BH_Mdot / Edd_Mdot), a=spin_parameter)

nu = sagn.nu_obs  # Frequency grid.
axis.plot(nu, nu * sagn.get_DiscComponent(rel=True), color='teal', lw=3, label='$\mathrm{Disc}$')
axis.loglog(nu, nu * sagn.get_WarmComponent(rel=True), color='plum', lw=3, label='$\mathrm{Warm}$')
axis.loglog(nu, nu * sagn.get_HotComponent(rel=True), color='dodgerblue', lw=3, label='$\mathrm{Hot}$')
axis.loglog(nu, nu * sagn.get_totSED(rel=True), color='k', lw=3, label='$\mathrm{Total}$')

# When using the original RELAGN remove log_rout=np.log10(r_out.value) to truncate the disc at r_sg and set method='NT'.
sagn = relagn(viscosity, rotation, bh_mass_dmsnless, bh_mdot_dmsnless, M=bh_mass_astro, method='NT',
              log_mdot=np.log10(BH_Mdot / Edd_Mdot), a=spin_parameter)

nu = sagn.nu_obs  # Frequency grid.
axis.plot(nu, nu * sagn.get_DiscComponent(rel=True), color='teal', alpha=0.5, ls='dashed', lw=3)
axis.loglog(nu, nu * sagn.get_WarmComponent(rel=True), color='plum', alpha=0.5, ls='dashed', lw=3)
axis.loglog(nu, nu * sagn.get_HotComponent(rel=True), color='dodgerblue', alpha=0.5, ls='dashed', lw=3)
axis.loglog(nu, nu * sagn.get_totSED(rel=True), color='k', alpha=0.5, ls='dashed', lw=3)

# Customise tick labels. #
start, stop = np.ceil(np.log10(axis.get_ylim()[0])), np.floor(np.log10(axis.get_ylim()[1]))
axis.set_yticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
axis.set_yticklabels(labels=['$10^{%.3s}$' % name for name in np.log10(axis.get_yticks()).astype(int)])
start, stop = np.ceil(np.log10(axis.get_xlim()[0])), np.floor(np.log10(axis.get_xlim()[1]))
axis.set_xticks(ticks=np.logspace(start=start, stop=stop, num=int(np.ceil(stop - start)) + 1))
axis.set_xticklabels(labels=['$10^{%.2s}$' % name for name in np.log10(axis.get_xticks()).astype(int)])

# Create the legend, save and close the figure, exit the script. #
axis.legend(ncol=2, loc='upper center', fontsize=25, frameon=False)
plt.savefig('sed' + '_' + time.strftime('%d|%m|%y|%H:%M:%S') + '.pdf', bbox_inches='tight')
plt.close()
