/*
 This file is part of a public black hole accretion model.
 Copyright (c) 2025. Dimitrios Irodotou (di.irodotou@gmail.com) and contributing co-authors.
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// Include existing library files. //
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../proto.h"
#include "../../../allvars.h"


/**
 * @brief Evaluate the local structure of the accretion disc before transferring to global (i.e. P[pt5]) properties.
 * @param pt5 Index of an active black hole particle.
 * @param bh_mass Black hole mass in gram.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param GES_radial_bins Valid 'radial_bins' for the 'GES' regime.
 * @param RES_radial_bins Valid 'radial_bins' for the 'RES' regime.
 * @param GFF_radial_bins Valid 'radial_bins' for the 'GFF' regime.
 * @param in_function Name of function that called global_tests.
 * @param in_script Name of script that contains the function that called global_tests.
 **/
void *local_tests(int pt5, double bh_mass, double spin, const char *rotation, double GES_radial_bins[NUM_BINS],
                  double RES_radial_bins[NUM_BINS], double GFF_radial_bins[NUM_BINS], const char *in_function,
                  const char *in_script) {
    
    // Declare and initialise the 'radial_bins' array with zeros. //
    double radial_bins[NUM_BINS];
    memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost' with 'radial_bin_width' step. //
    (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
    double radial_bin_width = fabs(log10(r_outermost(pt5, bh_mass)) - log10(r_isco(bh_mass, spin, rotation)))
                              / (NUM_BINS - 1);
    
    
    // Sanity checks, more info in the print statements. //
    // Loop over 'NUM_BINS' radial bins. //
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        if (GES_radial_bins[bin] > 0.0 && RES_radial_bins[bin] > 0.0
            && fabs(log10(GES_radial_bins[bin]) - log10(RES_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES_radial_bins' and 'RES_radial_bins' overlap", in_function, in_script);
            
        } else if (GES_radial_bins[bin] > 0.0 && GFF_radial_bins[bin] > 0.0
                   && fabs(log10(GES_radial_bins[bin]) - log10(GFF_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES_radial_bins' and 'GFF_radial_bins' overlap", in_function, in_script);
            
        } else if (RES_radial_bins[bin] > 0.0 && GFF_radial_bins[bin] > 0.0
                   && fabs(log10(RES_radial_bins[bin]) - log10(GFF_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'RES_radial_bins' and 'GFF_radial_bins' overlap", in_function, in_script);
        }
    }
    
    // Sanity checks, more info in the print statements. //
    // Loop over 'NUM_BINS' radial bins excluding the first and the last bin. //
    for (int bin = 1; bin < NUM_BINS - 1; ++bin) {
        if ((GES_radial_bins[bin - 1] > 0.0 && GES_radial_bins[bin - 1] == 0.0 && GES_radial_bins[bin + 1] > 0.0)
            || (GES_radial_bins[bin - 1] == 0.0 && GES_radial_bins[bin - 1] > 0.0 && GES_radial_bins[bin + 1] == 0.0)) {
            terminating_scope("'GES_radial_bins' are not continuous", in_function, in_script);
            
        } else if ((RES_radial_bins[bin - 1] > 0.0 && RES_radial_bins[bin - 1] == 0.0 && RES_radial_bins[bin + 1] > 0.0)
                   || (RES_radial_bins[bin - 1] == 0.0 && RES_radial_bins[bin - 1] > 0.0
                       && RES_radial_bins[bin + 1] == 0.0)) {
            terminating_scope("'RES_radial_bins' are not continuous", in_function, in_script);
            
        } else if ((GFF_radial_bins[bin - 1] > 0.0 && GFF_radial_bins[bin - 1] == 0.0 && GFF_radial_bins[bin + 1] > 0.0)
                   || (GFF_radial_bins[bin - 1] == 0.0 && GFF_radial_bins[bin - 1] > 0.0
                       && GFF_radial_bins[bin + 1] == 0.0)) {
            terminating_scope("'GFF_radial_bins' are not continuous", in_function, in_script);
        }
    }
    
    // Sanity checks, more info in the print statements. //
    double radial_bins_sum = 0.0;
    P[pt5].GES_sum = 0.0, P[pt5].RES_sum = 0.0, P[pt5].GFF_sum = 0.0;
    // Loop over 'NUM_BINS' radial bins. //
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        // Add up all valid 'radial_bins' for each regime. //
        if (GES_radial_bins[bin] > 0.0) { P[pt5].GES_sum += GES_radial_bins[bin]; }
        if (RES_radial_bins[bin] > 0.0) { P[pt5].RES_sum += RES_radial_bins[bin]; }
        if (GFF_radial_bins[bin] > 0.0) { P[pt5].GFF_sum += GFF_radial_bins[bin]; }
        radial_bins_sum += radial_bins[bin]; // Add up all 'radial_bins'.
    }
    
    if (GES_radial_bins[0] == 0.0 && RES_radial_bins[0] == 0.0 && GFF_radial_bins[0] == 0.0) {
        terminating_scope("None of the 'regimes' are valid at the ISCO", in_function, in_script);
        
    } else if (GES_radial_bins[NUM_BINS - 1] == 0.0 && RES_radial_bins[NUM_BINS - 1] == 0.0
               && GFF_radial_bins[NUM_BINS - 1] == 0.0) {
        terminating_scope("None of the 'regimes' are valid at the outermost radius.", in_function, in_script);
        
    } else if (P[pt5].GES_sum == 0.0 && P[pt5].RES_sum == 0.0 && P[pt5].GFF_sum == 0.0) {
        terminating_scope("None of the 'regimes' are valid anywhere", in_function, in_script);
        
    } else if (fabs(log10(P[pt5].GES_sum + P[pt5].RES_sum + P[pt5].GFF_sum) - log10(radial_bins_sum))
               > radial_bin_width / NUM_BINS) {
        terminating_scope("The 'regimes' do not fill all 'radial_bins'", in_function, in_script);
    }
    
} // local_tests(...)


/**
 * @brief Evaluate the global structure of properties (i.e. P[pt5]).
 * @param pt5 Index of an active black hole particle.
 * @param radial_bin_width Logarithmic step of the bins.
 * @param in_function Name of function that called global_tests.
 * @param in_script Name of script that contains the function that called global_tests.
 **/
void global_tests(int pt5, double radial_bin_width, const char *in_function, const char *in_script) {
    
    // Sanity check, more info in the print statement. //
    for (int bin = 0; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_starts[bin] > 0.0 && P[pt5].GES_stops[bin] > 0.0 &&
            fabs(log10(P[pt5].GES_starts[bin]) - log10(P[pt5].GES_stops[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES' starts and stops at the same radius", in_function, in_script);
        }
        if (P[pt5].RES_starts[bin] > 0.0 && P[pt5].RES_stops[bin] > 0.0 &&
            fabs(log10(P[pt5].RES_starts[bin]) - log10(P[pt5].RES_stops[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'RES' starts and stops at the same radius", in_function, in_script);
        }
        if (P[pt5].GFF_starts[bin] > 0.0 && P[pt5].GFF_stops[bin] > 0.0 &&
            fabs(log10(P[pt5].GFF_starts[bin]) - log10(P[pt5].GFF_stops[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GFF' starts and stops at the same radius", in_function, in_script);
        }
        if (P[pt5].GES_starts[bin] > 0.0 && P[pt5].RES_starts[bin] > 0.0 &&
            fabs(log10(P[pt5].GES_starts[bin]) - log10(P[pt5].RES_starts[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES' and 'RES' start at the same radius", in_function, in_script);
        }
        if (P[pt5].GES_starts[bin] > 0.0 && P[pt5].GFF_starts[bin] > 0.0 &&
            fabs(log10(P[pt5].GES_starts[bin]) - log10(P[pt5].GFF_starts[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES' and 'GFF' start at the same radius", in_function, in_script);
        }
        if (P[pt5].GES_stops[bin] > 0.0 && P[pt5].RES_stops[bin] > 0.0 &&
            fabs(log10(P[pt5].GES_stops[bin]) - log10(P[pt5].RES_starts[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES' and 'RES' stop at the same radius", in_function, in_script);
        }
        if (P[pt5].GES_stops[bin] > 0.0 && P[pt5].GFF_stops[bin] > 0.0 &&
            fabs(log10(P[pt5].GES_stops[bin]) - log10(P[pt5].GFF_starts[bin])) < radial_bin_width / NUM_BINS) {
            terminating_scope("'GES' and 'GFF' stop at the same radius", in_function, in_script);
        }
    }
    
    /* Check if every regime does not start and stop equal number of times (also prevents regimes from being valid only
     * in one 'radial_bin'). If yes, terminate. */
    int GES_i = 0, RES_i = 0, GFF_i = 0;
    for (int bin = 0; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_starts[bin] > 0.0) { GES_i++; }
        if (P[pt5].GES_stops[bin] > 0.0) { GES_i--; }
        if (P[pt5].RES_starts[bin] > 0.0) { RES_i++; }
        if (P[pt5].RES_stops[bin] > 0.0) { RES_i--; }
        if (P[pt5].GFF_starts[bin] > 0.0) { GFF_i++; }
        if (P[pt5].GFF_stops[bin] > 0.0) { GFF_i--; }
    }
    if (GES_i != 0) {terminating_scope("'GES' does not start and stop equal number of times", in_function, in_script); }
    if (RES_i != 0) {terminating_scope("'RES' does not start and stop equal number of times", in_function, in_script); }
    if (GFF_i != 0) {terminating_scope("'GFF' does not start and stop equal number of times", in_function, in_script); }
    
} // global_tests(...)


/**
 * @brief Compare the radii returned by C functions with pre-calculated values from independent Python and Mathematica
 * scripts.
 * @param pt5 Index of an active black hole particle.
 **/
void radii_unit_tests(int pt5) {
    
    // Declare local variables. These values should always remain the same. //
    P[pt5].Extent_Factor = 1.0;
    char *rotation = "co-rotation";
    double viscosity = 0.2, mass_dmnless = 1, accretion_rate_dmnless = 1, spin = 0.5;
    double bh_mass = mass_dmnless * (3.0 * 1.99e33), accretion_rate = accretion_rate_dmnless * 1.0e17;
    
    // Check if radii are different from their pre-calculated values +/- some tolerance. If yes, terminate. //
    if (fabs(r_photon(bh_mass, spin, rotation) - 1038545.6237750975) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_photon'");
    }
    if (fabs(r_isco(bh_mass, spin, rotation) - 1872863.7491740503) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_isco'");
    }
    if (fabs(r_grav(bh_mass) - 442443.33333333326) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'r_grav'"); }
    if (fabs(r_outermost(pt5, bh_mass) - 915581489.3011377) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_outermost'");
    }
    
    // Calculate the range of validity of different regimes, which combined form the accretion disc. //
    (void) ranges_of_validity(pt5, viscosity, bh_mass, accretion_rate, spin, rotation);
    
    // Check if ranges are different from their pre-calculated values +/- some tolerance. If yes, terminate. //
    if (fabs(P[pt5].GES_starts[0] - 1872863.7491740503) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'GES_starts[0]'");
    }
    if (fabs(P[pt5].GES_stops[0] - 2767545.6544088884) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GES_stops[0]'"); }
    if (fabs(P[pt5].RES_starts[0] - 2784752.9750331785) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'RES_starts[0]'");
    }
    if (fabs(P[pt5].RES_stops[0] - 11442895.526202) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'RES_stops[0]'"); }
    if (P[pt5].GFF_starts[0] != 0.0) {terminating("Invalid 'GFF_starts[0]'"); }
    if (P[pt5].GFF_stops[0] != 0.0) {terminating("Invalid 'GFF_stops[0]'"); }
    
    if (fabs(P[pt5].GES_starts[1] - 11514042.165418575) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'GES_starts[1]'");
    }
    if (fabs(P[pt5].GES_stops[1] - 915581489.3011377) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GES_stops[1]'"); }
    if (P[pt5].RES_starts[1] != 0.0) {terminating("Invalid 'RES_starts[1]'"); }
    if (P[pt5].RES_stops[1] != 0.0) {terminating("Invalid 'RES_stops[1]'"); }
    if (P[pt5].GFF_starts[1] != 0.0) {terminating("Invalid 'GFF_starts[1]'"); }
    if (P[pt5].GFF_stops[1] != 0.0) {terminating("Invalid 'GFF_stops[1]'"); }
    
    // Loop over the remaining radial bins (i.e. for 'bin = 2') and check if they are not empty. If yes, terminate. //
    for (int bin = 2; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_starts[bin] != 0.0) {terminating("Invalid 'GES_starts[bin]'"); }
        if (P[pt5].GES_stops[bin] != 0.0) {terminating("Invalid 'GES_stops[bin]'"); }
        if (P[pt5].RES_starts[bin] != 0.0) {terminating("Invalid 'RES_starts[bin]'"); }
        if (P[pt5].RES_stops[bin] != 0.0) {terminating("Invalid 'RES_stops[bin]'"); }
        if (P[pt5].GFF_starts[bin] != 0.0) {terminating("Invalid 'GFF_starts[bin]'"); }
        if (P[pt5].GFF_stops[bin] != 0.0) {terminating("Invalid 'GFF_stops[bin]'"); }
    }
    
    // Perform a second test with a different set of parameters. //
    // Declare local variables. These values should always remain the same. //
    P[pt5].Extent_Factor = 10.0;
    rotation = "counter-rotation";
    viscosity = 0.1, mass_dmnless = 1e7, accretion_rate_dmnless = 1e5, spin = 0.9;
    bh_mass = mass_dmnless * (3.0 * 1.99e33), accretion_rate = accretion_rate_dmnless * 1.0e17;
    
    // Check if radii are different from their pre-calculated values +/- some tolerance. If yes, terminate. //
    if (fabs(r_photon(bh_mass, spin, rotation) - 17300719812032.111) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_photon'");
    }
    if (fabs(r_isco(bh_mass, spin, rotation) - 38569344004300.266) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_isco'");
    }
    if (fabs(r_grav(bh_mass) - 4424433333333.333) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'r_grav'"); }
    if (fabs(r_outermost(pt5, bh_mass) - 3644995561168333.0) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'r_outermost'");
    }
    
    // Calculate the range of validity of different regimes, which combined form the accretion disc. //
    (void) ranges_of_validity(pt5, viscosity, bh_mass, accretion_rate, spin, rotation);
    
    // Check if ranges are different from their pre-calculated values +/- some tolerance. If yes, terminate. //
    if (fabs(P[pt5].GES_starts[0] - 45854786460043.1) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GES_starts[0]'"); }
    if (fabs(P[pt5].GES_stops[0] - 561012291213253.3) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GES_stops[0]'"); }
    if (P[pt5].RES_starts[0] != 0.0) {terminating("Invalid 'RES_starts[0]'"); }
    if (P[pt5].RES_stops[0] != 0.0) {terminating("Invalid 'RES_stops[0]'"); }
    if (fabs(P[pt5].GFF_starts[0] - 38569344004300.266) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'GFF_starts[0]'");
    }
    if (fabs(P[pt5].GFF_stops[0] - 45646474780109.47) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GFF_stops[0]'"); }
    
    // Check the next radial bin (i.e. for 'bin = 1'). //
    if (P[pt5].GES_starts[1] != 0.0) {terminating("Invalid 'GES_starts[1]'"); }
    if (P[pt5].GES_stops[1] != 0.0) {terminating("Invalid 'GES_stops[1]'"); }
    if (P[pt5].RES_starts[1] != 0.0) {terminating("Invalid 'RES_starts[1]'"); }
    if (P[pt5].RES_stops[1] != 0.0) {terminating("Invalid 'RES_stops[1]'"); }
    if (fabs(P[pt5].GFF_starts[1] - 563572519870756.9) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GFF_starts[1]'"); }
    if (fabs(P[pt5].GFF_stops[1] - 3644995561168333.0) > UNIT_TEST_TOLERANCE) {terminating("Invalid 'GFF_stops[1]'"); }
    
    // Loop over the remaining radial bins (i.e. for 'bin = 2') and check if they are not empty. If yes, terminate. //
    for (int bin = 2; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_starts[bin] != 0.0) {terminating("Invalid 'GES_starts[bin]'"); }
        if (P[pt5].GES_stops[bin] != 0.0) {terminating("Invalid 'GES_stops[bin]'"); }
        if (P[pt5].RES_starts[bin] != 0.0) {terminating("Invalid 'RES_starts[bin]'"); }
        if (P[pt5].RES_stops[bin] != 0.0) {terminating("Invalid 'RES_stops[bin]'"); }
        if (P[pt5].GFF_starts[bin] != 0.0) {terminating("Invalid 'GFF_starts[bin]'"); }
        if (P[pt5].GFF_stops[bin] != 0.0) {terminating("Invalid 'GFF_stops[bin]'"); }
    }
    
    // Reset the black hole properties. //
    P[pt5].Extent_Factor = 1.0;
    P[pt5].GES_sum = NUM_BINS, P[pt5].RES_sum = NUM_BINS, P[pt5].GFF_sum = NUM_BINS;
    memset(P[pt5].GES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_stops, 0.0, NUM_VALID * sizeof(double));
    
} // radii_unit_tests(...)


/**
 * @brief Compare the structure returned by C functions with pre-calculated values from independent Python and
 * Mathematica scripts.
 * @param pt5 Index of an active black hole particle.
 **/
void accretion_disc_unit_tests(int pt5) {
    
    // Declare local variables. These values should always remain the same. //
    P[pt5].Extent_Factor = 1.0;
    char *rotation = "co-rotation";
    double viscosity = 0.2, mass_dmnless = 1, accretion_rate_dmnless = 1, spin = 0.5;
    double bh_mass = mass_dmnless * (3.0 * 1.99e33), accretion_rate = accretion_rate_dmnless * 1.0e17;
    
    // Calculate the range of validity of different regimes, which combined form the accretion disc. //
    (void) ranges_of_validity(pt5, viscosity, bh_mass, accretion_rate, spin, rotation);
    
    // Declare and initialise the 'radial_bins' arrays with zeros. //
    double radial_bins_GES[NUM_BINS], radial_bins_RES[NUM_BINS], radial_bins_GES_1[NUM_BINS],
            radial_bins_ISCO[NUM_BINS];
    memset(radial_bins_GES, 0.0, NUM_BINS * sizeof(double)), memset(radial_bins_RES, 0.0, NUM_BINS * sizeof(double));
    memset(radial_bins_GES_1, 0.0, NUM_BINS * sizeof(double)), memset(radial_bins_ISCO, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced radial bins from each regime. //
    (void) logspace(log10(P[pt5].GES_starts[0]), log10(P[pt5].GES_stops[0]), NUM_BINS, radial_bins_GES);
    (void) logspace(log10(P[pt5].RES_starts[0]), log10(P[pt5].RES_stops[0]), NUM_BINS, radial_bins_RES);
    (void) logspace(log10(P[pt5].GES_starts[1]), log10(P[pt5].GES_stops[1]), NUM_BINS, radial_bins_GES_1);
    (void) logspace(log10(r_photon(bh_mass, spin, rotation)), log10(r_isco(bh_mass, spin, rotation)), NUM_BINS,
                    radial_bins_ISCO);
    
    // Calculate the surface density of each valid regime. //
    double total_surface_density_GES = 0.0, total_surface_density_RES = 0.0, total_surface_density_ISCO = 0.0;
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        total_surface_density_GES += surface_densities(pt5, radial_bins_GES[bin], viscosity, bh_mass, accretion_rate,
                                                       spin, rotation, "GES")
                                     + surface_densities(pt5, radial_bins_GES_1[bin], viscosity, bh_mass,
                                                         accretion_rate, spin, rotation, "GES");
        total_surface_density_RES += surface_densities(pt5, radial_bins_RES[bin], viscosity, bh_mass, accretion_rate,
                                                       spin, rotation, "RES");
    }
    
    // Loop over 'NUM_BINS' radial bins excluding the first (i.e. 'r_photon') and the last (i.e. 'r_isco') bins. //
    for (int bin = 1; bin < NUM_BINS - 1; ++bin) {
        total_surface_density_ISCO += surface_densities(pt5, radial_bins_ISCO[bin], viscosity, bh_mass, accretion_rate,
                                                        spin, rotation, "ISCO");
    }
    
    /* Check if surface densities are different from their pre-calculated values +/- some tolerance. If yes, terminate.
     */
    if (fabs(total_surface_density_GES - 13641673.27379524) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_GES'");
    }
    if (fabs(total_surface_density_RES - 10415344.435010917) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_GFF'");
    }
    if (fabs(total_surface_density_ISCO - 24712.258310225498) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_ISCO'");
    }
    
    // Perform a second test with a different set of parameters. //
    // Declare local variables. These values should always remain the same. //
    P[pt5].Extent_Factor = 10.0;
    rotation = "counter-rotation";
    viscosity = 0.1, mass_dmnless = 1e7, accretion_rate_dmnless = 1e5, spin = 0.9;
    bh_mass = mass_dmnless * (3.0 * 1.99e33), accretion_rate = accretion_rate_dmnless * 1.0e17;
    
    // Calculate the range of validity of different regimes, which combined form the accretion disc. //
    (void) ranges_of_validity(pt5, viscosity, bh_mass, accretion_rate, spin, rotation);
    
    // Reset the 'radial_bins' arrays to zero. //
    double radial_bins_GFF[NUM_BINS], radial_bins_GFF_1[NUM_BINS];
    memset(radial_bins_GES, 0.0, NUM_BINS * sizeof(double)), memset(radial_bins_ISCO, 0.0, NUM_BINS * sizeof(double));
    memset(radial_bins_GFF, 0.0, NUM_BINS * sizeof(double)), memset(radial_bins_GFF_1, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced radial bins from each regime. //
    (void) logspace(log10(P[pt5].GES_starts[0]), log10(P[pt5].GES_stops[0]), NUM_BINS, radial_bins_GES);
    (void) logspace(log10(P[pt5].GFF_starts[0]), log10(P[pt5].GFF_stops[0]), NUM_BINS, radial_bins_GFF);
    (void) logspace(log10(P[pt5].GFF_starts[1]), log10(P[pt5].GFF_stops[1]), NUM_BINS, radial_bins_GFF_1);
    (void) logspace(log10(r_photon(bh_mass, spin, rotation)), log10(r_isco(bh_mass, spin, rotation)), NUM_BINS,
                    radial_bins_ISCO);
    
    // Calculate the surface density of each valid regime. //
    double total_surface_density_GFF = 0.0;
    total_surface_density_GES = 0.0, total_surface_density_ISCO = 0.0;
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        total_surface_density_GES += surface_densities(pt5, radial_bins_GES[bin], viscosity, bh_mass, accretion_rate,
                                                       spin, rotation, "GES");
        total_surface_density_GFF += surface_densities(pt5, radial_bins_GFF[bin], viscosity, bh_mass, accretion_rate,
                                                       spin, rotation, "GFF")
                                     + surface_densities(pt5, radial_bins_GFF_1[bin], viscosity, bh_mass,
                                                         accretion_rate, spin, rotation, "GFF");
    }
    
    // Loop over 'NUM_BINS' radial bins excluding the first (i.e. 'r_photon') and the last (i.e. 'r_isco') bins. //
    for (int bin = 1; bin < NUM_BINS - 1; ++bin) {
        total_surface_density_ISCO += surface_densities(pt5, radial_bins_ISCO[bin], viscosity, bh_mass, accretion_rate,
                                                        spin, rotation, "ISCO");
    }
    
    /* Check if surface densities are different from their pre-calculated values +/- some tolerance. If yes, terminate.
     */
    if (fabs(total_surface_density_GES - 26104440.18451707) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_GES'");
    }
    if (fabs(total_surface_density_GFF - 16372132.648912407) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_GFF'");
    }
    if (fabs(total_surface_density_ISCO - 300.3822907121095) > UNIT_TEST_TOLERANCE) {
        terminating("Invalid 'total_surface_density_ISCO'");
    }
    
    // Reset the black hole properties. //
    P[pt5].Extent_Factor = 1.0;
    P[pt5].GES_sum = NUM_BINS, P[pt5].RES_sum = NUM_BINS, P[pt5].GFF_sum = NUM_BINS;
    memset(P[pt5].GES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_stops, 0.0, NUM_VALID * sizeof(double));
    
} // accretion_disc_unit_tests(...)