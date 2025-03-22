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
#include <stdlib.h>
#include <string.h>

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../proto.h"


/**
 * @brief Calculate the radial extent of the accretion disc.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 **/
void accretion_disc_extent(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin,
                           const char *rotation) {
    
    double alignment; // Declare local variable.
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Declare local variables. //
    int iterations = 1;
    const char *outermost_regime;
    double extent = 0.0, innermost_radial_bin = 0.0, innermost_toomre_parameter = 0.0;
    
    /* Increase the 'Extent_Factor' by a factor of 'NUM_STEPS' in every iteration, as long as a self-gravitation radius
     * has be defined. */
    do {
        // Calculate the range of validity of different regimes, which combined form the accretion disc. //
        (void) ranges_of_validity(pt5, viscosity, bh_mass, accretion_rate, spin, rotation);
        
        // Find the maximum starting radius for each regime. //
        double GES_max = maximum_array_value(NUM_VALID, P[pt5].GES_starts);
        double RES_max = maximum_array_value(NUM_VALID, P[pt5].RES_starts);
        double GFF_max = maximum_array_value(NUM_VALID, P[pt5].GFF_starts);
        
        /* Check if 'extent' is positive. If yes, find the regime that is valid at the outskirts of the accretion
        * disc. Else, 'outermost_regime' has been calculated in the previous do-while iteration as the second to last
        * 'outermost_regime' . */
        if (extent >= 0.0) {
            if (RES_max > GES_max && RES_max > GFF_max) { outermost_regime = "RES"; }
            else if (GES_max > RES_max && GES_max > GFF_max) { outermost_regime = "GES"; }
            else if (GFF_max > GES_max && GFF_max > RES_max) { outermost_regime = "GFF"; }
            else {terminating("Invalid 'outermost_regime'"); } // Sanity check, more info in the print statement.
        }
        
        // Declare and initialise the 'radial_bins' array with zeros. //
        double radial_bins[NUM_BINS];
        memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
        
        /* Calculate evenly log-spaced 'radial_bins' from the starting radius of the 'outermost_regime' to 'r_outermost'
         * with 'radial_bin_width' step . */
        if (strcmp(outermost_regime, "RES") == 0) {
            (void) logspace(log10(RES_max), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
        } else if (strcmp(outermost_regime, "GES") == 0) {
            (void) logspace(log10(GES_max), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
        } else if (strcmp(outermost_regime, "GFF") == 0) {
            (void) logspace(log10(GFF_max), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
        }
        
        // Declare and initialise the 'toomre_parameters' array with zeros. //
        double toomre_parameters[NUM_BINS];
        memset(toomre_parameters, 0.0, NUM_BINS * sizeof(double));
        
        /* Check if the 'toomre_parameter' at the outskirts of the accretion disk (i.e. for r=radial_bins[NUM_BINS])
         * is less than one. If yes, then calculate the 'toomre_parameters' for all 'radial_bins' because instability
         * will occur in these 'radial_bins'. */
        double outermost_toomre_parameter = toomre_parameter(pt5, radial_bins[NUM_BINS - 1], viscosity, bh_mass,
                                                             accretion_rate, spin, rotation, outermost_regime);
        
        // Sanity check, more info in the print statement. //
        // Loop over 'NUM_BINS' radial bins excluding the first bin. //
        for (int bin = 0; bin < NUM_BINS; ++bin) {
            if (toomre_parameter(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                 outermost_regime) < 1.0) {
                // Loop over the remaining radial bins. //
                for (int remaining_bin = bin + 1; remaining_bin < NUM_BINS; ++remaining_bin) {
                    if (toomre_parameter(pt5, radial_bins[remaining_bin], viscosity, bh_mass, accretion_rate, spin,
                                         rotation, outermost_regime)
                        > toomre_parameter(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                           outermost_regime)) {
                        terminating("'toomre_parameter' is not a decreasing function of 'radial_bins'");
                    }
                }
                break;
            }
        }
        
        if (outermost_toomre_parameter < 1.0) {
            for (int bin = 0; bin < NUM_BINS; ++bin) {
                toomre_parameters[bin] = toomre_parameter(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate,
                                                          spin, rotation, outermost_regime);
                
                // Check if the local 'toomre_parameters' is below one. //
                if (toomre_parameters[bin] < 1.0) {
                    /* Check if the first 'radial_bin' has a 'toomre_parameters' below one. If yes, that means the
                     * 'Extent_Factor' has increased more than enough in the previous iteration, so find the regime that
                     * is valid right before the 'outermost_regime'. Else, interpolate between the two 'radial_bins' to
                     * get the 'extent'. */
                    if (bin == 0) {
                        /* Find the regime that is valid before the 'outermost_regime' and set it as the new
                         * 'outermost_regime'. Also, restrict 'Extent_Factor' so 'r_outermost' is (1 - 1/'NUM_BINS') of
                         * the current 'outermost_regime', which will prevent the current 'outermost_regime' to appear
                         * again. */
                        double current_max;
                        if (strcmp(outermost_regime, "RES") == 0) {
                            current_max = RES_max;
                            if (GES_max > GFF_max) { outermost_regime = "GES"; }
                            else if (GFF_max > GES_max) { outermost_regime = "GFF"; }
                        } else if (strcmp(outermost_regime, "GES") == 0) {
                            current_max = GES_max;
                            if (RES_max > GFF_max) { outermost_regime = "RES"; }
                            else if (GFF_max > RES_max) { outermost_regime = "GFF"; }
                        } else if (strcmp(outermost_regime, "GFF") == 0) {
                            current_max = GFF_max;
                            if (RES_max > GES_max) { outermost_regime = "RES"; }
                            else if (GES_max > RES_max) { outermost_regime = "GES"; }
                        }
                        
                        /* Set the 'Extent_Factor' so that the next 'r_outermost' is 99% of the current, hence the
                         * current 'outermost_regime' does not exist when calculating the new 'ranges_of_validity' in
                         * the next do-while iteration. */
                        P[pt5].Extent_Factor = 1.0;
                        P[pt5].Extent_Factor = 0.99 * current_max / r_outermost(pt5, bh_mass);
                        extent = -1.0; // Set to -1 so 'Extent_Factor' will not update again in this do-while iteration.
                        innermost_radial_bin = radial_bins[bin];
                        innermost_toomre_parameter = toomre_parameters[bin];
                        break; // Break the 'radial_bins[bin]' loop.
                        
                    } else if (bin > 0 && toomre_parameters[bin - 1] > 1.0) {
                        extent = interpolate(radial_bins[bin - 1], toomre_parameters[bin - 1], radial_bins[bin],
                                             toomre_parameters[bin]);
                        break; // Break the 'radial_bins[bin]' loop.
                    }
                }
            }
        }
        
        /* Check if the current 'outermost_regime' is everywhere stable but previously 'outermost_regime' was not.
         * If yes, the instability happened between radial bins so 'interpolate' the solution. */
        if (outermost_toomre_parameter > 1.0 && extent < 0.0) {
            extent = interpolate(radial_bins[NUM_BINS - 1], outermost_toomre_parameter, innermost_radial_bin,
                                 innermost_toomre_parameter);
        }
        
        /* Check if a self-gravitation radius has not been found. If yes, increase 'Extent_Factor' by a factor of
         * 'NUM_STEPS'. Else, increase 'Extent_Factor' so 'r_outermost' results in 'extent'. */
        if (extent == 0.0) {
            P[pt5].Extent_Factor *= NUM_STEPS;
        } else if (extent > 0.0) {
            // Set 'Extent_Factor' to 1.0 and then re-set it appropriately based on the 'extent/r_outermost' ratio. //
            P[pt5].Extent_Factor = 1.0;
            P[pt5].Extent_Factor = extent / r_outermost(pt5, bh_mass);
            break; // Break the 'Extent_Factor' loop.
        }
        
        // Sanity check, more info in the print statement. //
        if (iterations++ == MAX_ITERATIONS) {terminating("MAX_ITERATIONS has been reached"); }
    } while (extent < 0.0 || fmod(P[pt5].Extent_Factor, 10.0) == 0.0 || fmod(P[pt5].Extent_Factor, 10.0) == 1.0);
    
} // accretion_disc_extent(...)


/**
 * @brief Calculate the range of validity of the local solutions, this information helps build the global structure of
 * an accretion disc as a combination of locally valid regimes.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 **/
void ranges_of_validity(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin,
                        const char *rotation) {
    
    // Declare and initialise the 'radial_bins' array with zeros. //
    double radial_bins[NUM_BINS];
    memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost' with 'radial_bin_width' step. //
    (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
    double radial_bin_width = fabs(log10(r_outermost(pt5, bh_mass)) - log10(r_isco(bh_mass, spin, rotation)))
                              / (NUM_BINS - 1);
    
    // Declare arrays to store the ranges of validity and initialise them with zeros. //
    double GES_radial_bins[NUM_BINS], RES_radial_bins[NUM_BINS], GFF_radial_bins[NUM_BINS];
    memset(GES_radial_bins, 0.0, NUM_BINS * sizeof(double));
    memset(RES_radial_bins, 0.0, NUM_BINS * sizeof(double));
    memset(GFF_radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    /* Check for every 'bin' and every 'regime' if its 'pressure_ratio' and 'opacity_ratio' are both less than one
     * (i.e. if the 'regime' is valid for a given 'bin'). If yes, then store the corresponding 'radial_bins'. */
    // Loop over 'NUM_BINS' radial bins. //
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        if (pressure_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GES") < 1.0
            && opacity_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GES") < 1.0) {
            if (bin == 0) { GES_radial_bins[bin] = r_isco(bh_mass, spin, rotation); }
            else if (bin == NUM_BINS - 1) { GES_radial_bins[bin] = r_outermost(pt5, bh_mass); }
            else { GES_radial_bins[bin] = radial_bins[bin]; }
        }
        if (pressure_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "RES") < 1.0
            && opacity_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "RES") < 1.0) {
            if (bin == 0) { RES_radial_bins[bin] = r_isco(bh_mass, spin, rotation); }
            else if (bin == NUM_BINS - 1) { RES_radial_bins[bin] = r_outermost(pt5, bh_mass); }
            else { RES_radial_bins[bin] = radial_bins[bin]; }
        }
        if (pressure_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GFF") < 1.0
            && opacity_ratio(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GFF") < 1.0) {
            if (bin == 0) { GFF_radial_bins[bin] = r_isco(bh_mass, spin, rotation); }
            else if (bin == NUM_BINS - 1) { GFF_radial_bins[bin] = r_outermost(pt5, bh_mass); }
            else { GFF_radial_bins[bin] = radial_bins[bin]; }
        }
    }
    
    /* 'resolve_overlapping_regimes' and 'resolve_lone_regimes' will make small internal changes to the structure of
     * the accretion disc, which -- given the connection between the extent of the accretion disc and the exact range
     * of validity of each regime -- will slightly move it away from the gravitational stability established in
     * 'accretion_disc_extent'. However, re-calculating the 'accretion_disc_extent' every time even a single radial bin
     * has been changed might lead to an endless repetition to find the 'accretion_disc_extent',
     * 'resolve_overlapping_regimes', 'resolve_lone_regimes', and pass the 'local_tests' and 'global_tests'. */
    
    // If two or more regimes are valid in the same radial bin, split them based on their 'pressure_ratio'. //
    (void) resolve_overlapping_regimes(pt5, viscosity, bh_mass, accretion_rate, spin, rotation, GES_radial_bins,
                                       RES_radial_bins, GFF_radial_bins);
    
    // Before 'resolve_lone_regimes', 'resolve_empty_regimes' (i.e. assign empty bins to adjacent regimes). //
    (void) resolve_empty_regimes(pt5, viscosity, bh_mass, accretion_rate, spin, rotation, GES_radial_bins,
                                 RES_radial_bins, GFF_radial_bins);
    
    // If a regime is only valid in one radial bin, replace it with an adjacent valid regime. //
    (void) resolve_lone_regimes(pt5, viscosity, bh_mass, accretion_rate, spin, rotation, GES_radial_bins,
                                RES_radial_bins, GFF_radial_bins);
    
    // Sanity checks, more info in their print statements. //
    (void) local_tests(pt5, bh_mass, spin, rotation, GES_radial_bins, RES_radial_bins, GFF_radial_bins, __FUNCTION__,
                       __FILE__);
    
    // Set the black hole properties that store the starting and stopping radii to zero. //
    memset(P[pt5].GES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].RES_stops, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_starts, 0.0, NUM_VALID * sizeof(double));
    memset(P[pt5].GFF_stops, 0.0, NUM_VALID * sizeof(double));
    
    // Check if the 'GES' regime is valid somewhere (i.e. if 'GES_sum' is positive). //
    if (P[pt5].GES_sum > 0.0) {
        /* Check if the first entry and 'r_isco' are valid in the same 'radial_bin' (i.e. have difference smaller
         * than the 'radial_bin_width / NUM_BINS'). If yes, then 'GES' starts at 'r_isco', so record that radius and
         * increase 'start' by one. */
        int start = 0, stop = 0; // Declare indices to store starting and stopping radii in P[pt5] properties. //
        if (GES_radial_bins[0] > 0.0
            && fabs(log10(GES_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation))) < radial_bin_width / NUM_BINS) {
            P[pt5].Regime_ISCO = "GES";
            P[pt5].GES_starts[start] = r_isco(bh_mass, spin, rotation), start++;
        }
        
        // Loop over 'NUM_BINS' radial bins excluding the first bin. //
        for (int bin = 1; bin < NUM_BINS; ++bin) {
            /* Check if the previous radial bin is larger than the current and not the ISCO, and the current is zero.
             * If yes, then 'GES' has stopped in the previous radial bin, so record that radius and increase 'stop' by
             * one. */
            if (GES_radial_bins[bin - 1] > GES_radial_bins[bin] && GES_radial_bins[bin] == 0.0
                && fabs(log10(GES_radial_bins[bin - 1]) - log10(r_isco(bh_mass, spin, rotation)))
                   > radial_bin_width / NUM_BINS) {
                P[pt5].GES_stops[stop] = GES_radial_bins[bin - 1], stop++;
            }
            /* Check if the previous radial bin is zero, and the current is positive. If yes, then 'GES' has started
             * again in the current radial bin, so record that radius and increase 'start' by one. */
            if (GES_radial_bins[bin - 1] == 0.0 && GES_radial_bins[bin] > 0.0) {
                P[pt5].GES_starts[start] = GES_radial_bins[bin], start++;
            }
        }
        
        /* Check if the last entry and 'r_outermost' are valid in the same 'radial_bin' (i.e. have difference
         * smaller than the 'radial_bin_width / NUM_BINS'). If yes, then 'GES' ends at the outermost radius. */
        if (GES_radial_bins[NUM_BINS - 1] > 0.0
            && fabs(log10(GES_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass)))
               < radial_bin_width / NUM_BINS) {
            P[pt5].GES_stops[stop] = r_outermost(pt5, bh_mass);
        }
        
        // Sanity check, more info in the print statement.
        if (start > NUM_VALID) {terminating("'start' is larger than 'NUM_VALID' for 'GES"); }
        if (stop > NUM_VALID) {terminating("'stop' is larger than 'NUM_VALID' for 'GES"); }
    }
    
    // Check if the 'RES' regime is valid somewhere (i.e. if 'RES_sum' is positive). //
    if (P[pt5].RES_sum > 0.0) {
        /* Check if the first entry and 'r_isco' are valid in the same 'radial_bin' (i.e. have difference smaller
         * than the 'radial_bin_width / NUM_BINS'). If yes, then 'RES' starts at 'r_isco', so record that radius and
         * increase 'start' by one. */
        int start = 0, stop = 0; // Declare indices to store starting and stopping radii in P[pt5] properties. //
        if (RES_radial_bins[0] > 0.0
            && fabs(log10(RES_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation))) < radial_bin_width / NUM_BINS) {
            P[pt5].Regime_ISCO = "RES";
            P[pt5].RES_starts[start] = r_isco(bh_mass, spin, rotation), start++;
        }
        
        // Loop over 'NUM_BINS' radial bins excluding the first bin. //
        for (int bin = 1; bin < NUM_BINS; ++bin) {
            /* Check if the previous radial bin is larger than the current and not the ISCO, and the current is zero.
             * If yes, then 'RES' has stopped in the previous radial bin, so record that radius and increase 'stop' by
             * one. */
            if (RES_radial_bins[bin - 1] > RES_radial_bins[bin] && RES_radial_bins[bin] == 0.0
                && fabs(log10(RES_radial_bins[bin - 1]) - log10(r_isco(bh_mass, spin, rotation)))
                   > radial_bin_width / NUM_BINS) {
                P[pt5].RES_stops[stop] = RES_radial_bins[bin - 1], stop++;
            }
            /* Check if the previous radial bin is zero, and the current is positive. If yes, then 'RES' has started
             * again in the current radial bin, so record that radius and increase 'start' by one. */
            if (RES_radial_bins[bin - 1] == 0.0 && RES_radial_bins[bin] > 0.0) {
                P[pt5].RES_starts[start] = RES_radial_bins[bin], start++;
            }
        }
        
        /* Check if the last entry and 'r_outermost' are valid in the same 'radial_bin' (i.e. have difference
         * smaller than the 'radial_bin_width / NUM_BINS'). If yes, then 'RES' ends at the outermost radius. */
        if (RES_radial_bins[NUM_BINS - 1] > 0.0
            && fabs(log10(RES_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass)))
               < radial_bin_width / NUM_BINS) {
            P[pt5].RES_stops[stop] = r_outermost(pt5, bh_mass);
        }
        
        // Sanity check, more info in the print statement.
        if (start > NUM_VALID) {terminating("'start' is larger than 'NUM_VALID' for 'RES"); }
        if (stop > NUM_VALID) {terminating("'stop' is larger than 'NUM_VALID' for 'RES"); }
    }
    
    // Check if the 'GFF' regime is valid somewhere (i.e. if 'GFF_sum' is positive). //
    if (P[pt5].GFF_sum > 0.0) {
        /* Check if the first entry and 'r_isco' are valid in the same 'radial_bin' (i.e. have difference smaller
         * than the 'radial_bin_width / NUM_BINS'). If yes, then 'GFF' starts at 'r_isco', so record that radius and
         * increase 'start' by one. */
        int start = 0, stop = 0; // Declare indices to store starting and stopping radii in P[pt5] properties. //
        if (GFF_radial_bins[0] > 0.0
            && fabs(log10(GFF_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation))) < radial_bin_width / NUM_BINS) {
            P[pt5].Regime_ISCO = "GFF";
            P[pt5].GFF_starts[start] = r_isco(bh_mass, spin, rotation), start++;
        }
        
        // Loop over 'NUM_BINS' radial bins excluding the first bin. //
        for (int bin = 1; bin < NUM_BINS; ++bin) {
            /* Check if the previous radial bin is larger than the current and not the ISCO, and the current is zero.
             * If yes, then 'GFF' has stopped in the previous radial bin, so record that radius and increase 'stop' by
             * one. */
            if (GFF_radial_bins[bin - 1] > GFF_radial_bins[bin] && GFF_radial_bins[bin] == 0.0
                && fabs(log10(GFF_radial_bins[bin - 1]) - log10(r_isco(bh_mass, spin, rotation)))
                   > radial_bin_width / NUM_BINS) {
                P[pt5].GFF_stops[stop] = GFF_radial_bins[bin - 1], stop++;
            }
            /* Check if the previous radial bin is zero, and the current is positive. If yes, then 'GFF' has started
             * again in the current radial bin, so record that radius and increase 'start' by one. */
            if (GFF_radial_bins[bin - 1] == 0.0 && GFF_radial_bins[bin] > 0.0) {
                P[pt5].GFF_starts[start] = GFF_radial_bins[bin], start++;
            }
        }
        
        /* Check if the last entry and 'r_outermost' are valid in the same 'radial_bin' (i.e. have difference
         * smaller than the 'radial_bin_width / NUM_BINS'). If yes, then 'GFF' ends at the outermost radius. */
        if (GFF_radial_bins[NUM_BINS - 1] > 0.0
            && fabs(log10(GFF_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass)))
               < radial_bin_width / NUM_BINS) {
            P[pt5].GFF_stops[stop] = r_outermost(pt5, bh_mass);
        }
        
        // Sanity check, more info in the print statement.
        if (start > NUM_VALID) {terminating("'start' is larger than 'NUM_VALID' for 'GFF"); }
        if (stop > NUM_VALID) {terminating("'stop' is larger than 'NUM_VALID' for 'GFF"); }
    }
    
    // Sanity checks, more info in their print statements. //
    (void) global_tests(pt5, radial_bin_width, __FUNCTION__, __FILE__);
    
} // ranges_of_validity(...)


/**
 * @brief Calculate the gravitational radius.
 * @param bh_mass Black hole mass in gram.
 * @return Gravitational radius in centimeter.
 **/
double r_grav(double bh_mass) {
    
    return RA_G_CGS * bh_mass * pow(RA_C_CGS, -2.0);
    
} // r_grav(...)


/**
 * @brief Calculate the radius of the outermost edge of the accretion disc.
 * @param pt5 Index of an active black hole particle.
 * @param bh_mass Black hole mass in gram.
 * @return Outermost radius in centimeter.
 **/
double r_outermost(int pt5, double bh_mass) {
    
    // Fitting formula from Morgan+10 (2010ApJ...712.1129M) multiplied by 'Extent_Factor'. //
    return P[pt5].Extent_Factor * pow(10.0, 15.78) * pow(bh_mass / (1.0e9 * RA_SOLAR_MASS_CGS), 0.8);
    
} // r_outermost(...)


/**
 * @brief Calculate the radius of the innermost stable circular orbit (ISCO).
 * @param bh_mass Black hole mass in gram.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return ISCO radius in centimeter.
 **/
double r_isco(double bh_mass, double spin, const char *rotation) {
    
    double alignment; // Declare local variable.
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Calculate the spin-depended corrections from Bardeen+72 (1972ApJ...178..347B). //
    double z_1 = 1.0 + pow(1.0 - pow(spin, 2.0), 1.0 / 3.0) * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double z_2 = sqrt(3.0 * pow(spin, 2.0) + pow(z_1, 2.0));
    double z = 3.0 + z_2 - alignment * sqrt((3.0 - z_1) * (3.0 + z_1 + 2.0 * z_2));
    
    return r_grav(bh_mass) * z;
    
} // r_isco(...)


/**
 * @brief Calculate the radius of the photon orbit based on Eq. 2.18 from 1972ApJ...178..347B.
 * @param bh_mass Black hole mass in gram.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Radius in centimeter.
 **/
double r_photon(double bh_mass, double spin, const char *rotation) {
    
    double alignment; // Declare local variable.
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    return 2.0 * r_grav(bh_mass) * (1.0 + cos(2.0 / 3.0 * acos(-1.0 * alignment * spin)));
    
} // r_photon(...)


/**
 * @brief Calculate the warp radius of an accretion disc that is misaligned with respect to the black hole's equator.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Warp radius in centimeter.
 **/
void r_warp(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin, const char *rotation) {
    
    // Check if the accretion disc exists for at least one time-step. If yes, calculate its warp radius. //
    if (P[pt5].AD_Lifetime > 0.0) {
        double t_precession = 0.0, t_bardeen_petterson = 0.0; // Declare local variables.
        
        // Declare and initialise the 'radial_bins' array with zeros. //
        double radial_bins[NUM_BINS];
        memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
        
        // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost'. //
        (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS,
                        radial_bins);
        
        // Calculate the black hole's time-step (i.e. time since it was last active) and convert to cgs units. //
        double dt = (P[pt5].TimeBin ? (1 << P[pt5].TimeBin) : 0.0) * All.Timebase_interval / hubble_parameter();
        double dt_cgs = convert_code_to_cgs_units("time", dt);
        
        // Loop over 'NUM_BINS' radial bins excluding the first bin (integral is not valid at the ISCO). //
        for (int bin = 1; bin < NUM_BINS; ++bin) {
            // Calculate the precession and Bardeen-Petterson timescales for a given radial bin. //
            t_precession = precession_timescale(pt5, radial_bins[bin]);
            t_bardeen_petterson = bardeen_petterson_timescale(pt5, radial_bins[bin], viscosity, bh_mass, accretion_rate,
                                                              spin, rotation);
            
            /* Alignment occurs in radii where the precession timescale is smaller than the Bardeen-Petterson timescale,
             * and only if both timescales are smaller than the lifetime of the accretion disc. */
            if (radial_bins[bin] >= P[pt5].R_Warp && t_bardeen_petterson <= P[pt5].AD_Lifetime
                && t_precession <= P[pt5].AD_Lifetime) {
                // Check if the Lifetime is the same as the current timestep with some tolerance. If yes, then this is
                // the first timestep for this AD so the whole AD should be integrated as misaligned (i.e. R_Warp = 0).
                if (fabs(dt_cgs - P[pt5].AD_Lifetime) > dt_cgs / UNIT_TEST_TOLERANCE) { P[pt5].R_Warp = 0.0; }
                else if (t_precession <= t_bardeen_petterson) { P[pt5].R_Warp = radial_bins[bin]; }
            }
        }
    }
    
} // r_warp(...)