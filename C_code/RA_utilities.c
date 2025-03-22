/*
 This file is part of a public black hole accretion model.
` Copyright (c) 2025. Dimitrios Irodotou (di.irodotou@gmail.com) and contributing co-authors.
 
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
 * @brief Perform unit tests and initialise the black hole and accretion disc properties.
 **/
void initialise_setup(void) {
    
    // Loop over all active particles and check if the active particle is a black hole with positive mass. //
    for (int pt5 = FirstActiveParticle; pt5 >= 0; pt5 = NextActiveParticle[pt5]) {
        if (P[pt5].Type == 5 && P[pt5].BH_Mass > 0.0) {
            // Perform unit tests. //
            (void) radii_unit_tests(pt5);
            (void) accretion_disc_unit_tests(pt5);
            
            // Initialise the black hole spin with a random value. //
            (void) initialise_black_hole_spin(pt5);
            
            // Create a log file per black hole and log the initial values of some properties. //
            (void) initialise_log_files(pt5);
            (void) log_details(pt5, "ICs");
            
            // Set the accretion disc properties to their default values. //
            (void) reset_accretion_disc(pt5);
        }
    }
    
} // initialise_setup(...)



/**
 * @brief Initialise a directory and text file to keep the 'relativistic_accretion_logs' from all CPUs.
 * @param pt5 Index of an active black hole particle.
 **/
void initialise_log_files(int pt5) {
    
    char mode[2], buf[1000]; // Declare local variables.
    
    /* Check if the simulation starts from initial conditions (instead of from restart or snapshot files). If yes,
     * mode 'w' is needed to create an empty file for writing. Else, mode 'a' is needed to append to an existing file.
     * 'RestartFlag' is passed as an argument from command line when running the code. 0 is for starting from initial
     * conditions, 1 is for resuming from restart files, while 2 is for resuming from a snapshot file. */
    if (RestartFlag == 0) { strcpy(mode, "w"); }
    else { strcpy(mode, "a"); }
    
    /* Create a directory 'relativistic_accretion_logs' and a file 'ID_*.txt' inside the directory 'All.OutputDir'
     * for each black hole. */
    sprintf(buf, "%srelativistic_accretion_logs", All.OutputDir);
    mkdir(buf, 02755); // 02755 permission means readable and executable by everyone but writable only by the owner.
    
    // Sanity check, more info in the print statement. //
    sprintf(buf, "%srelativistic_accretion_logs/ID_%d.txt", All.OutputDir, P[pt5].ID);
    if (!(FdAccretionDiscsAll = fopen(buf, mode))) {
        terminating("Cannot open 'relativistic_accretion_logs/ID_.txt' file");
    }
    
    // Sanity check, more info in the print statement. //
    sprintf(buf, "%srelativistic_accretion_logs/capture.txt", All.OutputDir);
    if (!(FdHoyleLyttletonAll = fopen(buf, mode))) {
        terminating("Cannot open 'relativistic_accretion_logs/r_hl.txt' file");
    }
    
} // initialise_log_files(...)


/**
 * @brief Append information in 'relativistic_accretion_logs/ID_*.txt' from all CPUs that have an active black hole.
 * @param pt5 Index of an active black hole particle.
 * @param io Either "ICS", "in", or "out".
 **/
void log_details(int pt5, const char *io) {
    
    char buf[1000]; // Declare local variable.
    
    // Sanity check, more info in the print statement. //
    sprintf(buf, "%srelativistic_accretion_logs/ID_%d.txt", All.OutputDir, P[pt5].ID);
    if (!(FdAccretionDiscsAll = fopen(buf, "a"))) {
        terminating("Cannot open 'relativistic_accretion_logs/ID_.txt' file");
    }
    
    // Calculate the black hole's time-step (i.e. time since it was last active) in code units. //
    double dt = (P[pt5].TimeBin ? (1 << P[pt5].TimeBin) : 0.0) * All.Timebase_interval / hubble_parameter();
    
    fprintf(FdAccretionDiscsAll,
            "Sync-Point:%d, Time:%.20e, dt:%.20e, Task:%d, io:%s, rotation:%s, "
            "ID:%d, Mass:%.20e, BH_Mass:%.20e, AD_Mass:%.20e, BH_accreted_Mass:%.20e, BH_Mdot:%.20e, AD_Mdot:%.20e, "
            "Edd_Mdot:%.20e, BHL_Mdot:%.20e, Epsilon_r:%.20e, Epsilon_spin:%.20e, AD_Energy:%.20e, "
            "spin_parameter:%.20e, P[pt5].Spin[0]:%.20e, P[pt5].Spin[1]:%.20e, P[pt5].Spin[2]:%.20e, "
            "P[pt5].AD_AngMomentum[0]:%.20e, P[pt5].AD_AngMomentum[1]:%.20e, P[pt5].AD_AngMomentum[2]:%.20e \n ",
            All.NumCurrentTiStep, All.Time, dt, ThisTask, io, assess_alignment(pt5),
            P[pt5].ID, P[pt5].Mass, P[pt5].BH_Mass, P[pt5].AD_Mass, P[pt5].b4.dBH_accreted_Mass, P[pt5].BH_Mdot,
            P[pt5].AD_Mdot, P[pt5].Edd_Mdot, bondi_hoyle_accretion_rate(pt5), P[pt5].Epsilon_r, P[pt5].Epsilon_spin,
            P[pt5].AD_Energy, calculate_spin_parameter(pt5), P[pt5].Spin[0], P[pt5].Spin[1], P[pt5].Spin[2],
            P[pt5].AD_AngMomentum[0], P[pt5].AD_AngMomentum[1], P[pt5].AD_AngMomentum[2]);
    
    fclose(FdAccretionDiscsAll); // Write the output buffer (i.e. flush) and close the file.
    
} // log_details(...)


/**
 * @brief Re-set the 'Extent_Factor' to one and the rest of the accretion disc properties to zero.
 * @param pt5 Index of an active black hole particle.
 **/
void reset_accretion_disc(int pt5) {
    
    P[pt5].Regime_ISCO = "";
    P[pt5].Extent_Factor = 1.0;
    P[pt5].Drain_Flag = 0, P[pt5].Form_Flag = 0;
    P[pt5].Epsilon_r = 0.0, P[pt5].Epsilon_spin = 0.0;
    P[pt5].GES_sum = 0.0, P[pt5].RES_sum = 0.0, P[pt5].GFF_sum = 0.0;
    P[pt5].T_Viscous = 0.0, P[pt5].T_Deplete = 0.0, P[pt5].AD_Lifetime = 0.0;
    P[pt5].AD_AngMomentum[0] = 0.0, P[pt5].AD_AngMomentum[1] = 0.0, P[pt5].AD_AngMomentum[2] = 0.0;
    P[pt5].AD_AngMomentumDir[0] = 0.0, P[pt5].AD_AngMomentumDir[1] = 0.0, P[pt5].AD_AngMomentumDir[2] = 0.0;
    P[pt5].AD_Mass = 0.0, P[pt5].R_Warp = 0.0, P[pt5].AD_Mdot = 0.0, P[pt5].AD_Energy = 0.0, P[pt5].Edd_Mdot = 0.0;
    
} // reset_accretion_disc(...)


/**
 * @brief Convert standard GADGET code units (length: 1 kpc, mass: 1e10 Msun/h, time:9.8e8 yr/h) to cgs.
 * @param property Select between 'length', 'mass', 'time', 'mass_rate', or 'angmomentum'.
 * @param value Value of the property.
 * @return Value in cgs units.
 **/
double convert_code_to_cgs_units(const char *property, double value) {
    
    if ((strcmp(property, "length") == 0)) {
        double centimeter = (All.UnitLength_in_cm / All.HubbleParam);
        return value * centimeter;
        
    } else if (strcmp(property, "mass") == 0) {
        double gram = (All.UnitMass_in_g / All.HubbleParam);
        return value * gram;
        
    } else if (strcmp(property, "time") == 0) {
        double second = (All.UnitTime_in_s / All.HubbleParam);
        return value * second;
        
    } else if (strcmp(property, "mass_rate") == 0) {
        double gram = (All.UnitMass_in_g / All.HubbleParam);
        double second = (All.UnitTime_in_s / All.HubbleParam);
        return value * gram / second;
        
    } else if (strcmp(property, "angmomentum") == 0) {
        double velocity = All.UnitVelocity_in_cm_per_s;
        double gram = (All.UnitMass_in_g / All.HubbleParam);
        double centimeter = (All.UnitLength_in_cm / All.HubbleParam);
        return value * gram * velocity * centimeter;
        
    } else {terminating("Invalid 'property'"); }
    
} // convert_code_to_cgs_units(...)


/**
 * @brief Convert cgs to standard GADGET code units (length: 1 kpc, mass: 1e10 Msun/h, time:9.8e8 yr/h).
 * @param property Select between 'length', 'mass', 'time', 'mass_rate', 'angmomentum', 'gravity', or 'light'.
 * @param value Value of the property.
 * @return Value in code units
 **/
double convert_cgs_to_code_units(const char *property, double value) {
    
    if ((strcmp(property, "length") == 0)) {
        double code_unit_length = (All.HubbleParam / All.UnitLength_in_cm);
        return value * code_unit_length;
        
    } else if (strcmp(property, "mass") == 0) {
        double code_unit_mass = (All.HubbleParam / All.UnitMass_in_g);
        return value * code_unit_mass;
        
    } else if (strcmp(property, "time") == 0) {
        double code_unit_time = (All.HubbleParam / All.UnitTime_in_s);
        return value * code_unit_time;
        
    } else if (strcmp(property, "mass_rate") == 0) {
        double code_unit_time = (All.HubbleParam / All.UnitTime_in_s);
        double code_unit_mass = (All.HubbleParam / All.UnitMass_in_g);
        return value * code_unit_mass / code_unit_time;
        
    } else if (strcmp(property, "angmomentum") == 0) {
        double code_unit_velocity = 1.0 / All.UnitVelocity_in_cm_per_s;
        double code_unit_mass = (All.HubbleParam / All.UnitMass_in_g);
        double code_unit_length = (All.HubbleParam / All.UnitLength_in_cm);
        return value * code_unit_mass * code_unit_velocity * code_unit_length;
        
    } else if (strcmp(property, "gravity") == 0) {
        return value / pow(All.UnitLength_in_cm, 3.0) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2.0);
        
    } else if (strcmp(property, "light") == 0) {
        return value / All.UnitVelocity_in_cm_per_s;
        
    } else {terminating("Invalid 'property'"); }
    
} // convert_cgs_to_code_units(...)


/**
 * @brief Create an array with evenly spaced numbers in log space.
 * @param start Starting exponent value.
 * @param stop Stopping exponent value.
 * @param num Length of array.
 * @param array Array to store data.
 **/
void *logspace(double start, double stop, int num, double array[]) {
    
    double step = (stop - start) / (num - 1.0); // Calculate the step size.
    
    // Loop over 'num' elements excluding the last element and populate them with evenly spaced numbers in log space. //
    for (int element = 0; element < num - 1; ++element) { array[element] = pow(10.0, start + element * step); }
    array[num - 1] = pow(10.0, stop); // Set the value of the last element to ten to the power of 'stop'.
    
} // logspace(...)


/**
 * @brief Find the maximum value inside an array of values.
 * @param size Size of the array.
 * @param array Array to find the maximum value in.
 * @return Maximum value of the array.
 **/
double maximum_array_value(int size, const double array[]) {
    
    double max_value = array[0]; // Initialize maximum element.
    
    // Loop over array elements and compare every element with the latest 'max_value'. //
    for (int element = 1; element < size; ++element) { if (array[element] > max_value) { max_value = array[element]; }}
    
    return max_value;
    
} // maximum_array_value(...)


/**
 * @brief Linearly interpolate the x-coordinate of a given y-coordinate between two points.
 * @param x1 x-coordinate of first point.
 * @param y1 y-coordinate of first point.
 * @param x2 x-coordinate of second point.
 * @param y2 y-coordinate of second point.
 * @return Interpolated x-coordinate.
 **/
double interpolate(double x1, double y1, double x2, double y2) {
    
    double y = 1.0; // Represents the marginally stable solution with Toomre parameter equal to one.
    return x1 + (y - y1) * (x2 - x1) / (y2 - y1);
    
} // interpolate(...)


/**
 * @brief If two or more regimes are valid in the same radial bin, assign the bins accordingly to the adjacent regimes.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param GES_radial_bins Valid 'radial_bins' for the 'GES' regime.
 * @param RES_radial_bins Valid 'radial_bins' for the 'RES' regime.
 * @param GFF_radial_bins Valid 'radial_bins' for the 'GFF' regime.
 **/
void *resolve_overlapping_regimes(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin,
                                  const char *rotation, double GES_radial_bins[NUM_BINS],
                                  double RES_radial_bins[NUM_BINS], double GFF_radial_bins[NUM_BINS]) {
    
    // Declare and initialise the 'radial_bins' array with zeros. //
    double radial_bins[NUM_BINS];
    memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost' with 'radial_bin_width' step. //
    (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
    double radial_bin_width = fabs(log10(r_outermost(pt5, bh_mass)) - log10(r_isco(bh_mass, spin, rotation)))
                              / (NUM_BINS - 1);
    
    // Loop over 'NUM_BINS' radial bins. //
    for (int bin = 0; bin < NUM_BINS; ++bin) {
        /* Check if two regimes are valid in the same 'radial_bins' (i.e. have difference smaller than
         * 'radial_bin_width / NUM_BINS'). If yes, then store the 'bin' this happens and how many times. */
        if (GES_radial_bins[bin] > 0.0 && RES_radial_bins[bin] > 0.0
            && fabs(log10(GES_radial_bins[bin]) - log10(RES_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            
            if (pressure_ratio(pt5, GES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GES")
                <
                pressure_ratio(pt5, RES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "RES")) {
                RES_radial_bins[bin] = 0.0;
                
            } else if (pressure_ratio(pt5, RES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                      "RES")
                       < pressure_ratio(pt5, GES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                        "GES")) {
                GES_radial_bins[bin] = 0.0;
            }
        }
        if (GES_radial_bins[bin] > 0.0 && GFF_radial_bins[bin] > 0.0
            && fabs(log10(GES_radial_bins[bin]) - log10(GFF_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            
            if (pressure_ratio(pt5, GES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GES")
                <
                pressure_ratio(pt5, GFF_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GFF")) {
                GFF_radial_bins[bin] = 0.0;
                
            } else if (pressure_ratio(pt5, GFF_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                      "GFF")
                       < pressure_ratio(pt5, GES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                        "GES")) {
                GES_radial_bins[bin] = 0.0;
            }
        }
        if (RES_radial_bins[bin] > 0.0 && GFF_radial_bins[bin] > 0.0
            && fabs(log10(RES_radial_bins[bin]) - log10(GFF_radial_bins[bin])) < radial_bin_width / NUM_BINS) {
            
            if (pressure_ratio(pt5, RES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "RES")
                <
                pressure_ratio(pt5, GFF_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation, "GFF")) {
                GFF_radial_bins[bin] = 0.0;
                
            } else if (pressure_ratio(pt5, GFF_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                      "GFF")
                       < pressure_ratio(pt5, RES_radial_bins[bin], viscosity, bh_mass, accretion_rate, spin, rotation,
                                        "RES")) {
                RES_radial_bins[bin] = 0.0;
            }
        }
    }
    
} // resolve_overlapping_regimes(...)


/**
 * @brief If a radial bin is empty, assign it to the adjacent regime with the lowest 'pressure_ratio'.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param GES_radial_bins Valid 'radial_bins' for the 'GES' regime.
 * @param RES_radial_bins Valid 'radial_bins' for the 'RES' regime.
 * @param GFF_radial_bins Valid 'radial_bins' for the 'GFF' regime.
 **/
void *resolve_empty_regimes(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin,
                            const char *rotation, double GES_radial_bins[NUM_BINS],
                            double RES_radial_bins[NUM_BINS], double GFF_radial_bins[NUM_BINS]) {
    
    // Declare and initialise the 'radial_bins' array with zeros. //
    double radial_bins[NUM_BINS];
    memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost' with 'radial_bin_width' step. //
    (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
    
    /* Check if none of the 'regimes' are valid at the ISCO. If yes, then find the valid regime that has the minimum
     * 'pressure_ratio' in the bin right after the ISCO and assign it to 'r_isco'. */
    if (GES_radial_bins[0] == 0.0 && RES_radial_bins[0] == 0.0 && GFF_radial_bins[0] == 0.0) {
        double GES_pressure_ratio_after = 0.0, RES_pressure_ratio_after = 0.0, GFF_pressure_ratio_after = 0.0;
        
        if (GES_radial_bins[1] > 0.0) {
            GES_pressure_ratio_after = pressure_ratio(pt5, GES_radial_bins[1], viscosity, bh_mass, accretion_rate, spin,
                                                      rotation, "GES");
        }
        if (RES_radial_bins[1] > 0.0) {
            RES_pressure_ratio_after = pressure_ratio(pt5, RES_radial_bins[1], viscosity, bh_mass, accretion_rate, spin,
                                                      rotation, "RES");
        }
        if (GFF_radial_bins[1] > 0.0) {
            GFF_pressure_ratio_after = pressure_ratio(pt5, GFF_radial_bins[1], viscosity, bh_mass, accretion_rate, spin,
                                                      rotation, "GFF");
        }
        
        if (GES_pressure_ratio_after > 0.0 && GES_pressure_ratio_after < RES_pressure_ratio_after
            && GES_pressure_ratio_after < GFF_pressure_ratio_after) {
            GES_radial_bins[0] = r_isco(bh_mass, spin, rotation);
        }
        
        if (RES_pressure_ratio_after > 0.0 && RES_pressure_ratio_after < GES_pressure_ratio_after
            && RES_pressure_ratio_after < GFF_pressure_ratio_after) {
            RES_radial_bins[0] = r_isco(bh_mass, spin, rotation);
        }
        
        if (GFF_pressure_ratio_after > 0.0 && GFF_pressure_ratio_after < GES_pressure_ratio_after
            && GFF_pressure_ratio_after < RES_pressure_ratio_after) {
            GFF_radial_bins[0] = r_isco(bh_mass, spin, rotation);
        }
    }
    
    // Loop over 'NUM_BINS' radial bins excluding the first and the last bins. //
    for (int bin = 1; bin < NUM_BINS - 1; ++bin) {
        /* Check if none of the 'regimes' are valid at a given 'bin'. If yes, then find the regime that has the minimum
         * 'pressure_ratio' between all valid regimes in the bin just before and right after and assign it to the
         * previously empty bin. */
        if (GES_radial_bins[bin] == 0.0 && RES_radial_bins[bin] == 0.0 && GFF_radial_bins[bin] == 0.0) {
            
            double GES_pressure_ratio_before = 1.0, RES_pressure_ratio_before = 1.0, GFF_pressure_ratio_before = 1.0;
            double GES_pressure_ratio_after = 1.0, RES_pressure_ratio_after = 1.0, GFF_pressure_ratio_after = 1.0;
            
            if (GES_radial_bins[bin - 1] > 0.0) {
                GES_pressure_ratio_before = pressure_ratio(pt5, GES_radial_bins[bin - 1], viscosity, bh_mass,
                                                           accretion_rate, spin, rotation, "GES");
            }
            if (RES_radial_bins[bin - 1] > 0.0) {
                RES_pressure_ratio_before = pressure_ratio(pt5, RES_radial_bins[bin - 1], viscosity, bh_mass,
                                                           accretion_rate, spin, rotation, "RES");
            }
            if (GFF_radial_bins[bin - 1] > 0.0) {
                GFF_pressure_ratio_before = pressure_ratio(pt5, GFF_radial_bins[bin - 1], viscosity, bh_mass,
                                                           accretion_rate, spin, rotation, "GFF");
            }
            
            if (GES_radial_bins[bin + 1] > 0.0) {
                GES_pressure_ratio_after = pressure_ratio(pt5, GES_radial_bins[bin + 1], viscosity, bh_mass,
                                                          accretion_rate, spin, rotation, "GES");
            }
            if (RES_radial_bins[bin + 1] > 0.0) {
                RES_pressure_ratio_after = pressure_ratio(pt5, RES_radial_bins[bin + 1], viscosity, bh_mass,
                                                          accretion_rate, spin, rotation, "RES");
            }
            if (GFF_radial_bins[bin + 1] > 0.0) {
                GFF_pressure_ratio_after = pressure_ratio(pt5, GFF_radial_bins[bin + 1], viscosity, bh_mass,
                                                          accretion_rate, spin, rotation, "GFF");
            }
            
            if (GES_pressure_ratio_before < RES_pressure_ratio_before
                && GES_pressure_ratio_before < GFF_pressure_ratio_before
                && GES_pressure_ratio_before < GES_pressure_ratio_after
                && GES_pressure_ratio_before < RES_pressure_ratio_after
                && GES_pressure_ratio_before < GFF_pressure_ratio_after) { GES_radial_bins[bin] = radial_bins[bin]; }
            
            if (RES_pressure_ratio_before < GES_pressure_ratio_before
                && RES_pressure_ratio_before < GFF_pressure_ratio_before
                && RES_pressure_ratio_before < GES_pressure_ratio_after
                && RES_pressure_ratio_before < RES_pressure_ratio_after
                && RES_pressure_ratio_before < GFF_pressure_ratio_after) { RES_radial_bins[bin] = radial_bins[bin]; }
            
            if (GFF_pressure_ratio_before < GES_pressure_ratio_before
                && GFF_pressure_ratio_before < RES_pressure_ratio_before
                && GFF_pressure_ratio_before < GES_pressure_ratio_after
                && GFF_pressure_ratio_before < RES_pressure_ratio_after
                && GFF_pressure_ratio_before < GFF_pressure_ratio_after) { GFF_radial_bins[bin] = radial_bins[bin]; }
            
            if (GES_pressure_ratio_after < RES_pressure_ratio_after
                && GES_pressure_ratio_after < GFF_pressure_ratio_after
                && GES_pressure_ratio_after < GES_pressure_ratio_before
                && GES_pressure_ratio_after < RES_pressure_ratio_before
                && GES_pressure_ratio_after < GFF_pressure_ratio_before) { GES_radial_bins[bin] = radial_bins[bin]; }
            
            if (RES_pressure_ratio_after < GES_pressure_ratio_after
                && RES_pressure_ratio_after < GFF_pressure_ratio_after
                && RES_pressure_ratio_after < GES_pressure_ratio_before
                && RES_pressure_ratio_after < RES_pressure_ratio_before
                && RES_pressure_ratio_after < GFF_pressure_ratio_before) { RES_radial_bins[bin] = radial_bins[bin]; }
            
            if (GFF_pressure_ratio_after < GES_pressure_ratio_after
                && GFF_pressure_ratio_after < RES_pressure_ratio_after
                && GFF_pressure_ratio_after < GES_pressure_ratio_before
                && GFF_pressure_ratio_after < RES_pressure_ratio_before
                && GFF_pressure_ratio_after < GFF_pressure_ratio_before) { GFF_radial_bins[bin] = radial_bins[bin]; }
        }
    }
    
    /* Check if none of the 'regimes' are valid at the outermost edge. If yes, then find the valid regime that has the
     * minimum 'pressure_ratio' in the bin just before the edge and assign it to 'r_outermost'. */
    if (GES_radial_bins[NUM_BINS - 1] == 0.0 && RES_radial_bins[NUM_BINS - 1] == 0.0
        && GFF_radial_bins[NUM_BINS - 1] == 0.0) {
        double GES_pressure_ratio_before = 0.0, RES_pressure_ratio_before = 0.0, GFF_pressure_ratio_before = 0.0;
        
        if (GES_radial_bins[NUM_BINS - 2] > 0.0) {
            GES_pressure_ratio_before = pressure_ratio(pt5, GES_radial_bins[NUM_BINS - 2], viscosity, bh_mass,
                                                       accretion_rate, spin, rotation, "GES");
        }
        if (RES_radial_bins[NUM_BINS - 2] > 0.0) {
            RES_pressure_ratio_before = pressure_ratio(pt5, RES_radial_bins[NUM_BINS - 2], viscosity, bh_mass,
                                                       accretion_rate, spin, rotation, "RES");
        }
        if (GFF_radial_bins[NUM_BINS - 2] > 0.0) {
            GFF_pressure_ratio_before = pressure_ratio(pt5, GFF_radial_bins[NUM_BINS - 2], viscosity, bh_mass,
                                                       accretion_rate, spin, rotation, "GFF");
        }
        
        if (GES_pressure_ratio_before > 0.0 && GES_pressure_ratio_before < RES_pressure_ratio_before
            && GES_pressure_ratio_before < GFF_pressure_ratio_before) {
            GES_radial_bins[NUM_BINS - 1] = r_outermost(pt5, bh_mass);
        }
        
        if (RES_pressure_ratio_before > 0.0 && RES_pressure_ratio_before < GES_pressure_ratio_before
            && RES_pressure_ratio_before < GFF_pressure_ratio_before) {
            RES_radial_bins[NUM_BINS - 1] = r_outermost(pt5, bh_mass);
        }
        
        if (GFF_pressure_ratio_before > 0.0 && GFF_pressure_ratio_before < GES_pressure_ratio_before
            && GFF_pressure_ratio_before < RES_pressure_ratio_before) {
            GFF_radial_bins[NUM_BINS - 1] = r_outermost(pt5, bh_mass);
        }
    }
    
} // resolve_empty_regimes(...)


/**
 * @brief If a regime is valid in only one radial bin, assign the bin accordingly to an adjacent regimes.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param GES_radial_bins Valid 'radial_bins' for the 'GES' regime.
 * @param RES_radial_bins Valid 'radial_bins' for the 'RES' regime.
 * @param GFF_radial_bins Valid 'radial_bins' for the 'GFF' regime.
 **/
void *resolve_lone_regimes(int pt5, double viscosity, double bh_mass, double accretion_rate, double spin,
                           const char *rotation, double GES_radial_bins[NUM_BINS],
                           double RES_radial_bins[NUM_BINS], double GFF_radial_bins[NUM_BINS]) {
    
    // Declare and initialise the 'radial_bins' array with zeros. //
    double radial_bins[NUM_BINS];
    memset(radial_bins, 0.0, NUM_BINS * sizeof(double));
    
    // Calculate evenly log-spaced 'radial_bins' from 'r_isco' to 'r_outermost' with 'radial_bin_width' step. //
    (void) logspace(log10(r_isco(bh_mass, spin, rotation)), log10(r_outermost(pt5, bh_mass)), NUM_BINS, radial_bins);
    double radial_bin_width = fabs(log10(r_outermost(pt5, bh_mass)) - log10(r_isco(bh_mass, spin, rotation)))
                              / (NUM_BINS - 1);
    
    /* Check if a regime is valid at the ISCO but not in the second bin. If yes, then extend the adjacent valid regime
     * to 'r_isco' to avoid shrinking the accretion disc. */
    if (fabs(log10(GES_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation)))
        < radial_bin_width / NUM_BINS && GES_radial_bins[1] == 0) {
        
        if (RES_radial_bins[1] > 0.0 && GFF_radial_bins[1] == 0.0) {
            RES_radial_bins[0] = GES_radial_bins[0];
            GES_radial_bins[0] = 0.0;
        } else if (GFF_radial_bins[1] > 0.0 && RES_radial_bins[1] == 0.0) {
            GFF_radial_bins[0] = GES_radial_bins[0];
            GES_radial_bins[0] = 0.0;
        } else {terminating("Neither 'GES' nor 'GFF' are valid in the bin after the ISCO"); }
    }
    
    if (fabs(log10(RES_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation)))
        < radial_bin_width / NUM_BINS && RES_radial_bins[1] == 0) {
        
        if (GES_radial_bins[1] > 0.0 && GFF_radial_bins[1] == 0.0) {
            GES_radial_bins[0] = RES_radial_bins[0];
            RES_radial_bins[0] = 0.0;
        } else if (GFF_radial_bins[1] > 0.0 && GES_radial_bins[1] == 0.0) {
            GFF_radial_bins[0] = RES_radial_bins[0];
            RES_radial_bins[0] = 0.0;
        } else {terminating("Neither 'RES' nor 'GFF' are valid in the bin after the ISCO"); }
    }
    
    if (fabs(log10(GFF_radial_bins[0]) - log10(r_isco(bh_mass, spin, rotation)))
        < radial_bin_width / NUM_BINS && GFF_radial_bins[1] == 0) {
        
        if (GES_radial_bins[1] > 0.0 && RES_radial_bins[1] == 0.0) {
            GES_radial_bins[0] = GFF_radial_bins[0];
            GFF_radial_bins[0] = 0.0;
        } else if (RES_radial_bins[1] > 0.0 && GES_radial_bins[1] == 0.0) {
            RES_radial_bins[0] = GFF_radial_bins[0];
            GFF_radial_bins[0] = 0.0;
        } else {terminating("Neither 'GES' nor 'RES' are valid in the bin after the ISCO"); }
    }
    
    /* Check if a regime is valid in the last but not the second to last bin. If yes, then extend the adjacent valid
     * regime to 'r_outermost' to avoid shrinking the accretion disc. */
    if (fabs(log10(GES_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass))) < radial_bin_width / NUM_BINS
        && GES_radial_bins[NUM_BINS - 2] == 0) {
        
        if (RES_radial_bins[NUM_BINS - 2] > 0.0 && GFF_radial_bins[NUM_BINS - 2] == 0.0) {
            RES_radial_bins[NUM_BINS - 1] = GES_radial_bins[NUM_BINS - 1];
            GES_radial_bins[NUM_BINS - 1] = 0.0;
        } else if (GFF_radial_bins[NUM_BINS - 2] > 0.0 && RES_radial_bins[NUM_BINS - 2] == 0.0) {
            GFF_radial_bins[NUM_BINS - 1] = GES_radial_bins[NUM_BINS - 1];
            GES_radial_bins[NUM_BINS - 1] = 0.0;
        } else {terminating("Neither 'RES' nor 'GFF' are valid in the second to last bin."); }
    }
    
    if (fabs(log10(RES_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass))) < radial_bin_width / NUM_BINS
        && RES_radial_bins[NUM_BINS - 2] == 0) {
        
        if (GES_radial_bins[NUM_BINS - 2] > 0.0 && GFF_radial_bins[NUM_BINS - 2] == 0.0) {
            GES_radial_bins[NUM_BINS - 1] = RES_radial_bins[NUM_BINS - 1];
            RES_radial_bins[NUM_BINS - 1] = 0.0;
        } else if (GFF_radial_bins[NUM_BINS - 2] > 0.0 && GES_radial_bins[NUM_BINS - 2] == 0.0) {
            GFF_radial_bins[NUM_BINS - 1] = RES_radial_bins[NUM_BINS - 1];
            RES_radial_bins[NUM_BINS - 1] = 0.0;
        } else {terminating("Neither 'GES' nor 'GFF' are valid in the second to last bin."); }
    }
    
    if (fabs(log10(GFF_radial_bins[NUM_BINS - 1]) - log10(r_outermost(pt5, bh_mass))) < radial_bin_width / NUM_BINS
        && GFF_radial_bins[NUM_BINS - 2] == 0) {
        
        if (GES_radial_bins[NUM_BINS - 2] > 0.0 && RES_radial_bins[NUM_BINS - 2] == 0.0) {
            GES_radial_bins[NUM_BINS - 1] = GFF_radial_bins[NUM_BINS - 1];
            GFF_radial_bins[NUM_BINS - 1] = 0.0;
        } else if (RES_radial_bins[NUM_BINS - 2] > 0.0 && GES_radial_bins[NUM_BINS - 2] == 0.0) {
            RES_radial_bins[NUM_BINS - 1] = GFF_radial_bins[NUM_BINS - 1];
            GFF_radial_bins[NUM_BINS - 1] = 0.0;
        } else {terminating("Neither 'GES' nor 'RES' are valid in the second to last bin."); }
    }
    
} // resolve_lone_regimes(...)


/**
 * @brief Integrate the different profiles for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param mode Select between 'spin', 't_viscous', or 'angmomentum'.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param start Lower limit of the integration.
 * @param stop Upper limit of the integration.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Angular momentum in centimeter^2 gram second^-1.
 * @return
 **/
double integrate(int pt5, const char *mode, double viscosity, double bh_mass, double accretion_rate, double spin,
                 double start, double stop, const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double result, error;
    double epsabs = 0.0, epsrel = 1.0e-7; // The desired absolute and relative error limits.
    struct integrand_parameters custom_parameters; // Declare a structure for storing and passing 'custom_parameters'.
    
    // Allocate a workspace 'w' to hold 'MAX_INTERVALS'. //
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAX_INTERVALS);
    
    // Declare and initialize the integrand function and its parameters. //
    gsl_function F;
    if (strcmp(mode, "spin") == 0) { F.function = &spin_integrand; }
    else if (strcmp(mode, "t_viscous") == 0) { F.function = &t_viscous_integrand; }
    else if (strcmp(mode, "angmomentum") == 0) { F.function = &angmomentum_integrand; }
    else {terminating("Invalid 'mode'"); }
    
    F.params = &custom_parameters;
    custom_parameters.pt5 = pt5;
    custom_parameters.viscosity = viscosity;
    custom_parameters.bh_mass = bh_mass;
    custom_parameters.accretion_rate = accretion_rate;
    custom_parameters.spin = spin;
    custom_parameters.rotation = rotation;
    custom_parameters.regime = regime;
    
    // Calculate the integral and free the memory associated with the workspace w. //
    gsl_integration_qags(&F, start, stop, epsabs, epsrel, MAX_INTERVALS, w, &result, &error);
    gsl_integration_workspace_free(w);
    
    return result;
    
} // (integrate)


/**
 * @brief Apply additional limits to the code's time step.
 * @param pt5 Index of an active black hole particle.
 * @param time_step Current timestep.
 * @return Updated timestep in code units.
 **/
double timestep_limiter(int pt5, double time_step) {
    
    /* Check if either the viscous or the depletion timescales exist. If yes, then if allowed by the criteria, limit
     * the 'time_step'. */
    if (P[pt5].T_Viscous > 0.0 || P[pt5].T_Deplete > 0.0) {
        
        double dt_viscous = 1.0e-1 * P[pt5].T_Viscous;
        double dt_deplete = 1.0e-1 * P[pt5].T_Deplete;
        double kilo_year = 1.0e3 * SEC_PER_YEAR / All.UnitTime_in_s;
        
        if (dt_viscous < time_step) { time_step = dt_viscous; }
        if (dt_deplete < time_step) { time_step = dt_deplete; }
        if (time_step < kilo_year) { time_step = kilo_year; }
        
        /* Re-set the accretion disc's depletion timescale and viscous timescale to zero, so any future time-steps
         * without an accretion discs are not restricted by the last accretion episode. */
        P[pt5].T_Viscous = P[pt5].T_Deplete = 0.0;
    }
    
    return time_step;
    
} // timestep_limiter(...)


/**
 * @brief Sort the starting and stopping radii of all valid regimes.
 * @param pt5 Index of an active black hole particle.
 * @return Arrays of sorted starting and stopping radii.
 **/
void *sort_ranges_of_validity(int pt5, double starts[], double stops[]) {
    
    // Loop over all starting radii and combine them in the 'starts' array. //
    for (int bin = 0, all_bin = 0; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_starts[bin] > 0.0) { starts[all_bin] = P[pt5].GES_starts[bin], ++all_bin; }
        if (P[pt5].RES_starts[bin] > 0.0) { starts[all_bin] = P[pt5].RES_starts[bin], ++all_bin; }
        if (P[pt5].GFF_starts[bin] > 0.0) { starts[all_bin] = P[pt5].GFF_starts[bin], ++all_bin; }
    }
    
    // Loop over all stopping radii and combine them in the 'stops' array. //
    for (int bin = 0, all_bin = 0; bin < NUM_VALID; ++bin) {
        if (P[pt5].GES_stops[bin] > 0.0) { stops[all_bin] = P[pt5].GES_stops[bin], ++all_bin; }
        if (P[pt5].RES_stops[bin] > 0.0) { stops[all_bin] = P[pt5].RES_stops[bin], ++all_bin; }
        if (P[pt5].GFF_stops[bin] > 0.0) { stops[all_bin] = P[pt5].GFF_stops[bin], ++all_bin; }
    }
    
    /* Loop over all 'starts' and 'stops' excluding the fist starting and stopping radii, respectively and sort them in
    * ascending order. */
    for (int bin = 1; bin < 3 * NUM_VALID; ++bin) {
        double tmp_start = 0.0, tmp_stop = 0.0;
        if (starts[bin] > 0.0 && starts[bin - 1] > 0.0 && starts[bin] < starts[bin - 1]) {
            tmp_start = starts[bin];
            starts[bin] = starts[bin - 1];
            starts[bin - 1] = tmp_start;
        }
        if (stops[bin] > 0.0 && stops[bin - 1] > 0.0 && stops[bin] < stops[bin - 1]) {
            tmp_stop = stops[bin];
            stops[bin] = stops[bin - 1];
            stops[bin - 1] = tmp_stop;
        }
    }
    
} // sort_ranges_of_validity(...)


/**
 * @brief Calculate the Hubble parameter from the first Friedmann equation written as a function of the present day
 * density values and a scale factor of one.
 * @return Hubble parameter for a given scale factor.
 **/
double hubble_parameter() {
    
    // Check if co-moving integration is on. If yes, 'All.Time' is the scale factor. Else, return one. //
    if (All.ComovingIntegrationOn == 1) {
        double omega_k = 1.0 - All.Omega0 - All.OmegaLambda; // 'Omega0' represents the present day matter density.
        double H = All.Hubble
                   * sqrt(All.Omega0 * pow(All.Time, -3.0) + omega_k * pow(All.Time, -2.0) + All.OmegaLambda);
        
        // Check if the Hubble parameter is positive. If yes, return it, else terminate. //
        if (H > 0.0) { return H; }
        else {terminating("Invalid 'H'")}
    } else { return 1.0; }
    
} // hubble_parameter(...)