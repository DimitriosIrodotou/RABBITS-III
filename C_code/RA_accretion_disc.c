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
#include <gsl/gsl_integration.h>

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../allvars.h"


/**
 * @brief Calculate the radially integrated profiles.
 * @param pt5 Index of an active black hole particle.
 * @param mode Select between 'spin', 'angmomentum', or 't_viscous'.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Black hole accretion rate in gram second^-1.
 **/
double integrated_profiles(int pt5, const char *mode, double viscosity, double bh_mass, double accretion_rate,
                           double spin, const char *rotation) {
    
    /* In the following three blocks of code, check for every regime if it is anywhere valid. If yes, then integrate
     * the corresponding 'mode' profile. */
    
    double GES = 0.0;
    // Check if the 'GES' regime is valid somewhere. //
    if (P[pt5].GES_sum > 0.0) {
        // Loop over 'NUM_VALID' radial bins. //
        for (int bin = 0; bin < NUM_VALID; ++bin) {
            /* Check if there are valid starting and stopping radii. If yes, then integrate corresponding 'mode'
             * profile of the 'GES' regime. */
            if (P[pt5].GES_starts[bin] > 0.0 && P[pt5].GES_stops[bin] > 0.0) {
                /* When integrating for the Bardeen-Petterson effect (i.e. mode = 'spin'), stop at 'R_Warp' if it is
                 * between 'GES_starts' and 'GES_stops' */
                if (strcmp(mode, "spin") == 0 && P[pt5].R_Warp > P[pt5].GES_starts[bin]
                    && P[pt5].R_Warp < P[pt5].GES_stops[bin]) {
                    GES += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].R_Warp,
                                     P[pt5].GES_stops[bin], rotation, "GES");
                } else {
                    GES += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].GES_starts[bin],
                                     P[pt5].GES_stops[bin], rotation, "GES");
                }
            }
        }
    }
    
    double RES = 0.0;
    // Check if the 'RES' regime is valid somewhere. //
    if (P[pt5].RES_sum > 0.0) {
        // Loop over 'NUM_VALID' radial bins. //
        for (int bin = 0; bin < NUM_VALID; ++bin) {
            /* Check if there are valid starting and stopping radii. If yes, then integrate corresponding 'mode'
             * profile of the 'RES' regime. */
            if (P[pt5].RES_starts[bin] > 0.0 && P[pt5].RES_stops[bin] > 0.0) {
                /* When integrating for the Bardeen-Petterson effect (i.e. mode = 'spin'), stop at 'R_Warp' if it is
                 * between 'RES_starts' and 'RES_stops' */
                if (strcmp(mode, "spin") == 0 && P[pt5].R_Warp > P[pt5].RES_starts[bin]
                    && P[pt5].R_Warp < P[pt5].RES_stops[bin]) {
                    RES += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].R_Warp,
                                     P[pt5].RES_stops[bin], rotation, "RES");
                } else {
                    RES += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].RES_starts[bin],
                                     P[pt5].RES_stops[bin], rotation, "RES");
                }
            }
        }
    }
    
    double GFF = 0.0;
    // Check if the 'GFF' regime is valid somewhere. //
    if (P[pt5].GFF_sum > 0.0) {
        // Loop over 'NUM_VALID' radial bins. //
        for (int bin = 0; bin < NUM_VALID; ++bin) {
            /* Check if there are valid starting and stopping radii. If yes, then integrate corresponding 'mode'
             * profile of the 'GFF' regime. */
            if (P[pt5].GFF_starts[bin] > 0.0 && P[pt5].GFF_stops[bin] > 0.0) {
                /* When integrating for the Bardeen-Petterson effect (i.e. mode = 'spin'), stop at 'R_Warp' if it is
                 * between 'GFF_starts' and 'GFF_stops' */
                if (strcmp(mode, "spin") == 0 && P[pt5].R_Warp >= P[pt5].GFF_starts[bin]
                    && P[pt5].R_Warp <= P[pt5].GFF_stops[bin]) {
                    GFF += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].R_Warp,
                                     P[pt5].GFF_stops[bin], rotation, "GFF");
                } else {
                    GFF += integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, P[pt5].GFF_starts[bin],
                                     P[pt5].GFF_stops[bin], rotation, "GFF");
                }
            }
        }
    }
    
    /* Integrate corresponding 'mode' profile of the intra-ISCO regime from bin right after 'r_photon' to the bin just
     * before 'r_isco' (i.e. to 'radial_bins_ISCO[NUM_BINS - 2]') */
    double radial_bins_ISCO[NUM_BINS], intraISCO;
    (void) logspace(log10(r_photon(bh_mass, spin, rotation)), log10(r_isco(bh_mass, spin, rotation)), NUM_BINS,
                    radial_bins_ISCO);
    
    intraISCO = integrate(pt5, mode, viscosity, bh_mass, accretion_rate, spin, radial_bins_ISCO[1],
                          radial_bins_ISCO[NUM_BINS - 2], rotation, "ISCO");
    
    double total = GES + RES + GFF + intraISCO; // Integrated profile.
    
    if (strcmp(mode, "angmomentum") == 0) {
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            // Update the 'AD_AngMomentum' vector in code units. //
            P[pt5].AD_AngMomentum[direction] += convert_cgs_to_code_units("angmomentum", total)
                                                * P[pt5].AD_AngMomentumDir[direction];
        }
    }
    
    // Sanity checks, more info in the print statement. //
    if (GES < 0.0 && RES > 0.0 && GFF > 0.0 && intraISCO > 0.0) {terminating("GES_starts > GES_stops in 'integral'"); }
    if (RES < 0.0 && GES > 0.0 && GFF > 0.0 && intraISCO > 0.0) {terminating("RES_starts > RES_stops in 'integral'"); }
    if (GFF < 0.0 && GES > 0.0 && RES > 0.0 && intraISCO > 0.0) {terminating("GFF_starts > GFF_stops in 'integral'"); }
    if (intraISCO < 0.0 && GES > 0.0 && RES > 0.0 && GFF > 0.0) {
        terminating("r_photon() > radial_bins_ISCO[NUM_BINS - 2] in 'integral'");
    }
    
    return total;
    
} // integrated_profiles(...)


/**
 * @brief Calculate the black hole accretion rate.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param ad_mass Accretion disc mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Black hole accretion rate in gram second^-1.
 **/
double accretion_rate(int pt5, double viscosity, double bh_mass, double ad_mass, double accretion_rate, double spin,
                      const char *rotation) {
    
    // Calculate the accretion disc's angular momentum 'integrated_angmomentum' and the black hole spin 'S_BH'.
    double integrated_angmomentum = integrated_profiles(pt5, "angmomentum", viscosity, bh_mass, accretion_rate, spin,
                                                        rotation);
    double S_BH = RA_G_CGS * pow(bh_mass, 2.0) * spin / RA_C_CGS;
    
    // Calculate the total radiative efficiency 'Epsilon_r' based on which 'regime' is valid at the ISCO. //
    P[pt5].Epsilon_spin = radiative_efficiency_spin(bh_mass, spin, rotation);
    P[pt5].Epsilon_r = P[pt5].Epsilon_spin + radiative_efficiency_isco(viscosity, bh_mass, accretion_rate, spin,
                                                                       rotation, P[pt5].Regime_ISCO);
    
    // Calculate the black hole accretion rate as a fraction of the Eddington accretion rate. //
    double f_edd = 0.76 * (P[pt5].Epsilon_r / 0.1) * pow(ad_mass / (1.0e4 * RA_SOLAR_MASS_CGS), 5.0)
                   * pow(bh_mass / (1.0e6 * RA_SOLAR_MASS_CGS), -47.0 / 7.0)
                   * pow((spin * fabs(integrated_angmomentum) / S_BH) / 3.0, -25.0 / 7.0);
    
    double mdot_edd = (4.0 * RA_PI * RA_G_CGS * RA_PROTON_MASS_CGS / (RA_THOMPSON_CGS * RA_C_CGS))
                      * (bh_mass / P[pt5].Epsilon_r);
    double min = (f_edd < 1.0) ? f_edd : 1.0, max = (0.01 > min) ? 0.01 : min;
    
    // Convert the Eddington accretion rate to code units to export it in the 'log_details'. //
    P[pt5].Edd_Mdot = convert_cgs_to_code_units("mass_rate", mdot_edd);
    
    return mdot_edd * max;
    
} // accretion_rate(...)


/**
 * @brief Calculate the angular momentum integrand for different regimes.
 * @param r Radial coordinate.
 * @param pt5 Index of an active black hole particle.
 * @param passed_parameters Structure of parameters.
 * @return Angular momentum integrand in centimeter gram second^-1.
 **/
double angmomentum_integrand(double r, void *passed_parameters) {
    
    // Declare a structure for 'local_parameters' and a local variable. //
    struct integrand_parameters *local_parameters;
    local_parameters = (struct integrand_parameters *) passed_parameters;
    int pt5 = local_parameters->pt5;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, local_parameters->viscosity, local_parameters->bh_mass,
                                local_parameters->accretion_rate,
                                local_parameters->spin, local_parameters->rotation, local_parameters->regime,
                                corrections);
    
    // Calculate the angular momentum integrand. //
    if (strcmp(local_parameters->regime, "GES") == 0 || strcmp(local_parameters->regime, "RES") == 0 ||
        strcmp(local_parameters->regime, "GFF") == 0 || strcmp(local_parameters->regime, "ISCO") == 0) {
        double L_tilde = sqrt(RA_G_CGS * local_parameters->bh_mass * r) * pow(corrections[0], -1.0 / 2.0)
                         * corrections[4];
        return 2.0 * RA_PI * r * L_tilde * surface_densities(pt5, r, local_parameters->viscosity,
                                                             local_parameters->bh_mass,
                                                             local_parameters->accretion_rate, local_parameters->spin,
                                                             local_parameters->rotation, local_parameters->regime);
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // angmomentum_integrand(...)


/**
 * @brief Calculate the viscous timescale integrand for different regimes.
 * @param r Radial coordinate.
 * @param pt5 Index of an active black hole particle.
 * @param passed_parameters Structure of parameters.
 * @return Inverse radial velocity in centimeter^-1 second.
 **/
double t_viscous_integrand(double r, void *passed_parameters) {
    
    // Declare a structure for 'local_parameters' and a local variable. //
    struct integrand_parameters *local_parameters;
    local_parameters = (struct integrand_parameters *) passed_parameters;
    int pt5 = local_parameters->pt5;
    
    // Calculate the viscous timescale integrand. //
    if (strcmp(local_parameters->regime, "GES") == 0 || strcmp(local_parameters->regime, "RES") == 0 ||
        strcmp(local_parameters->regime, "GFF") == 0 || strcmp(local_parameters->regime, "ISCO") == 0) {
        return 1.0 / radial_velocities(pt5, r, local_parameters->viscosity, local_parameters->bh_mass,
                                       local_parameters->accretion_rate, local_parameters->spin,
                                       local_parameters->rotation, local_parameters->regime);
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // t_viscous_integrand(...)


/**
 * @brief Calculate the spin integrand for different regimes.
 * @param r Radial coordinate.
 * @param passed_parameters Structure of parameters.
 * @return Spin integrand in gram^2 second^-2.
 **/
double spin_integrand(double r, void *passed_parameters) {
    
    // Declare a structure for 'local_parameters' and a local variable. //
    struct integrand_parameters *local_parameters;
    local_parameters = (struct integrand_parameters *) passed_parameters;
    int pt5 = local_parameters->pt5;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, local_parameters->viscosity, local_parameters->bh_mass,
                                local_parameters->accretion_rate, local_parameters->spin, local_parameters->rotation,
                                local_parameters->regime, corrections);
    
    // Calculate the spin integrand. //
    if (strcmp(local_parameters->regime, "GES") == 0 || strcmp(local_parameters->regime, "RES") == 0 ||
        strcmp(local_parameters->regime, "GFF") == 0 || strcmp(local_parameters->regime, "ISCO") == 0) {
        double L_tilde = sqrt(RA_G_CGS * local_parameters->bh_mass * r) * pow(corrections[0], -1.0 / 2.0)
                         * corrections[4];
        
        // Calculate the sine of the angle between 'AD_AngMomentum' and 'Spin' by using the cross product method.
        double cross_0 = P[pt5].AD_AngMomentum[1] * P[pt5].Spin[2] - P[pt5].AD_AngMomentum[2] * P[pt5].Spin[1];
        double cross_1 = -(P[pt5].AD_AngMomentum[0] * P[pt5].Spin[2] - P[pt5].AD_AngMomentum[2] * P[pt5].Spin[0]);
        double cross_2 = P[pt5].AD_AngMomentum[0] * P[pt5].Spin[1] - P[pt5].AD_AngMomentum[1] * P[pt5].Spin[0];
        
        double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
        double ad_angmomentum_magnitude = sqrt(pow(P[pt5].AD_AngMomentum[0], 2.0) + pow(P[pt5].AD_AngMomentum[1], 2.0)
                                               + pow(P[pt5].AD_AngMomentum[2], 2.0));
        double cross_product_magnitude = sqrt(pow(cross_0, 2) + pow(cross_1, 2) + pow(cross_2, 2));
        double sin_theta = cross_product_magnitude / (ad_angmomentum_magnitude * spin_magnitude);
        
        double spin_magnitude_cgs = convert_code_to_cgs_units("angmomentum", spin_magnitude);
        return surface_densities(pt5, r, local_parameters->viscosity, local_parameters->bh_mass,
                                 local_parameters->accretion_rate, local_parameters->spin,
                                 local_parameters->rotation, local_parameters->regime)
               * L_tilde * spin_magnitude_cgs * sin_theta / pow(r, 2);
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // spin_integrand(...)