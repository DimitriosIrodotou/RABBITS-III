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
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../proto.h"
#include "../../../allvars.h"


/**
 * @brief Assign an initial, random spin to a black hole particle, and calculate its spin parameter.
 * @param pt5 Index of an active black hole particle.
 **/
void initialise_black_hole_spin(int pt5) {
    
    // Set and seed the desired random number generator (uniform in the [0.0,1.0) range). //
    gsl_rng *random_generator = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(random_generator, RNG_SEED);
    
    /* If 'spin_parameter' is higher than 0.998 (a physically motivated maximum value), then keep drawing random numbers
     * until either it is less than or equal to 0.998, or 'MAX_ITERATIONS' has been reached. */
    int iterations = 1;
    double spin_parameter;
    do {
        spin_parameter = gsl_rng_uniform(random_generator);
        if (iterations++ == MAX_ITERATIONS) {terminating("MAX_ITERATIONS has been reached"); }
    } while (spin_parameter > 0.998);
    
    // Calculate the black hole's spin vector magnitude in code units. //
    double spin = RA_G_CODE * pow(P[pt5].BH_Mass, 2.0) * spin_parameter / RA_C_CODE;
    
    // Assign random values to the black hole's spin vector components. //
    double phi = 2.0 * RA_PI * gsl_rng_uniform(random_generator);
    double theta = acos(2.0 * gsl_rng_uniform(random_generator) - 1.0);
    
    P[pt5].Spin[0] = spin * sin(theta) * cos(phi);
    P[pt5].Spin[1] = spin * sin(theta) * sin(phi);
    P[pt5].Spin[2] = spin * cos(theta);
    
    gsl_rng_free(random_generator); // Free the associated memory.
    
} // initialise_black_hole_spin(...)


/**
 * @brief Update the angular momenta of the black hole and of the accretion disc based on environmental processes and/or
 * based on their interactions.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param spin Dimensionless spin parameter.
 * @param accreting_mass Mass captured from ISM
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 **/
void exchange_angular_momenta(int pt5, double viscosity, double spin, double accreting_mass, const char *rotation) {
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment' and 'zero' when 'AD_Mass' is 0. //
    double alignment;
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else if (strcmp(rotation, "zero") == 0) { alignment = 0.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Calculate the relativistic specific and total angular momentum at the ISCO in code units. //
    double J_isco_code = isco_angmomentum(pt5, spin, rotation);
    double accreting_angmomentum = accreting_mass * J_isco_code;
    
    /* Calculate the 'Spin' and 'AD_AngMomentum' magnitudes and directions (i.e. unit vectors) and update their
     * magnitudes by 'accreting_angmomentum', which can be negative in 'counter-rotation' cases. */
    double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
    double ad_angmomentum_magnitude = sqrt(pow(P[pt5].AD_AngMomentum[0], 2.0) + pow(P[pt5].AD_AngMomentum[1], 2.0)
                                           + pow(P[pt5].AD_AngMomentum[2], 2.0));
    
    for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
        double spin_direction = P[pt5].Spin[direction] / spin_magnitude;
        double ad_angmomentum_direction = P[pt5].AD_AngMomentum[direction] / ad_angmomentum_magnitude;
        
        /* The magnitudes of the 'Spin' and 'AD_AngMomentum' vectors are updated whilst their directions remain the
         * same, since each vector component is changed by the same amount. */
        P[pt5].Spin[direction] = (spin_magnitude + accreting_angmomentum) * spin_direction;
        P[pt5].AD_AngMomentum[direction] = (ad_angmomentum_magnitude - accreting_angmomentum)
                                           * ad_angmomentum_direction;
    }
    
    // Convert properties to cgs units. //
    double bh_mass_cgs = convert_code_to_cgs_units("mass", P[pt5].BH_Mass);
    double ad_mass_cgs = convert_code_to_cgs_units("mass", P[pt5].AD_Mass);
    double accretion_rate_cgs;
    if (P[pt5].Form_Flag == 1) { accretion_rate_cgs = convert_code_to_cgs_units("mass_rate", P[pt5].AD_Mdot); }
    else { accretion_rate_cgs = convert_code_to_cgs_units("mass_rate", P[pt5].BH_Mdot); }
    
    /* Exchange angular momentum between the accretion disc and the black hole due to the Bardeen-Petterson effect
     * (i.e. the torque between the black hole and the accretion disc), whilst conserving the vector magnitudes. */
    (void) Bardeen_Petterson_effect(pt5, viscosity, bh_mass_cgs, ad_mass_cgs, accretion_rate_cgs, spin, rotation);
    
    /* Check if the total 'AD_AngMomentum' is lower than the angular momentum at the ISCO. If yes, move 'AD_Mass' to
     * 'BH_Mass' and set the accretion disc's properties to zero (i.e. dissolve the accretion disc). */
    ad_angmomentum_magnitude = sqrt(pow(P[pt5].AD_AngMomentum[0], 2.0) + pow(P[pt5].AD_AngMomentum[1], 2.0)
                                    + pow(P[pt5].AD_AngMomentum[2], 2.0));
    if (ad_angmomentum_magnitude < J_isco_code * P[pt5].AD_Mass) {
        P[pt5].BH_Mass += P[pt5].AD_Mass;
        (void) reset_accretion_disc(pt5);
    }
    
} // exchange_angular_momenta(...)


/**
 * @brief Update the angular momenta of the black hole and of the accretion disc based on environmental processes and/or
 * based on their interactions.
 * @param pt5 Index of an active black hole particle.
 * @param mode Select in RA_vars.c between 'chaotic' or 'coherent'.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 **/
void accrete_angular_momenta(int pt5, const char *mode, double spin) {
    
    /* Check if the 'mode' is chaotic. If yes, then the accretion disc's angular momentum vector gets assigned a random
     * direction. Else, the angular momentum vector from the captured gas particles gets added to the accretion disc's
     * angular momentum vector. */
    if ((strcmp(mode, "chaotic") == 0)) {
        
        // Set and seed the desired random number generator (uniform in the [0.0,1.0) range). //
        gsl_rng *random_generator = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(random_generator, RNG_SEED * All.NumCurrentTiStep);
        
        // Assign random values to the accretion disc's angular momentum vector components. //
        double phi = 2.0 * RA_PI * gsl_rng_uniform(random_generator);
        double theta = acos(2.0 * gsl_rng_uniform(random_generator) - 1.0);
        
        P[pt5].AD_AngMomentumDir[0] = sin(theta) * cos(phi);
        P[pt5].AD_AngMomentumDir[1] = sin(theta) * sin(phi);
        P[pt5].AD_AngMomentumDir[2] = cos(theta);
        
        gsl_rng_free(random_generator); // Free the associated memory.
        
    } else if (strcmp(mode, "coherent") == 0) {
        
        // Declare and initialise local variables with zeros. //
        double accreted_angmomentum_direction[3], dot_product = 0.0;
        memset(accreted_angmomentum_direction, 0.0, 3 * sizeof(double));
        
        // Calculate the 'Spin' and 'AD_AngMomentum' directions (i.e. unit vectors) and their 'dot_product'. //
        double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
        double accreted_angmomentum_magnitude = sqrt(pow(P[pt5].b0.dBH_accreted_angmomentum[0], 2.0)
                                                     + pow(P[pt5].b0.dBH_accreted_angmomentum[1], 2.0)
                                                     + pow(P[pt5].b0.dBH_accreted_angmomentum[2], 2.0));
        
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            double spin_direction = P[pt5].Spin[direction] / spin_magnitude;
            accreted_angmomentum_direction[direction] = P[pt5].b0.dBH_accreted_angmomentum[direction]
                                                        / accreted_angmomentum_magnitude;
            dot_product += spin_direction * accreted_angmomentum_direction[direction];
        }
        
        // Positive/negative 'dot_product' corresponds to positive/negative 'alignment'. //
        double alignment;
        if (dot_product >= 0.0) { alignment = 1.0; }
        else if (dot_product < 0.0) { alignment = -1.0; }
        else {terminating("Invalid 'rotation'"); }
        
        // Calculate the relativistic corrections 'c_correction' and 'f_correction' at 'x_outermost'. //
        double mass_cgs = convert_code_to_cgs_units("mass", P[pt5].BH_Mass + P[pt5].AD_Mass);
        double x_outermost = sqrt(r_outermost(pt5, mass_cgs) / r_grav(mass_cgs)); // Dimensionless outermost radius.
        double c_correction = 1.0 - 3.0 * pow(x_outermost, -2.0) + alignment * 2.0 * spin * pow(x_outermost, -3.0);
        double f_correction = alignment * (1.0 - alignment * 2.0 * spin * pow(x_outermost, -3.0)
                                           + pow(spin, 2.0) * pow(x_outermost, -4.0));
        
        // Calculate the relativistic corrected angular momentum to be added to the accretion disc. //
        double radius_code = convert_cgs_to_code_units("length", r_outermost(pt5, mass_cgs));
        double angmomentum_magnitude = sqrt(RA_G_CODE * (P[pt5].BH_Mass + P[pt5].AD_Mass) * radius_code)
                                       * pow(c_correction, -1.0 / 2.0) * f_correction * P[pt5].b4.dBH_accreted_Mass;
        
        /* The absolute value of 'angmomentum_magnitude' should be used because 'f_correction' is negative when
         * 'dot_product' is negative, which flips the sign of each 'AD_AngMomentum' component and makes it always
         * align with 'Spin' (i.e. theta < pi/2). */
        P[pt5].AD_AngMomentum[0] += fabs(angmomentum_magnitude) * accreted_angmomentum_direction[0];
        P[pt5].AD_AngMomentum[1] += fabs(angmomentum_magnitude) * accreted_angmomentum_direction[1];
        P[pt5].AD_AngMomentum[2] += fabs(angmomentum_magnitude) * accreted_angmomentum_direction[2];
        
    } else {terminating("Invalid 'mode' for 'add'"); }
    
} // accrete_angular_momenta(...)


/**
 * @brief Assess the alignment/counter-alignment between the accretion disc and the black hole.
 * @param pt5 Index of an active black hole particle.
 * @return Either 'co-rotation' or 'counter-rotation'.
 **/
const char *assess_alignment(int pt5) {
    
    double dot_product = 0.0; // Declare local variables.
    
    // Calculate the 'Spin' and 'AD_AngMomentum' directions (i.e. unit vectors) and their 'dot_product'. //
    double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
    double ad_angmomentum_magnitude = sqrt(pow(P[pt5].AD_AngMomentum[0], 2.0) + pow(P[pt5].AD_AngMomentum[1], 2.0)
                                           + pow(P[pt5].AD_AngMomentum[2], 2.0));
    
    /* Check if the accretion disc already has angular momentum (i.e. it exists). If yes, use the King et al. 2005
     * (2005MNRAS.363...49K) criterion. Else, use the 'dot_product' between the black hole spin and the accretion disc
     * angular momenta. */
    if (ad_angmomentum_magnitude > 0.0) {
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            double spin_direction = P[pt5].Spin[direction] / spin_magnitude;
            double ad_angmomentum_direction = P[pt5].AD_AngMomentum[direction] / ad_angmomentum_magnitude;
            dot_product += spin_direction * ad_angmomentum_direction;
        }
        
        // Assess the alignment/counter-alignment. //
        double spin_ratio = -ad_angmomentum_magnitude / (2.0 * spin_magnitude);
        
        if (dot_product >= spin_ratio) { return "co-rotation"; }
        else if (dot_product < spin_ratio) { return "counter-rotation"; }
        else {terminating("Invalid 'rotation'"); }
        
    } else {
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            double spin_direction = P[pt5].Spin[direction] / spin_magnitude;
            dot_product += spin_direction * P[pt5].AD_AngMomentumDir[direction];
        }
        
        // Positive/negative 'dot_product' corresponds to positive/negative 'alignment'. //
        if (dot_product >= 0.0) { return "co-rotation"; }
        else if (dot_product < 0.0) { return "counter-rotation"; }
        else {terminating("Invalid 'rotation'"); }
    }
    
} // assess_alignment(...)


/**
 * @brief Calculate the dimensionless spin parameter from a black hole's spin vector.
 * @param pt5 Index of an active black hole particle.
 * @return Dimensionless 'spin_parameter'.
 **/
double calculate_spin_parameter(int pt5) {
    
    // Calculate the black hole's spin vector magnitude. //
    double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
    
    // Check if the 'spin_parameter' is between 0 and 0.998. If yes, then return it. Else, terminate. //
    double spin_parameter = spin_magnitude * RA_C_CODE / (RA_G_CODE * pow(P[pt5].BH_Mass, 2.0));
    if (spin_parameter >= 0.0 && spin_parameter <= 0.998) { return spin_parameter; }
    else {terminating("Invalid 'spin_parameter'"); }
    
} // calculate_spin_parameter(...)


/**
 * @brief Calculate the evolution of the black hole spin due to the Bardeen_Petterson_effect.
 * @param pt5 Index of an active black hole particle.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param ad_mass Accretion disc mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 **/
void Bardeen_Petterson_effect(int pt5, double viscosity, double bh_mass, double ad_mass, double accretion_rate,
                              double spin, const char *rotation) {
    
    // Calculate the 'Spin' and 'AD_AngMomentum' magnitudes and directions (i.e. unit vectors). //
    double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
    double ad_angmomentum_magnitude = sqrt(pow(P[pt5].AD_AngMomentum[0], 2.0) + pow(P[pt5].AD_AngMomentum[1], 2.0)
                                           + pow(P[pt5].AD_AngMomentum[2], 2.0));
    
    double spin_direction[3], ad_angmomentum_direction[3];
    for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
        spin_direction[direction] = P[pt5].Spin[direction] / spin_magnitude;
        ad_angmomentum_direction[direction] = P[pt5].AD_AngMomentum[direction] / ad_angmomentum_magnitude;
    }
    
    // Calculate the cross product between the unit vectors of 'AD_AngMomentum' and 'Spin' //
    double BP_spin_direction[3];
    BP_spin_direction[0] = ad_angmomentum_direction[1] * spin_direction[2]
                           - ad_angmomentum_direction[2] * spin_direction[1];
    BP_spin_direction[1] = -(ad_angmomentum_direction[0] * spin_direction[2]
                             - ad_angmomentum_direction[2] * spin_direction[0]);
    BP_spin_direction[2] = ad_angmomentum_direction[0] * spin_direction[1]
                           - ad_angmomentum_direction[1] * spin_direction[0];
    
    // Convert 'BP_spin_direction[]' into a unit vector's directions (i.e. divide them by the vector magnitude). //
    double BP_spin_magnitude = sqrt(pow(BP_spin_direction[0], 2.0) + pow(BP_spin_direction[1], 2.0)
                                    + pow(BP_spin_direction[2], 2.0));
    
    /* Check if full alignment has not been accomplished (i.e. the 'BP_spin_magnitude' is higher than a fraction of the
     * 'ad_angmomentum_magnitude') */
    if (BP_spin_magnitude > ad_angmomentum_magnitude / pow(UNIT_TEST_TOLERANCE, 2)) {
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            BP_spin_direction[direction] /= BP_spin_magnitude;
        }
        
        // Calculate the black hole's time-step (i.e. time since it was last active) and convert to cgs units. //
        double dt = (P[pt5].TimeBin ? (1 << P[pt5].TimeBin) : 0.0) * All.Timebase_interval / hubble_parameter();
        double dt_cgs = convert_code_to_cgs_units("time", dt);
        
        // Calculate the warp radius of an accretion disc that's misaligned with respect to the black hole's equator.
        (void) r_warp(pt5, VISCOSITY, bh_mass, accretion_rate, spin, rotation);
        
        /* Calculate the magnitude of the Bardeen-Petterson torque (i.e. rate of change of angular momentum) and convert to
         * code units. */
        double integrated_spin = integrated_profiles(pt5, "spin", viscosity, bh_mass, accretion_rate, spin, rotation);
        double BP_spin_magnitude_cgs = 4.0 * RA_PI * RA_G_CGS / pow(RA_C_CGS, 2) * integrated_spin * dt_cgs;
        double BP_spin_magnitude_code = convert_cgs_to_code_units("angmomentum", BP_spin_magnitude_cgs);
        
        /* The absolute value of 'BP_spin_magnitude_code' should be used because 'L_tilde' in 'spin_integrand' is negative
         * when 'alignment' is negative, due to 'correction[4]' (i.e. calligraphic F relativistic correction). */
        double new_spin_magnitude = sqrt(pow(fabs(BP_spin_magnitude_code) / spin_magnitude * BP_spin_direction[0], 2.0)
                                         +
                                         pow(fabs(BP_spin_magnitude_code) / spin_magnitude * BP_spin_direction[1], 2.0)
                                         + pow(fabs(BP_spin_magnitude_code) / spin_magnitude * BP_spin_direction[2],
                                               2.0));
        
        double new_spin_direction[3];
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            new_spin_direction[direction] = fabs(BP_spin_magnitude_code) / spin_magnitude
                                            * BP_spin_direction[direction] / new_spin_magnitude;
        }
        
        double updated_magnitude = sqrt(pow(spin_direction[0] + new_spin_direction[0], 2.0)
                                        + pow(spin_direction[1] + new_spin_direction[1], 2.0)
                                        + pow(spin_direction[2] + new_spin_direction[2], 2.0));
        
        for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
            P[pt5].Spin[direction] = spin_magnitude * (spin_direction[direction] + new_spin_direction[direction])
                                     / updated_magnitude;
            
            // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
            double alignment;
            if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
            else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
            else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
            
            /* If the accretion disc is counter-rotating wrt the black hole's spin, then the accretion disc will
             * counter-align.*/
            P[pt5].AD_AngMomentum[direction] = ad_angmomentum_magnitude * alignment
                                               * (spin_direction[direction] + new_spin_direction[direction])
                                               / updated_magnitude;
        }
    }
} // Bardeen_Petterson_effect(...)


/**
 * @brief Calculate the relativistic specific angular momentum at the ISCO.
 * @param pt5 Index of an active black hole particle.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return  Angular momentum integrand in code units.
 **/
double isco_angmomentum(int pt5, double spin, const char *rotation) {
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    double alignment;
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Calculate the spin-depended corrections from Bardeen+72 (1972ApJ...178..347B). //
    double z_1 = 1 + pow(1.0 - pow(spin, 2.0), 1.0 / 3.0)
                     * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double z_2 = sqrt(3.0 * pow(spin, 2.0) + pow(z_1, 2.0));
    double z = 3.0 + z_2 - alignment * sqrt((3.0 - z_1) * (3.0 + z_1 + 2.0 * z_2));
    
    return alignment * RA_G_CODE * P[pt5].BH_Mass / RA_C_CODE
           * (pow(z, 2.0) - alignment * 2.0 * spin * sqrt(z) + pow(spin, 2.0))
           / (z * sqrt(z - 3.0 + alignment * 2.0 * spin * pow(z, -1.0 / 2.0)));
    
} // isco_angmomentum(...)