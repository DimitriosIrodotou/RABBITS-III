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
 * @brief Main function to call the functions defined in RA_proto.h and declared in the other 'RA_*.c' scripts.
 **/
int relativistic_accretion(void) {
    
    /* Check if this is the first time step. If yes, then initialise black hole and accretion disc properties. Else,
     * evolve black holes and their accretion discs. */
    if (All.NumCurrentTiStep == 0) { (void) initialise_setup(); }
    else {
        // All CPUs search for gas particles that meet the capture criteria and store their properties. //
        (void) capture_SPH_particles();
        
        // Loop over all active particles and check if an active particle is a black hole with positive mass. //
        for (int pt5 = FirstActiveParticle; pt5 >= 0; pt5 = NextActiveParticle[pt5]) {
            if (P[pt5].Type == 5 && P[pt5].BH_Mass > 0.0) {
                
                // If a binary has formed, exit this function and move to the 2211.11788 binary accretion model. //
                if (P[pt5].isBinary) { return 1; }
                
                printf("----------RELATIVISTIC_ACCRETION---------- \n"), fflush(stdout);
                
                // Calculate the 'spin_parameter' due to initialisation or previous accretion episode. //
                double spin_parameter = calculate_spin_parameter(pt5);
                
                /* Check if mass has been accreted in 'capture_SPH_particles'. If yes, update black hole and accretion
                 * disc properties. */
                if (P[pt5].b4.dBH_accreted_Mass > 0.0) {
                    // Update the black hole's velocity vector components first and then its dynamical mass. //
                    for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
                        P[pt5].Vel[direction] = (P[pt5].Vel[direction] * P[pt5].Mass
                                                 + P[pt5].b6.dBH_accreted_momentum[direction])
                                                / (P[pt5].Mass + P[pt5].b4.dBH_accreted_Mass);
                    }
                    P[pt5].Mass += P[pt5].b4.dBH_accreted_Mass;
                    
                    // Check if the 'AD_Mass' and 'Form_Flag' are zero. If yes, then assign an accretion rate. //
                    if (P[pt5].AD_Mass == 0.0 && P[pt5].Form_Flag == 0) {
                        P[pt5].Form_Flag = 1;
                        P[pt5].AD_Mdot = P[pt5].b4.dBH_accreted_Mass / bondi_hoyle_timescale(pt5);
                    }
                    
                    // Update the accretion disc's angular momentum first and then its mass. //
                    (void) accrete_angular_momenta(pt5, MODE, spin_parameter);
                    P[pt5].AD_Mass += P[pt5].b4.dBH_accreted_Mass;
                }
                
                /* Check if a black hole and an accretion disc exist. If yes, calculate the appropriate accretion rate
                 * and use the structure of the accretion disc to evolve black hole properties. */
                if (P[pt5].BH_Mass > 0.0 && P[pt5].AD_Mass > 0.0) {
                    (void) log_details(pt5, "in"); // Log information before updating any properties.
                    
                    // Calculate the black hole's time-step (i.e. time since it was last active) in code units. //
                    double dt = (P[pt5].TimeBin ? (1 << P[pt5].TimeBin) : 0.0)
                                * All.Timebase_interval / hubble_parameter();
                    
                    /* Check if an accretion disc has been formed in this time-step. If yes, calculate the black hole's
                     * accretion rate based on the SPH mass captured over the Bondi-Hoyle timescale. Else, use the
                     * accretion rate calculated by evolving an existing accretion disc and increase the accretion
                     * disc's lifetime by a time-step. */
                    double accretion_rate_cgs;
                    if (P[pt5].Form_Flag == 1) {
                        accretion_rate_cgs = convert_code_to_cgs_units("mass_rate", P[pt5].AD_Mdot);
                    } else {
                        P[pt5].AD_Lifetime += convert_code_to_cgs_units("time", dt);
                        accretion_rate_cgs = convert_code_to_cgs_units("mass_rate", P[pt5].BH_Mdot);
                    }
                    
                    // Assess the alignment/counter-alignment between the accretion disc and the black hole. //
                    const char *rotation = assess_alignment(pt5);
                    
                    // Calculate the radial extent of the accretion disc. //
                    double bh_mass_cgs = convert_code_to_cgs_units("mass", P[pt5].BH_Mass);
                    (void) accretion_disc_extent(pt5, VISCOSITY, bh_mass_cgs, accretion_rate_cgs, spin_parameter,
                                                 rotation);
                    
                    // Calculate the range of validity of different regimes, which combined form the accretion disc. //
                    (void) ranges_of_validity(pt5, VISCOSITY, bh_mass_cgs, accretion_rate_cgs, spin_parameter,
                                              rotation);
                    
                    // Calculate the black hole's accretion rate. //
                    double ad_mass_cgs = convert_code_to_cgs_units("mass", P[pt5].AD_Mass);
                    double mdot = accretion_rate(pt5, VISCOSITY, bh_mass_cgs, ad_mass_cgs, accretion_rate_cgs,
                                                 spin_parameter, rotation);
                    
                    // Convert the accretion rate back into code units and update 'BH_Mdot'. //
                    P[pt5].BH_Mdot = convert_cgs_to_code_units("mass_rate", mdot);
                    
                    // Calculate the accretion disc's depletion and viscous timescale for the time-step limiter. //
                    P[pt5].T_Deplete = P[pt5].AD_Mass / P[pt5].BH_Mdot;
                    double t_viscous_cgs = integrated_profiles(pt5, "t_viscous", VISCOSITY, bh_mass_cgs, mdot,
                                                               spin_parameter, rotation);
                    P[pt5].T_Viscous = convert_cgs_to_code_units("time", t_viscous_cgs);
                    
                    /* Check if the mass to be accreted by the black hole is more than the accretion disc's mass. If
                     * yes, the black hole acquires the whole accretion disc. Else, a fraction of the accretion disc's
                     * mass and angular momentum is transferred to the black hole. */
                    double accreting_mass = (1.0 - P[pt5].Epsilon_r) * P[pt5].BH_Mdot * dt;
                    if (accreting_mass >= P[pt5].AD_Mass) {
                        // Transfer all the accretion disc's angular momentum and mass to the black hole. //
                        (void) exchange_angular_momenta(pt5, VISCOSITY, spin_parameter, P[pt5].AD_Mass, rotation);
                        P[pt5].BH_Mass += P[pt5].AD_Mass;
                        
                        /* Update the black hole's accretion rate and set the 'Drain_Flag' to one (i.e. dissolve the
                         * accretion disc at the end of the time-step). */
                        P[pt5].BH_Mdot = P[pt5].AD_Mass / ((1.0 - P[pt5].Epsilon_r) * dt);
                        P[pt5].Drain_Flag = 1;
                    } else {
                        // Exchange angular momentum and mass between the black hole and the accretion disc. //
                        (void) exchange_angular_momenta(pt5, VISCOSITY, spin_parameter, accreting_mass, rotation);
                        P[pt5].BH_Mass += accreting_mass, P[pt5].AD_Mass -= accreting_mass;
                    }
                    
                    (void) log_details(pt5, "out"); // Log information after updating the properties.
                }
                
                // Re-set the 'Form_Flag', 'accreted_mass', 'accreted_momentum', and 'accreted_angmomentum' to zero. //
                P[pt5].Form_Flag = 0;
                P[pt5].b4.dBH_accreted_Mass = 0.0;
                for (int direction = 0; direction < 3; ++direction) { // Loop over x, y, and z direction.
                    P[pt5].b6.dBH_accreted_momentum[direction] = 0.0;
                    P[pt5].b0.dBH_accreted_angmomentum[direction] = 0.0;
                }
                printf("----------RELATIVISTIC_ACCRETION---------- \n"), fflush(stdout);
            }
        }
        
        // Calculate the released thermal energy and distribute it to the neighbouring gas particles. //
        thermal_feedback();
        
        /* Loop over all active particles and check if an active particle is a black
         * hole with positive mass whose accretion disc has been dissolved. If yes, then re-set the accretion disc's
         * properties. */
        for (int pt5 = FirstActiveParticle; pt5 >= 0; pt5 = NextActiveParticle[pt5]) {
            if (P[pt5].Type == 5 && P[pt5].BH_Mass > 0.0 && P[pt5].Drain_Flag == 1) {
                (void) reset_accretion_disc(pt5);
            }
        }
    }
    return 0;
} // relativistic_accretion(...)