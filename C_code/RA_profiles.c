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

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../allvars.h"


/**
 * @brief Calculate the relativistic correction functions for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @param corrections Array to store correction functions.
 **/
void *correction_functions(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                           const char *rotation, const char *regime, double corrections[]) {
    
    // Declare local variables. //
    double h_0, alignment;
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double x_isco = sqrt(r_isco(bh_mass, spin, rotation) / r_grav(bh_mass)); // Dimensionless radius of the ISCO.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Correction functions and values at the ISCO (see in-line comments for their names). //
    corrections[0] = 1.0 - 3.0 * pow(x, -2.0) + alignment * 2.0 * spin * pow(x, -3.0); // C
    corrections[1] = 1.0 - 3.0 * pow(x_isco, -2.0) + alignment * 2.0 * spin * pow(x_isco, -3.0); // C_isco
    corrections[2] = 1.0 - 2.0 * pow(x, -2.0) + pow(spin, 2.0) * pow(x, -4.0); // D
    corrections[3] = 1.0 - 2.0 * pow(x_isco, -2.0) + pow(spin, 2.0) * pow(x_isco, -4.0); // D_isco
    corrections[4] = alignment * (1.0 - alignment * 2.0 * spin * pow(x, -3.0) + pow(spin, 2.0) * pow(x, -4.0)); // F
    corrections[5] = alignment * (1.0 - alignment * 2.0 * spin * pow(x_isco, -3.0) + pow(spin, 2.0)
                                                                                     * pow(x_isco, -4.0)); // F_isco
    corrections[6] = 1.0 - 2.0 * pow(x, -2.0) + alignment * spin * pow(x, -3.0); // G
    corrections[7] = 1.0 - 2.0 * pow(x_isco, -2.0) + alignment * spin * pow(x_isco, -3.0); // G_isco
    corrections[8] = pow(corrections[0], -1.0) * pow(corrections[4], 2.0)
                     - pow(spin, 2.0) * pow(x, -2.0)
                       * (pow(corrections[0], -1.0 / 2.0) * corrections[6] - 1.0); // R
    corrections[9] = pow(corrections[1], -1.0) * pow(corrections[5], 2.0)
                     - pow(spin, 2.0) * pow(x_isco, -2.0)
                       * (pow(corrections[1], -1.0 / 2.0) * corrections[7] - 1.0); // R_isco
    
    if (strcmp(regime, "GES") == 0) {
        h_0 = 0.00184101952815808 * pow(viscosity, 1.0 / 8.0) * pow(mass_dmnsless, -3.0 / 8.0)
              * pow(accretion_rate_dmnsless, 1.0 / 4.0) * pow(x_isco, 1.0 / 8.0) * pow(corrections[1], -1.0 / 8.0)
              * pow(corrections[9], -1.0 / 2.0);
        
    } else if (strcmp(regime, "RES") == 0) {
        h_0 = 0.002;
        
    } else if (strcmp(regime, "GFF") == 0) {
        h_0 = 0.0013166295323149043 * pow(viscosity, 1.0 / 17.0) * pow(mass_dmnsless, -5.0 / 17.0)
              * pow(accretion_rate_dmnsless, 3.0 / 17.0) * pow(x_isco, 5.0 / 17.0) * pow(corrections[1], -1.0 / 17.0)
              * pow(corrections[3], -1.0 / 34.0) * pow(corrections[9], -8.0 / 17.0);
        
    } else if (strcmp(regime, "ISCO") == 0) {
        if (strcmp(P[pt5].Regime_ISCO, "GES") == 0) {
            h_0 = 0.00184101952815808 * pow(viscosity, 1.0 / 8.0) * pow(mass_dmnsless, -3.0 / 8.0)
                  * pow(accretion_rate_dmnsless, 1.0 / 4.0) * pow(x_isco, 1.0 / 8.0) * pow(corrections[1], -1.0 / 8.0)
                  * pow(corrections[9], -1.0 / 2.0);
            
        } else if (strcmp(P[pt5].Regime_ISCO, "RES") == 0) {
            h_0 = 0.002;
            
        } else if (strcmp(P[pt5].Regime_ISCO, "GFF") == 0) {
            h_0 = 0.0013166295323149043 * pow(viscosity, 1.0 / 17.0) * pow(mass_dmnsless, -5.0 / 17.0)
                  * pow(accretion_rate_dmnsless, 3.0 / 17.0) * pow(x_isco, 5.0 / 17.0) *
                  pow(corrections[1], -1.0 / 17.0)
                  * pow(corrections[3], -1.0 / 34.0) * pow(corrections[9], -8.0 / 17.0);
            
        } else {terminating("Invalid 'P[pt5].Regime_ISCO'"); } // Sanity check, more info in the print statement.
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
    corrections[10] = pow(2.0, -1.0 / 2.0) * viscosity * x_isco * h_0 * pow(corrections[3], 1.0 / 2.0)
                      * sqrt(corrections[9]); // P_isco
    
    double x_1 = 2.0 * cos(acos(alignment * spin) / 3.0 - RA_PI / 3.0);
    double x_2 = 2.0 * cos(acos(alignment * spin) / 3.0 + RA_PI / 3.0);
    double x_3 = -2.0 * cos(acos(alignment * spin) / 3.0);
    
    corrections[11] = corrections[10] + x - x_isco - alignment * (3.0 * spin / 2.0) * log(x / x_isco)
                      - (3.0 * pow(x_1 - alignment * spin, 2.0)) / (x_1 * (x_1 - x_2) * (x_1 - x_3))
                        * log((x - x_1) / (x_isco - x_1))
                      - (3.0 * pow(x_2 - alignment * spin, 2.0)) / (x_2 * (x_2 - x_3) * (x_2 - x_1))
                        * log((x - x_2) / (x_isco - x_2))
                      - (3.0 * pow(x_3 - alignment * spin, 2.0)) / (x_3 * (x_3 - x_1) * (x_3 - x_2))
                        * log((x - x_3) / (x_isco - x_3)); // P
    
} // correction_functions(...)


/**
 * @brief Calculate the ratio between the gas and radiation pressure or vice versa for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Dimensionless pressure ratio.
 **/
double pressure_ratio(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                      const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Calculate the pressure ratio for each 'regime'. //
    if (strcmp(regime, "GES") == 0) {
        return 68.98449619493478 * pow(viscosity, 1.0 / 10.0) * pow(mass_dmnsless, -7.0 / 10.0)
               * pow(accretion_rate_dmnsless, 4.0 / 5.0) * pow(x, -29.0 / 10.0)
               * pow(corrections[0], -9.0 / 10.0) * pow(corrections[2], 1.0 / 10.0)
               * pow(corrections[11], 4.0 / 5.0) * pow(corrections[8], -1.0 / 2.0);
    } else if (strcmp(regime, "RES") == 0) {
        return 0.000025300042906895175 * pow(viscosity, -1.0 / 4.0) * pow(mass_dmnsless, 7.0 / 4.0)
               * pow(accretion_rate_dmnsless, -2.0) * pow(x, 29.0 / 4.0) * pow(corrections[0], 9.0 / 4.0)
               * pow(corrections[2], -1.0 / 4.0) * pow(corrections[11], -2.0)
               * pow(corrections[8], 5.0 / 4.0);
    } else if (strcmp(regime, "GFF") == 0) {
        return 0.26699048406463 * pow(viscosity, -1.0 / 10.0) * pow(mass_dmnsless, -1.0 / 4.0)
               * pow(accretion_rate_dmnsless, 7.0 / 20.0) * pow(x, -11.0 / 10.0)
               * pow(corrections[0], -9.0 / 20.0)
               * pow(corrections[2], 1.0 / 10.0) * pow(corrections[11], 7.0 / 20.0)
               * pow(corrections[8], -11.0 / 40.0);
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // pressure_ratio(...)


/**
 * @brief Calculate the ratio between the free-free and electron scattering opacity or vice versa for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Dimensionless opacity ratio.
 **/
double opacity_ratio(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                     const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Calculate the opacity ratio for each 'regime'. //
    if (strcmp(regime, "GES") == 0) {
        return 4.359450454931594e-6 * mass_dmnsless * pow(accretion_rate_dmnsless, -1.0) * pow(x, 4.0)
               * corrections[0] * pow(corrections[11], -1.0) * sqrt(corrections[8]);
    } else if (strcmp(regime, "RES") == 0) {
        return 2.1927664368169804e-8 * pow(mass_dmnsless, 15.0 / 8.0) * pow(viscosity, -1.0 / 8.0)
               * pow(accretion_rate_dmnsless, -2.0) * pow(x, 61.0 / 8.0) * pow(corrections[0], 17.0 / 8.0)
               * pow(corrections[2], -1.0 / 8.0) * pow(corrections[11], -2.0) * pow(corrections[8], 9.0 / 8.0);
    } else if (strcmp(regime, "GFF") == 0) {
        return 478.9433271561746 * pow(mass_dmnsless, -1.0 / 2.0) * pow(accretion_rate_dmnsless, 1.0 / 2.0)
               * pow(x, -2.0) * pow(corrections[0], -1.0 / 2.0) * pow(corrections[11], 1.0 / 2.0)
               * pow(corrections[8], -1.0 / 4.0);
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // opacity_ratio(...)


/**
 * @brief Calculate the surface density profiles for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Surface density in gram centimeter^-2.
 **/
double surface_densities(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                         const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Surface density profiles for each regime. //
    if (strcmp(regime, "GES") == 0) {
        return 49610.3 * pow(viscosity, -4.0 / 5.0) * pow(mass_dmnsless, -2.0 / 5.0)
               * pow(accretion_rate_dmnsless, 3.0 / 5.0) * pow(x, -9.0 / 5.0) * pow(corrections[0], 1.0 / 5.0)
               * pow(corrections[2], -4.0 / 5.0) * pow(corrections[11], 3.0 / 5.0);
        
    } else if (strcmp(regime, "RES") == 0) {
        return 10.4248 * pow(viscosity, -1.0) * mass_dmnsless * pow(accretion_rate_dmnsless, -1.0) * pow(x, 4.0)
               * pow(corrections[0], 2.0) * pow(corrections[2], -1.0) * pow(corrections[11], -1.0) * corrections[8];
        
    } else if (strcmp(regime, "GFF") == 0) {
        return 170462 * pow(viscosity, -4.0 / 5.0) * pow(mass_dmnsless, -1.0 / 2.0)
               * pow(accretion_rate_dmnsless, 7.0 / 10.0) * pow(x, -11.0 / 5.0) * pow(corrections[0], 1.0 / 10.0)
               * pow(corrections[2], -4.0 / 5.0) * pow(corrections[11], 7.0 / 10.0) * pow(corrections[8], -1.0 / 20.0);
        
    } else if (strcmp(regime, "ISCO") == 0) {
        // Get 'correction_functions' at the ISCO and store them in the 'local_corrections' array. //
        double local_corrections[12];
        (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, P[pt5].Regime_ISCO,
                                    local_corrections);
        
        double surface_density_isco, radial_velocity_isco; // Declare local variable.
        double x_isco = sqrt(r_isco(bh_mass, spin, rotation) / r_grav(bh_mass)); // Dimensionless radius of the ISCO.
        
        // Surface density profiles at the ISCO for each regime. //
        if (strcmp(P[pt5].Regime_ISCO, "GES") == 0) {
            surface_density_isco = 49610.3 * pow(viscosity, -4.0 / 5.0) * pow(mass_dmnsless, -2.0 / 5.0)
                                   * pow(accretion_rate_dmnsless, 3.0 / 5.0) * pow(x_isco, -9.0 / 5.0)
                                   * pow(corrections[1], 1.0 / 5.0) * pow(corrections[3], -4.0 / 5.0)
                                   * pow(corrections[10], 3.0 / 5.0);
            
        } else if (strcmp(P[pt5].Regime_ISCO, "RES") == 0) {
            surface_density_isco = 10.4248 * pow(viscosity, -1.0) * mass_dmnsless * pow(accretion_rate_dmnsless, -1.0)
                                   * pow(x_isco, 4.0) * pow(corrections[1], 2.0) * pow(corrections[3], -1.0)
                                   * pow(corrections[10], -1.0) * corrections[9];
            
        } else if (strcmp(P[pt5].Regime_ISCO, "GFF") == 0) {
            surface_density_isco = 170462 * pow(viscosity, -4.0 / 5.0) * pow(mass_dmnsless, -1.0 / 2.0)
                                   * pow(accretion_rate_dmnsless, 7.0 / 10.0) * pow(x_isco, -11.0 / 5.0)
                                   * pow(corrections[1], 1.0 / 10.0) * pow(corrections[3], -4.0 / 5.0)
                                   * pow(corrections[10], 7.0 / 10.0) * pow(corrections[9], -1.0 / 20.0);
            
        } else {terminating("Invalid 'P[pt5].Regime_ISCO'"); } // Sanity check, more info in the print statement.
        
        if (strcmp(P[pt5].Regime_ISCO, "GES") == 0) {
            radial_velocity_isco = fabs(-725088.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -3.0 / 5.0)
                                   * pow(accretion_rate_dmnsless, 2.0 / 5.0) * pow(x_isco, -1.0 / 5.0)
                                   * pow(corrections[1], -1.0 / 5.0) * pow(corrections[3], 4.0 / 5.0)
                                   * pow(corrections[10], -3.0 / 5.0);
            
        } else if (strcmp(P[pt5].Regime_ISCO, "RES") == 0) {
            radial_velocity_isco = fabs(-3.45059e9) * viscosity * pow(mass_dmnsless, -2.0)
                                   * pow(accretion_rate_dmnsless, 2.0) * pow(x_isco, -6.0) * pow(corrections[0], -2.0)
                                   * corrections[3] * pow(corrections[9], -1.0) * corrections[10];
            
        } else if (strcmp(P[pt5].Regime_ISCO, "GFF") == 0) {
            radial_velocity_isco = fabs(-211025.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -1.0 / 2.0)
                                   * pow(accretion_rate_dmnsless, 3.0 / 10.0) * pow(x_isco, 1.0 / 5.0)
                                   * pow(corrections[1], -1.0 / 10.0) * pow(corrections[3], 4.0 / 5.0)
                                   * pow(corrections[10], -7.0 / 10.0) * pow(corrections[9], 1.0 / 20.0);
            
        } else {terminating("Invalid 'P[pt5].Regime_ISCO'"); } // Sanity check, more info in the print statement.
        
        // Calculate the surface density between the photon radius and the ISCO. //
        double epsilon = radial_velocity_isco / RA_C_CGS
                         * sqrt(3.0 * r_isco(bh_mass, spin, rotation) / (2.0 * r_grav(bh_mass)));
        
        return surface_density_isco * r_isco(bh_mass, spin, rotation) / r
               * pow(pow(epsilon, -1.0) * pow(r_isco(bh_mass, spin, rotation) / r - 1.0, 3.0 / 2.0) + 1.0, -1.0);
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // surface_densities(...)


/**
 * @brief Calculate the pressure profiles for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Pressure in dyn centimeter^-2.
 **/
double pressures(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                 const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Pressure profiles for each regime. //
    if (strcmp(regime, "GES") == 0) {
        return 1.75407e17 * pow(viscosity, -9.0 / 10.0) * pow(mass_dmnsless, -17.0 / 10.0)
               * pow(accretion_rate_dmnsless, 4.0 / 5.0) * pow(x, -59.0 / 10.0) * pow(corrections[0], 1.0 / 10.0)
               * pow(corrections[2], -9.0 / 10.0) * pow(corrections[11], 4.0 / 5.0) * pow(corrections[8], 1.0 / 2.0);
        
    } else if (strcmp(regime, "RES") == 0) {
        return 2.5427e15 * pow(viscosity, -1.0) * pow(mass_dmnsless, -1.0) * pow(x, -3.0)
               * corrections[0] * pow(corrections[2], -1.0) * corrections[8];
        
    } else if (strcmp(regime, "GFF") == 0) {
        return 3.251432754512535e17 * pow(viscosity, -9.0 / 10.0) * pow(mass_dmnsless, -7.0 / 4.0)
               * pow(accretion_rate_dmnsless, 17.0 / 20.0) * pow(x, -61.0 / 10.0) * pow(corrections[0], 1.0 / 20.0)
               * pow(corrections[2], -9.0 / 10.0) * pow(corrections[11], 17.0 / 20.0)
               * pow(corrections[8], 19.0 / 40.0);
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // pressures(...)


/**
 * @brief Calculate the density profiles for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Density in gram centimeter^-3.
 **/
double densities(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                 const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Density profiles for each regime. //
    if (strcmp(regime, "GES") == 0) {
        return 8.0637 * pow(viscosity, -7.0 / 10.0) * pow(mass_dmnsless, -11.0 / 10.0)
               * pow(accretion_rate_dmnsless, 2.0 / 5.0) * pow(x, -37.0 / 10.0) * pow(corrections[0], 3.0 / 10.0)
               * pow(corrections[2], -7.0 / 10.0) * pow(corrections[11], 2.0 / 5.0) * pow(corrections[8], 1.0 / 2.0);
        
    } else if (strcmp(regime, "RES") == 0) {
        return 0.0000245629 * pow(viscosity, -1.0) * mass_dmnsless * pow(accretion_rate_dmnsless, -2.0) * pow(x, 5.0)
               * pow(corrections[0], 3.0) * pow(corrections[2], -1.0) * pow(corrections[11], -2.0)
               * pow(corrections[8], 2.0);
        
    } else if (strcmp(regime, "GFF") == 0) {
        return 51.3593663069338 * pow(viscosity, -7.0 / 10.0) * pow(mass_dmnsless, -5.0 / 4.0)
               * pow(accretion_rate_dmnsless, 11.0 / 20.0) * pow(x, -43.0 / 10.0) * pow(corrections[0], 3.0 / 20.0)
               * pow(corrections[2], -7.0 / 10.0) * pow(corrections[11], 11.0 / 20.0)
               * pow(corrections[8], 17.0 / 40.0);
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // densities(...)


/**
 * @brief Calculate the radial velocity profiles for different regimes.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Velocity in centimeter second^-1.
 **/
double radial_velocities(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                         const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double x = sqrt(r / r_grav(bh_mass)); // Dimensionless radial coordinate.
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
    
    // Get 'correction_functions' and store them in the 'corrections' array. //
    double corrections[12];
    (void) correction_functions(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime, corrections);
    
    // Radial velocity profiles for each regime. //
    if (strcmp(regime, "GES") == 0) {
        return fabs(-725088.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -3.0 / 5.0)
               * pow(accretion_rate_dmnsless, 2.0 / 5.0) * pow(x, -1.0 / 5.0) * pow(corrections[0], -1.0 / 5.0)
               * pow(corrections[2], 4.0 / 5.0) * pow(corrections[11], -3.0 / 5.0);
        
    } else if (strcmp(regime, "RES") == 0) {
        return fabs(-3.45059e9) * viscosity * pow(mass_dmnsless, -2.0) * pow(accretion_rate_dmnsless, 2.0)
               * pow(x, -6.0) * pow(corrections[0], -2.0) * corrections[2] * pow(corrections[8], -1.0)
               * corrections[11];
        
    } else if (strcmp(regime, "GFF") == 0) {
        return fabs(-211025.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -1.0 / 2.0)
               * pow(accretion_rate_dmnsless, 3.0 / 10.0) * pow(x, 1.0 / 5.0) * pow(corrections[0], -1.0 / 10.0)
               * pow(corrections[2], 4.0 / 5.0) * pow(corrections[11], -7.0 / 10.0) * pow(corrections[8], 1.0 / 20.0);
        
    } else if (strcmp(regime, "ISCO") == 0) {
        
        double radial_velocity_isco; // Declare local variable.
        double x_isco = sqrt(r_isco(bh_mass, spin, rotation) / r_grav(bh_mass)); // Dimensionless radius of the ISCO.
        
        // Radial velocity profiles at the ISCO for each regime. //
        if (strcmp(P[pt5].Regime_ISCO, "GES") == 0) {
            radial_velocity_isco = fabs(-725088.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -3.0 / 5.0)
                                   * pow(accretion_rate_dmnsless, 2.0 / 5.0) * pow(x_isco, -1.0 / 5.0)
                                   * pow(corrections[1], -1.0 / 5.0) * pow(corrections[3], 4.0 / 5.0)
                                   * pow(corrections[10], -3.0 / 5.0);
            
        } else if (strcmp(P[pt5].Regime_ISCO, "RES") == 0) {
            radial_velocity_isco = fabs(-3.45059e9) * viscosity * pow(mass_dmnsless, -2.0)
                                   * pow(accretion_rate_dmnsless, 2.0) * pow(x_isco, -6.0) * pow(corrections[0], -2.0)
                                   * corrections[3] * pow(corrections[9], -1.0) * corrections[10];
            
        } else if (strcmp(P[pt5].Regime_ISCO, "GFF") == 0) {
            radial_velocity_isco = fabs(-211025.0) * pow(viscosity, 4.0 / 5.0) * pow(mass_dmnsless, -1.0 / 2.0)
                                   * pow(accretion_rate_dmnsless, 3.0 / 10.0) * pow(x_isco, 1.0 / 5.0)
                                   * pow(corrections[1], -1.0 / 10.0) * pow(corrections[3], 4.0 / 5.0)
                                   * pow(corrections[10], -7.0 / 10.0) * pow(corrections[9], 1.0 / 20.0);
            
        } else {terminating("Invalid 'P[pt5].Regime_ISCO'"); } // Sanity check, more info in the print statement.
        
        // Calculate the surface density between the photon radius and the ISCO. //
        double epsilon = radial_velocity_isco / RA_C_CGS
                         * sqrt(3.0 * r_isco(bh_mass, spin, rotation) / (2.0 * r_grav(bh_mass)));
        
        return radial_velocity_isco
               * (pow(epsilon, -1.0) * pow(r_isco(bh_mass, spin, rotation) / r - 1.0, 3.0 / 2.0) + 1.0);
        
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
} // radial_velocities(...)


/**
 * @brief Calculate the local relativistic Toomre Q parameter.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Dimensionless Toomre parameter.
 **/
double toomre_parameter(int pt5, double r, double viscosity, double bh_mass, double accretion_rate, double spin,
                        const char *rotation, const char *regime) {
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    double alignment;
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Calculate the relativistic angular frequency. //
    double omega = pow(RA_C_CGS, 3.0)
                   / ((RA_G_CGS * bh_mass) * (pow(sqrt(r / r_grav(bh_mass)), 3.0) + alignment * spin));
    
    // Calculate the relativistic epicyclic frequency. //
    double kappa = omega * sqrt(1.0 - 6.0 * pow(r, -2.0) + alignment * 8.0 * spin * pow(r, -3.0)
                                - 3.0 * pow(spin, 2.0) * pow(r, -4.0));
    
    // Calculate the relativistic sound speed . //
    double sound_speed = sqrt(pressures(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime)
                              / densities(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime));
    
    return kappa * sound_speed
           / (RA_PI * RA_G_CGS * surface_densities(pt5, r, viscosity, bh_mass, accretion_rate, spin, rotation, regime));
    
} // toomre_parameter(...)


/**
 * @brief Calculate the viscosity profile.
 * @param r Radial coordinate.
 * @param bh_mass Black hole mass in gram.
 * @param spin Dimensionless spin parameter.
 * @return Dimensionless viscosity parameter.
 **/
double viscosities(double r, double bh_mass, double spin) {
    
    // Calculate the viscosity. //
    double q = 3.0 / 2.0
               * (1 - 2 * pow(r / r_grav(bh_mass), -1.0) + pow(spin, 2.0) * pow(r / r_grav(bh_mass), -2.0))
               / (1 - 3 * pow(r / r_grav(bh_mass), -1.0) + 2 * spin * pow(r / r_grav(bh_mass), -3.0 / 2.0));
    
    return 0.025 * pow(q / (3.0 / 2.0), 6.0);
    
} // viscosities(...)


/**
 * @brief Calculate the precession timescale due to the Lense-Thirring effect.
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @return Timescale in seconds.
 **/
double precession_timescale(int pt5, double r) {
    
    // Calculate the black hole's spin vector magnitude and convert to cgs units. //
    double spin_magnitude = sqrt(pow(P[pt5].Spin[0], 2.0) + pow(P[pt5].Spin[1], 2.0) + pow(P[pt5].Spin[2], 2.0));
    double spin_magnitude_cgs = convert_code_to_cgs_units("angmomentum", spin_magnitude);
    
    // Calculate the angular velocity and precession timescale. //
    double omega = (2.0 * RA_G_CGS * spin_magnitude_cgs) / (pow(RA_C_CGS, 2.0) * pow(r, 3.0));
    double t_precession = (2.0 * RA_PI) / omega;
    
    // Check if the precession timescale is positive. If yes, return it, else terminate. //
    if (t_precession > 0.0) { return t_precession; }
    else {terminating("Invalid 't_precession'")}
    
} // precession_timescale(...)


/**
 * @brief Calculate the Bardeen-Petterson timescale .
 * @param pt5 Index of an active black hole particle.
 * @param r Radial coordinate.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Timescale in seconds.
 **/
double bardeen_petterson_timescale(int pt5, double r, double viscosity, double bh_mass, double accretion_rate,
                                   double spin, const char *rotation) {
    
    double cumulative_timescale = 0.0;
    // Loop over 'NUM_VALID' radial bins. //
    for (int bin = 0; bin < NUM_VALID; ++bin) {
        // Find the regime 'r' belongs to and integrate the . //
        if (r >= P[pt5].GES_starts[bin] && r < P[pt5].GES_stops[bin] && P[pt5].GES_starts[bin] > 0.0
            && P[pt5].GES_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].GES_starts[bin], r, rotation, "GES");
        }
        if (r >= P[pt5].RES_starts[bin] && r < P[pt5].RES_stops[bin] && P[pt5].RES_starts[bin] > 0.0
            && P[pt5].RES_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].RES_starts[bin], r, rotation, "RES");
        }
        if (r >= P[pt5].GFF_starts[bin] && r < P[pt5].GFF_stops[bin] && P[pt5].GFF_starts[bin] > 0.0
            && P[pt5].GFF_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].GFF_starts[bin], r, rotation, "GFF");
        }
        if (r > P[pt5].GES_starts[bin] && r >= P[pt5].GES_stops[bin] && P[pt5].GES_starts[bin] > 0.0
            && P[pt5].GES_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].GES_starts[bin], P[pt5].GES_stops[bin], rotation, "GES");
        }
        if (r > P[pt5].RES_starts[bin] && r >= P[pt5].RES_stops[bin] && P[pt5].RES_starts[bin] > 0.0
            && P[pt5].RES_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].RES_starts[bin], P[pt5].RES_stops[bin], rotation, "RES");
        }
        if (r > P[pt5].GFF_starts[bin] && r >= P[pt5].GFF_stops[bin] && P[pt5].GFF_starts[bin] > 0.0
            && P[pt5].GFF_stops[bin] > 0.0) {
            cumulative_timescale += integrate(pt5, "t_viscous", viscosity, bh_mass, accretion_rate, spin,
                                              P[pt5].GFF_starts[bin], P[pt5].GFF_stops[bin], rotation, "GFF");
        }
    }
    
    double t_bardeen_petterson = 4.0 * (1.0 + 7.0 * pow(viscosity, 2.0)) * cumulative_timescale
                                 / (2.0 * pow(viscosity, 2.0) * (4.0 + pow(viscosity, 2.0)));
    
    // Check if the Bardeen-Petterson timescale is positive. If yes, return it, else terminate. //
    if (t_bardeen_petterson > 0.0) { return t_bardeen_petterson; }
    else {terminating("Invalid 't_bardeen_petterson'")}
} // bardeen_petterson_timescale(...)