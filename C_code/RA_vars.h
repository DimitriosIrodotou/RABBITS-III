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

#ifndef RA_VARS_H
#define RA_VARS_H

// Include existing library files. //
#include <stdio.h>
#include <stdbool.h>

// Declare constants in cgs units. //
#define RA_PI               3.14159265358979323846
#define RA_C_CGS            3.0e10 // In centimeter second^-1.
#define RA_G_CGS            6.67e-8 // In centimeter gram^-1 centimeter^2 second^-2.
#define RA_THOMPSON_CGS     6.65e-25 // In centimeter^2.
#define RA_SOLAR_MASS_CGS   1.99e33 // In gram.
#define RA_PROTON_MASS_CGS  1.67e-24 // In gram.

// Convert constants from cgs to code units. //
#define RA_C_CODE           convert_cgs_to_code_units("light", RA_C_CGS)
#define RA_G_CODE           convert_cgs_to_code_units("gravity", RA_G_CGS)

// Declare macros. //
#define  terminating(x) { char termbuf[2000]; \
sprintf(termbuf, "'Terminating! %s in function:'%s()', file:'%s, task:%d' \n", x, __FUNCTION__, __FILE__, ThisTask); \
printf("%s", termbuf); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1); exit(EXIT_FAILURE); }

#define  terminating_scope(text, function, file) { char termbuf[2000]; \
sprintf(termbuf, "'Terminating! %s in function:'%s()', file:'%s, task:%d' \n", text, function, file, ThisTask); \
printf("%s", termbuf); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1); exit(EXIT_FAILURE); }

#define NUM_VALID 10 // Number of valid starting/stopping radii for regimes.

// Declare global variables. //
extern const char *MODE; // Select in RA_vars.c between 'chaotic' or 'coherent' for 'accrete_angular_momenta()'.

extern const double VISCOSITY; // Dimensionless viscosity parameter.
extern const double RNG_SEED; // Seed value for the random number generator (RNG).
extern const double UNIT_TEST_TOLERANCE; // In centimeter for distance comparisons, otherwise dimensionless.

extern FILE *FdHoyleLyttletonAll; // File handle for 'relativistic_accretion_logs' per CPU.
extern FILE *FdAccretionDiscsAll; // File handle for 'relativistic_accretion_logs' for all CPUs.

extern const int NUM_STEPS; // Factor to adjust the accretion disc's extent when calculating Toomre Q parameter.
extern const int NUM_BINS; // Number of radial bins the intra-ISCO and extra-ISCO regions are divided into.
extern const int MAX_ITERATIONS; // Number of maximum iterations allowed.
extern const int MAX_INTERVALS; // Maximum number of sub-intervals the integration region is divided into.

// A structure for the parameters of the integrand function. //
struct integrand_parameters {
    
    int pt5;
    double viscosity;
    double bh_mass;
    double accretion_rate;
    double spin;
    const char *rotation;
    const char *regime;
    
};

#endif // RA_VARS_H