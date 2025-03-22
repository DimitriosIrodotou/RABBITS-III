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

// Include user-defined header files. //
#include "RA_vars.h"


// Declare global variables. //
const char *MODE = "coherent"; // Select between 'chaotic' or 'coherent' for 'accrete_angular_momenta()'.

const double VISCOSITY = 0.1; // Dimensionless viscosity parameter.
const double RNG_SEED = 1312.0; // Seed value for the random number generator (RNG).
const double UNIT_TEST_TOLERANCE = 100.0; // In centimeter for distance comparisons, otherwise dimensionless.

FILE *FdHoyleLyttletonAll; // File handle for 'relativistic_accretion_logs' per CPU.
FILE *FdAccretionDiscsAll; // File handle for 'relativistic_accretion_logs' for all CPUs.

const int NUM_STEPS = 10; // Factor to adjust the accretion disc's extent when calculating Toomre Q parameter.
const int NUM_BINS = 1000; // Number of radial bins the intra-ISCO and extra-ISCO regions are divided into.
const int MAX_ITERATIONS = 10; // Number of maximum iterations allowed.
const int MAX_INTERVALS = 1000; // Maximum number of sub-intervals the integration region is divided into.

struct integrand_parameters *integrand_parameters;