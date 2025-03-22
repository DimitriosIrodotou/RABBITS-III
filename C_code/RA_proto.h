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

#ifndef RA_PROTO_H
#define RA_PROTO_H

// Declare prototype functions for RA_main.c. //
int relativistic_accretion(void);

// Declare prototype functions for RA_accretion_disc.c. //
double integrated_profiles(int, const char *, double, double, double, double, const char *);

double accretion_rate(int, double, double, double, double, double, const char *);

double angmomentum_integrand(double, void *);

double t_viscous_integrand(double, void *);

double spin_integrand(double, void *);

// Declare prototype functions for RA_profiles.c. //
void *correction_functions(int, double, double, double, double, double, const char *, const char *, double *);

double pressure_ratio(int, double, double, double, double, double, const char *, const char *);

double opacity_ratio(int, double, double, double, double, double, const char *, const char *);

double surface_densities(int, double, double, double, double, double, const char *, const char *);

double pressures(int, double, double, double, double, double, const char *, const char *);

double densities(int, double, double, double, double, double, const char *, const char *);

double radial_velocities(int, double, double, double, double, double, const char *, const char *);

double toomre_parameter(int, double, double, double, double, double, const char *, const char *);

double viscosities(double, double, double);

double precession_timescale(int, double);

double bardeen_petterson_timescale(int, double, double, double, double, double, const char *);

// Declare prototype functions for RA_radii.c. //
void accretion_disc_extent(int, double, double, double, double, const char *);

void ranges_of_validity(int, double, double, double, double, const char *);

double r_grav(double);

double r_outermost(int, double);

double r_isco(double, double, const char *);

double r_photon(double, double, const char *);

void r_warp(int, double, double, double, double, const char *);

// Declare prototype functions for RA_feeding.c. //
void capture_SPH_particles(void);

int find_neighbours(double [3], double, int, int *, int, int *, int *, int);

double bondi_hoyle_accretion_rate(int);

double bondi_hoyle_timescale(int);

// Declare prototype functions for RA_feedback.c. //
double radiative_efficiency_spin(double, double, const char *);

double radiative_efficiency_isco(double, double, double, double, const char *, const char *);

void thermal_feedback(void);

// Declare prototype functions for RA_spins.c. //
void initialise_black_hole_spin(int);

void exchange_angular_momenta(int, double, double, double, const char *);

void accrete_angular_momenta(int, const char *, double);

const char *assess_alignment(int);

double calculate_spin_parameter(int);

void Bardeen_Petterson_effect(int, double, double, double, double, double, const char *);

double isco_angmomentum(int, double, const char *);

// Declare prototype functions for RA_utilities.c. //
void initialise_setup(void);

void initialise_log_files(int);

void log_details(int, const char *);

void reset_accretion_disc(int);

double convert_code_to_cgs_units(const char *, double);

double convert_cgs_to_code_units(const char *, double);

void *logspace(double, double, int, double *);

double maximum_array_value(int, const double *);

double interpolate(double, double, double, double);

void *resolve_overlapping_regimes(int, double, double, double, double, const char *, double *, double *, double *);

void *resolve_lone_regimes(int, double, double, double, double, const char *, double *, double *, double *);

void *resolve_empty_regimes(int, double, double, double, double, const char *, double *, double *, double *);

double integrate(int, const char *, double, double, double, double, double, double, const char *, const char *);

double timestep_limiter(int, double);

void *sort_ranges_of_validity(int, double *, double *);

double hubble_parameter();

// Declare prototype functions for RA_tests.c. //

void *local_tests(int, double, double, const char *, double *, double *, double *, const char *, const char *);

void global_tests(int, double, const char *, const char *);

void radii_unit_tests(int);

void accretion_disc_unit_tests(int);

#endif // RA_PROTO_H