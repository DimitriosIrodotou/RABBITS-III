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
#include <stdlib.h>

// Include user-defined header files. //
#include "RA_vars.h"
#include "RA_proto.h"
#include "../../../proto.h"
#include "../../../kernel.h"
#include "../../../allvars.h"


// Declare structures to hold data to be sent ('data_in') and received ('data_out') between CPUs. //
struct data_in {
    
    double Dt;
    double Hsml;
    double Pos[3];
    double BH_Mdot;
    double Drain_Flag;
    double Epsilon_r;
    double BH_Density;
    int NodeList[NODELISTLENGTH];
#ifdef TIMESTEP_LIMITER_FB
    int TimeStep;
#endif // TIMESTEP_LIMITER_FB

};

struct data_out {
    double Energy;
};


// Declare the internal function. //
static int perform_thermal_feedback(int target, int mode, int *nexport, int *nSend_local, struct data_in *DataGet,
                                    struct data_out *DataResult, int use_part_type);


/**
 * @brief Calculate the spin-dependent radiative efficiency.
 * @param bh_mass Black hole mass in gram.
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @return Spin-dependent radiative efficiency.
 **/
double radiative_efficiency_spin(double bh_mass, double spin, const char *rotation) {
    
    // Sanity check, more info in the print statement. //
    if (strcmp(rotation, "co-rotation") != 0 && strcmp(rotation, "counter-rotation") != 0) {
        terminating("Invalid 'rotation'");
    }
    
    return 1.0 - sqrt(1.0 - 2.0 * r_grav(bh_mass) / (3.0 * r_isco(bh_mass, spin, rotation)));
    
} // radiative_efficiency_spin(...)


/**
 * @brief Calculate the relativistic radiative efficiency at the ISCO.
 * @param viscosity Dimensionless viscosity parameter.
 * @param bh_mass Black hole mass in gram.
 * @param accretion_rate Black hole accretion rate in gram second^-1 .
 * @param spin Dimensionless spin parameter.
 * @param rotation Select between 'co-rotation' or 'counter-rotation'.
 * @param regime Select between 'Gas-ES', 'Rad-ES', or 'Gas-FF'.
 * @return Spin-dependent radiative efficiency.
 **/
double radiative_efficiency_isco(double viscosity, double bh_mass, double accretion_rate, double spin,
                                 const char *rotation, const char *regime) {
    
    // Declare local variables. //
    double h_0, alignment;
    
    // 'co-rotation'/'counter-rotation' corresponds to positive/negative 'alignment'. //
    if (strcmp(rotation, "co-rotation") == 0) { alignment = 1.0; }
    else if (strcmp(rotation, "counter-rotation") == 0) { alignment = -1.0; }
    else {terminating("Invalid 'rotation'"); } // Sanity check, more info in the print statement.
    
    // Calculate the 'corrections' functions. //
    double corrections[12], x_isco = sqrt(r_isco(bh_mass, spin, rotation) / r_grav(bh_mass));
    corrections[1] = 1.0 - 3.0 * pow(x_isco, -2.0) + alignment * 2.0 * spin * pow(x_isco, -3.0); // C_isco
    corrections[3] = 1.0 - 2.0 * pow(x_isco, -2.0) + pow(spin, 2.0) * pow(x_isco, -4.0); // D_isco
    corrections[5] = alignment * (1.0 - alignment * 2.0 * spin * pow(x_isco, -3.0)
                                  + pow(spin, 2.0) * pow(x_isco, -4.0)); // F_isco
    corrections[7] = 1.0 - 2.0 * pow(x_isco, -2.0) + alignment * spin * pow(x_isco, -3.0); // G_isco
    corrections[9] = pow(corrections[1], -1.0) * pow(corrections[5], 2.0)
                     - pow(spin, 2.0) * pow(x_isco, -2.0)
                       * (pow(corrections[1], -1.0 / 2.0) * corrections[7] - 1.0); // R_isco
    
    // Calculate the dimensionless mass and accretion rate. //
    double mass_dmnsless = bh_mass / (3.0 * RA_SOLAR_MASS_CGS), accretion_rate_dmnsless = accretion_rate / 1.0e17;
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
    } else {terminating("Invalid 'regime'"); } // Sanity check, more info in the print statement.
    
    corrections[10] = pow(2.0, -1.0 / 2.0) * viscosity * x_isco * h_0
                      * sqrt(corrections[3]) * sqrt(corrections[9]); // P_isco
    
    // Relativistic angular velocity at the ISCO. //
    double omega_isco = alignment * sqrt(RA_G_CGS * bh_mass) / (pow(r_isco(bh_mass, spin, rotation), 3.0 / 2.0)
                                                                + alignment * spin * pow(r_grav(bh_mass), 3.0 / 2.0));
    
    // Calculate the spin-depended corrections from Bardeen+72 (1972ApJ...178..347B). //
    double z_1 = 1 + pow(1.0 - pow(spin, 2.0), 1.0 / 3.0) * (pow(1.0 + spin, 1.0 / 3.0) + pow(1.0 - spin, 1.0 / 3.0));
    double z_2 = sqrt(3.0 * pow(spin, 2.0) + pow(z_1, 2.0));
    double z = 3.0 + z_2 - alignment * sqrt((3.0 - z_1) * (3.0 + z_1 + 2.0 * z_2));
    
    // Calculate the relativistic specific angular momentum at the ISCO. //
    double J_isco = alignment * RA_G_CGS * bh_mass / RA_C_CGS
                    * (pow(z, 2.0) - alignment * 2.0 * spin * sqrt(z) + pow(spin, 2.0))
                    / (z * sqrt(z - 3.0 + alignment * 2.0 * spin * pow(z, -1.0 / 2.0)));
    
    // Additional torque at the ISCO. //
    double E_isco = sqrt(1.0 - 2.0 * r_grav(bh_mass) / (3.0 * r_isco(bh_mass, spin, rotation)));
    double g_isco = bh_mass * accretion_rate * corrections[10]
                    * pow(pow(RA_C_CGS, 2.0) * E_isco - omega_isco * J_isco, -1.0);
    
    return RA_G_CGS * g_isco * omega_isco / (RA_C_CGS * accretion_rate);
    
} // radiative_efficiency_isco(...)


/**
 * @brief Calculate the released energy in the form of thermal feedback and distribute it by following a kernel
 * weighted method to the neighbouring gas particles that reside within the black hole's smoothing length by increasing t
 * heir internal energy.
 **/
void thermal_feedback(void) {
    
    // Dynamically allocate memory to hold indices of neighbouring particles equal to 'NumPart' in the local CPU.
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    
    // Calculate the number of particles that can fit into a memory buffer of size 'BufferSize'. //
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) /
                           (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct data_in) +
                            sizeof(struct data_out) + sizemax(sizeof(struct data_in), sizeof(struct data_out))));
    
    // Dynamically allocate memory to hold the indices of the CPUs the exported particles belong to. //
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    // Declare local variables and pointers to the structures used in the MPI communication. //
    struct data_in *DataIn, *DataGet;
    struct data_out *DataResult, *DataOut;
    int nexport, nimport, ndone_flag, ndone;
    double hubble_a = All.ComovingIntegrationOn ? hubble_function(All.Time) : 1.0;
    
    MPI_Status status; // Get access to additional information when receiving data from other CPUs.
    
    /* Start from the 'FirstActiveParticle' in this time step and loop over CPUs to add the thermal energy from each
     * active BH particle to its gas neighbours. */
    int i = FirstActiveParticle;
    do {
        // Loop over 'NTask' number of CPUs and initialise 'Send_count' to 0 and 'Exportflag' to -1
        for (int j = 0; j < NTask; j++) {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        
        /* Loop over local (i.e. belonging to this CPU) active (i.e. the ones that need a force update at the current
         * time step) particles and 'perform_capture'. */
        for (nexport = 0; i >= 0; i = NextActiveParticle[i]) {
            if (P[i].Type == 5 && P[i].BH_Mass > 0.0) {
                int use_part_type = 1; // = 2^0 since only gas particles will be captured.
                
                /* Since mode=0, declare and use "dummy" structure variables because this is the local CPU so there will
                 * be no exchange of information between CPUs */
                int mode = 0;
                struct data_in dummy1;
                struct data_out dummy2;
                if (perform_thermal_feedback(i, mode, &nexport, Send_count, &dummy1, &dummy2, use_part_type) < 0) {
                    break; // Break the 'NextActiveParticle' loop if no neighbours were found.
                }
            }
        }
        
        // Sort the indices of CPUs. //
        MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
        
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        
        nimport = 0;
        Recv_offset[0] = 0, Send_offset[0] = 0;
        for (int j = 0; j < NTask; j++) {
            nimport += Recv_count[j];
            if (j > 0) {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        // Dynamically allocate memory to hold the structure of data to be sent ('data_in') between CPUs. //
        DataGet = (struct data_in *) mymalloc("DataGet", nimport * sizeof(struct data_in));
        DataIn = (struct data_in *) mymalloc("DataIn", nexport * sizeof(struct data_in));
        
        /* Fill the 'DataIn' structure with black hole information to be sent and used by other CPUs. The properties
         * below should much the ones in 'struct data_in' */
        for (int j = 0; j < nexport; j++) {
            int place = DataIndexTable[j].Index;
            
            DataIn[j].Hsml = P[place].Hsml;
            DataIn[j].BH_Mdot = P[place].BH_Mdot;
            DataIn[j].Epsilon_r = P[place].Epsilon_r;
            DataIn[j].BH_Density = P[place].b1.BH_Density;
            DataIn[j].Drain_Flag = P[place].Drain_Flag;
            DataIn[j].Dt = (P[place].TimeBin ? (1 << P[place].TimeBin) : 0.0) * All.Timebase_interval / hubble_a;
            for (int k = 0; k < 3; k++) { DataIn[j].Pos[k] = P[place].Pos[k]; }
            
            memcpy(DataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#ifdef TIMESTEP_LIMITER_FB
            DataIn[j].TimeStep = P[place].TimeStep;
#endif // TIMESTEP_LIMITER_FB
        }
        
        // Exchange particle data
        for (int ngrp = 1; ngrp < (1 << PTask); ngrp++) {
            int recvTask = ThisTask ^ ngrp;
            
            if (recvTask < NTask) {
                if (Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {
                    MPI_Sendrecv(&DataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct data_in), MPI_BYTE, recvTask, TAG_DENS_A,
                                 &DataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct data_in), MPI_BYTE, recvTask, TAG_DENS_A,
                                 MPI_COMM_WORLD, &status);
                }
            }
        }
        
        myfree(DataIn); // Free the dynamically allocated memory.
        
        DataResult = (struct data_out *) mymalloc("DataResult", nimport * sizeof(struct data_out));
        DataOut = (struct data_out *) mymalloc("DataOut", nexport * sizeof(struct data_out));
        
        // Now 'perform_thermal_feedback' for the particles that were sent from the other CPUs. //
        for (int j = 0; j < nimport; j++) {
            int use_part_type = 1; // = 2^0 since only gas particles will be captured.
            int mode = 1;
            int dummy;
            perform_thermal_feedback(j, mode, &dummy, &dummy, DataGet, DataResult, use_part_type);
        }
        
        if (i < 0) { ndone_flag = 1; }
        else { ndone_flag = 0; }
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        // Get the result
        for (int ngrp = 1; ngrp < (1 << PTask); ngrp++) {
            int recvTask = ThisTask ^ ngrp;
            
            if (recvTask < NTask) {
                if (Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {
                    MPI_Sendrecv(&DataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct data_out), MPI_BYTE, recvTask, TAG_DENS_B,
                                 &DataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct data_out), MPI_BYTE, recvTask, TAG_DENS_B,
                                 MPI_COMM_WORLD, &status);
                }
            }
        }
        /* Collect the results in the b4 and b6 unions, so you can transfer the 'accreted_mass', 'accreted_momentum',
        * and 'accreted_angmomentum' to the black hole particle and update its properties in 'RA_main.c'. */
        for (int j = 0; j < nexport; j++) {
            int place = DataIndexTable[j].Index;
            P[place].AD_Energy += DataOut[j].Energy;
        }
        
        // Free the dynamically allocated memory. //
        myfree(DataOut);
        myfree(DataResult);
        myfree(DataGet);
    } while (ndone < NTask);
    
    // Free the dynamically allocated memory. //
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
}


/**
 * @brief
 * @param target
 * @param mode
 * @param nexport
 * @param nSend_local
 * @param DataGet
 * @param DataResult
 * @param use_part_type
 * @param pt5 Index of an active black hole particle.
 * @return
 **/
static int perform_thermal_feedback(int target, int mode, int *nexport, int *nSend_local, struct data_in *DataGet,
                                    struct data_out *DataResult, int use_part_type) {
    
    // Declare and initialise local variables. //
    int drain_flag;
    double *pos, rho, mdot, dt, h_i, epsilon_r, energy = 0.0;
    double hubble_a = All.ComovingIntegrationOn ? hubble_function(All.Time) : 1.0;
#ifdef TIMESTEP_LIMITER_FB
    int timestep;
#endif // TIMESTEP_LIMITER_FB
    
    /* Check if this is the CPU that contains the black hole (i.e. mode == 0). If yes, extract its information directly.
     * Else, use the 'DataGet' structure to extract information sent from other CPUs */
    if (mode == 0) {
        pos = P[target].Pos;
        h_i = P[target].Hsml;
        mdot = P[target].BH_Mdot;
        rho = P[target].b1.BH_Density;
        epsilon_r = P[target].Epsilon_r;
        drain_flag = P[target].Drain_Flag;
        dt = (P[target].TimeBin ? (1 << P[target].TimeBin) : 0) * All.Timebase_interval / hubble_a;
#ifdef TIMESTEP_LIMITER_FB
        timestep = P[target].TimeStep;
#endif // TIMESTEP_LIMITER_FB
    } else {
        dt = DataGet[target].Dt;
        pos = DataGet[target].Pos;
        h_i = DataGet[target].Hsml;
        mdot = DataGet[target].BH_Mdot;
        rho = DataGet[target].BH_Density;
        drain_flag = DataGet[target].Drain_Flag;
        epsilon_r = DataGet[target].Epsilon_r;
#ifdef TIMESTEP_LIMITER_FB
        timestep = DataGet[target].TimeStep;
#endif // TIMESTEP_LIMITER_FB
    }
    
    /* The index convention for accessing tree nodes is the following: the indices 0...NumPart-1 reference
     * single particles, the indices All.MaxPart...All.MaxPart+nodes-1 reference tree nodes. `Nodes' is shifted
     * such that Nodes[All.MaxPart] gives the first tree node. */
    int start_node, list_index = 0;
    if (mode == 0) { start_node = All.MaxPart; } // Start from the root node.
    else {
        start_node = DataGet[target].NodeList[0];
        start_node = Nodes[start_node].u.d.nextnode; // Open the start node.
    }
    
    // Loop over oct-tree nodes and
    while (start_node >= 0) {
        // Loop over oct-tree nodes and
        while (start_node >= 0) {
            // Search for neighbouring SPH particles (i.e. for use_part_type = 1). If none, return -1. //
            int number_of_ngbrs =
                    find_neighbours(pos, h_i, target, &start_node, mode, nexport, nSend_local, use_part_type);
            if (number_of_ngbrs < 0) { return -1; }
            
            // Loop over neighbouring SPH particles (i.e. with Type = 0). //
            for (int n = 0; n < number_of_ngbrs; n++) {
                int pt0 = Ngblist[n];
                if (P[pt0].Type == 0 && P[pt0].Mass > 0.0) {
                    // Calculate the spatial difference between gas particles and the black hole. //
                    double dx = pos[0] - P[pt0].Pos[0], dy = pos[1] - P[pt0].Pos[1], dz = pos[2] - P[pt0].Pos[2];

#ifdef PERIODIC
                    dx = NEAREST_X(dx), dy = NEAREST_Y(dy), dz = NEAREST_Z(dz);
#endif // PERIODIC
                    double r_i = sqrt(dx * dx + dy * dy + dz * dz);
                    
                    // Check if a gas particle is within the black hole's smoothing length. If yes,
                    if (r_i < h_i) {
#ifdef TIMESTEP_LIMITER_FB
                        /* Check if an inactive neighbour (i.e. whose 'TimeBinActive' is 0) has a longer timestep
                         * than this active particle. If yes, set the timestep for the inactive particle equal to
                         * that of this active particle */
                        int bin = get_timestep_bin(P[pt0].TimeStep); // Time bin this neighbour's timestep is in.
                        if (TimeBinActive[bin] == 0) {
                            if (P[pt0].TimeStep > timestep) { P[pt0].TimeStep = timestep; }
                        }
#endif // TIMESTEP_LIMITER_FB
                        
                        // Calculate the SPH kerner weights. //
                        double hinv, hinv3, hinv4, wk, dwk;
                        
                        kernel_hinv(h_i, &hinv, &hinv3, &hinv4);
                        double u = r_i * hinv;
                        kernel_main(u, hinv3, hinv4, &wk, &dwk, 0);
                        
                        // Calculate the thermal feedback energy. //
                        energy = epsilon_r * All.KetjuData.options.bh_thermal_fb_efficiency * mdot * dt
                                 * pow(RA_C_CODE, 2.0);
                        
                        /* Check if radiative efficiency has not been calculated in 'accretion_rate' (i.e. if there has been an accretion
                            * episode happened in this time step). If yes, return and don't 'perform_thermal_feedback'. */
                        if (drain_flag == 1) { energy = 0.0; }
                        
                        if (rho > 0) { SphP[pt0].i.dInjected_BH_Energy += FLT(energy * P[pt0].Mass * wk / rho); }
                    }
                }
            }
        }
        
        if (mode == 1) {
            list_index++;
            if (list_index < NODELISTLENGTH) {
                start_node = DataGet[target].NodeList[list_index];
                if (start_node >= 0) {
                    start_node = Nodes[start_node].u.d.nextnode; // Open the start node.
                }
            }
        }
    }
    
    // Collect the results, so you can transfer the 'energy' and export it in 'log_details'.
    if (mode == 0) { P[target].AD_Energy = energy; }
    else { DataResult[target].Energy = energy; }
    
    return 0;
    
}

