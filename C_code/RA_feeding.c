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
#include "../../../proto.h"
#include "../../../kernel.h"
#include "../../../allvars.h"


// Declare structures to hold data to be sent ('data_in') and received ('data_out') between CPUs. //
struct data_in {
    
    double Hsml;
    double Mass;
    double Pos[3];
    double Vel[3];
    double BH_Mass;
    int NodeList[NODELISTLENGTH];
#ifdef TIMESTEP_LIMITER_FB
    int TimeStep;
#endif // TIMESTEP_LIMITER_FB

};

struct data_out {
    
    double Accreted_Mass;
    double Accreted_Momentum[3];
    double Accreted_AngMomentum[3];
    
};


// Declare the internal function. //
static int perform_capture(int target, int mode, int *nexport, int *nSend_local, struct data_in *DataGet,
                           struct data_out *DataResult, int use_part_type);


/**
 * @brief Ask all CPUs to find and return gas particles that belong to them and meet the capture criteria.
 **/
void capture_SPH_particles(void) {
    
    // Dynamically allocate memory to hold indices of neighbouring particles equal to 'NumPart' in the local CPU. //
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    
    // Calculate the number of particles that can fit into a memory buffer of size 'BufferSize'. //
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024)
                           / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct data_in)
                              + sizeof(struct data_out) + sizemax(sizeof(struct data_in), sizeof(struct data_out))));
    
    // Dynamically allocate memory to hold the indices of the CPUs the exported particles belong to. //
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    // Declare local variables and pointers to the structures used in the MPI communication. //
    struct data_in *DataIn, *DataGet;
    struct data_out *DataResult, *DataOut;
    int nexport, nimport, ndone_flag, ndone;
    
    MPI_Status status; // Get access to additional information when receiving data from other CPUs.
    
    /* Start from the 'FirstActiveParticle' in this time step and loop over CPUs to get the mass of gas particles
     * around each active black hole that are within the Hoyle-Lyttleton radius and have angular momenta less than the
     * circular. */
    int i = FirstActiveParticle;
    do {
        // Loop over 'NTask' number of CPUs and initialise 'Send_count' to 0 and 'Exportflag' to -1. //
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
                if (perform_capture(i, mode, &nexport, Send_count, &dummy1, &dummy2, use_part_type) < 0) {
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
         * below should much the ones in 'struct data_in'. */
        for (int j = 0; j < nexport; j++) {
            int place = DataIndexTable[j].Index;
            
            DataIn[j].Hsml = P[place].Hsml;
            DataIn[j].Mass = P[place].Mass;
            DataIn[j].BH_Mass = P[place].BH_Mass;
            for (int k = 0; k < 3; k++) { DataIn[j].Pos[k] = P[place].Pos[k]; }
            for (int k = 0; k < 3; k++) { DataIn[j].Vel[k] = P[place].Vel[k]; }
            
            memcpy(DataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#ifdef TIMESTEP_LIMITER_FB
            DataIn[j].TimeStep = P[place].TimeStep;
#endif // TIMESTEP_LIMITER_FB
        }
        
        // Exchange particle data. //
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
        
        // Now 'perform_capture' for the particles that were sent from the other CPUs. //
        for (int j = 0; j < nimport; j++) {
            int use_part_type = 1; // = 2^0 since only gas particles will be captured.
            int mode = 1;
            int dummy;
            perform_capture(j, mode, &dummy, &dummy, DataGet, DataResult, use_part_type);
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
            
            P[place].b4.dBH_accreted_Mass += DataOut[j].Accreted_Mass;
            for (int k = 0; k < 3; k++) {
                P[place].b6.dBH_accreted_momentum[k] += DataOut[j].Accreted_Momentum[k];
                P[place].b0.dBH_accreted_angmomentum[k] += DataOut[j].Accreted_AngMomentum[k];
            }
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
    
} // capture_SPH_particles(...)


/**
 * @brief Capture gas particles and store their properties to the black hole's structure.
 * @param target
 * @param mode
 * @param nexport
 * @param nSend_local
 * @param DataGet
 * @param DataResult
 * @param use_part_type
 * @return
 **/
static int perform_capture(int target, int mode, int *nexport, int *nSend_local, struct data_in *DataGet,
                           struct data_out *DataResult, int use_part_type) {
    
    // Declare and initialise local variables. //
    double *pos, *vel, mass, bh_mass, h_i, accreted_mass, accreted_momentum[3], accreted_angmomentum[3];
    accreted_angmomentum[0] = 0.0, accreted_angmomentum[1] = 0.0, accreted_angmomentum[2] = 0.0;
    accreted_mass = 0.0, accreted_momentum[0] = 0.0, accreted_momentum[1] = 0.0, accreted_momentum[2] = 0.0;

#ifdef TIMESTEP_LIMITER_FB
    int timestep;
#endif // TIMESTEP_LIMITER_FB
    
    /* Check if this is the CPU that contains the black hole (i.e. mode == 0). If yes, extract its information directly.
     * Else, use the 'DataGet' structure to extract information sent from other CPUs */
    if (mode == 0) {
        pos = P[target].Pos;
        vel = P[target].Vel;
        h_i = P[target].Hsml;
        mass = P[target].Mass;
        bh_mass = P[target].BH_Mass;
#ifdef TIMESTEP_LIMITER_FB
        timestep = P[target].TimeStep;
#endif // TIMESTEP_LIMITER_FB
    } else {
        pos = DataGet[target].Pos;
        vel = DataGet[target].Vel;
        h_i = DataGet[target].Hsml;
        mass = DataGet[target].Mass;
        bh_mass = DataGet[target].BH_Mass;
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
                        
                        // Calculate the Hoyle & Lyttleton 1939 radius. //
                        double dvx, dvy, dvz;
                        dvx = vel[0] - P[pt0].Vel[0], dvy = vel[1] - P[pt0].Vel[1], dvz = vel[2] - P[pt0].Vel[2];
                        double relative_velocity = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
                        
                        double sound_speed =
                                sqrt(GAMMA * SphP[pt0].Entropy * pow(SphP[pt0].d.Density, GAMMA_MINUS1));
                        double r_hl = 2.0 * RA_G_CODE * bh_mass
                                      / (pow(relative_velocity, 2.0) + pow(sound_speed, 2.0));
                        
                        /* Calculate the capture radius as the maximum between the Hoyle & Lyttleton 1939 radius and
                         * the gas' softening length. */
                        double r_capture = (All.SofteningGas > r_hl) ? All.SofteningGas : r_hl;
                        
                        /* Calculate the gas particle's circular velocity based on the black hole + accretion disc mass
                         * (i.e. the dynamical mass of the black hole particle) */
                        double circular_velocity = sqrt(RA_G_CODE * mass / r_i);
                        
                        /* Check if a gas particle is within the capture radius and has relative velocity smaller than
                         * that circular velocity. If yes, add its mass in the 'accreted_mass' reservoir. */
                        if (r_i < r_capture && relative_velocity < circular_velocity) {
                            
                            // Log in 'capture.txt' the distance and Hoyle-Lyttleton radius of the captured gas particle.
                            char buf[1000]; // Declare local variable.
                            sprintf(buf, "%srelativistic_accretion_logs/capture.txt", All.OutputDir);
                            if (!(FdHoyleLyttletonAll = fopen(buf, "a"))) {
                                terminating("Cannot open 'relativistic_accretion_logs/capture.txt' file");
                            }
                            fprintf(FdHoyleLyttletonAll,
                                    " Task:%d, Time:%.20e, Sync-Point:%d, r_i:%.20e, r_hl:%.20e, r_capture:%.20e, "
                                    "sound_speed:%.20e, relative_velocity:%.20e, circular_velocity:%.20e \n", ThisTask,
                                    All.Time, All.NumCurrentTiStep, r_i, r_hl, r_capture, sound_speed,
                                    relative_velocity, circular_velocity);
                            fclose(FdHoyleLyttletonAll);
                            
                            accreted_mass += P[pt0].Mass;
                            accreted_angmomentum[0] += P[pt0].Mass * (dy * dvz - dvy * dz);
                            accreted_angmomentum[1] += P[pt0].Mass * (dvx * dz - dx * dvz);
                            accreted_angmomentum[2] += P[pt0].Mass * (dx * dvy - dvx * dy);
                            for (int k = 0; k < 3; k++) { accreted_momentum[k] += P[pt0].Mass * P[pt0].Vel[k]; }
                            
                            // Turn the captured gas particles into a ghost particles. //
                            P[pt0].Mass = 0.0;
#ifdef CS_MODEL
#ifdef CS_WSS_N
                            for (int k = 0; k < 12; k++) {P[pt0].Zm[k] = 0;}
#endif // CS_WSS_N
#endif // CS_MODEL
                        }
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
    
    /* Collect the results in the b4 and b6 unions, so you can transfer the 'accreted_mass', 'accreted_momentum', and
     * 'accreted_angmomentum' to the black hole particle and update its properties in 'RA_main.c'. */
    if (mode == 0) {
        P[target].b4.dBH_accreted_Mass = accreted_mass;
        for (int k = 0; k < 3; k++) {
            P[target].b6.dBH_accreted_momentum[k] = accreted_momentum[k];
            P[target].b0.dBH_accreted_angmomentum[k] = accreted_angmomentum[k];
        }
    } else {
        DataResult[target].Accreted_Mass = accreted_mass;
        for (int k = 0; k < 3; k++) {
            DataResult[target].Accreted_Momentum[k] = accreted_momentum[k];
            DataResult[target].Accreted_AngMomentum[k] = accreted_angmomentum[k];
        }
    }
    
    return 0;
    
}


///// Function to search particles with desired particle types within a given search center and a given radius
// This function originally comes from src/blackhole.c. The main modification is to add an additional parameter
// use_part_type. This parameter should be given as use_part_type = \sum 2^partType. For example, if we want search
// gas and dark matter particles within a given center and a given radius, we should set use_part_type = 2^0 + 2^1 = 3.
// The list of the desired particles will be given in Ngblist.
///*! This routine finds all neighbours `j' that can interact with the
//*  particle 'i' in the communication buffer.
/**
 * @brief
 * @param searchcenter
 * @param hsml
 * @param target
 * @param start_node
 * @param mode
 * @param nexport
 * @param nsend_local
 * @param use_part_type
 * @return
 **/
int find_neighbours(double searchcenter[3], double hsml, int target, int *start_node, int mode, int *nexport,
                    int *nsend_local, int use_part_type) {
    
    int number_of_ngbrs, no, p, task, nexport_save;
    struct NODE *current;
    double dx, dy, dz, dist;

#ifdef PERIODIC
    double xtmp;
#endif // PERIODIC
    
    nexport_save = *nexport;
    number_of_ngbrs = 0;
    no = *start_node;
    
    while (no >= 0) {
        if (no < All.MaxPart) {
            p = no;
            no = Nextnode[no];
            
            if (!((1 << P[p].Type) & use_part_type)) { continue; }
            
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
            if (dx > dist) { continue; }
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
            if (dy > dist) { continue; }
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
            if (dz > dist) { continue; }
            if (dx * dx + dy * dy + dz * dz > dist * dist) { continue; }
            
            Ngblist[number_of_ngbrs++] = p;
        } else {
            if (no >= All.MaxPart + MaxNodes) {
                if (mode % 2 == 1) {
                    printf("KETJU_BH, mode: %d, no: %d, All.MaxPart: %d, MaxNodes: %d \n",
                           mode, no, All.MaxPart, MaxNodes);
                    fflush(stdout);
                    endrun(12312);
                }
                
                if (Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target) {
                    Exportflag[task] = target;
                    Exportnodecount[task] = NODELISTLENGTH;
                }
                
                if (Exportnodecount[task] == NODELISTLENGTH) {
                    if (*nexport >= All.BunchSize) {
                        *nexport = nexport_save;
                        for (task = 0; task < NTask; task++) {
                            nsend_local[task] = 0;
                        }
                        for (no = 0; no < nexport_save; no++) {
                            nsend_local[DataIndexTable[no].Task]++;
                        }
                        return -1;
                    }
                    Exportnodecount[task] = 0;
                    Exportindex[task] = *nexport;
                    DataIndexTable[*nexport].Task = task;
                    DataIndexTable[*nexport].Index = target;
                    DataIndexTable[*nexport].IndexGet = *nexport;
                    *nexport = *nexport + 1;
                    nsend_local[task]++;
                }
                
                DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                
                if (Exportnodecount[task] < NODELISTLENGTH) {
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }
                
                no = Nextnode[no - MaxNodes];
                continue;
            }
            
            current = &Nodes[no];
            
            if (mode % 2 == 1) {
                if (current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)) {
                    *start_node = -1;
                    return number_of_ngbrs;
                }
            }
            
            no = current->u.d.sibling;
            
            dist = hsml + 0.5 * current->len;
            dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
            if (dx > dist) { continue; }
            dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
            if (dy > dist) { continue; }
            dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
            if (dz > dist) { continue; }
            
            dist += FACT1 * current->len;
            if (dx * dx + dy * dy + dz * dz > dist * dist) { continue; }
            
            no = current->u.d.nextnode;
        }
    }
    
    *start_node = -1;
    return number_of_ngbrs;
    
}


/**
 * @brief Calculate the Bondi-Hoyle accretion rate based on Eq. from
 * @param pt5 Index of an active black hole particle.
 * @return Bondi-Hoyle accretion rate in code units.
 **/
double bondi_hoyle_accretion_rate(int pt5) {
    
    // Calculate the black hole - gas relative velocity. //
    double dvx = P[pt5].Vel[0] - P[pt5].b3.BH_SurroundingGasVel[0];
    double dvy = P[pt5].Vel[1] - P[pt5].b3.BH_SurroundingGasVel[1];
    double dvz = P[pt5].Vel[2] - P[pt5].b3.BH_SurroundingGasVel[2];
    double relative_velocity = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
    
    // Check if co-moving integration is on. If yes, 'All.Time' is the scale factor to correct the properties. //
    if (All.ComovingIntegrationOn == 1) {
        relative_velocity /= All.Time;
        P[pt5].b1.BH_Density /= pow(All.Time, 3.0);
    }
    
    // Calculate the sound speed. //
    double sound_speed = sqrt(GAMMA * P[pt5].b2.BH_Entropy * pow(P[pt5].b1.BH_Density, GAMMA_MINUS1));
    
    // Calculate the Bondi-Hoyle accretion rate. //
    double mdot = 4.0 * RA_PI * pow(All.G, 2.0) * pow(P[pt5].BH_Mass, 2.0) * P[pt5].b1.BH_Density
                  / pow((pow(sound_speed, 2.0) + pow(relative_velocity, 2.0)), 3.0 / 2.0);
    
    // Check if the Bondi-Hoyle accretion rate is positive. If yes, return it, else terminate. //
    if (mdot > 0.0) { return mdot; }
    else {terminating("Invalid 'mdot'")}
    
} // bondi_hoyle_accretion_rate(...)


/**
 * @brief Calculate the Bondi-Hoyle timescale based on Eq. from
 * @param pt5 Index of an active black hole particle.
 * @return Bondi-Hoyle time in code units.
 **/
double bondi_hoyle_timescale(int pt5) {
    
    // Calculate the black hole-gas relative velocity. //
    double dvx = P[pt5].Vel[0] - P[pt5].b3.BH_SurroundingGasVel[0];
    double dvy = P[pt5].Vel[1] - P[pt5].b3.BH_SurroundingGasVel[1];
    double dvz = P[pt5].Vel[2] - P[pt5].b3.BH_SurroundingGasVel[2];
    double relative_velocity = sqrt(pow(dvx, 2.0) + pow(dvy, 2.0) + pow(dvz, 2.0));
    
    // Check if co-moving integration is on. If yes, 'All.Time' is the scale factor to correct the properties. //
    if (All.ComovingIntegrationOn == 1) {
        relative_velocity /= All.Time;
        P[pt5].b1.BH_Density /= pow(All.Time, 3.0);
    }
    
    // Calculate the sound speed. //
    double sound_speed = sqrt(GAMMA * P[pt5].b2.BH_Entropy * pow(P[pt5].b1.BH_Density, GAMMA_MINUS1));
    
    // Calculate the Bondi-Hoyle accretion rate. //
    double mach_number = relative_velocity / sound_speed;
    double r_bh = RA_G_CODE * P[pt5].BH_Mass / (pow(sound_speed, 2.0) * (1 + pow(mach_number, 2.0)));
    double t_bh = r_bh / (sound_speed * sqrt(1 + pow(mach_number, 2.0)));
    
    // Check if the Bondi-Hoyle timescale is positive. If yes, return it, else terminate. //
    if (t_bh > 0.0) { return t_bh; }
    else {terminating("Invalid 't_bh'")}
    
} // bondi_hoyle_timescale(...)