#ifdef FEWBODY
#define ASSERT(x) assert(x)

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "stdio.h"
#include <vector>
#include <algorithm>
#include "../global.h"
#include "../def.h"


#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"
#include "GR_energy_loss.hpp"

// #define SEVN
#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
#endif

void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);
void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);

void NewFBInitialization3(Group* group);
void FBTermination2(Group* group);


// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp

// true: integrated normally, false: terminated by stellar merger, TDE, GW merger, etc.
// If Intererrupt_mode != none, then bin_termination = true;
void Group::ARIntegration(double next_time){

    for (int dim=0; dim<Dim; dim++) {
        sym_int.particles.cm.Position[dim] = groupCM->Position[dim];
        sym_int.particles.cm.Velocity[dim] = groupCM->Velocity[dim];
        for (int j=0; j<HERMITE_ORDER; j++)
            sym_int.particles.cm.a_irr[dim][j] = groupCM->a_irr[dim][j];
    }
    
    sym_int.particles.cm.NumberOfNeighbor = groupCM->NumberOfNeighbor;
    for (int i=0; i<groupCM->NumberOfNeighbor; i++)
        sym_int.particles.cm.Neighbors[i] = groupCM->Neighbors[i];

    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);

// /* PN corrections
    // if (PNon && bin_interrupt.status == AR::InterruptStatus::none) { // Only particles with BH
    if (bin_interrupt.status == AR::InterruptStatus::none) { // Every bound orbit
        
        auto& bin_root = sym_int.info.getBinaryTreeRoot();

        if (bin_root.getMemberN() == 2) {
            bin_root.calcOrbit(double(1.0));
            if (bin_root.ecc < 1) {
                GR_energy_loss(bin_interrupt, bin_root, CurrentTime, next_time);    
            }
        }
        else {
            GR_energy_loss_iter(bin_interrupt, bin_root, CurrentTime, next_time);
        }
        if (bin_interrupt.status == AR::InterruptStatus::none)
            sym_int.initialIntegration(next_time*EnzoTimeStep); // Eunwoo: this should be fixed later // Eunwoo: I don't think so!
    }    
// */

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        isMerger = true;

        if (sym_int.particles.getSize() == 2) {

            double pos[Dim], vel[Dim];

            groupCM->predictParticleSecondOrder(bin_interrupt.time_now/EnzoTimeStep - CurrentTime, pos, vel);
            // This might be changed later because changing Pos & Vel during Irregular Acceleration calculation is not good
            // But if SDAR integration is done after Irregular Acceleration calculation, this is fine
            // (Query) by EW 2025.1.6
            for (int dim=0; dim<Dim; dim++) {
                groupCM->Position[dim] = pos[dim];
                groupCM->Velocity[dim] = vel[dim];
            }

            // I don't use the following code anymore because it might take a long time by EW 2025.1.6
            /*
            for (int dim=0; dim<Dim; dim++) {
                groupCM->Position[dim] = pos[dim];
                groupCM->Velocity[dim] = vel[dim];
                sym_int.particles.cm.Position[dim] = pos[dim];
                sym_int.particles.cm.Velocity[dim] = vel[dim];
            }
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();
            */

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<Dim; dim++) {
                    particles[members->ParticleIndex].Position[dim] = groupCM->Position[dim] + members->Position[dim];
                    particles[members->ParticleIndex].Velocity[dim] = groupCM->Velocity[dim] + members->Velocity[dim];
                }
                particles[members->ParticleIndex].Mass = members->Mass;
            }

            isTerminate = true;

            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;

            FBTermination2(this);

            return;   
        }
        else {

            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;

            // I don't use the following code anymore because it might take a long time by EW 2025.1.6
            /*
            for (int dim=0; dim<Dim; dim++) {
                sym_int.particles.cm.Position[dim] = groupCM->Position[dim];
                sym_int.particles.cm.Velocity[dim] = groupCM->Velocity[dim];
            }
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();
            */

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<Dim; dim++) {
                    particles[members->ParticleIndex].Position[dim] = groupCM->Position[dim] + members->Position[dim];
                    particles[members->ParticleIndex].Velocity[dim] = groupCM->Velocity[dim] + members->Velocity[dim];
                }
                particles[members->ParticleIndex].Mass = members->Mass;
            }

            NewFBInitialization3(this);
            return;
        }
    }
#ifdef SEVN
    if (useSEVN) {

        bool breakEvolution = false;
        bool evolved = false;
        while (EvolutionTime + EvolutionTimeStep < next_time*EnzoTimeStep*1e4) {

            evolved = true;

            // if (groupCM->PID == -10917) {
            //     fprintf(stderr, "Stellar evolution!\n");
            //     fflush(stderr);
            // }

            EvolutionTime += EvolutionTimeStep;

            REAL dt_evolve_next = NUMERIC_FLOAT_MAX; // Myr

            for (int i=0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];
                if (members->star == nullptr || members->star->amiremnant())
                    continue;
                members->star->sync_with(EvolutionTimeStep);
                members->EvolutionTime += members->star->getp(Timestep::ID);
                // fprintf(SEVNout, "INT. PID: %d. GT: %e, ET: %e\n", members->PID, CurrentTime*EnzoTimeStep*1e4, members->EvolutionTime);
                // fflush(SEVNout);
                members->star->evolve();
                // fprintf(SEVNout, "INT. After. PID: %d, ET: %e, WT: %e\n", members->PID, members->EvolutionTime, members->star->getp(Worldtime::ID));
                // fflush(SEVNout);
                UpdateEvolution(members);
                // fprintf(SEVNout, "Int. PID: %d, Mass: %e, Radius: %e, EvolutionTime: %e\n", members->PID, members->Mass*mass_unit, members->radius*position_unit, members->EvolutionTime);
                // fflush(SEVNout);
                if (members->star->amiempty() || members->star->vkick[3] > 0.0) {
                    breakEvolution = true;
                    break;
                }
                if (members->star->amiBH())
                    PNon = true;
                if (!members->star->amiremnant() && members->star->getp(Timestep::ID) < dt_evolve_next)
                    dt_evolve_next = members->star->getp(Timestep::ID);
            }
            EvolutionTimeStep = dt_evolve_next;
            // fprintf(SEVNout, "INT. EvolutionTimeStep: %e\n", EvolutionTimeStep);
            // fflush(SEVNout);
        }

        useSEVN = false;
        for (int i=0; i < sym_int.particles.getSize(); i++) {
            Particle* members = &sym_int.particles[i];
            if (members->star != nullptr && !members->star->amiremnant()) {
                useSEVN = true;
                break;
            }
        }
        if (breakEvolution) { // break SDAR
            sym_int.particles.shiftToOriginFrame();
            sym_int.particles.template writeBackMemberAll<Particle>();

            FBTermination2(groupCM, next_time, particle);
            return;
        }
        if (evolved) { // Eunwoo: orbital parameters should be re-calculated due to mass changes during stellar evolution!
            sym_int.particles.shiftToOriginFrame();
            sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);
            sym_int.initialIntegration(next_time*EnzoTimeStep);
        }
    }
#endif

    // for write_out_group function by EW 2025.1.6
    assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
    for (int i = 0; i < sym_int.particles.getSize(); i++) {
        Particle* members = &sym_int.particles[i];

        for (int dim=0; dim<Dim; dim++) {
            particles[members->ParticleIndex].Position[dim] = groupCM->NewPosition[dim] + members->Position[dim];
            particles[members->ParticleIndex].Velocity[dim] = groupCM->NewVelocity[dim] + members->Velocity[dim];
        }
        particles[members->ParticleIndex].Mass = members->Mass;
    }

    CurrentTime = next_time;
    return; 
}
#endif