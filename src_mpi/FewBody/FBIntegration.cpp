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


#ifdef SEVN
void UpdateEvolution(Particle* ptcl);
void Mix(Star* star1, Star* star2);
#endif


void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);
void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);
void remnantSpinMass(Particle* p1, Particle* p2);
void recoilKick(Particle* p1, Particle* p2);


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


#ifdef SEVN
    bool evolved = false;
    bool kicked = false;
    for (int i=0; i < sym_int.particles.getSize(); i++) {
        Particle* members = &sym_int.particles[i];
        if (members->Mass != particles[members->ParticleIndex].Mass) {
            members->Mass = particles[members->ParticleIndex].Mass;
            evolved = true;
            if (particles[members->ParticleIndex].getBinaryInterruptState() == BinaryInterruptState::kicked) {
                kicked = true;
                break;
            }
        }
    }
    if (!kicked && evolved) { // Eunwoo: orbital parameters should be re-calculated due to mass changes during stellar evolution!
        sym_int.particles.shiftToOriginFrame();
        sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);
        sym_int.initialIntegration(next_time*EnzoTimeStep);
    }
    if (kicked) {
        for (int i = 0; i < sym_int.particles.getSize(); i++) {
            Particle* members = &sym_int.particles[i];
            particles[members->ParticleIndex].CurrentTimeIrr = CurrentTime;
        }
        isTerminate = true;
        return; 
    }
#endif

/*
    for (int i=0; i < sym_int.particles.getSize(); i++) {
		Particle* members = &sym_int.particles[i];
		fprintf(stderr, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(stderr, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(stderr, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
    }
    
    auto& bin_root = sym_int.info.getBinaryTreeRoot();
    std::cerr<<"  Binary: "
                <<"  i1="<<bin_root.getMemberIndex(0)
                <<"  i2="<<bin_root.getMemberIndex(1)
                <<"  m1="<<bin_root.m1*mass_unit
                <<"  m2="<<bin_root.m2*mass_unit
                <<"  sep="<<bin_root.r*position_unit
                <<"  semi= "<<bin_root.semi*position_unit
                <<"  ecc= "<<bin_root.ecc
                <<"  period= "<<bin_root.period*1e4
                <<"  stab= "<<bin_root.stab
                <<"  SD= "<<bin_root.slowdown.getSlowDownFactor()
                <<"  SD_org= "<<bin_root.slowdown.getSlowDownFactorOrigin()
                <<"  Tscale= "<<bin_root.slowdown.timescale
                <<"  pert_in= "<<bin_root.slowdown.pert_in
                <<"  pert_out= "<<bin_root.slowdown.pert_out;
    std::cerr<<std::endl;

    std::cerr << "Current Time: " << CurrentTime*EnzoTimeStep*1e4 << "Myr, " << "Next Time: " << next_time*EnzoTimeStep*1e4 << "Myr" << std::endl;
*/  
    // fprintf(stderr, "Before integration\n");
    // fprintf(stdout, "ARInt. current time: %e Myr, next_time: %e Myr\n", CurrentTime*EnzoTimeStep*1e4, next_time*EnzoTimeStep*1e4);
    // for (int i=0; i < sym_int.particles.getSize(); i++) {
	// 	Particle* members = &sym_int.particles[i];
    //     fprintf(stderr, "CM. Position (pc) - x:%e, y:%e, z:%e, \n", sym_int.particles.cm.Position[0]*position_unit, sym_int.particles.cm.Position[1]*position_unit, sym_int.particles.cm.Position[2]*position_unit);
	// 	fprintf(stderr, "CM. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", sym_int.particles.cm.Velocity[0]*velocity_unit/yr*pc/1e5, sym_int.particles.cm.Velocity[1]*velocity_unit/yr*pc/1e5, sym_int.particles.cm.Velocity[2]*velocity_unit/yr*pc/1e5);
	// 	fprintf(stderr, "CM. Mass (Msol) - %e, \n", sym_int.particles.cm.Mass*mass_unit);
	// 	fprintf(stderr, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
	// 	fprintf(stderr, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
	// 	fprintf(stderr, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
    // }
    
    assert(next_time > CurrentTime);
    auto bin_interrupt = sym_int.integrateToTime(next_time*EnzoTimeStep);

    // fprintf(stderr, "After integration\n");
    // fprintf(stderr, "bin_interrupt.time_now: %e Myr, bin_interrupt.time_end: %e Myr\n", bin_interrupt.time_now*1e4, bin_interrupt.time_end*1e4);
    // for (int i=0; i < sym_int.particles.getSize(); i++) {
	// 	Particle* members = &sym_int.particles[i];
    //     fprintf(stderr, "CM. Position (pc) - x:%e, y:%e, z:%e, \n", sym_int.particles.cm.Position[0]*position_unit, sym_int.particles.cm.Position[1]*position_unit, sym_int.particles.cm.Position[2]*position_unit);
	// 	fprintf(stderr, "CM. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", sym_int.particles.cm.Velocity[0]*velocity_unit/yr*pc/1e5, sym_int.particles.cm.Velocity[1]*velocity_unit/yr*pc/1e5, sym_int.particles.cm.Velocity[2]*velocity_unit/yr*pc/1e5);
	// 	fprintf(stderr, "CM. Mass (Msol) - %e, \n", sym_int.particles.cm.Mass*mass_unit);
	// 	fprintf(stderr, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
	// 	fprintf(stderr, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
	// 	fprintf(stderr, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
    // }
/*
    for (int i=0; i < sym_int.particles.getSize(); i++) {
		Particle* members = &sym_int.particles[i];
		fprintf(stderr, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(stderr, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(stderr, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
    }

    bin_root = sym_int.info.getBinaryTreeRoot();
    std::cerr<<"  Binary: "
                <<"  i1="<<bin_root.getMemberIndex(0)
                <<"  i2="<<bin_root.getMemberIndex(1)
                <<"  m1="<<bin_root.m1*mass_unit
                <<"  m2="<<bin_root.m2*mass_unit
                <<"  sep="<<bin_root.r*position_unit
                <<"  semi= "<<bin_root.semi*position_unit
                <<"  ecc= "<<bin_root.ecc
                <<"  period= "<<bin_root.period*1e4
                <<"  stab= "<<bin_root.stab
                <<"  SD= "<<bin_root.slowdown.getSlowDownFactor()
                <<"  SD_org= "<<bin_root.slowdown.getSlowDownFactorOrigin()
                <<"  Tscale= "<<bin_root.slowdown.timescale
                <<"  pert_in= "<<bin_root.slowdown.pert_in
                <<"  pert_out= "<<bin_root.slowdown.pert_out;
    std::cerr<<std::endl;
*/
// /* PN corrections
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
        groupCM->setBinaryInterruptState(BinaryInterruptState::merger);

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
            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;
            groupCM->CurrentTimeIrr = CurrentTime;

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<Dim; dim++) {
                    particles[members->ParticleIndex].Position[dim] = groupCM->Position[dim] + members->Position[dim];
                    particles[members->ParticleIndex].Velocity[dim] = groupCM->Velocity[dim] + members->Velocity[dim];
                }
                particles[members->ParticleIndex].Mass = members->Mass;
                particles[members->ParticleIndex].binary_state = members->binary_state;
                particles[members->ParticleIndex].CurrentTimeIrr = CurrentTime;
            }

            isTerminate = true;

            return;   
        }
        else {

            CurrentTime = bin_interrupt.time_now/EnzoTimeStep;

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<Dim; dim++) {
                    particles[members->ParticleIndex].Position[dim] = groupCM->Position[dim] + members->Position[dim];
                    particles[members->ParticleIndex].Velocity[dim] = groupCM->Velocity[dim] + members->Velocity[dim];
                }
                particles[members->ParticleIndex].Mass = members->Mass;
                particles[members->ParticleIndex].binary_state = members->binary_state;
            }

            // NewFBInitialization3(this);
            return;
        }
    }

    // for write_out_group function by EW 2025.1.6
    assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
    for (int i = 0; i < sym_int.particles.getSize(); i++) {
        Particle* members = &sym_int.particles[i];

        for (int dim=0; dim<Dim; dim++) {
            particles[members->ParticleIndex].Position[dim] = groupCM->NewPosition[dim] + members->Position[dim];
            particles[members->ParticleIndex].Velocity[dim] = groupCM->NewVelocity[dim] + members->Velocity[dim];
        }
        particles[members->ParticleIndex].Mass = members->Mass;
        particles[members->ParticleIndex].CurrentTimeIrr = next_time;
    }
    
    CurrentTime = next_time;
    return; 
}


void Merge(Particle* p1, Particle* p2) { // Stellar merger

    if (p1->Mass < p2->Mass) 
        std::swap(p1, p2); // p1 should have the larger mass than p2 (p1->Mass > p2->Mass)

    p1->setBinaryInterruptState(BinaryInterruptState::none);
    p2->setBinaryInterruptState(BinaryInterruptState::none);

    double radius;

    if (p1->ParticleType == (Blackhole+SingleStar) && p2->ParticleType == (Blackhole+SingleStar)) {

        radius = (p1->radius > p2->radius) ? 3*p1->radius : 3*p2->radius; // r_ISCO == 3 * Schwartzschild radius
        fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
        // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(mergerout, "r_ISCO: %e pc\n", radius*position_unit);

        fprintf(mergerout, "GW driven merger happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
        fprintf(mergerout, "Time: %e Myr\n", p1->CurrentTimeIrr*EnzoTimeStep*1e4);
        // fprintf(mergerout, "In center-of-mass frame...\n");
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
        fprintf(mergerout, "PID: %d. Dimensionless spin - %e, %e, %e\n", p1->PID, p1->a_spin[0], p1->a_spin[1], p1->a_spin[2]);
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);
        fprintf(mergerout, "PID: %d. Dimensionless spin - %e, %e, %e\n", p2->PID, p2->a_spin[0], p2->a_spin[1], p2->a_spin[2]);


        double mcm = p1->Mass + p2->Mass;
        for (int k=0; k<3; k++) {
            p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
            p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
            p2->Position[k] = 0.0;
            p2->Velocity[k] = 0.0;
        }
        recoilKick(p1, p2);
        remnantSpinMass(p1, p2);
        p1->radius = 2*p1->Mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius

        p2->Mass = 0.0;
        fprintf(mergerout, "---------------Merger remnant properties---------------\n");
        fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
        fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
    }
    else if ((p1->ParticleType == (Blackhole+SingleStar) && p2->ParticleType == (NormalStar+SingleStar)) ||
            (p1->ParticleType == (NormalStar+SingleStar) && p2->ParticleType == (Blackhole+SingleStar))) {

        if (p2->ParticleType == (Blackhole+SingleStar)) 
            std::swap(p1, p2); // p1 should be BH

        radius = 1.3*pow((p1->Mass + p2->Mass)/p2->Mass, 1./3)*p2->radius; // TDE radius

        fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
        // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(mergerout, "r_TDE: %e pc\n", radius*position_unit);

        fprintf(mergerout, "TDE happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
        fprintf(mergerout, "Time: %e Myr\n", p1->CurrentTimeIrr*EnzoTimeStep*1e4);
        // fprintf(mergerout, "In center-of-mass frame...\n");
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);


        double mcm = p1->Mass + p2->Mass * 0.5; // If TDE happens, the half of the mass of star is accreted to a BH.
        for (int k=0; k<3; k++) {
            p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
            p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
        }

        p1->radius = 2*p1->Mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius

        p1->dm += 0.5 * p2->Mass;
        p1->Mass = mcm;
        p2->Mass = 0.0;
        fprintf(mergerout, "---------------Merger remnant properties---------------\n");
        fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
        fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
    }
    else if (p1->ParticleType == (NormalStar+SingleStar) && p2->ParticleType == (NormalStar+SingleStar)) {

        radius = p1->radius + p2->radius; // Sum of two stellar radius
        // fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
        // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(mergerout, "r1 + r2: %e pc\n", radius*position_unit);

        fprintf(mergerout, "Stellar merger happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
        fprintf(mergerout, "Time: %e Myr\n", p1->CurrentTimeIrr*EnzoTimeStep*1e4);
        // fprintf(mergerout, "In center-of-mass frame...\n");
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
        fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
        fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);

        double mcm = p1->Mass + p2->Mass;
        for (int k=0; k<3; k++) {
            p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
            p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
            p2->Position[k] = p1->Position[k];
            p2->Velocity[k] = p1->Velocity[k];
            
        }

#ifdef SEVN
        if (p1->StellarEvolution == nullptr && p2->StellarEvolution == nullptr) {

            // p1->dm = mcm - p1->Mass;
            // p2->dm = -p2->Mass;
            p1->Mass = mcm;
            p2->Mass = 0.0;

            p1->radius = 2.25461e-8/position_unit*pow(p1->Mass*mass_unit, 1./3); // stellar radius in code unit

            if (mcm*mass_unit > 2.2 && mcm*mass_unit < 600) {

                std::vector<std::string> args = {"empty", // Not used
                                            // "-myself", "/data/vinicius/NbodyPlus/SEVN",
                                            "-tables", "/data/vinicius/mpi/NbodyPlues/SEVN/tables/SEVNtracks_parsec_ov04_AGB", 
                                            //  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
                                            // "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
                                            // "-turn_WR_to_pureHe", "false",
                                            "-snmode", "delayed",
                                            "-Z", "0.0002",
                                            "-spin", "0.0",
                                            "-tini", "zams", 
                                            "-tf", "end",
                                            // "-tf", "0.000122",
                                            "-dtout", "events",
                                            "-xspinmode", "geneva"};
                std::vector<char*> c_args;
                for (auto& arg : args) {
                    c_args.push_back(&arg[0]);
                }

                IO* sevnio; // Eunwoo: global variable -> We can initialize Star and Binstar class anywhere.
                sevnio = new IO;
                sevnio->load(c_args.size(), c_args.data());

                std::vector<std::string> init_params{std::to_string(double(p1->Mass*mass_unit)), "0.0002", "0.0", "delayed", "zams", "end", "events"};
                size_t id = p1->PID;
                p1->StellarEvolution = new Star(sevnio, init_params, id, false);

                p1->FormationTime = p1->CurrentTimeIrr*EnzoTimeStep*1e4;
                p1->WorldTime = p1->CurrentTimeIrr*EnzoTimeStep*1e4;
                p1->Mzams = p1->StellarEvolution->get_zams();
                p1->radius = p1->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
                fprintf(SEVNout, "New Star class made!\n");
                fprintf(SEVNout, "PID: %d. Mass: %e Msol, Radius: %e pc\n", p1->PID, p1->Mass*mass_unit, p1->radius*position_unit);
            }
            fprintf(mergerout, "---------------Merger remnant properties---------------\n");
            fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
            // fprintf(mergerout, "Type: %d, \n", int(p1->StellarEvolution->getp(Phase::ID)));
            fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if (p1->StellarEvolution != nullptr && p2->StellarEvolution != nullptr) {

            Mix(p1->StellarEvolution, p2->StellarEvolution);
            UpdateEvolution(p1);
            UpdateEvolution(p2);
            fprintf(SEVNout, "Mix done!\n");

            if (p1->StellarEvolution->amiempty() && !p2->StellarEvolution->amiempty()) {
                p1->Mass = 0.0;
                p2->Mzams = p1->Mzams + p2->Mzams;
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
                fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(mergerout, "Mass (Msol) - %e, \n", p2->Mass*mass_unit);
                // if (p2->StellarEvolution->amiremnant())
                //     fprintf(mergerout, "Type: %d, \n", 8 + int(p2->StellarEvolution->getp(RemnantType::ID)));
                // else
                //     fprintf(mergerout, "Type: %d, \n", int(p2->StellarEvolution->getp(Phase::ID)));
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (!p1->StellarEvolution->amiempty() && p2->StellarEvolution->amiempty()) {
                p2->Mass = 0.0;
                p1->Mzams = p1->Mzams + p2->Mzams;
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
                // if (p1->StellarEvolution->amiremnant())
                //     fprintf(mergerout, "Type: %d, \n", 8 + int(p1->StellarEvolution->getp(RemnantType::ID)));
                // else
                //     fprintf(mergerout, "Type: %d, \n", int(p1->StellarEvolution->getp(Phase::ID)));
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (p1->StellarEvolution->amiempty() && p2->StellarEvolution->amiempty()) { // Type Ia supernova
                p1->dm += p1->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
                p1->Mass = 0.0;
                p2->dm += p2->Mass;
                p2->Mass = 0.0;
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Type Ia Supernova event! Both of the stars becomes empty!\n");
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else
                throw std::runtime_error("None of stars are empty: Something wrong in stellar merger!");
        }
        else if (p1->StellarEvolution == nullptr && p2->StellarEvolution != nullptr) {

            p2->StellarEvolution->update_from_binary(Mass::ID, p1->Mass*mass_unit);
            p2->StellarEvolution->update_from_binary(dMcumul_binary::ID, p1->Mass*mass_unit);
            if (p2->StellarEvolution->aminakedhelium())
                p2->StellarEvolution->jump_to_normal_tracks();
            else
                p2->StellarEvolution->find_new_track_after_merger();

            UpdateEvolution(p2);
            p1->Mass = 0.0;
            p2->Mzams = p1->Mzams + p2->Mzams;

            fprintf(SEVNout, "Mix with no done!\n");
            fprintf(mergerout, "---------------Merger remnant properties---------------\n");
            fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
            fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "Mass (Msol) - %e, \n", p2->Mass*mass_unit);
            // if (p2->StellarEvolution->amiremnant())
            //     fprintf(mergerout, "Type: %d, \n", 8 + int(p2->StellarEvolution->getp(RemnantType::ID)));
            // else
            //     fprintf(mergerout, "Type: %d, \n", int(p2->StellarEvolution->getp(Phase::ID)));
            fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if (p1->StellarEvolution != nullptr && p2->StellarEvolution == nullptr) {

            p1->StellarEvolution->update_from_binary(Mass::ID, p2->Mass*mass_unit);
            p1->StellarEvolution->update_from_binary(dMcumul_binary::ID, p2->Mass*mass_unit);
            if (p1->StellarEvolution->aminakedhelium())
                p1->StellarEvolution->jump_to_normal_tracks();
            else
                p1->StellarEvolution->find_new_track_after_merger();

            UpdateEvolution(p1);
            p2->Mass = 0.0;
            p1->Mzams = p1->Mzams + p2->Mzams;

            fprintf(SEVNout, "Mix with no done!\n");
            fprintf(mergerout, "---------------Merger remnant properties---------------\n");
            fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
            fflush(mergerout);
            // if (p1->StellarEvolution->amiremnant())
            //     fprintf(mergerout, "Type: %d, \n", 8 + int(p1->StellarEvolution->getp(RemnantType::ID)));
            // else
            //     fprintf(mergerout, "Type: %d, \n", int(p1->StellarEvolution->getp(Phase::ID)));
            fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        }
    }
    fflush(mergerout);
    fflush(SEVNout);
#else
        p1->radius = 2.25461e-8/position_unit*pow(p1->Mass*mass_unit, 1./3); // stellar radius in code unit

        p1->dm = mcm - p1->Mass;
        p2->dm = -p2->Mass;
        p1->Mass = mcm;
        p2->Mass = 0.0;
    }
    fflush(mergerout);   
#endif
}

// Reference for remnant spin: Hofmann et al. (2016) (https://iopscience.iop.org/article/10.3847/2041-8205/825/2/L19/pdf)
// Using eq (2)-(6), (13)-(16)
// n_M = 3, n_J = 4
// Reference for remnant mass: Barausse et al. (2012) (https://iopscience.iop.org/article/10.1088/0004-637X/758/1/63/pdf)
// Using eq (1)-(5), (12), (15)-(18)
void remnantSpinMass(Particle* p1, Particle* p2) {

    double k[4][5] = {{-5.9, 3.39221, 4.48865, -5.77101, -13.0459},
                    {35.1287, -72.9336, -86.0036, 93.7371, 200.975},
                    {-146.822, 387.184, 447.009, -467.383, -884.339},
                    {223.911, -648.502, -697.177, 753.738, 1166.89}};
    double ksi = 0.474046;

    double a1 = sqrt(mag(p1->a_spin));
    double a2 = sqrt(mag(p2->a_spin));
    double Mtot = p1->Mass + p2->Mass;
    double q = p2->Mass/p1->Mass;
    assert(q <= 1);
    double nu = q/(1+q)/(1+q);
    const double c = 299752.458 / (velocity_unit / yr * pc / 1e5);

    double pos_rel[3];
    double vel_rel[3];
    double L_ang[3]; // specific angular momentum (r_rel x v_rel)

    for (int dim=0; dim<Dim; dim++) {
        pos_rel[dim] = p1->Position[dim] - p2->Position[dim];
        vel_rel[dim] = p1->Velocity[dim] - p2->Velocity[dim];
    }
    L_ang[0] = pos_rel[1] * vel_rel[2] - pos_rel[2] * vel_rel[1];
    L_ang[1] = pos_rel[2] * vel_rel[0] - pos_rel[0] * vel_rel[2];
    L_ang[2] = pos_rel[0] * vel_rel[1] - pos_rel[1] * vel_rel[0];

    fprintf(mergerout, "L_orbit: (%e, %e, %e)\n", (p1->Mass*p2->Mass/Mtot) * L_ang[0],
                                                    (p1->Mass*p2->Mass/Mtot) * L_ang[1],
                                                    (p1->Mass*p2->Mass/Mtot) * L_ang[2]);
    fprintf(mergerout, "S_1: (%e, %e, %e)\n", (p1->Mass*p1->Mass/c) * p1->a_spin[0], 
                                                (p1->Mass*p1->Mass/c) * p1->a_spin[1], 
                                                (p1->Mass*p1->Mass/c) * p1->a_spin[2]);
    fprintf(mergerout, "S_2: (%e, %e, %e)\n", (p2->Mass*p2->Mass/c) * p2->a_spin[0], 
                                                (p2->Mass*p2->Mass/c) * p2->a_spin[1], 
                                                (p2->Mass*p2->Mass/c) * p2->a_spin[2]);
    fprintf(mergerout, "In code unit!\n");

    auto cosine = [&](double a[3], double b[3]) {
        if (mag(a) == 0 || mag(b) == 0)
            return 0.;
        else
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/sqrt(mag(a))/sqrt(mag(b));
    };

    double cosa = cosine(p1->a_spin, p2->a_spin); // cos(alpha): This angle should be changed to the initial value. I will change this later.
    double cosb = cosine(L_ang, p1->a_spin);      // cos(beta)
    double cosg = cosine(L_ang, p2->a_spin);      // cos(gamma)

    double atot = (a1*cosb + a2*cosg*q*q)/(1 + q)/(1 + q); // atilde in Barausse et al. (2012)
    double aeff = atot + ksi*nu*(a1*cosb + a2*cosg);

    auto Z1 = [&](double a) {
        return 1 + pow(1 - a*a, 1./3)*(pow(1 + a, 1./3) + pow(1 - a, 1./3));
    };
    auto Z2 = [&](double a, double Z1) {
        return sqrt(3*a*a + Z1*Z1);
    };
    auto r_ISCO = [&](double a, double Z1, double Z2) {
        if (a > 0)
            return 3 + Z2 - sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
        else if (a < 0)
            return 3 + Z2 + sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
        else
            return 3 + Z2;
    };
    auto E_ISCO = [&](double r) {
        return sqrt(1 - 2./3/r);
    };

    double Z1_aeff = Z1(aeff);
    double Z2_aeff = Z2(aeff, Z1_aeff);
    double r_ISCO_aeff = r_ISCO(aeff, Z1_aeff, Z2_aeff);
    double E_ISCO_aeff = E_ISCO(r_ISCO_aeff);
    double L_ISCO_aeff = 2/3/sqrt(3)*(1 + 2*sqrt(3*r_ISCO_aeff - 2));

    double l = L_ISCO_aeff - 2 * atot * (E_ISCO_aeff - 1);
    for (int i = 0; i <= 3; ++i) {      // n_M = 3
        for (int j = 0; j <= 4; ++j) {  // n_J = 4
            l += k[i][j] * pow(nu, 1 + i) * pow(aeff, j);
        }
    }
    l = abs(l);

    double afin = 1./pow(1 + q, 2)*sqrt(a1*a1 + a2*a2*pow(q, 4) + 2*a1*a2*q*q*cosa + 2*(a1*cosb + a2*q*q*cosg)*l*q + l*l*q*q);

    double J_tot[3];
    for (int i=0; i<3; i++)
        J_tot[i] = (p1->Mass*p2->Mass/Mtot) * L_ang[i] + (p1->Mass*p1->Mass/c) * p1->a_spin[i] + (p2->Mass*p2->Mass/c) * p2->a_spin[i];

    double J_tot_norm[3];
    for (int i=0; i<3; i++)
        J_tot_norm[i] = J_tot[i]/sqrt(mag(J_tot));

    for (int i=0; i<3; i++)
        p1->a_spin[i] = afin*J_tot_norm[i];

    fprintf(mergerout, "Dimensionless spin of remnant BH: (%e, %e, %e)\n", p1->a_spin[0], p1->a_spin[1], p1->a_spin[2]);

    double Z1_atot = Z1(atot);
    double Z2_atot = Z2(atot, Z1_atot);
    double r_ISCO_atot = r_ISCO(atot, Z1_atot, Z2_atot);
    double E_ISCO_atot = E_ISCO(r_ISCO_atot);

    double Erad = (1 - E_ISCO_atot)*nu + 4*nu*nu*(4*0.04827 + 16*0.01707*atot*(atot+1) + E_ISCO_atot - 1);

    p1->Mass = (1 - Erad) * Mtot;

    fprintf(mergerout, "Mass of remnant BH: %e Msol\n", p1->Mass*mass_unit);
    fflush(mergerout);
}

// Reference: Arca Sedda et al. (2020) (https://iopscience.iop.org/article/10.3847/1538-4357/ab88b2/pdf)
void recoilKick(Particle* p1, Particle* p2) {
    
    double A      = 1.2e4; // km/s
    double B      = -0.93;
    double H      = 6.9e3; // km/s
    double ksi    = 145 * M_PI / 180; // rad
    double V11    = 3677.76; // km/s
    double VA     = 2.481e3; // km/s
    double VB     = 1.793e3; // km/s
    double VC     = 1.507e3; // km/s

    auto cosine = [&](double a[3], double b[3]) {
        if (mag(a) == 0 || mag(b) == 0)
            return 0.;
        else
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/sqrt(mag(a))/sqrt(mag(b));
    };

    double a1 = sqrt(mag(p1->a_spin));
    double a2 = sqrt(mag(p2->a_spin));
    double q = p2->Mass/p1->Mass;
    // fprintf(mergerout, "q: %e\n", q);
    assert(q <= 1);
    double nu = q/(1+q)/(1+q);
    // fprintf(mergerout, "nu: %e\n", nu);

    double pos_rel[3];
    double vel_rel[3];
    double L_ang[3]; // specific angular momentum (r_rel x v_rel)

    for (int dim=0; dim<Dim; dim++) {
        pos_rel[dim] = p1->Position[dim] - p2->Position[dim];
        vel_rel[dim] = p1->Velocity[dim] - p2->Velocity[dim];
    }
    L_ang[0] = pos_rel[1] * vel_rel[2] - pos_rel[2] * vel_rel[1];
    L_ang[1] = pos_rel[2] * vel_rel[0] - pos_rel[0] * vel_rel[2];
    L_ang[2] = pos_rel[0] * vel_rel[1] - pos_rel[1] * vel_rel[0];


    double cosa = cosine(p1->a_spin, p2->a_spin); // cos(alpha): This angle should be changed to the initial value. I will change this later.
    double sina = sin(acos(cosa));
    double cosb = cosine(L_ang, p1->a_spin);      // cos(beta)
    double sinb = sin(acos(cosb));
    double cosg = cosine(L_ang, p2->a_spin);      // cos(gamma)
    double sing = sin(acos(cosg));
    // fprintf(mergerout, "cosa: %e, sina: %e, cosb: %e, sinb: %e, cosg: %e, sing: %e\n", cosa, sina, cosb, sinb, cosg, sing);

    double a2par  = a2 * cosg;
    // fprintf(mergerout, "a2par: %e\n", a2par);
    double a2per1 = a2 * sing;
    // fprintf(mergerout, "a2per1, %e\n", a2per1);
    double a2per2 = 0.0;
    // fprintf(mergerout, "a2per2, %e\n", a2per2);
    
    double a1par  = a1 * cosb;
    // fprintf(mergerout, "a1par: %e\n", a1par);
    double a1per1 = a1 * sinb*cosa;
    // fprintf(mergerout, "a1per1: %e\n", a1per1);
    double a1per2 = a1 * sinb*sina;
    // fprintf(mergerout, "a1per2: %e\n", a1per2);

    double KSIpar = 2 * (a2par + q*q*a1par) / (1 + q) / (1 + q);
    // fprintf(mergerout, "KSIpar: %e\n", KSIpar);
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    double phi = 2 * M_PI * distr(mt); // phi_Delta - phi_1
    // fprintf(mergerout, "phi: %e\n", phi);

    double vm = A*nu*nu*sqrt(1 - 4*nu) * (1 + B*nu);
    // fprintf(mergerout, "vm: %e\n", vm);
    double vper = H*nu*nu / (1 + q) * (a2par - q * a1par);
    // fprintf(mergerout, "vper: %e\n", vper);
    double vpar = 16*nu*nu / (1 + q) * (V11 + VA*KSIpar + VB*KSIpar*KSIpar + VC*KSIpar*KSIpar*KSIpar);
    vpar *= sqrt((a2per1 - q * a1per1) * (a2per1 - q * a1per1) + (a2per2 - q * a1per2) * (a2per2 - q * a1per2)) * cos(phi);
    // fprintf(mergerout, "vpar: %e\n", vpar);

    double e_par[3]   = {L_ang[0]/sqrt(mag(L_ang)), L_ang[1]/sqrt(mag(L_ang)), L_ang[2]/sqrt(mag(L_ang))};
    double e_per1[3]  = {p2->a_spin[0] - p2->a_spin[0]*cosg, p2->a_spin[1] - p2->a_spin[1]*cosg, p2->a_spin[2] - p2->a_spin[2]*cosg};
    double norm1      = sqrt(mag(e_per1));
    if (norm1 != 0) {
        for (int i=0; i<3; i++)
            e_per1[i]   /= norm1; 
    }
            
    double e_per2[3]  =   {e_par[1] * e_per1[2] - e_par[2] * e_per1[1], 
                        e_par[2] * e_per1[0] - e_par[0] * e_per1[2], 
                        e_par[0] * e_per1[1] - e_par[1] * e_per1[0]}; // cross product: e2 = e3 x e1
    // fprintf(mergerout, "e1: %e, %e, %e\n", e_per1[0], e_per1[1], e_per1[2]);
    // fprintf(mergerout, "e2: %e, %e, %e\n", e_per2[0], e_per2[1], e_per2[2]);
    // fprintf(mergerout, "e3: %e, %e, %e\n", e_par[0], e_par[1], e_par[2]);

    double vkick[3];
    for (int i=0; i<3; i++)
        vkick[i] = (vm + vper*cos(ksi)) * e_per1[i] + vper*sin(ksi) * e_per2[i] + vpar * e_par[i];

    fprintf(mergerout, "GW recoil kick: (%e, %e, %e) km/s\n", vkick[0], vkick[1], vkick[2]);
    fprintf(mergerout, "\t magnitude: %e km/s\n", sqrt(mag(vkick)));

    for (int i=0; i<3; i++)
        p1->Velocity[i] += vkick[i]/(velocity_unit/yr*pc/1e5); // km/s to code unit

    // fprintf(mergerout, "Remnant velocity: (%e, %e, %e) km/s\n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
    fflush(mergerout);
}



#endif