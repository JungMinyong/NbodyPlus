// #include <mpi.h>
#include <iostream>
#include <math.h> //for sqrt;
#include "global.h"
#include "defs.h"


Particle* FirstParticleInEnzo = nullptr;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;

void InitializeParticle(Particle* newParticle, std::vector<Particle*> &particle);

void PointPotential(std::vector<Particle*> &particle, double M_point){
	double BackgroundAcceleration[Dim], dx[Dim];
	double r2, r3;

	M_point /= mass_unit;
	for (Particle *ptcl: particle){
		r2 = 0;
		for (int dim=0; dim<Dim; dim++){
		dx[dim] = ptcl->PredPosition[dim];
		r2 += dx[dim]*dx[dim];
		}
		r3 = r2 * std::sqrt(r2);
		for (int dim=0; dim<Dim; dim++){
		// BackgroundAcceleration[dim] = - dx[dim] / r3;
		ptcl->BackgroundAcceleration[dim] = - dx[dim] / r3;
		}
		
	}
}
int Communication_test(std::vector<Particle*> &particle) {
	double M_point = 1000; //1000Msun BH
	PointPotential(particle, M_point);
	return true;
}



