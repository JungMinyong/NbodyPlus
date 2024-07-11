// #include <mpi.h>
#include <iostream>
#include <math.h> //for sqrt;
#include "global.h"
#include "defs.h"


Particle* FirstParticleInEnzo = nullptr;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;

void InitializeParticle(Particle* newParticle, std::vector<Particle*> &particle);
void GetCenterOfMass(std::vector<Particle*> &particle, double x_com[], double v_com[]);

void PointPotential(std::vector<Particle*> &particle, double M_point, double eps2){
	double BackgroundAcceleration[Dim], dx[Dim];
	double r2, r3;

	fprintf(stdout, "Applying Point Potential with M_point = %e Msun\n", M_point);
	M_point /= mass_unit;

	for (Particle *ptcl: particle){
		r2 = 0.0;
		for (int dim=0; dim<Dim; dim++){
			dx[dim] = ptcl->Position[dim]; //PredPosition?
			r2 += dx[dim]*dx[dim];
		}
		r2 += eps2;
		r3 = r2 * sqrt(r2);
		for (int dim=0; dim<Dim; dim++){
			ptcl->BackgroundAcceleration[dim] = - M_point * dx[dim] / r3;
		}
		
	}
}

void computeTidalTensor(std::vector<Particle*> &particles, double M_point, double eps2) {
    double r2, r5, Adot;
    double dx[Dim];
	double TidalTensor[Dim][Dim];

    for (Particle *ptcl : particles) {
        r2 = 0.0;
        for (int dim = 0; dim < Dim; dim++) {
            dx[dim] = ptcl->Position[dim];
            r2 += dx[dim] * dx[dim];
        }
        r5 = (r2 + eps2);
		r5 = r5 * r5 * sqrt(r5);

        for (int i = 0; i < Dim; i++) {
            for (int j = 0; j < Dim; j++) {
                if (i == j) {
                    TidalTensor[i][j] = M_point * (r2 - 3.0 * dx[i] * dx[j]) / r5;
                } else {
                    TidalTensor[i][j] = -3.0 * M_point * dx[i] * dx[j] / r5;
                }
            }
        }


        for (int i = 0; i < Dim; i++) {
			Adot = 0.;
            for (int j = 0; j < Dim; j++) {
            	Adot += TidalTensor[i][j] * ptcl->Velocity[j];
            }
			ptcl->BackgroundAccelerationDot[i] = Adot; //ignore the partial derivative of acceleration wrt time
        }
    }
}

int Communication_test(std::vector<Particle*> &particle) {
	double ClusterAcceleration[Dim], ClusterPosition[Dim], ClusterVelocity[Dim];
	double M_point = 1e5;  // 9.56543900e+10: point-mass of galaxy [Msun]
	double eps2;
	
	eps2 = 0.8/position_unit; //0.8 pc of softening length
	eps2 *= eps2;

	PointPotential(particle, M_point, eps2);
	#ifdef TT
	computeTidalTensor(particle, M_point, eps2);
	#endif
	
	#ifdef acc_debug
	GetCenterOfMass(particle, ClusterPosition, ClusterVelocity);
	ClusterAcceleration[0] = 0.;
	ClusterAcceleration[1] = 0.;
	ClusterAcceleration[2] = 0.;
	double total_mass=0.;
	for (Particle *ptcl: particle) {
		for (int dim=0; dim<Dim; dim++) {
			ClusterAcceleration[dim] += ptcl->Mass*ptcl->BackgroundAcceleration[dim];
		}
		total_mass += ptcl->Mass;
	}
	for (int dim=0; dim<Dim; dim++) 
			ClusterAcceleration[dim] /= total_mass;
	#endif

	return true;
}




void GetCenterOfMass(std::vector<Particle*> &particle, double x_com[], double v_com[]) {
	double total_mass=0.;
	for (int dim=0; dim<Dim; dim++) {
		x_com[dim] = 0.;
		v_com[dim] = 0.;
	}

	for (Particle *ptcl: particle) {
		for (int dim=0; dim<Dim; dim++) {
			x_com[dim] += ptcl->Mass * ptcl->Position[dim];
			v_com[dim] += ptcl->Mass * ptcl->Velocity[dim];
		}
		total_mass += ptcl->Mass;
	}

	for (int dim=0; dim<Dim; dim++) {
		x_com[dim] /= total_mass;
		v_com[dim] /= total_mass;
	}
}