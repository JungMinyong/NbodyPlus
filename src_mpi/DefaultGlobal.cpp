#include <iostream>
#include <stdio.h>
#include "global.h"
#include "GlobalVariable.h"
#include <mpi.h>



Particle *particles_original;
Particle *particles;

MPI_Win win;
MPI_Win win2;

MPI_Comm shared_comm;
int MyRank;
int NumberOfProcessor;
int NumberOfWorker;
int NumberOfCommunication;


int LoadBalanceParticle;


GlobalVariable *global_variable;
GlobalVariable *global_variable_original;
int NumberOfParticle;
int NewPID;

// Task
int Task[NumberOfTask];

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
int time_block;
double time_step;
ULL block_max;
double outputTimeStep;
double endTime;

double binary_time;
double binary_time_prev;
ULL binary_block;

// Enzo to Nbody
Particle* FirstEnzoParticle;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoTimeStep;



// i/o
char* fname;
double inputTime;
bool restart;
char* foutput;
bool IsOutput;
double outputTime;
int outNum;

FILE* binout;
FILE* mergerout;

#ifdef PerformanceTrace
Performance performance;
#endif

void DefaultGlobal() {

	/* Task initialization */
	//int Task[NumberOfTask];
	for (int i=0;i<NumberOfTask; i++) {
		Task[i] = i;
	}

	// Eunwoo: I think these lines are redundant
	/* Number of particles */
	// NumberOfParticle = 1000;
	// NumberOfSingle = 1000;
	// Eunwoo: I think these lines are redundant

	NumberOfCommunication = 0;

	LoadBalanceParticle = 10;

	/* Timesteps */
	endTime = 1;
	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	time_block = -30;
	block_max = static_cast<ULL>(pow(2, -time_block));
	time_step = std::pow(2,time_block);

	inputTime = 0.0;
	endTime = 0.0;
	outputTimeStep = 0.;

}
