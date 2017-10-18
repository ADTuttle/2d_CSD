#include "constants.h"
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>
//Function to initialize the state
void init(struct SimState *);
//Solve until steady state (mostly to update gating_vars)
void initialize_data(struct SimState*,struct GateType*,struct ConstVars*,struct Solver*,struct FluxData*);

//Set the paramaters based on the constants
void set_params(struct SimState *,struct ConstVars*,struct GateType*,struct FluxData*);

//Linear current-voltage flux relation
void mclin(struct FluxData *,int,double,int,double,double,double,int);
//GHK Relation 
void mcGoldman(struct FluxData *,int,double,int,double,double,double,int);
//Conductance for potassium inward rectifier
double inwardrect(double,double,double);
//Returns of c_i*z_i
double cz(double *,const int *,int,int,int); 
//Computes diffusion coef
void diff_coef(double *,double *,double);
//Calculate the ion fluxes and derivatives
void ionmflux(struct FluxData *,struct SimState *,struct SimState *,struct GateType *, struct ExctType *,struct ConstVars *);
//Water flow
void wflowm(struct FluxData *,struct SimState *,struct ConstVars *);

//Update the gating variables
void gatevars_update(struct GateType *,struct SimState *,double,int);
//Initial excitation
void excitation(struct ExctType *,double);

// Indexing functions
int c_index(int,int,int,int);
int phi_index(int,int,int);
int al_index(int,int,int);
int xy_index(int,int);
int Ind_1(int,int,int,int);

//Create the sparse row and column vectors
void Assemble_Index(PetscInt *,PetscInt *);
//Create Petsc Structures
PetscErrorCode initialize_petsc(struct Solver *,int, char **);

//Newton Solver
void newton_solve(struct SimState *, double, struct GateType *, struct ExctType *, struct ConstVars *,struct Solver *,struct FluxData*); 

//Find abs. max value of an array
double array_max(double *,size_t);
// abs. max, but for difference
double array_diff_max(double *,double *,size_t);



#endif