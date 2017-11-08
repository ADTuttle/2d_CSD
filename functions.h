#include "constants.h"
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>
#include "petscsys.h" 
//Function to initialize the state
void init(struct SimState *);
//Solve until steady state (mostly to update gating_vars)
void initialize_data(struct SimState*,struct SimState *,struct GateType*,struct ConstVars*,struct Solver*,struct FluxData*);

//Set the paramaters based on the constants
void set_params(struct SimState *,struct ConstVars*,struct GateType*,struct FluxData*);

//Linear current-voltage flux relation
void mclin(struct FluxData *,int,double,int,double,double,double,int);
//GHK Relation 
void mcGoldman(struct FluxData *,int,double,int,double,double,double,int);
//Conductance for potassium inward rectifier
double inwardrect(double,double,double);
//Returns of c_i*z_i
double cz(double *,const PetscInt *,int,int,int); 
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
PetscInt c_index(int,int,int,int);
PetscInt phi_index(int,int,int);
PetscInt al_index(int,int,int);
PetscInt xy_index(int,int);
PetscInt Ind_1(int,int,int,int);


//Create Petsc Structures
PetscErrorCode initialize_petsc(struct Solver *,int, char **);
void Get_Nonzero_in_Rows(int *);
PetscErrorCode initialize_jacobian(Mat);

//Newton Solver
PetscErrorCode newton_solve(struct SimState *,struct SimState *, double, struct GateType *, struct ExctType *, struct ConstVars *,struct Solver *,struct FluxData*); 
//Calculate residual
PetscErrorCode calc_residual(Vec,struct SimState *,struct SimState *,double,double *,double *,struct FluxData *,struct ConstVars *);
PetscErrorCode calc_jacobian(struct Solver *,struct SimState *,struct SimState *,double,double *,double *,struct FluxData *,struct ConstVars *);

//Find abs. max value of an array
double array_max(double *,size_t);
// abs. max, but for difference
double array_diff_max(double *,double *,size_t);
//Calc l2_norm of one array
double l2_norm(double *,size_t);

void print_all(double *, double *,struct ConstVars *,struct FluxData *, struct GateType *, struct SimState*, struct Solver *);



#endif