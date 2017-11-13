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
double cz(const double *,const PetscInt *,int,int,int);
//Computes diffusion coef
void diff_coef(double *,const PetscReal *,PetscReal);
//Calculate the ion fluxes and derivatives
void ionmflux(struct FluxData *,struct SimState *,struct SimState *,struct GateType *, struct ExctType *,struct ConstVars *);
//Water flow
void wflowm(struct FluxData *,struct SimState *,struct ConstVars *);

//Update the gating variables
void gatevars_update(struct GateType *,struct SimState *,PetscReal,int);
//Initial excitation
void excitation(struct ExctType *,PetscReal);

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
PetscErrorCode newton_solve(struct SimState *,struct SimState *, PetscReal, struct GateType *, struct ExctType *, struct ConstVars *,struct Solver *,struct FluxData*);
//Calculate residual
PetscErrorCode calc_residual(Vec,struct SimState *,struct SimState *,PetscReal,PetscReal *,PetscReal *,struct FluxData *,struct ConstVars *);
PetscErrorCode calc_jacobian(Mat, struct SimState *, struct SimState *, PetscReal, PetscReal *, PetscReal *, struct FluxData *,
                             struct ConstVars *, int iter);

//Find abs. max value of an array
PetscReal array_max(PetscReal *,size_t);
// abs. max, but for difference
PetscReal array_diff_max(PetscReal *,PetscReal *,size_t);
//Calc l2_norm of one array
PetscReal l2_norm(PetscReal *,size_t);

void print_all(PetscReal *, PetscReal *,struct ConstVars *,struct FluxData *, struct GateType *, struct SimState*, struct Solver *);
const char* getfield(char* , int );
void find_print(int, int, double, int iter);
void compare_res(double *, int );


#endif