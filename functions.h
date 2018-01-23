#include "constants.h"
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <petsctime.h>
#include <petscsnes.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petsclog.h>
//Function to initialize the state
void init(Vec,struct SimState *);
//Solve until steady state (mostly to update gating_vars)
void initialize_data(Vec,struct AppCtx*);

//Set the paramaters based on the constants
void set_params(Vec,struct SimState *,struct ConstVars*,struct GateType*,struct FluxData*);

//Data management functions
PetscErrorCode init_simstate(Vec,struct SimState*);
PetscErrorCode extract_subarray(Vec,struct SimState*);
PetscErrorCode restore_subarray(Vec,struct SimState*);
PetscErrorCode copy_simstate(Vec,struct SimState *state_vars_past);

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
//Update Volume fraction
void volume_update(struct SimState*,struct SimState*,struct AppCtx*);
//Initial excitation
void excitation(struct ExctType *,PetscReal);

// Indexing functions
PetscInt c_index(PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt phi_index(PetscInt,PetscInt,PetscInt);
PetscInt al_index(PetscInt,PetscInt,PetscInt);
PetscInt xy_index(PetscInt,PetscInt);
PetscInt Ind_1(PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt Ind_nx(PetscInt,PetscInt ,PetscInt ,PetscInt , PetscInt );


//Create Petsc Structures

PetscErrorCode initialize_petsc(struct Solver *,int, char **,struct AppCtx*);
void Get_Nonzero_in_Rows(int *);
PetscErrorCode initialize_jacobian(Mat);

//Multigrid functions
PetscErrorCode Create_Restriction(Mat R,PetscInt nx, PetscInt ny);
PetscErrorCode Create_Interpolation(Mat R,PetscInt nx, PetscInt ny);
PetscErrorCode Initialize_PCMG(PC pc,Mat A);

//Newton Solver
PetscErrorCode newton_solve(Vec,struct AppCtx*);
//Calculate residual
PetscErrorCode calc_residual(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
PetscErrorCode calc_residual_no_vol(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_no_vol(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx

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
void write_data(FILE *,struct SimState *,int );
void init_events(struct AppCtx *);


#endif