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
void init(Vec,struct SimState *,struct AppCtx*);
//Solve until steady state (mostly to update gating_vars)
void initialize_data(Vec,struct AppCtx*);

//Set the paramaters based on the constants
void set_params(Vec,struct SimState *,struct ConstVars*,struct GateType*,struct FluxData*,struct AppCtx*);

//Data management functions
void init_arrays(struct AppCtx*);
PetscErrorCode init_simstate(Vec,struct SimState*,struct AppCtx*);
PetscErrorCode extract_subarray(Vec,struct SimState*, int);
PetscErrorCode restore_subarray(Vec,struct SimState*, int);
PetscErrorCode copy_simstate(Vec,struct SimState *state_vars_past,int);

//Linear current-voltage flux relation
void mclin(struct FluxData *,int,double,int,double,double,double,int);
//GHK Relation 
void mcGoldman(struct FluxData *,int,double,int,double,double,double,int);
//Conductance for potassium inward rectifier
double inwardrect(double,double,double);
//Returns of c_i*z_i
double cz(const double *,const PetscInt *,int,int,int,struct AppCtx*);
//Computes diffusion coef
void diff_coef(double *,const PetscReal *,PetscReal,struct AppCtx*);
//Calculate the ion fluxes and derivatives
void ionmflux(struct AppCtx*);
//Water flow
void wflowm(struct AppCtx*);

//Update the gating variables
void gatevars_update(struct GateType *,struct SimState *,PetscReal,struct AppCtx *,int);
//Update Volume fraction
void volume_update(struct SimState*,struct SimState*,struct AppCtx*);
//Initial excitation
void excitation(struct AppCtx*,PetscReal);

// Indexing functions
PetscInt c_index(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt phi_index(PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt al_index(PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt xy_index(PetscInt,PetscInt,PetscInt);
PetscInt Ind_1(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
PetscInt Ind_nx(PetscInt,PetscInt ,PetscInt ,PetscInt , PetscInt );


//Create Petsc Structures

PetscErrorCode initialize_petsc(struct Solver *,int, char **,struct AppCtx*);
PetscErrorCode initialize_fast_petsc(struct Solver *,int, char **,struct AppCtx*);
void Get_Nonzero_in_Rows(int *,struct AppCtx*);
PetscErrorCode initialize_jacobian(Mat,struct AppCtx*);

//Multigrid functions
PetscErrorCode Create_Restriction(Mat R,PetscInt nx, PetscInt ny);
PetscErrorCode Create_Interpolation(Mat R,PetscInt nx, PetscInt ny);
PetscErrorCode Initialize_PCMG(PC pc,Mat A,struct AppCtx*);

//Newton Solver
PetscErrorCode newton_solve(Vec,struct AppCtx*);
PetscErrorCode update_fast_vars(FILE*,struct SimState *,struct SimState *, struct AppCtx *,PetscReal);
//Calculate residuals and jacobians

//Nonlinear discretizations
// Derivative of CC with volume
PetscErrorCode calc_residual(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
//Derivative of CC with no volume
PetscErrorCode calc_residual_no_vol(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_no_vol(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
//Algebraic CC with volume
PetscErrorCode calc_residual_algebraic(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_algebraic(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
//Algebraic CC with no volume
PetscErrorCode calc_residual_algebraic_no_vol(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_algebraic_no_vol(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx

//Linear discretizations
// Linear system Algebraic CC with no volume
PetscErrorCode calc_residual_linear_algebraic(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_linear_algebraic(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
PetscErrorCode calc_residual_linear_deriv(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_linear_deriv(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx

//Multirate discretizations
//Derivative of CC with no volume
//Fast eqns
PetscErrorCode calc_residual_fast(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_fast(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
//Slow eqns
PetscErrorCode calc_residual_slow(SNES,Vec,Vec,void*); //void is masked AppCtx
PetscErrorCode calc_jacobian_slow(SNES, Vec, Mat,Mat, void*); //void is masked AppCtx
PetscErrorCode point_jacobian_fast(Mat ,int ,int, void *ctx);//point-wise solve
PetscErrorCode point_residual_fast(Vec ,int ,int, void *ctx);//point-wise solve
PetscErrorCode Point_Solve(struct SimState* ,struct SimState *,struct AppCtx*,PetscInt, PetscInt,PetscReal);

void point_ionmflux(PetscInt , PetscInt ,struct AppCtx*);
void gatevars_update_point(struct GateType *,struct SimState *,PetscReal ,struct AppCtx *,PetscInt , PetscInt );
void excitation_point(struct AppCtx* ,PetscReal ,PetscInt ,PetscInt );

//Find abs. max value of an array
PetscReal array_max(PetscReal *,size_t);
// abs. max, but for difference
PetscReal array_diff_max(PetscReal *,PetscReal *,size_t);
//Calc l2_norm of one array
PetscReal l2_norm(PetscReal *,size_t);

void print_all(struct AppCtx*);
const char* getfield(char* , int );
void find_print(int, int, double, int iter);
void compare_res(double *, int );
void write_data(FILE *,struct AppCtx *,PetscInt,int );
void write_point(FILE *,struct AppCtx *,PetscInt,int );
void write_fast(FILE *,struct AppCtx *,PetscInt,int );
void record_fast_timestep(FILE *,struct AppCtx *,PetscInt,int );
void measure_flux(FILE *,struct AppCtx *,PetscInt,int );
void init_events(struct AppCtx *);

void recombine(struct SimState*, struct AppCtx*);

#endif