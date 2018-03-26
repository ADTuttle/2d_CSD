#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsnes.h>


//set global parameters here (constants)

#define use_en_deriv 1 //if true, will use the derivative of the electroneutrality condition for the system of equations
#define separate_vol 0 //if true, will solve c,phi separate from alpha.
#define details 0 //if true, will show how many iterations were necessary for each newton solve, and the residual
#define mid_points_exct 1
#define one_point_exct 0 //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
#define Profiling_on 0 //Turns timing of functions on/off.
#define Linear_Diffusion 0 //Changes to a linear discretization of electrodiffusion.
#define trecordstep 0.01//0.5 //determines how often to record
#define save_one_var 0 //Instead of saving all 14 vars, save just 1 (specified in write_data)

//basic ion extern constants
#define   Ni  3           //number of ion species (Na, K, Cl)
//extern const PetscInt z[3]; //valences of ion species
//extern const PetscReal D[3]; //diffusion coefficients in cm^2/sec
static const   PetscInt z[3] = {1,1,-1};//valences of ion species
static const   PetscReal D[3] = {1.33e-5, 1.96e-5, 2.03e-5};      //diffusion coefficients in cm^2/sec

//grid parameters
#define Time 5.0   //total simulated time in seconds
//#define  Time  60.0//2e-2
#define   Nc 3           //number of compartments
//#define Lx 0.32        //width of domain in cm (x)
//#define Ly 0.32         //length of domain in cm (y)
#define Lx 0.5        //width of domain in cm (x)
#define Ly 0.5         //length of domain in cm (y)
//#define  Nx  50     //number of grid points in the x direction
//#define  Ny  50     //number of grid points in the y direction
//#define dx (Lx/Nx)       //grid size in x direction (in cm)
//#define dy (Ly/Ny)        //grid size in y direction (in cm)

#define width_size  1

//number of variables to be solved for at each grid point
//#define  Nv  ((Ni+2)*Nc-1) //version if volume is included
//#define  Nv  ((Ni+1)*Nc) //version if volume is excluded
#define Nv  (((Ni+2)*Nc-1)*(!separate_vol)+((Ni+1)*Nc)*separate_vol)  //combining the above two with the separate_vol
//#define NA  (Nx*Ny*Nv)     //total number of unknowns
//#define Nz  (Ni*Nc*(4*(Nx-1)*Ny+4*(Ny-1)*Nx+2*Nx*Ny)+Ni*(Nc-1)*6*Nx*Ny+(Nc*Ni+1)*Nx*Ny+(Nc-1)*(6*Nx*Ny+Nx*Ny*(Nc-2)+Ni*2*Nx*Ny)) //number of nonzeros in Jacobian

//Newton solve parameters
#define  itermax  20      //maximum Newton iterations allowed
#define  reltol  1e-11    //relative tolerance


//physical extern constants
#define   Rc  8.314472e6    //gas constant in nJ/K/mmol
#define   T (273.15+37)    //absolute temperature in K
#define  FC  9.64853399e7 //Faraday constant in muC/mmol
#define RTFC (Rc*T/FC)
#define pi 3.14159265358979323846 //M_PI is defined in some C implementations, but isn't in gen.
//So we define pi here.

//Bath variables
//extern PetscReal cbath[3]; //Na, K, and Cl
static const PetscReal cbath[3]={140*1e-3,3.4*1e-3,120*1e-3}; //Na, K, and Cl
#define Batheps 1.0 //Bath diffusion multiplier
#define phibath (-0/RTFC) //Voltage of outside bath

//excitation parameters
#define pmax  (1e-1/3)          //max value for excitation
//#define pmax  (1e1)          //max value for excitation
//#define pmax  50          //max value for excitation
//#define texct 2         //time for excitation
//#define texct 0.5
#define texct 0.05         //time for excitation
#define Lexct 0.05          //Length of region for excitation in each direction
//#define Lexct 0.025          //Length of region for excitation in each direction

//initial state setup
#define rest_state  1        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
#define spatially_uniform  0 //if true, all initial values and parameters are spatially uniform


//set "relaxed" volume fractions and initial volume fractions
#define  alphaon 0.5         //base neuronal volume fraction
#define alphaog 0.300000       //base glial volume fraction
//extern const PetscReal alphao[2];
//extern const PetscReal alpha0[2];
static const PetscReal alphao[2]={alphaon,alphaog};
static const PetscReal alpha0[2]={alphaon,alphaog};
// membrane parameters
#define  cmt  0.75e-3            //membrane capacitance in mF/cm^2
 #define sa  1.586e-5          //membrane area in cm^2
 #define voli  2.16e-9         //intracellular volume in cm^3
#define vole (0.15*voli)
#define ell (voli+vole/sa)    //average membrane separation in cm
//extern const PetscReal cm[2];     //membrane capacitance in mF/cm^2 converted to mmol/cm^3
static const PetscReal cm[2] ={cmt*RTFC/FC/ell,cmt*RTFC/FC/ell};     //membrane capacitance in mF/cm^2 converted to mmol/cm^3


//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
#define pNaT  0                 //1e-4%0%1e-3%if set to 0, recovery possible
#define pNaP  2e-5
#define pKDR  1e-3
#define pKA  1e-4

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
#define pKLeak  (7e-2*RTFC/FC)     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
#define pClLeak  (10e-2*RTFC/FC)   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
#define pClLeakg  (5e-2*RTFC/FC)

//Glial KIR from Steinberg et al 2005
#define pKIR  (.13*RTFC/FC)        //conductance in mS/cm^2 converted to mmol/cm^2/s
#define pKLeakadjust  1.0       //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
#define mK  2e-3                 //pump current constant in mmol/cm^3=mol/l
#define mNa  7.7e-3              //pump current constant in mmol/cm^3=mol/l
#define glpump  1.0              //multiplier to change glial pump rate
#define npump  1.0


// Data Structures
struct SimState{
    PetscScalar *c; //Variable arrays
    PetscScalar *phi;
    PetscScalar *alpha;
    Vec v;  // Full variable vec
    Vec c_vec;  //Inidivual variable vecs
    Vec phi_vec;
    Vec al_vec;
    IS c_ind;    //Indices for c/phi/al.
    IS al_ind;
    IS phi_ind;

};

struct FluxData{
	PetscReal *mflux;
	PetscReal *dfdci;
	PetscReal *dfdce;
	PetscReal *dfdphim;
	PetscReal *wflow;
	PetscReal *dwdpi;
	PetscReal *dwdal;
};


struct GateType{
	PetscReal *mNaT;
	PetscReal *hNaT;
	PetscReal *gNaT;
	PetscReal *mNaP;
	PetscReal *hNaP;
	PetscReal *gNaP;
	PetscReal *mKDR;
	PetscReal *gKDR;
	PetscReal *mKA;
	PetscReal *hKA;
	PetscReal *gKA;
};

struct ExctType{
	PetscReal *pNa;
	PetscReal *pK;
	PetscReal *pCl;
};

struct ConstVars{
	PetscReal pNaKCl;
	PetscReal Imax;
	PetscReal pNaLeak;
	PetscReal Imaxg;
	PetscReal pNaLeakg;
	PetscReal *ao;
	PetscReal *zo;
	PetscReal kappa;
	PetscReal *zeta1;
	int S; //boolean
	PetscReal *zetaalpha;
};
struct Solver{
	Vec Q;      /* Update*/
	Vec Res;
  	Mat A;            /* linear system matrix */
    SNES snes;       /*Nonlinear solver context*/
  	KSP ksp;         /* linear solver context */
  	PC pc;           /* preconditioner context */
  	PetscMPIInt   size;

    PetscInt NA;
};


struct AppCtx{
    struct SimState *state_vars;
    struct SimState *state_vars_past;
    struct SimState *grid_vars;
    struct SimState *grid_vars_past;
    struct Solver *slvr;
    struct Solver *grid_slvr;
    struct FluxData *flux;
    struct GateType *gate_vars;
    struct GateType *gate_vars_past;
    struct GateType *grid_gate_vars;
    struct ExctType *gexct;
    struct ConstVars *con_vars;
    PetscReal *Dcs;
    PetscReal *Dcb;
    PetscReal dt;
    PetscReal dx;
    PetscReal dy;
    PetscInt Nx;
    PetscInt Ny;
    PetscInt Nz;
    PetscReal *dt_space;
};
PetscLogEvent event[10];

#endif