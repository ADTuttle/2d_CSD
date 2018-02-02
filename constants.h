#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsnes.h>

//set global parameters here (extern constants)

extern const int use_en_deriv ; //if true, will use the derivative of the electroneutrality condition for the system of equations
#define separate_vol 1 //if true, will solve c,phi separate from alpha.
extern const int details; //if true, will show how many iterations were necessary for each newton solve, and the residual
extern const int two_points_exct;   //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
extern const int Profiling_on; //Turns timing of functions on/off.
extern const int Linear_Diffusion; //Changes to a linear discretization of electrodiffusion.
//extern int krecordfreq = 10; //determines how many time steps to run before recording the state variables
#define krecordfreq 1 //determines how many time steps to run before recording the state variables
extern const int save_one_var;
//basic ion extern constants
#define   Ni  3           //number of ion species (Na, K, Cl)
extern  const  PetscInt z[3];//valences of ion species
extern const   PetscReal D[3];      //diffusion coefficients in cm^2/sec

//grid parameters
// extern   PetscReal dt = 1e-2 ;        //time step (in s)
//extern 	PetscReal dt = 0.001;
//extern   PetscReal Time = 3e-2;
// extern   PetscReal Time = 1e-1;
//extern   PetscReal Time = 1;
#define Time 10.0
//extern  PetscReal   Time = 24;
//extern  PetscReal  Time = 60;//2e-2        //total simulated time in seconds
// extern    Time=2e-2
#define   Nc 3           //number of compartments
  // extern PetscInt  Nx;         //number of grid points in the x direction
  // extern PetscInt   Ny;      //number of grid points in the y direction
// extern  PetscInt  Nx = 100;
// extern  PetscInt  Ny = 100;
#define  Nx  64
#define  Ny  64
#define dx 0.01       //grid size in x direction (in cm)
#define dy 0.01        //grid size in y direction (in cm)
//extern PetscReal  dx = 0.02;        //grid size in x direction (in cm)
//extern PetscReal   dy = 0.02;        //grid size in y direction (in cm)
extern PetscReal  Lx;          //width of domain in cm (x)
extern PetscReal  Ly;         //length of domain in cm (y)
// extern PetscInt  Nt;     //total number of time steps

//Newton solve parameters
//number of variables to be solved for at each grid point
//extern PetscInt  Nv = (Ni+2)*Nc-1; //version if volume is included
//extern PetscInt  Nv = (Ni+1)*Nc; //version if volume is excluded
extern  PetscInt  Nv;   //combining the above two with the separate_vol
extern PetscInt  NA;     //total number of unknowns
extern PetscInt Nz;

 extern PetscInt  itermax;      //maximum Newton iterations allowed
//extern PetscInt itermax = 4;
extern PetscReal  reltol;    //relative tolerance


//physical extern constants
#define   Rc  8.314472e6    //gas constant in nJ/K/mmol
#define   T (273.15+37)    //absolute temperature in K
#define  FC  9.64853399e7 //Faraday constant in muC/mmol
extern const  PetscReal  RTFC;
#define pi 3.14159265358979323846 //M_PI is defined in some C implementations, but isn't in gen.
//So we define pi here.

//Bath variables
extern PetscReal cbath[3]; //Na, K, and Cl
extern PetscReal Batheps; //Bath diffusion multiplier
extern PetscReal phibath;

//excitation parameters
extern  PetscReal  pmax;          //max value for excitation
//extern  PetscReal  pmax  = 25;          //max value for excitation
extern  PetscReal  texct;          //time for excitation
extern  PetscInt  Nexct;          //number of grid points for excitation in each direction

//initial state setup
extern const PetscInt rest_state;        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
extern const PetscInt spatially_uniform ; //if true, all initial values and parameters are spatially uniform



//set "relaxed" volume fractions and initial volume fractions
#define  alphaon 0.5         //base neuronal volume fraction
#define alphaog 0.300000       //base glial volume fraction
      // alphaotemp = zeros(Nx,Ny,Nc-1)
      // alphaotemp[:,:,1] = alphaon
      // alphaotemp[:,:,2] = alphaog
// extern alphao = alphaotemp
      // alpha0 = copy(alphao) //initial volume fractions equal to "relaxed" volume fractions
extern const PetscReal alphao[2];
//extern PetscReal alphao[2]={0.5,0.3};
extern const PetscReal alpha0[2];
// membrane parameters
 #define   cmt  0.75e-3            //membrane capacitance in mF/cm^2
    // cmt = repmat([cmt],Nc-1)
 #define sa  1.586e-5          //membrane area in cm^2
 #define voli  2.16e-9         //intracellular volume in cm^3
extern const PetscReal vole;        //extracellular volume
extern const PetscReal ell ;    //average membrane separation in cm
extern const PetscReal cm[2];     //membrane capacitance in mF/cm^2 converted to mmol/cm^3

//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
extern PetscReal pNaT;                  //1e-4%0%1e-3%if set to 0, recovery possible
extern PetscReal pNaP;
extern PetscReal pKDR;
extern PetscReal pKA;

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
extern PetscReal pKLeak;     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
extern PetscReal pClLeak;   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
extern PetscReal pClLeakg;

//Glial KIR from Steinberg et al 2005
extern PetscReal pKIR;        //conductance in mS/cm^2 converted to mmol/cm^2/s
extern PetscReal pKLeakadjust;        //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
extern PetscReal mK;                 //pump current extern constant in mmol/cm^3=mol/l
extern PetscReal mNa;              //pump current extern constant in mmol/cm^3=mol/l
extern PetscReal glpump;              //multiplier to change glial pump rate
extern PetscReal npump;               //multiplier to change neuronal pump rate


// extern struct Simstate;
// extern struct FluxData;
// extern struct AppCtx;
// extern struct GateType;
// extern struct ConstVars;
// extern struct ExctType;
// extern struct Solver;

//struct SimState{
//	PetscReal c[Nx*Ny*Nc*Ni];
//	PetscReal phi[Nx*Ny*Nc];
//	PetscReal alpha[Nx*Ny*(Nc-1)];
//};
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
struct GatePoint{
	PetscReal mNaT;
	PetscReal hNaT;
	PetscReal gNaT;
	PetscReal mNaP;
	PetscReal hNaP;
	PetscReal gNaP;
	PetscReal mKDR;
	PetscReal gKDR;
	PetscReal mKA;
	PetscReal hKA;
	PetscReal gKA;
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
};

struct AppCtx{
    struct SimState *state_vars;
    struct SimState *state_vars_past;
    struct Solver *slvr;
    struct FluxData *flux;
    struct GateType *gate_vars;
    struct ExctType *gexct;
    struct ConstVars *con_vars;
    PetscReal *Dcs;
    PetscReal *Dcb;
    PetscReal dt;
};
PetscLogEvent event[10];
/*
struct ConstVars{
	PetscReal pNaKCl;
	PetscReal Imax[Nx*Ny];
	PetscReal pNaLeak[Nx*Ny];
	PetscReal Imaxg[Nx*Ny];
	PetscReal pNaLeakg[Nx*Ny];
	PetscReal ao[Nx*Ny];
	PetscReal zo[Nx*Ny];
	PetscReal kappa[Nx*Ny];
	PetscReal zeta1[Nx*Ny];
	PetscInt S; //boolean
	PetscReal zetaalpha[Nx*Ny];
};
*/
#endif