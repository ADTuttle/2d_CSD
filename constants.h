#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsnes.h>

//set global parameters here (extern constants)

extern const int use_en_deriv ; //if true, will use the derivative of the electroneutrality condition for the system of equations
extern const PetscInt separate_vol; //if true, will solve c,phi separate from alpha.
extern const int details; //if true, will show how many iterations were necessary for each newton solve, and the residual
extern const int two_points_exct;   //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
extern const int Profiling_on; //Turns timing of functions on/off.
extern const int Linear_Diffusion; //Changes to a linear discretization of electrodiffusion.
//extern const int krecordfreq = 10; //determines how many time steps to run before recording the state variables
extern const int krecordfreq; //determines how many time steps to run before recording the state variables
extern const int save_one_var;
//basic ion extern constants
static const   PetscInt Ni = 3;            //number of ion species (Na, K, Cl)
extern const   PetscInt z[3];//valences of ion species
extern const   PetscReal D[3];      //diffusion coefficients in cm^2/sec

//grid parameters
// extern const   PetscReal dt = 1e-2 ;        //time step (in s)
extern const 	PetscReal dt;
//extern const 	PetscReal dt = 0.001;
//extern const   PetscReal Time = 3e-2;
// extern const   PetscReal Time = 1e-1;
//extern const   PetscReal Time = 1;
 extern const  PetscReal   Time;
//extern const  PetscReal   Time = 24;
//extern const  PetscReal  Time = 60;//2e-2        //total simulated time in seconds
// extern const    Time=2e-2
static const  PetscInt  Nc = 3;            //number of compartments
  // extern const PetscInt  Nx;         //number of grid points in the x direction
  // extern const PetscInt   Ny;      //number of grid points in the y direction
// extern const  PetscInt  Nx = 100;
// extern const  PetscInt  Ny = 100;
static const PetscInt   Nx = 32;
static const PetscInt  Ny = 32;
extern const PetscReal  dx;        //grid size in x direction (in cm)
extern const PetscReal   dy;        //grid size in y direction (in cm)
//extern const PetscReal  dx = 0.02;        //grid size in x direction (in cm)
//extern const PetscReal   dy = 0.02;        //grid size in y direction (in cm)
extern const PetscReal  Lx;          //width of domain in cm (x)
extern const PetscReal  Ly;         //length of domain in cm (y)
extern const PetscInt  Nt;     //total number of time steps

extern const PetscReal numrecords;

//Newton solve parameters
//number of variables to be solved for at each grid point
//extern const PetscInt  Nv = (Ni+2)*Nc-1; //version if volume is included
//extern const PetscInt  Nv = (Ni+1)*Nc; //version if volume is excluded
extern const PetscInt  Nv;   //combining the above two with the separate_vol
extern const PetscInt  NA;     //total number of unknowns
extern const PetscInt Nz;

 extern const PetscInt  itermax;      //maximum Newton iterations allowed
//extern const PetscInt itermax = 4;
extern const PetscReal  reltol;    //relative tolerance


//physical extern constants
extern const  PetscReal  R;    //gas extern constant in nJ/K/mmol
extern const  PetscReal  T;     //absolute temperature in K
extern const  PetscReal  FC; //Faraday extern constant in muC/mmol
extern const  PetscReal  RTFC;
#define pi 3.14159265358979323846 //M_PI is defined in some C implementations, but isn't in gen.
//So we define pi here.

//Bath variables
extern const PetscReal cbath[3]; //Na, K, and Cl
extern const PetscReal Batheps; //Bath diffusion multiplier
extern const PetscReal phibath;

//excitation parameters
extern const  PetscReal  pmax;          //max value for excitation
//extern const  PetscReal  pmax  = 25;          //max value for excitation
extern const  PetscReal  texct;          //time for excitation
extern const  PetscInt  Nexct;          //number of grid points for excitation in each direction

//initial state setup
extern const PetscInt rest_state;        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
extern const PetscInt spatially_uniform ; //if true, all initial values and parameters are spatially uniform



//set "relaxed" volume fractions and initial volume fractions
extern const long double alphaon;         //base neuronal volume fraction
extern const long double alphaog;        //base glial volume fraction
      // alphaotemp = zeros(Nx,Ny,Nc-1)
      // alphaotemp[:,:,1] = alphaon
      // alphaotemp[:,:,2] = alphaog
// extern const alphao = alphaotemp
      // alpha0 = copy(alphao) //initial volume fractions equal to "relaxed" volume fractions
extern const PetscReal alphao[2];
//extern const PetscReal alphao[2]={0.5,0.3};
extern const PetscReal alpha0[2];
// membrane parameters
extern const    PetscReal cmt;            //membrane capacitance in mF/cm^2
      // cmt = repmat([cmt],Nc-1)
extern const PetscReal sa;           //membrane area in cm^2
extern const PetscReal voli;          //intracellular volume in cm^3
extern const PetscReal vole;        //extracellular volume
extern const PetscReal ell ;    //average membrane separation in cm
extern const PetscReal cm[2];     //membrane capacitance in mF/cm^2 converted to mmol/cm^3

//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
extern const PetscReal pNaT;                  //1e-4%0%1e-3%if set to 0, recovery possible
extern const PetscReal pNaP;
extern const PetscReal pKDR;
extern const PetscReal pKA;

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
extern const PetscReal pKLeak;     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
extern const PetscReal pClLeak;   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
extern const PetscReal pClLeakg;

//Glial KIR from Steinberg et al 2005
extern const PetscReal pKIR;        //conductance in mS/cm^2 converted to mmol/cm^2/s
extern const PetscReal pKLeakadjust;        //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
extern const PetscReal mK;                 //pump current extern constant in mmol/cm^3=mol/l
extern const PetscReal mNa;              //pump current extern constant in mmol/cm^3=mol/l
extern const PetscReal glpump;              //multiplier to change glial pump rate
extern const PetscReal npump;               //multiplier to change neuronal pump rate


// extern const struct Simstate;
// extern const struct FluxData;
// extern const struct AppCtx;
// extern const struct GateType;
// extern const struct ConstVars;
// extern const struct ExctType;
// extern const struct Solver;

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
	PetscReal mflux[Nx*Ny*Ni*Nc];
	PetscReal dfdci[Nx*Ny*Ni*Nc];
	PetscReal dfdce[Nx*Ny*Ni*Nc];
	PetscReal dfdphim[Nx*Ny*Ni*Nc];
	PetscReal wflow[Nx*Ny*(Nc-1)];
	PetscReal dwdpi[Nx*Ny*(Nc-1)];
	PetscReal dwdal[Nx*Ny*(Nc-1)];
};


struct GateType{
	PetscReal mNaT[Nx*Ny];
	PetscReal hNaT[Nx*Ny];
	PetscReal gNaT[Nx*Ny];
	PetscReal mNaP[Nx*Ny];
	PetscReal hNaP[Nx*Ny];
	PetscReal gNaP[Nx*Ny];
	PetscReal mKDR[Nx*Ny];
	PetscReal gKDR[Nx*Ny];
	PetscReal mKA[Nx*Ny];
	PetscReal hKA[Nx*Ny];
	PetscReal gKA[Nx*Ny];
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
	PetscReal pNa[Nx*Ny];
	PetscReal pK[Nx*Ny];
	PetscReal pCl[Nx*Ny];
};

struct ConstVars{
	PetscReal pNaKCl;
	PetscReal Imax;
	PetscReal pNaLeak;
	PetscReal Imaxg;
	PetscReal pNaLeakg;
	PetscReal ao[Nc];
	PetscReal zo[Nc];
	PetscReal kappa;
	PetscReal zeta1[Nc-1];
	int S; //boolean
	PetscReal zetaalpha[Nc-1];
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
    PetscReal Dcs[Nx*Ny*Ni*Nc*2];
    PetscReal Dcb[Nx*Ny*Ni*Nc*2];
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