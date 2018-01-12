#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>

//set global parameters here (static constants)

static const int use_en_deriv = 0; //if true, will use the derivative of the electroneutrality condition for the system of equations
static const int use_direct_solve = 0; //if true, will use gmres instead of direct solve
static const int details = 0; //if true, will show how many iterations were necessary for each newton solve, and the residual
static const int krecordfreq = 10; //determines how many time steps to run before recording the state variables
static const int two_points_exct = 0;   //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
static const int savefreq = 500;

//basic ion static constants
static const   PetscInt Ni = 3;            //number of ion species (Na, K, Cl)
static const   PetscInt z[3] = {1,1,-1};//valences of ion species
static const   PetscReal D[3] = {1.33e-5, 1.96e-5, 2.03e-5};      //diffusion coefficients in cm^2/sec

//grid parameters
// static const   PetscReal dt = 1e-2 ;        //time step (in s)
static const 	PetscReal dt = 0.01;
//static const   PetscReal Time = 3e-2;
// static const   PetscReal Time = 1e-1;
//static const   PetscReal Time = 1;
 static const  PetscReal   Time=10;
//static const    Time = 60//2e-2        //total simulated time in seconds
// static const    Time=2e-2
static const  PetscInt  Nc = 3;            //number of compartments
  static const PetscInt  Nx = 64;         //number of grid points in the x direction
  static const PetscInt   Ny = 64;      //number of grid points in the y direction
// static const  PetscInt  Nx = 100;
// static const  PetscInt  Ny = 100;
//static const PetscInt   Nx = 9;
//static const PetscInt  Ny = 8;
static const PetscReal  dx = 0.01;        //grid size in x direction (in cm)
static const PetscReal   dy = 0.01;        //grid size in y direction (in cm)
static const PetscReal  Lx = Nx*dx;          //width of domain in cm (x)
static const PetscReal  Ly = Ny*dy;         //length of domain in cm (y)
static const PetscInt  Nt = Time/dt;     //total number of time steps

static const PetscReal numrecords = Time/(dt*krecordfreq);

//Newton solve parameters
static const PetscInt  Nv = (Ni+2)*Nc-1;  //number of variables to be solved for at each grid point
static const PetscInt  NA = Nx*Ny*Nv;     //total number of unknowns
static const PetscInt Nz = Ni*Nc*(4*(Nx-1)*Ny+4*(Ny-1)*Nx+2*Nx*Ny)+Ni*(Nc-1)*6*Nx*Ny+(Nc*Ni+1)*Nx*Ny+(Nc-1)*(6*Nx*Ny+Nx*Ny*(Nc-2)+Ni*2*Nx*Ny);

 static const PetscInt  itermax = 10;      //maximum Newton iterations allowed
//static const PetscInt itermax = 4;
static const PetscReal  reltol = 1e-11;    //relative tolerance


//physical static constants
static const  PetscReal  R = 8.314472e6;    //gas static constant in nJ/K/mmol
static const  PetscReal  T = 273.15+37;     //absolute temperature in K
static const  PetscReal  FC = 9.64853399e7; //Faraday static constant in muC/mmol
static const  PetscReal  RTFC = R*T/FC;
#define pi 3.14159265358979323846 //M_PI is defined in some C implementations, but isn't in gen.
//So we define pi here.

//Bath variables
static const PetscReal cbath[3]={140*1e-3,3.4*1e-3,120*1e-3}; //Na, K, and Cl
static const PetscReal Batheps=1; //Bath diffusion multiplier
static const PetscReal phibath=-0/RTFC;

//excitation parameters
static const  PetscReal  pmax  = 5;          //max value for excitation
static const  PetscReal  texct = 2;          //time for excitation
static const  PetscInt  Nexct = 5;          //number of grid points for excitation in each directionexport pmax, texct, Nexct

//initial state setup
static const PetscInt rest_state = 1;        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
static const PetscInt spatially_uniform = 0; //if true, all initial values and parameters are spatially uniform



//set "relaxed" volume fractions and initial volume fractions
static const long double alphaon = .5;         //base neuronal volume fraction
static const long double alphaog = .300000;        //base glial volume fraction
      // alphaotemp = zeros(Nx,Ny,Nc-1)
      // alphaotemp[:,:,1] = alphaon
      // alphaotemp[:,:,2] = alphaog
// static const alphao = alphaotemp
      // alpha0 = copy(alphao) //initial volume fractions equal to "relaxed" volume fractions
static const PetscReal alphao[2]={alphaon,alphaog};
//static const PetscReal alphao[2]={0.5,0.3};
static const PetscReal alpha0[2]={alphaon,alphaog};
// membrane parameters
static const    PetscReal cmt = 0.75e-3;            //membrane capacitance in mF/cm^2
      // cmt = repmat([cmt],Nc-1)
static const PetscReal sa = 1.586e-5;           //membrane area in cm^2
static const PetscReal voli = 2.16e-9;          //intracellular volume in cm^3
static const PetscReal vole = 0.15*voli;        //extracellular volume
static const PetscReal ell = (voli+vole)/sa;    //average membrane separation in cm
static const PetscReal cm[2] ={cmt*RTFC/FC/ell,cmt*RTFC/FC/ell};     //membrane capacitance in mF/cm^2 converted to mmol/cm^3

//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
static const PetscReal pNaT = 0;                  //1e-4%0%1e-3%if set to 0, recovery possible
static const PetscReal pNaP = 2e-5;
static const PetscReal pKDR = 1e-3;
static const PetscReal pKA = 1e-4;

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
static const PetscReal pKLeak = 7e-2*RTFC/FC;     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
static const PetscReal pClLeak = 10e-2*RTFC/FC;   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
static const PetscReal pClLeakg = 5e-2*RTFC/FC;

//Glial KIR from Steinberg et al 2005
static const PetscReal pKIR = .13*RTFC/FC;        //conductance in mS/cm^2 converted to mmol/cm^2/s
static const PetscReal pKLeakadjust = 1.0;        //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
static const PetscReal mK = 2e-3;                 //pump current static constant in mmol/cm^3=mol/l
static const PetscReal mNa = 7.7e-3;              //pump current static constant in mmol/cm^3=mol/l
static const PetscReal glpump = 1.0;              //multiplier to change glial pump rate
static const PetscReal npump = 1.0;               //multiplier to change neuronal pump rate



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
  	KSP ksp;         /* linear solver context */
  	PC pc;           /* preconditioner context */
  	PetscMPIInt   size;
};

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