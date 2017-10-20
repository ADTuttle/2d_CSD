#ifndef __CONSTANTS__
#define __CONSTANTS__
#include <petsctime.h>
#include <petscksp.h>
#include <petscmat.h>

//set global parameters here (static constants)

static const int use_en_deriv = 0; //if true, will use the derivative of the electroneutrality condition for the system of equations
static const int use_direct_solve = 0; //if true, will use gmres instead of direct solve
static const int details = 1; //if true, will show how many iterations were necessary for each newton solve, and the residual
static const int krecordfreq = 10; //determines how many time steps to run before recording the state variables
static const int two_points_exct = 0;   //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
static const int savefreq = 500;

//basic ion static constants
static const   int Ni = 3;            //number of ion species (Na, K, Cl)
static const   int z[3] = {1,1,-1};//valences of ion species
static const   double D[3] = {1.33e-5, 1.96e-5, 2.03e-5};      //diffusion coefficients in cm^2/sec

//grid parameters
static const   double dt = 1e-2 ;        //time step (in s)
static const   double Time = 2e-2;
// static const    Time = 1e-1
// static const     Time=10
//static const    Time = 60//2e-2        //total simulated time in seconds
// static const    Time=2e-2
static const  int  Nc = 3;            //number of compartments
 // static const int  Nx = 100;         //number of grid points in the x direction
 // static const int   Ny = 100;      //number of grid points in the y direction
// static const  int  Nx = 50;
// static const  int  Ny = 50;
static const int   Nx = 5;
static const int  Ny = 5;
static const double  dx = 0.01;        //grid size in x direction (in cm)
static const double   dy = 0.01;        //grid size in y direction (in cm)
static const double  Lx = Nx*dx;          //width of domain in cm (x)
static const double  Ly = Ny*dy;         //length of domain in cm (y)
static const int  Nt = (int)Time/dt;     //total number of time steps

//Newton solve parameters
static const int  Nv = (Ni+2)*Nc-1;  //number of variables to be solved for at each grid point
static const int  NA = Nx*Ny*Nv;     //total number of unknowns
static const int Nz = Ni*Nc*(4*(Nx-1)*Ny+4*(Ny-1)*Nx+2*Nx*Ny)+Ni*(Nc-1)*6*Nx*Ny+(Nc*Ni+1)*Nx*Ny+(Nc-1)*(6*Nx*Ny+Nx*Ny*(Nc-2)+Ni*2*Nx*Ny);

static const int  itermax = 10;      //maximum Newton iterations allowed
static const double  reltol = 1e-11;    //relative tolerance


//physical static constants
static const  double  R = 8.314472e6;    //gas static constant in nJ/K/mmol
static const  double  T = 273.15+37;     //absolute temperature in K
static const  double  FC = 9.64853399e7; //Faraday static constant in muC/mmol
static const  double  RTFC = R*T/FC;
#define pi 3.14159265358979323846 //M_PI is defined in some C implementations, but isn't in gen.
//So we define pi here.

//Bath variables
static const double cbath[3]={140*1e-3,3.4*1e-3,120*1e-3}; //Na, K, and Cl
static const double Batheps=1; //Bath diffusion multiplier
static const double phibath=-0/RTFC;

//excitation parameters
static const  double  pmax  = 5;          //max value for excitation
static const  double  texct = 2;          //time for excitation
static const  int  Nexct = 5;          //number of grid points for excitation in each directionexport pmax, texct, Nexct

//initial state setup
static const int rest_state = 1;        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
static const int spatially_uniform = 0; //if true, all initial values and parameters are spatially uniform



//set "relaxed" volume fractions and initial volume fractions
static const double alphaon = .5;         //base neuronal volume fraction
static const double alphaog = .3;        //base glial volume fraction
      // alphaotemp = zeros(Nx,Ny,Nc-1)
      // alphaotemp[:,:,1] = alphaon
      // alphaotemp[:,:,2] = alphaog
// static const alphao = alphaotemp
      // alpha0 = copy(alphao) //initial volume fractions equal to "relaxed" volume fractions
static const double alphao[2]={alphaon,alphaog};
static const double alpha0[2]={alphaon,alphaog};
// membrane parameters
static const    double cmt = 0.75e-3;            //membrane capacitance in mF/cm^2
      // cmt = repmat([cmt],Nc-1)
static const double sa = 1.586e-5;           //membrane area in cm^2
static const double voli = 2.16e-9;          //intracellular volume in cm^3
static const double vole = 0.15*voli;        //extracellular volume
static const double ell = (voli+vole)/sa;    //average membrane separation in cm
static const double cm[2] ={cmt*RTFC/FC/ell,cmt*RTFC/FC/ell};     //membrane capacitance in mF/cm^2 converted to mmol/cm^3

//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
static const double pNaT = 0;                  //1e-4%0%1e-3%if set to 0, recovery possible
static const double pNaP = 2e-5;
static const double pKDR = 1e-3;
static const double pKA = 1e-4;

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
static const double pKLeak = 7e-2*RTFC/FC;     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
static const double pClLeak = 10e-2*RTFC/FC;   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
static const double pClLeakg = 5e-2*RTFC/FC;

//Glial KIR from Steinberg et al 2005
static const double pKIR = .13*RTFC/FC;        //conductance in mS/cm^2 converted to mmol/cm^2/s
static const double pKLeakadjust = 1.0;        //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
static const double mK = 2e-3;                 //pump current static constant in mmol/cm^3=mol/l
static const double mNa = 7.7e-3;              //pump current static constant in mmol/cm^3=mol/l
static const double glpump = 1.0;              //multiplier to change glial pump rate
static const double npump = 1.0;               //multiplier to change neuronal pump rate



struct SimState{
	double c[Nx*Ny*Nc*Ni];
	double phi[Nx*Ny*Nc];
	double alpha[Nx*Ny*(Nc-1)];
};

struct FluxData{
	double mflux[Nx*Ny*Ni*Nc];
	double dfdci[Nx*Ny*Ni*Nc];
	double dfdce[Nx*Ny*Ni*Nc];
	double dfdphim[Nx*Ny*Ni*Nc];
	double wflow[Nx*Ny*(Nc-1)];
	double dwdpi[Nx*Ny*(Nc-1)];
	double dwdal[Nx*Ny*(Nc-1)];
};


struct GateType{
	double mNaT[Nx*Ny];
	double hNaT[Nx*Ny];
	double gNaT[Nx*Ny];
	double mNaP[Nx*Ny];
	double hNaP[Nx*Ny];
	double gNaP[Nx*Ny];
	double mKDR[Nx*Ny];
	double gKDR[Nx*Ny];
	double mKA[Nx*Ny];
	double hKA[Nx*Ny];
	double gKA[Nx*Ny];
};
struct GatePoint{
	double mNaT;
	double hNaT;
	double gNaT;
	double mNaP;
	double hNaP;
	double gNaP;
	double mKDR;
	double gKDR;
	double mKA;
	double hKA;
	double gKA;
};

struct ExctType{
	double pNa[Nx*Ny];
	double pK[Nx*Ny];
	double pCl[Nx*Ny];
};

struct ConstVars{
	double pNaKCl;
	double Imax;
	double pNaLeak;
	double Imaxg;
	double pNaLeakg;
	double ao[Nc];
	double zo[Nc];
	double kappa;
	double zeta1[Nc-1];
	int S; //boolean
	double zetaalpha[Nc-1];
};
struct Solver{
	PetscInt row[Nz];
	PetscInt col[Nz];
	Vec Q;      /* Update*/
	Vec Res;
  	Mat A;            /* linear system matrix */
  	KSP ksp;         /* linear solver context */
  	PC pc;           /* preconditioner context */
  	PetscMPIInt    size;
};

/*
struct ConstVars{
	double pNaKCl;
	double Imax[Nx*Ny];
	double pNaLeak[Nx*Ny];
	double Imaxg[Nx*Ny];
	double pNaLeakg[Nx*Ny];
	double ao[Nx*Ny];
	double zo[Nx*Ny];
	double kappa[Nx*Ny];
	double zeta1[Nx*Ny];
	int S; //boolean
	double zetaalpha[Nx*Ny];
};
*/





#endif