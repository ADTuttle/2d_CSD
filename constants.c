#include "constants.h"
#include "functions.h"

//set global parameters here (constants)

 const int use_en_deriv = 1; //if true, will use the derivative of the electroneutrality condition for the system of equations
    // const PetscInt separate_vol = 1; //if true, will solve c,phi separate from alpha.
    const int details = 0; //if true, will show how many iterations were necessary for each newton solve, and the residual
    const int two_points_exct = 0;   //if true, triggers SD at origin and (Nx/2,1) (halfway along x-axis)
    const int Profiling_on = 1; //Turns timing of functions on/off.
    const int Linear_Diffusion = 1; //Changes to a linear discretization of electrodiffusion.
//const int krecordfreq = 10; //determines how many time steps to run before recording the state variables
    // const int krecordfreq = 1; //determines how many time steps to run before recording the state variables
    const int save_one_var = 0;
//basic ion constants
           //number of ion species (Na, K, Cl)
    const   PetscInt z[3] = {1,1,-1};//valences of ion species
    const   PetscReal D[3] = {1.33e-5, 1.96e-5, 2.03e-5};      //diffusion coefficients in cm^2/sec

//grid parameters
// const   PetscReal dt = 1e-2 ;        //time step (in s)
    // const 	PetscReal dt = 0.01;
//const 	PetscReal dt = 0.001;
//const   PetscReal Time = 3e-2;
// const   PetscReal Time = 1e-1;
//const   PetscReal Time = 1;
    // const  PetscReal   Time = 10;
//const  PetscReal   Time = 24;
//const  PetscReal  Time = 60;//2e-2        //total simulated time in seconds
// const    Time=2e-2           //number of compartments
    // const PetscInt  Nx = 64;         //number of grid points in the x direction
    // const PetscInt   Ny = 64;      //number of grid points in the y direction
// const  PetscInt  Nx = 100;
// const  PetscInt  Ny = 100;
//const PetscInt   Nx = 32;
//const PetscInt  Ny = 32;
    // const PetscReal  dx = 0.01;        //grid size in x direction (in cm)
    // const PetscReal   dy = 0.01;        //grid size in y direction (in cm)
//const PetscReal  dx = 0.02;        //grid size in x direction (in cm)
//const PetscReal   dy = 0.02;        //grid size in y direction (in cm)
     PetscReal  Lx = Nx*dx;          //width of domain in cm (x)
     PetscReal  Ly = Ny*dy;         //length of domain in cm (y)
     // PetscInt  Nt = Time/dt;     //total number of time steps


//Newton solve parameters
//number of variables to be solved for at each grid point
//const PetscInt  Nv = (Ni+2)*Nc-1; //version if volume is included
//const PetscInt  Nv = (Ni+1)*Nc; //version if volume is excluded
     PetscInt  Nv = ((Ni+2)*Nc-1)*(!separate_vol)+((Ni+1)*Nc)*separate_vol;  //combining the above two with the separate_vol
     // PetscInt  NA = Nx*Ny*Nv;     //total number of unknowns
     PetscInt NA = Nx*Ny*(((Ni+2)*Nc-1)*(!separate_vol)+((Ni+1)*Nc)*separate_vol);
     PetscInt Nz = Ni*Nc*(4*(Nx-1)*Ny+4*(Ny-1)*Nx+2*Nx*Ny)+Ni*(Nc-1)*6*Nx*Ny+(Nc*Ni+1)*Nx*Ny+(Nc-1)*(6*Nx*Ny+Nx*Ny*(Nc-2)+Ni*2*Nx*Ny);

     PetscInt  itermax = 10;      //maximum Newton iterations allowed
//const PetscInt itermax = 4;
     PetscReal  reltol = 1e-11;    //relative tolerance


//physical constants
    // const  PetscReal  R = 8.314472e6;    //gas constant in nJ/K/mmol
    // const  PetscReal  T = 273.15+37;     //absolute temperature in K
    // const  PetscReal  FC = 9.64853399e7; //Faraday constant in muC/mmol
    const  PetscReal  RTFC = Rc*T/FC;
//So we define pi here.

//Bath variables
     PetscReal cbath[3]={140*1e-3,3.4*1e-3,120*1e-3}; //Na, K, and Cl
     PetscReal Batheps=1; //Bath diffusion multiplier
     PetscReal phibath=-0/RTFC;

//excitation parameters
      PetscReal  pmax  = 5;          //max value for excitation
//const  PetscReal  pmax  = 25;          //max value for excitation
      PetscReal  texct = 2;          //time for excitation
      PetscInt  Nexct = 5;          //number of grid points for excitation in each direction

//initial state setup
    const  PetscInt rest_state = 1;        //if true, membrane parameters are set so that the initial voltages and concentrations are at a rest state
    const PetscInt spatially_uniform = 0; //if true, all initial values and parameters are spatially uniform



//set "relaxed" volume fractions and initial volume fractions
    // const long double alphaon = .5;         //base neuronal volume fraction
    // const long double alphaog = .300000;        //base glial volume fraction
    // alphaotemp = zeros(Nx,Ny,Nc-1)
    // alphaotemp[:,:,1] = alphaon
    // alphaotemp[:,:,2] = alphaog
// const alphao = alphaotemp
    // alpha0 = copy(alphao) //initial volume fractions equal to "relaxed" volume fractions
    const PetscReal alphao[2]={alphaon,alphaog};
//const PetscReal alphao[2]={0.5,0.3};
    const PetscReal alpha0[2]={alphaon,alphaog};
// membrane parameters
    // const    PetscReal cmt = 0.75e-3;            //membrane capacitance in mF/cm^2
    // cmt = repmat([cmt],Nc-1)
    // const PetscReal sa = 1.586e-5;           //membrane area in cm^2
    // const PetscReal voli = 2.16e-9;          //intracellular volume in cm^3
    const PetscReal vole = 0.15*voli;        //extracellular volume
    const PetscReal ell = (voli+vole)/sa;    //average membrane separation in cm
    const PetscReal cm[2] ={cmt*RTFC/FC/ell,cmt*RTFC/FC/ell};     //membrane capacitance in mF/cm^2 converted to mmol/cm^3

//data for ion channel currents
//permeabilities in cm/s from Kager, 2000 and Yao, Huang, Miura, 2011.
     PetscReal pNaT = 0;                  //1e-4%0%1e-3%if set to 0, recovery possible
     PetscReal pNaP = 2e-5;
     PetscReal pKDR = 1e-3;
     PetscReal pKA = 1e-4;

//Leak conductances in mS/cm^2 from Kager, 2000 or Yao, Huang, Miura, 2011.
     PetscReal pKLeak = 7e-2*RTFC/FC;     //Kager:10e-2,Miura:7e-2%K Leak conductance in mS/cm^2 converted to mmol/cm^2/s
     PetscReal pClLeak = 10e-2*RTFC/FC;   //Kager:10e-2,Miura:20e-2%Cl Leak conductance in mS/cm^2 converted to mmol/cm^2/s
     PetscReal pClLeakg = 5e-2*RTFC/FC;

//Glial KIR from Steinberg et al 2005
     PetscReal pKIR = .13*RTFC/FC;        //conductance in mS/cm^2 converted to mmol/cm^2/s
     PetscReal pKLeakadjust = 1.0;        //parameter for varying glial permeability

//pump current, parameters from Yao, Huang, Miura, 2011
     PetscReal mK = 2e-3;                 //pump current constant in mmol/cm^3=mol/l
     PetscReal mNa = 7.7e-3;              //pump current constant in mmol/cm^3=mol/l
     PetscReal glpump = 1.0;              //multiplier to change glial pump rate
     PetscReal npump = 1.0;



