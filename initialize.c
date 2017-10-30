#include "constants.h"
#include "functions.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

PetscInt c_index(PetscInt x,PetscInt y,PetscInt comp,PetscInt ion)
{
	return Nc*Ni* (Nx * y + x) + comp*Ni+ion; 
}
PetscInt phi_index(PetscInt x,PetscInt y,PetscInt comp)
{
	return Nc* (Nx * y + x) + comp; 
}
PetscInt al_index(PetscInt x,PetscInt y,PetscInt comp)
{
	return (Nc-1)* (Nx * y + x) + comp; 
}
PetscInt xy_index(PetscInt x,PetscInt y)
{
	return Nx*y+x;
}
PetscInt Ind_1(PetscInt x,PetscInt y,PetscInt ion,PetscInt comp)
{
  return Nv*(Nx*y+x)+ion*Nc+comp;
}


void init(struct SimState *state_vars)
{
	for(PetscInt x=0;x<Nx;x++)
	{
		for(PetscInt y=0;y<Ny;y++)
		{
			//initial volume fractions
			state_vars->alpha[al_index(x,y,0)]=alphao[0];
			state_vars->alpha[al_index(x,y,1)]=alphao[1]; 
			//initial voltages (dimensionless)
			state_vars->phi[phi_index(x,y,0)] = -70/RTFC; //neuronal voltage
			state_vars->phi[phi_index(x,y,1)] = -85/RTFC; //glial voltage
			state_vars->phi[phi_index(x,y,2)] = -0/RTFC; //extracell voltage
			//initial concentrations in mmol/cm^3=1e-3 mmol/l
			state_vars->c[c_index(x,y,0,0)] = 10e-3;     //neuronal Na concentration
			state_vars->c[c_index(x,y,1,0)] = 10e-3;      //glial Na concentration
			state_vars->c[c_index(x,y,2,0)] = 140e-3;     //extracellular Na concentration
			state_vars->c[c_index(x,y,0,1)] = 130e-3;     //neuronal K concentration
			state_vars->c[c_index(x,y,1,1)] = 130e-3;     //glial K concentration
			state_vars->c[c_index(x,y,2,1)] = 3.4e-3;     //extracellular K concentration
			state_vars->c[c_index(x,y,0,2)] = 10e-3;       //neuronal Cl concentration
			state_vars->c[c_index(x,y,1,2)] = 10e-3; 		//glial Cl concentraion
			state_vars->c[c_index(x,y,2,2)] = 120e-3;       //143.5e-3%extracellular Cl

		}
	}
}

void set_params(struct SimState* state_vars,struct ConstVars* con_vars,struct GateType* gate_vars,struct FluxData *flux)
{
	//Everything that follows will asume spatially uniform
	//At rest state
	double c[Ni*Nc];
	double phi[Nc];
	double alpha[Nc];
	for(PetscInt comp=0;comp<Nc;comp++)
	{
		for(PetscInt ion=0;ion<Ni;ion++)
		{
			c[c_index(0,0,comp,ion)]=state_vars->c[c_index(0,0,comp,ion)];
		}
		phi[phi_index(0,0,comp)]=state_vars->phi[phi_index(0,0,comp)];
		if(comp<Nc-1)
		{
			alpha[al_index(0,0,comp)]=state_vars->alpha[al_index(0,0,comp)];
		}
		else
		{
			alpha[al_index(0,0,comp)]=1; //Adding in the extracellular vol
			for(PetscInt i=0;i<Nc-1;i++)
			{
				alpha[al_index(0,0,comp)]-=alpha[al_index(0,0,i)];
			}
		}
	}
	double vm = phi[phi_index(0,0,0)]-phi[phi_index(0,0,2)]; //neuronal membrane potential
    double vmg = phi[phi_index(0,0,1)]-phi[phi_index(0,0,2)]; //glial membrane potential

    //compute neuronal Cl concentration (since neuron has only leak conductance, must be at reversal potential for Cl)

    c[c_index(0,0,0,2)] = c[c_index(0,0,2,2)]*exp(vm);
    //set glial Cl concentration equal to neuronal Cl concentration
    c[c_index(0,0,1,2)] = c[c_index(0,0,0,2)];

    // struct FluxPoPetscInt *flux;
    // flux = (struct FluxPoint*)malloc(sizeof(struct FluxPoint));

    //compute cotransporter permeability so that glial Cl is at rest
    mclin(flux,c_index(0,0,1,2),pClLeakg,-1,c[c_index(0,0,1,2)],c[c_index(0,0,2,2)],vmg,0);
    con_vars->pNaKCl = -flux->mflux[c_index(0,0,1,2)]/2/log(c[c_index(0,0,1,0)]*c[c_index(0,0,1,1)]*c[c_index(0,0,1,2)]*c[c_index(0,0,1,2)]/(c[c_index(0,0,2,0)]*c[c_index(0,0,2,1)]*c[c_index(0,0,2,2)]*c[c_index(0,0,2,2)]));
    double NaKCl = -flux->mflux[c_index(0,0,1,2)]/2;

    //compute gating variables
    gatevars_update(gate_vars,state_vars,0,1);

    //compute K channel currents (neuron)
    double pKGHK = pKDR*gate_vars->gKDR[0]+pKA*gate_vars->gKA[0];
    //Initialize the KGHK flux
    mcGoldman(flux,c_index(0,0,0,1),pKGHK,1,c[c_index(0,0,0,1)],c[c_index(0,0,Nc-1,1)],vm,0);
    //Add the KLeak flux to it
    mclin(flux,c_index(0,0,0,1),pKLeak,1,c[c_index(0,0,0,1)],c[c_index(0,0,Nc-1,1)],vm,1);

    //compute neuronal ATPase value
    con_vars->Imax = flux->mflux[c_index(0,0,0,1)]*(pow(1+mK/c[c_index(0,0,Nc-1,1)],2)*pow(1+mNa/c[c_index(0,0,0,0)],3))/2;


    //compute neuronal sodium currents and leak permeability value
    double pNaGHK = pNaT*gate_vars->gNaT[0]+pNaP*gate_vars->gNaP[0];
    mcGoldman(flux,c_index(0,0,0,0),pNaGHK,1,c[c_index(0,0,0,0)],c[c_index(0,0,Nc-1,0)],vm,0);
    double Ipump = npump*con_vars->Imax/(pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,0,0)]),3));
    con_vars->pNaLeak = (-flux->mflux[c_index(0,0,0,0)]-3*Ipump)/(log(c[c_index(0,0,0,0)]/c[c_index(0,0,Nc-1,0)])+vm);

    //compute K channel currents (glial)
    double pKLinG = pKIR*inwardrect(c[c_index(0,0,1,1)],c[c_index(0,0,Nc-1,1)],vmg)*pKLeakadjust;
    mclin(flux,c_index(0,0,1,1),pKLinG,1,c[c_index(0,0,1,1)],c[c_index(0,0,Nc-1,1)],vmg,0);
    flux->mflux[c_index(0,0,1,1)] += NaKCl;

    //compute glial ATPase value
    con_vars->Imaxg = flux->mflux[c_index(0,0,1,1)]*pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,1,0)]),3)/2;

    //compute glial sodium current and leak permeability value
    double Ipumpg = glpump*con_vars->Imaxg/(pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,1,0)]),3));
    con_vars->pNaLeakg = (-NaKCl-3*Ipumpg)/(log(c[c_index(0,0,1,0)]/c[c_index(0,0,Nc-1,0)])+vmg);

    //Compute resting organic anion amounts and average valences
    //set extracellular organic anion amounts and valence to ensure electroneutrality
    con_vars->ao[Nc-1] = 5e-4;
    double cmphi[Nc];
    double osmotic;
    for(PetscInt k=0;k<Nc-1;k++)
    {
      cmphi[k] = cm[k]*(phi[k]-phi[Nc-1]);
      cmphi[Nc-1] += cmphi[k];
      //set intracellular organic anion amounts to ensure osmotic pressure balance
      osmotic=0;
      for(PetscInt ion=0;ion<Ni;ion++)
      {
      	osmotic += c[c_index(0,0,Nc-1,ion)]-c[c_index(0,0,k,ion)];
      }
      con_vars->ao[k] = alpha[k]*(con_vars->ao[Nc-1]/alpha[Nc-1]+osmotic);
      //set average valence to ensure electroneutrality
      con_vars->zo[k] = (-cz(c,z,0,0,k)*alpha[k]+cmphi[k])/con_vars->ao[k];
    }
    con_vars->zo[Nc-1] = (-cz(c,z,0,0,Nc-1)*alpha[Nc-1]-cmphi[Nc-1])/con_vars->ao[Nc-1];
    //Copy the point data to vectors.
    //Only needed for uniform data
    for(PetscInt x=0;x<Nx;x++)
    {
    	for(PetscInt y=0;y<Ny;y++)
    	{
    		//Gating variables (already set in gatevars_update)
    		//We changed c_index(0,0,0/1,2), neuronal/glial Cl.
    		state_vars->c[c_index(x,y,0,2)] = c[c_index(0,0,0,2)];
    		state_vars->c[c_index(x,y,1,2)] = c[c_index(0,0,1,2)];
    	}
    }

    //Set kappa to 0 for no flow
    con_vars->kappa = 0;

 	//parameters for osmotic water flow
  	
	 double zetaadjust = 1; //modify glial permeability  
  	for(PetscInt comp=0;comp<Nc-1;comp++)
  	{
  		//based on B.E. Shapiro dissertation (2000) 
  		con_vars->zeta1[comp] = 5.4e-5;  //hydraulic permeability in cm/sec/(mmol/cm^3)
	  	con_vars->zeta1[comp] /= ell;  //conversion to 1/sec/(mmol/cm^3)
	    //based on Strieter, Stephenson, Palmer,
	  	//Weinstein, Journal or General Physiology, 1990.
	  	//zeta=7e-8%6e-10%hydraulic permeability in cm/sec/mmHg
	  	//zeta=zeta*7.501e-6%conversion to cm/sec/mPa
	  	//zeta=zeta*R*T%conversion to cm/sec/(mmol/cm^3)
	  	//zeta=zeta/ell%conversion to 1/sec/(mmol/cm^3)
	  	if(comp==1)
	  	{          //parameter for varying glial hydraulic permeability
	  		con_vars->zeta1[comp] *= zetaadjust; //adjust glial hydraulic permeability
	  	}
	  	con_vars->zetaalpha[comp] = 0;  //stiffness constant or 1/stiffness constant
  	}

    con_vars->S = 1;  //Indicates whether zetaalpha is the stiffness (true) or 1/stiffness (false)



    free(flux);
    return;
}

void initialize_data(struct SimState *state_vars,struct GateType* gate_vars,struct ConstVars* con_vars,struct Solver *slvr,struct FluxData *flux)
{
	double reltol = 1e-11;
	double tol = reltol*array_max(state_vars->c,(size_t)Nx*Ny*Nc*Ni);
  	double rsd = 1.0;
  	double *cp;
  	cp = (double *)malloc(sizeof(double)*Nx*Ny*Ni*Nc);
  	
    //Compute Gating variables
    //compute gating variables
    gatevars_update(gate_vars,state_vars,0,1);

  	//Initialize and comput the excitation (it's zeros here)
  	struct ExctType *gexct;
  	gexct = (struct ExctType*)malloc(sizeof(struct ExctType));
  	excitation(gexct,texct+1);
  	PetscInt k = 0;
  	double dt_temp = 0.1;
    // double dt_temp = 0.01;
  	
  	while(rsd>tol && dt_temp*k<10)
  	{
    	memcpy(cp,state_vars->c,sizeof(double)*Nx*Ny*Ni*Nc);
    	newton_solve(state_vars, dt_temp, gate_vars, gexct,con_vars,slvr,flux);
    	gatevars_update(gate_vars,state_vars,dt_temp*1e3,0);
    	rsd = array_diff_max(state_vars->c,cp,(size_t)Nx*Ny*Nc*Ni)/dt_temp;
        printf("Init_Data rsd: %f, Tol: %f\n",1e6*rsd,1e6*tol);
    	k++;
	}
  	
  	free(cp);
  	free(gexct);
	if(rsd>1e-7)
  	{
    	fprintf(stderr, "Did not converge! Aborting...\n");
    	exit(EXIT_FAILURE); /* indicate failure.*/
	}	
  	else
  	{
    	return;
	}
	
}


PetscErrorCode initialize_petsc(struct Solver *slvr,int argc, char **argv)
{
    PetscErrorCode ierr;
	//Init Petsc
	PetscInitialize(&argc,&argv,(char*)0,NULL);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&slvr->size);CHKERRQ(ierr);

  	//Create Vectors
    /*
  	ierr = VecCreate(PETSC_COMM_WORLD,&slvr->Q);CHKERRQ(ierr);
  	ierr = PetscObjectSetName((PetscObject) slvr->Q, "Solution");CHKERRQ(ierr);
  	ierr = VecSetSizes(slvr->Q,PETSC_DECIDE,NA);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(slvr->Q);CHKERRQ(ierr);
  	ierr = VecDuplicate(slvr->Q,&slvr->Res);CHKERRQ(ierr);
    */
    ierr = VecCreateSeq(PETSC_COMM_WORLD,NA,&slvr->Q);CHKERRQ(ierr);
    ierr = VecSetFromOptions(slvr->Q);CHKERRQ(ierr);
    ierr = VecDuplicate(slvr->Q,&slvr->Res);CHKERRQ(ierr);

  	//Create Matrix
  	ierr = MatCreate(PETSC_COMM_WORLD,&slvr->A);CHKERRQ(ierr);
  	ierr = MatSetType(slvr->A,MATSEQAIJ);CHKERRQ(ierr);
  	ierr = MatSetSizes(slvr->A,PETSC_DECIDE,PETSC_DECIDE,NA,NA);CHKERRQ(ierr);
  	ierr = MatSeqAIJSetPreallocation(slvr->A,3*Nv,NULL);CHKERRQ(ierr);
  	ierr = MatSetFromOptions(slvr->A);CHKERRQ(ierr);
  	ierr = MatSetUp(slvr->A);CHKERRQ(ierr);

  	//Create Solver Contexts
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&slvr->ksp);CHKERRQ(ierr);
    /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
    */
    ierr = KSPSetOperators(slvr->ksp,slvr->A,slvr->A);CHKERRQ(ierr);
    ierr = KSPSetType(slvr->ksp,KSPBCGS);CHKERRQ(ierr);
    // ILU Precond
    ierr = KSPGetPC(slvr->ksp,&slvr->pc);CHKERRQ(ierr);
    ierr = PCSetType(slvr->pc,PCILU);CHKERRQ(ierr);
    ierr = PCFactorSetFill(slvr->pc,3.0);CHKERRQ(ierr);
    ierr = PCFactorSetLevels(slvr->pc,1);CHKERRQ(ierr);
    ierr = PCFactorSetAllowDiagonalFill(slvr->pc,PETSC_TRUE);CHKERRQ(ierr);
    // ierr = PCFactorSetUseInPlace(slvr->pc,PETSC_TRUE);CHKERRQ(ierr);

    PetscReal div_tol = 1e12;
    PetscReal abs_tol = 1e-15;
    ierr = KSPSetTolerances(slvr->ksp,1e-10,abs_tol,div_tol,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetNormType(slvr->ksp,KSP_NORM_UNPRECONDITIONED);CHKERRQ(ierr);
    /*
        Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */
    ierr = KSPSetFromOptions(slvr->ksp);CHKERRQ(ierr);
    ierr = PCSetFromOptions(slvr->pc);CHKERRQ(ierr);

  return ierr;
}



