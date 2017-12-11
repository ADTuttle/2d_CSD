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
PetscInt Ind_nx(PetscInt x,PetscInt y,PetscInt ion,PetscInt comp, PetscInt nx)
{
    return Nv*(nx*y+x)+ion*Nc+comp;
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
	PetscReal c[Ni*Nc];
	PetscReal phi[Nc];
	PetscReal alpha[Nc];
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
	PetscReal vm = phi[phi_index(0,0,0)]-phi[phi_index(0,0,2)]; //neuronal membrane potential
    PetscReal vmg = phi[phi_index(0,0,1)]-phi[phi_index(0,0,2)]; //glial membrane potential

    //compute neuronal Cl concentration (since neuron has only leak conductance, must be at reversal potential for Cl)

    c[c_index(0,0,0,2)] = c[c_index(0,0,2,2)]*exp(vm);
    //set glial Cl concentration equal to neuronal Cl concentration
    c[c_index(0,0,1,2)] = c[c_index(0,0,0,2)];

    // struct FluxPoPetscInt *flux;
    // flux = (struct FluxPoint*)malloc(sizeof(struct FluxPoint));

    //compute cotransporter permeability so that glial Cl is at rest
    mclin(flux,c_index(0,0,1,2),pClLeakg,-1,c[c_index(0,0,1,2)],c[c_index(0,0,2,2)],vmg,0);
    con_vars->pNaKCl = -flux->mflux[c_index(0,0,1,2)]/2/log(c[c_index(0,0,1,0)]*c[c_index(0,0,1,1)]*c[c_index(0,0,1,2)]*c[c_index(0,0,1,2)]/(c[c_index(0,0,2,0)]*c[c_index(0,0,2,1)]*c[c_index(0,0,2,2)]*c[c_index(0,0,2,2)]));
    PetscReal NaKCl = -flux->mflux[c_index(0,0,1,2)]/2;

    //compute gating variables
    gatevars_update(gate_vars,state_vars,0,1);

    //compute K channel currents (neuron)
    PetscReal pKGHK = pKDR*gate_vars->gKDR[0]+pKA*gate_vars->gKA[0];
    //Initialize the KGHK flux
    mcGoldman(flux,c_index(0,0,0,1),pKGHK,1,c[c_index(0,0,0,1)],c[c_index(0,0,Nc-1,1)],vm,0);
    //Add the KLeak flux to it
    mclin(flux,c_index(0,0,0,1),pKLeak,1,c[c_index(0,0,0,1)],c[c_index(0,0,Nc-1,1)],vm,1);

    //compute neuronal ATPase value
    con_vars->Imax = flux->mflux[c_index(0,0,0,1)]*(pow(1+mK/c[c_index(0,0,Nc-1,1)],2)*pow(1+mNa/c[c_index(0,0,0,0)],3))/2;


    //compute neuronal sodium currents and leak permeability value
    PetscReal pNaGHK = pNaT*gate_vars->gNaT[0]+pNaP*gate_vars->gNaP[0];
    mcGoldman(flux,c_index(0,0,0,0),pNaGHK,1,c[c_index(0,0,0,0)],c[c_index(0,0,Nc-1,0)],vm,0);
    PetscReal Ipump = npump*con_vars->Imax/(pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,0,0)]),3));
    con_vars->pNaLeak = (-flux->mflux[c_index(0,0,0,0)]-3*Ipump)/(log(c[c_index(0,0,0,0)]/c[c_index(0,0,Nc-1,0)])+vm);

    //compute K channel currents (glial)
    PetscReal pKLinG = pKIR*inwardrect(c[c_index(0,0,1,1)],c[c_index(0,0,Nc-1,1)],vmg)*pKLeakadjust;
    mclin(flux,c_index(0,0,1,1),pKLinG,1,c[c_index(0,0,1,1)],c[c_index(0,0,Nc-1,1)],vmg,0);
    flux->mflux[c_index(0,0,1,1)] += NaKCl;

    //compute glial ATPase value
    con_vars->Imaxg = flux->mflux[c_index(0,0,1,1)]*pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,1,0)]),3)/2;

    //compute glial sodium current and leak permeability value
    PetscReal Ipumpg = glpump*con_vars->Imaxg/(pow((1+mK/c[c_index(0,0,Nc-1,1)]),2)*pow((1+mNa/c[c_index(0,0,1,0)]),3));
    con_vars->pNaLeakg = (-NaKCl-3*Ipumpg)/(log(c[c_index(0,0,1,0)]/c[c_index(0,0,Nc-1,0)])+vmg);

    //Compute resting organic anion amounts and average valences
    //set extracellular organic anion amounts and valence to ensure electroneutrality
    con_vars->ao[Nc-1] = 5e-4;
    PetscReal cmphi[Nc];
    PetscReal osmotic;
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
  	
	 PetscReal zetaadjust = 1; //modify glial permeability  
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



    return;
}

void initialize_data(struct SimState *state_vars,struct SimState *state_vars_past, struct GateType* gate_vars,struct ConstVars* con_vars,struct Solver *slvr,struct FluxData *flux)
{
	PetscReal reltol = 1e-11;
	PetscReal tol = reltol*array_max(state_vars->c,(size_t)Nx*Ny*Nc*Ni);
  	PetscReal rsd = 1.0;
  	PetscReal *cp;
  	cp = (PetscReal *)malloc(sizeof(PetscReal)*Nx*Ny*Ni*Nc);
  	
    //Compute Gating variables
    //compute gating variables
    gatevars_update(gate_vars,state_vars,0,1);

  	//Initialize and comput the excitation (it's zeros here)
  	struct ExctType *gexct;
  	gexct = (struct ExctType*)malloc(sizeof(struct ExctType));
  	excitation(gexct,texct+1);
  	PetscInt k = 0;
  	PetscReal dt_temp = 0.1;
    // PetscReal dt_temp = 0.01;
  	
  	while(rsd>tol && dt_temp*k<10)
  	{
    	memcpy(cp,state_vars->c,sizeof(PetscReal)*Nx*Ny*Ni*Nc);
    	newton_solve(state_vars,state_vars_past, dt_temp, gate_vars, gexct,con_vars,slvr,flux);
    	gatevars_update(gate_vars,state_vars,dt_temp*1e3,0);
    	rsd = array_diff_max(state_vars->c,cp,(size_t)Nx*Ny*Nc*Ni)/dt_temp;
        printf("Init_Data rsd: %.10e, Tol: %.10e\n",rsd,tol);
    	k++;
	}
  	
  	free(cp);
  	free(gexct);
	if(rsd>1e-7)
  	{
    	fprintf(stderr, "Did not converge! Aborting...\n");
    	exit(EXIT_FAILURE); /* indicate failure.*/
	} else
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
    ierr = VecCreate(PETSC_COMM_WORLD,&slvr->Q);CHKERRQ(ierr);
    ierr = VecSetType(slvr->Q,VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(slvr->Q,PETSC_DECIDE,NA);CHKERRQ(ierr);
    ierr = VecDuplicate(slvr->Q,&slvr->Res);CHKERRQ(ierr);

  	//Create Matrix
    //Get number of nonzeros in each row
    int nnz[NA];
    Get_Nonzero_in_Rows(nnz);
    //Construct matrix using that
  	// ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,NA,NA,5*Nv,nnz,&slvr->A);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&slvr->A);CHKERRQ(ierr);
  	ierr = MatSetType(slvr->A,MATSEQAIJ);CHKERRQ(ierr);
  	ierr = MatSetSizes(slvr->A,PETSC_DECIDE,PETSC_DECIDE,NA,NA);CHKERRQ(ierr);
  	ierr = MatSeqAIJSetPreallocation(slvr->A,5*Nv,nnz);CHKERRQ(ierr);
  	// ierr = MatSetFromOptions(slvr->A);CHKERRQ(ierr);
  	ierr = MatSetUp(slvr->A);CHKERRQ(ierr);

    //Initialize Space
    ierr = initialize_jacobian(slvr->A); CHKERRQ(ierr);
    ierr = MatSetOption(slvr->A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  	//Create Solver Contexts

    ierr = KSPCreate(PETSC_COMM_WORLD,&slvr->ksp);CHKERRQ(ierr);

    //Iterative Solver Types
//    ierr = KSPSetType(slvr->ksp,KSPPREONLY);CHKERRQ(ierr);
//    ierr = KSPSetType(slvr->ksp,KSPBCGS);CHKERRQ(ierr);
//    ierr = KSPSetType(slvr->ksp,KSPRICHARDSON);CHKERRQ(ierr);
     ierr = KSPSetType(slvr->ksp,KSPGMRES);CHKERRQ(ierr);

    ierr = KSPGetPC(slvr->ksp,&slvr->pc);CHKERRQ(ierr);


//    /*
    // Multigrid Precond

    ierr = Initialize_PCMG(slvr->pc,slvr->A); CHKERRQ(ierr);
//    /*
//    ierr = PCSetType(slvr->pc,PCLU);CHKERRQ(ierr);
//    ierr = PCFactorSetMatSolverPackage(slvr->pc, MATSOLVERSUPERLU); CHKERRQ(ierr);
//     */
    // ILU Precond
    /*
    ierr = PCSetType(slvr->pc,PCILU);CHKERRQ(ierr);
    ierr = PCFactorSetFill(slvr->pc,3.0);CHKERRQ(ierr);
    ierr = PCFactorSetLevels(slvr->pc,1);CHKERRQ(ierr);
    ierr = PCFactorSetAllowDiagonalFill(slvr->pc,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(slvr->pc,MATORDERINGRCM); CHKERRQ(ierr);
    ierr = PCFactorSetReuseOrdering(slvr->pc,PETSC_TRUE); CHKERRQ(ierr);
    */

    PetscReal div_tol = 1e12;
    // PetscReal abs_tol = 1e-13;
    // PetscReal rel_tol = 1e-10;
    PetscReal abs_tol = 1e-12;
    PetscReal rel_tol = 1e-8;
    ierr = KSPSetTolerances(slvr->ksp,rel_tol,abs_tol,div_tol,PETSC_DEFAULT);CHKERRQ(ierr);
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

void Get_Nonzero_in_Rows(int *nnz)
{
    //Make sure nnz is initialized to zero
    for(int i=0;i<NA;i++)
    {
        nnz[i]=0;
    }
    int ind = 0;
    int x,y,comp,ion;
    //Ionic concentration equations
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    //Electrodiffusion contributions
                    if(x<Nx-1)
                    {
                        nnz[Ind_1(x+1,y,ion,comp)]++; //Ind_1(x,y,ion,comp)
                        ind++;
                        //Right c with left phi (-Fph0x)
                        nnz[Ind_1(x+1,y,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                        ind++;
                    }
                    if(x>0)
                    {
                        //left c with right c (-Fc1x)
                        nnz[Ind_1(x-1,y,ion,comp)]++;//Ind_1(x,y,ion,comp)
                        ind++;
                        //Left c with right phi (-Fph1x)
                        nnz[Ind_1(x-1,y,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                        ind++;
                    }
                    if(y<Ny-1)
                    {
                        // Upper c with lower c (-Fc0y)
                        nnz[Ind_1(x,y+1,ion,comp)]++;//Ind_1(x,y,ion,comp);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        nnz[Ind_1(x,y+1,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                        ind++;
                    }
                    if(y>0)
                    {
                        //Lower c with Upper c (-Fc1y)
                        nnz[Ind_1(x,y-1,ion,comp)]++;//Ind_1(x,y,ion,comp)
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        nnz[Ind_1(x,y-1,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                        ind++;
                    }
      
                    //membrane current contributions
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    nnz[Ind_1(x,y,ion,Nc-1)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                    // C Intra with C Extra
                    nnz[Ind_1(x,y,ion,comp)]++;//Ind_1(x,y,ion,Nc-1)
                    ind++;
                    // C Extracellular with Phi Inside
                    nnz[Ind_1(x,y,ion,Nc-1)]++;//Ind_1(x,y,Ni,comp)
                    ind++;
                    // C Intra with Phi Extra
                    nnz[Ind_1(x,y,ion,comp)]++;//Ind_1(x,y,Ni,Nc-1)
                    ind++;
                    //Volume terms
                    //C extra with intra alpha
                    nnz[Ind_1(x,y,ion,Nc-1)]++;//Ind_1(x,y,Ni+1,comp)
                    ind++;
                    //C intra with intra alpha
                    nnz[Ind_1(x,y,ion,comp)]++;//Ind_1(x,y,Ni+1,comp)
                    ind++;
                    //Same compartment terms
                    // c with c
                    nnz[Ind_1(x,y,ion,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                     // c with phi
                    nnz[Ind_1(x,y,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                if(x<Nx-1)
                {
                    // Right c with left c (-Fc0x)
                    nnz[Ind_1(x+1,y,ion,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                    //Right c with left phi (-Fph0x)
                    nnz[Ind_1(x+1,y,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                    ind++;
                }
                if(x>0)
                {
                    //left c with right c (-Fc1x)
                    nnz[Ind_1(x-1,y,ion,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                    //Left c with right phi (-Fph1x)
                    nnz[Ind_1(x-1,y,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                    ind++;
                }
                if(y<Ny-1)
                {
                    // Upper c with lower c (-Fc0y)
                    nnz[Ind_1(x,y+1,ion,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    nnz[Ind_1(x,y+1,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                    ind++;
                }
                if(y>0)
                {
                    //Lower c with Upper c (-Fc1y)
                    nnz[Ind_1(x,y-1,ion,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    nnz[Ind_1(x,y-1,ion,comp)]++;//Ind_1(x,y,Ni,comp)
                    ind++;
                }
            
                //Membrane current contribution
                //Add bath contributions
                //Insert extracell to extracell parts
                // c with c
                nnz[Ind_1(x,y,ion,Nc-1)]++;//Ind_1(x,y,ion,Nc-1)
                ind++;
                 // c with phi
                nnz[Ind_1(x,y,ion,Nc-1)]++;//Ind_1(x,y,Ni,Nc-1)
                ind++;
            }
        }
    }
    //Electroneutrality charge-capcitance condition
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            //electroneutral-concentration entries
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    //Phi with C entries
                    nnz[Ind_1(x,y,Ni,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                }
                //Phi with C extracellular one
                comp = Nc-1;
                nnz[Ind_1(x,y,Ni,comp)]++;//Ind_1(x,y,ion,comp)
                ind++;

            }
            //electroneutrality-voltage entries
            //extraphi with extra phi
            nnz[Ind_1(x,y,Ni,Nc-1)]++;//Ind_1(x,y,Ni,Nc-1)
            ind++;
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extra phi with intra phi
                nnz[Ind_1(x,y,Ni,Nc-1)]++;//Ind_1(x,y,Ni,comp)
                ind++;
                // Intra phi with Extraphi
                nnz[Ind_1(x,y,Ni,comp)]++;//Ind_1(x,y,Ni,Nc-1)
                ind++;
                //Intra phi with Intra phi
                nnz[Ind_1(x,y,Ni,comp)]++;//Ind_1(x,y,Ni,comp)
                ind++;
                //Extra phi with intra-Volume
                nnz[Ind_1(x,y,Ni,Nc-1)]++;//Ind_1(x,y,Ni+1,comp)
                ind++;
                //Intra phi with Intra Vol
                nnz[Ind_1(x,y,Ni,comp)]++;//Ind_1(x,y,Ni+1,comp)
                ind++;
            }
        }
    }
    //water flow
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(comp=0;comp<Nc-1;comp++)
            {
                //Water flow volume fraction entries
                //Volume to Volume
                nnz[Ind_1(x,y,Ni+1,comp)]++;//Ind_1(x,y,Ni+1,comp)
                ind++;
                //Off diagonal (from aNc=1-sum(ak))
                for (PetscInt l=0; l<comp; l++)
                {
                    nnz[Ind_1(x,y,Ni+1,comp)]++;//Ind_1(x,y,Ni+1,l)
                    ind++; 
                }
                for (PetscInt l=comp+1; l<Nc-1; l++)
                {
                    nnz[Ind_1(x,y,Ni+1,comp)]++;//Ind_1(x,y,Ni+1,l)
                    ind++;
                }
                for (ion=0; ion<Ni; ion++)
                {
                  //Volume to extra c
                    nnz[Ind_1(x,y,Ni+1,comp)]++;//Ind_1(x,y,ion,Nc-1)
                    ind++;
                  //Volume to intra c
                    nnz[Ind_1(x,y,Ni+1,comp)]++;//Ind_1(x,y,ion,comp)
                    ind++;
                }
            }
        }
    }
    return;
}


PetscErrorCode initialize_jacobian(Mat Jac) {
    printf("Initializing Jacobian Memory\n");
    PetscErrorCode ierr;
    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    //Ionic concentration equations
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    //Electrodiffusion contributions
                    if(x<Nx-1)
                    {
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0)
                    {
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1)
                    {
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0)
                    {
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    //membrane current contributions
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,ion,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume terms
                    //C extra with intra alpha
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni+1,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //C intra with intra alpha
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni+1,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                if(x<Nx-1)
                {
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0)
                {
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1)
                {
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0)
                {
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution

                //Membrane current contribution
                //Add bath contributions
                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,ion,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
        }
    }

    //Electroneutrality charge-capcitance condition
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            //electroneutral-concentration entries
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    //Phi with C entries
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                }
                //Phi with C extracellular one
                comp = Nc-1;
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES); CHKERRQ(ierr);
                ind++;

            }
            //electroneutrality-voltage entries
            //extraphi with extra phi
            ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
            ind++;
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extra phi with intra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // Intra phi with Extraphi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra phi with Intra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Extra phi with intra-Volume
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni+1,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra phi with Intra Vol
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni+1,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
        }
    }
    //water flow
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(comp=0;comp<Nc-1;comp++)
            {
                //Water flow volume fraction entries
                //Volume to Volume
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Off diagonal (from aNc=1-sum(ak))
                for (PetscInt l=0; l<comp; l++)
                {
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,l),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (PetscInt l=comp+1; l<Nc-1; l++)
                {
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,l),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (ion=0; ion<Ni; ion++)
                {
                    //Volume to extra c
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,ion,Nc-1),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume to intra c
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,ion,comp),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
            }
        }
    }

//    for (int i = 0; i < NA; i++)
//    {
//        ierr = MatSetValue(Jac,i,i,0,INSERT_VALUES); CHKERRQ(ierr);
//    }
    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return ierr;
}

PetscErrorCode Create_Restriction(Mat R,PetscInt nx, PetscInt ny)
{
    PetscErrorCode  ierr;
    int x,y,ion,comp;
    for(x=1;x<nx/2-1;x++) {
        for (y = 1; y < ny / 2-1; y++) {
            //Restriction for concentrations
            for (ion = 0; ion < Ni; ion++) {
                for (comp = 0; comp < Nc; comp++) {
                    //Center point
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 4,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);

                    //Up/down/left/right
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, ion, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, ion, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, ion, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, ion, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);

                    //Four diagonals
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y - 1, ion, comp, nx),
                                       1.0 / 16, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y + 1, ion, comp, nx),
                                       1.0 / 16, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y - 1, ion, comp, nx),
                                       1.0 / 16, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, ion, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y + 1, ion, comp, nx),
                                       1.0 / 16, INSERT_VALUES);
                    CHKERRQ(ierr);

                }
            }
            //Restriction for Voltage
            for (comp = 0; comp < Nc; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 4,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni, comp, nx), 1.0 / 8,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni, comp, nx), 1.0 / 8,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni, comp, nx), 1.0 / 8,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni, comp, nx), 1.0 / 8,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y - 1, Ni, comp, nx),
                                   1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y + 1, Ni, comp, nx),
                                   1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y - 1, Ni, comp, nx),
                                   1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y + 1, Ni, comp, nx),
                                   1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);

            }
            //Restriction for Volume
            for (comp = 0; comp < Nc-1; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 4,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                                   1.0 / 8, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 8, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                                   1.0 / 8, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 8, INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 16, INSERT_VALUES);
                CHKERRQ(ierr);

            }

        }
    }
    x=0;
    for(y=1;y<ny/2-1;y++) {
        //Restriction for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                //Center point
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y-1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y+1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

                //Two diagonals
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y-1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y+1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Restriction for Voltage
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y-1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y+1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }
        //Restriction for Volume
        for(comp=0;comp<Nc-1;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y-1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y+1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }

    }
    x=nx/2-1;
    for(y=1;y<ny/2-1;y++) {
        //Restriction for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                //Center point
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y-1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y+1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

                //Two diagonals
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y-1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y+1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Restriction for Voltage
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y-1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y+1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }
        //Restriction for Volume
        for(comp=0;comp<Nc-1;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y-1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y+1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }

    }
    y=0;
    for(x=1;x<nx/2-1;x++) {
        //Restriction for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                //Center point
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y+1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

                //Two diagonals
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y+1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y+1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Restriction for Voltage
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y+1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }
        //Restriction for Volume
        for(comp=0;comp<Nc-1;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y+1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }

    }
    y=nx/2-1;
    for(x=1;x<nx/2-1;x++) {
        //Restriction for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                //Center point
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y-1,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y,ion,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

                //Two diagonals
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y-1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y-1,ion,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Restriction for Voltage
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y-1,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y,Ni,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }
        //Restriction for Volume
        for(comp=0;comp<Nc-1;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),1.0/3,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y-1,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y,Ni+1,comp,nx),1.0/6,INSERT_VALUES); CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni+1,comp,nx),1.0/12,INSERT_VALUES); CHKERRQ(ierr);

        }

    }
    x=0;y=0;
    //Restriction for concentrations
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y+1,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Two diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y+1,ion,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

        }
    }
    //Restriction for Voltage
    for(comp=0;comp<Nc;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y+1,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    //Restriction for Volume
    for(comp=0;comp<Nc-1;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y+1,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y+1,Ni+1,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    x=nx/2-1;y=0;
    //Restriction for concentrations
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y+1,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Two diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y+1,ion,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

        }
    }
    //Restriction for Voltage
    for(comp=0;comp<Nc;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y+1,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    //Restriction for Volume
    for(comp=0;comp<Nc-1;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y+1,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y+1,Ni+1,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    x=0;y=ny/2-1;
    //Restriction for concentrations
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y-1,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Two diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x+1,2*y-1,ion,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

        }
    }
    //Restriction for Voltage
    for(comp=0;comp<Nc;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y-1,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    //Restriction for Volume
    for(comp=0;comp<Nc-1;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y-1,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x+1,2*y-1,Ni+1,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    x=nx/2-1;y=ny/2-1;
    //Restriction for concentrations
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc;comp++) {
            //Center point
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x,2*y-1,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y,ion,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

            //Two diagonals
            ierr = MatSetValue(R,Ind_nx(x,y,ion,comp,nx/2),Ind_nx(2*x-1,2*y-1,ion,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

        }
    }
    //Restriction for Voltage
    for(comp=0;comp<Nc;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x,2*y-1,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y,Ni,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    //Restriction for Volume
    for(comp=0;comp<Nc-1;comp++) {
        //Center point
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y,ion,comp,nx),4.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Up/down/left/right
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x,2*y-1,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y,Ni+1,comp,nx),2.0/9,INSERT_VALUES); CHKERRQ(ierr);

        //Four diagonals
        ierr = MatSetValue(R,Ind_nx(x,y,Ni+1,comp,nx/2),Ind_nx(2*x-1,2*y-1,Ni+1,comp,nx),1.0/9,INSERT_VALUES); CHKERRQ(ierr);

    }
    ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return ierr;
}
PetscErrorCode Create_Interpolation(Mat R,PetscInt nx, PetscInt ny)
{
    PetscErrorCode  ierr;
    int x,y,ion,comp;
    for(x=0;x<nx-1;x++)
    {
        for(y=0;y<ny-1;y++)
        {
            //Interpolation for concentrations
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc;comp++) {
                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x+1,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                    ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y+1,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x+1,y,ion,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y+1,ion,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                    ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x+1,y+1,ion,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);

                    ierr = MatSetValue(R,Ind_nx(2*x,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);

                }
            }
            //Interpolation for Voltage
            for(comp=0;comp<Nc;comp++) {
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x+1,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y+1,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x+1,y,Ni,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y+1,Ni,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x+1,y+1,Ni,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
            }
            //Interpolation for Volume
            for(comp=0;comp<Nc-1;comp++) {

                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x+1,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y+1,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x+1,y,Ni+1,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y+1,Ni+1,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x+1,y+1,Ni+1,comp,nx),0.25,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
            }

        }
    }
    x=nx-1;
    for(y=0;y<ny-1;y++) {
        //Interpolation for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y+1,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y+1,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Interpolation for Voltage
        for(comp=0;comp<Nc;comp++) {
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y+1,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y+1,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
        }
        //Interpolation for Volume
        for(comp=0;comp<Nc-1;comp++) {

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y+1,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y+1,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    y=ny-1;
    for(x=0;x<nx-1;x++) {
        //Interpolation for concentrations
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc;comp++) {
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x+1,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
                ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x+1,y,ion,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

                ierr = MatSetValue(R,Ind_nx(2*x,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);

            }
        }
        //Interpolation for Voltage
        for(comp=0;comp<Nc;comp++) {
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x+1,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x+1,y,Ni,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
        }
        //Interpolation for Volume
        for(comp=0;comp<Nc-1;comp++) {

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x+1,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x+1,y,Ni+1,comp,nx),0.5,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    x=nx-1;y=ny-1;
    //Interpolation for concentrations
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc;comp++) {
            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

            ierr = MatSetValue(R,Ind_nx(2*x,2*y,ion,comp,2*nx),Ind_nx(x,y,ion,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);

        }
    }
    //Interpolation for Voltage
    for(comp=0;comp<Nc;comp++) {
        ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni,comp,2*nx),Ind_nx(x,y,Ni,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
    }
    //Interpolation for Volume
    for(comp=0;comp<Nc-1;comp++) {

        ierr = MatSetValue(R,Ind_nx(2*x+1,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x+1,2*y+1,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1.0,INSERT_VALUES); CHKERRQ(ierr);

        ierr = MatSetValue(R,Ind_nx(2*x,2*y,Ni+1,comp,2*nx),Ind_nx(x,y,Ni+1,comp,nx),1,INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//    MatView(R,PETSC_VIEWER_STDOUT_SELF);

    return ierr;
}

PetscErrorCode Initialize_PCMG(PC pc,Mat A)
{
    PetscErrorCode  ierr;

    KSP coarse_ksp,sksp;
    PC coarse_pc,spc;

    PetscInt nlevels = 3;
    PetscInt nx = Nx;
    PetscInt ny = Ny;
    Mat R,P;
    ierr = PCSetType(pc,PCMG); CHKERRQ(ierr);
    ierr = PCSetOperators(pc,A,A);CHKERRQ(ierr);
    ierr = PCMGSetType(pc,PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH); CHKERRQ(ierr);
    PCMGSetLevels(pc,nlevels,PETSC_NULL);
    ierr = PCMGSetCycleType(pc,	PC_MG_CYCLE_V); CHKERRQ(ierr);

    //Set coarse solve
    ierr = PCMGGetCoarseSolve(pc,&coarse_ksp); CHKERRQ(ierr);
    ierr = KSPSetType(coarse_ksp,KSPPREONLY);CHKERRQ(ierr);
    ierr = KSPGetPC(coarse_ksp,&coarse_pc);CHKERRQ(ierr);
    ierr = PCSetType(coarse_pc,PCLU); CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(coarse_pc, MATSOLVERSUPERLU); CHKERRQ(ierr);


/*
    //Make Multigrid operators at each level
    ierr = PCMGGetSmoother(pc,0,&sksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(sksp,A,A); CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD,&Aprev); CHKERRQ(ierr);
    ierr = MatDuplicate(A,MAT_COPY_VALUES,&Aprev); CHKERRQ(ierr);
    for (int i=1; i<nlevels; i++) {

        //Make restriction operator
        ierr = MatCreate(PETSC_COMM_WORLD,&R);CHKERRQ(ierr);
        ierr = MatSetType(R,MATSEQAIJ);CHKERRQ(ierr);
        ierr = MatSetSizes(R,PETSC_DECIDE,PETSC_DECIDE,nx/2*ny/2*Nv,nx*ny*Nv);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(R,9,NULL);CHKERRQ(ierr);
        ierr = MatSetOption(R,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);

        ierr = Create_Restriction(R,nx,ny); CHKERRQ(ierr);

        ierr = PCMGSetRestriction(pc,i,R); CHKERRQ(ierr);


        //Make interpolation opperator
        ierr = MatCreate(PETSC_COMM_WORLD,&P);CHKERRQ(ierr);
        ierr = MatSetType(P,MATSEQAIJ);CHKERRQ(ierr);
        ierr = MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,nx*ny*Nv,nx/2*ny/2*Nv);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(P,4,NULL);CHKERRQ(ierr);
        ierr = MatSetOption(P,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);

        ierr = Create_Interpolation(P,nx/2,ny/2); CHKERRQ(ierr);

        ierr = PCMGSetInterpolation(pc,i,P); CHKERRQ(ierr);

        ierr = PCMGGetSmoother(pc,i,&sksp); CHKERRQ(ierr);

        ierr = MatCreate(PETSC_COMM_WORLD,&Ai);CHKERRQ(ierr);
        ierr = MatSetType(Ai,MATSEQAIJ);CHKERRQ(ierr);
        ierr = MatSetSizes(Ai,PETSC_DECIDE,PETSC_DECIDE,nx/2*ny/2*Nv,nx/2*ny/2*Nv);CHKERRQ(ierr);
        ierr = MatSetUp(Ai); CHKERRQ(ierr);

        ierr = MatMatMatMult(R,Aprev,P,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Ai); CHKERRQ(ierr);
        ierr = KSPSetOperators(sksp,Ai,Ai); CHKERRQ(ierr);
        ierr = MatDuplicate(Ai,MAT_COPY_VALUES,&Aprev); CHKERRQ(ierr);

        ierr = MatDestroy(&Ai);CHKERRQ(ierr);
        ierr = MatDestroy(&R); CHKERRQ(ierr);
        ierr = MatDestroy(&P); CHKERRQ(ierr);

        nx = nx/2;
        ny = ny/2;
    }
    */



    //Make restriction operators
    for (int i=nlevels-1; i>0; i--) {
        ierr = MatCreate(PETSC_COMM_WORLD,&R);CHKERRQ(ierr);
        ierr = MatSetType(R,MATSEQAIJ);CHKERRQ(ierr);
        ierr = MatSetSizes(R,PETSC_DECIDE,PETSC_DECIDE,nx/2*ny/2*Nv,nx*ny*Nv);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(R,9,NULL);CHKERRQ(ierr);
        ierr = MatSetOption(R,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);

        ierr = Create_Restriction(R,nx,ny); CHKERRQ(ierr);

        ierr = PCMGSetRestriction(pc,i,R); CHKERRQ(ierr);
//        ierr = PCMGSetResidual(pc,i,PCMGResidualDefault,A); CHKERRQ(ierr);

        ierr = MatDestroy(&R); CHKERRQ(ierr);

        nx = nx/2;
        ny = ny/2;
    }

    // Make interpolation Ops
    for (int i=1; i<nlevels; i++) {
        ierr = MatCreate(PETSC_COMM_WORLD,&P);CHKERRQ(ierr);
        ierr = MatSetType(P,MATSEQAIJ);CHKERRQ(ierr);
        ierr = MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,nx*2*ny*2*Nv,nx*ny*Nv);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(P,4,NULL);CHKERRQ(ierr);
        ierr = MatSetOption(P,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);

        ierr = Create_Interpolation(P,nx,ny); CHKERRQ(ierr);

        ierr = PCMGSetInterpolation(pc,i,P); CHKERRQ(ierr);

        ierr = MatDestroy(&P); CHKERRQ(ierr);

        nx = nx*2;
        ny = ny*2;
    }


    //Modify the smoother (default KSP is chebyshev with SOR)
    for(int i=1;i<nlevels;i++){
        ierr = PCMGGetSmoother(pc,i,&sksp); CHKERRQ(ierr);
        ierr = KSPGetPC(sksp,&spc); CHKERRQ(ierr);


        //Smoother KSP
//        ierr = KSPSetType(sksp,KSPRICHARDSON); CHKERRQ(ierr);
//        ierr = KSPRichardsonSetScale(sksp,0.5); CHKERRQ(ierr);
        ierr = KSPSetType(sksp,KSPBCGS); CHKERRQ(ierr);
//        ierr = KSPSetType(sksp,KSPPREONLY); CHKERRQ(ierr);
        //Smoother Precond
        /*
        ierr = PCSetType(spc,PCSOR); CHKERRQ(ierr);
        ierr = PCSORSetSymmetric(spc,SOR_LOCAL_BACKWARD_SWEEP); CHKERRQ(ierr);
        ierr = PCSORSetIterations(spc,10,10); CHKERRQ(ierr);
        ierr = PCSORSetOmega(spc,1.999);
         */
//        ierr = PCSetType(spc, PCJACOBI);CHKERRQ(ierr);
//        ierr = PCJacobiSetType(spc,PC_JACOBI_ROWMAX); CHKERRQ(ierr);

//        ierr = PCSetType(spc, PCASM); CHKERRQ(ierr);

//        /*
        ierr = PCSetType(spc,PCILU);CHKERRQ(ierr);
        ierr = PCFactorSetFill(spc,3.0);CHKERRQ(ierr);
        ierr = PCFactorSetLevels(spc,1);CHKERRQ(ierr);
        ierr = PCFactorSetAllowDiagonalFill(spc,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PCFactorSetMatOrderingType(spc,MATORDERINGRCM); CHKERRQ(ierr);
        ierr = PCFactorSetReuseOrdering(spc,PETSC_TRUE); CHKERRQ(ierr);
//        */
    }

    return ierr;
}