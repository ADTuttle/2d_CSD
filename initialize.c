#include "constants.h"
#include "functions.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

PetscInt c_index(PetscInt x,PetscInt y,PetscInt comp,PetscInt ion,PetscInt Nx)
{
	return Nc*Ni* (Nx * y + x) + comp*Ni+ion; 
}
PetscInt phi_index(PetscInt x,PetscInt y,PetscInt comp,PetscInt Nx)
{
	return Nc* (Nx * y + x) + comp; 
}
PetscInt al_index(PetscInt x,PetscInt y,PetscInt comp,PetscInt Nx)
{
	return (Nc-1)* (Nx * y + x) + comp; 
}
PetscInt xy_index(PetscInt x,PetscInt y,PetscInt Nx)
{
	return Nx*y+x;
}
PetscInt Ind_1(PetscInt x,PetscInt y,PetscInt ion,PetscInt comp,PetscInt Nx)
{
    return Nv*(Nx*y+x)+ion*Nc+comp;
}

PetscInt Ind_nx(PetscInt x,PetscInt y,PetscInt ion,PetscInt comp, PetscInt nx)
{
    return Nv*(nx*y+x)+ion*Nc+comp;
}

PetscErrorCode init_simstate(Vec state,struct SimState *state_vars,struct AppCtx *user)
{
    PetscErrorCode ierr;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    //Setup indices
    int x,y,comp,ion;
    PetscInt c_ind[Nx*Ny*Nc*Ni];
    PetscInt phi_ind[Nx*Ny*Nc];
    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){
            for(comp=0;comp<Nc;comp++)
            {
                for(ion=0;ion<Ni;ion++)
                {
                    c_ind[c_index(x,y,comp,ion,Nx)] = Ind_1(x,y,ion,comp,Nx);
                }
                phi_ind[phi_index(x,y,comp,Nx)] = Ind_1(x,y,Ni,comp,Nx);
            }
        }
    }
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,Nx*Ny*Ni*Nc,c_ind,PETSC_COPY_VALUES,&state_vars->c_ind); CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,Nx*Ny*Nc,phi_ind,PETSC_COPY_VALUES,&state_vars->phi_ind); CHKERRQ(ierr);

    if(!separate_vol) {
        PetscInt al_ind[Nx*Ny*(Nc-1)];
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc - 1; comp++) {
                    al_ind[al_index(x, y, comp,Nx)] = Ind_1(x, y, Ni + 1, comp,Nx);
                }
            }
        }
        ierr = ISCreateGeneral(PETSC_COMM_WORLD, Nx * Ny * (Nc - 1), al_ind, PETSC_COPY_VALUES, &state_vars->al_ind);
        CHKERRQ(ierr);
    }
    else{
        state_vars->alpha = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*(Nc-1));
    }
    extract_subarray(state,state_vars,1);
    return ierr;
}
PetscErrorCode extract_subarray(Vec state,struct SimState *state_vars,int option)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[2], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;

    //If option=2 then we are extracting a slow variable
    //If option=1 then we are extraction a fast variable
    if(option==2) {
        //Get Slow vars
        ierr = VecGetSubVector(state_vars->v, state_vars->c_ind, &state_vars->c_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->c_vec, &state_vars->c);CHKERRQ(ierr);

        ierr = VecGetSubVector(state_vars->v, state_vars->phi_ind, &state_vars->phi_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->phi_vec, &state_vars->phi);CHKERRQ(ierr);
        //Get Fast variables
        ierr = VecGetSubVector(state, state_vars->c_ind, &state_vars->c_fast_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->c_fast_vec, &state_vars->c_fast);CHKERRQ(ierr);

        ierr = VecGetSubVector(state, state_vars->phi_ind, &state_vars->phi_fast_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->phi_fast_vec, &state_vars->phi_fast);CHKERRQ(ierr);
    } else {
        //Get Slow vars
        ierr = VecGetSubVector(state, state_vars->c_ind, &state_vars->c_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->c_vec, &state_vars->c);CHKERRQ(ierr);

        ierr = VecGetSubVector(state, state_vars->phi_ind, &state_vars->phi_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->phi_vec, &state_vars->phi);CHKERRQ(ierr);
        //Get Fast variables
        ierr = VecGetSubVector(state_vars->v_fast, state_vars->c_ind, &state_vars->c_fast_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->c_fast_vec, &state_vars->c_fast);CHKERRQ(ierr);

        ierr = VecGetSubVector(state_vars->v_fast, state_vars->phi_ind, &state_vars->phi_fast_vec);CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->phi_fast_vec, &state_vars->phi_fast);CHKERRQ(ierr);
    }

    if(!separate_vol) {
        ierr = VecGetSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[2], 0, 0, 0, 0);
    }

    return ierr;

}

PetscErrorCode restore_subarray(Vec state,struct SimState *state_vars,int option)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[3], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;

    //If option=2 then we are extracting a slow variable
    //If option=1 then we are extraction a fast variable
    if(option==2) {
        ierr = VecRestoreArray(state_vars->c_vec, &state_vars->c);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state_vars->v, state_vars->c_ind, &state_vars->c_vec);CHKERRQ(ierr);


        ierr = VecRestoreArray(state_vars->phi_vec, &state_vars->phi);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state_vars->v, state_vars->phi_ind, &state_vars->phi_vec);CHKERRQ(ierr);

        ierr = VecRestoreArray(state_vars->c_fast_vec, &state_vars->c_fast);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->c_ind, &state_vars->c_fast_vec);CHKERRQ(ierr);


        ierr = VecRestoreArray(state_vars->phi_fast_vec, &state_vars->phi_fast);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->phi_ind, &state_vars->phi_fast_vec);CHKERRQ(ierr);
    } else{
        ierr = VecRestoreArray(state_vars->c_vec, &state_vars->c);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->c_ind, &state_vars->c_vec);CHKERRQ(ierr);


        ierr = VecRestoreArray(state_vars->phi_vec, &state_vars->phi);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->phi_ind, &state_vars->phi_vec);CHKERRQ(ierr);

        ierr = VecRestoreArray(state_vars->c_fast_vec, &state_vars->c_fast);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state_vars->v_fast, state_vars->c_ind, &state_vars->c_fast_vec);CHKERRQ(ierr);


        ierr = VecRestoreArray(state_vars->phi_fast_vec, &state_vars->phi_fast);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state_vars->v_fast, state_vars->phi_ind, &state_vars->phi_fast_vec);CHKERRQ(ierr);
    }

    if(!separate_vol) {
        ierr = VecRestoreArray(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        state_vars->alpha = NULL;
    }

    state_vars->c = NULL;
    state_vars->phi = NULL;
    state_vars->c_fast = NULL;
    state_vars->phi_fast = NULL;
    if(Profiling_on) {
        PetscLogEventEnd(event[3], 0, 0, 0, 0);
    }

    return ierr;

}
PetscErrorCode copy_simstate(Vec current_state,struct SimState *state_vars_past,int option)
{
    PetscErrorCode ierr;
    if(option==2){
        ierr = VecCopy(current_state, state_vars_past->v_fast);CHKERRQ(ierr);
        ierr = extract_subarray(state_vars_past->v_fast, state_vars_past,option);CHKERRQ(ierr);
    } else {
        ierr = VecCopy(current_state, state_vars_past->v);CHKERRQ(ierr);
        ierr = extract_subarray(state_vars_past->v, state_vars_past,option);CHKERRQ(ierr);
    }
    return ierr;
}

void init(Vec state,struct SimState *state_vars,struct AppCtx*user)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    extract_subarray(state,state_vars,1);
	for(PetscInt x=0;x<Nx;x++)
	{
		for(PetscInt y=0;y<Ny;y++)
		{
			//initial volume fractions
			state_vars->alpha[al_index(x,y,0,Nx)]=alphao[0];
			state_vars->alpha[al_index(x,y,1,Nx)]=alphao[1]; 
			//initial voltages (dimensionless)
			state_vars->phi[phi_index(x,y,0,Nx)] = -70/RTFC; //neuronal voltage
			state_vars->phi[phi_index(x,y,1,Nx)] = -85/RTFC; //glial voltage
			state_vars->phi[phi_index(x,y,2,Nx)] = -0/RTFC; //extracell voltage
			//initial concentrations in mmol/cm^3=1e-3 mmol/l
			state_vars->c[c_index(x,y,0,0,Nx)] = 10e-3;     //neuronal Na concentration
			state_vars->c[c_index(x,y,1,0,Nx)] = 10e-3;      //glial Na concentration
			state_vars->c[c_index(x,y,2,0,Nx)] = 140e-3;     //extracellular Na concentration
			state_vars->c[c_index(x,y,0,1,Nx)] = 130e-3;     //neuronal K concentration
			state_vars->c[c_index(x,y,1,1,Nx)] = 130e-3;     //glial K concentration
			state_vars->c[c_index(x,y,2,1,Nx)] = 3.4e-3;     //extracellular K concentration
			state_vars->c[c_index(x,y,0,2,Nx)] = 10e-3;       //neuronal Cl concentration
			state_vars->c[c_index(x,y,1,2,Nx)] = 10e-3; 		//glial Cl concentraion
			state_vars->c[c_index(x,y,2,2,Nx)] = 120e-3;       //143.5e-3%extracellular Cl

		}
	}
    restore_subarray(state,state_vars,1);
    VecSet(user->state_vars->v_fast,0);
    VecSet(user->state_vars_past->v_fast,0);
}

void init_arrays(struct AppCtx*user)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    //Flux quantities
    user->flux->mflux = (PetscReal*) malloc(Nx*Ny*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdci = (PetscReal*) malloc(Nx*Ny*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdce = (PetscReal*) malloc(Nx*Ny*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdphim = (PetscReal*) malloc(Nx*Ny*Ni*Nc*sizeof(PetscReal));
    user->flux->wflow = (PetscReal*) malloc(Nx*Ny*(Nc-1)*sizeof(PetscReal));
    user->flux->dwdpi = (PetscReal*) malloc(Nx*Ny*(Nc-1)*sizeof(PetscReal));
    user->flux->dwdal = (PetscReal*) malloc(Nx*Ny*(Nc-1)*sizeof(PetscReal));

    //Gating variables
    user->gate_vars->mNaT = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->hNaT = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->gNaT = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->mNaP = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->hNaP = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->gNaP = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->mKDR = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->gKDR = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->mKA = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->hKA = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gate_vars->gKA = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));


    //Excitation
    user->gexct->pNa = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gexct->pK = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    user->gexct->pCl = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));

    //Constant params
    user->con_vars->ao = (PetscReal*) malloc(Nc*sizeof(PetscReal));
    user->con_vars->zo = (PetscReal*) malloc(Nc*sizeof(PetscReal));
    user->con_vars->zeta1 = (PetscReal*) malloc((Nc-1)*sizeof(PetscReal));
    user->con_vars->zetaalpha = (PetscReal*) malloc((Nc-1)*sizeof(PetscReal));

    //Diffusion in ctx
    user->Dcs = (PetscReal*) malloc(Nx*Ny*Ni*Nc*2*sizeof(PetscReal));
    user->Dcb = (PetscReal*) malloc(Nx*Ny*Ni*Nc*2*sizeof(PetscReal));
}
void set_params(Vec state,struct SimState* state_vars,struct ConstVars* con_vars,struct GateType* gate_vars,struct FluxData *flux,struct AppCtx*user)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    extract_subarray(state,state_vars,1);
	//Everything that follows will asume spatially uniform
	//At rest state
	PetscReal c[Ni*Nc];
	PetscReal phi[Nc];
	PetscReal alpha[Nc];
	for(PetscInt comp=0;comp<Nc;comp++)
	{
		for(PetscInt ion=0;ion<Ni;ion++)
		{
			c[c_index(0,0,comp,ion,Nx)]=state_vars->c[c_index(0,0,comp,ion,Nx)];
		}
		phi[phi_index(0,0,comp,Nx)]=state_vars->phi[phi_index(0,0,comp,Nx)];
		if(comp<Nc-1)
		{
			alpha[al_index(0,0,comp,Nx)]=state_vars->alpha[al_index(0,0,comp,Nx)];
		}
		else
		{
			alpha[al_index(0,0,comp,Nx)]=1; //Adding in the extracellular vol
			for(PetscInt i=0;i<Nc-1;i++)
			{
				alpha[al_index(0,0,comp,Nx)]-=alpha[al_index(0,0,i,Nx)];
			}
		}
	}
	PetscReal vm = phi[phi_index(0,0,0,Nx)]-phi[phi_index(0,0,2,Nx)]; //neuronal membrane potential
    PetscReal vmg = phi[phi_index(0,0,1,Nx)]-phi[phi_index(0,0,2,Nx)]; //glial membrane potential

    //compute neuronal Cl concentration (since neuron has only leak conductance, must be at reversal potential for Cl)

    c[c_index(0,0,0,2,Nx)] = c[c_index(0,0,2,2,Nx)]*exp(vm);
    //set glial Cl concentration equal to neuronal Cl concentration
    c[c_index(0,0,1,2,Nx)] = c[c_index(0,0,0,2,Nx)];

    // struct FluxPoPetscInt *flux;
    // flux = (struct FluxPoint*)malloc(sizeof(struct FluxPoint));

    //compute cotransporter permeability so that glial Cl is at rest
    mclin(flux,c_index(0,0,1,2,Nx),pClLeakg,-1,c[c_index(0,0,1,2,Nx)],c[c_index(0,0,2,2,Nx)],vmg,0);
    con_vars->pNaKCl = -flux->mflux[c_index(0,0,1,2,Nx)]/2/log(c[c_index(0,0,1,0,Nx)]*c[c_index(0,0,1,1,Nx)]*c[c_index(0,0,1,2,Nx)]*c[c_index(0,0,1,2,Nx)]/(c[c_index(0,0,2,0,Nx)]*c[c_index(0,0,2,1,Nx)]*c[c_index(0,0,2,2,Nx)]*c[c_index(0,0,2,2,Nx)]));
    PetscReal NaKCl = -flux->mflux[c_index(0,0,1,2,Nx)]/2;

    //compute gating variables
    gatevars_update(gate_vars,state_vars,0,user,1);

    //compute K channel currents (neuron)
    PetscReal pKGHK = pKDR*gate_vars->gKDR[0]+pKA*gate_vars->gKA[0];
    //Initialize the KGHK flux
    mcGoldman(flux,c_index(0,0,0,1,Nx),pKGHK,1,c[c_index(0,0,0,1,Nx)],c[c_index(0,0,Nc-1,1,Nx)],vm,0);
    //Add the KLeak flux to it
    mclin(flux,c_index(0,0,0,1,Nx),pKLeak,1,c[c_index(0,0,0,1,Nx)],c[c_index(0,0,Nc-1,1,Nx)],vm,1);

    //compute neuronal ATPase value
    con_vars->Imax = flux->mflux[c_index(0,0,0,1,Nx)]*(pow(1+mK/c[c_index(0,0,Nc-1,1,Nx)],2)*pow(1+mNa/c[c_index(0,0,0,0,Nx)],3))/2;


    //compute neuronal sodium currents and leak permeability value
    PetscReal pNaGHK = pNaT*gate_vars->gNaT[0]+pNaP*gate_vars->gNaP[0];
    mcGoldman(flux,c_index(0,0,0,0,Nx),pNaGHK,1,c[c_index(0,0,0,0,Nx)],c[c_index(0,0,Nc-1,0,Nx)],vm,0);
    PetscReal Ipump = npump*con_vars->Imax/(pow((1+mK/c[c_index(0,0,Nc-1,1,Nx)]),2)*pow((1+mNa/c[c_index(0,0,0,0,Nx)]),3));
    con_vars->pNaLeak = (-flux->mflux[c_index(0,0,0,0,Nx)]-3*Ipump)/(log(c[c_index(0,0,0,0,Nx)]/c[c_index(0,0,Nc-1,0,Nx)])+vm);

    //compute K channel currents (glial)
    PetscReal pKLinG = pKIR*inwardrect(c[c_index(0,0,1,1,Nx)],c[c_index(0,0,Nc-1,1,Nx)],vmg)*pKLeakadjust;
    mclin(flux,c_index(0,0,1,1,Nx),pKLinG,1,c[c_index(0,0,1,1,Nx)],c[c_index(0,0,Nc-1,1,Nx)],vmg,0);
    flux->mflux[c_index(0,0,1,1,Nx)] += NaKCl;

    //compute glial ATPase value
    con_vars->Imaxg = flux->mflux[c_index(0,0,1,1,Nx)]*pow((1+mK/c[c_index(0,0,Nc-1,1,Nx)]),2)*pow((1+mNa/c[c_index(0,0,1,0,Nx)]),3)/2;

    //compute glial sodium current and leak permeability value
    PetscReal Ipumpg = glpump*con_vars->Imaxg/(pow((1+mK/c[c_index(0,0,Nc-1,1,Nx)]),2)*pow((1+mNa/c[c_index(0,0,1,0,Nx)]),3));
    con_vars->pNaLeakg = (-NaKCl-3*Ipumpg)/(log(c[c_index(0,0,1,0,Nx)]/c[c_index(0,0,Nc-1,0,Nx)])+vmg);

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
      	osmotic += c[c_index(0,0,Nc-1,ion,Nx)]-c[c_index(0,0,k,ion,Nx)];
      }
      con_vars->ao[k] = alpha[k]*(con_vars->ao[Nc-1]/alpha[Nc-1]+osmotic);
      //set average valence to ensure electroneutrality
      con_vars->zo[k] = (-cz(c,z,0,0,k,user)*alpha[k]+cmphi[k])/con_vars->ao[k];
    }
    con_vars->zo[Nc-1] = (-cz(c,z,0,0,Nc-1,user)*alpha[Nc-1]-cmphi[Nc-1])/con_vars->ao[Nc-1];
    //Copy the point data to vectors.
    //Only needed for uniform data
    for(PetscInt x=0;x<Nx;x++)
    {
    	for(PetscInt y=0;y<Ny;y++)
    	{
    		//Gating variables (already set in gatevars_update)
    		//We changed c_index(0,0,0/1,2,Nx), neuronal/glial Cl.
    		state_vars->c[c_index(x,y,0,2,Nx)] = c[c_index(0,0,0,2,Nx)];
    		state_vars->c[c_index(x,y,1,2,Nx)] = c[c_index(0,0,1,2,Nx)];
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

    restore_subarray(state,state_vars,1);


    return;
}

void initialize_data(Vec current_state,struct AppCtx *user)
{

    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

	PetscReal convtol = 1e-11;
    extract_subarray(current_state,user->state_vars,1);
	PetscReal tol = convtol*array_max(user->state_vars->c,(size_t)Nx*Ny*Nc*Ni);
  	PetscReal rsd = 1.0;
  	PetscReal *cp;
  	cp = (PetscReal *)malloc(sizeof(PetscReal)*Nx*Ny*Ni*Nc);
    //Compute Gating variables
    //compute gating variables
    gatevars_update(user->gate_vars,user->state_vars,0,user,1);
    restore_subarray(current_state,user->state_vars,1);

  	//Initialize and compute the excitation (it's zeros here)
  	excitation(user,texct+1);
  	PetscInt k = 0;
    user->dt = 0.1;
//  	PetscReal dt_temp = 0.1;
    // PetscReal dt_temp = 0.01;
  	
  	while(rsd>tol && user->dt*k<10)
  	{
        extract_subarray(current_state,user->state_vars,1);
    	memcpy(cp,user->state_vars->c,sizeof(PetscReal)*Nx*Ny*Ni*Nc);
        //Save the "current" aka past state
        restore_subarray(user->state_vars_past->v,user->state_vars_past,1);
        copy_simstate(user->state_vars->v_fast,user->state_vars_past,2);
        restore_subarray(user->state_vars_past->v,user->state_vars_past,2);
        copy_simstate(current_state,user->state_vars_past,1);
        if(separate_vol) {
            //Update volume
            volume_update(user->state_vars, user->state_vars_past, user);
        }
        //compute diffusion coefficients
        diff_coef(user->Dcs,user->state_vars->alpha,1,user);
        //Bath diffusion
        diff_coef(user->Dcb,user->state_vars->alpha,Batheps,user);
        restore_subarray(current_state,user->state_vars,1);

//    	newton_solve(current_state,user);
        //Solve fast var
        SNESSolve(user->fast_slvr->snes,NULL,user->state_vars->v_fast);
        //Solve slow var
        SNESSolve(user->slvr->snes,NULL,current_state);


        //Update gating variables
        extract_subarray(current_state,user->state_vars,1);
        gatevars_update(user->gate_vars,user->state_vars,user->dt*1e3,user,0);

        //Update Excitation
    	rsd = array_diff_max(user->state_vars->c,cp,(size_t)Nx*Ny*Nc*Ni)/user->dt;
        restore_subarray(current_state,user->state_vars,1);
        printf("Init_Data rsd: %.10e, Tol: %.10e\n",rsd,tol);
    	k++;
	}
  	
  	free(cp);
	if(rsd>1e-7)
  	{
    	fprintf(stderr, "Did not converge! Aborting...\n");
    	exit(EXIT_FAILURE); /* indicate failure.*/
	} else
  	{
    	return;
	}
	
}


PetscErrorCode initialize_petsc(struct Solver *slvr,int argc, char **argv,struct AppCtx *user)
{
    PetscErrorCode ierr;
	//Init Petsc
	PetscInitialize(&argc,&argv,(char*)0,NULL);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&slvr->size);CHKERRQ(ierr);
    //    Get Nx, Ny, and dt from options if possible

    user->Nx = 50;
    user->Ny = 50;
    user->dt =0.1;
//    user->dt = 1.0/1000.0/2.0;
    PetscOptionsGetInt(NULL,NULL,"-Nx",&user->Nx,NULL);
    PetscOptionsGetInt(NULL,NULL,"-Ny",&user->Ny,NULL);
    PetscOptionsGetReal(NULL,NULL,"-dt",&user->dt,NULL);


    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    user->dx = Lx/Nx;
    user->dy = Ly/Ny;
    user->NA = Nx*Ny*Nv;//total number of unknowns
    user->Nz = (Ni*Nc*(4*(Nx-1)*Ny+4*(Ny-1)*Nx+2*Nx*Ny)+Ni*(Nc-1)*6*Nx*Ny+(Nc*Ni+1)*Nx*Ny+(Nc-1)*(6*Nx*Ny+Nx*Ny*(Nc-2)+Ni*2*Nx*Ny)); //number of nonzeros in Jacobian

    PetscInt NA = user->NA;
    PetscInt Nz = user->Nz;


  	//Create Vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&slvr->Q);CHKERRQ(ierr);
    ierr = VecSetType(slvr->Q,VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(slvr->Q,PETSC_DECIDE,NA);CHKERRQ(ierr);
    ierr = VecDuplicate(slvr->Q,&slvr->Res);CHKERRQ(ierr);

  	//Create Matrix
    //Get number of nonzeros in each row
    int *nnz = (int*) malloc(sizeof(int)*NA);
    Get_Nonzero_in_Rows(nnz,user);
    //Construct matrix using that
  	// ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,NA,NA,5*Nv,nnz,&slvr->A);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&slvr->A);CHKERRQ(ierr);
  	ierr = MatSetType(slvr->A,MATSEQAIJ);CHKERRQ(ierr);
  	ierr = MatSetSizes(slvr->A,PETSC_DECIDE,PETSC_DECIDE,NA,NA);CHKERRQ(ierr);
  	ierr = MatSeqAIJSetPreallocation(slvr->A,5*Nv,nnz);CHKERRQ(ierr);
  	// ierr = MatSetFromOptions(slvr->A);CHKERRQ(ierr);
  	ierr = MatSetUp(slvr->A);CHKERRQ(ierr);


    //Initialize Space

    ierr = initialize_jacobian(slvr->A,user); CHKERRQ(ierr);
    ierr = MatSetOption(slvr->A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  	//Create Solver Contexts
    ierr = SNESCreate(PETSC_COMM_WORLD,&slvr->snes); CHKERRQ(ierr);


//    ierr = KSPCreate(PETSC_COMM_WORLD,&slvr->ksp);CHKERRQ(ierr);
    ierr = SNESGetKSP(slvr->snes,&slvr->ksp); CHKERRQ(ierr);
    ierr = KSPGetPC(slvr->ksp,&slvr->pc);CHKERRQ(ierr);

    if(Linear_Diffusion){
        if(use_en_deriv){
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_linear_deriv, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_linear_deriv, user);
            CHKERRQ(ierr);
        } else {
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_linear_algebraic, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_linear_algebraic, user);
            CHKERRQ(ierr);
        }
    }else {
        if (separate_vol && use_en_deriv) {
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_no_vol, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_no_vol, user);
            CHKERRQ(ierr);
        } else if (!separate_vol && !use_en_deriv) {
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_algebraic, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_algebraic, user);
            CHKERRQ(ierr);
        } else if (separate_vol && !use_en_deriv) {
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_algebraic_no_vol, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_algebraic_no_vol, user);
            CHKERRQ(ierr);
        } else {
            //Set Function eval
            ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual, user);
            CHKERRQ(ierr);
            //Set Jacobian eval
            ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian, user);
            CHKERRQ(ierr);
        }
    }

    ierr = SNESSetFunction(slvr->snes,slvr->Res,calc_residual_slow,user); CHKERRQ(ierr);
    ierr = SNESSetJacobian(slvr->snes,slvr->A,slvr->A,calc_jacobian_slow,user); CHKERRQ(ierr);
    //Set SNES types
    ierr = SNESSetType(slvr->snes,SNESNEWTONLS); CHKERRQ(ierr);
//    ierr = SNESSetType(slvr->snes,SNESNEWTONTR); CHKERRQ(ierr);


//    ierr = KSPSetOperators(slvr->ksp,slvr->A,slvr->A);CHKERRQ(ierr);
//    ierr = KSPSetType(slvr->ksp,KSPPREONLY);CHKERRQ(ierr);
//     ierr = KSPSetType(slvr->ksp,KSPBCGS);CHKERRQ(ierr);

    //Gmres type methods
//     ierr = KSPSetType(slvr->ksp,KSPGMRES);CHKERRQ(ierr);


    if( Nx%2==0 && Nx>64){
        printf("Using FGMRES+Multigrid(Richardson+SOR)\n");
        ierr = KSPSetType(slvr->ksp,KSPFGMRES);CHKERRQ(ierr);

        //Multigrid precond
        ierr = Initialize_PCMG(slvr->pc,slvr->A,user); CHKERRQ(ierr);
    }else {
        printf("Using DGMRES+ILU\n");
        //DGmres KSP
        ierr = KSPSetType(slvr->ksp, KSPDGMRES);CHKERRQ(ierr);
        ierr = KSPGMRESSetRestart(slvr->ksp, 40);CHKERRQ(ierr);
        ierr = PetscOptionsSetValue(NULL, "-ksp_dgmres_eigen", "10");CHKERRQ(ierr);
        ierr = PetscOptionsSetValue(NULL, "-ksp_dgmres_max_eigen", "100");CHKERRQ(ierr);
        ierr = PetscOptionsSetValue(NULL, "-ksp_dgmres_force", "");CHKERRQ(ierr);

        // ILU Precond
        ierr = PCSetType(slvr->pc, PCILU);CHKERRQ(ierr);
        ierr = PCFactorSetFill(slvr->pc, 3.0);CHKERRQ(ierr);
        ierr = PCFactorSetLevels(slvr->pc, 1);CHKERRQ(ierr);
        ierr = PCFactorSetAllowDiagonalFill(slvr->pc, PETSC_TRUE);CHKERRQ(ierr);
        ierr = PCFactorSetMatOrderingType(slvr->pc, MATORDERINGNATURAL);CHKERRQ(ierr);
    }




    //LU Direct solve
    /*
    ierr = PCSetType(slvr->pc,PCLU);CHKERRQ(ierr);
    ierr = KSPSetPC(slvr->ksp,slvr->pc);CHKERRQ(ierr);
     ierr = PCFactorSetMatSolverPackage(slvr->pc, MATSOLVERSUPERLU); CHKERRQ(ierr);
    */

//     ierr = PCFactorSetUseInPlace(slvr->pc,PETSC_TRUE);CHKERRQ(ierr);
    /*
    PetscReal div_tol = 1e12;
    PetscReal abs_tol = 1e-16;
    PetscReal rel_tol = 1e-14;
//    PetscReal abs_tol = 1e-12;
//    PetscReal rel_tol = 1e-8;
    ierr = KSPSetTolerances(slvr->ksp,rel_tol,abs_tol,div_tol,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetNormType(slvr->ksp,KSP_NORM_UNPRECONDITIONED);CHKERRQ(ierr);
    */

    /*
        Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
    */
    ierr = SNESSetFromOptions(slvr->snes);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(slvr->ksp);CHKERRQ(ierr);
    ierr = PCSetFromOptions(slvr->pc);CHKERRQ(ierr);


    free(nnz);

    return ierr;
}

PetscErrorCode initialize_fast_petsc(struct Solver *slvr,int argc, char **argv,struct AppCtx *user)
{
    PetscErrorCode ierr;

    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;


    //Create Vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&slvr->Q);CHKERRQ(ierr);
    ierr = VecSetType(slvr->Q,VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(slvr->Q,PETSC_DECIDE,Nx*Ny*Nc*(Ni+1));CHKERRQ(ierr);
    ierr = VecDuplicate(slvr->Q,&slvr->Res);CHKERRQ(ierr);

    //Create Matrix
    ierr = MatCreate(PETSC_COMM_WORLD,&slvr->A);CHKERRQ(ierr);
    ierr = MatSetType(slvr->A,MATSEQAIJ);CHKERRQ(ierr);
    ierr = MatSetSizes(slvr->A,PETSC_DECIDE,PETSC_DECIDE,Nx*Ny*Nc*(Ni+1),Nx*Ny*Nc*(Ni+1));CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(slvr->A,12,NULL);CHKERRQ(ierr);
    // ierr = MatSetFromOptions(slvr->A);CHKERRQ(ierr);
    ierr = MatSetUp(slvr->A);CHKERRQ(ierr);


    //Create Solver Contexts
    ierr = SNESCreate(PETSC_COMM_WORLD,&slvr->snes); CHKERRQ(ierr);

    ierr = SNESGetKSP(slvr->snes,&slvr->ksp); CHKERRQ(ierr);
    ierr = KSPGetPC(slvr->ksp,&slvr->pc);CHKERRQ(ierr);

    //Set Function eval
    ierr = SNESSetFunction(slvr->snes, slvr->Res, calc_residual_fast, user);CHKERRQ(ierr);
    //Set Jacobian eval
    ierr = SNESSetJacobian(slvr->snes, slvr->A, slvr->A, calc_jacobian_fast, user);CHKERRQ(ierr);

    //Set SNES types
    ierr = SNESSetType(slvr->snes,SNESNEWTONLS); CHKERRQ(ierr);
    ierr = KSPSetType(slvr->ksp,KSPPREONLY);CHKERRQ(ierr);
//    ierr = KSPSetType(slvr->ksp,KSPGMRES);CHKERRQ(ierr);
    //LU Direct solve

//    /*
    ierr = PCSetType(slvr->pc,PCLU);CHKERRQ(ierr);
    ierr = KSPSetPC(slvr->ksp,slvr->pc);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(slvr->pc, MATSOLVERSUPERLU); CHKERRQ(ierr);
//    */
    /*
    // ILU Precond
    ierr = PCSetType(slvr->pc, PCILU);CHKERRQ(ierr);
    ierr = PCFactorSetFill(slvr->pc, 3.0);CHKERRQ(ierr);
    ierr = PCFactorSetLevels(slvr->pc, 1);CHKERRQ(ierr);
    ierr = PCFactorSetAllowDiagonalFill(slvr->pc, PETSC_TRUE);CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(slvr->pc, MATORDERINGNATURAL);CHKERRQ(ierr);
 */
//     ierr = PCFactorSetUseInPlace(slvr->pc,PETSC_TRUE);CHKERRQ(ierr);
    /*
    PetscReal div_tol = 1e12;
    PetscReal abs_tol = 1e-16;
    PetscReal rel_tol = 1e-14;
//    PetscReal abs_tol = 1e-12;
//    PetscReal rel_tol = 1e-8;
    ierr = KSPSetTolerances(slvr->ksp,rel_tol,abs_tol,div_tol,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetNormType(slvr->ksp,KSP_NORM_UNPRECONDITIONED);CHKERRQ(ierr);
    */
    ierr = SNESSetFromOptions(slvr->snes);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(slvr->ksp);CHKERRQ(ierr);
    ierr = PCSetFromOptions(slvr->pc);CHKERRQ(ierr);

    return ierr;
}
void Get_Nonzero_in_Rows(int *nnz,struct AppCtx *user)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt NA = user->NA;
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
                        nnz[Ind_1(x+1,y,ion,comp,Nx)]++; //Ind_1(x,y,ion,comp,Nx)
                        ind++;
                        //Right c with left phi (-Fph0x)
                        nnz[Ind_1(x+1,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                        ind++;
                        nnz[Ind_1(x+1,y,Ni,comp,Nx)]++;
                        ind++;
                    }
                    if(x>0)
                    {
                        //left c with right c (-Fc1x)
                        nnz[Ind_1(x-1,y,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                        ind++;
                        //Left c with right phi (-Fph1x)
                        nnz[Ind_1(x-1,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                        ind++;
                        nnz[Ind_1(x-1,y,Ni,comp,Nx)]++;
                        ind++;
                    }
                    if(y<Ny-1)
                    {
                        // Upper c with lower c (-Fc0y)
                        nnz[Ind_1(x,y+1,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        nnz[Ind_1(x,y+1,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                        ind++;
                        nnz[Ind_1(x,y+1,Ni,comp,Nx)]++;
                        ind++;
                    }
                    if(y>0)
                    {
                        //Lower c with Upper c (-Fc1y)
                        nnz[Ind_1(x,y-1,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        nnz[Ind_1(x,y-1,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                        ind++;
                        nnz[Ind_1(x,y-1,Ni,comp,Nx)]++;
                        ind++;
                    }
      
                    //membrane current contributions
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    nnz[Ind_1(x,y,ion,Nc-1,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                    // C Intra with C Extra
                    nnz[Ind_1(x,y,ion,comp,Nx)]++;//Ind_1(x,y,ion,Nc-1,Nx)
                    ind++;
                    // C Extracellular with Phi Inside
                    nnz[Ind_1(x,y,ion,Nc-1,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                    // C Intra with Phi Extra
                    nnz[Ind_1(x,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,Nc-1,Nx)
                    ind++;
                    if(!separate_vol) {
                        //Volume terms
                        //C extra with intra alpha
                        nnz[Ind_1(x, y, ion, Nc - 1,Nx)]++;//Ind_1(x,y,Ni+1,comp,Nx)
                        ind++;
                        //C intra with intra alpha
                        nnz[Ind_1(x, y, ion, comp,Nx)]++;//Ind_1(x,y,Ni+1,comp,Nx)
                        ind++;
                    }
                    //Same compartment terms
                    // c with c
                    nnz[Ind_1(x,y,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                     // c with phi
                    nnz[Ind_1(x,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    nnz[Ind_1(x,y,Ni,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    //IntraPhi with c extra(volt eqn)
                    nnz[Ind_1(x,y,Ni,comp,Nx)]++;//Ind_1(x,y,ion,Nc-1,Nx)

                    //Extra-Phi with intra-c (voltage eqn)
                    nnz[Ind_1(x,y,Ni,Nc-1,Nx)]++;//Ind_1(x,y,ion,comp,Nx)


                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                if(x<Nx-1)
                {
                    // Right c with left c (-Fc0x)
                    nnz[Ind_1(x+1,y,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                    //Right c with left phi (-Fph0x)
                    nnz[Ind_1(x+1,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                    nnz[Ind_1(x+1,y,Ni,comp,Nx)]++;
                    ind++;
                }
                if(x>0)
                {
                    //left c with right c (-Fc1x)
                    nnz[Ind_1(x-1,y,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                    //Left c with right phi (-Fph1x)
                    nnz[Ind_1(x-1,y,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                    nnz[Ind_1(x-1,y,Ni,comp,Nx)]++;
                    ind++;
                }
                if(y<Ny-1)
                {
                    // Upper c with lower c (-Fc0y)
                    nnz[Ind_1(x,y+1,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    nnz[Ind_1(x,y+1,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                    nnz[Ind_1(x,y+1,Ni,comp,Nx)]++;
                    ind++;
                }
                if(y>0)
                {
                    //Lower c with Upper c (-Fc1y)
                    nnz[Ind_1(x,y-1,ion,comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    nnz[Ind_1(x,y-1,ion,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                    nnz[Ind_1(x,y-1,Ni,comp,Nx)]++;
                    ind++;
                }
            
                //Membrane current contribution
                //Add bath contributions
                //Insert extracell to extracell parts
                // c with c
                nnz[Ind_1(x,y,ion,Nc-1,Nx)]++;//Ind_1(x,y,ion,Nc-1,Nx)
                ind++;
                 // c with phi
                nnz[Ind_1(x,y,ion,Nc-1,Nx)]++;//Ind_1(x,y,Ni,Nc-1,Nx)
                ind++;
                //Extra phi with c (volt eqn)
                nnz[Ind_1(x,y,Ni,Nc-1,Nx)]++;
                ind++;
            }
            //Derivative of charge-capacitance
            for(comp=0;comp<Nc-1;comp++){
                if(x<Nx-1)
                {
                    //Right phi with left phi (-Fph0x)
                    nnz[Ind_1(x+1,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                }
                if(x>0)
                {
                    //Left phi with right phi (-Fph1x)
                    nnz[Ind_1(x-1,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                }
                if(y<Ny-1)
                {
                    //Upper phi with lower phi (-Fph0y)
                    nnz[Ind_1(x,y+1,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                }
                if(y>0)
                {
                    //Lower phi with upper phi (-Fph1y)
                    nnz[Ind_1(x,y-1,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                    ind++;
                }
                //Intra-phi with Intra-phi
                nnz[Ind_1(x,y,Ni,comp,Nx)]++;
                ind++;
                //Intra-phi with extra-phi
                nnz[Ind_1(x,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,Nc-1,Nx)
                ind++;
            }
            //Extracellular terms
            comp = Nc-1;
            if(x<Nx-1)
            {
                //Right phi with left phi (-Fph0x)
                nnz[Ind_1(x+1,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                ind++;
            }
            if(x>0)
            {
                //Left phi with right phi (-Fph1x)
                nnz[Ind_1(x-1,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                ind++;
            }
            if(y<Ny-1)
            {
                //Upper phi with lower phi (-Fph0y)
                nnz[Ind_1(x,y+1,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                ind++;
            }
            if(y>0)
            {
                //Lower phi with upper phi (-Fph1y)
                nnz[Ind_1(x,y-1,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
                ind++;
            }
            for(int k=0;k<Nc-1;k++){

                //Extra-phi with Intra-phi
                nnz[Ind_1(x,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,k,Nx)
                ind++;
            }
            //extra-phi with extra-phi
            nnz[Ind_1(x,y,Ni,comp,Nx)]++;//Ind_1(x,y,Ni,comp,Nx)
            ind++;
        }
    }
    if(!separate_vol) {
        //water flow
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Water flow volume fraction entries
                    //Volume to Volume
                    nnz[Ind_1(x, y, Ni + 1, comp,Nx)]++;//Ind_1(x,y,Ni+1,comp,Nx)
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for (PetscInt l = 0; l < comp; l++) {
                        nnz[Ind_1(x, y, Ni + 1, comp,Nx)]++;//Ind_1(x,y,Ni+1,l,Nx)
                        ind++;
                    }
                    for (PetscInt l = comp + 1; l < Nc - 1; l++) {
                        nnz[Ind_1(x, y, Ni + 1, comp,Nx)]++;//Ind_1(x,y,Ni+1,l,Nx)
                        ind++;
                    }
                    for (ion = 0; ion < Ni; ion++) {
                        //Volume to extra c
                        nnz[Ind_1(x, y, Ni + 1, comp,Nx)]++;//Ind_1(x,y,ion,Nc-1,Nx)
                        ind++;
                        //Volume to intra c
                        nnz[Ind_1(x, y, Ni + 1, comp,Nx)]++;//Ind_1(x,y,ion,comp,Nx)
                        ind++;
                    }
                }
            }
        }
    }
    printf("Nz: %d, ind: %d\n",user->Nz,ind);
    return;
}


PetscErrorCode initialize_jacobian(Mat Jac,struct AppCtx *user) {
    printf("Initializing Jacobian Memory\n");
    PetscErrorCode ierr;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions

                    if(x<Nx-1)
                    {
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv) {
                            //Right phi with left c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }


                    }
                    if(x>0)
                    {
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv) {
                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }
                    if(y<Ny-1)
                    {
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv) {
                            //Upper phi with lower c in voltage eqn
                            ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }
                    if(y>0)
                    {
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv) {
                            //Lower phi with upper c in voltage eqn
                            ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }

                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,ion,Nc-1,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,Ni,Nc-1,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if(!separate_vol) {
                        //Volume terms
                        //C extra with intra alpha
                        ierr = MatSetValue(Jac, Ind_1(x, y, ion, Nc - 1,Nx), Ind_1(x, y, Ni + 1, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //C intra with intra alpha
                        ierr = MatSetValue(Jac, Ind_1(x, y, ion, comp,Nx), Ind_1(x, y, Ni + 1, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv) {
                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, ion, Nc - 1,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                if(x<Nx-1)
                {
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv) {
                        // left Phi with right c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(x>0)
                {
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv) {
                        // left Phi with right c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(y<Ny-1)
                {
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv) {
                        // Upper Phi with lower c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(y>0)
                {
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv) {
                        // Lower Phi with upper c (voltage eqn)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,ion,Nc-1,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,Ni,Nc-1,Nx),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                if (use_en_deriv) {
                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1,Nx), Ind_1(x, y, ion, Nc - 1,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
            }
            if (use_en_deriv) {
                //Derivative of charge-capacitance
                for (comp = 0; comp < Nc - 1; comp++) {
                    if (x < Nx - 1) {
                        //Right phi with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (x > 0) {
                        //Left phi with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y < Ny - 1) {
                        //Upper phi with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y > 0) {
                        //Lower phi with upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, Nc - 1,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc - 1;
                if (x < Nx - 1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                for (int k = 0; k < Nc - 1; k++) {
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, k,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //extra-phi with extra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }

        }
    }
    if(!use_en_deriv) {
        //Electroneutrality charge-capcitance condition
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                //electroneutral-concentration entries
                for (ion = 0; ion < Ni; ion++) {
                    for (comp = 0; comp < Nc - 1; comp++) {
                        //Phi with C entries
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc - 1;
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //electroneutrality-voltage entries

                //extraphi with extra phi
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1,Nx), Ind_1(x, y, Ni, Nc - 1,Nx), 0, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Extra phi with intra phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // Intra phi with Extraphi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, Nc - 1,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    if(!separate_vol) {
                        //Extra phi with intra-Volume
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1,Nx), Ind_1(x, y, Ni + 1, comp,Nx), 0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Intra phi with Intra Vol
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp,Nx), Ind_1(x, y, Ni + 1, comp,Nx), 0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
            }
        }
    }
    if(!separate_vol) {
        //water flow
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Water flow volume fraction entries
                    //Volume to Volume
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni + 1, comp,Nx), Ind_1(x, y, Ni + 1, comp,Nx), 0, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for (PetscInt l = 0; l < comp; l++) {
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni + 1, comp,Nx), Ind_1(x, y, Ni + 1, l,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for (PetscInt l = comp + 1; l < Nc - 1; l++) {
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni + 1, comp,Nx), Ind_1(x, y, Ni + 1, l,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for (ion = 0; ion < Ni; ion++) {
                        //Volume to extra c
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni + 1, comp,Nx), Ind_1(x, y, ion, Nc - 1,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume to intra c
                        ierr = MatSetValue(Jac, Ind_1(x, y, Ni + 1, comp,Nx), Ind_1(x, y, ion, comp,Nx), 0, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
            }
        }
    }
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
            if(!separate_vol) {
                //Restriction for Volume
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Center point
                    ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx),
                                       1.0 / 4,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);

                    //Up/down/left/right
                    ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                       Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                       Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                       Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                                       1.0 / 8, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                       Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
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
        if(!separate_vol) {
            //Restriction for Volume
            for (comp = 0; comp < Nc - 1; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 3,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);

            }
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
        if(!separate_vol) {
            //Restriction for Volume
            for (comp = 0; comp < Nc - 1; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 3,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);

            }
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
        if(!separate_vol) {
            //Restriction for Volume
            for (comp = 0; comp < Nc - 1; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 3,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);

            }
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
        if(!separate_vol) {
            //Restriction for Volume
            for (comp = 0; comp < Nc - 1; comp++) {
                //Center point
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 1.0 / 3,
                                   INSERT_VALUES);
                CHKERRQ(ierr);

                //Up/down/left/right
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                                   1.0 / 6, INSERT_VALUES);
                CHKERRQ(ierr);

                //Four diagonals
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x - 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2),
                                   Ind_nx(2 * x + 1, 2 * y - 1, Ni + 1, comp, nx), 1.0 / 12, INSERT_VALUES);
                CHKERRQ(ierr);

            }
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
    if(!separate_vol) {
        //Restriction for Volume
        for (comp = 0; comp < Nc - 1; comp++) {
            //Center point
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 4.0 / 9,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, nx),
                               1.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

        }
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
    if(!separate_vol) {
        //Restriction for Volume
        for (comp = 0; comp < Nc - 1; comp++) {
            //Center point
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 4.0 / 9,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y + 1, Ni + 1, comp, nx),
                               1.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

        }
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
    if(!separate_vol) {
        //Restriction for Volume
        for (comp = 0; comp < Nc - 1; comp++) {
            //Center point
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 4.0 / 9,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x + 1, 2 * y - 1, Ni + 1, comp, nx),
                               1.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

        }
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
    if(!separate_vol) {
        //Restriction for Volume
        for (comp = 0; comp < Nc - 1; comp++) {
            //Center point
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y, ion, comp, nx), 4.0 / 9,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            //Up/down/left/right
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x, 2 * y - 1, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y, Ni + 1, comp, nx),
                               2.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

            //Four diagonals
            ierr = MatSetValue(R, Ind_nx(x, y, Ni + 1, comp, nx / 2), Ind_nx(2 * x - 1, 2 * y - 1, Ni + 1, comp, nx),
                               1.0 / 9, INSERT_VALUES);
            CHKERRQ(ierr);

        }
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
            if(!separate_vol) {
                //Interpolation for Volume
                for (comp = 0; comp < Nc - 1; comp++) {

                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x + 1, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                    CHKERRQ(ierr);

                    ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x, y + 1, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                    CHKERRQ(ierr);

                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x, y, Ni + 1, comp, nx), 0.25, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x + 1, y, Ni + 1, comp, nx), 0.25, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x, y + 1, Ni + 1, comp, nx), 0.25, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                       Ind_nx(x + 1, y + 1, Ni + 1, comp, nx), 0.25, INSERT_VALUES);
                    CHKERRQ(ierr);

                    ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                }
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
        if(!separate_vol) {
            //Interpolation for Volume
            for (comp = 0; comp < Nc - 1; comp++) {

                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx),
                                   1.0, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx),
                                   0.5, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x, y + 1, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x, y + 1, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
            }
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
        if(!separate_vol) {
            //Interpolation for Volume
            for (comp = 0; comp < Nc - 1; comp++) {

                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx),
                                   0.5, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x + 1, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx),
                                   1.0, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);
                ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx),
                                   Ind_nx(x + 1, y, Ni + 1, comp, nx), 0.5, INSERT_VALUES);
                CHKERRQ(ierr);

                ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
            }
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
    if(!separate_vol) {
        //Interpolation for Volume
        for (comp = 0; comp < Nc - 1; comp++) {

            ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1.0,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y + 1, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1.0,
                               INSERT_VALUES);
            CHKERRQ(ierr);

            ierr = MatSetValue(R, Ind_nx(2 * x + 1, 2 * y + 1, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx),
                               1.0, INSERT_VALUES);
            CHKERRQ(ierr);

            ierr = MatSetValue(R, Ind_nx(2 * x, 2 * y, Ni + 1, comp, 2 * nx), Ind_nx(x, y, Ni + 1, comp, nx), 1,
                               INSERT_VALUES);
            CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//    MatView(R,PETSC_VIEWER_STDOUT_SELF);

    return ierr;
}

PetscErrorCode Initialize_PCMG(PC pc,Mat A,struct AppCtx*user)
{
    PetscErrorCode  ierr;

    KSP coarse_ksp,sksp;
    PC coarse_pc,spc;


    PetscInt nx = user->Nx;
    PetscInt ny = user->Ny;
    PetscInt nlevels = (PetscInt) log2(nx);
    Mat R,P;

    ierr = PCSetType(pc,PCMG); CHKERRQ(ierr);
    ierr = PCSetOperators(pc,A,A);CHKERRQ(ierr);
    ierr = PCMGSetType(pc,PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
//    ierr = PCMGSetType(pc,PC_MG_KASKADE); CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH); CHKERRQ(ierr);
    PCMGSetLevels(pc,nlevels,PETSC_NULL);
//    ierr = PCMGSetCycleType(pc,	PC_MG_CYCLE_V); CHKERRQ(ierr);
    ierr = PCMGSetCycleType(pc,PC_MG_CYCLE_W); CHKERRQ(ierr);



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
        ierr = KSPSetType(sksp,KSPRICHARDSON); CHKERRQ(ierr);
        ierr = KSPRichardsonSetScale(sksp,1.0); CHKERRQ(ierr);
//        ierr = KSPSetType(sksp,KSPBCGS); CHKERRQ(ierr);
//        ierr = KSPSetType(sksp,KSPGMRES); CHKERRQ(ierr);
//        ierr = KSPSetType(sksp,KSPPREONLY); CHKERRQ(ierr);
        //Smoother Precond
//        /*
        ierr = PCSetType(spc,PCSOR); CHKERRQ(ierr);
//        ierr = PCSORSetSymmetric(spc,SOR_LOCAL_BACKWARD_SWEEP); CHKERRQ(ierr);
//        ierr = PCSORSetSymmetric(spc,SOR_LOCAL_SYMMETRIC_SWEEP); CHKERRQ(ierr);
        ierr = PCSORSetSymmetric(spc,SOR_LOCAL_FORWARD_SWEEP); CHKERRQ(ierr);
        ierr = PCSORSetIterations(spc,2,2); CHKERRQ(ierr);
        ierr = PCSORSetOmega(spc,1.0);
//         */
//        ierr = PCSetType(spc, PCJACOBI);CHKERRQ(ierr);
//        ierr = PCJacobiSetType(spc,PC_JACOBI_ROWMAX); CHKERRQ(ierr);
        /*
        ierr = PCSetType(spc, PCASM); CHKERRQ(ierr);
        ierr = PCASMSetType(spc,PC_ASM_BASIC); CHKERRQ(ierr);
        ierr = PCASMSetLocalType(spc,PC_COMPOSITE_ADDITIVE); CHKERRQ(ierr);
//        ierr = PCASMSetLocalType(spc,PC_COMPOSITE_MULTIPLICATIVE); CHKERRQ(ierr);
         */

        /*
        ierr = PCSetType(spc,PCILU);CHKERRQ(ierr);
        ierr = PCFactorSetFill(spc,3.0);CHKERRQ(ierr);
        ierr = PCFactorSetLevels(spc,1);CHKERRQ(ierr);
        ierr = PCFactorSetAllowDiagonalFill(spc,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PCFactorSetMatOrderingType(spc,MATORDERINGRCM); CHKERRQ(ierr);
//        ierr = PCFactorSetReuseOrdering(spc,PETSC_TRUE); CHKERRQ(ierr);
        */
    }

    return ierr;
}