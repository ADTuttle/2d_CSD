
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"



int main(int argc, char **argv)
{
    printf("\n\n\nGrid size: %dx%d, with %d ions, and %d compartments.\n",Nx,Ny,Ni,Nc);
    PetscLogDouble tic,toc,full_tic,full_toc;
    //Create state_variables struct
    struct SimState *state_vars;
    state_vars=(struct SimState*)malloc(sizeof(struct SimState));
    struct SimState *state_vars_past;
    state_vars_past=(struct SimState*)malloc(sizeof(struct SimState));
    //Initialize
    init(state_vars);

    printf("Init Value: c: %f,ph: %f,al: %f\n",state_vars->c[0],state_vars->phi[10],state_vars->alpha[25]);

    //Create the constant ion channel vars
    struct ConstVars *con_vars;
    con_vars=(struct ConstVars*)malloc(sizeof(struct ConstVars));

    //Create the gating variables
    struct GateType *gate_vars;
    gate_vars = (struct GateType*) malloc(sizeof(struct GateType));
    //Create the flux structure
    struct FluxData *flux;
    flux = (struct FluxData*) malloc(sizeof(struct FluxData));

    //Set the constant variables
    set_params(state_vars,con_vars,gate_vars,flux);

    printf("Initialization Routine\n");

  	struct Solver *slvr;
  	slvr = (struct Solver*)malloc(sizeof(struct Solver));
  	PetscErrorCode ierr=0;
  	ierr = initialize_petsc(slvr,argc,argv);CHKERRQ(ierr);

	//Run Initialization routine to get to steady state
	initialize_data(state_vars,state_vars_past,gate_vars,con_vars,slvr,flux);

    printf("Beginning Main Routine \n");
    printf("\n\n\n");
    //Open file to write to
    FILE *fp;
    fp = fopen("data_csd.txt","w");
    write_data(fp,state_vars,1);
    //Create the excitation
    struct ExctType *gexct;
    gexct = (struct ExctType*)malloc(sizeof(struct ExctType));
    excitation(gexct,0);
    int count = 0;
    PetscTime(&full_tic);
    for(PetscReal t=dt;t<=Time;t+=dt)
    {
        PetscTime(&tic);
        newton_solve(state_vars,state_vars_past,dt,gate_vars,gexct,con_vars,slvr,flux);
        PetscTime(&toc);
        printf("Newton time: %f\n",toc-tic);
        //Update gating variables
        gatevars_update(gate_vars,state_vars,dt*1e3,0);
        //Update Excitation
        excitation(gexct,t);
        count++;
        if(count%krecordfreq==0) {
            printf("Time: %f, Netwon SolveTime %f\n",t,toc-tic);
            write_data(fp, state_vars, 0);
        }

    }
    PetscTime(&full_toc);
    //Close
    fclose(fp);
    printf("Finished Running. Full solve time: %.10e\n",full_toc-full_tic);
    //Free memory
    free(state_vars);free(con_vars);free(gate_vars);
    VecDestroy(&slvr->Q); VecDestroy(&slvr->Res); MatDestroy(&slvr->A);
    KSPDestroy(&slvr->ksp);
    free(slvr);
    return 0;
}

