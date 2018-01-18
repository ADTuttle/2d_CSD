
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"



int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    //Petsc Initialize
    struct Solver *slvr;
    slvr = (struct Solver*)malloc(sizeof(struct Solver));
    struct AppCtx *user;
    user = (struct AppCtx*)malloc(sizeof(struct AppCtx));
    ierr = initialize_petsc(slvr,argc,argv,user);CHKERRQ(ierr);

    printf("\n\n\nGrid size: %dx%d, with %d ions, and %d compartments.\n",Nx,Ny,Ni,Nc);
    PetscLogDouble tic,toc,full_tic,full_toc;
    //Create state_variables struct
    struct SimState *state_vars;
    state_vars=(struct SimState*)malloc(sizeof(struct SimState));
    Vec current_state;
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&current_state);CHKERRQ(ierr);
    ierr = VecSetType(current_state,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(current_state,PETSC_DECIDE,NA);CHKERRQ(ierr);

    struct SimState *state_vars_past;
    state_vars_past=(struct SimState*)malloc(sizeof(struct SimState));
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&state_vars_past->v);CHKERRQ(ierr);
    ierr = VecSetType(state_vars_past->v,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(state_vars_past->v,PETSC_DECIDE,NA);CHKERRQ(ierr);
    //Initialize
    printf("Initialize Data Routines\n");



    //Data struct creation
    ierr = init_simstate(current_state,state_vars); CHKERRQ(ierr);
    ierr = init_simstate(state_vars_past->v, state_vars_past); CHKERRQ(ierr);
    //In order to nicely copy into the past variable we leave this here.
    //Variable initiation
    init(current_state,state_vars);

    ierr = extract_subarray(current_state,state_vars); CHKERRQ(ierr);
    printf("Init Value: c: %f,ph: %f,al: %f\n",state_vars->c[0],state_vars->phi[10],state_vars->alpha[25]);
    ierr = restore_subarray(current_state,state_vars); CHKERRQ(ierr);
    //Create the constant ion channel vars
    struct ConstVars *con_vars;
    con_vars=(struct ConstVars*)malloc(sizeof(struct ConstVars));

    //Create the gating variables
    struct GateType *gate_vars;
    gate_vars = (struct GateType*) malloc(sizeof(struct GateType));
    //Create the flux structure
    struct FluxData *flux;
    flux = (struct FluxData*) malloc(sizeof(struct FluxData));
    //Create Excitation
    struct ExctType *gexct;
    gexct = (struct ExctType*)malloc(sizeof(struct ExctType));

    //Set the constant variables
    set_params(current_state,state_vars,con_vars,gate_vars,flux);

    //Pass data structs over to AppCtx

    user->slvr = slvr;
    user->con_vars = con_vars;
    user->gate_vars=gate_vars;
    user->flux=flux;
    user->gexct=gexct;
    user->state_vars_past=state_vars_past;
    user->state_vars=state_vars;

    printf("Steady State Routine\n");

	//Run Initialization routine to get to steady state
	initialize_data(current_state,user);

    printf("Beginning Main Routine \n");
    printf("\n\n\n");
    //Open file to write to
    FILE *fp;
    fp = fopen("data_csd.txt","w");

    FILE *fptime;
    fptime = fopen("timing.txt","a");
    extract_subarray(current_state,state_vars);
    write_data(fp,state_vars,1);

    //Reset time step
    user->dt = dt;
    //Create the excitation
    excitation(user->gexct,0);
    int count = 0;
    PetscInt num_iters;
    PetscTime(&full_tic);
    for(PetscReal t=dt;t<=Time;t+=dt)
    {
        //Save the "current" aka past state
        ierr = restore_subarray(user->state_vars_past->v,user->state_vars_past); CHKERRQ(ierr);
        ierr = copy_simstate(current_state,user->state_vars_past); CHKERRQ(ierr);
        //Update diffusion with past
        //compute diffusion coefficients
        diff_coef(user->Dcs,state_vars->alpha,1);
        //Bath diffusion
        diff_coef(user->Dcb,state_vars->alpha,Batheps);
        restore_subarray(current_state,state_vars);
        //Newton update
        PetscTime(&tic);
//        newton_solve(current_state,user);
        SNESSolve(user->slvr->snes,NULL,current_state);
        PetscTime(&toc);


        if(details) {
            SNESGetIterationNumber(user->slvr->snes,&num_iters);

            printf("Newton time: %f,iters:%d\n", toc - tic,num_iters);
        }
        //Update gating variables
        extract_subarray(current_state,user->state_vars);
        gatevars_update(user->gate_vars,state_vars,user->dt*1e3,0);
        //Update Excitation
        excitation(user->gexct,t);
        count++;
        if(count%krecordfreq==0) {
            printf("Time: %f, Netwon SolveTime %f\n",t,toc-tic);
            write_data(fp, state_vars, 0);
        }

    }
    PetscTime(&full_toc);
    //Close
    fclose(fp);
    fprintf(fptime,"%d,%f\n",count,full_toc-full_tic);
    fclose(fptime);
    printf("Finished Running. Full solve time: %.10e\n",full_toc-full_tic);
    //Free memory
    free(state_vars);free(con_vars);free(gate_vars);
    VecDestroy(&slvr->Q); VecDestroy(&slvr->Res); MatDestroy(&slvr->A);
    KSPDestroy(&slvr->ksp);
    free(slvr);
    return 0;
}

