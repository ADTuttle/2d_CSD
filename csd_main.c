
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
    PetscReal dt = user->dt;
    PetscInt Nt = (PetscInt) floor(Time/dt);
    PetscInt numrecords = (PetscInt)floor(Time/trecordstep);
    PetscInt krecordfreq = (PetscInt)floor(trecordstep/dt);

    if(Profiling_on) {
        PetscLogStage stage1;
        PetscLogStageRegister("Initialize", &stage1);
        //Start events
        init_events(user);
        PetscLogStagePush(stage1);
    }

    printf("\n\n\nGrid size: %dx%d, with %d ions, and %d compartments. For %f sec at step %f\n",user->Nx,user->Ny,Ni,Nc,Time,dt);
    PetscLogDouble tic,toc,full_tic,full_toc;
    //Create state_variables struct
    struct SimState *state_vars;
    state_vars=(struct SimState*)malloc(sizeof(struct SimState));
    Vec current_state;
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&current_state);CHKERRQ(ierr);
    ierr = VecSetType(current_state,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(current_state,PETSC_DECIDE,user->NA);CHKERRQ(ierr);

    struct SimState *state_vars_past;
    state_vars_past=(struct SimState*)malloc(sizeof(struct SimState));
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&state_vars_past->v);CHKERRQ(ierr);
    ierr = VecSetType(state_vars_past->v,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(state_vars_past->v,PETSC_DECIDE,user->NA);CHKERRQ(ierr);
    //Initialize
    printf("Initialize Data Routines\n");



    //Data struct creation
    ierr = init_simstate(current_state,state_vars,user); CHKERRQ(ierr);
    ierr = init_simstate(state_vars_past->v, state_vars_past,user); CHKERRQ(ierr);
    //In order to nicely copy into the past variable we leave this here.
    //Variable initiation
    init(current_state,state_vars,user);

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
    struct BathType *bath;
    bath = (struct BathType*)malloc(sizeof(struct BathType));
    //Pass data structs over to AppCtx

    user->slvr = slvr;
    user->con_vars = con_vars;
    user->gate_vars=gate_vars;
    user->flux=flux;
    user->gexct=gexct;
    user->state_vars_past=state_vars_past;
    user->state_vars=state_vars;
    user->bath = bath;

    //Init misc. array sizes
    init_arrays(user);

    //Set the constant variables
    set_params(current_state,state_vars,con_vars,gate_vars,flux,user);



    printf("Steady State Routine\n");

	//Run Initialization routine to get to steady state
	initialize_data(current_state,user);

    if(Profiling_on) {
        PetscLogStage stage2;
        PetscLogStageRegister("Main Loop", &stage2);
        PetscLogStagePop();
        PetscLogStagePush(stage2);
        init_events(user);
    }
    printf("Beginning Main Routine \n");
    printf("\n\n\n");
    //Open file to write to
    FILE *fp;
    fp = fopen("data_csd.txt","w");
    FILE *fptime;
    fptime = fopen("timing.txt","a");
    extract_subarray(current_state,state_vars);
    write_data(fp,user,numrecords,1);
//    write_point(fp,user,numrecords,1);
    //Reset time step
    user->dt = dt;
    //Create the excitation
    excitation(user,texct+1);
    bath_excitation(user,0);
    int count = 0;
    PetscInt num_iters;
    SNESConvergedReason reason;
    PetscTime(&full_tic);
    for(PetscReal t=dt;t<=Time;t+=dt)
    {
        //Save the "current" aka past state
        ierr = restore_subarray(user->state_vars_past->v,user->state_vars_past); CHKERRQ(ierr);
        ierr = copy_simstate(current_state,user->state_vars_past); CHKERRQ(ierr);
        if(separate_vol) {
            //Update volume(uses past c values for wflow)
            volume_update(user->state_vars, user->state_vars_past, user);
        }
        //Update diffusion with past
        //compute diffusion coefficients
        diff_coef(user->Dcs,state_vars_past->alpha,1,user);
        //Bath diffusion
        diff_coef(user->Dcb,state_vars_past->alpha,Batheps,user);
        restore_subarray(current_state,state_vars);
        //Newton update
        PetscTime(&tic);
//        newton_solve(current_state,user);
        SNESSolve(user->slvr->snes,NULL,current_state);
        PetscTime(&toc);


        if(details) {
            SNESGetIterationNumber(user->slvr->snes,&num_iters);
            SNESGetConvergedReason(user->slvr->snes,&reason);
            printf("Newton time: %f,iters:%d, Reason: %d\n", toc - tic,num_iters,reason);
        }
        //Update gating variables
        extract_subarray(current_state,user->state_vars);

        gatevars_update(user->gate_vars,user->state_vars,user->dt*1e3,user,0);
        if(separate_vol) {
            //Update volume (this uses new c values for wflow)
//            volume_update(user->state_vars, user->state_vars_past, user);
        }
        //Update Excitation
//        excitation(user,t);
        bath_excitation(user,t);
        count++;
        if(count%krecordfreq==0) {
            SNESGetIterationNumber(user->slvr->snes,&num_iters);
            SNESGetConvergedReason(user->slvr->snes,&reason);
            printf("Time: %f,Newton time: %f,iters:%d, Reason: %d\n",t, toc - tic,num_iters,reason);
//            write_point(fp, user,numrecords, 0);
            write_data(fp, user,numrecords, 0);
        }
        SNESGetConvergedReason(user->slvr->snes,&reason);
        if(reason<0){
            // Failure Close
            PetscTime(&full_toc);
            fclose(fp);
            fprintf(fptime,"%d,%d,%d,%d,%f,%f\n",0,count,user->Nx,user->Ny,user->dt,full_toc-full_tic);
            fclose(fptime);
            fprintf(stderr, "Netwon Iteration did not converge! Stopping at %f...\n",t);
            exit(EXIT_FAILURE); /* indicate failure.*/}

    }
    PetscTime(&full_toc);
    //Close
    fclose(fp);
    fprintf(fptime,"%d,%d,%d,%d,%f,%f\n",1,count,user->Nx,user->Ny,user->dt,full_toc-full_tic);
    fclose(fptime);
    printf("Finished Running. Full solve time: %.10e\n",full_toc-full_tic);

    if(Profiling_on) {
        PetscLogStagePop();
        PetscLogView(PETSC_VIEWER_STDOUT_SELF);
    }
    //Free memory
    VecDestroy(&current_state); VecDestroy(&state_vars_past->v);
    free(state_vars);free(con_vars);free(gate_vars);
    VecDestroy(&slvr->Q); VecDestroy(&slvr->Res); MatDestroy(&slvr->A);
    KSPDestroy(&slvr->ksp);
    free(slvr);
    PetscFinalize();
    return 0;
}

