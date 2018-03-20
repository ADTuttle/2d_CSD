
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char **argv)
{

    PetscErrorCode ierr;
    //Petsc Initialize
    struct Solver *slvr = (struct Solver*)malloc(sizeof(struct Solver));
    struct AppCtx *user = (struct AppCtx*)malloc(sizeof(struct AppCtx));
    ierr = initialize_petsc(slvr,argc,argv,user);CHKERRQ(ierr);
    user->fast_slvr = (struct Solver*)malloc(sizeof(struct Solver));
    ierr = initialize_fast_petsc(user->fast_slvr,argc,argv,user); CHKERRQ(ierr);
    PetscReal dt = user->dt;
    user->dtf = dt/Nfast;
    PetscInt Nt = (PetscInt) floor(Time/dt);
    PetscInt numrecords = (PetscInt)floor(Time/trecordstep);
    PetscInt krecordfreq = (PetscInt)floor(trecordstep/dt);

//    PetscReal atol;
//    atol=1e-19;
//    SNESSetTolerances(user->fast_slvr->snes,atol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

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
    struct SimState *state_vars = (struct SimState*)malloc(sizeof(struct SimState));
    Vec current_state;
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&current_state);CHKERRQ(ierr);
    ierr = VecSetType(current_state,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(current_state,PETSC_DECIDE,user->NA);CHKERRQ(ierr);

    struct SimState *state_vars_past = (struct SimState*)malloc(sizeof(struct SimState));
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&state_vars_past->v);CHKERRQ(ierr);
    ierr = VecSetType(state_vars_past->v,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(state_vars_past->v,PETSC_DECIDE,user->NA);CHKERRQ(ierr);

    //Init fast variables
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&state_vars->v_fast);CHKERRQ(ierr);
    ierr = VecSetType(state_vars->v_fast,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(state_vars->v_fast,PETSC_DECIDE,user->Nx*user->Ny*Nc*(Ni+1));CHKERRQ(ierr);
    ierr = VecDuplicate(state_vars->v_fast,&state_vars_past->v_fast); CHKERRQ(ierr);

    user->state_vars_past=state_vars_past;
    user->state_vars=state_vars;
    user->state_vars->v = current_state;
    //Initialize
    printf("Initialize Data Routines\n");



    //Data struct creation
    ierr = init_simstate(current_state,state_vars,user); CHKERRQ(ierr);
    ierr = init_simstate(state_vars_past->v, state_vars_past,user); CHKERRQ(ierr);
    //In order to nicely copy into the past variable we leave this here.
    //Variable initiation
    init(current_state,state_vars,user);

    ierr = extract_subarray(current_state,state_vars,1); CHKERRQ(ierr);
    printf("Init Value: c: %f,ph: %f,al: %f\n",state_vars->c[0],state_vars->phi[10],state_vars->alpha[25]);
    ierr = restore_subarray(current_state,state_vars,1); CHKERRQ(ierr);
    //Create the constant ion channel vars
    struct ConstVars *con_vars = (struct ConstVars*)malloc(sizeof(struct ConstVars));

    //Create the gating variables
    struct GateType *gate_vars = (struct GateType*) malloc(sizeof(struct GateType));
    //Create the flux structure
    struct FluxData *flux = (struct FluxData*) malloc(sizeof(struct FluxData));
    //Create Excitation
    struct ExctType *gexct = (struct ExctType*)malloc(sizeof(struct ExctType));
    //Pass data structs over to AppCtx

    user->slvr = slvr;
    user->con_vars = con_vars;
    user->gate_vars=gate_vars;
    user->flux=flux;
    user->gexct=gexct;

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
    extract_subarray(current_state,state_vars,1);
    write_data(fp,user,numrecords,1);

    FILE *ftime;
    ftime = fopen("time_step.txt","w");
    record_fast_timestep(ftime,user,numrecords,1);


//    write_point(fp,user,numrecords,1);
    FILE *fpflux;
    fpflux = fopen("flux_csd.txt","w");
    FILE *fp_fast = fopen("data_fastcsd.txt","w");
    write_fast(fp_fast,user,numrecords*Nfast-(Nfast-1),1);

    measure_flux(fpflux,user,numrecords,1);

    //Reset time step
    user->dt = dt;
    //Create the excitation
    excitation(user,0);
    int count = 0;
    PetscInt num_iters,ksp_iters_old,ksp_iters_new;
    SNESConvergedReason reason;
    PetscTime(&full_tic);
    for(PetscReal t=dt;t<=Time;t+=dt)
    {
        //Save the "current" aka past state
        ierr = restore_subarray(state_vars_past->v,state_vars_past,1); CHKERRQ(ierr);
        ierr = copy_simstate(current_state,state_vars_past,1); CHKERRQ(ierr);

        //Update diffusion with past
        //compute diffusion coefficients
        diff_coef(user->Dcs,state_vars_past->alpha,1,user);
        //Bath diffusion
        diff_coef(user->Dcb,state_vars_past->alpha,Batheps,user);

        //update the fast variables
        update_fast_vars(fp_fast,state_vars,state_vars_past,user,t);
        //record time step
        record_fast_timestep(ftime,user,numrecords,0);

        if(separate_vol) {
            //Update volume(uses past c values for wflow)
            extract_subarray(current_state,state_vars,1);
            extract_subarray(state_vars_past->v,state_vars_past,1);
            volume_update(user->state_vars, user->state_vars_past, user);
            restore_subarray(current_state,state_vars,1);
        }


        //Newton update
        PetscTime(&tic);
        SNESSolve(user->slvr->snes,NULL,current_state);
        PetscTime(&toc);

        SNESGetIterationNumber(user->slvr->snes,&num_iters);
        SNESGetConvergedReason(user->slvr->snes,&reason);
        KSPGetTotalIterations(user->slvr->ksp,&ksp_iters_new);
        if(details) {
            printf("Newton time: %f,SNesiters:%d, Reason: %d, KSPIters: %d\n", toc - tic,num_iters,reason,ksp_iters_new-ksp_iters_old);

        }
        //Reconfigured split, setting phi_fast to zero
        recombine(state_vars,user);

        count++;
        if(count%krecordfreq==0) {
            printf("\nTime: %f,Newton time: %f,iters:%d, Reason: %d,KSPIters: %d\n\n",t, toc - tic,num_iters,reason,ksp_iters_new-ksp_iters_old);
            extract_subarray(state_vars->v,state_vars,1);
//            write_point(fp, user,numrecords, 0);
            write_data(fp, user,numrecords, 0);
            measure_flux(fpflux,user,numrecords,0);
            restore_subarray(state_vars->v,state_vars,1);
        }
        ksp_iters_old = ksp_iters_new;

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
    fclose(fp_fast);
    fclose(ftime);
    fprintf(fptime,"%d,%d,%d,%d,%f,%f\n",1,count,user->Nx,user->Ny,user->dt,full_toc-full_tic);
    fclose(fptime);
    fclose(fpflux);
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
PetscErrorCode update_fast_vars(FILE *fp_fast,struct SimState *state_vars,struct SimState *state_vars_past, struct AppCtx *user,PetscReal t)
{
    PetscErrorCode ierr=0;
    PetscLogDouble tic,toc;

    PetscInt x,y;


        PetscTime(&tic);

        extract_subarray(state_vars_past->v,state_vars_past,1);
        extract_subarray(state_vars->v,state_vars,1);
        for(x=0;x<user->Nx;x++) {
            for (y = 0; y < user->Ny; y++) {

                Point_Solve(state_vars, state_vars_past, user,x,y,t);
            }
        }

        PetscTime(&toc);

        //Update gating variables
//        printf("Nfast: %d,Solve time: %f\n",nfast,toc - tic);
//        printf("phi:%.10e, c: %.10e\n",state_vars->phi_fast[phi_index(user->Nx/2,user->Ny/2,0,user->Nx)],state_vars->c[c_index(user->Nx/2,user->Ny/2,0,0,user->Nx)]);

        write_fast(fp_fast, user,-10, 0);

        restore_subarray(state_vars->v,state_vars,1);
        restore_subarray(state_vars_past->v,state_vars_past,1);


    printf("Fast Solve total: %f\n",toc - tic);

    return ierr;
}
/*
PetscErrorCode update_fast_vars(FILE *fp_fast,struct SimState *state_vars,struct SimState *state_vars_past, struct AppCtx *user,PetscReal t)
{
    PetscErrorCode ierr=0;
    PetscInt num_iters;
    SNESConvergedReason reason;
    PetscLogDouble tic,toc,fulltic,fulltoc;

    PetscTime(&fulltic);
    //Fast time steps
    for(PetscInt nfast=1;nfast<=Nfast;nfast++){

        //Save the "current" aka past state
        ierr = restore_subarray(user->state_vars_past->v_fast,user->state_vars_past,2); CHKERRQ(ierr);
        ierr = copy_simstate(state_vars->v_fast,state_vars_past,2); CHKERRQ(ierr);


        PetscTime(&tic);
//        SNESSolve(user->fast_slvr->snes,NULL,state_vars->v_fast);
        Point_Solve(state_vars,state_vars_past,user);
        PetscTime(&toc);
        SNESGetIterationNumber(user->fast_slvr->snes,&num_iters);
        SNESGetConvergedReason(user->fast_slvr->snes,&reason);

        recombine(state_vars,user);
        //Update gating variables
        extract_subarray(state_vars->v,state_vars,1);
//        printf("Nfast: %d,Solve time: %f\n",nfast,toc - tic);
//        printf("phi:%.10e, c: %.10e\n",state_vars->phi_fast[phi_index(user->Nx/2,user->Ny/2,0,user->Nx)],state_vars->c[c_index(user->Nx/2,user->Ny/2,0,0,user->Nx)]);

        gatevars_update(user->gate_vars,user->state_vars,user->dtf*1e3,user,0);

        write_fast(fp_fast, user,-10, 0);

        restore_subarray(state_vars->v,state_vars,1);

        //update the excitation
        excitation(user,t-user->dt+user->dtf*nfast);

        recombine(state_vars,user);

        if(reason<0){
            extract_subarray(state_vars->v,state_vars,1);
            printf("Fast Newton time: %f,iters:%d, Reason: %d\n",toc - tic,num_iters,reason);
            printf("phi:%.10e, c: %.10e\n",state_vars->phi_fast[phi_index(user->Nx/2,user->Ny/2,0,user->Nx)],state_vars->c[c_index(user->Nx/2,user->Ny/2,0,0,user->Nx)]);

            // Failure Close
            printf("Fast Netwon Solve did not converge! Stopping at %f...\n",t-user->dt+user->dtf*nfast);
            fprintf(stderr, "Fast Netwon Solve did not converge! Stopping at %f...\n",t-user->dt+user->dtf*nfast);
            exit(EXIT_FAILURE); }

    }
    PetscTime(&fulltoc);
    printf("Fast Solve total: %f,Newton time: %f,iters:%d, Reason: %d\n",fulltoc-fulltic,toc - tic,num_iters,reason);

    return ierr;
}
*/