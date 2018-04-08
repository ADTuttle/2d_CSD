
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


int main(int argc, char **argv)
{

    PetscErrorCode ierr;
    //Petsc Initialize
    struct Solver *slvr = (struct Solver*)malloc(sizeof(struct Solver));
    struct Solver *grid_slvr = (struct Solver*)malloc(sizeof(struct Solver));
    struct AppCtx *user = (struct AppCtx*)malloc(sizeof(struct AppCtx));
    user->slvr = slvr;
    user->grid_slvr = grid_slvr;
    ierr = initialize_petsc(slvr,argc,argv,user);CHKERRQ(ierr);
    if(Predictor) {
        ierr = initialize_grid_slvr(grid_slvr, argc, argv, user);CHKERRQ(ierr);
    }

    PetscReal dt = user->dt;
    PetscInt Nt = (PetscInt) floor(Time/dt);
    PetscInt numrecords = (PetscInt)floor(Time/trecordstep);
    PetscInt krecordfreq = (PetscInt)floor(trecordstep/dt);
    PetscInt x,y,comp,ion;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscReal dx = Lx/Nx;
    PetscReal dy = Ly/Ny;

    if(Profiling_on) {
        PetscLogStage stage1;
        PetscLogStageRegister("Initialize", &stage1);
        //Start events
        init_events(user);
        PetscLogStagePush(stage1);
    }

    printf("\n\n\nGrid size: %dx%d, with %d ions, and %d compartments. For %f sec at step %f\n",user->Nx,user->Ny,Ni,Nc,Time,dt);
    PetscLogDouble tic,toc,full_tic,full_toc,grid_tic,grid_toc;
    //Create state_variables struct
    struct SimState *state_vars = (struct SimState*)malloc(sizeof(struct SimState));
    Vec current_state;
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&current_state);CHKERRQ(ierr);
    ierr = VecSetType(current_state,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(current_state,PETSC_DECIDE,user->Nx*user->Ny*Nv);CHKERRQ(ierr);

    struct SimState *state_vars_past = (struct SimState*)malloc(sizeof(struct SimState));
    //Create Vector
    ierr = VecCreate(PETSC_COMM_WORLD,&state_vars_past->v);CHKERRQ(ierr);
    ierr = VecSetType(state_vars_past->v,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(state_vars_past->v,PETSC_DECIDE,user->Nx*user->Ny*Nv);CHKERRQ(ierr);

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
    struct ConstVars *con_vars = (struct ConstVars*)malloc(sizeof(struct ConstVars));

    //Create the gating variables
    struct GateType *gate_vars = (struct GateType*) malloc(sizeof(struct GateType));
    struct GateType *gate_vars_past = (struct GateType*) malloc(sizeof(struct GateType));
    //Create the flux structure
    struct FluxData *flux = (struct FluxData*) malloc(sizeof(struct FluxData));
    //Create Excitation
    struct ExctType *gexct = (struct ExctType*)malloc(sizeof(struct ExctType));

    //Create small grid variables
    struct SimState *grid_vars = (struct SimState*)malloc(sizeof(struct SimState));
    struct SimState *grid_vars_past = (struct SimState*)malloc(sizeof(struct SimState));
    struct GateType *grid_gate_vars = (struct GateType*) malloc(sizeof(struct GateType));

    //Pass data structs over to AppCtx
    user->con_vars = con_vars;
    user->gate_vars=gate_vars;
    user->gate_vars_past=gate_vars_past;
    user->flux=flux;
    user->gexct=gexct;
    user->state_vars_past=state_vars_past;
    user->state_vars=state_vars;
    user->grid_gate_vars = grid_gate_vars;
    user->grid_vars_past = grid_vars_past;
    user->grid_vars = grid_vars;

    //Init misc. array sizes
    init_arrays(user);
    parameter_dependence(user);

    //Set the constant variables
    set_params(current_state,state_vars,con_vars,gate_vars,flux,user);



    printf("Steady State Routine\n");

    //Run Initialization routine to get to steady state
    initialize_data(current_state,user);

    /*
    for(x=0;x<user->Nx/2;x++){
        for(y=4;y<6;y++){
            user->con_vars->pNaP[xy_index(x,y,user->Nx)]=0;
            user->con_vars->pKDR[xy_index(x,y,user->Nx)]=0;
            user->con_vars->pKA[xy_index(x,y,user->Nx)]=0;
        }
    }
     */
    PetscReal rad1,rad2;
    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){
//            rad1 = sqrt(pow((x+0.5)*dx-Lx/2,2)+pow((y+0.5)*dy,2));
//            rad2 = sqrt(pow((x+0.5)*dx-Lx,2)+pow((y+0.5)*dy,2));
            rad1 = sqrt(pow((x+0.5)*dx-Lx/2,2)+pow((y+0.5)*dy-Ly/2,2));
            if(rad1<Lx/4 || (x>3*Nx/8 && x<5*Nx/8 && y<Nx/4)){
//            if(rad1<3*Lx/4){
//            if(rad2<3*Lx/4&&rad1>=Lx/4){
                user->con_vars->pNaP[xy_index(x, y, user->Nx)] = 0;
                user->con_vars->pKDR[xy_index(x, y, user->Nx)] = 0;
                user->con_vars->pKA[xy_index(x, y, user->Nx)] = 0;
                if(y==Ny-1){
                    printf("x\n");
                }else{
                    printf("x");
                }
            } else {
                user->con_vars->pNaP[xy_index(x, y, user->Nx)] = basepNaP;
                user->con_vars->pKDR[xy_index(x, y, user->Nx)] = basepKDR;
                user->con_vars->pKA[xy_index(x, y, user->Nx)] = basepKA;
                if(y==Ny-1){
                    printf(".\n");
                }else{
                    printf(".");
                }
            }
        }
    }


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

    FILE *fdt;
    if(Predictor) {
        fdt = fopen("csd_dt.txt", "w");
        save_timestep(fdt,user,numrecords-1,1);
    }

    FILE *fptime;
    fptime = fopen("timing.txt","a");
    extract_subarray(current_state,state_vars);
    write_data(fp,user,numrecords,1);
//    write_point(fp,user,numrecords,1);
    FILE *fpflux;
    fpflux = fopen("flux_csd.txt","w");
    measure_flux(fpflux,user,numrecords,1);

    //Reset time step
    user->dt = dt;
    int count = 0;
    PetscInt num_iters,ksp_iters_old,ksp_iters_new,grid_ksp_old;
    PetscInt total_newton = 0;
    int refinement;
    SNESConvergedReason reason;
    PetscTime(&full_tic);
    for(PetscReal t=dt;t<=Time;t+=dt) {
        //This makes a spiral
        /*
        if(t>100.00 && t<101.00) {
            for (x = 0; x < Nx; x++) {
                for (y = 0; y < Ny; y++) {
                    rad1 = sqrt(pow((x + 0.5) * dx - Lx / 2, 2) + pow((y + 0.5) * dy - Ly / 2, 2));
                    if (rad1 < Lx / 4 || (x > 3 * Nx / 8 && x < 5 * Nx / 8 && y < Nx / 4)) {
                        user->con_vars->pNaP[xy_index(x, y, user->Nx)] = basepNaP;
                        user->con_vars->pKDR[xy_index(x, y, user->Nx)] = basepKDR;
                        user->con_vars->pKA[xy_index(x, y, user->Nx)] = basepKA;
                    }
                }
            }
        }*/
        //This makes a wave along a circle
        if(t>100.00 && t<101.00) {
            for (x = 0; x < Nx; x++) {
                for (y = 0; y < Ny; y++) {
                    rad1 = sqrt(pow((x + 0.5) * dx - Lx / 2, 2) + pow((y + 0.5) * dy - Ly / 2, 2));
                    if (rad1 < Lx / 4) {
                        user->con_vars->pNaP[xy_index(x, y, user->Nx)] = 0;
                        user->con_vars->pKDR[xy_index(x, y, user->Nx)] = 0;
                        user->con_vars->pKA[xy_index(x, y, user->Nx)] = 0;
                    } else{
                        user->con_vars->pNaP[xy_index(x, y, user->Nx)] = basepNaP;
                        user->con_vars->pKDR[xy_index(x, y, user->Nx)] = basepKDR;
                        user->con_vars->pKA[xy_index(x, y, user->Nx)] = basepKA;
                    }
                }
            }
        }
        count++;
        //Save the "current" aka past state
        ierr = restore_subarray(user->state_vars_past->v, user->state_vars_past);CHKERRQ(ierr);
        ierr = copy_simstate(current_state, user->state_vars_past);CHKERRQ(ierr);
        if (separate_vol) {
            memcpy(user->state_vars_past->alpha, user->state_vars->alpha,sizeof(PetscReal) * user->Nx * user->Ny * (Nc - 1));
        }

        //Predict if chosen
        if(Predictor) {
            PetscTime(&grid_tic);
            Update_Solution(current_state, t, user);
            PetscTime(&grid_toc);

            if(count%krecordfreq==0) {
                refinement=0;
                for(int z=0;z<user->Nx*user->Ny;z++){
                    if(user->dt_space[z]<user->dt){refinement++;}
                }
                KSPGetTotalIterations(user->grid_slvr->ksp,&ksp_iters_new);
                printf("Grid Time: %f, Refined %d, AvgKspIters: %.2f\n", grid_toc - grid_tic,refinement,((double)ksp_iters_new-grid_ksp_old)/(user->Nx*user->Ny));
                save_timestep(fdt,user,numrecords,0);
                grid_ksp_old = ksp_iters_new;
            }
        }
        if(separate_vol) {
            //Update volume(uses past c values for wflow)
            volume_update(user->state_vars, user->state_vars_past, user);
        }
        //Update diffusion with past
//        compute diffusion coefficients
        diff_coef(user->Dcs,state_vars_past->alpha,1,user);
//        Bath diffusion
        diff_coef(user->Dcb,state_vars_past->alpha,Batheps,user);
        restore_subarray(current_state,state_vars);

        //Update Excitation
        excitation(user,t-dt);
        PetscTime(&tic);
        SNESSolve(user->slvr->snes,NULL,current_state);
        PetscTime(&toc);


        SNESGetIterationNumber(user->slvr->snes,&num_iters);
        SNESGetConvergedReason(user->slvr->snes,&reason);
        KSPGetTotalIterations(user->slvr->ksp,&ksp_iters_new);
        total_newton+=num_iters;
        if(details) {
            printf("Newton time: %f,SNesiters:%d, Reason: %d, KSPIters: %d\n", toc - tic,num_iters,reason,ksp_iters_new-ksp_iters_old);

        }


        //Update gating variables
        extract_subarray(current_state,user->state_vars);

        gatevars_update(user->gate_vars,user->gate_vars_past,user->state_vars,user->dt*1e3,user,0);

        //Copy old gating variables
        //Save the gating variables
        memcpy(user->gate_vars_past->mNaT,user->gate_vars->mNaT,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->hNaT,user->gate_vars->hNaT,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->gNaT,user->gate_vars->gNaT,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->mNaP,user->gate_vars->mNaP,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->hNaP,user->gate_vars->hNaP,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->gNaP,user->gate_vars->gNaP,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->gKA,user->gate_vars->gKA,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->hKA,user->gate_vars->hKA,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->mKA,user->gate_vars->mKA,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->mKDR,user->gate_vars->mKDR,sizeof(PetscReal)*user->Nx*user->Ny);
        memcpy(user->gate_vars_past->gKDR,user->gate_vars->gKDR,sizeof(PetscReal)*user->Nx*user->Ny);
        //Update the past membrane voltage
        if(Predictor){
            for(x=0;x<user->Nx;x++){
                for(y=0;y<user->Ny;y++){
                    user->vm_past[xy_index(x,y,user->Nx)]=(user->state_vars_past->phi[phi_index(x,y,0,user->Nx)]-user->state_vars_past->phi[phi_index(x,y,Nc-1,user->Nx)])*RTFC;
                }
            }
        }


        if(count%krecordfreq==0) {
            printf("Time: %f,Newton time: %f,iters:%d, Reason: %d,KSPIters: %d\n",t,toc - tic,num_iters,reason,ksp_iters_new-ksp_iters_old);
//            write_point(fp, user,numrecords, 0);
            write_data(fp, user,numrecords, 0);
            measure_flux(fpflux,user,numrecords,0);
            if(count%1000){
                fclose(fp);
                fp = fopen("data_csd.txt","a");
            }
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
    fprintf(fptime,"%d,%d,%d,%d,%f,%f\n",1,count,user->Nx,user->Ny,user->dt,full_toc-full_tic);
    fclose(fptime);
    fclose(fpflux);
    printf("Finished Running. Full solve time: %.10e\n",full_toc-full_tic);
    printf("Total newton iterations:%d\n",total_newton);

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

