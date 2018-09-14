#include "constants.h"
#include "functions.h"


PetscErrorCode newton_solve(Vec current_state,struct Solver *slvr,struct AppCtx *user)
{

    PetscReal rsd;
    PetscErrorCode ierr = 0;
    PetscReal const *temp;
    PetscInt num_iter,comp,ion,x,y;
    PetscReal rnorm;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;


    PetscLogDouble tic,toc;

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at even positions (0,2,4,...)
    //y at odd (1,3,5,...)
    //still use c_index(x,y,comp,ion,Nx), but with ind*2 or ind*2+1

    extract_subarray(current_state,user->state_vars);

    PetscReal tol = reltol*array_max(user->state_vars->c,(size_t)Nx*Ny*Ni*Nc);
    restore_subarray(current_state,user->state_vars);
    rsd = tol+1;

    for(PetscInt iter=0;iter<1;iter++)
    {
        if(separate_vol){
            if(use_en_deriv){
                ierr = calc_residual_no_vol(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);

            } else{
                ierr = calc_residual_algebraic_no_vol(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);
            }

        }else{
            if(use_en_deriv){
                ierr = calc_residual(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);

            } else{
                ierr = calc_residual_algebraic(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);
            }

        }

        ierr = VecNorm(slvr->Res,NORM_MAX,&rsd);CHKERRQ(ierr);
//        printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
        if(rsd<tol)
        {
            if(details)
            {
                printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
            }
            return ierr;
        }
        if(separate_vol){
            if(use_en_deriv){
                ierr = calc_jacobian_no_vol(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }else{
                ierr = calc_jacobian_algebraic_no_vol(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }

        }else{
            if(use_en_deriv){
                ierr = calc_jacobian(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }else{
                ierr = calc_jacobian_algebraic(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }
        }
        //Set the new operator
        ierr = KSPSetOperators(slvr->ksp,slvr->A,slvr->A);CHKERRQ(ierr);

        //Solve
        PetscTime(&tic);
        ierr = KSPSolve(slvr->ksp,slvr->Res,slvr->Q);CHKERRQ(ierr);
        PetscTime(&toc);

        ierr = KSPGetIterationNumber(user->slvr->ksp,&num_iter); CHKERRQ(ierr);
        ierr =  KSPGetResidualNorm(user->slvr->ksp,&rnorm); CHKERRQ(ierr);


        if(details) {
            printf("KSP Solve time: %f, iter num:%d, norm: %.10e\n",toc-tic,num_iter,rnorm);
        }

//        PetscTime(&tic);
        ierr = VecGetArrayRead(slvr->Q,&temp); CHKERRQ(ierr);
        extract_subarray(current_state,user->state_vars);
        for(x=0;x<Nx;x++){
            for(y=0;y<Ny;y++){
                for(comp=0;comp<Nc;comp++){
                    for(ion=0;ion<Ni;ion++){
                        user->state_vars->c[c_index(x, y, 0, comp, ion, Nx, 0)]-=temp[Ind_1(x, y, 0, ion, comp, Nx, 0)];
                    }
                    user->state_vars->phi[phi_index(x, y, 0, comp, Nx, 0)]-=temp[Ind_1(x, y, 0, Ni, comp, Nx, 0)];
                }
                if(!separate_vol){
                    for(comp=0;comp<Nc-1;comp++) {
                        user->state_vars->alpha[al_index(x, y, 0, comp, Nx, 0)] -= temp[Ind_1(x, y, 0, Ni + 1, comp, Nx,
                                                                                              0)];
                    }
                }
            }
        }

        ierr = VecRestoreArrayRead(slvr->Q,&temp);
        restore_subarray(current_state,user->state_vars);

        if(details)
        {
            printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
        }
    }
/*
    if(rsd>tol)
    {
        fprintf(stderr, "Netwon Iteration did not converge! Stopping...\n");
        exit(EXIT_FAILURE);
    }
    */
    return ierr;
}


PetscErrorCode calc_residual(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    //Volume is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);

    //Compute membrane water flow related quantities
    wflowm(user);


    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc], RphxRight[Nc], RphyUp[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            //Init voltage tracking to zero
            for(comp=0;comp<Nc;comp++) {
                Rphx[comp]=0;
                Rphy[comp]=0;
                RphxRight[comp]=0;
                RphyUp[comp]=0;
            }
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0)]) / 2;
                        Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x + 1, y,
                                                                                                               0, comp,
                                                                                                               ion, Nx,
                                                                                                               0)]) / 2;
                        RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                                comp,
                                                                                                                ion, Nx,
                                                                                                                0)]) + z_charge[ion] * (phi[phi_index(
                                x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                                 Nx, 0)] + cp[c_index(x,
                                                                                                                      y,
                                                                                                                      0,
                                                                                                                      comp,
                                                                                                                      ion,
                                                                                                                      Nx,
                                                                                                                      0)]) / 2;
                        Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                               0)] + cp[c_index(x,
                                                                                                                y + 1, 0,
                                                                                                                comp,
                                                                                                                ion,
                                                                                                                Nx,
                                                                                                                0)]) / 2;
                        RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) + z_charge[ion] * (phi[phi_index(
                                x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    Resc = al[al_index(x, y, 0, comp, Nx, 0)] * c[c_index(x, y, 0, comp, ion, Nx, 0)] - alp[al_index(x,
                                                                                                                     y,
                                                                                                                     0,
                                                                                                                     comp,
                                                                                                                     Nx,
                                                                                                                     0)] * cp[c_index(
                            x, y, 0, comp, ion, Nx, 0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;

                    ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp]+=z_charge[ion]*Rcvx;
                    Rphy[comp]+=z_charge[ion]*Rcvy;
                    RphxRight[comp]+=z_charge[ion]*RcvxRight;
                    RphyUp[comp]+=z_charge[ion]*RcvyUp;

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)];
                alpNc = 1 - alp[al_index(x, y, 0, 0, Nx, 0)] - alp[al_index(x, y, 0, 1, Nx, 0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                         0)] + cp[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) / 2;
                    Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                            comp, ion,
                                                                                                            Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp,
                                                                                                              ion, Nx,
                                                                                                              0)]) / 2;
                    Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                      ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                Resc = alNc*c[c_index(x, y, 0, comp, ion, Nx, 0)] - alpNc * cp[c_index(x, y, 0, comp, ion, Nx, 0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, comp, ion,
                                                                                                   Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
                ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z_charge[ion]*Rcvx;
                Rphy[comp]+=z_charge[ion]*Rcvy;
                RphxRight[comp]+=z_charge[ion]*RcvxRight;
                RphyUp[comp]+=z_charge[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[comp]*(phi[phi_index(x, y, 0, comp, Nx, 0)] - phi[phi_index(x, y, 0, Nc - 1, Nx, 0)]) - cm[comp] * (phip[phi_index(
                        x,
                        y,
                        0,
                        comp,
                        Nx,
                        0)] - phip[phi_index(
                        x, y, 0, Nc - 1, Nx, 0)]);
                for(ion=0;ion<Ni;ion++){
                    //Ion channel
                    Resph += z_charge[ion] * flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add the terms shared with extracell
                ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                Resph += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resph, INSERT_VALUES); CHKERRQ(ierr);
            }

            //Finish adding extracell
            comp = Nc-1;
            //Add bath contribution
            for(ion=0;ion<Ni;ion++){

                ResphN -= z_charge[ion] * sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), ResphN, INSERT_VALUES); CHKERRQ(ierr);
        }
    }


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            //Residual for water flow
            //Plus modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++) {
                //Water flow
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), al[al_index(x, y, 0, comp, Nx, 0)] - alp[al_index(x, y,
                                                                                                                               0,
                                                                                                                               comp,
                                                                                                                               Nx,
                                                                                                                               0)] + flux->wflow[al_index(
                        x, y, 0, comp, Nx, 0)] * dt, INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Volume is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    PetscReal Fphph0x[Nc],Fphph1x[Nc];
    PetscReal Fphph0y[Nc],Fphph1y[Nc];

    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(comp=0;comp<Nc;comp++){
                Fphph0x[comp]=0;
                Fphph1x[comp]=0;
                Fphph0y[comp]=0;
                Fphph1y[comp]=0;
            }
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions
                    Ftmpx = 0;
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Ftmpy = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if(x<Nx-1) {
                        Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                                x + 1, y, 0,
                                comp, ion, Nx,
                                0)]) / 2 / dx * dt / dx;
                        Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0) {
                        Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dx * dt / dx;
                        Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x, y + 1,
                                                                                                               0,
                                                                                                               comp,
                                                                                                               ion,
                                                                                                               Nx,
                                                                                                               0)]) / 2 / dy * dt / dy;
                        Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)] + cp[c_index(
                                x, y, 0,
                                comp,
                                ion, Nx,
                                0)]) / 2 / dy * dt / dy;
                        Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x, y, 0, comp, Nx, 0)] + Fc0x + Fc1x + Fc0y + Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Add up terms for voltage eqns
                    Fphph0x[comp]+=z_charge[ion]*Fph0x;
                    Fphph1x[comp]+=z_charge[ion]*Fph1x;
                    Fphph0y[comp]+=z_charge[ion]*Fph0y;
                    Fphph1y[comp]+=z_charge[ion]*Fph1y;

                    //membrane current contributions
                    Ac+= flux->dfdci[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi+= flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dfdci[c_index(x, y, 0,
                                                                                                                                       comp,
                                                                                                                                       ion,
                                                                                                                                       Nx,
                                                                                                                                       0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dfdce[c_index(x, y, 0,
                                                                                                                                      comp,
                                                                                                                                      ion,
                                                                                                                                      Nx,
                                                                                                                                      0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume terms
                    //C extra with intra alpha
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), -c[c_index(x, y, 0, Nc - 1,
                                                                                                                                ion, Nx, 0)], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //C intra with intra alpha
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), c[c_index(x, y, 0, comp,
                                                                                                                             ion,
                                                                                                                             Nx, 0)], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * (Fc0x + Fc1x + Fc0y + Fc1y + flux->dfdci[c_index(
                            x, y, 0, comp, ion, Nx, 0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //IntraPhi with c extra(volt eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), z_charge[ion] * (flux->dfdce[c_index(x,
                                                                                                                                               y,
                                                                                                                                               0,
                                                                                                                                               comp,
                                                                                                                                               ion,
                                                                                                                                               Nx,
                                                                                                                                               0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //Extra-Phi with intra-c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * (flux->dfdci[c_index(
                            x,
                            y,
                            0,
                            comp,
                            ion,
                            Nx,
                            0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                Ftmpx = 0;
                Fc0x = 0;
                Fc1x = 0;
                Fph0x = 0;
                Fph1x = 0;
                Ftmpy = 0;
                Fc0y = 0;
                Fc1y = 0;
                Fph0y = 0;
                Fph1y = 0;
                if(x<Nx-1) {
                    Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0, comp,
                            ion, Nx, 0)]) / 2 / dx * dt / dx;
                    Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0x = z_charge[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Right Phi with left c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                          0)] + cp[c_index(x, y, 0,
                                                                                                           comp,
                                                                                                           ion, Nx,
                                                                                                           0)]) / 2 / dx * dt / dx;
                    Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1x = z_charge[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // left Phi with right c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx, 0)]) / 2 / dy * dt / dy;
                    Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0y = z_charge[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Upper Phi with lower c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dy * dt / dy;
                    Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1y = z_charge[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Lower Phi with upper c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) + Fc0x + Fc1x + Fc0y + Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                Avolt = z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y);

                //Add up terms for voltage eqns
                Fphph0x[comp]+=z_charge[ion]*Fph0x;
                Fphph1x[comp]+=z_charge[ion]*Fph1x;
                Fphph0y[comp]+=z_charge[ion]*Fph0y;
                Fphph1y[comp]+=z_charge[ion]*Fph1y;

                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++) {
                    Ac -= flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Avolt -= z_charge[ion] * flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1, ion,
                                                                                                   Nx, 0) * 2 + 1], 2));
                Ac -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y, 0, Nc - 1,
                                                                                                       ion, Nx, 0)]) * dt;
                Aphi -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;

                Avolt -= z_charge[ion] * Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y,
                                                                                                                   0,
                                                                                                                   Nc -
                                                                                                                   1,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) * dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for(comp=0;comp<Nc-1;comp++) {
                if(x<Nx-1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0x[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1x[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0y[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1y[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                Avolt = cm[comp]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];
                AvoltN = -cm[comp];
                for(ion=0;ion<Ni;ion++) {
                    Avolt+= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    AvoltN-= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), AvoltN, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc-1;
            if(x<Nx-1) {
                //Right phi with left phi (-Fph0x)
                ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0x[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(x>0) {
                //Left phi with right phi (-Fph1x)
                ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1x[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y<Ny-1) {
                //Upper phi with lower phi (-Fph0y)
                ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0y[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y>0) {
                //Lower phi with upper phi (-Fph1y)
                ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1y[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            AvoltN = 0;

            for(int k=0;k<Nc-1;k++) {
                AvoltN += cm[k];
                Avolt = -cm[k];
                for(ion=0;ion<Ni;ion++) {
                    Avolt-= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, k, ion, Nx, 0)] * dt;
                    AvoltN+= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, k, ion, Nx, 0)] * dt;
                }
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, k, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }

            AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];

            //Bath terms
            for(ion=0;ion<Ni;ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1,
                                                                                                     ion, Nx, 0) * 2 + 1], 2));
                AvoltN -= z_charge[ion] * Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), AvoltN, INSERT_VALUES);CHKERRQ(ierr);
            ind++;

        }
    }
    //water flow
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(comp=0;comp<Nc-1;comp++) {
                //Water flow volume fraction entries
                //Volume to Volume
                Ac=1+ (flux->dwdpi[al_index(x, y, 0, comp, Nx, 0)] * (con_vars->ao[phi_index(0, 0, 0, Nc - 1, Nx, 0)] / (pow(1 - al[al_index(
                        x, y, 0, 0, Nx, 0)] - al[al_index(
                        x, y, 0, 1, Nx, 0)], 2)) + con_vars->ao[phi_index(
                        0, 0, 0, comp, Nx, 0)] / pow(al[al_index(x, y, 0, comp, Nx, 0)], 2)) + flux->dwdal[al_index(x,
                                                                                                                    y, 0,
                                                                                                                    comp,
                                                                                                                    Nx,
                                                                                                                    0)]) * dt;
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Off diagonal (from aNc=1-sum(ak))
                for (PetscInt l=0; l<comp; l++) {
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, l, Nx, 0), flux->dwdpi[al_index(x, y, 0,
                                                                                                                                        comp,
                                                                                                                                        Nx,
                                                                                                                                        0)] * con_vars->ao[phi_index(
                            0, 0, 0, Nc - 1, Nx, 0)] / pow(1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1,
                                                                                                             Nx,
                                                                                                             0)], 2) * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (PetscInt l=comp+1; l<Nc-1; l++) {
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, l, Nx, 0), flux->dwdpi[al_index(x, y, 0,
                                                                                                                                        comp,
                                                                                                                                        Nx,
                                                                                                                                        0)] * con_vars->ao[phi_index(
                            0, 0, 0, Nc - 1, Nx, 0)] / ((1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1,
                                                                                                           Nx, 0)]) * (1 - al[al_index(
                            x, y, 0, 0, Nx, 0)] - al[al_index(
                            x, y, 0, 1, Nx, 0)])) * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (ion=0; ion<Ni; ion++) {
                    //Volume to extra c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dwdpi[al_index(x, y,
                                                                                                                                          0,
                                                                                                                                          comp,
                                                                                                                                          Nx,
                                                                                                                                          0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume to intra c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dwdpi[al_index(x, y,
                                                                                                                                         0,
                                                                                                                                         comp,
                                                                                                                                         Nx,
                                                                                                                                         0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
            }
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

void volume_update(struct SimState *state_vars,struct SimState *state_vars_past, struct AppCtx *user)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[7], 0, 0, 0, 0);
    }
    int x,y,comp;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    for(int n=0;n<1;n++) {

//        memcpy(state_vars_past->alpha, state_vars->alpha, sizeof(PetscReal) * Nx * Ny * (Nc - 1));
        //Forward Euler update
/*
    wflowm(user->flux,user->state_vars_past,user->con_vars);
    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){
            for(comp=0;comp<Nc-1;comp++) {
                state_vars->alpha[al_index(x, y, comp,Nx)] = state_vars_past->alpha[al_index(x,y,comp,Nx)]+user->dt*user->flux->wflow[al_index(x,y,comp,Nx)];
            }
        }
    }
*/
//    /*
        //Backward Euler update
        PetscReal res, Func, Deriv, max_res;
        for (int iter = 0; iter < 10; iter++) {
            max_res = 0;
            wflowm(user);
            for (x = 0; x < Nx; x++) {
                for (y = 0; y < Ny; y++) {
                    for (comp = 0; comp < Nc - 1; comp++) {

                        Func = state_vars->alpha[al_index(x, y, 0, comp, Nx, 0)] - state_vars_past->alpha[al_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   Nx,
                                                                                                                   0)] +
                               dt * user->flux->wflow[al_index(x, y, 0, comp, Nx, 0)];

                        Deriv = 1 + dt * user->flux->dwdal[al_index(x, y, 0, comp, Nx, 0)];

                        res = -Func / Deriv;

                        state_vars->alpha[al_index(x, y, 0, comp, Nx, 0)] += res;

                        if (fabs(res) > max_res) { max_res = fabs(res); }
                    }
                }
            }
            if (max_res < reltol) {
                if(Profiling_on) {
                    PetscLogEventEnd(event[7], 0, 0, 0, 0);
                }
                return; }
        }
    }
    fprintf(stderr,"Volume failed to update!\n");
    if(Profiling_on) {
        PetscLogEventEnd(event[7], 0, 0, 0, 0);
    }
    exit(EXIT_FAILURE); /* indicate failure.*/
}

PetscErrorCode calc_residual_no_vol(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);

    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc], RphxRight[Nc], RphyUp[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            //Init voltage tracking to zero
            for(comp=0;comp<Nc;comp++) {
                Rphx[comp]=0;
                Rphy[comp]=0;
                RphxRight[comp]=0;
                RphyUp[comp]=0;
            }
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0)]) / 2;
                        Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x + 1, y,
                                                                                                               0, comp,
                                                                                                               ion, Nx,
                                                                                                               0)]) / 2;
                        RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                                comp,
                                                                                                                ion, Nx,
                                                                                                                0)]) + z_charge[ion] * (phi[phi_index(
                                x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                                 Nx, 0)] + cp[c_index(x,
                                                                                                                      y,
                                                                                                                      0,
                                                                                                                      comp,
                                                                                                                      ion,
                                                                                                                      Nx,
                                                                                                                      0)]) / 2;
                        Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                               0)] + cp[c_index(x,
                                                                                                                y + 1, 0,
                                                                                                                comp,
                                                                                                                ion,
                                                                                                                Nx,
                                                                                                                0)]) / 2;
                        RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) + z_charge[ion] * (phi[phi_index(
                                x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    Resc = al[al_index(x, y, 0, comp, Nx, 0)] * c[c_index(x, y, 0, comp, ion, Nx, 0)] - alp[al_index(x,
                                                                                                                     y,
                                                                                                                     0,
                                                                                                                     comp,
                                                                                                                     Nx,
                                                                                                                     0)] * cp[c_index(
                            x, y, 0, comp, ion, Nx, 0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;

                    ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp]+=z_charge[ion]*Rcvx;
                    Rphy[comp]+=z_charge[ion]*Rcvy;
                    RphxRight[comp]+=z_charge[ion]*RcvxRight;
                    RphyUp[comp]+=z_charge[ion]*RcvyUp;

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)];
                alpNc = 1 - alp[al_index(x, y, 0, 0, Nx, 0)] - alp[al_index(x, y, 0, 1, Nx, 0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                         0)] + cp[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) / 2;
                    Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                            comp, ion,
                                                                                                            Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp,
                                                                                                              ion, Nx,
                                                                                                              0)]) / 2;
                    Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                      ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                Resc = alNc*c[c_index(x, y, 0, comp, ion, Nx, 0)] - alpNc * cp[c_index(x, y, 0, comp, ion, Nx, 0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, comp, ion,
                                                                                                   Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
                ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z_charge[ion]*Rcvx;
                Rphy[comp]+=z_charge[ion]*Rcvy;
                RphxRight[comp]+=z_charge[ion]*RcvxRight;
                RphyUp[comp]+=z_charge[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[comp]*(phi[phi_index(x, y, 0, comp, Nx, 0)] - phi[phi_index(x, y, 0, Nc - 1, Nx, 0)]) - cm[comp] * (phip[phi_index(
                        x,
                        y,
                        0,
                        comp,
                        Nx,
                        0)] - phip[phi_index(
                        x, y, 0, Nc - 1, Nx, 0)]);
                for(ion=0;ion<Ni;ion++){
                    //Ion channel
                    Resph += z_charge[ion] * flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add the terms shared with extracell
                ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                Resph += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resph, INSERT_VALUES); CHKERRQ(ierr);
            }

            //Finish adding extracell
            comp = Nc-1;
            //Add bath contribution
            for(ion=0;ion<Ni;ion++){

                ResphN -= z_charge[ion] * sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), ResphN, INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_no_vol(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    PetscReal Fphph0x[Nc],Fphph1x[Nc];
    PetscReal Fphph0y[Nc],Fphph1y[Nc];

    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(comp=0;comp<Nc;comp++){
                Fphph0x[comp]=0;
                Fphph1x[comp]=0;
                Fphph0y[comp]=0;
                Fphph1y[comp]=0;
            }
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions
                    Ftmpx = 0;
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Ftmpy = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if(x<Nx-1) {
                        Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                                x + 1, y, 0,
                                comp, ion, Nx,
                                0)]) / 2 / dx * dt / dx;
                        Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0) {
                        Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dx * dt / dx;
                        Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x, y + 1,
                                                                                                               0,
                                                                                                               comp,
                                                                                                               ion,
                                                                                                               Nx,
                                                                                                               0)]) / 2 / dy * dt / dy;
                        Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)] + cp[c_index(
                                x, y, 0,
                                comp,
                                ion, Nx,
                                0)]) / 2 / dy * dt / dy;
                        Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x, y, 0, comp, Nx, 0)] + Fc0x + Fc1x + Fc0y + Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Add up terms for voltage eqns
                    Fphph0x[comp]+=z_charge[ion]*Fph0x;
                    Fphph1x[comp]+=z_charge[ion]*Fph1x;
                    Fphph0y[comp]+=z_charge[ion]*Fph0y;
                    Fphph1y[comp]+=z_charge[ion]*Fph1y;

                    //membrane current contributions
                    Ac+= flux->dfdci[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi+= flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dfdci[c_index(x, y, 0,
                                                                                                                                       comp,
                                                                                                                                       ion,
                                                                                                                                       Nx,
                                                                                                                                       0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dfdce[c_index(x, y, 0,
                                                                                                                                      comp,
                                                                                                                                      ion,
                                                                                                                                      Nx,
                                                                                                                                      0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * (Fc0x + Fc1x + Fc0y + Fc1y + flux->dfdci[c_index(
                            x, y, 0, comp, ion, Nx, 0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //IntraPhi with c extra(volt eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), z_charge[ion] * (flux->dfdce[c_index(x,
                                                                                                                                               y,
                                                                                                                                               0,
                                                                                                                                               comp,
                                                                                                                                               ion,
                                                                                                                                               Nx,
                                                                                                                                               0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //Extra-Phi with intra-c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * (flux->dfdci[c_index(
                            x,
                            y,
                            0,
                            comp,
                            ion,
                            Nx,
                            0)] * dt), INSERT_VALUES); CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                Ftmpx = 0;
                Fc0x = 0;
                Fc1x = 0;
                Fph0x = 0;
                Fph1x = 0;
                Ftmpy = 0;
                Fc0y = 0;
                Fc1y = 0;
                Fph0y = 0;
                Fph1y = 0;
                if(x<Nx-1) {
                    Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0, comp,
                            ion, Nx, 0)]) / 2 / dx * dt / dx;
                    Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0x = z_charge[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Right Phi with left c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                          0)] + cp[c_index(x, y, 0,
                                                                                                           comp,
                                                                                                           ion, Nx,
                                                                                                           0)]) / 2 / dx * dt / dx;
                    Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1x = z_charge[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // left Phi with right c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx, 0)]) / 2 / dy * dt / dy;
                    Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0y = z_charge[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Upper Phi with lower c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dy * dt / dy;
                    Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1y = z_charge[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    // Lower Phi with upper c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -z_charge[ion] * Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) + Fc0x + Fc1x + Fc0y + Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                Avolt = z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y);

                //Add up terms for voltage eqns
                Fphph0x[comp]+=z_charge[ion]*Fph0x;
                Fphph1x[comp]+=z_charge[ion]*Fph1x;
                Fphph0y[comp]+=z_charge[ion]*Fph0y;
                Fphph1y[comp]+=z_charge[ion]*Fph1y;

                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++) {
                    Ac -= flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Avolt -= z_charge[ion] * flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1, ion,
                                                                                                   Nx, 0) * 2 + 1], 2));
                Ac -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y, 0, Nc - 1,
                                                                                                       ion, Nx, 0)]) * dt;
                Aphi -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;

                Avolt -= z_charge[ion] * Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y,
                                                                                                                   0,
                                                                                                                   Nc -
                                                                                                                   1,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) * dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for(comp=0;comp<Nc-1;comp++) {
                if(x<Nx-1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0x[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1x[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0y[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1y[comp], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                Avolt = cm[comp]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];
                AvoltN = -cm[comp];
                for(ion=0;ion<Ni;ion++) {
                    Avolt+= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    AvoltN-= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), AvoltN, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc-1;
            if(x<Nx-1) {
                //Right phi with left phi (-Fph0x)
                ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0x[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(x>0) {
                //Left phi with right phi (-Fph1x)
                ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1x[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y<Ny-1) {
                //Upper phi with lower phi (-Fph0y)
                ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph0y[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y>0) {
                //Lower phi with upper phi (-Fph1y)
                ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fphph1y[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            AvoltN = 0;

            for(int k=0;k<Nc-1;k++) {
                AvoltN += cm[k];
                Avolt = -cm[k];
                for(ion=0;ion<Ni;ion++) {
                    Avolt-= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, k, ion, Nx, 0)] * dt;
                    AvoltN+= z_charge[ion] * flux->dfdphim[c_index(x, y, 0, k, ion, Nx, 0)] * dt;
                }
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, k, Nx, 0), Avolt, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }

            AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];

            //Bath terms
            for(ion=0;ion<Ni;ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1,
                                                                                                     ion, Nx, 0) * 2 + 1], 2));
                AvoltN -= z_charge[ion] * Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), AvoltN, INSERT_VALUES);CHKERRQ(ierr);
            ind++;

        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode calc_residual_algebraic(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using algebraic version of the charge-capacitance relation
    //Alpha is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);


    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0)]) / 2;
                        Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x + 1, y,
                                                                                                               0, comp,
                                                                                                               ion, Nx,
                                                                                                               0)]) / 2;
                        RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                                comp,
                                                                                                                ion, Nx,
                                                                                                                0)]) + z_charge[ion] * (phi[phi_index(
                                x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                                 Nx, 0)] + cp[c_index(x,
                                                                                                                      y,
                                                                                                                      0,
                                                                                                                      comp,
                                                                                                                      ion,
                                                                                                                      Nx,
                                                                                                                      0)]) / 2;
                        Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                               0)] + cp[c_index(x,
                                                                                                                y + 1, 0,
                                                                                                                comp,
                                                                                                                ion,
                                                                                                                Nx,
                                                                                                                0)]) / 2;
                        RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) + z_charge[ion] * (phi[phi_index(
                                x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    Resc = al[al_index(x, y, 0, comp, Nx, 0)] * c[c_index(x, y, 0, comp, ion, Nx, 0)] - alp[al_index(x,
                                                                                                                     y,
                                                                                                                     0,
                                                                                                                     comp,
                                                                                                                     Nx,
                                                                                                                     0)] * cp[c_index(
                            x, y, 0, comp, ion, Nx, 0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;

                    ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)];
                alpNc = 1 - alp[al_index(x, y, 0, 0, Nx, 0)] - alp[al_index(x, y, 0, 1, Nx, 0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                         0)] + cp[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) / 2;
                    Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                            comp, ion,
                                                                                                            Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp,
                                                                                                              ion, Nx,
                                                                                                              0)]) / 2;
                    Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                      ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                Resc = alNc*c[c_index(x, y, 0, comp, ion, Nx, 0)] - alpNc * cp[c_index(x, y, 0, comp, ion, Nx, 0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, comp, ion,
                                                                                                   Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
                ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {

            //Residual for electroneutrality condition
            for(comp=0;comp<Nc-1;comp++)
            {

                Resc = al[al_index(x, y, 0, comp, Nx, 0)] *cz(c,z_charge,x,y,0,Nx,0,comp,user)+user->con_vars->zo[phi_index(0,0,0,comp,
                                                                                                                            Nx,0)]*user->con_vars->ao[phi_index(
                        0, 0, 0, comp, Nx, 0)];
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resc, INSERT_VALUES); CHKERRQ(ierr);
            }
            //Extracellular term
            comp=Nc-1;
            Resc = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) *
                           cz(c,z_charge,x,y,0,Nx,0,comp,user)+user->con_vars->zo[phi_index(
                    0, 0, 0, comp, Nx, 0)]*user->con_vars->ao[phi_index(
                    0, 0, 0, comp, Nx, 0)];
            ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resc, INSERT_VALUES); CHKERRQ(ierr);

            //Residual for water flow
            //Plus modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++) {
                //Water flow
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), al[al_index(x, y, 0, comp, Nx, 0)] - alp[al_index(x, y,
                                                                                                                               0,
                                                                                                                               comp,
                                                                                                                               Nx,
                                                                                                                               0)] + flux->wflow[al_index(
                        x, y, 0, comp, Nx, 0)] * dt, INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            // Add Modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extracell voltage
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -cm[comp] * (phi[phi_index(x, y, 0, Nc - 1, Nx, 0)] - phi[phi_index(
                        x, y,
                        0,
                        comp,
                        Nx, 0)]), ADD_VALUES);CHKERRQ(ierr);
                //Intracell voltage mod
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), -cm[comp] * (phi[phi_index(x, y, 0, comp, Nx, 0)] - phi[phi_index(
                        x, y,
                        0,
                        Nc -
                        1, Nx,
                        0)]), ADD_VALUES);CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);



    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_algebraic(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using algebraic version of the charge-capacitance relation
    // Alpha is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi;


    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions
                    Ftmpx = 0;
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Ftmpy = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if(x<Nx-1) {
                        Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                                x + 1, y, 0,
                                comp, ion, Nx,
                                0)]) / 2 / dx * dt / dx;
                        Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                    }
                    if(x>0) {
                        Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dx * dt / dx;
                        Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x, y + 1,
                                                                                                               0,
                                                                                                               comp,
                                                                                                               ion,
                                                                                                               Nx,
                                                                                                               0)]) / 2 / dy * dt / dy;
                        Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)] + cp[c_index(
                                x, y, 0,
                                comp,
                                ion, Nx,
                                0)]) / 2 / dy * dt / dy;
                        Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x, y, 0, comp, Nx, 0)] + Fc0x + Fc1x + Fc0y + Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;


                    //membrane current contributions
                    Ac+= flux->dfdci[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi+= flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dfdci[c_index(x, y, 0,
                                                                                                                                       comp,
                                                                                                                                       ion,
                                                                                                                                       Nx,
                                                                                                                                       0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dfdce[c_index(x, y, 0,
                                                                                                                                      comp,
                                                                                                                                      ion,
                                                                                                                                      Nx,
                                                                                                                                      0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume terms
                    //C extra with intra alpha
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), -c[c_index(x, y, 0, Nc - 1,
                                                                                                                                ion, Nx, 0)], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //C intra with intra alpha
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), c[c_index(x, y, 0, comp,
                                                                                                                             ion,
                                                                                                                             Nx, 0)], INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                Ftmpx = 0;
                Fc0x = 0;
                Fc1x = 0;
                Fph0x = 0;
                Fph1x = 0;
                Ftmpy = 0;
                Fc0y = 0;
                Fc1y = 0;
                Fph0y = 0;
                Fph1y = 0;
                if(x<Nx-1) {
                    Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0, comp,
                            ion, Nx, 0)]) / 2 / dx * dt / dx;
                    Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0x = z_charge[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                          0)] + cp[c_index(x, y, 0,
                                                                                                           comp,
                                                                                                           ion, Nx,
                                                                                                           0)]) / 2 / dx * dt / dx;
                    Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1x = z_charge[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx, 0)]) / 2 / dy * dt / dy;
                    Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0y = z_charge[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dy * dt / dy;
                    Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1y = z_charge[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) + Fc0x + Fc1x + Fc0y + Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++) {
                    Ac -= flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1, ion,
                                                                                                   Nx, 0) * 2 + 1], 2));
                Ac -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y, 0, Nc - 1,
                                                                                                       ion, Nx, 0)]) * dt;
                Aphi -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
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
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * al[al_index(x, y, 0,
                                                                                                                                    comp, Nx, 0)], INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                }
                //Phi with C extracellular one
                comp = Nc-1;
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * (1 - al[al_index(x, y, 0,
                                                                                                                                     0,
                                                                                                                                     Nx, 0)] - al[al_index(
                        x, y, 0, 1, Nx, 0)]), INSERT_VALUES); CHKERRQ(ierr);
                ind++;

            }
            //electroneutrality-voltage entries
            Aphi = 0;
            for(comp=0;comp<Nc-1;comp++)
            {
                Aphi -= cm[comp];
            }
            //extraphi with extra phi
            ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
            ind++;
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extra phi with intra phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // Intra phi with Extraphi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra phi with Intra phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Extra phi with intra-Volume
                ierr = MatSetValue(Jac,Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0),Ind_1(x, y, 0, Ni + 1, comp, Nx, 0),-cz(c,z_charge,x,y,0,Nx,0,
                                                                                                                 Nc-
                                                                                                                 1,user),INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra phi with Intra Vol
                ierr = MatSetValue(Jac,Ind_1(x, y, 0, Ni, comp, Nx, 0),Ind_1(x, y, 0, Ni + 1, comp, Nx, 0),cz(c,z_charge,x,y,0,Nx,0,comp,user),INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
        }
    }
    //water flow
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(comp=0;comp<Nc-1;comp++) {
                //Water flow volume fraction entries
                //Volume to Volume
                Ac=1+ (flux->dwdpi[al_index(x, y, 0, comp, Nx, 0)] * (con_vars->ao[phi_index(0, 0, 0, Nc - 1, Nx, 0)] / (pow(1 - al[al_index(
                        x, y, 0, 0, Nx, 0)] - al[al_index(
                        x, y, 0, 1, Nx, 0)], 2)) + con_vars->ao[phi_index(
                        0, 0, 0, comp, Nx, 0)] / pow(al[al_index(x, y, 0, comp, Nx, 0)], 2)) + flux->dwdal[al_index(x,
                                                                                                                    y, 0,
                                                                                                                    comp,
                                                                                                                    Nx,
                                                                                                                    0)]) * dt;
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Off diagonal (from aNc=1-sum(ak))
                for (PetscInt l=0; l<comp; l++) {
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, l, Nx, 0), flux->dwdpi[al_index(x, y, 0,
                                                                                                                                        comp,
                                                                                                                                        Nx,
                                                                                                                                        0)] * con_vars->ao[phi_index(
                            0, 0, 0, Nc - 1, Nx, 0)] / pow(1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1,
                                                                                                             Nx,
                                                                                                             0)], 2) * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (PetscInt l=comp+1; l<Nc-1; l++) {
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, Ni + 1, l, Nx, 0), flux->dwdpi[al_index(x, y, 0,
                                                                                                                                        comp,
                                                                                                                                        Nx,
                                                                                                                                        0)] * con_vars->ao[phi_index(
                            0, 0, 0, Nc - 1, Nx, 0)] / ((1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1,
                                                                                                           Nx, 0)]) * (1 - al[al_index(
                            x, y, 0, 0, Nx, 0)] - al[al_index(
                            x, y, 0, 1, Nx, 0)])) * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                for (ion=0; ion<Ni; ion++) {
                    //Volume to extra c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dwdpi[al_index(x, y,
                                                                                                                                          0,
                                                                                                                                          comp,
                                                                                                                                          Nx,
                                                                                                                                          0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Volume to intra c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni + 1, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dwdpi[al_index(x, y,
                                                                                                                                         0,
                                                                                                                                         comp,
                                                                                                                                         Nx,
                                                                                                                                         0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
            }
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode calc_residual_algebraic_no_vol(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using algebraic version of the charge-capacitance relation
    //Alpha is not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);


    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp, ion,
                                                                                                              Nx, 0)]) / 2;
                        Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x + 1, y,
                                                                                                               0, comp,
                                                                                                               ion, Nx,
                                                                                                               0)]) / 2;
                        RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                                comp,
                                                                                                                ion, Nx,
                                                                                                                0)]) + z_charge[ion] * (phi[phi_index(
                                x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                                 Nx, 0)] + cp[c_index(x,
                                                                                                                      y,
                                                                                                                      0,
                                                                                                                      comp,
                                                                                                                      ion,
                                                                                                                      Nx,
                                                                                                                      0)]) / 2;
                        Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                                x, y, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                               0)] + cp[c_index(x,
                                                                                                                y + 1, 0,
                                                                                                                comp,
                                                                                                                ion,
                                                                                                                Nx,
                                                                                                                0)]) / 2;
                        RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) + z_charge[ion] * (phi[phi_index(
                                x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                                x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                    }
                    Resc = al[al_index(x, y, 0, comp, Nx, 0)] * c[c_index(x, y, 0, comp, ion, Nx, 0)] - alp[al_index(x,
                                                                                                                     y,
                                                                                                                     0,
                                                                                                                     comp,
                                                                                                                     Nx,
                                                                                                                     0)] * cp[c_index(
                            x, y, 0, comp, ion, Nx, 0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;

                    ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)];
                alpNc = 1 - alp[al_index(x, y, 0, 0, Nx, 0)] - alp[al_index(x, y, 0, 1, Nx, 0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                         0)] + cp[c_index(x, y, 0, comp,
                                                                                                          ion, Nx,
                                                                                                          0)]) / 2;
                    Rcvx = Rcvx * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x - 1, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvxRight = RcvxRight * (log(c[c_index(x + 1, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0,
                                                                                                            comp, ion,
                                                                                                            Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x + 1, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dx * dt / dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion, Nx,
                                                                                             0)] + cp[c_index(x, y, 0,
                                                                                                              comp,
                                                                                                              ion, Nx,
                                                                                                              0)]) / 2;
                    Rcvy = Rcvy * (log(c[c_index(x, y, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y - 1, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx,
                            0)]) / 2;
                    RcvyUp = RcvyUp * (log(c[c_index(x, y + 1, 0, comp, ion, Nx, 0)]) - log(c[c_index(x, y, 0, comp,
                                                                                                      ion, Nx, 0)]) + z_charge[ion] * (phi[phi_index(
                            x, y + 1, 0, comp, Nx, 0)] - phi[phi_index(
                            x, y, 0, comp, Nx, 0)])) / dy * dt / dy;
                }
                Resc = alNc*c[c_index(x, y, 0, comp, ion, Nx, 0)] - alpNc * cp[c_index(x, y, 0, comp, ion, Nx, 0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x, y, 0, comp, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, comp, ion,
                                                                                                   Nx, 0) * 2 + 1], 2)) * (cp[c_index(
                        x, y, 0, comp, ion, Nx, 0)] + cbath[ion]) / 2.0 * (log(c[c_index(
                        x, y, 0, comp, ion, Nx, 0)]) - log(cbath[ion]) + z_charge[ion] * phi[phi_index(x, y, 0, comp, Nx, 0)] - z_charge[ion] * phibath) * dt;
                ierr = VecSetValue(Res, Ind_1(x, y, 0, ion, comp, Nx, 0), Resc, INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {

            //Residual for electroneutrality condition
            for(comp=0;comp<Nc-1;comp++)
            {

                Resc = al[al_index(x, y, 0, comp, Nx, 0)] *cz(c,z_charge,x,y,0,Nx,0,comp,user)+user->con_vars->zo[phi_index(0,0,0,comp,
                                                                                                                            Nx,0)]*user->con_vars->ao[phi_index(
                        0, 0, 0, comp, Nx, 0)];
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resc, INSERT_VALUES); CHKERRQ(ierr);
            }
            //Extracellular term
            comp=Nc-1;
            Resc = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) *
                           cz(c,z_charge,x,y,0,Nx,0,comp,user)+user->con_vars->zo[phi_index(
                    0, 0, 0, comp, Nx, 0)]*user->con_vars->ao[phi_index(
                    0, 0, 0, comp, Nx, 0)];
            ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), Resc, INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            // Add Modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extracell voltage
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -cm[comp] * (phi[phi_index(x, y, 0, Nc - 1, Nx, 0)] - phi[phi_index(
                        x, y,
                        0,
                        comp,
                        Nx, 0)]), ADD_VALUES);CHKERRQ(ierr);
                //Intracell voltage mod
                ierr = VecSetValue(Res, Ind_1(x, y, 0, Ni, comp, Nx, 0), -cm[comp] * (phi[phi_index(x, y, 0, comp, Nx, 0)] - phi[phi_index(
                        x, y,
                        0,
                        Nc -
                        1, Nx,
                        0)]), ADD_VALUES);CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_algebraic_no_vol(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using algebraic version of the charge-capacitance relation
    // Alpha is not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi;


    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions
                    Ftmpx = 0;
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Ftmpy = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if(x<Nx-1) {
                        Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                                x + 1, y, 0,
                                comp, ion, Nx,
                                0)]) / 2 / dx * dt / dx;
                        Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                    }
                    if(x>0) {
                        Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dx * dt / dx;
                        Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx,
                                                                                              0)] + cp[c_index(x, y + 1,
                                                                                                               0,
                                                                                                               comp,
                                                                                                               ion,
                                                                                                               Nx,
                                                                                                               0)]) / 2 / dy * dt / dy;
                        Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp,
                                                                                                  ion, Nx, 0)] + cp[c_index(
                                x, y, 0,
                                comp,
                                ion, Nx,
                                0)]) / 2 / dy * dt / dy;
                        Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x, y, 0, comp, Nx, 0)] + Fc0x + Fc1x + Fc0y + Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;


                    //membrane current contributions
                    Ac+= flux->dfdci[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi+= flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -flux->dfdci[c_index(x, y, 0,
                                                                                                                                       comp,
                                                                                                                                       ion,
                                                                                                                                       Nx,
                                                                                                                                       0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), flux->dfdce[c_index(x, y, 0,
                                                                                                                                      comp,
                                                                                                                                      ion,
                                                                                                                                      Nx,
                                                                                                                                      0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), -flux->dfdphim[c_index(x, y,
                                                                                                                                        0,
                                                                                                                                        comp,
                                                                                                                                        ion,
                                                                                                                                        Nx,
                                                                                                                                        0)] * dt, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                Ftmpx = 0;
                Fc0x = 0;
                Fc1x = 0;
                Fph0x = 0;
                Fph1x = 0;
                Ftmpy = 0;
                Fc0y = 0;
                Fc1y = 0;
                Fph0y = 0;
                Fph1y = 0;
                if(x<Nx-1) {
                    Ftmpx = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x + 1, y, 0, comp,
                            ion, Nx, 0)]) / 2 / dx * dt / dx;
                    Fc0x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0x = z_charge[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    Ftmpx = Dcs[c_index(x - 1, y, 0, comp, ion, Nx, 0) * 2] * (cp[c_index(x - 1, y, 0, comp, ion, Nx,
                                                                                          0)] + cp[c_index(x, y, 0,
                                                                                                           comp,
                                                                                                           ion, Nx,
                                                                                                           0)]) / 2 / dx * dt / dx;
                    Fc1x = Ftmpx/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1x = z_charge[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1x, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    Ftmpy = Dcs[c_index(x, y, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y, 0, comp, ion, Nx, 0)] + cp[c_index(
                            x, y + 1, 0,
                            comp, ion, Nx, 0)]) / 2 / dy * dt / dy;
                    Fc0y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph0y = z_charge[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph0y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    Ftmpy = Dcs[c_index(x, y - 1, 0, comp, ion, Nx, 0) * 2 + 1] * (cp[c_index(x, y - 1, 0, comp, ion,
                                                                                              Nx, 0)] + cp[c_index(x, y,
                                                                                                                   0,
                                                                                                                   comp,
                                                                                                                   ion,
                                                                                                                   Nx,
                                                                                                                   0)]) / 2 / dy * dt / dy;
                    Fc1y = Ftmpy/c[c_index(x, y, 0, comp, ion, Nx, 0)];
                    Fph1y = z_charge[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), -Fc1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, 0, ion, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -Fph1y, INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1 - al[al_index(x, y, 0, 0, Nx, 0)] - al[al_index(x, y, 0, 1, Nx, 0)]) + Fc0x + Fc1x + Fc0y + Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++) {
                    Ac -= flux->dfdce[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, 0, comp, ion, Nx, 0)] * dt;
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x, y, 0, Nc - 1, ion, Nx, 0) * 2], 2) + pow(Dcb[c_index(x, y, 0, Nc - 1, ion,
                                                                                                   Nx, 0) * 2 + 1], 2));
                Ac -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) / (2 * c[c_index(x, y, 0, Nc - 1,
                                                                                                       ion, Nx, 0)]) * dt;
                Aphi -= Ftmpx * (cp[c_index(x, y, 0, Nc - 1, ion, Nx, 0)] + cbath[ion]) * z_charge[ion] / 2 * dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ac, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, ion, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }

        }
    }

    //Electroneutrality charge-capacitence condition
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
                    ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * al[al_index(x, y, 0,
                                                                                                                                    comp, Nx, 0)], INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                }
                //Phi with C extracellular one
                comp = Nc-1;
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, ion, comp, Nx, 0), z_charge[ion] * (1 - al[al_index(x, y, 0,
                                                                                                                                     0,
                                                                                                                                     Nx, 0)] - al[al_index(
                        x, y, 0, 1, Nx, 0)]), INSERT_VALUES); CHKERRQ(ierr);
                ind++;

            }
            //electroneutrality-voltage entries
            Aphi = 0;
            for(comp=0;comp<Nc-1;comp++)
            {
                Aphi -= cm[comp];
            }
            //extraphi with extra phi
            ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Aphi, INSERT_VALUES);CHKERRQ(ierr);
            ind++;
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extra phi with intra phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // Intra phi with Extraphi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, Nc - 1, Nx, 0), cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra phi with Intra phi
                ierr = MatSetValue(Jac, Ind_1(x, y, 0, Ni, comp, Nx, 0), Ind_1(x, y, 0, Ni, comp, Nx, 0), -cm[comp], INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}



