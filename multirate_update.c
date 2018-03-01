#include "constants.h"
#include "functions.h"

void recombine(struct SimState* state_vars,struct AppCtx *user){
    extract_subarray(state_vars->v,state_vars,1);

    PetscInt x,y,comp;
    PetscInt Nx = user->Nx;
    for(x=0;x<Nx;x++){
        for ( y = 0; y < user->Ny; y++) {
            for (comp = 0; comp < Nc; ++comp) {
                state_vars->phi[phi_index(x,y,comp,Nx)]+=state_vars->phi_fast[phi_index(x,y,comp,Nx)];
                state_vars->phi_fast[phi_index(x,y,comp,Nx)]=0;
            }

        }

    }

    restore_subarray(state_vars->v,state_vars,1);
}

PetscErrorCode calc_residual_slow(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray(current_state,user->state_vars,1); CHKERRQ(ierr);

    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *cfast = user->state_vars->c_fast;
    PetscReal *phifast = user->state_vars->phi_fast;

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
                        Rcvx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfast[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x-1,y,comp,ion,Nx)])
                                     +z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)]-phi[phi_index(x-1,y,comp,Nx)]-phifast[phi_index(x-1,y,comp,Nx)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                        RcvxRight += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])
                                               +z[ion]*(phi[phi_index(x+1,y,comp,Nx)]+phifast[phi_index(x+1,y,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]-phifast[phi_index(x,y,comp,Nx)]))/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfast[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y-1,comp,ion,Nx)])
                                     +z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)]-phi[phi_index(x,y-1,comp,Nx)]-phifast[phi_index(x,y-1,comp,Nx)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                        RcvyUp += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])
                                         +z[ion]*(phi[phi_index(x,y+1,comp,Nx)]+phifast[phi_index(x,y+1,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]-phifast[phi_index(x,y,comp,Nx)]))/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,comp,Nx)]*c[c_index(x,y,comp,ion,Nx)]-alp[al_index(x,y,comp,Nx)]*cp[c_index(x,y,comp,ion,Nx)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion,Nx)]*dt;

                    ierr = VecSetValue(Res,Ind_1(x,y,ion,comp,Nx),Resc,INSERT_VALUES);CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp]+=z[ion]*Rcvx;
                    Rphy[comp]+=z[ion]*Rcvy;
                    RphxRight[comp]+=z[ion]*RcvxRight;
                    RphyUp[comp]+=z[ion]*RcvyUp;

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x,y,0,Nx)] - al[al_index(x,y,1,Nx)];
                alpNc = 1 - alp[al_index(x,y,0,Nx)] - alp[al_index(x,y,1,Nx)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfast[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x-1,y,comp,ion,Nx)])
                                 +z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)]-phi[phi_index(x-1,y,comp,Nx)]-phifast[phi_index(x-1,y,comp,Nx)]))/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                    RcvxRight += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])
                                           +z[ion]*(phi[phi_index(x+1,y,comp,Nx)]+phifast[phi_index(x+1,y,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]-phifast[phi_index(x,y,comp,Nx)]))/dx*dt/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfast[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y-1,comp,ion,Nx)])
                                 +z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)]-phi[phi_index(x,y-1,comp,Nx)]-phifast[phi_index(x,y-1,comp,Nx)]))/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                    RcvyUp += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])
                                     +z[ion]*(phi[phi_index(x,y+1,comp,Nx)]+phifast[phi_index(x,y+1,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]-phifast[phi_index(x,y,comp,Nx)]))/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,comp,ion,Nx)]-alpNc*cp[c_index(x,y,comp,ion,Nx)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion,Nx)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0
                        *(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)])-z[ion]*phibath)*dt;
                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp,Nx),Resc,INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z[ion]*Rcvx;
                Rphy[comp]+=z[ion]*Rcvy;
                RphxRight[comp]+=z[ion]*RcvxRight;
                RphyUp[comp]+=z[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[comp]*(phi[phi_index(x,y,comp,Nx)]-phi[phi_index(x,y,Nc-1,Nx)])-cm[comp]*(phip[phi_index(x,y,comp,Nx)]-phip[phi_index(x,y,Nc-1,Nx)]);
                for(ion=0;ion<Ni;ion++){
                    //Ion channel
                    Resph +=z[ion]*flux->mflux[c_index(x,y,comp,ion,Nx)]*dt;
                }
                //Add the terms shared with extracell
                ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                Resph += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp,Nx),Resph,INSERT_VALUES); CHKERRQ(ierr);
            }

            //Finish adding extracell
            comp = Nc-1;
            //Add bath contribution
            for(ion=0;ion<Ni;ion++){

                ResphN -=z[ion]*sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0
                         *(log(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*(phi[phi_index(x,y,comp,Nx)]+phifast[phi_index(x,y,comp,Nx)])-z[ion]*phibath)*dt;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp,Nx),ResphN,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray(current_state,user->state_vars,1); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_slow(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray(current_state,user->state_vars,1); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *cfast = user->state_vars->c_fast;
    PetscReal *phifast = user->state_vars->phi_fast;

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
                        Ftmpx = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2/dx*dt/dx;
                        Ftmpx += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])/2/dx*dt/dx;
                        Fc0x = Ftmpx/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                        Fph0x = z[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0) {
                        Ftmpx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2/dx*dt/dx;
                        Ftmpx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfast[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                        Fph1x = z[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Ftmpy = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2/dy*dt/dy;
                        Ftmpy += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])/2/dy*dt/dy;
                        Fc0y = Ftmpy/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                        Fph0y = z[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Ftmpy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2/dy*dt/dy;
                        Ftmpy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfast[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2/dy*dt/dy;
                        Fc1y = Ftmpy/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                        Fph1y = z[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x,y,comp,Nx)]+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Add up terms for voltage eqns
                    Fphph0x[comp]+=z[ion]*Fph0x;
                    Fphph1x[comp]+=z[ion]*Fph1x;
                    Fphph0y[comp]+=z[ion]*Fph0y;
                    Fphph1y[comp]+=z[ion]*Fph1y;


                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),Ac,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),z[ion]*(Fc0x+Fc1x+Fc0y+Fc1y),INSERT_VALUES); CHKERRQ(ierr);
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
                    Ftmpx = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2/dx*dt/dx;
                    Ftmpx += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x+1,y,comp,ion,Nx)])/2/dx*dt/dx;
                    Fc0x = Ftmpx/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                    Fph0x = z[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right phi with left c in voltage eqn
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    Ftmpx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2/dx*dt/dx;
                    Ftmpx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfast[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2/dx*dt/dx;
                    Fc1x = Ftmpx/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                    Fph1x = z[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left phi with right c in voltage eqn
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    Ftmpy = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2/dy*dt/dy;
                    Ftmpy += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfast[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y+1,comp,ion,Nx)])/2/dy*dt/dy;
                    Fc0y = Ftmpy/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                    Fph0y = z[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper phi with lower c in voltage eqn
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    Ftmpy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2/dy*dt/dy;
                    Ftmpy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfast[c_index(x,y-1,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])/2/dy*dt/dy;
                    Fc1y = Ftmpy/(c[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)]);
                    Fph1y = z[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower phi with upper c in voltage eqn
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1-al[al_index(x,y,0,Nx)]-al[al_index(x,y,1,Nx)])+Fc0x+Fc1x+Fc0y+Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                Avolt = z[ion]*(Fc0x+Fc1x+Fc0y+Fc1y);

                //Add up terms for voltage eqns
                Fphph0x[comp]+=z[ion]*Fph0x;
                Fphph1x[comp]+=z[ion]*Fph1x;
                Fphph0y[comp]+=z[ion]*Fph0y;
                Fphph1y[comp]+=z[ion]*Fph1y;

                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2+1],2));
                Ac -= Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])/(2*(c[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]))*dt;
                Aphi -= Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])*z[ion]/2*dt;

                Avolt -=z[ion]*Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])/(2*(c[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]))*dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,ion,Nc-1,Nx),Ac,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,Ni,Nc-1,Nx),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1,Nx),Ind_1(x,y,ion,Nc-1,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for(comp=0;comp<Nc-1;comp++) {
                if(x<Nx-1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph0x[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph1x[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph0y[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph1y[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                Avolt = cm[comp]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];
                AvoltN = -cm[comp];

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,Nc-1,Nx),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc-1;
            if(x<Nx-1) {
                //Right phi with left phi (-Fph0x)
                ierr = MatSetValue(Jac,Ind_1(x+1,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph0x[comp],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(x>0) {
                //Left phi with right phi (-Fph1x)
                ierr = MatSetValue(Jac,Ind_1(x-1,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph1x[comp],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y<Ny-1) {
                //Upper phi with lower phi (-Fph0y)
                ierr = MatSetValue(Jac,Ind_1(x,y+1,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph0y[comp],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            if(y>0) {
                //Lower phi with upper phi (-Fph1y)
                ierr = MatSetValue(Jac,Ind_1(x,y-1,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),-Fphph1y[comp],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            AvoltN = 0;

            for(int k=0;k<Nc-1;k++) {
                AvoltN += cm[k];
                Avolt = -cm[k];
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,k,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }

            AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];

            //Bath terms
            for(ion=0;ion<Ni;ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2+1],2));
                AvoltN -= z[ion]*Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])*z[ion]/2*dt;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
            ind++;

        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray(current_state,user->state_vars,1); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}
PetscErrorCode calc_residual_fast(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray(current_state,user->state_vars,2); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);

    PetscReal *cfast = user->state_vars->c_fast;
    PetscReal *phifast = user->state_vars->phi_fast;
    PetscReal *cfastp = user->state_vars_past->c_fast;
    PetscReal *phifastp = user->state_vars_past->phi_fast;

    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dtf = user->dtf;
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
                        Rcvx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfastp[c_index(x-1,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvx = Rcvx*(log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])-log(cp[c_index(x-1,y,comp,ion,Nx)]+cfast[c_index(x-1,y,comp,ion,Nx)])
                                     +z[ion]*(phip[phi_index(x,y,comp,Nx)]+phifastp[phi_index(x,y,comp,Nx)]-phip[phi_index(x-1,y,comp,Nx)]-phifastp[phi_index(x-1,y,comp,Nx)]))/dx*dtf/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                        RcvxRight += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfastp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x+1,y,comp,ion,Nx)])/2;
                        RcvxRight = RcvxRight*(log(cp[c_index(x+1,y,comp,ion,Nx)]+cfastp[c_index(x+1,y,comp,ion,Nx)])-log(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])
                                               +z[ion]*(phip[phi_index(x+1,y,comp,Nx)]+phifastp[phi_index(x+1,y,comp,Nx)]-phip[phi_index(x,y,comp,Nx)]-phifastp[phi_index(x,y,comp,Nx)]))/dx*dtf/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfastp[c_index(x,y-1,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvy = Rcvy*(log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])-log(cp[c_index(x,y-1,comp,ion,Nx)]+cfastp[c_index(x,y-1,comp,ion,Nx)])
                                     +z[ion]*(phip[phi_index(x,y,comp,Nx)]+phifastp[phi_index(x,y,comp,Nx)]-phip[phi_index(x,y-1,comp,Nx)]-phifastp[phi_index(x,y-1,comp,Nx)]))/dy*dtf/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                        RcvyUp += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfastp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y+1,comp,ion,Nx)])/2;
                        RcvyUp = RcvyUp*(log(cp[c_index(x,y+1,comp,ion,Nx)]+cfastp[c_index(x,y+1,comp,ion,Nx)])-log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])
                                         +z[ion]*(phip[phi_index(x,y+1,comp,Nx)]+phifastp[phi_index(x,y+1,comp,Nx)]-phip[phi_index(x,y,comp,Nx)]-phifastp[phi_index(x,y,comp,Nx)]))/dy*dtf/dy;
                    }
                    Resc = alp[al_index(x,y,comp,Nx)]*cfast[c_index(x,y,comp,ion,Nx)]-alp[al_index(x,y,comp,Nx)]*cfastp[c_index(x,y,comp,ion,Nx)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion,Nx)]*dtf;

                    ierr = VecSetValue(Res,Ind_1(x,y,ion,comp,Nx),Resc,INSERT_VALUES);CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp]+=z[ion]*Rcvx;
                    Rphy[comp]+=z[ion]*Rcvy;
                    RphxRight[comp]+=z[ion]*RcvxRight;
                    RphyUp[comp]+=z[ion]*RcvyUp;

                }
                //Set Extracellular values
                alpNc = 1 - alp[al_index(x,y,0,Nx)] - alp[al_index(x,y,1,Nx)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cp[c_index(x-1,y,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvx += Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(cfastp[c_index(x-1,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvx = Rcvx*(log(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(cp[c_index(x-1,y,comp,ion,Nx)]+cfastp[c_index(x-1,y,comp,ion,Nx)])
                                 +z[ion]*(phip[phi_index(x,y,comp,Nx)]+phifastp[phi_index(x,y,comp,Nx)]-phip[phi_index(x-1,y,comp,Nx)]-phifastp[phi_index(x-1,y,comp,Nx)]))/dx*dtf/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                    RcvxRight += Dcs[c_index(x,y,comp,ion,Nx)*2]*(cfastp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x+1,y,comp,ion,Nx)])/2;
                    RcvxRight = RcvxRight*(log(cp[c_index(x+1,y,comp,ion,Nx)]+cfastp[c_index(x+1,y,comp,ion,Nx)])-log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])
                                           +z[ion]*(phip[phi_index(x+1,y,comp,Nx)]+phifastp[phi_index(x+1,y,comp,Nx)]-phip[phi_index(x,y,comp,Nx)]-phifastp[phi_index(x,y,comp,Nx)]))/dx*dtf/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvy += Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cfastp[c_index(x,y-1,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvy = Rcvy*(log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])-log(cp[c_index(x,y-1,comp,ion,Nx)]+cfastp[c_index(x,y-1,comp,ion,Nx)])
                                 +z[ion]*(phip[phi_index(x,y,comp,Nx)]+phifastp[phi_index(x,y,comp,Nx)]-phip[phi_index(x,y-1,comp,Nx)]-phifastp[phi_index(x,y-1,comp,Nx)]))/dy*dtf/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                    RcvyUp += Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cfastp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y+1,comp,ion,Nx)])/2;
                    RcvyUp = RcvyUp*(log(cp[c_index(x,y+1,comp,ion,Nx)]+cfastp[c_index(x,y+1,comp,ion,Nx)])-log(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)])
                                     +z[ion]*(phip[phi_index(x,y+1,comp,Nx)]+phifastp[phi_index(x,y+1,comp,Nx)]-phip[phi_index(x,y,comp,Nx)]-phifastp[phi_index(x,y,comp,Nx)]))/dy*dtf/dy;
                }
                Resc = alpNc*cfast[c_index(x,y,comp,ion,Nx)]-alpNc*cfastp[c_index(x,y,comp,ion,Nx)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion,Nx)]*dtf;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0
                        *(log(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*phip[phi_index(x,y,comp,Nx)]+z[ion]*phifast[phi_index(x,y,comp,Nx)]-z[ion]*phibath)*dtf;
                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp,Nx),Resc,INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z[ion]*Rcvx;
                Rphy[comp]+=z[ion]*Rcvy;
                RphxRight[comp]+=z[ion]*RcvxRight;
                RphyUp[comp]+=z[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[comp]*(phifast[phi_index(x,y,comp,Nx)]-phifast[phi_index(x,y,Nc-1,Nx)])-cm[comp]*(phifastp[phi_index(x,y,comp,Nx)]-phifastp[phi_index(x,y,Nc-1,Nx)]);
                for(ion=0;ion<Ni;ion++){
                    //Ion channel
                    Resph +=z[ion]*flux->mflux[c_index(x,y,comp,ion,Nx)]*dtf;
                }
                //Add the terms shared with extracell
                ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                Resph += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp,Nx),Resph,INSERT_VALUES); CHKERRQ(ierr);
            }

            //Finish adding extracell
            comp = Nc-1;
            //Add bath contribution
            for(ion=0;ion<Ni;ion++){

                ResphN -=z[ion]*sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cfastp[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0
                         *(log(cp[c_index(x,y,comp,ion,Nx)]+cfast[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*phip[phi_index(x,y,comp,Nx)]+z[ion]*phifast[phi_index(x,y,comp,Nx)]-z[ion]*phibath)*dtf;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp,Nx),ResphN,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray(current_state,user->state_vars,2); CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_fast(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray(current_state,user->state_vars,2); CHKERRQ(ierr);
    PetscReal *cfast = user->state_vars->c_fast;
    PetscReal *cfastp = user->state_vars_past->c_fast;


    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dtf = user->dtf;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Diagonal term contribution
                    Ac = alp[al_index(x,y,comp,Nx)];
                    Aphi = 0;

                    //membrane current contributions
                    Ac+=flux->dfdci[c_index(x,y,comp,ion,Nx)]*dtf;
                    Aphi+=flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,ion,comp,Nx),-flux->dfdci[c_index(x,y,comp,ion,Nx)]*dtf,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,ion,Nc-1,Nx),flux->dfdce[c_index(x,y,comp,ion,Nx)]*dtf,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,Ni,comp,Nx),-flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,Ni,Nc-1,Nx),-flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,ion,comp,Nx),Ac,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp,Nx),Ind_1(x,y,Ni,comp,Nx),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,ion,comp,Nx),z[ion]*(flux->dfdci[c_index(x,y,comp,ion,Nx)]*dtf),INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //IntraPhi with c extra(volt eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,ion,Nc-1,Nx),z[ion]*(flux->dfdce[c_index(x,y,comp,ion,Nx)]*dtf),INSERT_VALUES); CHKERRQ(ierr);
                    ind++;
                    //Extra-Phi with intra-c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1,Nx),Ind_1(x,y,ion,comp,Nx),-z[ion]*(flux->dfdci[c_index(x,y,comp,ion,Nx)]*dtf),INSERT_VALUES); CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                //Diagonal term contribution
                Ac = (1-alp[al_index(x,y,0,Nx)]-alp[al_index(x,y,1,Nx)]);
                Aphi = 0;

                Avolt = 0;

                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++) {
                    Ac -= flux->dfdce[c_index(x,y,comp,ion,Nx)]*dtf;
                    Aphi += flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf;
                    Avolt -=z[ion]*flux->dfdce[c_index(x,y,comp,ion,Nx)]*dtf;
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2+1],2));
                Ac -= Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfastp[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])/(2*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]))*dtf;
                Aphi -= Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfastp[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])*z[ion]/2*dtf;

                Avolt -=z[ion]*Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfastp[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])/(2*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfast[c_index(x,y,Nc-1,ion,Nx)]))*dtf;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,ion,Nc-1,Nx),Ac,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1,Nx),Ind_1(x,y,Ni,Nc-1,Nx),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1,Nx),Ind_1(x,y,ion,Nc-1,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for(comp=0;comp<Nc-1;comp++) {
                Avolt = cm[comp];
                AvoltN = -cm[comp];
                for(ion=0;ion<Ni;ion++) {
                    Avolt+=z[ion]*flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf;
                    AvoltN-=z[ion]*flux->dfdphim[c_index(x,y,comp,ion,Nx)]*dtf;
                }

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,Nc-1,Nx),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc-1;
            AvoltN = 0;

            for(int k=0;k<Nc-1;k++) {
                AvoltN += cm[k];
                Avolt = -cm[k];
                for(ion=0;ion<Ni;ion++) {
                    Avolt-=z[ion]*flux->dfdphim[c_index(x,y,k,ion,Nx)]*dtf;
                    AvoltN+=z[ion]*flux->dfdphim[c_index(x,y,k,ion,Nx)]*dtf;
                }
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,k,Nx),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
            }

            //Bath terms
            for(ion=0;ion<Ni;ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,Nc-1,ion,Nx)*2+1],2));
                AvoltN -= z[ion]*Ftmpx*(cp[c_index(x,y,Nc-1,ion,Nx)]+cfastp[c_index(x,y,Nc-1,ion,Nx)]+cbath[ion])*z[ion]/2*dtf;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp,Nx),Ind_1(x,y,Ni,comp,Nx),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
            ind++;

        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray(current_state,user->state_vars,2); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}