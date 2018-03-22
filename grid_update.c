#include "constants.h"
#include "functions.h"

void Load_Grid(struct AppCtx *user,PetscInt xi,PetscInt yi){

    struct SimState *state_vars = user->state_vars_past;
    struct SimState * grid_vars = user->grid_vars;
    struct GateType *gate_vars = user->gate_vars;
    struct GateType *grid_gate = user->grid_gate_vars;

    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt nx = 2*width_size+1;
    PetscInt ny = 2*width_size+1;
    PetscInt xind,yind,ion,comp,x,y;

    PetscInt neighbor_count;

    for( x=0;x<nx;x++) {
        for ( y = 0; y < ny; y++) {
            xind = x-width_size+xi;

            yind = y-width_size+yi;

            //If in interior just copy
            if(xind>0 && yind>0 && xind<Nx-1 && yind<Ny-1) {

                grid_gate->mNaT[xy_index(x, y, nx)] = gate_vars->mNaT[xy_index(xind, yind, Nx)];
                grid_gate->hNaT[xy_index(x, y, nx)] = gate_vars->hNaT[xy_index(xind, yind, Nx)];
                grid_gate->gNaT[xy_index(x, y, nx)] = gate_vars->gNaT[xy_index(xind, yind, Nx)];
                grid_gate->mNaP[xy_index(x, y, nx)] = gate_vars->mNaP[xy_index(xind, yind, Nx)];
                grid_gate->hNaP[xy_index(x, y, nx)] = gate_vars->hNaP[xy_index(xind, yind, Nx)];
                grid_gate->gNaP[xy_index(x, y, nx)] = gate_vars->gNaP[xy_index(xind, yind, Nx)];
                grid_gate->mKDR[xy_index(x, y, nx)] = gate_vars->mKDR[xy_index(xind, yind, Nx)];
                grid_gate->gKDR[xy_index(x, y, nx)] = gate_vars->gKDR[xy_index(xind, yind, Nx)];
                grid_gate->mKA[xy_index(x, y, nx)] = gate_vars->mKA[xy_index(xind, yind, Nx)];
                grid_gate->hKA[xy_index(x, y, nx)] = gate_vars->hKA[xy_index(xind, yind, Nx)];
                grid_gate->gKA[xy_index(x, y, nx)] = gate_vars->gKA[xy_index(xind, yind, Nx)];

                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        grid_vars->c[c_index(x, y, comp, ion, nx)] = state_vars->c[c_index(xind, yind, comp, ion, Nx)];
                    }
                    grid_vars->phi[phi_index(x, y, comp, nx)] = state_vars->phi[phi_index(xind, yind, comp, Nx)];
                }
                for (comp = 0; comp < Nc - 1; comp++) {
                    grid_vars->alpha[al_index(x, y, comp, nx)] = state_vars->alpha[al_index(xind, yind, comp, Nx)];
                }
            } else{
                //If not the interior try to average over neighbors
                neighbor_count = 0;
                grid_gate->mNaT[xy_index(x, y, nx)] = 0;
                grid_gate->hNaT[xy_index(x, y, nx)] = 0;
                grid_gate->gNaT[xy_index(x, y, nx)] = 0;
                grid_gate->mNaP[xy_index(x, y, nx)] = 0;
                grid_gate->hNaP[xy_index(x, y, nx)] = 0;
                grid_gate->gNaP[xy_index(x, y, nx)] = 0;
                grid_gate->mKDR[xy_index(x, y, nx)] = 0;
                grid_gate->gKDR[xy_index(x, y, nx)] = 0;
                grid_gate->mKA[xy_index(x, y, nx)] = 0;
                grid_gate->hKA[xy_index(x, y, nx)] = 0;
                grid_gate->gKA[xy_index(x, y, nx)] = 0;

                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        grid_vars->c[c_index(x, y, comp, ion, nx)] = 0;
                    }
                    grid_vars->phi[phi_index(x, y, comp, nx)] = 0;
                }
                for (comp = 0; comp < Nc - 1; comp++) {
                    grid_vars->alpha[al_index(x, y, comp, nx)] = 0;
                }
                if(xind>0){
                    neighbor_count++;
                    grid_gate->mNaT[xy_index(x, y, nx)] += gate_vars->mNaT[xy_index(xind-1, yind, Nx)];
                    grid_gate->hNaT[xy_index(x, y, nx)] += gate_vars->hNaT[xy_index(xind-1, yind, Nx)];
                    grid_gate->gNaT[xy_index(x, y, nx)] += gate_vars->gNaT[xy_index(xind-1, yind, Nx)];
                    grid_gate->mNaP[xy_index(x, y, nx)] += gate_vars->mNaP[xy_index(xind-1, yind, Nx)];
                    grid_gate->hNaP[xy_index(x, y, nx)] += gate_vars->hNaP[xy_index(xind-1, yind, Nx)];
                    grid_gate->gNaP[xy_index(x, y, nx)] += gate_vars->gNaP[xy_index(xind-1, yind, Nx)];
                    grid_gate->mKDR[xy_index(x, y, nx)] += gate_vars->mKDR[xy_index(xind-1, yind, Nx)];
                    grid_gate->gKDR[xy_index(x, y, nx)] += gate_vars->gKDR[xy_index(xind-1, yind, Nx)];
                    grid_gate->mKA[xy_index(x, y, nx)] += gate_vars->mKA[xy_index(xind-1, yind, Nx)];
                    grid_gate->hKA[xy_index(x, y, nx)] += gate_vars->hKA[xy_index(xind-1, yind, Nx)];
                    grid_gate->gKA[xy_index(x, y, nx)] += gate_vars->gKA[xy_index(xind-1, yind, Nx)];

                    for (comp = 0; comp < Nc; comp++) {
                        for (ion = 0; ion < Ni; ion++) {
                            grid_vars->c[c_index(x, y, comp, ion, nx)] += state_vars->c[c_index(xind-1, yind, comp, ion, Nx)];
                        }
                        grid_vars->phi[phi_index(x, y, comp, nx)] += state_vars->phi[phi_index(xind-1, yind, comp, Nx)];
                    }
                    for (comp = 0; comp < Nc - 1; comp++) {
                        grid_vars->alpha[al_index(x, y, comp, nx)] += state_vars->alpha[al_index(xind-1, yind, comp, Nx)];
                    }
                }
                if(yind>0){
                    neighbor_count++;
                    grid_gate->mNaT[xy_index(x, y, nx)] += gate_vars->mNaT[xy_index(xind, yind-1, Nx)];
                    grid_gate->hNaT[xy_index(x, y, nx)] += gate_vars->hNaT[xy_index(xind, yind-1, Nx)];
                    grid_gate->gNaT[xy_index(x, y, nx)] += gate_vars->gNaT[xy_index(xind, yind-1, Nx)];
                    grid_gate->mNaP[xy_index(x, y, nx)] += gate_vars->mNaP[xy_index(xind, yind-1, Nx)];
                    grid_gate->hNaP[xy_index(x, y, nx)] += gate_vars->hNaP[xy_index(xind, yind-1, Nx)];
                    grid_gate->gNaP[xy_index(x, y, nx)] += gate_vars->gNaP[xy_index(xind, yind-1, Nx)];
                    grid_gate->mKDR[xy_index(x, y, nx)] += gate_vars->mKDR[xy_index(xind, yind-1, Nx)];
                    grid_gate->gKDR[xy_index(x, y, nx)] += gate_vars->gKDR[xy_index(xind, yind-1, Nx)];
                    grid_gate->mKA[xy_index(x, y, nx)] += gate_vars->mKA[xy_index(xind, yind-1, Nx)];
                    grid_gate->hKA[xy_index(x, y, nx)] += gate_vars->hKA[xy_index(xind, yind-1, Nx)];
                    grid_gate->gKA[xy_index(x, y, nx)] += gate_vars->gKA[xy_index(xind, yind-1, Nx)];

                    for (comp = 0; comp < Nc; comp++) {
                        for (ion = 0; ion < Ni; ion++) {
                            grid_vars->c[c_index(x, y, comp, ion, nx)] += state_vars->c[c_index(xind, yind-1, comp, ion, Nx)];
                        }
                        grid_vars->phi[phi_index(x, y, comp, nx)] += state_vars->phi[phi_index(xind, yind-1, comp, Nx)];
                    }
                    for (comp = 0; comp < Nc - 1; comp++) {
                        grid_vars->alpha[al_index(x, y, comp, nx)] += state_vars->alpha[al_index(xind, yind-1, comp, Nx)];
                    }
                }
                if(xind<Nx-1){
                    neighbor_count++;
                    grid_gate->mNaT[xy_index(x, y, nx)] += gate_vars->mNaT[xy_index(xind+1, yind, Nx)];
                    grid_gate->hNaT[xy_index(x, y, nx)] += gate_vars->hNaT[xy_index(xind+1, yind, Nx)];
                    grid_gate->gNaT[xy_index(x, y, nx)] += gate_vars->gNaT[xy_index(xind+1, yind, Nx)];
                    grid_gate->mNaP[xy_index(x, y, nx)] += gate_vars->mNaP[xy_index(xind+1, yind, Nx)];
                    grid_gate->hNaP[xy_index(x, y, nx)] += gate_vars->hNaP[xy_index(xind+1, yind, Nx)];
                    grid_gate->gNaP[xy_index(x, y, nx)] += gate_vars->gNaP[xy_index(xind+1, yind, Nx)];
                    grid_gate->mKDR[xy_index(x, y, nx)] += gate_vars->mKDR[xy_index(xind+1, yind, Nx)];
                    grid_gate->gKDR[xy_index(x, y, nx)] += gate_vars->gKDR[xy_index(xind+1, yind, Nx)];
                    grid_gate->mKA[xy_index(x, y, nx)] += gate_vars->mKA[xy_index(xind+1, yind, Nx)];
                    grid_gate->hKA[xy_index(x, y, nx)] += gate_vars->hKA[xy_index(xind+1, yind, Nx)];
                    grid_gate->gKA[xy_index(x, y, nx)] += gate_vars->gKA[xy_index(xind+1, yind, Nx)];

                    for (comp = 0; comp < Nc; comp++) {
                        for (ion = 0; ion < Ni; ion++) {
                            grid_vars->c[c_index(x, y, comp, ion, nx)] += state_vars->c[c_index(xind+1, yind, comp, ion, Nx)];
                        }
                        grid_vars->phi[phi_index(x, y, comp, nx)] += state_vars->phi[phi_index(xind+1, yind, comp, Nx)];
                    }
                    for (comp = 0; comp < Nc - 1; comp++) {
                        grid_vars->alpha[al_index(x, y, comp, nx)] += state_vars->alpha[al_index(xind+1, yind, comp, Nx)];
                    }
                }
                if(yind<Ny-1){
                    neighbor_count++;
                    grid_gate->mNaT[xy_index(x, y, nx)] += gate_vars->mNaT[xy_index(xind, yind+1, Nx)];
                    grid_gate->hNaT[xy_index(x, y, nx)] += gate_vars->hNaT[xy_index(xind, yind+1, Nx)];
                    grid_gate->gNaT[xy_index(x, y, nx)] += gate_vars->gNaT[xy_index(xind, yind+1, Nx)];
                    grid_gate->mNaP[xy_index(x, y, nx)] += gate_vars->mNaP[xy_index(xind, yind+1, Nx)];
                    grid_gate->hNaP[xy_index(x, y, nx)] += gate_vars->hNaP[xy_index(xind, yind+1, Nx)];
                    grid_gate->gNaP[xy_index(x, y, nx)] += gate_vars->gNaP[xy_index(xind, yind+1, Nx)];
                    grid_gate->mKDR[xy_index(x, y, nx)] += gate_vars->mKDR[xy_index(xind, yind+1, Nx)];
                    grid_gate->gKDR[xy_index(x, y, nx)] += gate_vars->gKDR[xy_index(xind, yind+1, Nx)];
                    grid_gate->mKA[xy_index(x, y, nx)] += gate_vars->mKA[xy_index(xind, yind+1, Nx)];
                    grid_gate->hKA[xy_index(x, y, nx)] += gate_vars->hKA[xy_index(xind, yind+1, Nx)];
                    grid_gate->gKA[xy_index(x, y, nx)] += gate_vars->gKA[xy_index(xind, yind+1, Nx)];

                    for (comp = 0; comp < Nc; comp++) {
                        for (ion = 0; ion < Ni; ion++) {
                            grid_vars->c[c_index(x, y, comp, ion, nx)] += state_vars->c[c_index(xind, yind+1, comp, ion, Nx)];
                        }
                        grid_vars->phi[phi_index(x, y, comp, nx)] += state_vars->phi[phi_index(xind, yind+1, comp, Nx)];
                    }
                    for (comp = 0; comp < Nc - 1; comp++) {
                        grid_vars->alpha[al_index(x, y, comp, nx)] += state_vars->alpha[al_index(xind, yind+1, comp, Nx)];
                    }
                }

                //Now divide (we exclude the diagonals)
                grid_gate->mNaT[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->hNaT[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->gNaT[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->mNaP[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->hNaP[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->gNaP[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->mKDR[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->gKDR[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->mKA[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->hKA[xy_index(x, y, nx)] /= neighbor_count;
                grid_gate->gKA[xy_index(x, y, nx)] /= neighbor_count;

                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        grid_vars->c[c_index(x, y, comp, ion, nx)] /= neighbor_count;
                    }
                    grid_vars->phi[phi_index(x, y, comp, nx)] /= neighbor_count;
                }
                for (comp = 0; comp < Nc - 1; comp++) {
                    grid_vars->alpha[al_index(x, y, comp, nx)] /= neighbor_count;
                }

            }


        }
    }

}

PetscErrorCode Grid_Residual(Vec Res,PetscInt xi,PetscInt yi,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    //Compute membrane ionic flux relation quantitites
    grid_ionmflux(user);

    //Compute membrane water flow related quantities
    grid_wflowm(user);

    PetscReal *c = user->grid_vars->c;
    PetscReal *phi = user->grid_vars->phi;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;
    PetscReal *alp = user->grid_vars_past->alpha;
    PetscReal *phip = user->grid_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;

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
                        Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x-1,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y,comp,Nx)]-phi[phi_index(x-1,y,comp,Nx)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x+1,y,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]))/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                        Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x,y-1,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y,comp,Nx)]-phi[phi_index(x,y-1,comp,Nx)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y+1,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]))/dy*dt/dy;
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
                    Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x-1,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y,comp,Nx)]-phi[phi_index(x-1,y,comp,Nx)]))/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x,y,comp,ion,Nx)*2]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x+1,y,comp,ion,Nx)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x+1,y,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]))/dx*dt/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(cp[c_index(x,y-1,comp,ion,Nx)]+cp[c_index(x,y,comp,ion,Nx)])/2;
                    Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion,Nx)])-log(c[c_index(x,y-1,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y,comp,Nx)]-phi[phi_index(x,y-1,comp,Nx)]))/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(cp[c_index(x,y,comp,ion,Nx)]+cp[c_index(x,y+1,comp,ion,Nx)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion,Nx)])-log(c[c_index(x,y,comp,ion,Nx)])+z[ion]*(phi[phi_index(x,y+1,comp,Nx)]-phi[phi_index(x,y,comp,Nx)]))/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,comp,ion,Nx)]-alpNc*cp[c_index(x,y,comp,ion,Nx)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion,Nx)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp,Nx)]-z[ion]*phibath)*dt;
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

                ResphN -=z[ion]*sqrt(pow(Dcb[c_index(x,y,comp,ion,Nx)*2],2)+pow(Dcb[c_index(x,y,comp,ion,Nx)*2+1],2))*(cp[c_index(x,y,comp,ion,Nx)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion,Nx)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp,Nx)]-z[ion]*phibath)*dt;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp,Nx),ResphN,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode Grid_Jacobian(Mat Jac,PetscInt xi,PetscInt yi,void *ctx) {
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx *user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if (Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    PetscReal *c = user->grid_vars->c;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    struct ConstVars *con_vars = user->con_vars;

    PetscInt ind = 0;
    PetscInt x, y, ion, comp;

    PetscReal Ftmpx, Fc0x, Fc1x, Fph0x, Fph1x;
    PetscReal Ftmpy, Fc0y, Fc1y, Fph0y, Fph1y;
    PetscReal Ac, Aphi, Avolt, AvoltN;

    PetscReal Fphph0x[Nc], Fphph1x[Nc];
    PetscReal Fphph0y[Nc], Fphph1y[Nc];

    //Ionic concentration equations
    for (x = 0; x < Nx; x++) {
        for (y = 0; y < Ny; y++) {
            for (comp = 0; comp < Nc; comp++) {
                Fphph0x[comp] = 0;
                Fphph1x[comp] = 0;
                Fphph0y[comp] = 0;
                Fphph1y[comp] = 0;
            }
            for (ion = 0; ion < Ni; ion++) {
                for (comp = 0; comp < Nc - 1; comp++) {
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
                    if (x < Nx - 1) {
                        Ftmpx = Dcs[c_index(x, y, comp, ion, Nx) * 2] *
                                (cp[c_index(x, y, comp, ion, Nx)] + cp[c_index(x + 1, y, comp, ion, Nx)]) / 2 / dx *
                                dt / dx;
                        Fc0x = Ftmpx / c[c_index(x, y, comp, ion, Nx)];
                        Fph0x = z[ion] * Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx),
                                           -z[ion] * Fc0x, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (x > 0) {
                        Ftmpx = Dcs[c_index(x - 1, y, comp, ion, Nx) * 2] *
                                (cp[c_index(x - 1, y, comp, ion, Nx)] + cp[c_index(x, y, comp, ion, Nx)]) / 2 / dx *
                                dt / dx;
                        Fc1x = Ftmpx / c[c_index(x, y, comp, ion, Nx)];
                        Fph1x = z[ion] * Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx),
                                           -z[ion] * Fc1x, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y < Ny - 1) {
                        Ftmpy = Dcs[c_index(x, y, comp, ion, Nx) * 2 + 1] *
                                (cp[c_index(x, y, comp, ion, Nx)] + cp[c_index(x, y + 1, comp, ion, Nx)]) / 2 / dy *
                                dt / dy;
                        Fc0y = Ftmpy / c[c_index(x, y, comp, ion, Nx)];
                        Fph0y = z[ion] * Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx),
                                           -z[ion] * Fc0y, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y > 0) {
                        Ftmpy = Dcs[c_index(x, y - 1, comp, ion, Nx) * 2 + 1] *
                                (cp[c_index(x, y - 1, comp, ion, Nx)] + cp[c_index(x, y, comp, ion, Nx)]) / 2 / dy *
                                dt / dy;
                        Fc1y = Ftmpy / c[c_index(x, y, comp, ion, Nx)];
                        Fph1y = z[ion] * Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx),
                                           -z[ion] * Fc1y, INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x, y, comp, Nx)] + Fc0x + Fc1x + Fc0y + Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Add up terms for voltage eqns
                    Fphph0x[comp] += z[ion] * Fph0x;
                    Fphph1x[comp] += z[ion] * Fph1x;
                    Fphph0y[comp] += z[ion] * Fph0y;
                    Fphph1y[comp] += z[ion] * Fph1y;

                    //membrane current contributions
                    Ac += flux->dfdci[c_index(x, y, comp, ion, Nx)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, Nc - 1, Nx), Ind_1(x, y, ion, comp, Nx),
                                       -flux->dfdci[c_index(x, y, comp, ion, Nx)] * dt, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, comp, Nx), Ind_1(x, y, ion, Nc - 1, Nx),
                                       flux->dfdce[c_index(x, y, comp, ion, Nx)] * dt, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, Nc - 1, Nx), Ind_1(x, y, Ni, comp, Nx),
                                       -flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, comp, Nx), Ind_1(x, y, Ni, Nc - 1, Nx),
                                       -flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), Ac, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac, Ind_1(x, y, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), Aphi, INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx), z[ion] *
                                                                                                   (Fc0x + Fc1x + Fc0y +
                                                                                                    Fc1y +
                                                                                                    flux->dfdci[c_index(
                                                                                                            x, y, comp,
                                                                                                            ion, Nx)] *
                                                                                                    dt), INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //IntraPhi with c extra(volt eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, ion, Nc - 1, Nx),
                                       z[ion] * (flux->dfdce[c_index(x, y, comp, ion, Nx)] * dt), INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Extra-Phi with intra-c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1, Nx), Ind_1(x, y, ion, comp, Nx),
                                       -z[ion] * (flux->dfdci[c_index(x, y, comp, ion, Nx)] * dt), INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc - 1;
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
                if (x < Nx - 1) {
                    Ftmpx = Dcs[c_index(x, y, comp, ion, Nx) * 2] *
                            (cp[c_index(x, y, comp, ion, Nx)] + cp[c_index(x + 1, y, comp, ion, Nx)]) / 2 / dx * dt /
                            dx;
                    Fc0x = Ftmpx / c[c_index(x, y, comp, ion, Nx)];
                    Fph0x = z[ion] * Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Right Phi with left c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx), -z[ion] * Fc0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    Ftmpx = Dcs[c_index(x - 1, y, comp, ion, Nx) * 2] *
                            (cp[c_index(x - 1, y, comp, ion, Nx)] + cp[c_index(x, y, comp, ion, Nx)]) / 2 / dx * dt /
                            dx;
                    Fc1x = Ftmpx / c[c_index(x, y, comp, ion, Nx)];
                    Fph1x = z[ion] * Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // left Phi with right c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx), -z[ion] * Fc1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    Ftmpy = Dcs[c_index(x, y, comp, ion, Nx) * 2 + 1] *
                            (cp[c_index(x, y, comp, ion, Nx)] + cp[c_index(x, y + 1, comp, ion, Nx)]) / 2 / dy * dt /
                            dy;
                    Fc0y = Ftmpy / c[c_index(x, y, comp, ion, Nx)];
                    Fph0y = z[ion] * Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Upper Phi with lower c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx), -z[ion] * Fc0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    Ftmpy = Dcs[c_index(x, y - 1, comp, ion, Nx) * 2 + 1] *
                            (cp[c_index(x, y - 1, comp, ion, Nx)] + cp[c_index(x, y, comp, ion, Nx)]) / 2 / dy * dt /
                            dy;
                    Fc1y = Ftmpy / c[c_index(x, y, comp, ion, Nx)];
                    Fph1y = z[ion] * Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, ion, comp, Nx), Ind_1(x, y, ion, comp, Nx), -Fc1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, ion, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fph1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Lower Phi with upper c (voltage eqn)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp, Nx), Ind_1(x, y, ion, comp, Nx), -z[ion] * Fc1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1 - al[al_index(x, y, 0, Nx)] - al[al_index(x, y, 1, Nx)]) + Fc0x + Fc1x + Fc0y + Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                Avolt = z[ion] * (Fc0x + Fc1x + Fc0y + Fc1y);

                //Add up terms for voltage eqns
                Fphph0x[comp] += z[ion] * Fph0x;
                Fphph1x[comp] += z[ion] * Fph1x;
                Fphph0y[comp] += z[ion] * Fph0y;
                Fphph1y[comp] += z[ion] * Fph1y;

                //Membrane current contribution
                for (comp = 0; comp < Nc - 1; comp++) {
                    Ac -= flux->dfdce[c_index(x, y, comp, ion, Nx)] * dt;
                    Aphi += flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt;
                    Avolt -= z[ion] * flux->dfdce[c_index(x, y, comp, ion, Nx)] * dt;
                }
                //Add bath contributions
                Ftmpx = sqrt(pow(Dcb[c_index(x, y, Nc - 1, ion, Nx) * 2], 2) +
                             pow(Dcb[c_index(x, y, Nc - 1, ion, Nx) * 2 + 1], 2));
                Ac -= Ftmpx * (cp[c_index(x, y, Nc - 1, ion, Nx)] + cbath[ion]) /
                      (2 * c[c_index(x, y, Nc - 1, ion, Nx)]) * dt;
                Aphi -= Ftmpx * (cp[c_index(x, y, Nc - 1, ion, Nx)] + cbath[ion]) * z[ion] / 2 * dt;

                Avolt -= z[ion] * Ftmpx * (cp[c_index(x, y, Nc - 1, ion, Nx)] + cbath[ion]) /
                         (2 * c[c_index(x, y, Nc - 1, ion, Nx)]) * dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac, Ind_1(x, y, ion, Nc - 1, Nx), Ind_1(x, y, ion, Nc - 1, Nx), Ac, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac, Ind_1(x, y, ion, Nc - 1, Nx), Ind_1(x, y, Ni, Nc - 1, Nx), Aphi, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, Nc - 1, Nx), Ind_1(x, y, ion, Nc - 1, Nx), Avolt,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for (comp = 0; comp < Nc - 1; comp++) {
                if (x < Nx - 1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph0x[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph1x[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph0y[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph1y[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                Avolt = cm[comp] + Fphph0x[comp] + Fphph1x[comp] + Fphph0y[comp] + Fphph1y[comp];
                AvoltN = -cm[comp];
                for (ion = 0; ion < Ni; ion++) {
                    Avolt += z[ion] * flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt;
                    AvoltN -= z[ion] * flux->dfdphim[c_index(x, y, comp, ion, Nx)] * dt;
                }

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), Avolt, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, Ni, Nc - 1, Nx), AvoltN, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc - 1;
            if (x < Nx - 1) {
                //Right phi with left phi (-Fph0x)
                ierr = MatSetValue(Jac, Ind_1(x + 1, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph0x[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (x > 0) {
                //Left phi with right phi (-Fph1x)
                ierr = MatSetValue(Jac, Ind_1(x - 1, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph1x[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (y < Ny - 1) {
                //Upper phi with lower phi (-Fph0y)
                ierr = MatSetValue(Jac, Ind_1(x, y + 1, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph0y[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (y > 0) {
                //Lower phi with upper phi (-Fph1y)
                ierr = MatSetValue(Jac, Ind_1(x, y - 1, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), -Fphph1y[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            AvoltN = 0;

            for (int k = 0; k < Nc - 1; k++) {
                AvoltN += cm[k];
                Avolt = -cm[k];
                for (ion = 0; ion < Ni; ion++) {
                    Avolt -= z[ion] * flux->dfdphim[c_index(x, y, k, ion, Nx)] * dt;
                    AvoltN += z[ion] * flux->dfdphim[c_index(x, y, k, ion, Nx)] * dt;
                }
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, Ni, k, Nx), Avolt, INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }

            AvoltN += Fphph0x[comp] + Fphph1x[comp] + Fphph0y[comp] + Fphph1y[comp];

            //Bath terms
            for (ion = 0; ion < Ni; ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x, y, Nc - 1, ion, Nx) * 2], 2) +
                             pow(Dcb[c_index(x, y, Nc - 1, ion, Nx) * 2 + 1], 2));
                AvoltN -= z[ion] * Ftmpx * (cp[c_index(x, y, Nc - 1, ion, Nx)] + cbath[ion]) * z[ion] / 2 * dt;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac, Ind_1(x, y, Ni, comp, Nx), Ind_1(x, y, Ni, comp, Nx), AvoltN, INSERT_VALUES);
            CHKERRQ(ierr);
            ind++;

        }
    }

    ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    if (Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;

}

PetscErrorCode Newton_Solve_Grid(PetscInt xi, PetscInt yi,struct AppCtx *user) {


    PetscReal rsd;
    PetscErrorCode ierr = 0;
    PetscReal *temp;
    PetscInt num_iter;
    PetscReal rnorm;

    PetscInt x,y,comp,ion;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;


    PetscReal tol = reltol * array_max(user->grid_vars_past->c, (size_t) Nx * Ny * Ni * Nc);
    tol = reltol * tol;
    rsd = tol + 1;

    for (PetscInt iter = 0; iter < itermax; iter++) {

        ierr = Grid_Residual(user->slvr->Res, xi, yi, user);CHKERRQ(ierr);

        ierr = VecNorm(user->slvr->Res, NORM_MAX, &rsd);CHKERRQ(ierr);
        printf("Iteration: %d, Residual: %.10e\n", iter, rsd);
        if (rsd < tol) {
            if (details) {
                printf("Iteration: %d, Residual: %.10e\n", iter, rsd);
            }
            return ierr;
        }
        ierr = Grid_Jacobian(user->slvr->A, xi, yi, user);CHKERRQ(ierr);

        //Set the new operator
        ierr = KSPSetOperators(user->slvr->ksp, user->slvr->A, user->slvr->A);CHKERRQ(ierr);

        //Solve
        ierr = KSPSolve(user->slvr->ksp, user->slvr->Res, user->slvr->Q);CHKERRQ(ierr);

        ierr = KSPGetIterationNumber(user->slvr->ksp, &num_iter);CHKERRQ(ierr);
        ierr = KSPGetResidualNorm(user->slvr->ksp, &rnorm);CHKERRQ(ierr);

        // printf("KSP Solve time: %f, iter num:%d, norm: %.10e\n",toc-tic,num_iter,rnorm);
//        ierr = KSPView(slvr->ksp,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);


        if (details) {
            printf("iter num:%d, norm: %.10e\n", num_iter, rnorm);
        }

//        PetscTime(&tic);
        ierr = VecGetArray(user->slvr->Q, &temp);
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        user->grid_vars->c[c_index(x, y, comp, ion, Nx)] -= temp[Ind_1(x,y,ion,comp,Nx)];
                    }
                    user->grid_vars->phi[phi_index(x, y, comp, Nx)] -= temp[Ind_1(x,y,Ni,comp,Nx)];
                }
            }
        }
        ierr = VecRestoreArray(user->slvr->Q, &temp);


        if (details) {
            printf("Iteration: %d, Residual: %.10e\n", iter, rsd);
        }
    }

    if (rsd > tol) {
        fprintf(stderr, "Netwon Iteration did not converge! Stopping...\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }
    return ierr;
}
PetscErrorCode Update_Grid(PetscInt xi, PetscInt yi,PetscReal t,struct AppCtx *user)
{
    PetscErrorCode ierr;


    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscInt ion,comp,x,y;

    //Load current variable into past variable
    for( x=0;x<Nx;x++) {
        for ( y = 0; y < Ny; y++) {

            for (comp = 0; comp < Nc; comp++) {
                for (ion = 0; ion < Ni; ion++) {
                    user->grid_vars_past->c[c_index(x, y, comp, ion, Nx)] = user->grid_vars->c[c_index(x, y, comp, ion, Nx)] ;
                }
                user->grid_vars_past->phi[phi_index(x, y, comp, Nx)] = user->grid_vars->phi[phi_index(x, y, comp, Nx)];
            }
            for (comp = 0; comp < Nc - 1; comp++) {
                user->grid_vars_past->alpha[al_index(x, y, comp, Nx)] = user->grid_vars->alpha[al_index(x, y, comp, Nx)];
            }

        }
    }

    //Perform Newton Solve
    ierr = Newton_Solve_Grid(xi,yi,user);

    //Update Gating variable
    gatevars_update_grid(user->grid_gate_vars,user->grid_vars,user->dt*1e3,user);

    //Update Excitation
    excitation_grid(user,t-user->dt,xi,yi);


    return ierr;
}

PetscErrorCode Update_Solution(Vec current_state,PetscReal t,struct AppCtx *user)
{

    PetscErrorCode ierr = 0;
    PetscInt x,y,ion,comp;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    extract_subarray(current_state,user->state_vars);
    extract_subarray(user->state_vars_past->v,user->state_vars_past);

    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){

            //Load new gridpoint
            Load_Grid(user,x,y);
            //Update new grid
            Update_Grid(x,y,t,user);

            //Save the held variable
            for ( comp = 0; comp < Nc; comp++) {
                for ( ion = 0; ion < Ni; ion++) {
                    user->state_vars->c[c_index(x,y,comp,ion,Nx)]=user->grid_vars->c[c_index(x,y,comp,ion,Nx)];
                }
                user->state_vars->phi[phi_index(x,y,comp,Nx)]=user->grid_vars->phi[phi_index(x,y,comp,Nx)];
            }

        }
    }


    restore_subarray(current_state,user->state_vars);
    restore_subarray(user->state_vars_past->v,user->state_vars_past);

    return ierr;

}