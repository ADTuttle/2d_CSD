#include "functions.h"
#include "constants.h"
#include <string.h>
#include <math.h>


void print_all(struct AppCtx *user)
{
    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct ConstVars *con_vars = user->con_vars;
    struct FluxData *flux = user->flux;
    struct GateType *gvars = user->gate_vars;
    struct SimState *state_vars = user->state_vars;
    struct Solver *slvr = user->slvr;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    printf("ConstVars:\n");
    printf("%f,%f,%f\n",1e6*con_vars->pNaKCl,1e6*con_vars->Imax,1e6*con_vars->pNaLeak);
    printf("%f,%f\n",1e6*con_vars->Imaxg,1e6*con_vars->pNaLeakg);
    printf("%f,%f,%f\n",con_vars->ao[0],con_vars->ao[1],con_vars->ao[2]);
    printf("%f,%f,%f\n",con_vars->zo[0],con_vars->zo[1],con_vars->zo[2]);
    printf("%f,%f,%f\n",con_vars->kappa,con_vars->zeta1[0],con_vars->zeta1[1]);
    printf("%d,%f,%f\n",con_vars->S,1e6*con_vars->zetaalpha[0],1e6*con_vars->zetaalpha[1]);

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at even positions (0,2,4,...)
    //y at odd (1,3,5,...)
    //still use c_index(x,y,comp,ion,Nx), but with ind*2 or ind*2+1

    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Dcs: Ion %d, Comp %d ",ion,comp);
            printf("Dcs x: %f, Dcs y: %f\n",1e6*Dcs[c_index(0,0,comp,ion,Nx)*2],1e6*Dcs[c_index(0,4,comp,ion,Nx)*2+1]);
        }
    }
    printf("\n");

    //Bath diffusion

    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Dcb: Ion %d, Comp %d ",ion,comp);
            printf("Dcb x: %f, Dcb y: %f\n",1e6*Dcb[c_index(0,0,comp,ion,Nx)*2],1e6*Dcb[c_index(0,4,comp,ion,Nx)*2+1]);
        }
    }

    PetscInt x=0;PetscInt y=0;
    printf("\n");
    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Ion: %d, Comp %d, C: %.10e\n",ion,comp,state_vars->c[c_index(0,0,comp,ion,Nx)]);
        }
    }
    for(PetscInt comp=0;comp<Nc;comp++)
    {
        printf("Comp %d, Phi: %.10e\n",comp,state_vars->phi[phi_index(0,0,comp,Nx)]);
    }
    for(PetscInt comp=0;comp<Nc-1;comp++)
    {
        printf("Comp %d, alpha: %.10e\n",comp,state_vars->alpha[al_index(0,0,comp,Nx)]);
    }
    printf("Gvars:\n");
    printf("NaT :%f,%f,%f*1e-6\n",gvars->mNaT[0],gvars->hNaT[0],1e6*gvars->gNaT[0]);
    printf("NaP :%f,%f,%f\n",gvars->mNaP[0],gvars->hNaP[0],gvars->gNaP[0]);
    printf("KDR :%f,%f\n",gvars->mKDR[0],gvars->gKDR[0]);
    printf("KA :%f,%f,%f\n",gvars->mKA[0],gvars->hKA[0],gvars->gKA[0]);
    printf("\n");
    //Compute membrane ionic flux relation quantitites

    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Ion: %d, Comp %d\n",ion,comp);
            printf("Flux: %f*1e20, dfdci: %f, dfdce: %f, dfdphim: %f\n",1e20*flux->mflux[c_index(0,0,comp,ion,Nx)],flux->dfdci[c_index(0,0,comp,ion,Nx)],flux->dfdce[c_index(0,0,comp,ion,Nx)],flux->dfdphim[c_index(0,0,comp,ion,Nx)]);
        }
    }
    printf("\n");
    //Compute membrane water flow related quantities
    for(PetscInt comp=0;comp<Nc-1;comp++)
    {
        printf("Comp: %d\n",comp);
        printf("wFlux: %f,%f,%f\n",flux->wflow[al_index(x,y,comp,Nx)],flux->dwdpi[al_index(x,y,comp,Nx)],flux->dwdal[al_index(x,y,comp,Nx)]);
    }
    printf("\n");
    // VecView(slvr->Res,PETSC_VIEWER_STDOUT_SELF);

    // VecView(slvr->Q,PETSC_VIEWER_STDOUT_SELF);

    return;

}

const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
         tok && *tok;
         tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

void find_print(int rowC, int colC, double valC, int iter)
{

    if(iter!=-1){
        return;
    }
//    printf("%d,%d,%.10e\n",rowC,colC,valC);
//    return;
    if(rowC>70){
        return;
    }

    int row,col;
    double val;

    //Open Julia Mat file.
    FILE *fp;
    fp = fopen("../../../Julia_work/2d_modules/matrix.txt","r");
    if(fp==NULL)
    {
        fprintf(stderr, "File not found....\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

    //Begin searching for the right row and col.
    char line[1024];
    while (fgets(line, 1024, fp))
    {
        char* tmp = strdup(line);
        row = -1; col=-1; val=0;

        row = atoi(getfield(tmp,1));


        tmp = strdup(line);
        col = atoi(getfield(tmp,2));

        if(row==rowC && col==colC)
        {
            tmp = strdup(line);
            val = atof(getfield(tmp,3));
            if(val!=0) {
                if (fabs(val - valC)/fabs(val) > 5e-16) {
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16\n", row, col, val, valC,
                           1e16 * (val - valC) / val);
                }
            }
            /*
            if(val==0){
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16\n",row,col,val,valC,1e16*(val-valC));

            }
            else{
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16,Abs: %f\n",row,col,val,valC,1e16*(val-valC)/val,1e16*(val-valC));
            }
             */
            return;
        }
        free(tmp);
    }


    fclose(fp);
    return;

}

void compare_res(double *Res, int iter)
{
    if(iter!=-1){
        return;
    }

    FILE *fp;
    fp = fopen("../../../Julia_work/2d_modules/Res.txt","r");
    if(fp==NULL)
    {
        fprintf(stderr, "File not found....\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

    //Begin searching for the right row and col.
    char line[1024];
    int row=0;
    double val;
    while (fgets(line, 1024, fp))
    {
        char* tmp = strdup(line);
        val=0;
        val = atof(getfield(tmp,1));

        if(val==0){
            printf("Row: %d, J: %.10e, C: %.10e, diff: %.10e\n",row,val,Res[row],val-Res[row]);
        }else{
            printf("Row: %d, J: %.10e, C: %.10e,abs: %.10e, diff: %.10e\n",row,val,Res[row],val-Res[row],(val-Res[row])/val);
        }
        row++;
        if(row>70){
            break;
        }
        free(tmp);
    }


    fclose(fp);
    return;
}

void write_data(FILE *fp,struct AppCtx*user,PetscInt numrecords,int start)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[8], 0, 0, 0, 0);
    }
    struct SimState *state_vars = user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    if(!save_one_var) {
        if (start) {
            fprintf(fp, "%d,%d,%d,%d,%d\n", Nx, Ny, numrecords, Nc, Ni);
            write_data(fp, user, numrecords,0);
        } else {
            int ion, comp, x, y;
            for (ion = 0; ion < Ni; ion++) {
                for (comp = 0; comp < Nc; comp++) {
                    for (y = 0; y < Ny; y++) {
                        for (x = 0; x < Nx; x++) {
                            if (x == Nx - 1 & y == Ny - 1) {
                                fprintf(fp, "%f\n", state_vars->c[c_index(x, y, comp, ion,Nx)]);
                            } else {
                                fprintf(fp, "%f,", state_vars->c[c_index(x, y, comp, ion,Nx)]);
                            }
                        }
                    }
                }
            }
            for (comp = 0; comp < Nc; comp++) {
                for (y = 0; y < Ny; y++) {
                    for (x = 0; x < Nx; x++) {
                        if (x == Nx - 1 & y == Ny - 1) {
                            fprintf(fp, "%f\n", state_vars->phi[phi_index(x, y, comp,Nx)] * RTFC);
                        } else {
                            fprintf(fp, "%f,", state_vars->phi[phi_index(x, y, comp,Nx)] * RTFC);
                        }
                    }
                }
            }
            for (comp = 0; comp < Nc - 1; comp++) {
                for (y = 0; y < Ny; y++) {
                    for (x = 0; x < Nx; x++) {
                        if (x == Nx - 1 & y == Ny - 1) {
                            fprintf(fp, "%f\n", state_vars->alpha[al_index(x, y, comp,Nx)]);
                        } else {
                            fprintf(fp, "%f,", state_vars->alpha[al_index(x, y, comp,Nx)]);
                        }
                    }
                }
            }
        }
    } else{
        if (start) {
            fprintf(fp, "%d,%d,%d,%d,%d\n", Nx, Ny, (int) floor(numrecords), 0, 0);
            write_data(fp, user,numrecords, 0);
        } else {
            int ion, comp, x, y;
            comp = 0;
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    if (x == Nx - 1 & y == Ny - 1) {
//                        fprintf(fp, "%f\n", state_vars->phi[phi_index(x, y, Nc-1,Nx)] * RTFC);
                        fprintf(fp, "%f\n", (state_vars->phi[phi_index(x, y, comp,Nx)]-state_vars->phi[phi_index(x, y, Nc-1,Nx)]) * RTFC);
                    } else {
//                        fprintf(fp, "%f,", state_vars->phi[phi_index(x, y, Nc-1,Nx)] * RTFC);
                        fprintf(fp, "%f,", (state_vars->phi[phi_index(x, y, comp,Nx)]-state_vars->phi[phi_index(x, y, Nc-1,Nx)]) * RTFC);
                    }
                }
            }
        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[8], 0, 0, 0, 0);
    }
}
void write_point(FILE *fp,struct AppCtx* user,PetscInt numrecords,int start)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[8], 0, 0, 0, 0);
    }
    struct SimState *state_vars = user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    if (start) {
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d\n", Nx, Ny, numrecords, Nc, Ni,use_en_deriv,separate_vol,Linear_Diffusion);
        write_point(fp, user,numrecords, 0);
        } else {
            int ion, comp;
            int x =10;
            int y=10;
            for (ion = 0; ion < Ni; ion++) {
                for (comp = 0; comp < Nc; comp++) {
//                    fprintf(fp, "%f,", state_vars->c[c_index(x, y, comp, ion,Nx)]);
                    fprintf(fp, "%.10e,", state_vars->c[c_index(x, y, comp, ion,Nx)]);
                }
            }

            for (comp = 0; comp < Nc; comp++) {
//                fprintf(fp, "%f,", state_vars->phi[phi_index(x, y, comp,Nx)] * RTFC);
                fprintf(fp, "%.10e,", state_vars->phi[phi_index(x, y, comp,Nx)] * RTFC);
            }
            for (comp = 0; comp < Nc - 1; comp++) {
//                fprintf(fp, "%f,", state_vars->alpha[al_index(x, y, comp,Nx)]);
                fprintf(fp, "%.10e,", state_vars->alpha[al_index(x, y, comp,Nx)]);
            }
            fprintf(fp,"\n");
        }


}

void measure_flux(FILE *fp, struct AppCtx* user,PetscInt numrecords,int start)
{
    struct SimState *state_vars= user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    if (start) {
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d\n", Nx, Ny, numrecords, Nc, Ni, use_en_deriv, separate_vol,
                Linear_Diffusion);
        measure_flux(fp, user, numrecords, 0);
    } else{
        //compute diffusion coefficients
        diff_coef(user->Dcs,state_vars->alpha,1,user);
        //Bath diffusion
        diff_coef(user->Dcb,state_vars->alpha,Batheps,user);

        PetscReal *c = state_vars->c;
        PetscReal *phi = state_vars->phi;
        PetscReal *al = state_vars->alpha;
        PetscReal *Dcs = user->Dcs;
        PetscReal *Dcb = user->Dcb;
        PetscReal dx = user->dx;
        PetscReal dy = user->dy;

        PetscReal Ftmp;

        PetscReal Rcvx,RcvxRight,alNc;
        PetscReal Rcvy,RcvyUp;

        PetscReal Rphx,RphxRight;
        PetscReal Rphy,RphyUp;

        PetscReal RBath;

        PetscInt x,y,comp,ion;

        PetscReal Circ_Radius = 0.05;
        PetscReal radius;

        PetscReal Fluxc[Nc] = {0,0,0};
        PetscReal Fluxph[Nc] = {0,0,0};
        PetscReal Fluxbath[Nc] = {0,0,0};

        int num_points=0;

        for(x=0;x<Nx;x++){
            for(y=0;y<Ny;y++){
                radius = sqrt(pow((x + 0.5) * dx - Lx / 2, 2) + pow((y + 0.5) * dy - Lx / 2, 2));
                if(radius<=Circ_Radius) {
//                if(fabs((x+0.5)*dx-Lx/2)<=0.1 && fabs((y+0.5)*dy-Ly/2)<=0.1 ){
                    num_points++;
                    for(comp=0;comp<Nc;comp++) {
                        //Reset values
                        Rphx = 0;
                        RphxRight = 0;
                        Rphy = 0;
                        RphyUp = 0;
                        Rcvx = 0;
                        RcvxRight = 0;
                        Rcvy = 0;
                        RcvyUp = 0;
                        RBath = 0;
                        //Sum over all ions
                        for (ion = 0; ion < Ni; ion++) {
                            if (x > 0) {
                                //First difference term
                                Ftmp = z[ion] * Dcs[c_index(x - 1, y, comp, ion, Nx) * 2] / (dx * dx);
                                Rcvx += Ftmp * (c[c_index(x, y, comp, ion, Nx)] - c[c_index(x - 1, y, comp, ion, Nx)]);
                                Rphx += Ftmp * z[ion] * c[c_index(x - 1, y, comp, ion, Nx)] *
                                        (phi[phi_index(x, y, comp, Nx)] - phi[phi_index(x - 1, y, comp, Nx)]);
                            }
                            //Add Second right moving difference
                            if (x < Nx - 1) {
                                //Second difference term
                                Ftmp = z[ion] * Dcs[c_index(x, y, comp, ion, Nx) * 2] / (dx * dx);
                                RcvxRight +=
                                        Ftmp * (c[c_index(x + 1, y, comp, ion, Nx)] - c[c_index(x, y, comp, ion, Nx)]);
                                RphxRight += Ftmp * z[ion] * c[c_index(x, y, comp, ion, Nx)] *
                                             (phi[phi_index(x + 1, y, comp, Nx)] - phi[phi_index(x, y, comp, Nx)]);
                            }
                            if (y > 0) {
                                //Updown difference term
                                Ftmp = z[ion] * Dcs[c_index(x, y - 1, comp, ion, Nx) * 2 + 1] / (dy * dy);
                                Rcvy += Ftmp * (c[c_index(x, y, comp, ion, Nx)] - c[c_index(x, y - 1, comp, ion, Nx)]);
                                Rphy += Ftmp * z[ion] * c[c_index(x, y - 1, comp, ion, Nx)] *
                                        (phi[phi_index(x, y, comp, Nx)] - phi[phi_index(x, y - 1, comp, Nx)]);
                            }
                            //Next upward difference
                            if (y < Ny - 1) {
                                Ftmp = z[ion] * Dcs[c_index(x, y, comp, ion, Nx) * 2 + 1] / (dy * dy);
                                RcvyUp +=
                                        Ftmp * (c[c_index(x, y + 1, comp, ion, Nx)] - c[c_index(x, y, comp, ion, Nx)]);
                                RcvyUp += Ftmp * z[ion] * c[c_index(x, y, comp, ion, Nx)] *
                                          (phi[phi_index(x, y + 1, comp, Nx)] - phi[phi_index(x, y, comp, Nx)]);

                            }
                            if (comp == Nc - 1) {
                            Ftmp = z[ion] * sqrt(pow(Dcb[c_index(x, y, comp, ion, Nx) * 2], 2) +
                                                 pow(Dcb[c_index(x, y, comp, ion, Nx) * 2 + 1], 2));
                            RBath -= Ftmp * (c[c_index(x, y, comp, ion, Nx)] - cbath[ion]);
                            RBath -= Ftmp * (c[c_index(x, y, comp, ion, Nx)] + cbath[ion]) / 2.0 *
                                     (z[ion] * phi[phi_index(x, y, comp, Nx)] - z[ion] * phibath);
                            }
                        }

                        Fluxph[comp] += (Rphx - RphxRight + Rphy - RphyUp) * dx * dy;
                        Fluxc[comp] += (Rcvx - RcvxRight + Rcvy - RcvyUp) * dx * dy;
                        Fluxbath[comp] += RBath * dx * dy;
                    }
                }
            }
        }
/*
        for(comp=Nc-1;comp<Nc;comp++) {
            printf("Comp: %d\n",comp);
            printf("pts: %d,FluxC: %.10e, Fluxph: %.10e,FluxBath: %.10e\n", num_points, Fluxc[comp], Fluxph[comp], Fluxbath[comp]);
            printf("Flux nobath: %.10e, Flux Tot: %.10e\n", Fluxc[comp] + Fluxph[comp], Fluxbath[comp] + Fluxc[comp] + Fluxph[comp]);
        }
        */
        fprintf(fp,"%.10e,%.10e,%.10e\n",Fluxc[Nc-1], Fluxph[Nc-1], Fluxbath[Nc-1]);
        /*
        printf("All Comps:\n");
        printf("FluxC: %.10e, Fluxph: %.10e,FluxBath: %.10e\n", Fluxc[0]+Fluxc[1]+Fluxc[2], Fluxph[0]+Fluxph[1]+Fluxph[2], Fluxbath[2]);
        printf("Flux nobath: %.10e, Flux Tot: %.10e\n", Fluxc[0]+Fluxc[1]+Fluxc[2]+Fluxph[0]+Fluxph[1]+Fluxph[2], Fluxbath[2]+Fluxc[0]+Fluxc[1]+Fluxc[2]+Fluxph[0]+Fluxph[1]+Fluxph[2]);
    */
    }
}
void init_events(struct AppCtx *user)
{
    PetscLogDefaultBegin();


    PetscClassId id;
    PetscClassIdRegister("CSD",&id);

    PetscLogEventRegister("Jacobian",id,&event[0]);
    PetscLogEventRegister("Residual",id,&event[1]);
    PetscLogEventRegister("Extract Subarray",id,&event[2]);
    PetscLogEventRegister("Restore Subarray",id,&event[3]);
    PetscLogEventRegister("Gating Variables",id,&event[4]);
    PetscLogEventRegister("Ion Flux",id,&event[5]);
    PetscLogEventRegister("Water Flux",id,&event[6]);
    PetscLogEventRegister("Volume Update",id,&event[7]);
    PetscLogEventRegister("Write to File",id,&event[8]);

    //Deactivate Petsc tracking
    PetscLogEvent deactivate;

    PetscLogEventDeactivateClass(MAT_CLASSID);
    PetscLogEventDeactivateClass(KSP_CLASSID); /* includes PC and KSP */
    PetscLogEventDeactivateClass(VEC_CLASSID);
//    PetscLogEventDeactivateClass(SNES_CLASSID);

    //Some events are leftover somehow

    PetscLogEventGetId("PCApply",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("VecSet",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("MatAssemblyBegin",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("MatAssemblyEnd",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("SNESLineSearch",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("PCSetUp",&deactivate);
    PetscLogEventDeactivate(deactivate);

    PetscLogEventGetId("SNESFunctionEval",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("SNESJacobianEval",&deactivate);
    PetscLogEventDeactivate(deactivate);

}