#include "functions.h"
#include <string.h>
#include <math.h>


void print_all(PetscReal *Dcs, PetscReal *Dcb, struct ConstVars *con_vars, struct FluxData *flux, struct GateType *gvars,
               struct SimState *state_vars, struct Solver *slvr)
{
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
    //still use c_index(x,y,comp,ion), but with ind*2 or ind*2+1

    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Dcs: Ion %d, Comp %d ",ion,comp);
            printf("Dcs x: %f, Dcs y: %f\n",1e6*Dcs[c_index(0,0,comp,ion)*2],1e6*Dcs[c_index(0,4,comp,ion)*2+1]);
        }
    }
    printf("\n");

    //Bath diffusion

    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Dcb: Ion %d, Comp %d ",ion,comp);
            printf("Dcb x: %f, Dcb y: %f\n",1e6*Dcb[c_index(0,0,comp,ion)*2],1e6*Dcb[c_index(0,4,comp,ion)*2+1]);
        }
    }

    PetscInt x=0;PetscInt y=0;
    printf("\n");
    for(PetscInt ion=0;ion<Ni;ion++)
    {
        for(PetscInt comp=0;comp<Nc;comp++)
        {
            printf("Ion: %d, Comp %d, C: %.10e\n",ion,comp,state_vars->c[c_index(0,0,comp,ion)]);
        }
    }
    for(PetscInt comp=0;comp<Nc;comp++)
    {
        printf("Comp %d, Phi: %.10e\n",comp,state_vars->phi[phi_index(0,0,comp)]);
    }
    for(PetscInt comp=0;comp<Nc-1;comp++)
    {
        printf("Comp %d, alpha: %.10e\n",comp,state_vars->alpha[al_index(0,0,comp)]);
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
            printf("Flux: %f*1e20, dfdci: %f, dfdce: %f, dfdphim: %f\n",1e20*flux->mflux[c_index(0,0,comp,ion)],flux->dfdci[c_index(0,0,comp,ion)],flux->dfdce[c_index(0,0,comp,ion)],flux->dfdphim[c_index(0,0,comp,ion)]);
        }
    }
    printf("\n");
    //Compute membrane water flow related quantities
    for(PetscInt comp=0;comp<Nc-1;comp++)
    {
        printf("Comp: %d\n",comp);
        printf("wFlux: %f,%f,%f\n",flux->wflow[al_index(x,y,comp)],flux->dwdpi[al_index(x,y,comp)],flux->dwdal[al_index(x,y,comp)]);
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

void write_data(FILE *fp,struct SimState *state_vars,int start)
{
    PetscLogEventBegin(event[8],0,0,0,0);
    if(start) {

        fprintf(fp,"%d,%d,%d,%d,%d\n",Nx,Ny,(int)floor(numrecords),Nc,Ni);
        write_data(fp,state_vars,0);
    }else {
        int ion, comp, x, y;
        for (ion = 0; ion < Ni; ion++) {
            for (comp = 0; comp < Nc; comp++) {
                for (y = 0; y < Ny; y++) {
                    for (x = 0; x < Nx; x++) {
                        if(x==Nx-1 & y==Ny-1){
                            fprintf(fp, "%f\n", state_vars->c[c_index(x, y, comp, ion)]);
                        }else {
                            fprintf(fp, "%f,", state_vars->c[c_index(x, y, comp, ion)]);
                        }
                    }
                }
            }
        }
        for (comp = 0; comp < Nc; comp++) {
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    if(x==Nx-1 & y==Ny-1){
                        fprintf(fp, "%f\n", state_vars->phi[phi_index(x, y, comp)]);
                    } else {
                        fprintf(fp, "%f,", state_vars->phi[phi_index(x, y, comp)]);
                    }
                }
            }
        }
        for (comp = 0; comp < Nc - 1; comp++) {
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    if(x==Nx-1 & y==Ny-1){
                        fprintf(fp, "%f\n", state_vars->alpha[al_index(x, y, comp)]);
                    } else {
                        fprintf(fp, "%f,", state_vars->alpha[al_index(x, y, comp)]);
                    }
                }
            }
        }
    }
    PetscLogEventEnd(event[8],0,0,0,0);
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