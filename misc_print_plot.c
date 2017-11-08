#include "functions.h"
#include "constants.h"


void print_all(double *Dcs, double *Dcb, struct ConstVars *con_vars, struct FluxData *flux, struct GateType *gvars,
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
            printf("Ion: %d, Comp %d, C: %f\n",ion,comp,state_vars->c[c_index(0,0,comp,ion)]);
        }
    }
    for(PetscInt comp=0;comp<Nc;comp++)
    {
        printf("Comp %d, Phi: %f\n",comp,state_vars->phi[phi_index(0,0,comp)]);
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
            printf("Flux: %f*1e3, dfdci: %f, dfdce: %f, dfdphim: %f\n",1e3*flux->mflux[c_index(0,0,comp,ion)],flux->dfdci[c_index(0,0,comp,ion)],flux->dfdce[c_index(0,0,comp,ion)],flux->dfdphim[c_index(0,0,comp,ion)]);
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

