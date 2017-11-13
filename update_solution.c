#include "constants.h"
#include "functions.h"


PetscErrorCode newton_solve(struct SimState *state_vars,struct SimState *state_vars_past, PetscReal dt, struct GateType *gvars, struct ExctType *gexct, struct ConstVars *con_vars,struct Solver *slvr,struct FluxData *flux)
{

    PetscReal rsd;
    PetscErrorCode ierr = 0;
    PetscReal *temp;
    PetscInt num_iter;

    PetscInt index[NA];
    PetscReal before[NA];
    PetscReal after[NA];
    for(int i=0;i<NA;i++)
    {
        index[i]=i;
    }

    //Save the "current" aka past state
    memcpy(state_vars_past,state_vars,sizeof(struct SimState));

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at even positions (0,2,4,...)
    //y at odd (1,3,5,...)
    //still use c_index(x,y,comp,ion), but with ind*2 or ind*2+1
   	PetscReal Dcs[Nx*Ny*Ni*Nc*2];
   	PetscReal Dcb[Nx*Ny*Ni*Nc*2];
   	//compute diffusion coefficients
   	diff_coef(Dcs,state_vars->alpha,1);
   	//Bath diffusion
  	diff_coef(Dcb,state_vars->alpha,Batheps);

    PetscReal tol = reltol*array_max(state_vars->c,(size_t)Nx*Ny*Ni*Nc);
    rsd = tol+1;

    PetscInt x=0;PetscInt y=0;
    for(PetscInt iter=0;iter<itermax;iter++)
    {

        //Compute membrane ionic flux relation quantitites
        ionmflux(flux,state_vars,state_vars_past,gvars,gexct,con_vars);

        //Compute membrane water flow related quantities
        wflowm(flux,state_vars,con_vars);

        //Optional Print a bunch of info statement
        printf("\n\n Iteration: %d\n\n",iter);
        print_all(Dcs,Dcb,con_vars,flux,gvars,state_vars,slvr);

        // VecView(slvr->Res,PETSC_VIEWER_STDOUT_SELF);

//        if(iter==1){
//            printf("\n\n\n");
//            ierr = VecGetValues(slvr->Res,NA,index,before); CHKERRQ(ierr);
//        }

        PetscReal zero = 0;
        ierr = VecSet(slvr->Res,zero);CHKERRQ(ierr);
        ierr = VecSet(slvr->Q,zero);CHKERRQ(ierr);
        ierr = calc_residual(slvr->Res,state_vars_past,state_vars,dt,Dcs,Dcb,flux,con_vars);CHKERRQ(ierr);
//        if(iter==1){
//            ierr = VecGetValues(slvr->Res,NA,index,after); CHKERRQ(ierr);
//            compare_res(after,iter);
//        }

        ierr = VecNorm(slvr->Res,NORM_MAX,&rsd);CHKERRQ(ierr);
        PetscReal rsd2;
        ierr = VecNorm(slvr->Res,NORM_2,&rsd2);CHKERRQ(ierr);
//        printf("Rsd: %.10e,%.10e\n",rsd,rsd2);
        if(rsd<tol)
        {
            if(details)
            {
                printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
            }
            return ierr;
        }
        ierr = calc_jacobian(slvr->A, state_vars_past, state_vars, dt, Dcs, Dcb, flux, con_vars, iter);CHKERRQ(ierr);
        //Set the new operator
        // ierr = KSPSetOperators(slvr->ksp,NULL,NULL);CHKERRQ(ierr);
        ierr = KSPSetOperators(slvr->ksp,slvr->A,slvr->A);CHKERRQ(ierr);
        ierr = PCSetOperators(slvr->pc,slvr->A,slvr->A);CHKERRQ(ierr);
        ierr = KSPSetPC(slvr->ksp,slvr->pc);CHKERRQ(ierr);

        //Solve
//        ierr = KSPSetUp(slvr->ksp);CHKERRQ(ierr);
        ierr = KSPSolve(slvr->ksp,slvr->Res,slvr->Q);CHKERRQ(ierr);
//        ierr = KSPView(slvr->ksp,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = KSPGetIterationNumber(slvr->ksp,&num_iter); CHKERRQ(ierr);


//        if(iter==0) {
//            VecView(slvr->Res, PETSC_VIEWER_STDOUT_SELF);
//            VecView(slvr->Q, PETSC_VIEWER_STDOUT_SELF);
//        }
        printf("Number of KSP Iterations: %d\n",num_iter);
//         fprintf(stderr, "Performed first solve....\n");
//         exit(EXIT_FAILURE); /* indicate failure.*/
        ierr = VecGetArray(slvr->Q,&temp);CHKERRQ(ierr);
        for(x=0;x<Nx;x++)
        {
            for(y=0;y<Ny;y++)
            {
                for(PetscInt comp=0;comp<Nc;comp++)
                {
                    for(PetscInt ion = 0;ion<Ni;ion++)
                    {
                        // ierr = VecGetValues(slvr->Q,1,Ind_1(x,y,ion,comp),&temp); CHKERRQ(ierr);
                        state_vars->c[c_index(x,y,comp,ion)]-=temp[Ind_1(x,y,ion,comp)];
                    }
                    // ierr = VecGetValues(slvr->Q,1,Ind_1(x,y,Ni,comp),&temp); CHKERRQ(ierr);
                    state_vars->phi[phi_index(x,y,comp)]-=temp[Ind_1(x,y,Ni,comp)];
                }
                for(PetscInt comp=0;comp<Nc-1;comp++)
                {
                    // ierr = VecGetValues(slvr->Q,1,Ind_1(x,y,Ni+1,comp),&temp); CHKERRQ(ierr);
                    state_vars->alpha[al_index(x,y,comp)]-=temp[Ind_1(x,y,Ni+1,comp)];
                }
            }
        }
        ierr = VecRestoreArray(slvr->Q,&temp);CHKERRQ(ierr);


        if(details)
        {
            printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
        }
    }

    if(rsd>tol)
    {
        fprintf(stderr, "Netwon Iteration did not converge! Stopping...\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }
    return ierr;
}


PetscErrorCode calc_residual(Vec Res,struct SimState *state_vars_past,struct SimState *state_vars,PetscReal dt,PetscReal *Dcs,PetscReal *Dcb,struct FluxData *flux,struct ConstVars *con_vars)
{
    PetscReal *c = state_vars->c;
    PetscReal *phi = state_vars->phi;
    PetscReal *al = state_vars->alpha;
    PetscReal *cp = state_vars_past->c;
    PetscReal *alp = state_vars_past->alpha; 
    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;

    PetscErrorCode ierr;

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0)
                    {
                    //First difference term
                    Rcvx = Dcs[c_index(x-1,y,comp,ion)*2]*(cp[c_index(x-1,y,comp,ion)]+cp[c_index(x,y,comp,ion)])/2;
                    Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x-1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x-1,y,comp)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1)
                    {
                        RcvxRight = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x,y,comp)]))/dx*dt/dx;
                    } 
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0)
                    { 
                        Rcvy = Dcs[c_index(x,y-1,comp,ion)*2+1]*(cp[c_index(x,y-1,comp,ion)]+cp[c_index(x,y,comp,ion)])/2;
                        Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y-1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y-1,comp)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1)
                    {
                        RcvyUp = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y,comp)]))/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)]-alp[al_index(x,y,comp)]*cp[c_index(x,y,comp,ion)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;

                    ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
                alpNc = 1 - alp[al_index(x,y,0)] - alp[al_index(x,y,1)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0)
                {
                //First difference term
                Rcvx = Dcs[c_index(x-1,y,comp,ion)*2]*(cp[c_index(x-1,y,comp,ion)]+cp[c_index(x,y,comp,ion)])/2;
                Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x-1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x-1,y,comp)]))/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1)
                {
                    RcvxRight = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x,y,comp)]))/dx*dt/dx;
                } 
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0)
                { 
                    Rcvy = Dcs[c_index(x,y-1,comp,ion)*2+1]*(cp[c_index(x,y-1,comp,ion)]+cp[c_index(x,y,comp,ion)])/2;
                    Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y-1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y-1,comp)]))/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1)
                {
                    RcvyUp = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y,comp)]))/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,comp,ion)]-alpNc*cp[c_index(x,y,comp,ion)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,comp,ion)*2],2)+pow(Dcb[c_index(x,y,comp,ion)*2+1],2))*(cp[c_index(x,y,comp,ion)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp)]-z[ion]*phibath)*dt;
                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
    }
    

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {    
            //Residual for electroneutrality condition
            for(comp=0;comp<Nc-1;comp++)
            {

                Resc = al[al_index(x,y,comp)]*cz(c,z,x,y,comp)+con_vars->zo[phi_index(0,0,comp)]*con_vars->ao[phi_index(0,0,comp)];
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),Resc,INSERT_VALUES); CHKERRQ(ierr);
            }
            //Extracellular term
            comp=Nc-1;
            Resc = (1-al[al_index(x,y,0)]-al[al_index(x,y,1)])*cz(c,z,x,y,comp)+con_vars->zo[phi_index(0,0,comp)]*con_vars->ao[phi_index(0,0,comp)];
            ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),Resc,INSERT_VALUES); CHKERRQ(ierr);
            
            //Residual for water flow
            //Plus modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++)
            {
                //Water flow
                ierr = VecSetValue(Res,Ind_1(x,y,Ni+1,comp),al[al_index(x,y,comp)]-alp[al_index(x,y,comp)]+flux->wflow[al_index(x,y,comp)]*dt,INSERT_VALUES);CHKERRQ(ierr);

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
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,Nc-1),-cm[comp]*(phi[phi_index(x,y,Nc-1)]-phi[phi_index(x,y,comp)]),ADD_VALUES);CHKERRQ(ierr);
                //Intracell voltage mod
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),-cm[comp]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y,Nc-1)]),ADD_VALUES);CHKERRQ(ierr);
            }
        }
    }

    
//    for(int i=0;i<NA;i++)
//    {
//        ierr = VecSetValue(Res,i,1,INSERT_VALUES);CHKERRQ(ierr);
//    }
    
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    // VecView(Res,PETSC_VIEWER_STDOUT_SELF);


    return ierr;
}

PetscErrorCode
calc_jacobian(Mat Jac, struct SimState *state_vars_past, struct SimState *state_vars, PetscReal dt, PetscReal *Dcs,
              PetscReal *Dcb, struct FluxData *flux, struct ConstVars *con_vars, int iter)
{
    PetscErrorCode ierr;
    PetscReal *c = state_vars->c;
    PetscReal *al = state_vars->alpha;
    PetscReal *cp = state_vars_past->c;
    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi;

    //Ionic concentration equations
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
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
                    if(x<Nx-1)
                    {
                        Ftmpx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2/dx*dt/dx;
                        Fc0x = Ftmpx/c[c_index(x,y,comp,ion)];
                        Fph0x = z[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,ion,comp),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x + 1, y, ion, comp), Ind_1(x, y, ion, comp), -Fc0x, iter);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,Ni,comp),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x + 1, y, ion, comp), Ind_1(x, y, Ni, comp), -Fph0x, iter);
                        ind++;
                    }
                    if(x>0)
                    {
                        Ftmpx = Dcs[c_index(x-1,y,comp,ion)*2]*(cp[c_index(x-1,y,comp,ion)]+cp[c_index(x,y,comp,ion)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/c[c_index(x,y,comp,ion)];
                        Fph1x = z[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,ion,comp),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x - 1, y, ion, comp), Ind_1(x, y, ion, comp), -Fc1x, iter);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,Ni,comp),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x - 1, y, ion, comp), Ind_1(x, y, Ni, comp), -Fph1x, iter);
                        ind++;
                    }
                    if(y<Ny-1)
                    {
                        Ftmpy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2/dy*dt/dy;
                        Fc0y = Ftmpy/c[c_index(x,y,comp,ion)];
                        Fph0y = z[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,ion,comp),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x, y + 1, ion, comp), Ind_1(x, y, ion, comp), -Fc0y, iter);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,Ni,comp),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x, y + 1, ion, comp), Ind_1(x, y, Ni, comp), -Fph0y, iter);
                        ind++;
                    }
                    if(y>0)
                    {
                        Ftmpy = Dcs[c_index(x,y-1,comp,ion)*2+1]*(cp[c_index(x,y-1,comp,ion)]+cp[c_index(x,y,comp,ion)])/2/dy*dt/dy;
                        Fc1y = Ftmpy/c[c_index(x,y,comp,ion)];
                        Fph1y = z[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,ion,comp),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x, y - 1, ion, comp), Ind_1(x, y, ion, comp), -Fc1y, iter);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,Ni,comp),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                        find_print(Ind_1(x, y - 1, ion, comp), Ind_1(x, y, Ni, comp), -Fph1y, iter);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x,y,comp)]+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;
      
                    //membrane current contributions
                    Ac+=flux->dfdci[c_index(x,y,comp,ion)]*dt;
                    Aphi+=flux->dfdphim[c_index(x,y,comp,ion)]*dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,ion,comp),-flux->dfdci[c_index(x,y,comp,ion)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, Nc - 1), Ind_1(x, y, ion, comp),
                               -flux->dfdci[c_index(x, y, comp, ion)] * dt, iter);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,ion,Nc-1),flux->dfdce[c_index(x,y,comp,ion)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, comp), Ind_1(x, y, ion, Nc - 1),
                               flux->dfdce[c_index(x, y, comp, ion)] * dt, iter);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni,comp),-flux->dfdphim[c_index(x,y,comp,ion)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, Nc - 1), Ind_1(x, y, Ni, comp),
                               -flux->dfdphim[c_index(x, y, comp, ion)] * dt, iter);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni,Nc-1),-flux->dfdphim[c_index(x,y,comp,ion)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, comp), Ind_1(x, y, Ni, Nc - 1),
                               -flux->dfdphim[c_index(x, y, comp, ion)] * dt, iter);
                    ind++;
                    //Volume terms
                    //C extra with intra alpha
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni+1,comp),-c[c_index(x,y,Nc-1,ion)],INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, Nc - 1), Ind_1(x, y, Ni + 1, comp), -c[c_index(x, y, Nc - 1, ion)], iter);
                    ind++;
                    //C intra with intra alpha
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni+1,comp),c[c_index(x,y,comp,ion)],INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, comp), Ind_1(x, y, Ni + 1, comp), c[c_index(x, y, comp, ion)], iter);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,ion,comp),Ac,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, comp), Ind_1(x, y, ion, comp), Ac, iter);
                    ind++;
                     // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,ion,comp),Ind_1(x,y,Ni,comp),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, ion, comp), Ind_1(x, y, Ni, comp), Aphi, iter);
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
                if(x<Nx-1)
                {
                    Ftmpx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2/dx*dt/dx;
                    Fc0x = Ftmpx/c[c_index(x,y,comp,ion)];
                    Fph0x = z[ion]*Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,ion,comp),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x + 1, y, ion, comp), Ind_1(x, y, ion, comp), -Fc0x, iter);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,ion,comp),Ind_1(x,y,Ni,comp),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x + 1, y, ion, comp), Ind_1(x, y, Ni, comp), -Fph0x, iter);
                    ind++;
                }
                if(x>0)
                {
                    Ftmpx = Dcs[c_index(x-1,y,comp,ion)*2]*(cp[c_index(x-1,y,comp,ion)]+cp[c_index(x,y,comp,ion)])/2/dx*dt/dx;
                    Fc1x = Ftmpx/c[c_index(x,y,comp,ion)];
                    Fph1x = z[ion]*Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,ion,comp),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x - 1, y, ion, comp), Ind_1(x, y, ion, comp), -Fc1x, iter);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,ion,comp),Ind_1(x,y,Ni,comp),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x - 1, y, ion, comp), Ind_1(x, y, Ni, comp), -Fph1x, iter);
                    ind++;
                }
                if(y<Ny-1)
                {
                    Ftmpy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2/dy*dt/dy;
                    Fc0y = Ftmpy/c[c_index(x,y,comp,ion)];
                    Fph0y = z[ion]*Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,ion,comp),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y + 1, ion, comp), Ind_1(x, y, ion, comp), -Fc0y, iter);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,ion,comp),Ind_1(x,y,Ni,comp),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y + 1, ion, comp), Ind_1(x, y, Ni, comp), -Fph0y, iter);
                    ind++;
                }
                if(y>0)
                {
                    Ftmpy = Dcs[c_index(x,y-1,comp,ion)*2+1]*(cp[c_index(x,y-1,comp,ion)]+cp[c_index(x,y,comp,ion)])/2/dy*dt/dy;
                    Fc1y = Ftmpy/c[c_index(x,y,comp,ion)];
                    Fph1y = z[ion]*Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,ion,comp),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y - 1, ion, comp), Ind_1(x, y, ion, comp), -Fc1y, iter);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,ion,comp),Ind_1(x,y,Ni,comp),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y - 1, ion, comp), Ind_1(x, y, Ni, comp), -Fph1y, iter);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1-al[al_index(x,y,0)]-al[al_index(x,y,1)])+Fc0x+Fc1x+Fc0y+Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;
            
                //Membrane current contribution
                for(comp=0;comp<Nc-1;comp++)
                {
                    Ac -= flux->dfdce[c_index(x,y,comp,ion)]*dt;
                    Aphi += flux->dfdphim[c_index(x,y,comp,ion)]*dt; 
                }
                //Add bath contributions
                Ftmpx=sqrt(pow(Dcb[c_index(x,y,Nc-1,ion)*2],2)+pow(Dcb[c_index(x,y,Nc-1,ion)*2+1],2));
                Ac -= Ftmpx*(cp[c_index(x,y,Nc-1,ion)]+cbath[ion])/(2*c[c_index(x,y,Nc-1,ion)])*dt;
                Aphi -= Ftmpx*(cp[c_index(x,y,Nc-1,ion)]+cbath[ion])*z[ion]/2*dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,ion,Nc-1),Ac,INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, ion, Nc - 1), Ind_1(x, y, ion, Nc - 1), Ac, iter);
                ind++;
                 // c with phi
                ierr = MatSetValue(Jac,Ind_1(x,y,ion,Nc-1),Ind_1(x,y,Ni,Nc-1),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, ion, Nc - 1), Ind_1(x, y, Ni, Nc - 1), Aphi, iter);
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
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,ion,comp),z[ion]*al[al_index(x,y,comp)],INSERT_VALUES); CHKERRQ(ierr);
                    find_print(Ind_1(x, y, Ni, comp), Ind_1(x, y, ion, comp), z[ion] * al[al_index(x, y, comp)], iter);
                    ind++;
                }
                //Phi with C extracellular one
                comp = Nc-1;
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,ion,comp),z[ion]*(1-al[al_index(x,y,0)]-al[al_index(x,y,1)]),INSERT_VALUES); CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, comp), Ind_1(x, y, ion, comp),
                           z[ion] * (1 - al[al_index(x, y, 0)] - al[al_index(x, y, 1)]), iter);
                ind++;

            }
            //electroneutrality-voltage entries
            Aphi = 0;
            for(comp=0;comp<Nc-1;comp++)
            {
                Aphi -= cm[comp];
            }
            //extraphi with extra phi
            ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni,Nc-1),Aphi,INSERT_VALUES);CHKERRQ(ierr);
            find_print(Ind_1(x, y, Ni, Nc - 1), Ind_1(x, y, Ni, Nc - 1), Aphi, iter);
            ind++;
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extra phi with intra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni,comp),cm[comp],INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, Nc - 1), Ind_1(x, y, Ni, comp), cm[comp], iter);
                ind++;
                // Intra phi with Extraphi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni,Nc-1),cm[comp],INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, comp), Ind_1(x, y, Ni, Nc - 1), cm[comp], iter);
                ind++;
                //Intra phi with Intra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni,comp),-cm[comp],INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, comp), Ind_1(x, y, Ni, comp), -cm[comp], iter);
                ind++;
                //Extra phi with intra-Volume
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,Nc-1),Ind_1(x,y,Ni+1,comp),-cz(c,z,x,y,Nc-1),INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, Nc - 1), Ind_1(x, y, Ni + 1, comp), -cz(c, z, x, y, Nc - 1), iter);
                ind++;
                //Intra phi with Intra Vol
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni,comp),Ind_1(x,y,Ni+1,comp),cz(c,z,x,y,comp),INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni, comp), Ind_1(x, y, Ni + 1, comp), cz(c, z, x, y, comp), iter);
                ind++;
            }
        }
    }
    //water flow
    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            for(comp=0;comp<Nc-1;comp++)
            {
                //Water flow volume fraction entries
                //Volume to Volume
                Ac=1+(flux->dwdpi[al_index(x,y,comp)]*(con_vars->ao[phi_index(0,0,Nc-1)]/(pow(1-al[al_index(x,y,0)]-al[al_index(x,y,1)],2))+con_vars->ao[phi_index(0,0,comp)]/pow(al[al_index(x,y,comp)],2))+flux->dwdal[al_index(x,y,comp)])*dt;
                ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,comp),Ac,INSERT_VALUES);CHKERRQ(ierr);
                find_print(Ind_1(x, y, Ni + 1, comp), Ind_1(x, y, Ni + 1, comp), Ac, iter);
                ind++;
                //Off diagonal (from aNc=1-sum(ak))
                for (PetscInt l=0; l<comp; l++)
                {
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,l),flux->dwdpi[al_index(x,y,comp)]*con_vars->ao[phi_index(0,0,Nc-1)]/pow(1-al[al_index(x,y,0)]-al[al_index(x,y,1)],2)*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, Ni + 1, comp), Ind_1(x, y, Ni + 1, l),
                               flux->dwdpi[al_index(x, y, comp)] * con_vars->ao[phi_index(0, 0, Nc - 1)] /
                               pow(1 - al[al_index(x, y, 0)] - al[al_index(x, y, 1)], 2) * dt, iter);
                    ind++;
                }
                for (PetscInt l=comp+1; l<Nc-1; l++)
                {
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,Ni+1,l),flux->dwdpi[al_index(x,y,comp)]*con_vars->ao[phi_index(0,0,Nc-1)]/((1-al[al_index(x,y,0)]-al[al_index(x,y,1)])*(1-al[al_index(x,y,0)]-al[al_index(x,y,1)]))*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, Ni + 1, comp), Ind_1(x, y, Ni + 1, l),
                               flux->dwdpi[al_index(x, y, comp)] * con_vars->ao[phi_index(0, 0, Nc - 1)] /
                               ((1 - al[al_index(x, y, 0)] - al[al_index(x, y, 1)]) *
                                (1 - al[al_index(x, y, 0)] - al[al_index(x, y, 1)])) * dt, iter);
                    ind++;
                }
                for (ion=0; ion<Ni; ion++)
                {
                  //Volume to extra c
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,ion,Nc-1),flux->dwdpi[al_index(x,y,comp)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, Ni + 1, comp), Ind_1(x, y, ion, Nc - 1),
                               flux->dwdpi[al_index(x, y, comp)] * dt, iter);
                    ind++;
                  //Volume to intra c
                    ierr = MatSetValue(Jac,Ind_1(x,y,Ni+1,comp),Ind_1(x,y,ion,comp),-flux->dwdpi[al_index(x,y,comp)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                    find_print(Ind_1(x, y, Ni + 1, comp), Ind_1(x, y, ion, comp),
                               -flux->dwdpi[al_index(x, y, comp)] * dt, iter);
                    ind++;
                }
            }
        }
    }

    
//    for(int i=0;i<NA;i++)
//    {
//        ierr = MatSetValue(Jac,i,i,2,INSERT_VALUES);CHKERRQ(ierr);
//    }
    

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//    printf("Nz: %d, Ind: %d\n",Nz,ind);
    // MatView(Jac,PETSC_VIEWER_STDOUT_SELF);
//    PetscViewer viewer;
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",&viewer);
//    MatView(Jac,viewer);
    // MatView(Jac,PETSC_VIEWER_DRAW_WORLD);
//    PetscViewerDestroy(&viewer); /* indicate failure.*/
//     while(1){};
    return ierr;
}

