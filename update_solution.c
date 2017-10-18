#include "constants.h"
#include "functions.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void newton_solve(struct SimState *state_vars, double dt, struct GateType *gvars, struct ExctType *gexct, struct ConstVars *con_vars,struct Solver *slvr,struct FluxData *flux) 
{
    //Save the "current" aka past state
	struct SimState *state_vars_past;
    state_vars_past = (struct SimState*)malloc(sizeof(struct SimState));
    memcpy(state_vars_past,state_vars,sizeof(struct SimState)); 

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at even positions (0,2,4,...)
    //y at odd (1,3,5,...)
    //still use c_index(x,y,comp,ion), but with ind*2 or ind*2+1
   	double Dcs[Nx*Ny*Ni*Nc*2];
   	double Dcb[Nx*Ny*Ni*Nc*2];
   	//compute diffusion coefficients
   	diff_coef(Dcs,state_vars->alpha,1);
    for(int ion=0;ion<Ni;ion++)
    {
        for(int comp=0;comp<Nc;comp++)
        {
            printf("Dcs: Ion %d, Comp %d ",ion,comp);
            printf("Dcs x: %f, Dcs y: %f\n",1e6*Dcs[c_index(0,0,comp,ion)*2],1e6*Dcs[c_index(0,0,comp,ion)*2+1]);
        }
    }
   	//Bath diffusion
  	diff_coef(Dcb,state_vars->alpha,Batheps);
    for(int ion=0;ion<Ni;ion++)
    {
        for(int comp=0;comp<Nc;comp++)
        {
            printf("Dcb: Ion %d, Comp %d ",ion,comp);
            printf("Dcb x: %f, Dcb y: %f\n",1e6*Dcb[c_index(0,0,comp,ion)*2],1e6*Dcb[c_index(0,0,comp,ion)*2+1]);
        }
    }

    double tol = reltol*array_max(state_vars->c,(size_t)Nx*Ny*Ni*Nc);
    double rsd = tol+1;

    int x=0;int y=0;
    for(int iter=0;iter<itermax;iter++)
    {
        //Compute membrane ionic flux relation quantitites
        ionmflux(flux,state_vars,state_vars_past,gvars,gexct,con_vars);
        for(int ion=0;ion<Ni;ion++)
        {
            for(int comp=0;comp<Nc;comp++)
            {
                printf("Ion: %d, Comp %d\n",ion,comp);
                printf("Flux*1e6: %f, dfdci: %f, dfdce: %f, dfdphim: %f\n",1e6*flux->mflux[c_index(0,0,comp,ion)],flux->dfdci[c_index(0,0,comp,ion)],flux->dfdce[c_index(0,0,comp,ion)],flux->dfdphim[c_index(0,0,comp,ion)]);
            }
        }
        //Compute membrane water flow related quantities
        wflowm(flux,state_vars,con_vars);
        for(int comp=0;comp<Nc;comp++)
        {
            printf("Comp: %d\n",comp);
            printf("wFlux: %f,%f,%f\n",flux->wflow[al_index(x,y,comp)],flux->dwdpi[al_index(x,y,comp)],flux->dwdal[al_index(x,y,comp)]);
        }
        // VecView(slvr->Res,PETSC_VIEWER_STDOUT_SELF);
        calc_residual(slvr->Res,state_vars_past,state_vars,dt,Dcs,Dcb,flux,con_vars);
        // VecView(slvr->Res,PETSC_VIEWER_STDOUT_SELF);
        break;


    }

    if(rsd>tol)
    {
        fprintf(stderr, "Netwon Iteration did not converge! Stopping...\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

  	return;
}

/*
  al = fullalpha(state_vars.alpha)
  Dcs = diff_coeff(al) #compute diffusion coefficients
  Dcb=Batheps*diff_coeff(al)
  tol = reltol*maximum(abs,state_vars.c)
  rsd = tol+1
  iter = 0
  mflux = zeros(size(state_vars.c)) #initialize ion flux record
  for iter = 0:itermax
    Res = F(state_vars_past, state_vars, dt, Dcs, Dcb, fluxdata,con_vars)
    #if residual is small or iteration count is large, exit loop
    rsd = norm(Res,Inf)
#    write(diagnostics, "$(rsd)\n")
    if (rsd < tol)
        break
    end
    #println(iter)
    Jac = J(state_vars_past, state_vars, dt, Dcs, Dcb, fluxdata,con_vars)
    # println(full(Jac[1:2*Nv,1:2*Nv]))
    # Ka=0
    # K=Nx*Ny
    # println("Full Jacobian det: $(det(full(Jac[(1+Ka*Nv):K*Nv,(1+Ka*Nv):K*Nv])))")
    # println("Full cond:$(cond(full(Jac[(1+Ka*Nv):K*Nv,(1+Ka*Nv):K*Nv])))")
    # println("Full rank: $(rank(full(Jac[(1+Ka*Nv):K*Nv,(1+Ka*Nv):K*Nv])))")
    # println("Size: $(size(Jac[(1+Ka*Nv):K*Nv,(1+Ka*Nv):K*Nv]))")
#     if iter==1
#       global Jac1 = full(Jac)
#     end
    if !(use_direct_solve)
      #Q = -call_gmres(Jac, Res)
     #   @time Q = -iterative_solve(Jac,Res)
        # println("Preconditioner")
        P=Preconditioner_ILU(Jac)
        # b=rand(NA)
        # println("Precond norm: $(norm(P\b)/norm(b))")
        # @time P=Preconditioner_ADI(state_vars_past, state_vars, dt, Dcs,Dcb, fluxdata,con_vars)
        # println("Solve")
        # Q = -iterative_solve(Jac,Res,P,NA)
        println("Iterativesolve")
        @time Q = -iterative_solve(Jac,Res,P,NA)
        println("Iterative solve done")
    else
        @time Q = -Jac\Res#solve for correction
    end
    #   println("iteration: $(iter), residual: $(rsd)")
    #   state_vars.p=1
    #update c, phi, and alpha
    for k = 1:Nc
        for j=1:Ni
            state_vars.c[:,:,j,k] += reshape(Q[Ind_1(1,1,j,k):Nv:NA],Nx,Ny)
        end
        state_vars.phi[:,:,k] += reshape(Q[Ind_1(1,1,Ni+1,k):Nv:NA],Nx,Ny)
    end
    for k = 1:Nc-1
        state_vars.alpha[:,:,k] += reshape(Q[Ind_1(1,1,Ni+2,k):Nv:NA],Nx,Ny)
    end
    al = fullalpha(state_vars.alpha)
  end
  if rsd > tol
    println("REACHED MAXIMUM NEWTON SOLVE STEPS")
  end
  if details
    println("iteration: $(iter), residual: $(rsd)")
  end
  state_vars
end
*/
PetscErrorCode calc_residual(Vec Res,struct SimState *state_vars_past,struct SimState *state_vars,double dt,double *Dcs,double *Dcb,struct FluxData *flux,struct ConstVars *con_vars)
{
    double *c = state_vars->c;
    double *phi = state_vars->phi;
    double *al = state_vars->alpha;
    double *cp = state_vars_past->c;
    double *alp = state_vars_past->alpha; 
    //Residual for concentration equations
    double Rcvx,Rcvy,Resc;
    double RcvxRight,RcvyUp;

    double alNc,alpNc;
    int ion,comp,x,y;

    PetscErrorCode ierr;

    for(x=0;x<Nx-1;x++)
    {
        for(y=0;y<Ny-1;y++)
        {
            for(ion=0;ion<Ni;ion++)
            {
                for(comp=0;comp<Nc-1;comp++)
                {
                    //First difference term
                    Rcvx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
                    Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x+1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x+1,y,comp)]))/dx*dt/dx;
                    //Add Second right moving difference
                    RcvxRight = 0;
                    if(x<Nx-2)
                    {
                        RcvxRight = Dcs[c_index(x+1,y,comp,ion)*2]*(cp[c_index(x+1,y,comp,ion)]+cp[c_index(x+2,y,comp,ion)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x+2,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x+2,y,comp)]))/dx*dt/dx;
                    }  
                    //Up down difference 
                    Rcvy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
                    Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y+1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y+1,comp)]))/dy*dt/dy;
                    //Next upward difference
                    RcvyUp = 0;
                    if(y<Ny-2)
                    {
                        RcvyUp = Dcs[c_index(x,y+1,comp,ion)*2+1]*(cp[c_index(x,y+1,comp,ion)]+cp[c_index(x,y+2,comp,ion)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y+2,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y+2,comp)]))/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)]-alp[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;

                    ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);

                }
                //Set Extracellular values
                alNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
                alpNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
                comp = Nc-1;
                //First difference term
                Rcvx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
                Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x+1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x+1,y,comp)]))/dx*dt/dx;
                //Add Second right moving difference
                RcvxRight = 0;
                if(x<Nx-2)
                {
                    RcvxRight = Dcs[c_index(x+1,y,comp,ion)*2]*(cp[c_index(x+1,y,comp,ion)]+cp[c_index(x+2,y,comp,ion)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x+2,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x+2,y,comp)]))/dx*dt/dx;
                }  
                //Up down difference 
                Rcvy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
                Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y+1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y+1,comp)]))/dy*dt/dy;
                //Next upward difference
                RcvyUp = 0;
                if(y<Ny-2)
                {
                    RcvyUp = Dcs[c_index(x,y+1,comp,ion)*2+1]*(cp[c_index(x,y+1,comp,ion)]+cp[c_index(x,y+2,comp,ion)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y+2,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y+2,comp)]))/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,comp,ion)]-alpNc*c[c_index(x,y,comp,ion)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;
                //Add bath variables
                Resc += sqrt(pow(Dcb[c_index(x,y,comp,ion)*2],2)+pow(Dcb[c_index(x,y,comp,ion)*2+1],2))*(cp[c_index(x,y,comp,ion)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp)]-z[ion]*phibath)*dt;

                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
    }
    // Add the endpoints (x,Ny-1)
    y=Ny-1;
    for(x=0;x<Nx-1;x++)
    { 
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc-1;comp++)
            {
                //First difference term
                Rcvx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
                Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x+1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x+1,y,comp)]))/dx*dt/dx;
                //Add Second right moving difference
                RcvxRight = 0;
                if(x<Nx-2)
                {
                    RcvxRight = Dcs[c_index(x+1,y,comp,ion)*2]*(cp[c_index(x+1,y,comp,ion)]+cp[c_index(x+2,y,comp,ion)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x+2,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x+2,y,comp)]))/dx*dt/dx;
                }  
                //Up down difference 
                Rcvy = 0;
                //Next upward difference
                RcvyUp = 0;
                Resc = al[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)]-alp[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;

                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);

            }
            //Set Extracellular values
            alNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
            alpNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
            comp = Nc-1;
            //First difference term
            Rcvx = Dcs[c_index(x,y,comp,ion)*2]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x+1,y,comp,ion)])/2;
            Rcvx = Rcvx*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x+1,y,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x+1,y,comp)]))/dx*dt/dx;
            //Add Second right moving difference
            RcvxRight = 0;
            if(x<Nx-2)
            {
                RcvxRight = Dcs[c_index(x+1,y,comp,ion)*2]*(cp[c_index(x+1,y,comp,ion)]+cp[c_index(x+2,y,comp,ion)])/2;
                RcvxRight = RcvxRight*(log(c[c_index(x+1,y,comp,ion)])-log(c[c_index(x+2,y,comp,ion)])+z[ion]*(phi[phi_index(x+1,y,comp)]-phi[phi_index(x+2,y,comp)]))/dx*dt/dx;
            }  
            //Up down difference 
            Rcvy = 0;
            //Next upward difference
            RcvyUp = 0;

            Resc = alNc*c[c_index(x,y,comp,ion)]-alpNc*c[c_index(x,y,comp,ion)];
            Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;
            //Add bath variables
            Resc += sqrt(pow(Dcb[c_index(x,y,comp,ion)*2],2)+pow(Dcb[c_index(x,y,comp,ion)*2+1],2))*(cp[c_index(x,y,comp,ion)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp)]-z[ion]*phibath)*dt;

            ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);
        }

    }
    //Add the enpoints (Nx-1,y)
    x=Nx-1;
    for(y=0;y<Ny-1;y++)
    {
        for(ion=0;ion<Ni;ion++)
        {
            for(comp=0;comp<Nc-1;comp++)
            {
                //First difference term
                Rcvx = 0;
                //Add Second right moving difference
                RcvxRight = 0;
                //Up down difference 
                Rcvy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
                Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y+1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y+1,comp)]))/dy*dt/dy;
                //Next upward difference
                RcvyUp = 0;
                if(y<Ny-2)
                {
                    RcvyUp = Dcs[c_index(x,y+1,comp,ion)*2+1]*(cp[c_index(x,y+1,comp,ion)]+cp[c_index(x,y+2,comp,ion)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y+2,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y+2,comp)]))/dy*dt/dy;
                }
                Resc = al[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)]-alp[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;

                ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);

            }
            //Set Extracellular values
            alNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
            alpNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
            comp = Nc-1;
            //First difference term
            Rcvx = 0;
            //Add Second right moving difference
            RcvxRight = 0; 
            //Up down difference 
            Rcvy = Dcs[c_index(x,y,comp,ion)*2+1]*(cp[c_index(x,y,comp,ion)]+cp[c_index(x,y+1,comp,ion)])/2;
            Rcvy = Rcvy*(log(c[c_index(x,y,comp,ion)])-log(c[c_index(x,y+1,comp,ion)])+z[ion]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y+1,comp)]))/dy*dt/dy;
            //Next upward difference
            RcvyUp = 0;
            if(y<Ny-2)
            {
                RcvyUp = Dcs[c_index(x,y+1,comp,ion)*2+1]*(cp[c_index(x,y+1,comp,ion)]+cp[c_index(x,y+2,comp,ion)])/2;
                RcvyUp = RcvyUp*(log(c[c_index(x,y+1,comp,ion)])-log(c[c_index(x,y+2,comp,ion)])+z[ion]*(phi[phi_index(x,y+1,comp)]-phi[phi_index(x,y+2,comp)]))/dy*dt/dy;
            }
            Resc = alNc*c[c_index(x,y,comp,ion)]-alpNc*c[c_index(x,y,comp,ion)];
            Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;
            //Add bath variables
            Resc += sqrt(pow(Dcb[c_index(x,y,comp,ion)*2],2)+pow(Dcb[c_index(x,y,comp,ion)*2+1],2))*(cp[c_index(x,y,comp,ion)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp)]-z[ion]*phibath)*dt;

            ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);
        }
    }
    //Add the top right point: (Nx-1,Ny-1)
    x=Nx-1;y=Ny-1;
    for(ion=0;ion<Ni;ion++)
    {
        for(comp=0;comp<Nc-1;comp++)
        {
            //First difference term
            Rcvx = 0;
            //Add Second right moving difference
            RcvxRight = 0;  
            //Up down difference 
            Rcvy = 0;
            //Next upward difference
            RcvyUp = 0;

            Resc = al[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)]-alp[al_index(x,y,comp)]*c[c_index(x,y,comp,ion)];
            Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;

            ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);

        }
        //Set Extracellular values
        alNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
        alpNc = 1 - al[al_index(x,y,0)] - al[al_index(x,y,1)];
        comp = Nc-1;
        //First difference term
        Rcvx = 0;
        //Add Second right moving difference
        RcvxRight = 0; 
        //Up down difference 
        Rcvy = 0;
        //Next upward difference
        RcvyUp = 0;

        Resc = alNc*c[c_index(x,y,comp,ion)]-alpNc*c[c_index(x,y,comp,ion)];
        Resc += Rcvx - RcvxRight + Rcvy - RcvyUp + flux->mflux[c_index(x,y,comp,ion)]*dt;
        //Add bath variables
        Resc += sqrt(pow(Dcb[c_index(x,y,comp,ion)*2],2)+pow(Dcb[c_index(x,y,comp,ion)*2+1],2))*(cp[c_index(x,y,comp,ion)]+cbath[ion])/2.0*(log(c[c_index(x,y,comp,ion)])-log(cbath[ion])+z[ion]*phi[phi_index(x,y,comp)]-z[ion]*phibath)*dt;

        ierr = VecSetValue(Res,Ind_1(x,y,ion,comp),Resc,INSERT_VALUES);CHKERRQ(ierr);
    }

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {    
            //Residual for electroneutrality condition
            for(comp=0;comp<Nc-1;comp++)
            {

                Resc = al[al_index(x,y,comp)]*cz(c,z,x,y,comp)+con_vars->zo[phi_index(0,0,comp)]+con_vars->ao[phi_index(0,0,comp)];
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),Resc,INSERT_VALUES); CHKERRQ(ierr);
            }
            //Extracellular term
            comp=Nc-1;
            Resc = (1-al[al_index(x,y,0)]-al[al_index(x,y,1)])*cz(c,z,x,y,comp)+con_vars->zo[phi_index(0,0,comp)]+con_vars->ao[phi_index(0,0,comp)];
            ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),Resc,INSERT_VALUES); CHKERRQ(ierr);
            
            //Residual for water flow
            //Plus modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++)
            {
                //Water flow
                ierr = VecSetValue(Res,Ind_1(x,y,Ni+1,comp),al[al_index(x,y,comp)]-alp[al_index(x,y,comp)]+flux->wflow[al_index(x,y,comp)]*dt,INSERT_VALUES);

                //Extracell voltage
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,Nc-1),-cm[comp]*(phi[phi_index(x,y,Nc-1)]-phi[phi_index(x,y,comp)]),ADD_VALUES);
                //Intracell voltage mod
                ierr = VecSetValue(Res,Ind_1(x,y,Ni,comp),-cm[comp]*(phi[phi_index(x,y,comp)]-phi[phi_index(x,y,Nc-1)]),ADD_VALUES);
            }
        }
    }
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    // VecView(Res,PETSC_VIEWER_STDOUT_SELF);
    // printf("c: %f,%f,%f\n",al[al_index(0,0,0)],al[al_index(0,0,1)],al[al_index(3,0,0)]);



    return ierr;
}

