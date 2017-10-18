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
   	//Bath diffusion
  	diff_coef(Dcb,state_vars->alpha,Batheps);

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
        calc_residual(slvr.Res,state_vars_past,state_vars,dt,Dcs,Dcb,flux,con_vars);
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
    #compute membrane water flow related quantities
    piw = squeeze(sum(state_vars.c,3),3)+con_vars.ao./al
    wflow, dwdpi, dwdal = wflowm(piw, al)
    fluxdata = FluxData(mflux,dfdci,dfdce,dfdphim,wflow,dwdpi,dwdal)
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