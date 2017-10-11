/* Hello World program */

#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

int main(int argc, char **argv)
{
    printf("\n\n\nGrid size: %dx%d, with %d ions, and %d compartments.\n",Nx,Ny,Ni,Nc);
    
    //Create state_variables struct
    struct SimState *state_vars;
    state_vars=(struct SimState*)malloc(sizeof(struct SimState));

    //Initialize
    init(state_vars);

    printf("%f,%f,%f\n",state_vars->c[0],state_vars->phi[10],state_vars->alpha[25]);

    //Create the constant ion channel vars
    struct ConstVars *con_vars;
    con_vars=(struct ConstVars*)malloc(sizeof(struct ConstVars));
    
    //Set the constant variables
    set_params(state_vars,con_vars);





    return 0;
}

/*
state_vars = SimState(c0, phi0, alpha0)
(state_vars,con_vars) = set_params(state_vars)
krecord = false #true/false which determines whether or not to record values at any given time step
savefreq = 500
println("Initialization Routine")
global sIndex=sparseindex_1()
bIndex_x=bandindex_x()
bIndex_y=bandindex_y()
permute_x_to_y(Nx,Ny,Nv) #generate permutation matrix in Iterative_solver
F=numerical.F_1
J=numerical.J_1
initialdata = initdata
state_vars, gvars, stoptest = initialdata(state_vars, F, J)
gexct = exct(0)
numrecords = convert(Int64,floor(Nt/krecordfreq) + 1) #number of times to record values
#initialize arrays to save data
carray = zeros(Nx,Ny,Ni,Nc,numrecords)
phiarray = zeros(Nx,Ny,Nc,numrecords)
alphaarray = zeros(Nx,Ny,Nc-1,numrecords)
cpoints = zeros(2,Ni,Nc,Nt+1)
phipoints = zeros(2,Nc,Nt+1)
alphapoints = zeros(2,Nc-1,Nt+1)
if !fluid_flow
  state_vars_all = SimState(carray, phiarray, alphaarray)
  state_vars_point = SimState(cpoints,phipoints,alphapoints)
else
  parray = zeros(Nx,Ny,Nc,numrecords)
  uarray = zeros(Nx,Ny,Nc,numrecors)
  ppoints = zeros(2,Nc,Nt+1)
  upoints = zeros(2,Nc,Nt+1)
  state_vars_all = SimStateF(carray, phiarray, alphaarray, parray, uarray)
  state_vars_point = SimStateF(cpoints,phipoints,alphapoints,ppoints,upoints)
end
#save initial data
thirdsint=convert(Int64,floor(Nx*.3))
state_vars_all.c[:,:,:,:,1] = state_vars.c
state_vars_all.phi[:,:,:,1] = state_vars.phi
state_vars_all.alpha[:,:,:,1] = state_vars.alpha
state_vars_point.c[1,:,:,1] = state_vars.c[thirdsint,1,:,:]
state_vars_point.phi[1,:,1] = state_vars.phi[thirdsint,1,:]
state_vars_point.alpha[1,:,1] = state_vars.alpha[thirdsint,1,:]
state_vars_point.c[2,:,:,1] = state_vars.c[1,thirdsint,:,:]
state_vars_point.phi[2,:,1] = state_vars.phi[1,thirdsint,:]
state_vars_point.alpha[2,:,1] = state_vars.alpha[1,thirdsint,:]
if stoptest
  println("Initialization Routine Failed to Converge")
else

  println("Main Computation Routine")
  k = 0    #time step counter
  j = 1    #recording counter
  t = dt
  @time while t <= (Time + dt/2)
    println("Newton_solve:")
   @time state_vars = numerical.newton_solve(F, J, state_vars, dt, gvars, gexct,con_vars,sIndex)
   println("Newton_solve Done")
    gvars = gatevars(gvars,state_vars,false)
    gexct = exct(t)
    k += 1
    state_vars_point.c[1,:,:,k] = state_vars.c[thirdsint,1,:,:]
    state_vars_point.phi[1,:,k] = state_vars.phi[thirdsint,1,:]
    state_vars_point.alpha[1,:,k] = state_vars.alpha[thirdsint,1,:]
    state_vars_point.c[2,:,:,k] = state_vars.c[1,thirdsint,:,:]
    state_vars_point.phi[2,:,k] = state_vars.phi[1,thirdsint,:]
    state_vars_point.alpha[2,:,k] = state_vars.alpha[1,thirdsint,:]
    if mod(k,krecordfreq)==0
      krecord = true
    end
    if krecord
      j += 1
      state_vars_all.c[:,:,:,:,j] = state_vars.c
      state_vars_all.phi[:,:,:,j] = state_vars.phi
      state_vars_all.alpha[:,:,:,j] = state_vars.alpha
      krecord = false
      println("Time Step: $j\n")
    end
    # if mod(k,savefreq)==0
    #   if !two_points_exct
    #     matwrite("csd2dres.mat", ["state_vars" => state_vars_all, "gvars" => gvars, "gexct" => gexct, "t" => t, "points" => state_vars_point])
    #   else
    #     matwrite("csd2dres2points.mat", []"state_vars" => state_vars_all, "gvars" => gvars, "gexct"=> gexct, "t" => t, "points" => state_vars_point])
    #   end
    # end
    t += dt
  end

end
*/