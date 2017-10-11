#include "constants.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int c_index(x,y,comp,ion)
{
	return Nc*Ni* (Nx * y + x) + comp*Ni+ion; 
}
int phi_index(x,y,comp)
{
	return Nc* (Nx * y + x) + comp; 
}
int al_index(x,y,comp)
{
	return (Nc-1)* (Nx * y + x) + comp; 
}


void init(struct SimState *state_vars)
{
	for(int x=0;x<Nx;x++)
	{
		for(int y=0;y<Ny;y++)
		{
			//initial volume fractions
			state_vars->alpha[al_index(x,y,0)]=alphao[0];
			state_vars->alpha[al_index(x,y,1)]=alphao[1]; 
			//initial voltages (dimensionless)
			state_vars->phi[phi_index(x,y,0)] = -70/RTFC; //neuronal voltage
			state_vars->phi[phi_index(x,y,1)] = -85/RTFC; //glial voltage
			state_vars->phi[phi_index(x,y,2)] = -0/RTFC; //extracell voltage
			//initial concentrations in mmol/cm^3=1e-3 mmol/l
			state_vars->c[c_index(x,y,0,0)] = 10e-3;     //neuronal Na concentration
			state_vars->c[c_index(x,y,1,0)] = 10e-3;      //glial Na concentration
			state_vars->c[c_index(x,y,2,0)] = 140e-3;     //extracellular Na concentration
			state_vars->c[c_index(x,y,0,1)] = 130e-3;     //neuronal K concentration
			state_vars->c[c_index(x,y,1,1)] = 130e-3;     //glial K concentration
			state_vars->c[c_index(x,y,2,1)] = 3.4e-3;     //extracellular K concentration
			state_vars->c[c_index(x,y,0,2)] = 10e-3;       //neuronal Cl concentration
			state_vars->c[c_index(x,y,1,2)] = 10e-3; 		//glial Cl concentraion
			state_vars->c[c_index(x,y,2,2)] = 120e-3;       //143.5e-3%extracellular Cl

		}
	}
}

void set_params(struct SimState* state_vars,struct ConstVars* con_vars)
{
	//Everything that follows will asume spatially uniform
	//At rest state
	double c[Ni*Nc];
	double phi[Nc];
	double alpha[Nc];
	for(int comp=0;comp<Nc;comp++)
	{
		for(int ion=0;ion<Ni;ion++)
		{
			c[c_index(0,0,comp,ion)]=state_vars->c[c_index(0,0,comp,ion)];
		}
		phi[phi_index(0,0,comp)]=state_vars->phi[phi_index(0,0,comp)];
		if(comp<Nc-1)
		{
			alpha[al_index(0,0,comp)]=state_vars->alpha[al_index(0,0,comp)];
		}
		else
		{
			alpha[al_index(0,0,comp)]=1; //Adding in the extracellular vol
			for(int i=0;i<Nc-1;i++)
			{
				alpha[al_index(0,0,comp)]-=alpha[al_index(0,0,i)];
			}
		}
	}
	double vm = phi[phi_index(0,0,0)]-phi[phi_index(0,0,2)]; //neuronal membrane potential
    double vmg = phi[phi_index(0,0,1)]-phi[phi_index(0,0,2)]; //glial membrane potential

    //compute neuronal Cl concentration (since neuron has only leak conductance, must be at reversal potential for Cl)

    c[c_index(0,0,0,2)] = c[c_index(0,0,2,2)]*exp(vm);
    //set glial Cl concentration equal to neuronal Cl concentration
    c[c_index(0,0,1,2)] = c[c_index(0,0,0,2)];

    struct FluxPoint *flux;
    flux = (struct FluxPoint*)malloc(sizeof(struct FluxPoint));

    //compute cotransporter permeability so that glial Cl is at rest
    mclin(flux,pClLeakg,-1,c[c_index(0,0,1,2)],c[c_index(0,0,2,2)],vmg,c_index(0,0,1,2));
    con_vars->pNaKCl = -flux->mflux[c_index(0,0,1,2)]/2/log(c[c_index(0,0,1,0)]*c[c_index(0,0,2,2)]*c[c_index(0,0,1,2)]*c[c_index(0,0,1,2)]/(c[c_index(0,0,2,0)]*c[c_index(0,0,2,1)]*c[c_index(0,0,2,2)]*c[c_index(0,0,2,2)]));


    return;
}
/*
function set_params(state_vars)
  #Set remaining membrane parameters either arbitrarily or so that the system is at rest - need to turn this into a
  #function producing global constants and make sure it is called first in main routine
  if rest_state
    #if spatially uniform, reduce to one spatial point for computational simplicity
    if spatially_uniform
      phi = state_vars.phi[1,1,:]
      c = state_vars.c[1,1,:,:]
      alpha = state_vars.alpha[1,1,:]
      fluxdata = zeros(1,1,Ni,Nc-1,4)
    else
      phi = copy(state_vars.phi)
      c = copy(state_vars.c)
      alpha = copy(state_vars.alpha)
      fluxdata = zeros(Nx,Ny,Ni,Nc-1,4)
    end
    vm = phi[:,:,1]-phi[:,:,Nc] #neuronal membrane potential
    vmg = phi[:,:,2]-phi[:,:,Nc] #glial membrane potential

    #compute neuronal Cl concentration (since neuron has only leak conductance, must be at reversal potential for Cl)
    c[:,:,3,1] = c[:,:,3,3].*exp.(vm)
    #set glial Cl concentration equal to neuronal Cl concentration
    c[:,:,3,2] = c[:,:,3,1]

    #compute cotransporter permeability so that glial Cl is at rest
    fluxdata[:,:,3,2,:] = mclin(pClLeakg,-1,c[:,:,3,2],c[:,:,3,3],vmg)
    global const pNaKCl = -fluxdata[:,:,3,2,1]/2./log.(c[:,:,1,2].*c[:,:,2,2].*c[:,:,3,2].^2./(c[:,:,1,3].*c[:,:,2,3].*c[:,:,3,3].^2))
    NaKCl = pNaKCl.*log.(c[:,:,1,2].*c[:,:,2,2].*c[:,:,3,2].^2./(c[:,:,1,3].*c[:,:,2,3].*c[:,:,3,3].^2))

    #compute gating variables
    s = gatevars([],state_vars,true)
    #compute K channel currents (neuron)
    if spatially_uniform
      pKGHK = pKDR*s.gKDR[1,1]+pKA*s.gKA[1,1]
    else
      pKGHK = pKDR*s.gKDR+pKA*s.gKA
    end
    pKLin = pKLeak
    fluxdata[:,:,2,1,:] = mcGoldman(pKGHK,1,c[:,:,2,1],c[:,:,2,Nc],vm)+mclin(pKLin,1,c[:,:,2,1],c[:,:,2,Nc],vm)

    #compute neuronal ATPase value
    global const Imax = fluxdata[:,:,2,1,1].*((1+mK./c[:,:,2,Nc]).^2).*((1+mNa./c[:,:,1,1]).^3)./2

    #compute neuronal sodium currents and leak permeability value
    if spatially_uniform
      pNaGHK = pNaT*s.gNaT[1,1]+pNaP*s.gNaP[1,1]
    else
      pNaGHK = pNaT*s.gNaT+pNaP*s.gNaP
    end
    fluxdata[:,:,1,1,:] = mcGoldman(pNaGHK,1,c[:,:,1,1],c[:,:,1,Nc],vm)
    Ipump = Imax./(((1+mK./c[:,:,2,Nc]).^2).*((1+mNa./c[:,:,1,1]).^3))
    global const pNaLeak = (-fluxdata[:,:,1,1,1]-3*Ipump)./(log.(c[:,:,1,1]./c[:,:,1,3])+vm)

    #compute K channel currents (glial)
    pKLinG = pKIR*inwardrect(c[:,:,2,2],c[:,:,2,Nc],vmg)*pKLeakadjust
    fluxdata[:,:,2,2,:] = mclin(pKLinG,1,c[:,:,2,2],c[:,:,2,Nc],vmg)
    fluxdata[:,:,2,2,1] += NaKCl

    #compute glial ATPase value
    global const Imaxg = fluxdata[:,:,2,2,1].*((1+mK./c[:,:,2,Nc]).^2).*((1+mNa./c[:,:,1,2]).^3)./2

    #compute glial sodium current and leak permeability value
    Ipumpg = Imaxg./(((1+mK./c[:,:,2,Nc]).^2).*((1+mNa./c[:,:,1,2]).^3))
    global const pNaLeakg = (-NaKCl-3*Ipumpg)./(log.(c[:,:,1,2]./c[:,:,1,3])+vmg)


    #Compute resting organic anion amounts and average valences
    al = zeros(size(phi))
    al[:,:,1:Nc-1] = alpha
    al[:,:,Nc] = 1-sum(alpha,3)
    aot = zeros(size(phi))
    zot = copy(aot)
    #set extracellular organic anion amounts and valence to ensure electroneutrality
    aot[:,:,Nc] = 5e-4
    cmphi = zeros(size(phi))
    for k=1:Nc-1
      cmphi[:,:,k] = cm[k]*(phi[:,:,k]-phi[:,:,Nc])
      cmphi[:,:,Nc] += cmphi[:,:,k]
      #set intracellular organic anion amounts to ensure osmotic pressure balance
      aot[:,:,k] = al[:,:,k].*(aot[:,:,Nc]./al[:,:,Nc]+sum(c[:,:,:,Nc],3)-sum(c[:,:,:,k],3))
      #set average valence to ensure electroneutrality
      zot[:,:,k] = (-cz(c[:,:,:,k],z).*al[:,:,k]+cmphi[:,:,k])./aot[:,:,k]
    end
    zot[:,:,Nc] = (-cz(c[:,:,:,Nc],z).*al[:,:,Nc]-cmphi[:,:,Nc])./aot[:,:,Nc]
    if spatially_uniform #return to Nx by Ny instead of 1 by 1 if necessary - check whether repmat will do this more efficiently
      phitemp = copy(phi)
      ctemp = copy(c)
      alphatemp = copy(alpha)
      aotemp = copy(aot)
      zotemp = copy(zot)
      phi = zeros(Nx,Ny,Nc)
      c = zeros(Nx,Ny,Ni,Nc)
      alpha = zeros(Nx,Ny,Nc-1)
      aot = zeros(Nx,Ny,Nc)
      zot = zeros(Nx,Ny,Nc)
      for k = 1:Nc
        phi[:,:,k] = phitemp[1,1,k]
        aot[:,:,k] = aotemp[1,1,k]
        zot[:,:,k] = zotemp[1,1,k]
        for i=1:Ni
          c[:,:,i,k] = ctemp[1,1,i,k]
        end
      end
      for k=1:Nc-1
        alpha[:,:,k] = alphatemp[1,1,k]
      end
    end
    global const ao = aot
    global const zo = zot
  else
    #cotransporter parameters from Bennett 2008
    global const pNaKCl = 0.002*RTFC/FC      #Bennett: .002 in mS/cm^2 converted to mmol/cm^2/s
    global const Imax = 1.3e1/FC           #maximum pump current in muA/cm^2 converted to mmol/cm^2/s from Yao, Huang, Miura, 2011
    global const Imaxg = 1.3e1/FC
    global const pNaLeak = 2e-2*RTFC/FC    #Kager:1e-2,Miura:2e-2%Na Leak conductance in mS/cm^2 converted to mmol/cm^2/s
    global const pNaLeakg = 2e-2*RTFC/FC
    global const zo = repmat([-1.0],Nc,1) #valence of organic anion
    global const ao = repmat([5e-4],Nx,Ny,Nc)
    println("If you get here, some of the above constants need to be made into matrices")
  end
  if fluid_flow
  #hydraulic permeability constants from Basser 1992
    kappan = 5e-9 # 5e-9 is grey matter hydraulic permeability in cm^4/(dynes sec)
    #white matter hydraulic permeability 7.5e-9 cm^4/(dynes sec)
    kappag = 5e-9 #glial compartment value
    kappae = 5e-9
    kappan *= 1e-2 #conversion to cm^2/sec/mPa
    kappag *= 1e-2
    kappae *= 1e-2
    kappat = zeros(Nx,Ny,Nc)
    kappat[:,:,1] = kappan
    kappat[:,:,2] = kappag
    kappat[:,:,3] = kappae
    global const kappa = kappat*R*T #conversion to cm^2/sec/(mmol/cm^3)
  else
  #Set kappa to 0 for no flow
    global const kappa = zeros(Nx,Ny,Nc)
  end

  #parameters for osmotic water flow
  #based on B.E. Shapiro dissertation (2000)
  zetat = 5.4e-5    #hydraulic permeability in cm/sec/(mmol/cm^3)
  zetat /= ell  #conversion to 1/sec/(mmol/cm^3)
  #     %based on Strieter, Stephenson, Palmer,
  #     %Weinstein, Journal or General Physiology, 1990.
  #     zeta=7e-8%6e-10%hydraulic permeability in cm/sec/mmHg
  #     zeta=zeta*7.501e-6%conversion to cm/sec/mPa
  #     zeta=zeta*R*T%conversion to cm/sec/(mmol/cm^3)
  #     zeta=zeta/ell%conversion to 1/sec/(mmol/cm^3)
  zetat *= ones(Nx,Ny,Nc-1)
  zetaadjust = 1                       #parameter for varying glial hydraulic permeability
  zetat[:,:,2] *= zetaadjust #adjust glial hydraulic permeability
  global const zeta1 = zetat
  if fluid_flow     #if spatial water flow is possible
    zetaalphat = 0  #stiffness constant or 1/stiffness constant
    global const S = true  #Indicates whether zetaalpha is the stiffness (true) or 1/stiffness (false)
  else #constants may not be changed for no spatial water flow
    zetaalphat = 0  #stiffness constant or 1/stiffness constant
    global const S = true  #Indicates whether zetaalpha is the stiffness (true) or 1/stiffness (false)
  end
  global const zetaalpha = zetaalphat*ones(Nx,Ny,Nc-1)

  #set flow rate and pressure to initial values
  if fluid_flow
    u = zeros(Nx,Ny,Nc)
    #pressure in units of mol/l (mPa/R/T)
    p = zeros(Nx,Ny,Nc)
    #set intracellular pressure using stiffness constant
    if S
      p[:,:,1:Nc-1] = broadcast(+,zetaalpha.*(alpha-alphao),p[:,:,Nc:Nc])
    else
      p[:,:,1:Nc-1]=broadcast(+,1./zetaalpha.*(alpha-alphao),p[:,:,Nc:Nc])
    end
  else
    u = zeros(Nx,Ny,Nc)
    #pressure in units of mol/l (mPa/R/T)
    p = zeros(Nx,Ny,Nc)
  end

  global const pzero = 0
  global const pnx = 0
  if !fluid_flow
    (state_vars,con_vars)=[SimState(c,phi,alpha),ConstVars(pNaKCl,Imax,pNaLeak,Imaxg,pNaLeakg,ao,zo,kappa,zeta1,S,zetaalpha,pzero,pnx)]
  else
    (state_vars,con_vars)=[SimStateF(c,phi,alpha,p,u),ConstVars(pNaKCl,Imax,pNaLeak,Imaxg,pNaLeakg,ao,zo,kappa,zeta1,S,zetaalpha,u,p,pzero,pnx)]
  end
end
*/