#include "constants.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void mclin(struct FluxData *flux,PetscInt index,PetscReal pc,PetscInt zi,PetscReal ci,PetscReal ce,PetscReal phim,PetscInt ADD)
{
	//Returns the flux value by ref.
	// pc is the permeativity, zi is valence, ci/e intra/extra concentration
	//phim is membrane voltage, index is the index in the flux struct
	//compute value and derivatives of the function:
    //mflux=pc.*(log(ci./ce)+z*phim)
    //for trans-membrane flux of an ion obeying a linear current-voltage equation
    if(ADD) { //If add, we accumulate the result
    	flux->mflux[index] += pc*(log(ci/ce)+zi*phim);
	    flux->dfdci[index] += pc/ci;
		flux->dfdce[index] += -pc/ce;
	  	flux->dfdphim[index] += zi*pc;
    } else{ //If not add we reninitialize
		flux->mflux[index] = pc*(log(ci/ce)+zi*phim);
	    flux->dfdci[index] = pc/ci;
		flux->dfdce[index] = -pc/ce;
	  	flux->dfdphim[index] = zi*pc;
  	}
}
void mcGoldman(struct FluxData *flux,PetscInt index,PetscReal pc,PetscInt zi,PetscReal ci,PetscReal ce,PetscReal phim,PetscInt ADD)
{
    //compute value and derivatives of the function:
    //mflux=p.*(z*phim).*(ci.*exp(z*phim)-ce)./(exp(z*phim)-1)
    //for trans-membrane flux of an ion obeying the GHK equations
    PetscReal xi = zi*phim/2;
    PetscReal r = exp(xi);

    //compute s=x/sinh(x)
    //Watch out for division by zero
    PetscReal s;
    if(xi!=0) {
    	s = xi/sinh(xi);
    } else {
    	s = 1;
    }
    //compute dfdci,dfdce,mflux
    PetscReal dfdci = pc*s*r;
    PetscReal dfdce = -pc*s/r;
    PetscReal mi = ci*dfdci;
    PetscReal me = ce*dfdce;
    PetscReal mflux = mi+me;
    //compute w=(sinh(x)/x-cosh(x))/x
    PetscReal w;
    if(fabs(xi)<0.2) //use Taylor expan. for small values
    {
    	//w0 = -x0.^9/3991680-x0.^7/45360-x0.^5/840-x0.^3/30-x0/3
    	w = -pow(xi,9)/3991680-pow(xi,7)/45360-pow(xi,5)/840-pow(xi,3)/30-xi/3;
    }
    else //use formula for larger values
    {
    	w = (sinh(xi)/xi-cosh(xi))/xi;
    }
    //compute dfdphim
    PetscReal dfdphim = (PetscReal)zi/2*(mflux*s*w+mi-me);
    if(ADD) //If ADD, we accumulate this result
    {
    	flux->mflux[index] += mflux;
    	flux->dfdce[index] += dfdce;
    	flux->dfdci[index] += dfdci;
    	flux->dfdphim[index] += dfdphim;
    }
    else // If not ADD, we reset the values
    {
		flux->mflux[index] = mflux;
    	flux->dfdce[index] = dfdce;
    	flux->dfdci[index] = dfdci;
    	flux->dfdphim[index] = dfdphim;
    }
    return;
}
PetscReal xoverexpminusone(PetscReal v,PetscReal aa,PetscReal bb,PetscReal cc,PetscInt dd)
{
	//computes aa*(v+bb)/(exp(cc*(v+bb))-1) if dd==0
 	//computes aa*(v+bb)/(1-exp(-cc*(v+bb)) otherwise
  	//for computing gating variables
  	v+=bb;
  	if(v==0) {
  		return aa/cc;
  	}
  	if(dd==0) {
  		return aa*v/(2*sinh(cc/2*v))*exp(-cc/2*v);
  	}
  	else {
  		return aa*v/(2*sinh(cc/2*v))*exp(cc/2*v);
  	}
}
PetscReal inwardrect(PetscReal ci,PetscReal ce,PetscReal phim)
{
	//inwardrect determines the effective conductance for the inward-rectifying
  	//potassium channel - formula and constants from Steinberg et al 2005
  	PetscReal Enernst = RTFC*log(ce/ci);
  	PetscReal EKdef = -85.2; //#-85.2 mV is base reversal potential
  	PetscReal cKo = .003; //#3 mM is base extracellular potassium concentration
  	return sqrt(ce/cKo)*(1+exp(18.5/42.5))/(1+exp((RTFC*phim-Enernst+18.5)/42.5))*(1+exp((-118.6+EKdef)/44.1))/(1+exp((-118.6+RTFC*phim)/44.1));
}
PetscReal cz(const PetscReal *cmat,const PetscInt *zvals,PetscInt x,PetscInt y,PetscInt Nx,PetscInt comp,struct AppCtx *user)
{
	//function to compute sum over i of c_i*z_i
	PetscReal accumulate=0;
	for(PetscInt ion=0;ion<Ni;ion++) {
		accumulate += zvals[ion]*cmat[c_index(x,y,comp,ion,Nx)];
	}
	return accumulate;
}
void diff_coef(PetscReal *Dc,const PetscReal *alp,PetscReal scale,struct AppCtx* user)
{
  //diffusion coefficients at all points, for all ions, in all compartments, in both x and y directions
	PetscReal tortuosity=1.6;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
	PetscReal alNcL,alNcR,alNcU;
  	for(PetscInt x=0;x<Nx;x++) {
	  	for(PetscInt y=0;y<Ny;y++) {
	  		alNcL=1-alp[al_index(x,y,0,Nx)]-alp[al_index(x,y,1,Nx)]; //Left extracell
			alNcR = 0;
			if(x<Nx-1) {
				alNcR = 1 - alp[al_index(x + 1, y, 0,Nx)] - alp[al_index(x + 1, y, 1,Nx)]; //Right extracell
			}
            alNcU = 0;
			if(y<Ny-1) {
				alNcU = 1 - alp[al_index(x, y + 1, 0,Nx)] - alp[al_index(x, y + 1, 1,Nx)];
			}
		  	for(PetscInt ion = 0; ion<Ni;ion++) {
			    //diffusion coefficients in x direction
			    if(x==(Nx-1)) {
			    	//Boundary is zero
			    	Dc[c_index(x,y,Nc-1,ion,Nx)*2] = scale*D[ion]*(alNcL)/(tortuosity*tortuosity);
			    }
			    else {
					//diffusion coefficients in the extracellular space proportional to volume fraction
			    	Dc[c_index(x,y,Nc-1,ion,Nx)*2] = scale*D[ion]/2*(alNcL+alNcR)/(tortuosity*tortuosity);
			    }
			    //diffusion coefficients in neuronal compartment equal to 0
			    Dc[c_index(x,y,0,ion,Nx)*2] = 0.0 ;
			    //diffusion coefficients in glial compartment proportional to default volume fraction
			    Dc[c_index(x,y,1,ion,Nx)*2] = 0*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2);
//                Dc[c_index(x,y,1,ion,Nx)*2] = 0.25*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2); //0.25
			    //diffusion coefficients in y direction
			    if(y==(Ny-1)) {
			    	//Boundary is zero
			    	Dc[c_index(x,y,Nc-1,ion,Nx)*2+1] = scale*D[ion]*(alNcL)/pow(tortuosity,2);
			    }
			    else {
			    	//diffusion coefficients in the extracellular space proportional to volume fraction
			    	Dc[c_index(x,y,Nc-1,ion,Nx)*2+1] = scale*D[ion]/2*(alNcL+alNcU)/pow(tortuosity,2);

				}
				//diffusion coefficients in neuronal compartment equal to 0
			    Dc[c_index(x,y,0,ion,Nx)*2+1] = 0.0;
			    //diffusion coefficients in glial compartment proportional to default volume fraction
//			    Dc[c_index(x,y,1,ion,Nx)*2+1] = 0.25*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2); //0.25
                Dc[c_index(x,y,1,ion,Nx)*2+1] = 0.001*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2);

		  	}
		}
	}
}

void gatevars_update(struct GateType *gate_vars,struct GateType *gate_vars_past, struct SimState *state_vars,PetscReal dtms,struct AppCtx *user,PetscInt firstpass)
{
	if(Profiling_on) {
		PetscLogEventBegin(event[4], 0, 0, 0, 0);
	}
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
	if(firstpass) {
		//membrane potential in mV
		PetscReal v = (state_vars->phi[phi_index(0,0,0,Nx)]-state_vars->phi[phi_index(0,0,Nc-1,Nx)])*RTFC;
		//Iniitialize the poPetscInt gating variables
		//Cause we assume a uniform start
		PetscReal alpha,beta;

		//compute current NaT
  		//gating variables mNaT
  		alpha = xoverexpminusone(v,0.32,51.9,0.25,1); //0.32*(Vm+51.9)./(1-exp(-0.25*(Vm+51.9)))
  		beta = xoverexpminusone(v,0.28,24.89,0.2,0); //0.28*(Vm+24.89)./(exp(0.2*(Vm+24.89))-1)
    	gate_vars->mNaT[0] = alpha/(alpha+beta);


  		//gating variable hNaT
  		alpha = 0.128*exp(-(0.056*v+2.94));
  		beta = 4/(exp(-(0.2*v+6))+1);
    	gate_vars->hNaT[0] = alpha/(alpha+beta);

    	gate_vars->gNaT[0] = pow(gate_vars->mNaT[0],3)*gate_vars->hNaT[0];

  		//compute current NaP
  		//gating variable mNaP
  		alpha = 1/(1+exp(-(0.143*v+5.67)))/6;
  		beta = 1.0/6-alpha; //1./(1+exp(0.143*Vm+5.67))/6
    	gate_vars->mNaP[0] = alpha/(alpha+beta);

  		//gating variable hNaP
  		alpha = 5.12e-6*exp(-(0.056*v+2.94));
  		beta = 1.6e-4/(1+exp(-(0.2*v+8)));
    	gate_vars->hNaP[0] = alpha/(alpha+beta);

  		gate_vars->gNaP[0] = pow(gate_vars->mNaP[0],2)*gate_vars->hNaP[0];
  		//compute KDR current
  		//gating variable mKDR
  		alpha = xoverexpminusone(v,0.016,34.9,0.2,1); //0.016*(Vm+34.9)./(1-exp(-0.2*(Vm+34.9)))
  		beta = 0.25*exp(-(0.025*v+1.25));
    	gate_vars->mKDR[0] = alpha/(alpha+beta);

  		gate_vars->gKDR[0] = pow(gate_vars->mKDR[0],2);

  		//compute KA current
  		//gating variable mKA
  		alpha = xoverexpminusone(v,0.02,56.9,0.1,1); //0.02*(Vm+56.9)./(1-exp(-0.1*(Vm+56.9)))
  		beta = xoverexpminusone(v,0.0175,29.9,0.1,0); //0.0175*(Vm+29.9)./(exp(0.1*(Vm+29.9))-1)
    	gate_vars->mKA[0] = alpha/(alpha+beta);

 		//gating variable hKA
  		alpha = 0.016*exp(-(0.056*v+4.61));
  		beta = 0.5/(exp(-(0.2*v+11.98))+1);
    	gate_vars->hKA[0] = alpha/(alpha+beta);

  		gate_vars->gKA[0] = pow(gate_vars->mKA[0],2)*gate_vars->hKA[0];

  		//Copy them over the remaining Nx by Ny points.
  		for(PetscInt i=0;i<Nx*Ny;i++) {
  			gate_vars->mNaT[i] = gate_vars->mNaT[0];
  			gate_vars->hNaT[i] = gate_vars->hNaT[0];
  			gate_vars->gNaT[i] = gate_vars->gNaT[0];

  			gate_vars->mNaP[i] = gate_vars->mNaP[0];
  			gate_vars->hNaP[i] = gate_vars->hNaP[0];
  			gate_vars->gNaP[i] = gate_vars->gNaP[0];

  			gate_vars->mKDR[i] = gate_vars->mKDR[0];
  			gate_vars->gKDR[i] = gate_vars->gKDR[0];

  			gate_vars->mKA[i] = gate_vars->mKA[0];
  			gate_vars->hKA[i] = gate_vars->hKA[0];
  			gate_vars->gKA[i] = gate_vars->gKA[0];

  		}
    } else { //if it's not the firstpass, then we actually have values in v.
		PetscReal v, alpha,beta;
		for(PetscInt x=0;x<Nx;x++) {
			for(PetscInt y=0;y<Ny;y++) {
				//membrane potential in mV
				v = (state_vars->phi[phi_index(x,y,0,Nx)]-state_vars->phi[phi_index(x,y,Nc-1,Nx)])*RTFC;

				//compute current NaT
		  		//gating variables mNaT
		  		alpha = xoverexpminusone(v,0.32,51.9,0.25,1); //0.32*(Vm+51.9)./(1-exp(-0.25*(Vm+51.9)))
		  		beta = xoverexpminusone(v,0.28,24.89,0.2,0); //0.28*(Vm+24.89)./(exp(0.2*(Vm+24.89))-1)
		    	gate_vars->mNaT[xy_index(x,y,Nx)] = (gate_vars_past->mNaT[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		//gating variable hNaT
		  		alpha = 0.128*exp(-(0.056*v+2.94));
		  		beta = 4/(exp(-(0.2*v+6))+1);
		    	gate_vars->hNaT[xy_index(x,y,Nx)] = alpha/(alpha+beta);
		    	gate_vars->hNaT[xy_index(x,y,Nx)] = (gate_vars_past->hNaT[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		    	gate_vars->gNaT[xy_index(x,y,Nx)] = pow(gate_vars->mNaT[xy_index(x,y,Nx)],3)*gate_vars->hNaT[xy_index(x,y,Nx)];
		  		//compute current NaP
		  		//gating variable mNaP
		  		alpha = 1/(1+exp(-(0.143*v+5.67)))/6;
		  		beta = 1.0/6.0-alpha; //1./(1+exp(0.143*Vm+5.67))/6
		  		gate_vars->mNaP[xy_index(x,y,Nx)] = (gate_vars_past->mNaP[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		//gating variable hNaP
		  		alpha = 5.12e-6*exp(-(0.056*v+2.94));
		  		beta = 1.6e-4/(1+exp(-(0.2*v+8)));
		  		gate_vars->hNaP[xy_index(x,y,Nx)] = (gate_vars_past->hNaP[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gNaP[xy_index(x,y,Nx)] = pow(gate_vars->mNaP[xy_index(x,y,Nx)],2)*gate_vars->hNaP[xy_index(x,y,Nx)];

		  		//compute KDR current
		  		//gating variable mKDR
		  		alpha = xoverexpminusone(v,0.016,34.9,0.2,1); //0.016*(Vm+34.9)./(1-exp(-0.2*(Vm+34.9)))
		  		beta = 0.25*exp(-(0.025*v+1.25));
		    	gate_vars->mKDR[xy_index(x,y,Nx)] = (gate_vars_past->mKDR[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gKDR[xy_index(x,y,Nx)] = pow(gate_vars->mKDR[xy_index(x,y,Nx)],2);

		  		//compute KA current
		  		//gating variable mKA
		  		alpha = xoverexpminusone(v,0.02,56.9,0.1,1); //0.02*(Vm+56.9)./(1-exp(-0.1*(Vm+56.9)))
		  		beta = xoverexpminusone(v,0.0175,29.9,0.1,0); //0.0175*(Vm+29.9)./(exp(0.1*(Vm+29.9))-1)
		    	gate_vars->mKA[xy_index(x,y,Nx)] = (gate_vars_past->mKA[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		 		//gating variable hKA
		  		alpha = 0.016*exp(-(0.056*v+4.61));
		  		beta = 0.5/(exp(-(0.2*v+11.98))+1);
		    	gate_vars->hKA[xy_index(x,y,Nx)] = (gate_vars_past->hKA[xy_index(x,y,Nx)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gKA[xy_index(x,y,Nx)] = pow(gate_vars->mKA[xy_index(x,y,Nx)],2)*gate_vars->hKA[xy_index(x,y,Nx)];
  			}
  		}
  	}
	if(Profiling_on) {
		PetscLogEventEnd(event[4], 0, 0, 0, 0);
	}
}

void excitation(struct AppCtx* user,PetscReal t)
{
  //compute excitation conductance to trigger csd
  //Leak conductances in mS/cm^2
  //all units converted to mmol/cm^2/sec
  	PetscReal pexct,pany;
  	PetscReal xexct;
    PetscReal radius;
    struct ExctType *exct = user->gexct;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt num_points = 0;
    if(one_point_exct){
        for (PetscInt i = 0; i < Nx; i++) {
            for (PetscInt j = 0; j < Ny; j++) {

                if (t < texct && i==0 && j==0) {
                    num_points++;
                    pany = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
                    exct->pNa[xy_index(i, j,Nx)] = pany;
                    exct->pK[xy_index(i, j,Nx)] = pany;
                    exct->pCl[xy_index(i, j,Nx)] = pany;
                }else {
                    //pexct=0*RTFC/FC
                    exct->pNa[xy_index(i, j,Nx)] = 0;
                    exct->pK[xy_index(i, j,Nx)] = 0;
                    exct->pCl[xy_index(i, j,Nx)] = 0;
                }
            }
        }

    }else {
        for (PetscInt i = 0; i < Nx; i++) {
            for (PetscInt j = 0; j < Ny; j++) {
                if (mid_points_exct) {
                    radius = sqrt(pow((i + 0.5) * dx - Lx / 2, 2) + pow((j + 0.5) * dy - Lx / 2, 2));
                    if (t < texct && radius < Lexct) {
                        num_points++;
                        pexct = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
                        xexct = pow(cos(pi / 2 * (radius / Lexct)), 2);
                        pany = pexct * xexct;
                        exct->pNa[xy_index(i, j,Nx)] = pany;
                        exct->pK[xy_index(i, j,Nx)] = pany;
                        exct->pCl[xy_index(i, j,Nx)] = pany;
                    } else {
                        //pexct=0*RTFC/FC
                        exct->pNa[xy_index(i, j,Nx)] = 0;
                        exct->pK[xy_index(i, j,Nx)] = 0;
                        exct->pCl[xy_index(i, j,Nx)] = 0;
                    }
                } else {
                    radius = sqrt(pow((i + 0.5) * dx, 2) + pow((j + 0.5) * dy, 2));
                    if (t < texct && radius < Lexct) {
                        num_points++;
                        pexct = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
//	    		xexct=pow((cos(pi/2*(i+.5)/Nexct))*(cos(pi/2*(j+.5)/Nexct)),2);
                        xexct = pow(cos(pi / 2 * (radius / Lexct)), 2);
//				xexct=pow((cos(pi/2*((i+.5)*dx)/Lexct))*(cos(pi/2*((j+.5)*dy)/Lexct)),2);
                        pany = pexct * xexct;
                        exct->pNa[xy_index(i, j,Nx)] = pany;
                        exct->pK[xy_index(i, j,Nx)] = pany;
                        exct->pCl[xy_index(i, j,Nx)] = pany;
                    } else {
                        //pexct=0*RTFC/FC
                        exct->pNa[xy_index(i, j,Nx)] = 0;
                        exct->pK[xy_index(i, j,Nx)] = 0;
                        exct->pCl[xy_index(i, j,Nx)] = 0;
                    }
                }
            }
        }
    }
//    printf("Number of excited points: %d\n",num_points);
}

void ionmflux(struct AppCtx* user)
{
	if(Profiling_on) {
		PetscLogEventBegin(event[5], 0, 0, 0, 0);
	}
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    struct FluxData *flux = user->flux;
    struct SimState *state_vars=user->state_vars;
    struct SimState *state_vars_past = user->state_vars_past;
    struct GateType *gvars = user->gate_vars_past;
    struct ExctType *gexct = user->gexct;
    struct ConstVars *con_vars = user->con_vars;
    //Variables to save to for ease of notation
    PetscReal vm,vmg,vmgp;
    PetscReal ci,cg,ce,cgp,cep;
    PetscReal Na,K;//Variables for pump (so it's clear)

    //For calculationg permeabilities
    PetscReal pGHK,pLin;
    PetscReal Ipump,NaKCl;
    for(PetscInt x=0;x<Nx;x++) {
        for(PetscInt y=0;y<Ny;y++) {
            vm = state_vars->phi[phi_index(x,y,0,Nx)]-state_vars->phi[phi_index(x,y,Nc-1,Nx)];
            vmg = state_vars->phi[phi_index(x,y,1,Nx)]-state_vars->phi[phi_index(x,y,Nc-1,Nx)];
            vmgp = state_vars_past->phi[phi_index(x,y,1,Nx)]-state_vars_past->phi[phi_index(x,y,Nc-1,Nx)];

            //Compute Na Channel currents
            ci = state_vars->c[c_index(x,y,0,0,Nx)];
            cg = state_vars->c[c_index(x,y,1,0,Nx)];
            ce = state_vars->c[c_index(x,y,Nc-1,0,Nx)];

            //Neurons
            pGHK = pNaT*gvars->gNaT[xy_index(x,y,Nx)]+pNaP*gvars->gNaP[xy_index(x,y,Nx)];
            pLin = con_vars->pNaLeak + gexct->pNa[xy_index(x,y,Nx)]; //Add excitation
            //Initialize GHK Flux
            mcGoldman(flux,c_index(x,y,0,0,Nx),pGHK,1,ci,ce,vm,0);
            //Add leak current to that.
            mclin(flux,c_index(x,y,0,0,Nx),pLin,1,ci,ce,vm,1);
            //Glial NaLeak
            mclin(flux,c_index(x,y,1,0,Nx),con_vars->pNaLeakg,1,cg,ce,vmg,0);

            // Compute K channel Currents
            ci = state_vars->c[c_index(x,y,0,1,Nx)];
            cg = state_vars->c[c_index(x,y,1,1,Nx)];
            ce = state_vars->c[c_index(x,y,Nc-1,1,Nx)];

            //Neurons
            pGHK = pKDR*gvars->gKDR[xy_index(x,y,Nx)]+pKA*gvars->gKA[xy_index(x,y,Nx)];
            pLin = pKLeak+gexct->pK[xy_index(x,y,Nx)]; //add excitation
            mcGoldman(flux,c_index(x,y,0,1,Nx),pGHK,1,ci,ce,vm,0);
            mclin(flux,c_index(x,y,0,1,Nx),pLin,1,ci,ce,vm,1);

            // Glial K Leak(using past)
            cgp = state_vars_past->c[c_index(x,y,1,1,Nx)];
            cep = state_vars_past->c[c_index(x,y,Nc-1,1,Nx)];

            pLin = pKIR*inwardrect(cgp,cep,vmgp)*pKLeakadjust;
            mclin(flux,c_index(x,y,1,1,Nx),pLin,1,cg,ce,vmg,0);

            //Compute Cl Channel Current
            ci = state_vars->c[c_index(x,y,0,2,Nx)];
            cg = state_vars->c[c_index(x,y,1,2,Nx)];
            ce = state_vars->c[c_index(x,y,Nc-1,2,Nx)];

            //Neurons
            pLin = pClLeak+gexct->pCl[xy_index(x,y,Nx)]; //add excitation
            mclin(flux,c_index(x,y,0,2,Nx),pLin,-1,ci,ce,vm,0);

            //Glia
            mclin(flux,c_index(x,y,1,2,Nx),pClLeakg,-1,cg,ce,vmg,0);

            //Pump Currents(all past values)

            //Neurons
            Na = state_vars_past->c[c_index(x,y,0,0,Nx)];
            K = state_vars_past->c[c_index(x,y,Nc-1,1,Nx)];

            Ipump = npump*con_vars->Imax/(pow(1+mK/K,2)*pow(1+mNa/Na,3));

            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x,y,0,0,Nx)]+=3*Ipump; //Na part
            flux->mflux[c_index(x,y,0,1,Nx)]-=2*Ipump; //K part

            //Glia
            Na = state_vars_past->c[c_index(x,y,1,0,Nx)];
            //K is the same(extracellular)
            Ipump = glpump*con_vars->Imaxg/(pow(1+mK/K,2)*pow(1+mNa/Na,3));
            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x,y,1,0,Nx)]+=3*Ipump; //Na part
            flux->mflux[c_index(x,y,1,1,Nx)]-=2*Ipump; //K part

            //NaKCl Cotransporter
            //I'm going to reuse variables names...
            Na = state_vars_past->c[c_index(x,y,1,0,Nx)];//glia Na
            K = state_vars_past->c[c_index(x,y,1,1,Nx)]; // glia K.
            cgp = state_vars_past->c[c_index(x,y,1,2,Nx)]; //glia Cl

            cep = state_vars_past->c[c_index(x,y,Nc-1,0,Nx)];//Ext Na
            ce = state_vars_past->c[c_index(x,y,Nc-1,1,Nx)]; // Ext K.
            ci = state_vars_past->c[c_index(x,y,Nc-1,2,Nx)]; //Ext Cl

            NaKCl = con_vars->pNaKCl*log(Na*K*pow(cgp,2)/(cep*ce*pow(ci,2)));
            //Add to flux
            flux->mflux[c_index(x,y,1,0,Nx)]+=NaKCl; //Sodium
            flux->mflux[c_index(x,y,1,1,Nx)]+=NaKCl; //K
            flux->mflux[c_index(x,y,1,2,Nx)]+=2*NaKCl; //Cl

            //Change units of flux from mmol/cm^2 to mmol/cm^3/s
            for(PetscInt ion=0;ion<Ni;ion++) {
                flux->mflux[c_index(x,y,Nc-1,ion,Nx)] = 0;
                for(PetscInt comp = 0;comp<Nc-1;comp++) {
                    flux->mflux[c_index(x,y,comp,ion,Nx)]=flux->mflux[c_index(x,y,comp,ion,Nx)]/ell;
                    flux->dfdci[c_index(x,y,comp,ion,Nx)]=flux->dfdci[c_index(x,y,comp,ion,Nx)]/ell;
                    flux->dfdce[c_index(x,y,comp,ion,Nx)]=flux->dfdce[c_index(x,y,comp,ion,Nx)]/ell;
                    flux->dfdphim[c_index(x,y,comp,ion,Nx)]=flux->dfdphim[c_index(x,y,comp,ion,Nx)]/ell;

                    //And calculate the extracellular flux
                    flux->mflux[c_index(x,y,Nc-1,ion,Nx)] -= flux->mflux[c_index(x,y,comp,ion,Nx)];
                }
            }
        }
    }
	if(Profiling_on) {
		PetscLogEventEnd(event[5], 0, 0, 0, 0);
	}
}
void wflowm(struct AppCtx *user)
{
  //piw = sum of c over ions + ao/alpha
  // piw is the total number of ions in a compartment
  //outward transmembrane water flow seen as a function of
  //osmotic pressure and volume fraction or pressure.
	if(Profiling_on) {
		PetscLogEventBegin(event[6], 0, 0, 0, 0);
	}
    struct FluxData *flux = user->flux;
    struct SimState *state_vars = user->state_vars;
    struct ConstVars *con_vars = user->con_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;

    PetscReal dwdpi,dwdal,piw,piwNc;
    for(PetscInt x=0;x<Nx;x++) {
        for(PetscInt y=0;y<Ny;y++) {
            //Calculate the pi for extracellular
            piwNc = 0;
            for(PetscInt ion=0;ion<Ni;ion++) {
                piwNc +=state_vars->c[c_index(x,y,Nc-1,ion,Nx)];
            }
            piwNc +=con_vars->ao[al_index(0,0,Nc-1,Nx)]/(1-state_vars->alpha[al_index(x,y,0,Nx)]-state_vars->alpha[al_index(x,y,1,Nx)]);
            for(PetscInt comp = 0;comp<Nc-1;comp++) {
                piw = 0;
                for(PetscInt ion=0;ion<Ni;ion++) {
                    piw +=state_vars->c[c_index(x,y,comp,ion,Nx)];
                }
                piw +=con_vars->ao[al_index(0,0,comp,Nx)]/state_vars->alpha[al_index(x,y,comp,Nx)];
                //ao, zeta1, and zetalpha are currently constant in space
                dwdpi = con_vars->zeta1[al_index(0,0,comp,Nx)];
                dwdal = con_vars->zeta1[al_index(0,0,comp,Nx)]*con_vars->zetaalpha[al_index(0,0,comp,Nx)];

                flux->wflow[al_index(x,y,comp,Nx)] = dwdpi*(piwNc-piw)+dwdal*(state_vars->alpha[al_index(x,y,comp,Nx)]-alphao[comp]);
                flux->dwdpi[al_index(x,y,comp,Nx)] = dwdpi;
                flux->dwdal[al_index(x,y,comp,Nx)] = dwdal;
            }
        }
    }
	if(Profiling_on) {
		PetscLogEventEnd(event[6], 0, 0, 0, 0);
	}
}

void grid_wflowm(struct AppCtx *user)
{
    //piw = sum of c over ions + ao/alpha
    // piw is the total number of ions in a compartment
    //outward transmembrane water flow seen as a function of
    //osmotic pressure and volume fraction or pressure.
    if(Profiling_on) {
        PetscLogEventBegin(event[6], 0, 0, 0, 0);
    }
    struct FluxData *flux = user->flux;
    struct SimState *state_vars = user->grid_vars;
    struct ConstVars *con_vars = user->con_vars;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;

    PetscReal dwdpi,dwdal,piw,piwNc;
    for(PetscInt x=0;x<Nx;x++) {
        for (PetscInt y = 0; y < Ny; y++) {

            //Calculate the pi for extracellular
            piwNc = 0;
            for (PetscInt ion = 0; ion < Ni; ion++) {
                piwNc += state_vars->c[c_index(x, y, Nc - 1, ion, Nx)];
            }
            piwNc += con_vars->ao[al_index(0, 0, Nc - 1, Nx)] /
                     (1 - state_vars->alpha[al_index(x, y, 0, Nx)] - state_vars->alpha[al_index(x, y, 1, Nx)]);
            for (PetscInt comp = 0; comp < Nc - 1; comp++) {
                piw = 0;
                for (PetscInt ion = 0; ion < Ni; ion++) {
                    piw += state_vars->c[c_index(x, y, comp, ion, Nx)];
                }
                piw += con_vars->ao[al_index(0, 0, comp, Nx)] / state_vars->alpha[al_index(x, y, comp, Nx)];
                //ao, zeta1, and zetalpha are currently constant in space
                dwdpi = con_vars->zeta1[al_index(0, 0, comp, Nx)];
                dwdal = con_vars->zeta1[al_index(0, 0, comp, Nx)] * con_vars->zetaalpha[al_index(0, 0, comp, Nx)];

                flux->wflow[al_index(x, y, comp, Nx)] =
                        dwdpi * (piwNc - piw) + dwdal * (state_vars->alpha[al_index(x, y, comp, Nx)] - alphao[comp]);
                flux->dwdpi[al_index(x, y, comp, Nx)] = dwdpi;
                flux->dwdal[al_index(x, y, comp, Nx)] = dwdal;
            }
        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[6], 0, 0, 0, 0);
    }
}
void grid_ionmflux(struct AppCtx* user)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[5], 0, 0, 0, 0);
    }
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    struct FluxData *flux = user->flux;
    struct SimState *state_vars=user->grid_vars;
    struct SimState *state_vars_past = user->grid_vars_past;
    struct GateType *gvars = user->grid_gate_vars;
    struct ExctType *gexct = user->gexct;
    struct ConstVars *con_vars = user->con_vars;
    //Variables to save to for ease of notation
    PetscReal vm,vmg,vmgp;
    PetscReal ci,cg,ce,cgp,cep;
    PetscReal Na,K;//Variables for pump (so it's clear)

    //For calculationg permeabilities
    PetscReal pGHK,pLin;
    PetscReal Ipump,NaKCl;
    for(PetscInt x=0;x<Nx;x++){
        for(PetscInt y=0;y<Ny;y++) {
            vm = state_vars->phi[phi_index(x, y, 0, Nx)] - state_vars->phi[phi_index(x, y, Nc - 1, Nx)];
            vmg = state_vars->phi[phi_index(x, y, 1, Nx)] - state_vars->phi[phi_index(x, y, Nc - 1, Nx)];
            vmgp = state_vars_past->phi[phi_index(x, y, 1, Nx)] - state_vars_past->phi[phi_index(x, y, Nc - 1, Nx)];

            //Compute Na Channel currents
            ci = state_vars->c[c_index(x, y, 0, 0, Nx)];
            cg = state_vars->c[c_index(x, y, 1, 0, Nx)];
            ce = state_vars->c[c_index(x, y, Nc - 1, 0, Nx)];

            //Neurons
            pGHK = pNaT * gvars->gNaT[xy_index(x, y, Nx)] + pNaP * gvars->gNaP[xy_index(x, y, Nx)];
            pLin = con_vars->pNaLeak + gexct->pNa[xy_index(x, y, Nx)]; //Add excitation
            //Initialize GHK Flux
            mcGoldman(flux, c_index(x, y, 0, 0, Nx), pGHK, 1, ci, ce, vm, 0);
            //Add leak current to that.
            mclin(flux, c_index(x, y, 0, 0, Nx), pLin, 1, ci, ce, vm, 1);
            //Glial NaLeak
            mclin(flux, c_index(x, y, 1, 0, Nx), con_vars->pNaLeakg, 1, cg, ce, vmg, 0);

            // Compute K channel Currents
            ci = state_vars->c[c_index(x, y, 0, 1, Nx)];
            cg = state_vars->c[c_index(x, y, 1, 1, Nx)];
            ce = state_vars->c[c_index(x, y, Nc - 1, 1, Nx)];

            //Neurons
            pGHK = pKDR * gvars->gKDR[xy_index(x, y, Nx)] + pKA * gvars->gKA[xy_index(x, y, Nx)];
            pLin = pKLeak + gexct->pK[xy_index(x, y, Nx)]; //add excitation
            mcGoldman(flux, c_index(x, y, 0, 1, Nx), pGHK, 1, ci, ce, vm, 0);
            mclin(flux, c_index(x, y, 0, 1, Nx), pLin, 1, ci, ce, vm, 1);

            // Glial K Leak(using past)
            cgp = state_vars_past->c[c_index(x, y, 1, 1, Nx)];
            cep = state_vars_past->c[c_index(x, y, Nc - 1, 1, Nx)];

            pLin = pKIR * inwardrect(cgp, cep, vmgp) * pKLeakadjust;
            mclin(flux, c_index(x, y, 1, 1, Nx), pLin, 1, cg, ce, vmg, 0);

            //Compute Cl Channel Current
            ci = state_vars->c[c_index(x, y, 0, 2, Nx)];
            cg = state_vars->c[c_index(x, y, 1, 2, Nx)];
            ce = state_vars->c[c_index(x, y, Nc - 1, 2, Nx)];

            //Neurons
            pLin = pClLeak + gexct->pCl[xy_index(x, y, Nx)]; //add excitation
            mclin(flux, c_index(x, y, 0, 2, Nx), pLin, -1, ci, ce, vm, 0);

            //Glia
            mclin(flux, c_index(x, y, 1, 2, Nx), pClLeakg, -1, cg, ce, vmg, 0);

            //Pump Currents(all past values)

            //Neurons
            Na = state_vars_past->c[c_index(x, y, 0, 0, Nx)];
            K = state_vars_past->c[c_index(x, y, Nc - 1, 1, Nx)];

            Ipump = npump * con_vars->Imax / (pow(1 + mK / K, 2) * pow(1 + mNa / Na, 3));

            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x, y, 0, 0, Nx)] += 3 * Ipump; //Na part
            flux->mflux[c_index(x, y, 0, 1, Nx)] -= 2 * Ipump; //K part

            //Glia
            Na = state_vars_past->c[c_index(x, y, 1, 0, Nx)];
            //K is the same(extracellular)
            Ipump = glpump * con_vars->Imaxg / (pow(1 + mK / K, 2) * pow(1 + mNa / Na, 3));
            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x, y, 1, 0, Nx)] += 3 * Ipump; //Na part
            flux->mflux[c_index(x, y, 1, 1, Nx)] -= 2 * Ipump; //K part

            //NaKCl Cotransporter
            //I'm going to reuse variables names...
            Na = state_vars_past->c[c_index(x, y, 1, 0, Nx)];//glia Na
            K = state_vars_past->c[c_index(x, y, 1, 1, Nx)]; // glia K.
            cgp = state_vars_past->c[c_index(x, y, 1, 2, Nx)]; //glia Cl

            cep = state_vars_past->c[c_index(x, y, Nc - 1, 0, Nx)];//Ext Na
            ce = state_vars_past->c[c_index(x, y, Nc - 1, 1, Nx)]; // Ext K.
            ci = state_vars_past->c[c_index(x, y, Nc - 1, 2, Nx)]; //Ext Cl

            NaKCl = con_vars->pNaKCl * log(Na * K * pow(cgp, 2) / (cep * ce * pow(ci, 2)));
            //Add to flux
            flux->mflux[c_index(x, y, 1, 0, Nx)] += NaKCl; //Sodium
            flux->mflux[c_index(x, y, 1, 1, Nx)] += NaKCl; //K
            flux->mflux[c_index(x, y, 1, 2, Nx)] += 2 * NaKCl; //Cl

            //Change units of flux from mmol/cm^2 to mmol/cm^3/s
            for (PetscInt ion = 0; ion < Ni; ion++) {
                flux->mflux[c_index(x, y, Nc - 1, ion, Nx)] = 0;
                for (PetscInt comp = 0; comp < Nc - 1; comp++) {
                    flux->mflux[c_index(x, y, comp, ion, Nx)] = flux->mflux[c_index(x, y, comp, ion, Nx)] / ell;
                    flux->dfdci[c_index(x, y, comp, ion, Nx)] = flux->dfdci[c_index(x, y, comp, ion, Nx)] / ell;
                    flux->dfdce[c_index(x, y, comp, ion, Nx)] = flux->dfdce[c_index(x, y, comp, ion, Nx)] / ell;
                    flux->dfdphim[c_index(x, y, comp, ion, Nx)] = flux->dfdphim[c_index(x, y, comp, ion, Nx)] / ell;

                    //And calculate the extracellular flux
                    flux->mflux[c_index(x, y, Nc - 1, ion, Nx)] -= flux->mflux[c_index(x, y, comp, ion, Nx)];
                }

            }
        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[5], 0, 0, 0, 0);
    }
}

void gatevars_update_grid(struct GateType *gate_vars,struct SimState *state_vars,PetscReal dtms,struct AppCtx *user)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[4], 0, 0, 0, 0);
    }
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;

    PetscReal v, alpha,beta;

    for(PetscInt x=0;x<Nx;x++) {
        for (PetscInt y = 0; y < Ny; y++) {
            //membrane potential in mV
            v = (state_vars->phi[phi_index(x, y, 0, Nx)] - state_vars->phi[phi_index(x, y, Nc - 1, Nx)]) * RTFC;

            //compute current NaT
            //gating variables mNaT
            alpha = xoverexpminusone(v, 0.32, 51.9, 0.25, 1); //0.32*(Vm+51.9)./(1-exp(-0.25*(Vm+51.9)))
            beta = xoverexpminusone(v, 0.28, 24.89, 0.2, 0); //0.28*(Vm+24.89)./(exp(0.2*(Vm+24.89))-1)
            gate_vars->mNaT[xy_index(x, y, Nx)] =
                    (gate_vars->mNaT[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            //gating variable hNaT
            alpha = 0.128 * exp(-(0.056 * v + 2.94));
            beta = 4 / (exp(-(0.2 * v + 6)) + 1);
            gate_vars->hNaT[xy_index(x, y, Nx)] = alpha / (alpha + beta);
            gate_vars->hNaT[xy_index(x, y, Nx)] =
                    (gate_vars->hNaT[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            gate_vars->gNaT[xy_index(x, y, Nx)] =
                    pow(gate_vars->mNaT[xy_index(x, y, Nx)], 3) * gate_vars->hNaT[xy_index(x, y, Nx)];
            //compute current NaP
            //gating variable mNaP
            alpha = 1 / (1 + exp(-(0.143 * v + 5.67))) / 6;
            beta = 1.0 / 6.0 - alpha; //1./(1+exp(0.143*Vm+5.67))/6
            gate_vars->mNaP[xy_index(x, y, Nx)] =
                    (gate_vars->mNaP[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            //gating variable hNaP
            alpha = 5.12e-6 * exp(-(0.056 * v + 2.94));
            beta = 1.6e-4 / (1 + exp(-(0.2 * v + 8)));
            gate_vars->hNaP[xy_index(x, y, Nx)] =
                    (gate_vars->hNaP[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            gate_vars->gNaP[xy_index(x, y, Nx)] =
                    pow(gate_vars->mNaP[xy_index(x, y, Nx)], 2) * gate_vars->hNaP[xy_index(x, y, Nx)];

            //compute KDR current
            //gating variable mKDR
            alpha = xoverexpminusone(v, 0.016, 34.9, 0.2, 1); //0.016*(Vm+34.9)./(1-exp(-0.2*(Vm+34.9)))
            beta = 0.25 * exp(-(0.025 * v + 1.25));
            gate_vars->mKDR[xy_index(x, y, Nx)] =
                    (gate_vars->mKDR[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            gate_vars->gKDR[xy_index(x, y, Nx)] = pow(gate_vars->mKDR[xy_index(x, y, Nx)], 2);

            //compute KA current
            //gating variable mKA
            alpha = xoverexpminusone(v, 0.02, 56.9, 0.1, 1); //0.02*(Vm+56.9)./(1-exp(-0.1*(Vm+56.9)))
            beta = xoverexpminusone(v, 0.0175, 29.9, 0.1, 0); //0.0175*(Vm+29.9)./(exp(0.1*(Vm+29.9))-1)
            gate_vars->mKA[xy_index(x, y, Nx)] =
                    (gate_vars->mKA[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            //gating variable hKA
            alpha = 0.016 * exp(-(0.056 * v + 4.61));
            beta = 0.5 / (exp(-(0.2 * v + 11.98)) + 1);
            gate_vars->hKA[xy_index(x, y, Nx)] =
                    (gate_vars->hKA[xy_index(x, y, Nx)] + alpha * dtms) / (1 + (alpha + beta) * dtms);

            gate_vars->gKA[xy_index(x, y, Nx)] =
                    pow(gate_vars->mKA[xy_index(x, y, Nx)], 2) * gate_vars->hKA[xy_index(x, y, Nx)];

        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[4], 0, 0, 0, 0);
    }
}
void excitation_grid(struct AppCtx* user,PetscReal t,PetscInt xi,PetscInt yi)
{
    //compute excitation conductance to trigger csd
    //Leak conductances in mS/cm^2
    //all units converted to mmol/cm^2/sec
    PetscReal pexct,pany;
    PetscReal xexct;
    PetscReal radius;
    struct ExctType *exct = user->gexct;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    for(PetscInt i=0;i<Nx;i++){
        for(PetscInt j=0;j<Ny;j++) {
            if (one_point_exct) {
                if (t < texct && i+xi == 0 && j+yi == 0) {
                    pany = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
                    exct->pNa[xy_index(i, j, Nx)] = pany;
                    exct->pK[xy_index(i, j, Nx)] = pany;
                    exct->pCl[xy_index(i, j, Nx)] = pany;
                } else {
                    //pexct=0*RTFC/FC
                    exct->pNa[xy_index(i, j, Nx)] = 0;
                    exct->pK[xy_index(i, j, Nx)] = 0;
                    exct->pCl[xy_index(i, j, Nx)] = 0;
                }
            } else {
                if (mid_points_exct) {
                    radius = sqrt(pow((i+xi + 0.5) * dx - Lx / 2, 2) + pow((j+yi + 0.5) * dy - Lx / 2, 2));
                    if (t < texct && radius < Lexct) {
                        pexct = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
                        xexct = pow(cos(pi / 2 * (radius / Lexct)), 2);
                        pany = pexct * xexct;
                        exct->pNa[xy_index(i, j, Nx)] = pany;
                        exct->pK[xy_index(i, j, Nx)] = pany;
                        exct->pCl[xy_index(i, j, Nx)] = pany;
                    } else {
                        //pexct=0*RTFC/FC
                        exct->pNa[xy_index(i, j, Nx)] = 0;
                        exct->pK[xy_index(i, j, Nx)] = 0;
                        exct->pCl[xy_index(i, j, Nx)] = 0;
                    }
                } else {
                    radius = sqrt(pow((i+xi + 0.5) * dx, 2) + pow((j+yi + 0.5) * dy, 2));
                    if (t < texct && radius < Lexct) {
                        pexct = pmax * pow(sin(pi * t / texct), 2) * RTFC / FC;
//	    		xexct=pow((cos(pi/2*(i+xi+.5)/Nexct))*(cos(pi/2*(j+yi+.5)/Nexct)),2);
                        xexct = pow(cos(pi / 2 * (radius / Lexct)), 2);
//				xexct=pow((cos(pi/2*((i+xi+.5)*dx)/Lexct))*(cos(pi/2*((j+yi+.5)*dy)/Lexct)),2);
                        pany = pexct * xexct;
                        exct->pNa[xy_index(i, j, Nx)] = pany;
                        exct->pK[xy_index(i, j, Nx)] = pany;
                        exct->pCl[xy_index(i, j, Nx)] = pany;
                    } else {
                        //pexct=0*RTFC/FC
                        exct->pNa[xy_index(i, j, Nx)] = 0;
                        exct->pK[xy_index(i, j, Nx)] = 0;
                        exct->pCl[xy_index(i, j, Nx)] = 0;
                    }
                }
            }
        }
    }

//    printf("Number of excited points: %d\n",num_points);
}
void grid_diff_coef(PetscReal *Dc,const PetscReal *alp,PetscReal scale,struct AppCtx* user)
{
    //diffusion coefficients at all points, for all ions, in all compartments, in both x and y directions
    PetscReal tortuosity=1.6;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscReal alNcL,alNcR,alNcU;
    for(PetscInt x=0;x<Nx;x++) {
        for(PetscInt y=0;y<Ny;y++) {
            alNcL=1-alp[al_index(x,y,0,Nx)]-alp[al_index(x,y,1,Nx)]; //Left extracell
            alNcR = 0;
            if(x<Nx-1) {
                alNcR = 1 - alp[al_index(x + 1, y, 0,Nx)] - alp[al_index(x + 1, y, 1,Nx)]; //Right extracell
            }
            alNcU = 0;
            if(y<Ny-1) {
                alNcU = 1 - alp[al_index(x, y + 1, 0,Nx)] - alp[al_index(x, y + 1, 1,Nx)];
            }
            for(PetscInt ion = 0; ion<Ni;ion++) {
                //diffusion coefficients in x direction
                if(x==(Nx-1)) {
                    //Boundary is zero
                    Dc[c_index(x,y,Nc-1,ion,Nx)*2] = scale*D[ion]*(alNcL)/(tortuosity*tortuosity);
                }
                else {
                    //diffusion coefficients in the extracellular space proportional to volume fraction
                    Dc[c_index(x,y,Nc-1,ion,Nx)*2] = scale*D[ion]/2*(alNcL+alNcR)/(tortuosity*tortuosity);
                }
                //diffusion coefficients in neuronal compartment equal to 0
                Dc[c_index(x,y,0,ion,Nx)*2] = 0.0 ;
                //diffusion coefficients in glial compartment proportional to default volume fraction
                Dc[c_index(x,y,1,ion,Nx)*2] = 0*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2);
//                Dc[c_index(x,y,1,ion,Nx)*2] = 0.25*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2); //0.25
                //diffusion coefficients in y direction
                if(y==(Ny-1)) {
                    //Boundary is zero
                    Dc[c_index(x,y,Nc-1,ion,Nx)*2+1] = scale*D[ion]*(alNcL)/pow(tortuosity,2);
                }
                else {
                    //diffusion coefficients in the extracellular space proportional to volume fraction
                    Dc[c_index(x,y,Nc-1,ion,Nx)*2+1] = scale*D[ion]/2*(alNcL+alNcU)/pow(tortuosity,2);

                }
                //diffusion coefficients in neuronal compartment equal to 0
                Dc[c_index(x,y,0,ion,Nx)*2+1] = 0.0;
                //diffusion coefficients in glial compartment proportional to default volume fraction
//			    Dc[c_index(x,y,1,ion,Nx)*2+1] = 0.25*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2); //0.25
                Dc[c_index(x,y,1,ion,Nx)*2+1] = 0.001*scale*D[ion]*alphao[al_index(0,0,1,Nx)]/pow(tortuosity,2);

            }
        }
    }
}


