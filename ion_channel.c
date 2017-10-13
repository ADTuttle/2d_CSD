#include "constants.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void mclin(struct FluxPoint *flux,int index,double pc,int zi,double ci,double ce,double phim,int ADD)
{
	//Returns the flux value by ref.
	// pc is the permeativity, zi is valence, ci/e intra/extra concentration
	//phim is membrane voltage, index is the index in the flux struct
	//compute value and derivatives of the function:
    //mflux=pc.*(log(ci./ce)+z*phim)
    //for trans-membrane flux of an ion obeying a linear current-voltage equation
    if(ADD) //If add, we accumulate the result
    {
    	flux->mflux[index] += pc*(log(ci/ce)+zi*phim);
	    flux->dfdci[index] += pc/ci;
		flux->dfdce[index] += -pc/ce;
	  	flux->dfdphim[index] += zi*pc;
    }
    else //If not add we reninitialize
    {
		flux->mflux[index] = pc*(log(ci/ce)+zi*phim);
	    flux->dfdci[index] = pc/ci;
		flux->dfdce[index] = -pc/ce;
	  	flux->dfdphim[index] = zi*pc;
  	}

  	return;
}
void mcGoldman(struct FluxPoint *flux,int index,double pc,int zi,double ci,double ce,double phim,int ADD)
{
    //compute value and derivatives of the function:
    //mflux=p.*(z*phim).*(ci.*exp(z*phim)-ce)./(exp(z*phim)-1)
    //for trans-membrane flux of an ion obeying the GHK equations
    double xi = zi*phim/2;
    double r = exp(xi);

    //compute s=x/sinh(x)
    //Watch out for division by zero
    double s;
    if(xi!=0)
    {
    	s = xi/sinh(xi);
    }
    else
    {
    	s = 1;
    }
    //compute dfdci,dfdce,mflux
    double dfdci = pc*s*r;
    double dfdce = -pc*s/r;
    double mi = ci*dfdci;
    double me = ce*dfdce;
    double mflux = mi+me;
    //compute w=(sinh(x)/x-cosh(x))/x
    //p0 = abs.(xi).<0.2 #treat differently if small
    //p1 = !p0
    //x0 = xi[p0] #use Taylor expansion for small values
    //x1 = xi[p1] #use formula for larger values
    double w;
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
    double dfdphim = zi/2*(mflux*s*w+mi-me);

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
double xoverexpminusone(double v,double aa,double bb,double cc,int dd)
{
	//computes aa*(v+bb)/(exp(cc*(v+bb))-1) if dd==0
 	//computes aa*(v+bb)/(1-exp(-cc*(v+bb)) otherwise
  	//for computing gating variables
  	v+=bb;
  	if(v==0)
  	{
  		return aa/cc;
  	}
  	if(dd==0)
  	{
  		return aa*v/(2*sinh(cc/2*v))*exp(-cc/2*v);
  	}
  	else
  	{
  		return aa*v/(2*sinh(cc/2*v))*exp(cc/2*v);
  	}
}
double inwardrect(double ci,double ce,double phim)
{
	//inwardrect determines the effective conductance for the inward-rectifying
  	//potassium channel - formula and constants from Steinberg et al 2005
  	double Enernst = RTFC*log(ce/ci);
  	double EKdef = -85.2; //#-85.2 mV is base reversal potential
  	double cKo = .003; //#3 mM is base extracellular potassium concentration
  	return sqrt(ce/cKo)*(1+exp(18.5/42.5))/(1+exp((RTFC*phim-Enernst+18.5)/42.5))*(1+exp((-118.6+EKdef)/44.1))/(1+exp((-118.6+RTFC*phim)/44.1));
}
double cz(double *cmat,const int *zvals,int x,int y,int comp) 
{ 
	//function to compute sum over i of c_i*z_i
	double accumulate=0;
	for(int ion=0;ion<Ni;ion++)
	{
		accumulate += z[ion]*cmat[c_index(x,y,comp,ion)];
	}
	return accumulate;
}


void gatevars_update(struct GateType *gate_vars,struct SimState *state_vars,int firstpass)
{
	if(firstpass)
	{
		//membrane potential in mV
		double v = state_vars->phi[phi_index(0,0,0)]-state_vars->phi[phi_index(0,0,Nc-1)];
		//Iniitialize the point gating variables
		//Cause we assume a uniform start
		double alpha,beta;

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
  		for(int i=0;i<Nx*Ny;i++)
  		{
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
    }
    else //if it's not the firstpass, then we actually have values in v.
	{
		double v, alpha,beta;
		for(int x=0;x<Nx;x++)
		{
			for(int y=0;y<Ny;y++)
			{
				//membrane potential in mV
				v = state_vars->phi[phi_index(x,y,0)]-state_vars->phi[phi_index(x,y,Nc-1)];
				
				//compute current NaT
		  		//gating variables mNaT
		  		alpha = xoverexpminusone(v,0.32,51.9,0.25,1); //0.32*(Vm+51.9)./(1-exp(-0.25*(Vm+51.9)))
		  		beta = xoverexpminusone(v,0.28,24.89,0.2,0); //0.28*(Vm+24.89)./(exp(0.2*(Vm+24.89))-1)
		    	gate_vars->mNaT[xy_index(x,y)] = alpha/(alpha+beta);


		  		//gating variable hNaT
		  		alpha = 0.128*exp(-(0.056*v+2.94));
		  		beta = 4/(exp(-(0.2*v+6))+1);
		    	gate_vars->hNaT[xy_index(x,y)] = alpha/(alpha+beta);

		    	gate_vars->gNaT[xy_index(x,y)] = pow(gate_vars->mNaT[xy_index(x,y)],3)*gate_vars->hNaT[xy_index(x,y)];

		  		//compute current NaP
		  		//gating variable mNaP
		  		alpha = 1/(1+exp(-(0.143*v+5.67)))/6;
		  		beta = 1/6-alpha; //1./(1+exp(0.143*Vm+5.67))/6
		    	gate_vars->mNaP[xy_index(x,y)] = alpha/(alpha+beta);

		  		//gating variable hNaP
		  		alpha = 5.12e-6*exp(-(0.056*v+2.94));
		  		beta = 1.6e-4/(1+exp(-(0.2*v+8)));
		    	gate_vars->hNaP[xy_index(x,y)] = alpha/(alpha+beta);
		  	
		  		gate_vars->gNaP[xy_index(x,y)] = pow(gate_vars->mNaP[xy_index(x,y)],2)*gate_vars->hNaP[xy_index(x,y)];

		  		//compute KDR current
		  		//gating variable mKDR
		  		alpha = xoverexpminusone(v,0.016,34.9,0.2,1); //0.016*(Vm+34.9)./(1-exp(-0.2*(Vm+34.9)))
		  		beta = 0.25*exp(-(0.025*v+1.25));
		    	gate_vars->mKDR[xy_index(x,y)] = alpha/(alpha+beta);

		  		gate_vars->gKDR[xy_index(x,y)] = pow(gate_vars->mKDR[xy_index(x,y)],2);

		  		//compute KA current
		  		//gating variable mKA
		  		alpha = xoverexpminusone(v,0.02,56.9,0.1,1); //0.02*(Vm+56.9)./(1-exp(-0.1*(Vm+56.9)))
		  		beta = xoverexpminusone(v,0.0175,29.9,0.1,0); //0.0175*(Vm+29.9)./(exp(0.1*(Vm+29.9))-1)
		    	gate_vars->mKA[xy_index(x,y)] = alpha/(alpha+beta);

		 		//gating variable hKA
		  		alpha = 0.016*exp(-(0.056*v+4.61));
		  		beta = 0.5/(exp(-(0.2*v+11.98))+1);
		    	gate_vars->hKA[xy_index(x,y)] = alpha/(alpha+beta);

		  		gate_vars->gKA[xy_index(x,y)] = pow(gate_vars->mKA[xy_index(x,y)],2)*gate_vars->hKA[xy_index(x,y)];
  			}
  		}
  	}
}
