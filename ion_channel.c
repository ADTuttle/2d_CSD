#include "constants.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void mclin(struct FluxData *flux,PetscInt index,double pc,PetscInt zi,double ci,double ce,double phim,PetscInt ADD)
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
void mcGoldman(struct FluxData *flux,PetscInt index,double pc,PetscInt zi,double ci,double ce,double phim,PetscInt ADD)
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
    double dfdphim = (double)zi/2*(mflux*s*w+mi-me);
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
double xoverexpminusone(double v,double aa,double bb,double cc,PetscInt dd)
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
double cz(const double *cmat,const PetscInt *zvals,PetscInt x,PetscInt y,PetscInt comp)
{
	//function to compute sum over i of c_i*z_i
	double accumulate=0;
	for(PetscInt ion=0;ion<Ni;ion++)
	{
		accumulate += zvals[ion]*cmat[c_index(x,y,comp,ion)];
	}
	return accumulate;
}
void diff_coef(double *Dc,const double *alp,double scale)
{
  //diffusion coefficients at all points, for all ions, in all compartments, in both x and y directions
	double tortuosity=1.6;
	double alNcL,alNcR,alNcU;
  	for(PetscInt x=0;x<Nx;x++)
  	{
	  	for(PetscInt y=0;y<Ny;y++)
	  	{

	  		alNcL=1-alp[al_index(x,y,0)]-alp[al_index(x,y,1)]; //Left extracell
			alNcR = 0;
			if(x<Nx-1)
			{
				alNcR = 1 - alp[al_index(x + 1, y, 0)] - alp[al_index(x + 1, y, 1)]; //Right extracell
			}
            alNcU = 0;
			if(y<Ny-1)
			{
				alNcU = 1 - alp[al_index(x, y + 1, 0)] - alp[al_index(x, y + 1, 1)];
			}
		  	for(PetscInt ion = 0; ion<Ni;ion++)
		  	{
			    //diffusion coefficients in x direction
			    if(x==(Nx-1))
			    {
			    	//Boundary is zero
			    	Dc[c_index(x,y,Nc-1,ion)*2] = 0;
			    }
			    else
			    {
					//diffusion coefficients in the extracellular space proportional to volume fraction
			    	Dc[c_index(x,y,Nc-1,ion)*2] = scale*D[ion]/2*(alNcL+alNcR)/(pow(tortuosity,2));
			    }
			    //diffusion coefficients in neuronal compartment equal to 0
			    Dc[c_index(x,y,0,ion)*2] = 0.0 ;
			    //diffusion coefficients in glial compartment proportional to default volume fraction
			    Dc[c_index(x,y,1,ion)*2] = 0*scale*D[ion]*alphao[al_index(0,0,1)]/pow(tortuosity,2);
			    //diffusion coefficients in y direction
			    if(y==(Ny-1))
			    {
			    	//Boundary is zero
			    	Dc[c_index(x,y,Nc-1,ion)*2+1] = 0;
			    }
			    else
			    {
			    	//diffusion coefficients in the extracellular space proportional to volume fraction
			    	Dc[c_index(x,y,Nc-1,ion)*2+1] = scale*D[ion]/2*(alNcL+alNcU)/pow(tortuosity,2);

				}
				//diffusion coefficients in neuronal compartment equal to 0
			    Dc[c_index(x,y,0,ion)*2+1] = 0.0;
			    //diffusion coefficients in glial compartment proportional to default volume fraction
			    Dc[c_index(x,y,1,ion)*2+1] = 0.25*scale*D[ion]*alphao[al_index(0,0,1)]/pow(tortuosity,2);

		  	}
		}
	}
}

void gatevars_update(struct GateType *gate_vars,struct SimState *state_vars,double dtms,PetscInt firstpass)
{
	if(firstpass)
	{
		//membrane potential in mV
		double v = (state_vars->phi[phi_index(0,0,0)]-state_vars->phi[phi_index(0,0,Nc-1)])*RTFC;
		//Iniitialize the poPetscInt gating variables
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
  		for(PetscInt i=0;i<Nx*Ny;i++)
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
    } else //if it's not the firstpass, then we actually have values in v.
	{
		double v, alpha,beta;
		for(PetscInt x=0;x<Nx;x++)
		{
			for(PetscInt y=0;y<Ny;y++)
			{
				//membrane potential in mV
				v = (state_vars->phi[phi_index(x,y,0)]-state_vars->phi[phi_index(x,y,Nc-1)])*RTFC;

				//compute current NaT
		  		//gating variables mNaT
		  		alpha = xoverexpminusone(v,0.32,51.9,0.25,1); //0.32*(Vm+51.9)./(1-exp(-0.25*(Vm+51.9)))
		  		beta = xoverexpminusone(v,0.28,24.89,0.2,0); //0.28*(Vm+24.89)./(exp(0.2*(Vm+24.89))-1)
		    	gate_vars->mNaT[xy_index(x,y)] = (gate_vars->mNaT[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		//gating variable hNaT
		  		alpha = 0.128*exp(-(0.056*v+2.94));
		  		beta = 4/(exp(-(0.2*v+6))+1);
		    	gate_vars->hNaT[xy_index(x,y)] = alpha/(alpha+beta);
		    	gate_vars->hNaT[xy_index(x,y)] = (gate_vars->hNaT[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		    	gate_vars->gNaT[xy_index(x,y)] = pow(gate_vars->mNaT[xy_index(x,y)],3)*gate_vars->hNaT[xy_index(x,y)];
		  		//compute current NaP
		  		//gating variable mNaP
		  		alpha = 1/(1+exp(-(0.143*v+5.67)))/6;
		  		beta = 1.0/6.0-alpha; //1./(1+exp(0.143*Vm+5.67))/6
		  		gate_vars->mNaP[xy_index(x,y)] = (gate_vars->mNaP[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		//gating variable hNaP
		  		alpha = 5.12e-6*exp(-(0.056*v+2.94));
		  		beta = 1.6e-4/(1+exp(-(0.2*v+8)));
		  		gate_vars->hNaP[xy_index(x,y)] = (gate_vars->hNaP[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gNaP[xy_index(x,y)] = pow(gate_vars->mNaP[xy_index(x,y)],2)*gate_vars->hNaP[xy_index(x,y)];

		  		//compute KDR current
		  		//gating variable mKDR
		  		alpha = xoverexpminusone(v,0.016,34.9,0.2,1); //0.016*(Vm+34.9)./(1-exp(-0.2*(Vm+34.9)))
		  		beta = 0.25*exp(-(0.025*v+1.25));
		    	gate_vars->mKDR[xy_index(x,y)] = (gate_vars->mKDR[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gKDR[xy_index(x,y)] = pow(gate_vars->mKDR[xy_index(x,y)],2);

		  		//compute KA current
		  		//gating variable mKA
		  		alpha = xoverexpminusone(v,0.02,56.9,0.1,1); //0.02*(Vm+56.9)./(1-exp(-0.1*(Vm+56.9)))
		  		beta = xoverexpminusone(v,0.0175,29.9,0.1,0); //0.0175*(Vm+29.9)./(exp(0.1*(Vm+29.9))-1)
		    	gate_vars->mKA[xy_index(x,y)] = (gate_vars->mKA[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		 		//gating variable hKA
		  		alpha = 0.016*exp(-(0.056*v+4.61));
		  		beta = 0.5/(exp(-(0.2*v+11.98))+1);
		    	gate_vars->hKA[xy_index(x,y)] = (gate_vars->hKA[xy_index(x,y)] + alpha*dtms)/(1+(alpha+beta)*dtms);

		  		gate_vars->gKA[xy_index(x,y)] = pow(gate_vars->mKA[xy_index(x,y)],2)*gate_vars->hKA[xy_index(x,y)];
  			}
  		}
  	}
}

void excitation(struct ExctType *exct,double t)
{
  //compute excitation conductance to trigger csd
  //Leak conductances in mS/cm^2
  //all units converted to mmol/cm^2/sec
  	double pexct,pany;
  	double xexct;
  	for(PetscInt i=0;i<Nx;i++)
  	{
	  	for(PetscInt j=0;j<Ny;j++)
	  		{
	  		if(two_points_exct)
	  		{
	  			fprintf(stderr, "Two Points Exct is not implemented\n");
	    		exit(EXIT_FAILURE); /* indicate failure.*/
	    		// centerpointx=floor(Nx/2)
	    		// xexct[centerpointx+1:centerpointx+Nexct,1:Nexct]=xexct[1:Nexct,1:Nexct]
	    		// xexct[centerpointx+1-Nexct:centerpointx,1:Nexct]=xexct[Nexct:-1:1,1:Nexct]
	  		}
	  		if( t<texct && i<Nexct && j<Nexct)
	  		{
	    		pexct=pmax*(sin(pi*t/texct))*RTFC/FC;
	    		xexct=pow((cos(pi/2*(i+.5)/Nexct))*(cos(pi/2*(j+.5)/Nexct)),2);
	    		pany=pexct*xexct;
	    		exct->pNa[xy_index(i,j)]=pany;
	    		exct->pK[xy_index(i,j)]=pany;
	    		exct->pCl[xy_index(i,j)]=pany;
			}
	  		else
	  		{
	    		//pexct=0*RTFC/FC
	    		exct->pNa[xy_index(i,j)]=0;
	    		exct->pK[xy_index(i,j)]=0;
	    		exct->pCl[xy_index(i,j)]=0;
			}
		}
	}

}

void ionmflux(struct FluxData *flux,struct SimState *state_vars,struct SimState *state_vars_past,struct GateType *gvars, struct ExctType *gexct,struct ConstVars *con_vars)
{
    //Variables to save to for ease of notation
    double vm,vmg,vmgp;
    double ci,cg,ce,cgp,cep;
    double Na,K;//Variables for pump (so it's clear)

    //For calculationg permeabilities
    double pGHK,pLin;
    double Ipump,NaKCl;
    for(PetscInt x=0;x<Nx;x++)
    {
        for(PetscInt y=0;y<Ny;y++)
        {
            vm = state_vars->phi[phi_index(x,y,0)]-state_vars->phi[phi_index(x,y,Nc-1)];
            vmg = state_vars->phi[phi_index(x,y,1)]-state_vars->phi[phi_index(x,y,Nc-1)];
            vmgp = state_vars_past->phi[phi_index(x,y,1)]-state_vars_past->phi[phi_index(x,y,Nc-1)];

            //Compute Na Channel currents
            ci = state_vars->c[c_index(x,y,0,0)];
            cg = state_vars->c[c_index(x,y,1,0)];
            ce = state_vars->c[c_index(x,y,Nc-1,0)];

            //Neurons
            pGHK = pNaT*gvars->gNaT[xy_index(x,y)]+pNaP*gvars->gNaP[xy_index(x,y)];
            pLin = con_vars->pNaLeak + gexct->pNa[xy_index(x,y)]; //Add excitation
            //Initialize GHK Flux
            mcGoldman(flux,c_index(x,y,0,0),pGHK,1,ci,ce,vm,0);
            //Add leak current to that.
            mclin(flux,c_index(x,y,0,0),pLin,1,ci,ce,vm,1);
            //Glial NaLeak
            mclin(flux,c_index(x,y,1,0),con_vars->pNaLeakg,1,cg,ce,vmg,0);

            // Compute K channel Currents
            ci = state_vars->c[c_index(x,y,0,1)];
            cg = state_vars->c[c_index(x,y,1,1)];
            ce = state_vars->c[c_index(x,y,Nc-1,1)];

            //Neurons
            pGHK = pKDR*gvars->gKDR[xy_index(x,y)]+pKA*gvars->gKA[xy_index(x,y)];
            pLin = pKLeak+gexct->pK[xy_index(x,y)]; //add excitation
            mcGoldman(flux,c_index(x,y,0,1),pGHK,1,ci,ce,vm,0);
            mclin(flux,c_index(x,y,0,1),pLin,1,ci,ce,vm,1);

            // Glial K Leak(using past)
            cgp = state_vars_past->c[c_index(x,y,1,1)];
            cep = state_vars_past->c[c_index(x,y,Nc-1,1)];

            pLin = pKIR*inwardrect(cgp,cep,vmgp)*pKLeakadjust;
            mclin(flux,c_index(x,y,1,1),pLin,1,cg,ce,vmg,0);

            //Compute Cl Channel Current
            ci = state_vars->c[c_index(x,y,0,2)];
            cg = state_vars->c[c_index(x,y,1,2)];
            ce = state_vars->c[c_index(x,y,Nc-1,2)];

            //Neurons
            pLin = pClLeak+gexct->pCl[xy_index(x,y)]; //add excitation
            mclin(flux,c_index(x,y,0,2),pLin,-1,ci,ce,vm,0);

            //Glia
            mclin(flux,c_index(x,y,1,2),pClLeakg,-1,cg,ce,vmg,0);

            //Pump Currents(all past values)

            //Neurons
            Na = state_vars_past->c[c_index(x,y,0,0)];
            K = state_vars_past->c[c_index(x,y,Nc-1,1)];

            Ipump = npump*con_vars->Imax/(pow(1+mK/K,2)*pow(1+mNa/Na,3));

            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x,y,0,0)]+=3*Ipump; //Na part
            flux->mflux[c_index(x,y,0,1)]-=2*Ipump; //K part

            //Glia
            Na = state_vars_past->c[c_index(x,y,1,0)];
            //K is the same(extracellular)
            Ipump = glpump*con_vars->Imaxg/(pow(1+mK/K,2)*pow(1+mNa/Na,3));
            //Add to flux(it's explicit so no derivatives)
            flux->mflux[c_index(x,y,1,0)]+=3*Ipump; //Na part
            flux->mflux[c_index(x,y,1,1)]-=2*Ipump; //K part

            //NaKCl Cotransporter
            //I'm going to reuse variables names...
            Na = state_vars_past->c[c_index(x,y,1,0)];//glia Na
            K = state_vars_past->c[c_index(x,y,1,1)]; // glia K.
            cgp = state_vars_past->c[c_index(x,y,1,2)]; //glia Cl

            cep = state_vars_past->c[c_index(x,y,Nc-1,0)];//Ext Na
            ce = state_vars_past->c[c_index(x,y,Nc-1,1)]; // Ext K.
            ci = state_vars_past->c[c_index(x,y,Nc-1,2)]; //Ext Cl

            NaKCl = con_vars->pNaKCl*log(Na*K*pow(cgp,2)/(cep*ce*pow(ci,2)));
            //Add to flux
            flux->mflux[c_index(x,y,1,0)]+=NaKCl; //Sodium
            flux->mflux[c_index(x,y,1,1)]+=NaKCl; //K
            flux->mflux[c_index(x,y,1,2)]+=2*NaKCl; //Cl

            //Change units of flux from mmol/cm^2 to mmol/cm^3/s
            for(PetscInt ion=0;ion<Ni;ion++)
            {
                flux->mflux[c_index(x,y,Nc-1,ion)] = 0;
                for(PetscInt comp = 0;comp<Nc-1;comp++)
                {
                    flux->mflux[c_index(x,y,comp,ion)]=flux->mflux[c_index(x,y,comp,ion)]/ell;
                    flux->dfdci[c_index(x,y,comp,ion)]=flux->dfdci[c_index(x,y,comp,ion)]/ell;
                    flux->dfdce[c_index(x,y,comp,ion)]=flux->dfdce[c_index(x,y,comp,ion)]/ell;
                    flux->dfdphim[c_index(x,y,comp,ion)]=flux->dfdphim[c_index(x,y,comp,ion)]/ell;

                    //And calculate the extracellular flux
                    flux->mflux[c_index(x,y,Nc-1,ion)] -= flux->mflux[c_index(x,y,comp,ion)];
                }

            }
        }
    }
}
void wflowm(struct FluxData *flux,struct SimState *state_vars,struct ConstVars *con_vars)
{
  //piw = sum of c over ions + ao/alpha
  // piw is the total number of ions in a compartment
  //outward transmembrane water flow seen as a function of
  //osmotic pressure and volume fraction or pressure.
    double dwdpi,dwdal,piw,piwNc;
    for(PetscInt x=0;x<Nx;x++)
    {
        for(PetscInt y=0;y<Ny;y++)
        {
            //Calculate the pi for extracellular
            piwNc = 0;
            for(PetscInt ion=0;ion<Ni;ion++)
            {
                piwNc +=state_vars->c[c_index(x,y,Nc-1,ion)];
            }
            piwNc +=con_vars->ao[al_index(0,0,Nc-1)]/(1-state_vars->alpha[al_index(x,y,0)]-state_vars->alpha[al_index(x,y,1)]);
            for(PetscInt comp = 0;comp<Nc-1;comp++)
            {
                piw = 0;
                for(PetscInt ion=0;ion<Ni;ion++)
                {
                    piw +=state_vars->c[c_index(x,y,comp,ion)];
                }
                piw +=con_vars->ao[al_index(0,0,comp)]/state_vars->alpha[al_index(x,y,comp)];
                //ao, zeta1, and zetalpha are currently constant in space
                dwdpi = con_vars->zeta1[al_index(0,0,comp)];
                dwdal = con_vars->zeta1[al_index(0,0,comp)]*con_vars->zetaalpha[al_index(0,0,comp)];

                flux->wflow[al_index(x,y,comp)] = dwdpi*(piwNc-piw)+dwdal*(state_vars->alpha[al_index(x,y,comp)]-alphao[comp]);
                flux->dwdpi[al_index(x,y,comp)] = dwdpi;
                flux->dwdal[al_index(x,y,comp)] = dwdal;
            }
        }
    }
}


