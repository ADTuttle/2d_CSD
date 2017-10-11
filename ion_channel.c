#include "constants.h"
#include "functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void mclin(struct FluxPoint *flux,double pc,int zi,double ci,double ce,double phim,int index)
{
	//Returns the flux value by ref.
	// pc is the permeativity, zi is valence, ci/e intra/extra concentration
	//phim is membrane voltage, index is the index in the flux struct
	//compute value and derivatives of the function:
    //mflux=pc.*(log(ci./ce)+z*phim)
    //for trans-membrane flux of an ion obeying a linear current-voltage equation
	flux->mflux[index] = pc*(log(ci/ce)+zi*phim);
    flux->dfdci[index] = pc/ci;
	flux->dfdce[index] = -pc/ce;
  	flux->dfdphim[index] = zi*pc;

  	return;
}
