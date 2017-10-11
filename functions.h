#include "constants.h"
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

//Function to initialize the state
void init(struct SimState*);

//Set the paramaters based on the constants
void set_params(struct SimState*,struct ConstVars*);

//Linear current-voltage flux relation
void mclin(struct FluxPoint*,double,int,double,double,double,int);



#endif