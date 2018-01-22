# Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Future Fixes](#future)

## Introduction <a name="introduction"></a>

This code runs a simulation of cortical spreading depression. The model we use is a 2D electrodiffusion model. We simulate 4 main variables: ionic concentrations, voltages, cellular volume fractions, and gating variables. The model breaks up different cell types into "compartments" and models them differently using physiological parameters as well as ion channels. By default we have 3 compartments: neuronal, glial, and extracellular. We also model different ions separately. Currently we model: Sodium, Potassium, and Chloride ions. It is possible to add more with somewhat ease. But that would require an understanding of how the biology is implemented in the code.

All of the numerical analysis heavy-lifting is done using Petsc. So you must have an installation of it to run any variation of this code.

# Installation <a name="installation"></a>
Should be simple.
1. Clone/download this repo.
2. Install Petsc.
3. Modify the Makefile. Changing the Petsc_dir path at the top will do.
4. Run "make csd" or "make debug" to build the "csd" executable. 
5. ./csd and it's good to go. The run can be modified with any Petsc options menu things and will be taken into account. 


## Future Fixes <a name="future"></a>
Most parameters (most notably those stored in con_vars) are currently constant in space. However, the indexing on these still uses the indexing functions, just with x and y set to 0. If these were modified in the future, that needs to be updated.

However, some calls to con_vars->ao or zo, might have the wrong indexing function. It needs to be indexed by phi_index and not al_index.

A consequence of the transition is that c_index(x,y,comp,ion) does not agree with Ind_1(x,y,ion,comp). Ind_1 is used for the row and column ordering of the matrix. It has not been swapped because Ind_1 counts up the "ion" number by Nv and not Ni. This allows us to put in arguments like  Ind_1(x,y,Ni/Ni+1,comp) to represent the voltage and water flow equations.