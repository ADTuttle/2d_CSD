# 2D Cortical Spreading Depression Simulation

# Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [File Details](#files)
4. [Data Structures and Indexing](#data)
5. [Biology](#biology)
6. [Algorithm and Petsc](#algorithm)
7. [Future Fixes](#future)

<a name="introduction"></a>
## Introduction 

This code runs a simulation of cortical spreading depression. The model we use is a 2D electrodiffusion model. We simulate 4 main variables: ionic concentrations, voltages, cellular volume fractions, and gating variables. The model breaks up different cell types into "compartments" and models them differently using physiological parameters as well as ion channels. By default we have 3 compartments: neuronal, glial, and extracellular. We also model different ions separately. Currently we model: Sodium, Potassium, and Chloride ions. It is possible to add more with somewhat ease. But that would require an understanding of how the biology is implemented in the code.

All of the numerical analysis heavy-lifting is done using Petsc. So you must have an installation of it to run any variation of this code.

<a name="installation"></a>
# Installation 
Should be simple.
1. Clone/download this repo.
2. Install Petsc.
3. Modify the Makefile. Changing the Petsc_dir path at the top will do.
4. Run "make csd" or "make debug" to build the "csd" executable. 
5. ./csd and it's good to go. The run can be modified with any Petsc options menu things and will be taken into account. 

<a name="files"></a>
## File Details 

#### C-Files
* functions.h: Catch all header for every function declared in *.c files. 
* constants.h: Define all constants/parameters used throughout the code. Also includes data structure definitions. Change one constant here and it will effect everything.
* constants.c: Modify parameters to have spatial dependence. Also contains indexing routines. And allocates array memory.
* csd_main.c: Runs the main program. General order of initialize -> update solution -> record -> loop -> end.
* initialize.c: Initialize everything from structures, petsc, variables, parameters, initial steady state. 
* ion_channel.c: Functions to calculate all fluxes. Ion fluxes from ion channels and water fluxes. Plus all of their support functions: Calculate exponential quantities and permeabilities. Also contains update to gating variables (solved as individual ODES using backward euler).
* misc_print_plot.c: functions to print all state/flux variables whenever (unused). And the print to file function. 
* array_functions.c: Calculate maximum values, maximum difference, or l2 norm of arrays.
* update_solution.c: Jacobians and Residuals calculated here. Also my own Newton Solve is implemented here. 
* grid_update.c: Predictor function. Solves a smaller 3x3 grid problem as a predictor for update_solution.c
* linear_update.c: Jacoabians and Residuals for a linearized version of the discretization.

#### Others:
* makefile: makefile for code.
* CMakeLists.txt: Cmake file for code (I prefer the makefile).
* csd: Executable generated by the makefile (Cmake makes 2D_CSD).
* build_run.sh: Bash file to build and run code easily. Includes some petsc options.

#### Generated Files:
* data_csd.txt: Save file for each variable in space and time. Frequency of savinings defined by trecordstep in constants.h
* csd_dt.txt: Predictor adaptive time steping, time steps.
* flux_csd.txt: Total flux through a small circle in center of domain. Broken up by ionic diffusion and electro-diffusion parts.
* save_csd.txt: Saved state at end of file. Can be used to restart from a new run from the end of an old one.
* grad_field.txt: electro-diffusive gradient field of each ion and compartment.
* measures.txt: Any time varying measurement. For now this is whole grid potassium amount in glia. And percentage of each ion in glia.
* timing.txt: Saving Nx, Ny, dt, and total time a run took in a running tally file. Never get's overwritten.
* log.txt: Petsc Logging details if Profiling is set to on.

<a name="data"></a>
## Data Structures and Indexing
Before getting into this. This code was constructed initially from having c, phi, and alpha coded as arrays. So that legacy is leftover with how everything is calculated even though Petsc has better ways of doing this with vectors.

#### Simstate
Contains c,phi,alpha. These are put together into the Vector v. We use the function extract_subarray to do two things.

First: Extract a pointer to the c_vec, phi_vec, and al_vec portion of it (using what petsc calls Index sets: c_ind,phi_ind, and al_ind respectively). These are still Petsc Vectors.

Second: Get arrays from these subvectors. We use these to easily access data.

**IMPORTANT USAGE** 

The vector v is THE storage object. Everything else is a pointer. If you want to manually update the concentration you must call extract_subarray to get access to the data. Then restore_array to pass that data back to the vector.

If you want petsc to update the vector v (e.g what happens during a SNESSolve), you must restore_array before this is called.

#### FluxData
Contains all flux data. Even though the extracellular space is empty for most of these, we include that storage space for consistency sake.

* mflux: membrane flux value. Value for each coordinate, ion, and compartment.
* dfdci: derivative with respect to inside concentration. Value for each coordinate, ion, and either neuronal or glial compartment.
* dfdce: derivative with respect to outside concentration. Value for each coordinate, ion, and either neuronal or glial compartment.
*dfdphim: Derivative with respect to voltage. Value for each coordinate, ion, and either neuronal or glial compartment. Derivative of intracellular and extracellular are the same bit dfdphi_inside = -dfdphi_outside.
* wflow: Osmotic water flux. Value for each coordinate and intracellular compartment.
* dwdpi: Derivative of the osmotic pressure pi. (This is the term that appears in water flow derivative for concentrations, each ion has same constant)
*dwdal: Derivative of the osmotic pressure with respect to alpha.

#### GateType
Storage of each gating variable type at each coordinate.
* m: Activation for each channel.
* h: Inactivation for each channel.
* g: Overall gating variable openness.

List of ion channels:
* Neuronal transient sodium
* Neuronal perssitent sodium
* Neuronal potassium delayed rectifier
* Neuronal active potassium.

#### ExctType
The excitation type. To initiate CSD we introduce a flux of all 3 ions into the system. These are stored in this type. To change the excitation alter the "void excitation(struct ExctType *exct,PetscReal t)" function in ion_channels.c.

#### ConstVars
Container for various constants that are calculated to set the system in steady state. These get set in set_params. Contains:

* pNaKCl: Sodium, potassium, chloride cotransporter in glia.
* Imax/Imaxg: Sodium Potassium ATPase current modifier (neuron/glia respectively).
* pNaLeak/pNaLeakg: Sodium leak current (neuron/glia respectively).
* ao: Organic ion amount.
* zo: Average organic ion valence.
* kappa: water flow (since no fluid set to 0)
* zeta1: hydraulic permeability
* S: Indicates whether zetaalpha is the stiffness (true) or 1/stiffness (false)
* zetaalpha: Stiffness.

#### Solver:
Container for all Solver data.
* Q: Update vector (Is Ainv*Res)
* Res: Residual vector
* A: Jacobian matrix
* snes: Nonlinear solver context
* ksp: krylov solver context
* pc: Preconditioner context
* size: number of processors. 

#### AppCtx
User context that Petsc allows you to pass. It's a container for everything that isn't the function vector or matrix.
* state_vars: Current state
* state_vars_past: Previous state.
* slvr: Solver context
* flux: flux data
* gate_vars: gating variables
* gexct: excitation
* con_vars: constant variables.
* Dcs: Ionic diffusion array.
* Dcb: Ionic diffusion to bath.
* dt: Time step

### Indexing
Throughout the code we make use of index functions. The 6 currently in use are:
* c_index(x,y,comp,ion): tracks concentration variable/equation quantities. Also anything that varies by ion and compartment. E.g ion channel fluxes.
* phi_index(x,y,comp): Tracks voltage variable/equation quantities. Also anything that varies in all 3 compartments.
* al_index(x,y,comp): Tracks volume fraction variable/equation quantities. Also anything that varies in only intracellular spaces. 
* xy_index(x,y): Tracks xy coordinate indices when needed.
* Ind_1(x,y,"ion",comp): Tracks all variable quantities that have equations represented in the residual/jacobian. "ion" tracks ions if we are in the concentration equations. E.g 1,2,...,Ni-1 will be concentration equations. "ion"=Ni corresponds to voltage equations. And "ion"=Ni+1 corresponds to volume equations.
* Ind_nx(x,y,ion,comp,nx): Same as Ind_1 but allows variable number of x-y gridpoints. Needed in order to do multigrid.

To spell out more. For 3 ions and 3 compartments. There are 14 variables (Nv=14). These are broken into 3 groups.

Concentration equations:
Ind_1(x,y,0,0), Ind_1(x,y,0,1), Ind_1(x,y,0,2), Ind_1(x,y,1,0),...,Ind_1(x,y,2,0),..., Ind_1(x,y,2,2).

Voltage Equations:
Ind_1(x,y,Ni,0), Ind_1(x,y,Ni,1), Ind_1(x,y,Ni,2)

Volume Equations:
Ind_1(x,y,Ni+1,0), Ind_1(x,y,Ni+1,1)

<a name="biology"></a>
## Biology 
There are 3 parts where the underlying biology would be modified. First, the initial values of c, phi, and alpha are set in the function init in initialize.c.

Second, parameter values can be changed in constants.h (which has a huge list of constants). While others could be changed in the function set_params in initialize.c. set_params has values chosen so that given values from constants.h will start us at steady state. If new ion channels are added some of these formula will need to be checked.

Third, ion_channels.c. If a new ion channel was added you would have to modify it in this file most of all. The function ionmflux calculates the flux and derivatives from the gating variables. It uses two functions mclin (linear current voltage relationship) and mcgoldman (using GHK type relation), the exact functions for each are detailed in the comments in each function. Both these functions have the same format:
void mclin(struct FluxData *flux,PetscInt index,PetscReal pc,PetscInt zi,PetscReal ci,PetscReal ce,PetscReal phim,PetscInt ADD)

* flux: the whole flux structure being updated.
* index: index corresponding the ion we are updating. For example updating the neuronal sodium will at int x and int y will have index = c_index(x,y,0,0).
* pc: permeability
* zi: valence of ion
* ci: intracellular concentration of this ion.
* ce: extracellular concentration of this ion.
* phim: Membrane voltage difference.
* ADD: We update by passing flux by reference. ADD = 0 means we reset the current value. ADD = 1 means we add to the current value.

A new voltage gated ion channel being added would also need to modify gatevars_update. Time in this function is in ms (not seconds like the rest). Also, it has two parts a firstpass uniform initialization of values and an afterwards update using backwards euler (explicitly coded in).

Lastly, for different interesting excitations modify the values in excitation (can easily make them vary over the domain). Ionic diffusion is modified in diff_coef. And water flux in w_flux.

<a name="algorithm"></a>
## Algorithm and Petsc 
#### Algorithm
The algorithm we use can be summed up as:
1. Create all data structures (no new memory is allocated after)
2. Set initial values.
3. Iterate with a larger time step until steady state is reached.
4. Begin solve loop.
	a. Update c, phi, alpha (possibly alpha separately).
	b. Update gating variables.
	c. Record values every krecordfreq.

The update will always involve solving a nonlinear system via newton's method. Which in term has a matrix being "inverted". The details in this is determined by the petsc settings used. Additionally, my own newton_solve code is implemented and could be used if desired.

#### Petsc
All the Petsc settings are set early in the program by calling inialize_petsc. In here we:
1. Create vectors.
2. Create Matrix using Get_Nonzero_in_Rows to count non-zeros on each row.
3. Initialize the matrix values to zero (unnecessary, but allows us to set MAT_NEW_NONZERO_LOCATION_ERR for saftey).
4. Create SNES (nonlinear solver) and set the functions that calculate Jacobian and Function (aka residual for us) value.
5. Set the SNES solver.
6. Set the KSP (krylov subspace) solver.
7. Choose a preconditioner.
8. Set all of the above from options (allows the solver details to be modified from the command line)


Currently the matrix is set to be a sequential matrix. An improvement would be a block matrix, which should just require modification in this part of the code. The only other use of Petsc in the code is the SETVALUE calls in calc_residual and calc_jacobian and the SNESSolve command in csd_main.

<a name="future"></a>
## Future Fixes 
1. Most parameters (most notably those stored in con_vars) are currently constant in space. However, the indexing on these still uses the indexing functions, just with x and y set to 0. If these were modified in the future, that needs to be updated.

However, some calls to con_vars->ao or zo, might have the wrong indexing function. It needs to be indexed by phi_index and not al_index.

2. A consequence of the copying from Julia is that c_index(x,y,comp,ion) does not agree with Ind_1(x,y,ion,comp). Ind_1 is used for the row and column ordering of the matrix. It has not been swapped because Ind_1 counts up the "ion" number by Nv and not Ni. This allows us to put in arguments like  Ind_1(x,y,Ni/Ni+1,comp) to represent the voltage and water flow equations.

3. Vectorize sections. Everything relies on taking values from c/phi/alpha arrays. But now we have these as vectors. We could use Petsc vector operations to likely do these faster. This improvement is necessary in order to realistically run this with MPI enabled.

4. Use Petsc's Data Managment (DM) types. Specifically the DMDA for structured meshes. This will allow automatic use of multigrid and FAS.