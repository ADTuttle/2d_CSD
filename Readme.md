## Future Fixes
Most parameters (most notably those stored in con_vars) are currently constant in space. However, the indexing on these still uses the indexing functions, just with x and y set to 0. If these were modified in the future, that needs to be updated.

However, some calls to con_vars->ao or zo, might have the wrong indexing function. It needs to be indexed by phi_index and not al_index.

A consequence of the transition is that c_index(x,y,comp,ion) does not agree with Ind_1(x,y,ion,comp). Ind_1 is used for the row and column ordering of the matrix. It has not been swapped because Ind_1 counts up the "ion" number by Nv and not Ni. This allows us to put in arguments like  Ind_1(x,y,Ni/Ni+1,comp) to represent the voltage and water flow equations.