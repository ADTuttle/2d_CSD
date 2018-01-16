#!/bin/bash
rm csd
make csd
#make debug
#./csd -malloc_log -malloc_debug 
# ./csd
# ./csd -ksp_type dgmres -ksp_gmres_restart 40 -ksp_dgmres_eigen 10 -ksp_dgmres_max_eigen 100 -ksp_dgmres_force -pc_type ilu -pc_factor_levels 1 -pc_factor_fill 3 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type natural
./csd -ksp_type fgmres -ksp_gmres_restart 40 

#./csd -ksp_monitor 
#lldb ./csd

