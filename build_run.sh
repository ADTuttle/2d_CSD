#!/bin/bash
rm csd
make csd
#make debug
#./csd -malloc_log -malloc_debug 
#./csd
./csd -log_view
# ./csd -ksp_type gmres -ksp_gmres_restart 40  -pc_type ilu -pc_factor_levels 1 -pc_factor_fill 3 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type natural
# ./csd -ksp_view -ksp_type dgmres -ksp_gmres_restart 40 -ksp_dgmres_eigen 10 -ksp_dgmres_max_eigen 100 -ksp_dgmres_force -pc_type ilu -pc_factor_levels 1 -pc_factor_fill 3 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type natural
# ./csd -ksp_type fgmres -ksp_gmres_restart 40 
# ./csd -snes_type newtonls -snes_linesearch_type bt -snes_linesearch_order 2 -snes_linesearch_norms -snes_linesearch_alpha .1
# ./csd -snes_type newtontr -snes_trtol 1e-5 -snes_tr_eta 2 -snes_tr_delta0 1e-6
#./csd -ksp_monitor 
#lldb ./csd

