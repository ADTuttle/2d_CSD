#!/bin/bash

### ILU Preconditioning

# ./csd -ksp_type lgmres -ksp_gmres_restart 20 -ksp_lgmres_augment 10 -ksp_lgmres_constant -ksp_gmres_preallocate -pc_type ilu -pc_factor_levels 1 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type rcm

./csd -ksp_type bcgs -pc_type ilu -pc_factor_levels 1 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type rcm

# ./csd -ksp_type dgmres -ksp_gmres_restart 20 -ksp_dgmres_eigen 100 -ksp_dgmres_max_eigen 1000 -pc_type ilu -pc_factor_levels 1 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type rcm

# ./csd -ksp_type fbcgs -pc_type ilu -pc_factor_levels 1 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type rcm

# ./csd -ksp_type fbcgsr -pc_type ilu -pc_factor_levels 1 -pc_factor_reuse_ordering -pc_factor_mat_ordering_type rcm


#### ASM Preconditioning

# ./csd -ksp_type lgmres -ksp_gmres_restart 20 -ksp_lgmres_augment 10 -ksp_lgmres_constant -ksp_gmres_preallocate -pc_type asm -pc_asm_type none -pc_asm_local_type additive -pc_asm_blocks 1 -pc_asm_overlap 1 -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type bcgs -pc_type asm -pc_asm_type none -pc_asm_local_type additive -pc_asm_blocks 1 -pc_asm_overlap 1 -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type dgmres -ksp_gmres_restart 20 -ksp_dgmres_eigen 100 -ksp_dgmres_max_eigen 1000 -pc_type asm -pc_asm_type none -pc_asm_local_type additive -pc_asm_blocks 1 -pc_asm_overlap 1 -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type fbcgs -pc_type asm -pc_asm_type none -pc_asm_local_type additive -pc_asm_blocks 1 -pc_asm_overlap 1 -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type fbcgsr -pc_type asm -pc_asm_type none -pc_asm_local_type additive -pc_asm_blocks 1 -pc_asm_overlap 1 -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

#### GASM Preconditioning

# ./csd -ksp_type lgmres -ksp_gmres_restart 20 -ksp_lgmres_augment 10 -ksp_lgmres_constant -pc_type gasm -pc_gasm_type interpolate -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type bcgs -pc_type gasm -pc_gasm_type interpolate -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type dgmres -ksp_gmres_restart 20 -ksp_dgmres_eigen 100 -ksp_dgmres_max_eigen 1000 -pc_type gasm -pc_gasm_type interpolate -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type fbcgs -pc_type gasm -pc_gasm_type interpolate -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly

# ./csd -ksp_type fbcgsr -pc_type gasm -pc_gasm_type interpolate -sub_pctype ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly
