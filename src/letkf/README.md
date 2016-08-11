This directory contains the main letkf program.

letkf.f90 is the main program

letkf_tools.f90 is the primary working routine that loops through all model gridpoints

letkf_local.f90 performs the localization scheme
NOTE: this is currently a pointer to either letkf_local.f90.kdtree or letkf_local.f90.bruteforce

letkf_obs.f90 processes the observations

params_letkf.f90 is a module header that has the main parameters that may be updated at runtime using input.nml

input.nml contains changes to the default parameters at runtime

vars_letkf.f90 is a module header that contains the variables used throughout most of the letkf code

vars_obs.f90 contains the arrays for the observation data 
