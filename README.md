# Ocean-LETKF:  
Main repository for the Ocean-LETKF development

# Getting started with git:
https://gist.github.com/blackfalcon/8428401  
https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message  

# Configuration:  
All machine-specific details have been moved to 'config/'
For a new machine, create a set of files following the same format.

Code has been reorganized to merge support for multiple ocean models into one repository:
src/model_specific

# Parallel decomposition:  
The code has been updated to support the independent assimilation of subgrids of the global model grid.
To utilize this feature, the model must be compiled with the -DDynamic fortran pre-processor option.
This option will become the default in future builds.

Previously, the code read in the entire model grid and sent each grid point to a different processor, in an alternating fashion.
This is still done for each subgrid. Further load balancing could be applied.
(1) there may be a varying number of observations and ocean grid points in each subgrid
(2) the I/O for the model grid is the most expensive part (for large grids), so this cost has been greatly reduced (via parallelization),
(3) many models are already decomposed for executing the model integration, so this setup reduces intermediate processing required by the DA system
(4) The observation operator still produces global OMF files, so the subgrids still 'think' they are global

To improve on the load balancing, one might decide to:
(1) pre-process the observations to include only a halo around each subgrid. At the moment, all global observations are read in
(2) divy out panels to a node based on the number of observations in the panel (e.g. as a weight/cost function)

The observations are all read in and stored in a separate copy for each processor. This is fine with
few observations, but when the observation count gets large, this may become a memory constraint.
Instead, we should only have each cpu store enough observations that will be used by the gridpoints
assigned to that node (e.g. for example, within the local region plus a halo). So far, however,
it has been the model grids that have been the primary memory constraint as resolution increases.

# To correct merge conflict:
https://help.github.com/articles/resolving-a-merge-conflict-from-the-command-line/

