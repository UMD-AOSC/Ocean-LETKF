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
Currently, the code reads in the entire model grid and sends each grid point to a different processor, in an alternating fashion.
This is less efficient for the ocean, since   
(1) there are often very few observations,
(2) the I/O for the model grid is the most expensive part (for large grids),
(3) the MOM6 model is already decomposed into multiple files for each pe panel for the model forecast

For mom6, it would be better to read in the entire model grid in subgrid panels either:
(1) only if that core is using that grid point, or
(2) divy out panels based on the number of observations in the panel (e.g. as a weight/cost function)

The observations are all read in and stored in a separate copy for each processor. This is fine with
few observations, but when the observation count gets large, this is a major constraint.
Instead, we should only have each cpu store enough observations that will be used by the gridpoints
assigned to that node (e.g. for example, within the local region plus a halo).

* The current strategy is to perform an even load balance across all processors.
**The new strategy must be to minimize node-to-node communications.

# To correct merge conflict:
https://help.github.com/articles/resolving-a-merge-conflict-from-the-command-line/

