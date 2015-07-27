# Ocean-LETKF:  
Main repository for the Ocean-LETKF development

<<<<<<< HEAD
# Configuration:  
I am attempting to move all machine-specific details to 'config/'
For a new machine, create a set of files following the same format.


# Parallel decomposition:  
Currently, the code reads in the entire model grid and sends each grid point to a different processor, in an alternating fashion.
This is inefficient for the ocean, since   
(1) there are often very few observations,
(2) the I/O for the model grid is the most expensive part (for large grids),
(3) the MOM6 model is already decomposed into multiple files for each pe panel for the model forecast

For mom6, it would be better to read in each panel either:
(1) only if that core is using that grid point, or
(2) divy out panels based on the number of observations in the panel (e.g. as a weight/cost function)

The observations are all read in and stored in a separate copy for each processor. This is fine with
few observations, but when the observation count gets large, this is a major constraint.
Instead, we should only have each cpu store enough observations that will be used by the gridpoints
assigned to that node.

* The current strategy is to perform an even load balance across all processors.
**The new strategy must be to minimize node-to-node communications.
