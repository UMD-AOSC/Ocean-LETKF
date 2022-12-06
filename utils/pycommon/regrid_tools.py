import numpy as np
from numba import jit

# functions from https://github.com/NCAR/WOA_MOM6/blob/master/fill.py
# originally written by Gustavo Marques
@jit(nopython=True)
def iterative_fill_POP_core(var, fillmask, missing_value, tol=1.e-4, ltripole=True):

    done = False
    niter = 0
    nitermax = 100000
    nlat,nlon = var.shape

    work = np.empty((nlat, nlon))
    while not done:
        done = True
        niter += 1
        if niter > nitermax:
            break

        # assume bottom row is land, so skip it
        for j in range(1, nlat):
            jm1 = j - 1
            jp1 = j + 1

            for i in range(0, nlon):

                # assume periodic in x
                im1 = i - 1
                if i == 0:
                    im1 = nlon - 1
                ip1 = i + 1
                if i == nlon - 1:
                    ip1 = 0

                work[j, i] = var[j, i]

                if not fillmask[j, i]:
                    continue

                numer = 0.0
                denom = 0.0

                # East
                if var[j, ip1] != missing_value:
                    numer += var[j, ip1]
                    denom += 1.0

                # North
                if j < nlat - 1:
                    if var[jp1, i] != missing_value:
                        numer += var[jp1, i]
                        denom += 1.0

                else:
                    # assume only tripole has non-land top row
                    if ltripole:
                        if var[j, nlon - 1 - i] != missing_value:
                            numer += var[j, nlon - 1 - i]
                            denom += 1.0

                # West
                if var[j, im1] != missing_value:
                    numer += var[j, im1]
                    denom += 1.0

                # South
                if var[jm1, i] != missing_value:
                    numer += var[jm1, i]
                    denom += 1.0

                # self
                if var[j, i] != missing_value:
                    numer += denom * var[j, i]
                    denom *= 2.0

                if denom > 0.0:
                    work[j, i] = numer / denom
                    if var[j, i] == missing_value:
                        done = False
                    else:
                        delta = np.fabs(var[j, i] - work[j, i])
                        if delta > tol * np.abs(var[j, i]):
                            done = False

        var[1:nlat, :] = work[1:nlat, :]

        if niter%1000 == 0: 
           print("niter=",niter)

    print("niter_final=", niter)

