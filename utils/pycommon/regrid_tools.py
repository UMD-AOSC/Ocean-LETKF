import numpy as np
from numba import jit

@jit(nopython=True)
def fill_nan_grds_aux(v2d, aux2d, radius):
    nlat, nlon = v2d.shape
    wk2d = v2d.copy()

    n_replaced = 0

    for j in range(0,nlat):
        jpn = j + radius
        jmn = j - radius

        for i in range(0,nlon):
            ipn = i + radius
            imn = i - radius

            if not np.isnan( v2d[j,i] ):
                continue

            js = jmn if jmn>=0 else 0
            je = jpn if jpn<nlat else nlat-1
            n_total = 0
            n_nan = 0
            for j2 in range(js,je+1):
                for i2 in range(imn,ipn+1):
                    n_total += 1

                    if i2>=0 and i2<nlon:
                        i2_prdc = i2
                    elif i2<0:
                        i2_prdc = nlon + i2
                    elif i2>=nlon:
                        i2_prdc = i2 - nlon
                    else:
                        raise Exception("ERROR: fill_nan_grds_aux")

                    if np.isnan(v2d[j2,i2_prdc]):
                        n_nan += 1

            if n_nan == n_total:
                wk2d[j,i] = aux2d[j,i]
                n_replaced += 1


    print("n_replaced = ", n_replaced)
    return wk2d


# functions revised from https://github.com/NCAR/WOA_MOM6/blob/master/fill.py
# originally written by Gustavo Marques
# add an additional termination condition
@jit(nopython=True)
def iterative_fill_POP_core(var, fillmask, missing_value, tol=1.e-4, ltripole=True, nitermax = 10000, verbose=False):

    done = False
    niter = 0
    nlat,nlon = var.shape

    work = np.empty((nlat, nlon))
    while not done:
        done = True
        niter += 1
        if niter > nitermax:
            if verbose:
                print("[warning]: reach the maximum iteration: ", nitermax)
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

        if niter%3000 == 0 and verbose: 
           print("niter=",niter)

    print("niter_final=", niter)


@jit(nopython=True)
def iterative_fill_sor(nlat, nlon, var, fillmask, tol=5.0e-4,
            rc=1.6, max_iter=100, ltripole=False):
    """Iterative land fill algorithm via SOR solution of Laplace Equation."""

#    print("IN FOB version : _iterative_fill_sor")

    # Compute a zonal mean to use as a first guess
    # Apprarently jit doesn't like masked arrays so loop it out
    zoncnt = np.zeros(nlat)
    zonavg = np.zeros(nlat)
    for j in range(0,nlat) :
        zoncnt[j] = np.sum(np.where(fillmask[j,:],0,1))
        zonavg[j] = np.sum(np.where(fillmask[j,:],0,var[j,:]))
        if zoncnt[j] != 0 : zonavg[j] = zonavg[j]/zoncnt[j]

    # Fill missing zonal averages for rows that are entirely land
    for j in range(0,nlat-1) :   # northward pass
        if zoncnt[j] > 0 and zoncnt[j+1] == 0:
            zoncnt[j+1]=1
            zonavg[j+1] = zonavg[j]
    for j in range(0,nlat-1) :  # southward pass
        jrev = nlat-1-j
        if zoncnt[jrev] > 0 and zoncnt[jrev-1] == 0 :
            zoncnt[jrev-1]=1
            zonavg[jrev-1] = zonavg[jrev]

    # Replace the input array missing values with zonal average as first guess
    for j in range(0,nlat) :
        for i in range(0,nlon) :
            if fillmask[j,i] : var[j,i] = zonavg[j]

    # Now do the iterative 2D fill
    res = np.zeros((nlat,nlon))  # work array hold residuals
    res_max = tol
    iter0 = 0
    while iter0 < max_iter and res_max >= tol:
        res = res*0.0  # reset the residual to zero for this iteration

        # assume bottom row is all land, leave it set to zonal average
        # deal with top row separately below
        for j in range(1, nlat-1):
            jm1 = j - 1
            jp1 = j + 1

            for i in range(0, nlon):
                if fillmask[j, i]:
                    im1 = i - 1
                    if i == 0:                  # assume periodic in x
                        im1 = nlon - 1
                    ip1 = i + 1
                    if i == nlon - 1:
                        ip1 = 0

                    # this is SOR
                    res[j,i] = var[j,ip1] + var[j,im1] + var[jm1,i] + var[jp1,i] - 4.0*var[j,i]
                    var[j,i] = var[j,i] + rc*0.25*res[j,i]

        # do the top row if there was some valid data there in the input
        # otherwise leave it set to zonal average of northernmost row with valid data
        if  zoncnt[nlat-1] > 1 :
            j = nlat-1
            jm1 = j-1
            jp1 = j
            for i in range(0,nlon) :
                if fillmask[j,i] :
                    im1 = i-1
                    if i == 0:
                        im1 = nlon - 1
                    ip1 = i+1
                    if i == nlon - 1:
                        ip1 = 0
                    io = nlon-1-i

                    if ltripole :  # use cross-pole periodicity
                        res[j,i] = var[j,ip1] + var[j,im1] + var[jp1,io] + var[jm1,i] - 4.0*var[j,i]
                        var[j,i] = var[j,i] + rc*0.25*res[j,i]
                    else:          # do a 1D smooth on pole row
                        res[j,i] = var[j,ip1] + var[j,im1] - 2.0*var[j,i]
                        var[j,i] = var[j,i] + rc*0.5*res[j,i]


        res_max = np.max(np.abs(res))
        iter0 += 1
    
    print("SOR: niter, res_max=",iter0, res_max)
    return (iter0,res_max)
