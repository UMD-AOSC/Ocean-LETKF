MODULE vars_model

USE common,       ONLY: r_size
USE params_model, ONLY: nlon, nlat, nlev
USE params_model, ONLY: nv3d, nv2d

IMPLICIT NONE


PUBLIC

REAL(r_size),SAVE :: lon(nlon)
REAL(r_size),SAVE :: lat(nlat)
REAL(r_size),SAVE :: lev(nlev)                     !(OCEAN)
REAL(r_size),SAVE :: lon2d(nlon,nlat)              !(TRIPOLAR)
REAL(r_size),SAVE :: lat2d(nlon,nlat)              !(TRIPOLAR)

REAL(r_size),SAVE :: dx(nlon,nlat)
REAL(r_size),SAVE :: dy(nlon,nlat)
REAL(r_size),SAVE :: phi0(nlon,nlat)
REAL(r_size),SAVE :: kmt0(nlon,nlat)               !(OCEAN)

!STEVE: For generalized grid
REAL(r_size) :: lon0, lonf, lat0, latf
REAL(r_size) :: wrapgap

END MODULE vars_model
