MODULE params_letkf

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

!-------------------------------------------------------------------------------
! Ensemble Size
!-------------------------------------------------------------------------------
INTEGER,PARAMETER :: nbv=10


!-------------------------------------------------------------------------------
! From letkf_local.f90
!-------------------------------------------------------------------------------

!STEVE: Testing "Vertical Tube" localization:
!       i.e. the localization is not applied vertically
! This provides the benefit that 
! (1) the analysis only has to be computed once
! per horizontal gridpoint, thus providing a nlevX (40X) speedup
! (2) the altimetry, SST, SSH, and bottom pressure (GRACE) can now be applied
! as direct constraints on the water column.
!
! There is precedence for this as in the paper "Reconstructing the Ocean's
! Interior from Surface Data" Wang et al. (2013)
!
!STEVE: making namelist accessible:
LOGICAL :: DO_NO_VERT_LOC=.true. !STEVE: moved to letkf_local.f90 from letkf_tools.f90
INTEGER :: localization_method=1 !1 !(OCEAN) =0 for uniform radius (default), =1 for latitude-dependent

!-------------------------------------------------------------------------------
! From letkf_tools.f90
!-------------------------------------------------------------------------------

! > 0: globally constant covariance inflation
! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
REAL(r_size) :: cov_infl_mul = 1.0d0 !(NO INFLATION) !-1.0d0 => adaptive multiplicative inflation
REAL(r_size) :: sp_infl_add = 0.d0 !additive inflation
LOGICAL :: DO_INFL_RESET = .true. !(DO_SFCFLUXES) STEVE: for inflating only the surface forcing fields. Reset the 3D ocean interior to 0% inflation (rho = 1.0)


!-------------------------------------------------------------------------------
! From common_mom4.f90
!-------------------------------------------------------------------------------
! For (DRIFTERS)
LOGICAL :: DO_DRIFTERS=.false.
! For (SFCFLUX)
LOGICAL :: DO_SFCFLUX=.false.
! For (ALTIMETRY)
LOGICAL :: DO_ALTIMETRY=.false.

END MODULE params_letkf
