MODULE params_letkf

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

!-------------------------------------------------------------------------------
! Ensemble Size
!-------------------------------------------------------------------------------
INTEGER :: nbv=56

!-------------------------------------------------------------------------------
! Details for LETKF localization and quality control
!-------------------------------------------------------------------------------
!STEVE: making these namelist accessible:
INTEGER,SAVE :: nslots=5                  ! number of time slots for 4D-LETKF
INTEGER,SAVE :: nbslot=5 !1               !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
REAL(r_size),SAVE :: sigma_obs=720.0d3    !3x Rossby Radius of Deformation at Equ. according to Chelton
REAL(r_size),SAVE :: sigma_obs0=200.0d3   !20x Rossby Radius of Deformation at Pole, according to Chelton
REAL(r_size),SAVE :: sigma_obsv=-1.0d0    !STEVE: use a negative value with option "DO_NO_VERT_LOC"
REAL(r_size),SAVE :: sigma_obst=5.0d0     ! Not using this at the moment
REAL(r_size),SAVE :: gross_error=3.0d0    ! number of standard deviations
! REAL(r_size),PARAMETER :: gross_error=10.0d0 !3.0d0 ! number of standard deviations   (Use for OSSEs)
                                                      ! used to filter out observations

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

!-------------------------------------------------------------------------------
! From common_mom4.f90
!-------------------------------------------------------------------------------
! For (MLD, assimilate the mixed layer depth)
LOGICAL :: DO_MLD = .false.
! For (DRIFTERS)
LOGICAL :: DO_DRIFTERS=.false.
! For (ALTIMETRY)
LOGICAL :: DO_ALTIMETRY=.false.
! For (ALTIMETRY, Sea Level Anomalies)
LOGICAL :: DO_SLA=.false.
! For (ALTIMETRY, Absolute Dynamic Topography)
LOGICAL :: DO_ADT=.false.
! For (TRIPLOAR)
LOGICAL :: DO_TRIPOLAR=.true.
! For (IRREG_GRID) for non-orthogonal/non-rectilinear grids (not fully supported)
LOGICAL :: DO_IRREG_GRID=.false.

END MODULE params_letkf
