MODULE params_model

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC


! MOM4 ncep2012 tripolar converted to spherical
  INTEGER,PARAMETER :: nlon=720 !720
  INTEGER,PARAMETER :: nlat=410 !360
  INTEGER,PARAMETER :: nlev=40 !7  !40
! MOM4 NCEP_om3_core3 test case, tripolar converted to spherical
! INTEGER,PARAMETER :: nlon=360
! INTEGER,PARAMETER :: nlat=200
! INTEGER,PARAMETER :: nlev=50 !50 !STEVE: trying to reduce grid size to check memory issue
! MOM4 box1 test case
! INTEGER,PARAMETER :: nlon=24
! INTEGER,PARAMETER :: nlat=35
! INTEGER,PARAMETER :: nlev=18
! MOM4 mk3p51 test case, Global Spherical grid
! INTEGER,PARAMETER :: nlon=192
! INTEGER,PARAMETER :: nlat=189
! INTEGER,PARAMETER :: nlev=31
! MOM4 iom1 test case, Indian Ocean
! INTEGER,PARAMETER :: nlon=150
! INTEGER,PARAMETER :: nlat=150
! INTEGER,PARAMETER :: nlev=28
! MOM4 tripolar om3_core1 or om3_core3 test cases.  (spherical lat/lon grids use the same dimensions)
! INTEGER,PARAMETER :: nlon=360
! INTEGER,PARAMETER :: nlat=200
! INTEGER,PARAMETER :: nlev=50

  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
  INTEGER,PARAMETER :: nv4d=0 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
! INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss          !(OCEAN)
! INTEGER,PARAMETER :: nv2d=7 ! ssh/t/s, + sfc fluxes: taux,tauy,heat,freshwater
  INTEGER,PARAMETER :: nv2d=4 ! ssh,sst,sss,eta      !(OCEAN) !(ALTIMETRY)
  INTEGER,PARAMETER :: nvsfc=0 !14

  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall

  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                !(OCEAN)
                                               !          From ocean_sbc.res.nc:
  INTEGER,PARAMETER :: iv2d_ssh=1              !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_sst=2              !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model
  INTEGER,PARAMETER :: iv2d_sss=3              !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models
  INTEGER,PARAMETER :: iv2d_eta=4              !(OCEAN) ! eta sea surface perturbation from mom4's ocean_barotropic.res.nc restart file
  INTEGER,PARAMETER :: iv4d_x=1                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_y=2                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_z=3                !(OCEAN) (DRIFTERS)

  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'

  ! For grid_spec.nc data file:
  CHARACTER(12) :: gridfile = 'grid_spec.nc'

  ! Bounds checking (for output by common_mom4.f90::write_restart)
  LOGICAL :: do_physlimit=.true.
  REAL(r_size) :: max_t = 40.0d0 ! ÂC
  REAL(r_size) :: min_t = -4.0d0 ! ÂC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu

END MODULE params_model
