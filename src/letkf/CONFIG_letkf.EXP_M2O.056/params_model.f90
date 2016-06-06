MODULE params_model

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

  ! Now definable via namelist at runtime:

  ! MOM4 ncep2012 tripolar converted to spherical
  INTEGER,PARAMETER :: nlon=720
  INTEGER,PARAMETER :: nlat=410
  INTEGER,PARAMETER :: nlev=40

  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
! INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss          !(OCEAN)
! INTEGER,PARAMETER :: nv2d=7 ! ssh/t/s, + sfc fluxes: taux,tauy,heat,freshwater
  INTEGER,PARAMETER :: nv2d=5 ! ssh,sst,sss,eta,mld  !(OCEAN) !(ALTIMETRY)
  INTEGER,PARAMETER :: nv4d=0 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS

  ! Initialize via subroutine below:
  INTEGER,SAVE :: nij0           ! Number of gridpoints handled by myrank processor
  INTEGER,SAVE :: nlevall        ! Total number of variables and levels (3d + 2d)
  INTEGER,SAVE :: ngpv           ! Total number of gridpoints, including nij0*nlevall

  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                !(OCEAN)
                                               !          From ocean_sbc.res.nc:
  INTEGER,PARAMETER :: iv2d_ssh=1              !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_sst=2              !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model
  INTEGER,PARAMETER :: iv2d_sss=3              !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models
  INTEGER,PARAMETER :: iv2d_eta=4              !(OCEAN) ! eta sea surface perturbation from mom4's ocean_barotropic.res.nc restart file
  INTEGER,PARAMETER :: iv2d_mld=5              !(OCEAN) ! mixed layer depth
  INTEGER,PARAMETER :: iv4d_x=1                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_y=2                !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_z=3                !(OCEAN) (DRIFTERS)

  !
  ! Elements
  !
  CHARACTER(4) :: element(nv3d+nv2d+nv4d)
! element(iv3d_u) = 'U   '
! element(iv3d_v) = 'V   '
! element(iv3d_t) = 'T   '
! element(iv3d_s) = 'S   '               !(OCEAN)
! element(nv3d+iv2d_ssh) = 'SSH '        !(OCEAN)
! element(nv3d+iv2d_sst) = 'SST '        !(OCEAN)
! element(nv3d+iv2d_sss) = 'SSS '        !(OCEAN)
! if (DO_ALTIMETRY) then
!   element(nv3d+iv2d_eta) = 'eta '      !(OCEAN)
! endif
! if (DO_MLD) then
!   element(nv3d+iv2d_mld) = 'mld'       !(OCEAN)
! endif
! if (DO_DRIFTERS) then
!   element(nv3d+nv2d+iv4d_x) = 'X   '             !(OCEAN) (DRIFTERS)
!   element(nv3d+nv2d+iv4d_y) = 'Y   '             !(OCEAN) (DRIFTERS)
!   element(nv3d+nv2d+iv4d_z) = 'Z   '             !(OCEAN) (DRIFTERS)
! endif

  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
  CHARACTER(32) :: ts_basefile = 'ocean_temp_salt.res.nc'
  CHARACTER(32) :: uv_basefile = 'ocean_velocity.res.nc'
  CHARACTER(32) :: sf_basefile = 'ocean_sbc.res.nc'
  CHARACTER(32) :: sh_basefile = 'ocean_barotropic.res.nc'
  CHARACTER(32) :: hs_basefile = 'ocean_TS.nc'

  ! For grid_spec.nc data file:
  CHARACTER(12) :: gridfile = 'grid_spec.nc'

  ! Bounds checking (for output by common_mom4.f90::write_restart)
  LOGICAL :: do_physlimit=.true.
  REAL(r_size) :: max_t = 40.0d0 ! ÂC
  REAL(r_size) :: min_t = -4.0d0 ! ÂC
  REAL(r_size) :: max_s = 50.0d0 ! psu
  REAL(r_size) :: min_s =  0.0d0 ! psu

CONTAINS

SUBROUTINE initialize_params_model
  nij0=nlon*nlat
  nlevall=nlev*nv3d+nv2d
  ngpv=nij0*nlevall
END SUBROUTINE initialize_params_model

END MODULE params_model
