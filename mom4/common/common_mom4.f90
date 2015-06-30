MODULE common_mom4
!sfc=======================================================================
!
! [PURPOSE:] Common Information for MOM4
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!
!=======================================================================
!$USE OMP_LIB
  USE common
! use isa, ONLY: isnan, isnan4 !STEVE: for debugging (isnan, isinf)
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: slen=256
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
!
  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
  INTEGER,PARAMETER :: nv4d=3 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
! INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss          !(OCEAN)
! INTEGER,PARAMETER :: nv2d=7 ! ssh/t/s, + sfc fluxes: taux,tauy,heat,freshwater
  INTEGER,PARAMETER :: nv2d=3!4 ! ssh,sst,sss,eta      !(OCEAN) !(ALTIMETRY)
  INTEGER,PARAMETER :: nvsfc=0 !14
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                      !(OCEAN)
                                                     !          From ocean_sbc.res.nc:
  INTEGER,PARAMETER :: iv2d_ssh=1                    !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_sst=2                    !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model
  INTEGER,PARAMETER :: iv2d_sss=3                    !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models
  INTEGER,PARAMETER :: iv2d_eta=4                    !(OCEAN) ! eta sea surface perturbation from mom4's ocean_barotropic.res.nc restart file
  INTEGER,PARAMETER :: iv4d_x=1                      !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_y=2                      !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_z=3                      !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv2d_taux=5                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_tauy=6                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_tflx=7                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_qflx=8                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_u10=9                    !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_v10=10                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_t2m=11                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_q2m=12                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_pres=13                  !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_prate=14                 !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_dlw=15                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_dsw=16                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_nlw=17                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: iv2d_nsw=18                   !(OCEAN) (SFCFLUX)
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),SAVE :: lon(nlon)
  REAL(r_size),SAVE :: lat(nlat)
  REAL(r_size),SAVE :: lev(nlev)                     !(OCEAN)
  REAL(r_size),SAVE :: dx(nlon,nlat)
  REAL(r_size),SAVE :: dy(nlon,nlat)
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  REAL(r_size),SAVE :: kmt0(nlon,nlat)               !(OCEAN)
  REAL(r_size),SAVE :: wet(nlon,nlat)                !(OCEAN)
  REAL(r_size),SAVE :: area_t(nlon,nlat)             !(OCEAN)
  CHARACTER(4),SAVE :: element(nv3d+nv2d+nv4d)
  INTEGER, DIMENSION(nlon,nlat), SAVE     :: kmt=-1  !(OCEAN) STEVE: the bottom topography for mom4
  !STEVE: for Custom Localization
  INTEGER :: nobids(nij0*nlev)                       !(OCEAN)
  !STEVE: for verifying against input netcdf file
  INTEGER :: nlon0=0, nlat0=0, nlev0=0               !(OCEAN)
  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: vmax = 1.0e18
  !STEVE: For generalized grid
  REAL(r_size) :: lon0, lonf, lat0, latf
  REAL(r_size) :: wrapgap
  ! For (DRIFTERS)
  LOGICAL :: DO_DRIFTERS=.false.
  ! For (SFCFLUX)
  LOGICAL :: DO_SFCFLUX=.false.
  ! For (ALTIMETRY)
  LOGICAL :: DO_ALTIMETRY=.false.
  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
  REAL(r_size), DIMENSION(nlon,nlat) :: SSHclm_m
  REAL(r_sngl) :: buf4_2d(nlon,nlat)
  ! For grid_spec.nc data file:
  CHARACTER(12) :: gridfile = 'grid_spec.nc'
! For AMOC computation
  REAL(r_size) :: zb(nlev)
  REAL(r_size) :: dz(nlev)

! sfc fluxes: 
! taux: zonal surface stress
! tauy: meridional surface stress
! tflux: temperature flux
! qflux: freshwater flux
! u10: 10-meter zonal surface wind velocity
! v10: 10-meter meridional surface wind velocity
! t2m: 2 meter height temperature
! q2m: 2 meter height relative humidity
! prate: precipitation rate 
! pres: mean sea level pressure
! dlw: downward longwave radiation
! dsw: downward shortwave radiation
! longwv: net longwave radiation
! shrtwv: net shortwave radiation
  CHARACTER(slen), DIMENSION(14) :: sfc_infiles = &
(/'SFC_000_daily_TAUX.nc','SFC_000_daily_TAUY.nc','SFC_000_daily_TFLUX.nc','SFC_000_daily_QFLUX.nc','SFC_000_daily_U10.nc','SFC_000_daily_V10.nc', &
'SFC_000_daily_t2m.nc','SFC_000_daily_q2m.nc','SFC_000_daily_pres.nc','SFC_000_daily_PRATE.nc','SFC_000_daily_dlw.nc','SFC_000_daily_dsw.nc', &
'SFC_000_daily_LONGWV.nc','SFC_000_daily_SHRTWV.nc'/)
  CHARACTER(slen), DIMENSION(14) :: sfc_outfiles = &
(/'SFA_000_daily_TAUX.nc','SFA_000_daily_TAUY.nc','SFA_000_daily_TFLUX.nc','SFA_000_daily_QFLUX.nc','SFA_000_daily_U10.nc','SFA_000_daily_V10.nc', &
'SFA_000_daily_t2m.nc','SFA_000_daily_q2m.nc','SFA_000_daily_pres.nc','SFA_000_daily_PRATE.nc','SFA_000_daily_dlw.nc','SFA_000_daily_dsw.nc', &
'SFA_000_daily_LONGWV.nc','SFA_000_daily_SHRTWV.nc'/)
  CHARACTER(slen), DIMENSION(14) :: sfc_names = &
(/'uflx','vflx','tflux','qflux','u10','v10','t2m','q2m','PRESmsl','prate','dlw','dsw','longwv','shrtwv'/)
  CHARACTER(32) :: coeff_s2mfile = 'coeff_s2m.nc'
  CHARACTER(32) :: coeff_m2sfile = 'coeff_m2s.nc'
  INTEGER, PARAMETER :: nslot = 5
  INTEGER :: islot
  LOGICAL :: coeff_file_exists
  INTEGER, SAVE :: nslon=-1,nslat=-1
  INTEGER :: nmlon=nlon,nmlat=nlat
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: slon, slat, mlon, mlat
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: xc,yc
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: sfc_data
  INTEGER, DIMENSION(:), ALLOCATABLE :: xi,yi

  !STEVE: for debugging
  LOGICAL :: dodebug = .false.

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_mom4
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex

  WRITE(6,'(A)') 'Hello from set_common_mom4'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_s) = 'S   '             !(OCEAN)
  element(nv3d+iv2d_ssh) = 'SSH '      !(OCEAN)
  element(nv3d+iv2d_sst) = 'SST '      !(OCEAN)
  element(nv3d+iv2d_sss) = 'SSS '      !(OCEAN)
  if (DO_ALTIMETRY) then
    element(nv3d+iv2d_eta) = 'eta '      !(OCEAN)
  endif
  if (DO_SFCFLUX) then
    element(nv3d+iv2d_taux) = 'TAUX'   !(OCEAN)
    element(nv3d+iv2d_tauy) = 'TAUY'   !(OCEAN)
    element(nv3d+iv2d_tflx) = 'TFLX'  !(OCEAN)
    element(nv3d+iv2d_qflx) = 'QFLX'  !(OCEAN)
    element(nv3d+iv2d_u10) = 'U10 '     !(OCEAN)
    element(nv3d+iv2d_v10) = 'V10 '     !(OCEAN)
    element(nv3d+iv2d_t2m) = 'T2M '     !(OCEAN)
    element(nv3d+iv2d_q2m) = 'Q2M '     !(OCEAN)
    element(nv3d+iv2d_pres) = 'PRES'   !(OCEAN)
    element(nv3d+iv2d_prate) = 'PRAT' !(OCEAN)
    element(nv3d+iv2d_dlw) = 'DLW '     !(OCEAN)
    element(nv3d+iv2d_dsw) = 'DSW '     !(OCEAN)
    element(nv3d+iv2d_nlw) = 'NLWV'  !(OCEAN)
    element(nv3d+iv2d_nsw) = 'NSWV'  !(OCEAN)
  endif
  if (DO_DRIFTERS) then
    element(nv3d+nv2d+iv4d_x) = 'X   '             !(OCEAN) (DRIFTERS)
    element(nv3d+nv2d+iv4d_y) = 'Y   '             !(OCEAN) (DRIFTERS)
    element(nv3d+nv2d+iv4d_z) = 'Z   '             !(OCEAN) (DRIFTERS)
  endif
  if (DO_ALTIMETRY) then
    INQUIRE(FILE=trim(SSHclm_file),EXIST=ex)
    IF(ex) THEN
      ! Read in the model climatology
      CALL read_etaclm(SSHclm_file,SSHclm_m)
    ELSE
      WRITE(6,*) "The file does not exist: ", SSHclm_file
      WRITE(6,*) "Exiting common_mom4.f90..."
      STOP(1)
    ENDIF
  endif

  !
  ! Lon, Lat, f, orography
  !
!STEVE: this part adapted from ROMS, update from MOM4 netcdf files:
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile 
    WRITE(6,*) "Exiting common_mom4.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile
  call check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )
  call check( NF90_INQ_VARID(ncid,'grid_x_T',varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon) )
  WRITE(6,*) "lon(1) = ", lon(1)
  WRITE(6,*) "lon(nlon) = ", lon(nlon)
  call check( NF90_INQ_VARID(ncid,'grid_y_T',varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat) )
  WRITE(6,*) "lat(1) = ", lat(1)
  WRITE(6,*) "lat(nlat) = ", lat(nlat)
  call check( NF90_INQ_VARID(ncid,'zt',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)
! call check( NF90_INQ_VARID(ncid,'num_levels',varid) ) ! number of vertical levels
! call check( NF90_GET_VAR(ncid,varid,phi0) )
! WRITE(6,*) "ncid = ", ncid
! WRITE(6,*) "varid = ", varid
! WRITE(6,*) "phi0(1,1) = ", phi0(1,1)
! WRITE(6,*) "phi0(nlon,nlat) = ", phi0(nlon,nlat)
  !
  ! dx and dy
  !
  call check( NF90_INQ_VARID(ncid,'ds_01_21_T',varid) )    ! width of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dx) ) 
  call check( NF90_INQ_VARID(ncid,'ds_10_12_T',varid) )    ! height of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dy) ) 
  call check( NF90_INQ_VARID(ncid,'area_T',varid) )        ! area of T_cell
  call check( NF90_GET_VAR(ncid,varid,area_t) ) 
  WRITE(6,*) "common_mom4:: grid_spec.nc MIN(dx) = ", MINVAL(dx)
  WRITE(6,*) "common_mom4:: grid_spec.nc MAX(dx) = ", MAXVAL(dx)
  WRITE(6,*) "common_mom4:: grid_spec.nc MIN(dy) = ", MINVAL(dy)
  WRITE(6,*) "common_mom4:: grid_spec.nc MAX(dy) = ", MAXVAL(dy)
  WRITE(6,*) "common_mom4:: grid_spec.nc MIN(area_t) = ", MINVAL(area_t)
  WRITE(6,*) "common_mom4:: grid_spec.nc MAX(area_t) = ", MAXVAL(area_t)

  !
  ! kmt data
  !
  call check( NF90_INQ_VARID(ncid,'num_levels',varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,kmt0) )
  WRITE(6,*) "kmt0(1,1) = ", kmt0(1,1)
  WRITE(6,*) "kmt0(nlon,nlat) = ", kmt0(nlon,nlat)
  kmt = NINT(kmt0)
  call check( NF90_INQ_VARID(ncid,'wet',varid) )        ! land/sea flag (0=land) for T-cell
  call check( NF90_GET_VAR(ncid,varid,wet) )
  WRITE(6,*) "wet(1,1) = ", wet(1,1)
  WRITE(6,*) "wet(nlon,nlat) = ", wet(nlon,nlat)

  WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile
  WRITE(6,*) "dx(1,1) = ", dx(1,1)
  WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
  WRITE(6,*) "dy(1,1) = ", dy(1,1)
  WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)


  !STEVE: needed for computing the AMOC based on the streamfunction calculation:
  call check( NF90_INQ_VARID(ncid,'zb',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,zb) )
  WRITE(6,*) "zb(1) = ", zb(1)
  WRITE(6,*) "zb(nlev) = ", zb(nlev)

  ! Compute dz:
  dz(1) = zb(1)
  do k=2,nlev
    dz(k) = zb(k)-zb(k-1)
  enddo

  !
  ! Corioris parameter
  !
!$OMP PARALLEL WORKSHARE
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!$OMP END PARALLEL WORKSHARE

  ! Close the grid_spec.nc file:
  call check( NF90_CLOSE(ncid) )

  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

  RETURN
END SUBROUTINE set_common_mom4
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
SUBROUTINE read_etaclm(SSHclm_file,SSHclm_m)
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: SSHclm_file
  REAL(r_size), INTENT(OUT) :: SSHclm_m(nlon,nlat)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: i,j
  INTEGER :: ncid, varid

  ! read the model SSH climatology netcdf file
  ! read into: SSHclm_m
  call check( NF90_OPEN(SSHclm_file,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file

  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'ssh',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  DO j=1,nlat
    DO i=1,nlon
      !STEVE: Hopefully reading in meters here... (data might be in cm)
      SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
    END DO
  END DO

  call check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm

!STEVE: add this:
!-- Read a grid file in mom4 netcdf format ---------------------------------------------------

!-----------------------------------------------------------------------
! File I/O (netCDF) modified from ROMS
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(infile,v3d,v2d)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  CHARACTER(slen) :: tsfile,uvfile, sffile, bfile, drfile ! (TS) (UV) (SFC) (barotropic - eta) (DRIFTERS)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  !STEVE:
  REAL(r_size) :: meanSSH !STEVE: for temporary SSH estimate based on heat content
  REAL(r_size) :: videpth !STEVE: depth of vertically integrated heat content
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
  !LOGICAL :: dodebug = .true.
  CHARACTER(3) :: MEM3
  CHARACTER(32) :: sfc_infile

  tsfile = trim(infile)//'.ocean_temp_salt.res.nc'
  uvfile = trim(infile)//'.ocean_velocity.res.nc'
  sffile = trim(infile)//'.ocean_sbc.res.nc'
  bfile  = trim(infile)//'.ocean_barotropic.res.nc'
! read (FLUXES)

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the T/S netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_grd:: just opened file ", tsfile

  !!! t
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'temp',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_grd:: just got data for variable temp"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_t) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable temp"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_grd:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ",MAXVAL(v3d(:,:,k,iv3d_t))
    enddo
  endif
! !STEVE: end

  !!! s
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'salt',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_grd:: just got data for variable salt"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_s) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable salt"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_grd:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_s) = ", MAXVAL(v3d(:,:,k,iv3d_s))
    enddo 
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the U/V netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )
  IF(istat /= NF90_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR in read_grd for ',uvfile
    STOP(7)
  END IF
  WRITE(6,*) "read_grd:: just opened file ", uvfile

  !!! u
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_grd:: just got data for variable u"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_u) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable u"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_grd:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ",MAXVAL(v3d(:,:,k,iv3d_u))
    enddo
  endif
! !STEVE: end

  !!! v
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_grd:: just got data for variable v"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_v) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable v"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_grd:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_v) = ", MAXVAL(v3d(:,:,k,iv3d_v))
    enddo 
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set the SST, SSS, and SSH data for the SFC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.false.) then
    v2d(:,:,iv2d_sst) = v3d(:,:,ilev_sfc,iv3d_t)
    v2d(:,:,iv2d_sss) = v3d(:,:,ilev_sfc,iv3d_s)
  else
    call check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
    WRITE(6,*) "read_grd:: just opened file ", sffile

    !!! SST
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'t_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_grd:: just got data for variable sfc temp"
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. 1) v2d(i,j,iv2d_sst) = REAL(buf4(i,j,ilev_sfc),r_size) - t0c !kelvin
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable SST"

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SST"
      WRITE(6,*) "read_grd:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sst) = ",MAXVAL(v2d(:,:,iv2d_sst))
    endif
! !STEVE: end

    !!! SSS
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'s_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_grd:: just got data for variable sfc salt"
    DO j=1,nlat
      DO i=1,nlon
        v2d(i,j,iv2d_sss) = REAL(buf4(i,j,ilev_sfc),r_size)
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable SSS"

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SSS"
      WRITE(6,*) "read_grd:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sss) = ", MAXVAL(v2d(:,:,iv2d_sss))
    endif
! !STEVE: end

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Use the sbc SSH
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.false.) then
      !STEVE: for now, use vertically integrated heat content
      videpth = 300
      v2d(:,:,iv2d_ssh) = v3d(:,:,1,iv3d_t)*lev(1)
      do k=2,nlev
        if (lev(k) > videpth) EXIT
        v2d(:,:,iv2d_ssh) = v2d(:,:,iv2d_ssh) + v3d(:,:,k,iv3d_t)*(lev(k)-lev(k-1)) 
      enddo
  
      ! Divide out averages
      meanSSH = SUM(v2d(:,:,iv2d_ssh))/(nlon*nlat)
      where(v2d(:,:,iv2d_ssh) > 0) &
        v2d(:,:,iv2d_ssh) = v2d(:,:,iv2d_ssh)-meanSSH

    else !STEVE: use the ocean_sbc.res.nc file
      !!! SSH
      buf4=0.0
      call check( NF90_INQ_VARID(ncid,'sea_lev',varid) )
      call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
      if (dodebug) WRITE(6,*) "read_grd:: just got data for variable sea_lev"
      DO j=1,nlat
        DO i=1,nlon
          !STEVE: Hopefully reading in meters here... (data might be in cm)
          v2d(i,j,iv2d_ssh) = REAL(buf4(i,j,ilev_sfc),r_size)
        END DO
      END DO
      if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable SSH"

! !STEVE: debug
      if (dodebug) then
        WRITE(6,*) "POST-SSH"
        WRITE(6,*) "read_grd:: sffile = ", sffile
        WRITE(6,*) "max val for level v2d(:,:,iv2d_ssh) = ", MAXVAL(v2d(:,:,iv2d_ssh))
      endif
! !STEVE: end

    endif
    call check( NF90_CLOSE(ncid) )
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the ALTIMETRY netcdf restart file (eta)
  ! (These are the modeled sfc height perturbations used by GODAS)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  altimetry : if(DO_ALTIMETRY) then
    !STEVE: use the sea level perturbation from ocean_barotropic.res.nc
    call check( NF90_OPEN(bfile,NF90_NOWRITE,ncid) )
    WRITE(6,*) "read_grd:: just opened file ", bfile

    !!! SSH
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'eta_t',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_grd:: just got data for variable eta_t"
    DO j=1,nlat
      DO i=1,nlon
        !STEVE: Hopefully reading in meters here... (data might be in cm)
        v2d(i,j,iv2d_eta) = REAL(buf4(i,j,ilev_sfc),r_size)
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable SSH"

    ! Convert SSH eta stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,iv2d_eta) = v2d(:,:,iv2d_eta) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_grd:: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
      WRITE(6,*) "min val for level v2d(:,:,iv2d_eta) = ", MINVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  else
    WRITE(6,*) "read_grd:: Skipping SFC eta from: ", bfile
  endif altimetry

  ! For additional variables:
  ! E.g. surface fluxes, drifters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the surface fluxes netcdf file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !STEVE: debug test
  if (DO_SFCFLUX) then
    read( infile(3:4), '(i2)' )  islot
    MEM3 = infile(5:7)
    print *, "======================================================"
    print *, "read_grd :: DO_SFCFLUX ==============================="
    print *, "infile = ", infile
    print *, "islot = ", islot
    sfc_infile = sfc_infiles(1)
    sfc_infile(5:7) = MEM3
    print *, "sfc_infile = ", sfc_infile

    ! Read coefficients file
    INQUIRE(FILE=coeff_s2mfile, EXIST=coeff_file_exists)
    if (coeff_file_exists) then
      print *, "read_grd :: call read_sfc_grid ==============================="
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat) !STEVE: just need nslon and nslat
      print *, "read_grd :: call read_coeff ==============================="
      CALL read_coeff(coeff_s2mfile,xi,yi,xc,yc,mlon,mlat,nmlon,nmlat)
    else
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat)
      CALL coeff_sfc2model(slon,slat,lon,lat,xi,yi,xc,yc)
    endif

    ALLOCATE(sfc_data(nslon,nslat))

    ! Open and Read each of the SFC data files
    do i=1,nvsfc
      print *, "i,nvsfc = ", i,'/',nvsfc
      sfc_infile = sfc_infiles(i)
      sfc_infile(5:7) = MEM3
      print *, "sfc_infile = ", sfc_infile
      print *, "read_grd :: call read_sfc_data ==============================="
      CALL read_sfc_data(sfc_infile,sfc_names(i),islot+1,sfc_data) !STEVE: if using RIP, then use islot instead of islot+1

      ! Interpolate to the model grid using coefficients
      ! Convert sfc_data to model grid (on v2d)
      print *, "read_grd :: call itpl_sfc2model ==============================="
      CALL itpl_sfc2model(sfc_data,xi,yi,xc,yc,v2d(:,:,nv2d-nvsfc+i))
    enddo

    print *, "read_grd :: DEALLOCATE ==============================="
    DEALLOCATE(sfc_data)
    DEALLOCATE(xi,yi,xc,yc,slon,slat)
!   STOP(1)
  endif

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grd4(infile,v3d,v2d)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY: ocean_barotropic.res.nc, contains eta)
  INTEGER :: ncid,istat,varid
  !STEVE:
  INTEGER :: i,j,k
  REAL(r_size) :: meanSSH !STEVE: for temporary SSH estimate based on heat content
  REAL(r_size) :: videpth !STEVE: depth of vertically integrated heat content
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
! LOGICAL :: dodebug = .true.
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: model_data
  CHARACTER(3) :: MEM3
  CHARACTER(32) :: sfc_infile

  tsfile = trim(infile)//'.ocean_temp_salt.res.nc'
  uvfile = trim(infile)//'.ocean_velocity.res.nc'
  sffile = trim(infile)//'.ocean_sbc.res.nc'
  bfile  = trim(infile)//'.ocean_barotropic.res.nc'

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  call check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )

  !!! t
  call check( NF90_INQ_VARID(ncid,'temp',varid) )
  call check( NF90_GET_VAR(ncid,varid,v3d(:,:,:,iv3d_t)) )

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_grd4:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ", MAXVAL(v3d(:,:,k,iv3d_t))
    enddo 
  endif
! !STEVE: end

  !!! s
  call check( NF90_INQ_VARID(ncid,'salt',varid) )
  call check( NF90_GET_VAR(ncid,varid,v3d(:,:,:,iv3d_s)) )

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_grd4:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_s) = ", MAXVAL(v3d(:,:,k,iv3d_s))
    enddo 
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  call check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )

  !!! u
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  call check( NF90_GET_VAR(ncid,varid,v3d(:,:,:,iv3d_u)) )
! v3d(nlon,:,:,iv3d_u) = 0.0 !STEVE: why was this?

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_grd4:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ", MAXVAL(v3d(:,:,k,iv3d_u))
    enddo 
  endif
! !STEVE: end

  !!! v
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  call check( NF90_GET_VAR(ncid,varid,v3d(:,:,:,iv3d_v)) )
! v3d(:,nlat,:,iv3d_v) = 0.0 !STEVE: why was this?

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_grd4:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_v) = ", MAXVAL(v3d(:,:,k,iv3d_v))
    enddo 
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set the SST and SSS data for the SFC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.false.) then
    v2d(:,:,iv2d_sst) = v3d(:,:,ilev_sfc,iv3d_t)
    v2d(:,:,iv2d_sss) = v3d(:,:,ilev_sfc,iv3d_s)
  else
    call check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
    IF(istat /= NF90_NOERR) THEN
      WRITE(6,'(A)') 'netCDF OPEN ERROR in read_grd4 for ',sffile
      STOP(8)
    END IF
    WRITE(6,*) "read_grd4:: just opened file ", sffile

    !!! SST
    call check( NF90_INQ_VARID(ncid,'t_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,v2d(:,:,iv2d_sst)) )
    WHERE (kmt(:,:) .ge. 1) v2d(:,:,iv2d_sst) = v2d(:,:,iv2d_sst) - t0c ! kelvin to deg C

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SST"
      WRITE(6,*) "read_grd4:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sst) = ", MAXVAL(v2d(:,:,iv2d_sst))
    endif
! !STEVE: end

    !!! SSS
    call check( NF90_INQ_VARID(ncid,'s_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,v2d(:,:,iv2d_sss)) )
    !WRITE(6,*) "read_grd4:: just got data for variable temp"

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SSS"
      WRITE(6,*) "read_grd4:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sss) = ", MAXVAL(v2d(:,:,iv2d_sss))
    endif
! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the ALTIMETRY netcdf diagnostic file (SSH)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.false.) then
    !STEVE: for now, use vertically integrated heat content
    videpth = 300
    v2d(:,:,iv2d_ssh) = v3d(:,:,ilev_sfc,iv3d_t)*lev(ilev_sfc)
    do k=2,nlev
      if (lev(k) > videpth) EXIT
      v2d(:,:,iv2d_ssh) = v2d(:,:,iv2d_ssh) + v3d(:,:,k,iv3d_t)*(lev(k)-lev(k-1))
    enddo

    ! Divide out averages
    meanSSH = SUM(v2d(:,:,iv2d_ssh))/(nlon*nlat)
    where(v2d(:,:,iv2d_ssh) > 0) &
      v2d(:,:,iv2d_ssh) = v2d(:,:,iv2d_ssh)-meanSSH
  else
    call check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
    IF(istat /= NF90_NOERR) THEN
      WRITE(6,'(A)') 'netCDF OPEN ERROR in read_grd4 for ',sffile
      STOP(9)
    END IF
    WRITE(6,*) "read_grd4:: just read file ", sffile

    !!! SSH
    call check( NF90_INQ_VARID(ncid,'sea_lev',varid) )
    call check( NF90_GET_VAR(ncid,varid,v2d(:,:,iv2d_ssh)) )
    !WRITE(6,*) "read_grd4:: just got data for variable temp"

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SSH"
      WRITE(6,*) "read_grd4:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_ssh) = ", MAXVAL(v2d(:,:,iv2d_ssh))
    endif
! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  endif

  altimetry : if(DO_ALTIMETRY) then
    !STEVE: use the sea level perturbation from ocean_barotropic.res.nc
    call check( NF90_OPEN(bfile,NF90_NOWRITE,ncid) )
    WRITE(6,*) "read_grd:: just opened file ", bfile

    !!! SSH
    call check( NF90_INQ_VARID(ncid,'eta_t',varid) )
    call check( NF90_GET_VAR(ncid,varid,v2d(:,:,iv2d_eta)) )
    if (dodebug) WRITE(6,*) "read_grd:: just got data for variable eta_t"
    if (dodebug) WRITE(6,*) "read_grd:: finished processing data for variable SSH"

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,iv2d_eta) = v2d(:,:,iv2d_eta) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_grd:: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  endif altimetry


  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

! !STEVE: debug test
  if (.false.) then
    testfile = "test_read4.grd"
!   CALL write_bingrd4(trim(testfile),v3d,v2d)

    iunit=55
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=testfile,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

    WRITE(6,*) "Writing to", testfile
    irec=1
    DO n=1,nv3d
      DO k=1,nlev
        WRITE(6,*) "n, k, irec = ", n, k, irec
        WRITE(6,*) "max v3d(n) = ", MAXVAL(v3d(:,:,k,n))
        WRITE(iunit,REC=irec) v3d(:,:,k,n)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      WRITE(iunit,REC=irec) v2d(:,:,n)
      irec = irec + 1
    END DO
    CLOSE(iunit)

    WRITE(6,*) "Initially read from file: ", infile
    WRITE(6,*) "STOP 10"
    STOP(10)
  endif
! !STEVE: debug end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the surface fluxes netcdf file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !STEVE: debug test
  if (DO_SFCFLUX) then
    read( infile(3:4), '(i2)' )  islot
    MEM3 = infile(5:7)
    print *, " "
    print *, "======================================================="
    print *, "read_grd4 :: DO_SFCFLUX ==============================="
    print *, "infile = ", infile
    print *, "islot = ", islot
    sfc_infile = sfc_infiles(1)
    sfc_infile(5:7) = MEM3
    print *, "sfc_infile = ", sfc_infile

    ! Read coefficients file
    INQUIRE(FILE=coeff_s2mfile, EXIST=coeff_file_exists)
    if (coeff_file_exists) then
      print *, "read_grd4 :: call read_sfc_grid ==============================="
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat)
      print *, "read_grd4 :: call read_coeff ==============================="
      CALL read_coeff(coeff_s2mfile,xi,yi,xc,yc,mlon,mlat,nmlon,nmlat)
    else
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat)
      CALL coeff_sfc2model(slon,slat,lon,lat,xi,yi,xc,yc)
    endif

    ALLOCATE(sfc_data(nslon,nslat))
    ALLOCATE(model_data(nlon,nlat))

    ! Open and Read each of the SFC data files
    do i=1,nvsfc
      print *, "i,nvsfc = ", i,'/',nvsfc
      sfc_infile = sfc_infiles(i)
      sfc_infile(5:7) = MEM3
      print *, "sfc_infile = ", sfc_infile
      print *, "read_grd4 :: call read_sfc_data ==============================="
      CALL read_sfc_data(sfc_infile,sfc_names(i),islot+1,sfc_data) !STEVE: if using RIP, then use islot instead of islot+1

      ! Interpolate to the model grid using coefficients
      ! Convert sfc_data to model grid (on v2d)
      print *, "read_grd4 :: call itpl_sfc2model ==============================="
      CALL itpl_sfc2model(sfc_data,xi,yi,xc,yc,model_data)
      v2d(:,:,nv2d-nvsfc+i) = REAL(model_data,r_sngl)
    enddo

    print *, "read_grd4 :: DEALLOCATE", islot
    DEALLOCATE(sfc_data)
    DEALLOCATE(xi,yi,xc,yc,slon,slat)

!   STOP(2)
  endif

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  RETURN

END SUBROUTINE read_grd4

!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(outfile,v3d,v2d)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:,:) !(nlon,nlat,nlev)
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: mlon, mlat
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: model_data !(DO_SFCFLUX)
  REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
  CHARACTER(slen) :: sfc_outfile, sfc_infile
  CHARACTER(3) :: MEM3
  LOGICAL :: dodebug=.true.

  ALLOCATE(buf4(nlon,nlat,nlev))

  tsfile = trim(outfile)//'.ocean_temp_salt.res.nc'
  uvfile = trim(outfile)//'.ocean_velocity.res.nc'
  sffile = trim(outfile)//'.ocean_sbc.res.nc'
  bfile  = trim(outfile)//'.ocean_barotropic.res.nc'

  call check( NF90_OPEN(tsfile,NF90_WRITE,ncid) )
  WRITE(6,*) "write_grd:: just opened file ", tsfile

  !!! t
  call check( NF90_INQ_VARID(ncid,'temp',varid) )
  buf4=0.0
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. k) buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_t),r_sngl)
      END DO
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4) )

  !!! s
  call check( NF90_INQ_VARID(ncid,'salt',varid) )
  buf4=0.0
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. k) buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_s),r_sngl)
      END DO
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4) )

  call check( NF90_CLOSE(ncid) )

  call check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )
  WRITE(6,*) "write_grd:: just opened file ", uvfile

  !!! u
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  buf4=0.0
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. k) buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_u),r_sngl)
      END DO
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4) )

  !!! v
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  buf4=0.0
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. k) buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_v),r_sngl)
      END DO
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4) )

  call check( NF90_CLOSE(ncid) )

  call check( NF90_OPEN(sffile,NF90_WRITE,ncid) )
  WRITE(6,*) "write_grd:: just opened file ", sffile

  !!! SST
  call check( NF90_INQ_VARID(ncid,'t_surf',varid) )
  buf4=0.0
  DO j=1,nlat
    DO i=1,nlon
      if (kmt(i,j) .ge. 1) buf4(i,j,ilev_sfc) = REAL(v2d(i,j,iv2d_sst),r_sngl) + t0c ! deg C to kelvin
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )

  !!! SSS
  call check( NF90_INQ_VARID(ncid,'s_surf',varid) )
  buf4=0.0
  DO j=1,nlat
    DO i=1,nlon
      if (kmt(i,j) .ge. 1) buf4(i,j,ilev_sfc) = REAL(v2d(i,j,iv2d_sss),r_sngl)
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )

  ! Convert SLA stored in v2d back to SSH by adding back in model climatology
! if (DO_ALTIMETRY) then
!   v2d(:,:,iv2d_ssh) = v2d(:,:,iv2d_ssh) + SSHclm_m(:,:)
! endif

  !!! SSH
  call check( NF90_INQ_VARID(ncid,'sea_lev',varid) )
  buf4=0.0
  DO j=1,nlat
    DO i=1,nlon
      if (kmt(i,j) .ge. 1) buf4(i,j,ilev_sfc) = REAL(v2d(i,j,iv2d_ssh),r_sngl)
    END DO
  END DO
  call check( NF90_PUT_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the updated eta_t to analysis restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (DO_ALTIMETRY) then
    call check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    WRITE(6,*) "write_grd:: just opened file ", bfile

    call check( NF90_INQ_VARID(ncid,'eta_t',varid) )

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    buf4=0.0
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. 1) buf4(i,j,ilev_sfc) = REAL(v2d(i,j,iv2d_eta)+SSHclm_m(i,j),r_sngl)
      END DO
    END DO
    call check( NF90_PUT_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    call check( NF90_CLOSE(ncid) )
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the surface fluxes netcdf file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !STEVE: debug test
  if (DO_SFCFLUX) then
!   read( outfile(3:4), '(i2)' )  islot
    MEM3 = outfile(5:7)
    print *, " "
    print *, "======================================================="
    print *, "write_grd :: DO_SFCFLUX ==============================="
    print *, "outfile = ", outfile
    print *, "islot = ", islot
    !STEVE: NEED TO WRITE OUT FOR EACH MEMBER!, (not _000_)

    ! Read coefficients file
    INQUIRE(FILE=coeff_m2sfile, EXIST=coeff_file_exists)
    if (coeff_file_exists) then
      print *, "write_grd :: call read_coeff ==============================="
      CALL read_coeff(coeff_m2sfile,xi,yi,xc,yc,slon,slat,nslon,nslat)
    else
      sfc_infile = sfc_infiles(1)
      sfc_infile(5:7) = MEM3
      print *, "sfc_infile = ", sfc_infile
      print *, "write_grd :: call read_sfc_grid ==============================="
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat)
      print *, "write_grd :: call coeff_model2sfc ==============================="
      CALL coeff_model2sfc(lon,lat,slon,slat,xi,yi,xc,yc)
    endif
    ALLOCATE(model_data(nlon,nlat))
    ALLOCATE(sfc_data(nslon,nslat))
    ALLOCATE(data4D(nslon,nslat,1))

    ! Open and write each of the SFC data files
    do i=1,nvsfc
      print *, "i,nvsfc = ", i,'/',nvsfc
      sfc_outfile = sfc_outfiles(i)
      sfc_outfile(5:7) = MEM3
      print *, "sfc_outfile = ", sfc_outfile

      ! Interpolate to the model grid using coefficients
      ! Convert sfc_data to model grid (on v2d)
      print *, "write_grd :: call itpl_model2sfc ==============================="
      CALL itpl_model2sfc(v2d(:,:,nv2d-nvsfc+i),xi,yi,xc,yc,sfc_data)
      data4D(:,:,1) = sfc_data
      print *, "write_grd :: call write_sfc_data ==============================="
      CALL write_sfc_data(sfc_outfile,sfc_names(i),0,data4D,slon,slat,nslon,nslat) !STEVE: if using RIP, then use islot instead of islot+1
    enddo

    print *, "write_grd :: DEALLOCATE", islot
    DEALLOCATE(sfc_data)
    DEALLOCATE(xi,yi,xc,yc,slon,slat)

!   STOP(3)
  endif

  DEALLOCATE(buf4)

  RETURN
END SUBROUTINE write_grd

subroutine check(status)
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

SUBROUTINE write_grd4(outfile,v3d_in,v2d_in)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2d(:,:,:) !(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: t3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY)
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  LOGICAL, PARAMETER :: do_physlimit=.true.
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: mlon, mlat
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: model_data !(DO_SFCFLUX)
  REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
  CHARACTER(slen) :: sfc_outfile, sfc_infile
  CHARACTER(3) :: MEM3


  !STEVE: this is the routine that writes out the individual analysis files for
  !       each esnsemble member in netcdf format.

  tsfile = trim(outfile)//'.ocean_temp_salt.res.nc'
  uvfile = trim(outfile)//'.ocean_velocity.res.nc'
  sffile = trim(outfile)//'.ocean_sbc.res.nc'
  bfile  = trim(outfile)//'.ocean_barotropic.res.nc'

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  v3d = v3d_in
  v2d = v2d_in

  ! STEVE: for safety, clean up the variables for output:
  if (do_physlimit) then
  do k=1,nlev
    do j=1,nlat
      do i=1,nlon
!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_t) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
!         v3d(i,j,k,iv3d_t) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_s) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
!         v3d(i,j,k,iv3d_s) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_u) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
!         v3d(i,j,k,iv3d_u) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_v) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
!         v3d(i,j,k,iv3d_v) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (k .eq. 1 .and. kmt(i,j) .eq. 0) then 
!         if (v2d(i,j,iv2d_sst) .ne. 0.0 ) v2d(i,j,iv2d_sst) = 0.0 !NF90_FILL_FLOAT
!         if (v2d(i,j,iv2d_sss) .ne. 0.0 ) v2d(i,j,iv2d_sss) = 0.0 !NF90_FILL_FLOAT
!         if (v2d(i,j,iv2d_ssh) .ne. 0.0 ) v2d(i,j,iv2d_ssh) = 0.0 !NF90_FILL_FLOAT
!       endif

        if (v3d(i,j,k,iv3d_t) < -4) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = -4.0
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sst) < -4) v2d(i,j,iv2d_sst) = -4.0

        if (v3d(i,j,k,iv3d_t) > 40.0) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = 40.0
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sst) > 40.0) v2d(i,j,iv2d_sst) = 40.0

        if (v3d(i,j,k,iv3d_s) < 0 ) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = 0.0
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sss) < 0) v2d(i,j,iv2d_sss) = 0.0

        if (v3d(i,j,k,iv3d_s) > 50.0) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = 50.0
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sss) > 50.0) v2d(i,j,iv2d_sss) = 50.0

      enddo
    enddo
  enddo
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !STEVE: open temp/salt file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(tsfile,NF90_WRITE,ncid) )

  !!! t
  !STEVE: for debugging
  if (.false.) then
  do m=1,nv3d
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
!       if ( isnan( REAL(v3d(i,j,k,m),r_size) ) )then
!         WRITE(6,*) "common_mom4.f90::write_grd4:: ERROR: found NaN..."
!         WRITE(6,*) "v3d(i,j,k,m) contains NaN. i,j,k,m = ", i,j,k,m
!         STOP 1
!       endif
        enddo
      enddo
    enddo
  enddo
  endif

  call check( NF90_INQ_VARID(ncid,'temp',varid) )

  !STEVE: debug DEBUG
  !STEVE: switch out the data to see if this writes properly
  !ALLOCATE(t3d(nlon,nlat,nlev,nv3d))
  !t3d(:,:,:,iv3d_t) = 1.0
  !call check( NF90_PUT_VAR(ncid,varid,t3d(:,:,:,iv3d_t)) )
  !DEALLOCATE(t3d)
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_t)) )

  !!! s
  call check( NF90_INQ_VARID(ncid,'salt',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_s)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! uv file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )

  !!! u
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_u)) )

  !!! v
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_v)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! sfc file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(sffile,NF90_WRITE,ncid) )
  
  !!! SST
  call check( NF90_INQ_VARID(ncid,'t_surf',varid) )
  WHERE (kmt(:,:) .ge. 1) v2d(:,:,iv2d_sst) = v2d(:,:,iv2d_sst) + t0c ! kelvin
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d_sst)) )

  !!! SSS
  call check( NF90_INQ_VARID(ncid,'s_surf',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d_sss)) )

  !!! SSH
  call check( NF90_INQ_VARID(ncid,'sea_lev',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d_ssh)) )
  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the updated eta_t to analysis restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (DO_ALTIMETRY) then
    call check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    call check( NF90_INQ_VARID(ncid,'eta_t',varid) )

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,iv2d_eta) = v2d(:,:,iv2d_eta) + SSHclm_m(:,:)

    call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d_eta)) )
    call check( NF90_CLOSE(ncid) )
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the surface fluxes netcdf file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !STEVE: debug test
! !STEVE: the surface file output is only 1 timestamp, unlike one input for each timeslot
! !STEVE: Later, perhaps output the analysis weights Wa at the surface and apply
! to the next timestep, applying timefiltering to the weights so they are
! gradually applied as a spread-correction over time, and the analysis increment
! with time smoothing is applied to gradually bias-correct over time.
  if (DO_SFCFLUX) then
!   read( outfile(3:4), '(i2)' )  islot
    MEM3 = outfile(5:7)
    print *, " "
    print *, "======================================================="
    print *, "write_grd4 :: DO_SFCFLUX ==============================="
    print *, "outfile = ", outfile
!   print *, "islot = ", islot
    !STEVE: NEED TO WRITE OUT FOR EACH MEMBER!, (not _000_)
    !(should be fixed here, fix in write_grd also - 7/18/13)

    ! Read coefficients file
    INQUIRE(FILE=coeff_m2sfile, EXIST=coeff_file_exists)
    if (coeff_file_exists) then
      print *, "write_grd4 :: call read_coeff ==============================="
      CALL read_coeff(coeff_m2sfile,xi,yi,xc,yc,slon,slat,nslon,nslat)
    else
      sfc_infile = sfc_infiles(1)
      sfc_infile(5:7) = MEM3
      print *, "sfc_infile = ", sfc_infile
      print *, "write_grd4 :: call read_sfc_grid ==============================="
      CALL read_sfc_grid(sfc_infile,slon,slat,nslon,nslat)
      print *, "write_grd4 :: call coeff_model2sfc ==============================="
      CALL coeff_model2sfc(lon,lat,slon,slat,xi,yi,xc,yc)
    endif
    ALLOCATE(model_data(nlon,nlat))
    ALLOCATE(sfc_data(nslon,nslat))
    ALLOCATE(data4D(nslon,nslat,1))

    ! Open and write each of the SFC data files
    do i=1,nvsfc
      print *, "i,nvsfc = ", i,'/',nvsfc
      sfc_outfile = sfc_outfiles(i)
      sfc_outfile(5:7) = MEM3
      print *, "sfc_outfile = ", sfc_outfile

      ! Interpolate to the model grid using coefficients
      ! Convert sfc_data to model grid (on v2d)
      model_data = REAL(v2d(:,:,nv2d-nvsfc+i),r_size)
      print *, "write_grd4 :: call itpl_model2sfc ==============================="
      CALL itpl_model2sfc(model_data,xi,yi,xc,yc,sfc_data)

      print *, "write_grd4 :: call write_sfc_data ==============================="
      data4D(:,:,1) = sfc_data
      CALL write_sfc_data(sfc_outfile,sfc_names(i),0,data4D,slon,slat,nslon,nslat) !STEVE: if using RIP, then use islot instead of islot+1
    enddo

    print *, "write_grd4 :: DEALLOCATE" !, islot
    DEALLOCATE(sfc_data)
    DEALLOCATE(data4D)
    DEALLOCATE(xi,yi,xc,yc,slon,slat)

!   STOP(4)
  endif

  RETURN
END SUBROUTINE write_grd4

! Write a MOM4 XXX_increment.nc file for IAU
SUBROUTINE write_inc4(type,v3d,v2d)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: type
  CHARACTER(128) :: file
  ! These should contain (Analysis - Background)
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  CHARACTER(slen) :: outfile
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  INTEGER :: iv3d=0, iv2d=0
  ! For netcdf create:
  integer :: err,dimidx,dimidy,dimidz,dimidt,dd(4),varidx,varidy,varidz,varidt,myvar
  integer :: startA(1),start3D(4),countA(1),count3D(4),start2D(3),count2D(3)

  SELECT CASE(trim(type))
  CASE('temp')
    file = 'temp'
    iv3d = iv3d_t
  CASE('salt')
    file = 'salt'
    iv3d = iv3d_s
  CASE('u')
    file = 'u'
    iv3d = iv3d_u
  CASE('v')
    file = 'v'
    iv3d = iv3d_v
  CASE('ssh')
    ! STEVE: need this since mom4p1/Bluelink already coded it as 'eta'
    file = 'eta'
    iv2d = iv2d_ssh
  CASE('eta')
    file = 'eta'
    iv2d = iv2d_ssh
  CASE DEFAULT
    WRITE(6,*) "write_inc4:: unsupported type: ", type
    stop 1
  END SELECT

  outfile = trim(file)//'_increment.nc'

  !call check( NF90_OPEN(outfile,NF90_WRITE,ncid) )
  call check( NF90_CREATE(outfile,NF90_CLOBBER,ncid) )

  ! For NF90_CREATE:
  ! Define the dimensions
  call check( NF90_DEF_DIM(ncid,'xaxis_1',nlon,dimidx) )
  call check( NF90_DEF_DIM(ncid,'yaxis_1',nlat,dimidy) )
  call check( NF90_DEF_DIM(ncid,'zaxis_1',nlev,dimidz) )
  call check( NF90_DEF_DIM(ncid,'Time',NF90_UNLIMITED,dimidt) )

  !
  ! Define the variables (include 2 variables giving the dim. values)
  !
  call check( NF90_DEF_VAR(ncid,'xaxis_1',NF90_REAL,dimidx,varidx) )
  call check( NF90_DEF_VAR(ncid,'yaxis_1',NF90_REAL,dimidy,varidy) )
  call check( NF90_DEF_VAR(ncid,'zaxis_1',NF90_REAL,dimidz,varidz) )
  call check( NF90_DEF_VAR(ncid,'Time',NF90_REAL,dimidt,varidt) )

  dd(1)=dimidx
  dd(2)=dimidy
  dd(3)=dimidz
  dd(4)=dimidt
  call check( NF90_DEF_VAR(ncid,trim(file),NF90_REAL,dd(1:4),myvar) )

  !
  ! Define the Attributes:
  !

  ! x dimension attributes
  call check( NF90_PUT_ATT(ncid, varidx, 'units', 'degrees_east') )
  call check( NF90_PUT_ATT(ncid, varidx, 'cartesian_axis', 'X') )

  ! y dimension attributes
  call check( NF90_PUT_ATT(ncid, varidy, 'units', 'degrees_north') )
  call check( NF90_PUT_ATT(ncid, varidy, 'cartesian_axis', 'Y') )

  if (iv3d > 0) then
    ! z dimension attributes
    call check( NF90_PUT_ATT(ncid, varidz, 'units', 'meters') )
    call check( NF90_PUT_ATT(ncid, varidz, 'cartesian_axis', 'Z') )
  endif

  ! Time dimension attributes (required)
  ! TIME:cartesian_axis = "T" through ncatted.
  call check( NF90_PUT_ATT(ncid, varidt, 'cartesian_axis', 'T') )
  call check( NF90_PUT_ATT(ncid, varidt, 'units', 'seconds since 1900-01-01 00:00:00 UTC') )
! call check( NF90_PUT_ATT(ncid, varidt, 'calendar', 'GREGORIAN') )
  call check( NF90_PUT_ATT(ncid, varidt, 'calendar', 'julian') )

! print*,'x dim ID ',dimidx
! print*,'y dim ID ',dimidy
! print*,'z dim ID ',dimidz
! print*,'t dim ID ',dimidt
! print*,'x var ID ',varidx
! print*,'y var ID ',varidy
! print*,'z var ID ',varidz
! print*,'t var ID ',varidt
! print*,'main ID  ',myvar

  ! Change mode of netCDF operation
  call check( NF90_ENDDEF(ncid) )

  ! Output the values of the variables (include dimension variables)
  ! x dimension values
  call check( NF90_PUT_VAR(ncid,varidx,lon) )

  ! y dimension values
  call check( NF90_PUT_VAR(ncid,varidy,lat) )

  if (iv3d > 0) then
    ! z dimension values
    call check( NF90_PUT_VAR(ncid,varidz,lev) )
  endif

! ! t dimension values
! startA(1)=1
! countA(1)=1
! call check( NF90_PUT_VAR(ncid,varidt,1) )

  !!! t
  !STEVE: for debugging
  if (.false. .and. iv3d > 0) then
  do m=1,nv3d
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
        !if ( isnan( REAL(v3d(i,j,k,m),r_size) ) )then
        !  WRITE(6,*) "common_mom4.f90::write_inc4:: ERROR: found NaN..."
        !  WRITE(6,*) "v3d(i,j,k,m) contains NaN. i,j,k,m = ", i,j,k,m
        !  STOP(1)
        !endif
          if ( v3d(i,j,k,m) > 1 ) then
            WRITE(6,*) "i,j,k,m, v3d = ", i,j,k,m, v3d(i,j,k,m)
          endif
        enddo
      enddo
    enddo
  enddo
  endif

  call check( NF90_INQ_VARID(ncid,trim(file),varid) )
  IF(istat /= NF90_NOERR) THEN
    !STEVE: debugging
    WRITE(6,*) "common_mom4.f90::write_inc4:: ERROR: NF90_INQ_VARID failed"
    stop 1
  END IF

  ! Either output 3d or 2d data
  if (iv3d > 0) then 
    !call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d))
    istat=NF90_PUT_VAR(ncid,myvar,v3d(:,:,:,iv3d)) !start3D,count3D,
  elseif (iv2d > 0) then
    !call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d))
    istat=NF90_PUT_VAR(ncid,myvar,v2d(:,:,iv2d)) !start2D,count2D,
  endif

  IF(istat /= NF90_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (write_inc4:: temp)'
    WRITE(6,*) "call check( ", istat
    WRITE(6,*) "ncid = ", ncid
    WRITE(6,*) "varid = ", varid
    if (iv3d > 0) then 
      WRITE(6,*) "iv3d = ", iv3d
      WRITE(6,*) "MAXVAL(v3d(:,:,:,iv3d)) = ", MAXVAL(v3d(:,:,:,iv3d))
      WRITE(6,*) "MINVAL(v3d(:,:,:,iv3d)) = ", MINVAL(v3d(:,:,:,iv3d))
      WRITE(6,*) "v3d(:,:,:,iv3d) = ", v3d(:,:,:,iv3d)
    elseif (iv2d > 0) then
      WRITE(6,*) "iv2d = ", iv2d
      WRITE(6,*) "MAXVAL(v2d(:,:,iv2d)) = ", MAXVAL(v2d(:,:,iv2d))
      WRITE(6,*) "MINVAL(v2d(:,:,iv2d)) = ", MINVAL(v2d(:,:,iv2d))
      WRITE(6,*) "v2d(:,:,iv2d) = ", v2d(:,:,iv2d)
    endif
    STOP(11)
  END IF

  ! Close the file
  call check( NF90_CLOSE(ncid) )

  RETURN
END SUBROUTINE write_inc4

! Write a MOM4 XXX_increment.nc file for IAU
SUBROUTINE write_inc4_sf(type,v3d,v2d)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: type
  CHARACTER(128) :: file
  ! These should contain (Analysis - Background)
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  CHARACTER(slen) :: outfile
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  INTEGER :: iv3d=0, iv2d=0
  ! For netcdf create:
  integer :: err,dimidx,dimidy,dimidz,dimidt,dd(4),varidx,varidy,varidz,varidt,myvar
  integer :: startA(1),start3D(4),countA(1),count3D(4),start2D(3),count2D(3)

  SELECT CASE(trim(type))
  CASE('temp')
    file = 'temp'
    iv3d = iv3d_t
  CASE('salt')
    file = 'salt'
    iv3d = iv3d_s
  CASE('u')
    file = 'u'
    iv3d = iv3d_u
  CASE('v')
    file = 'v'
    iv3d = iv3d_v
  CASE('ssh')
    ! STEVE: need this since mom4p1/Bluelink already coded it as 'eta'
    file = 'eta'
    iv2d = iv2d_ssh
  CASE('eta')
    file = 'eta'
    iv2d = iv2d_ssh
  CASE DEFAULT
    WRITE(6,*) "write_inc4:: unsupported type: ", type
    stop 1
  END SELECT

  outfile = trim(file)//'_increment.nc'

  !call check( NF90_OPEN(outfile,NF90_WRITE,ncid) )
  call check( NF90_CREATE(outfile,NF90_CLOBBER,ncid) )

  ! For NF90_CREATE:
  ! Define the dimensions
  call check( NF90_DEF_DIM(ncid,'GRID_X_T',nlon,dimidx) )
  call check( NF90_DEF_DIM(ncid,'GRID_Y_T',nlat,dimidy) )
  call check( NF90_DEF_DIM(ncid,'ZT',nlev,dimidz) )
  call check( NF90_DEF_DIM(ncid,'TIME',NF90_UNLIMITED,dimidt) )

  !
  ! Define the variables (include 2 variables giving the dim. values)
  !
  call check( NF90_DEF_VAR(ncid,'GRID_X_T',NF90_REAL,dimidx,varidx) )
  call check( NF90_DEF_VAR(ncid,'GRID_Y_T',NF90_REAL,dimidy,varidy) )
  call check( NF90_DEF_VAR(ncid,'ZT',NF90_REAL,dimidz,varidz) )
  call check( NF90_DEF_VAR(ncid,'TIME',NF90_REAL,dimidt,varidt) )

  dd(1)=dimidx
  dd(2)=dimidy
  dd(3)=dimidz
  dd(4)=dimidt
  call check( NF90_DEF_VAR(ncid,trim(file),NF90_REAL,dd,myvar) )

  !
  ! Define the Attributes:
  !

  ! x dimension attributes
  call check( NF90_PUT_ATT(ncid, varidx, 'units', 'degrees_east') )
  call check( NF90_PUT_ATT(ncid, varidx, 'cartesian_axis', 'X') )

  ! y dimension attributes
  call check( NF90_PUT_ATT(ncid, varidy, 'units', 'degrees_north') )
  call check( NF90_PUT_ATT(ncid, varidy, 'cartesian_axis', 'Y') )

  if (iv3d > 0) then
    ! z dimension attributes
    call check( NF90_PUT_ATT(ncid, varidz, 'units', 'meters') )
    call check( NF90_PUT_ATT(ncid, varidz, 'cartesian_axis', 'Z') )
  endif

  ! Time dimension attributes (required)
  ! TIME:cartesian_axis = "T" through ncatted.
  call check( NF90_PUT_ATT(ncid, varidt, 'cartesian_axis', 'T') )
  call check( NF90_PUT_ATT(ncid, varidt, 'units', 'seconds since 1900-01-01 00:00:00 UTC') )
  call check( NF90_PUT_ATT(ncid, varidt, 'calendar_type', 'GREGORIAN') )
  call check( NF90_PUT_ATT(ncid, varidt, 'time_origin', '01-JAN-1900 00:00:00') )
  call check( NF90_PUT_ATT(ncid, varidt, 'modulo', ' ') )

! print*,'x dim ID ',dimidx
! print*,'y dim ID ',dimidy
! print*,'z dim ID ',dimidz
! print*,'t dim ID ',dimidt
! print*,'x var ID ',varidx
! print*,'y var ID ',varidy
! print*,'z var ID ',varidz
! print*,'t var ID ',varidt
! print*,'main ID  ',myvar

  ! Change mode of netCDF operation
  call check( NF90_ENDDEF(ncid) )

  ! Output the values of the variables (include dimension variables)
  ! x dimension values
  call check( NF90_PUT_VAR(ncid,varidx,lon) )

  ! y dimension values
  call check( NF90_PUT_VAR(ncid,varidy,lat) )

  if (iv3d > 0) then
    ! z dimension values
    call check( NF90_PUT_VAR(ncid,varidz,lev) )
  endif

! ! t dimension values
! call check( NF90_PUT_VAR(ncid,varidt,1) )

  !!! t
  !STEVE: for debugging
  if (.false. .and. iv3d > 0) then
  do m=1,nv3d
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
        !if ( isnan( REAL(v3d(i,j,k,m),r_size) ) )then
        !  WRITE(6,*) "common_mom4.f90::write_inc4:: ERROR: found NaN..."
        !  WRITE(6,*) "v3d(i,j,k,m) contains NaN. i,j,k,m = ", i,j,k,m
        !  STOP(1)
        !endif
          if ( v3d(i,j,k,m) > 1 ) then
            WRITE(6,*) "i,j,k,m, v3d = ", i,j,k,m, v3d(i,j,k,m)
          endif
        enddo
      enddo
    enddo
  enddo
  endif

  call check( NF90_INQ_VARID(ncid,trim(file),varid) )

  ! Either output 3d or 2d data
  if (iv3d > 0) then 
    !call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d)) )
    istat=NF90_PUT_VAR(ncid,myvar,v3d(:,:,:,iv3d))
  elseif (iv2d > 0) then
    !call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d)) )
    istat=NF90_PUT_VAR(ncid,myvar,v2d(:,:,iv2d))
  endif

  IF(istat /= NF90_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (write_inc4:: temp)'
    WRITE(6,*) "call check( ", istat
    WRITE(6,*) "ncid = ", ncid
    WRITE(6,*) "varid = ", varid
    if (iv3d > 0) then 
      WRITE(6,*) "iv3d = ", iv3d
      WRITE(6,*) "MAXVAL(v3d(:,:,:,iv3d)) = ", MAXVAL(v3d(:,:,:,iv3d))
      WRITE(6,*) "MINVAL(v3d(:,:,:,iv3d)) = ", MINVAL(v3d(:,:,:,iv3d))
      WRITE(6,*) "v3d(:,:,:,iv3d) = ", v3d(:,:,:,iv3d)
    elseif (iv2d > 0) then
      WRITE(6,*) "iv2d = ", iv2d
      WRITE(6,*) "MAXVAL(v2d(:,:,iv2d)) = ", MAXVAL(v2d(:,:,iv2d))
      WRITE(6,*) "MINVAL(v2d(:,:,iv2d)) = ", MINVAL(v2d(:,:,iv2d))
      WRITE(6,*) "v2d(:,:,iv2d) = ", v2d(:,:,iv2d)
    endif
    STOP(12)
  END IF

  ! Close the file
  call check( NF90_CLOSE(ncid) )

  RETURN
END SUBROUTINE write_inc4_sf


!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_bingrd(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:) !(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  ALLOCATE(buf4(nlon,nlat))

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      READ(iunit,REC=irec) buf4
      irec = irec + 1
      v3d(:,:,k,n) = REAL(buf4,r_size)
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  RETURN
END SUBROUTINE read_bingrd

SUBROUTINE read_bingrd4(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  RETURN
END SUBROUTINE read_bingrd4

!-- Write a grid file -------------------------------------------------
SUBROUTINE write_bingrd(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:) !(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec
  LOGICAL :: dodebug=.false.

  ALLOCATE(buf4(nlon,nlat))

  if (dodebug) print *, "write_bingrd:: open filename = ",filename
  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_bingrd:: nij0,iolength = ", nij0,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      buf4 = 0.0
      buf4 = REAL(v3d(:,:,k,n),r_sngl)
      if (dodebug) print *, "write_bingrd:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) buf4
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    buf4 = 0.0
    buf4 = REAL(v2d(:,:,n),r_sngl)
    if (dodebug) print *, "write_bingrd:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) buf4
    irec = irec + 1
  END DO

  CLOSE(iunit)

  DEALLOCATE(buf4)

  RETURN
END SUBROUTINE write_bingrd

SUBROUTINE write_bingrd4(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL :: dodebug=.false.

  if (dodebug) print *, "write_bingrd4:: open filename = ",filename
  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_bingrd4:: nij0,iolength = ", nij0,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      if (dodebug) print *, "write_bingrd4:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    if (dodebug) print *, "write_bingrd4:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_bingrd4
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

! ALLOCATE(v3dm(nij,nlev,nv3d),v2dm(nij,nv2d))

  DO n=1,nv3d
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STEVE: additional subroutines for OCEAN
!STEVE: all of these are still direct copies form mom2 version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE minkowski_flick(v3d0,v2d0,v3d,v2d)
  ! The new output v3d has been grown around the ocean perimeter by one grid point.
  ! This is a simple version of the general minkowski sum.
  REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d), INTENT(IN) :: v3d0
  REAL(r_size), DIMENSION(nlon,nlat,nv2d), INTENT(IN) :: v2d0
  REAL(r_size), INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size), INTENT(OUT) :: v2d(nlon,nlat,nv2d) !STEVE: not updating v2d for now
  REAL(r_size), ALLOCATABLE :: buf(:,:) !(nlon+2,nlat+2)
  REAL(r_size), ALLOCATABLE :: buf3d(:,:,:) !(nlon+2,nlat+2,nlev+2)
  REAL(r_size), ALLOCATABLE :: buf2d(:,:,:) !(nlon+2,nlat+2,3)
  REAL(r_size), ALLOCATABLE :: vcnt(:,:,:)   !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE :: mask(:,:,:,:) !(nlon,nlat,nlev,9) -> stores boundary location
  INTEGER, ALLOCATABLE :: kmt3d(:,:,:) !(nlon,nlat,nlev)
  INTEGER, ALLOCATABLE :: kmt2d(:,:) !(nlon,nlat)
  INTEGER :: nkmterr = 0
  INTEGER :: i,j,k,n
! LOGICAL, PARAMETER :: dodebug=.true.
  LOGICAL, PARAMETER :: do3d=.true.

  if (do3d) then
  WRITE(6,*) "In minkowski flick."

  ALLOCATE(buf(nlon+2,nlat+2))
  ALLOCATE(buf3d(nlon+2,nlat+2,nlev+2))
  ALLOCATE(buf2d(nlon+2,nlat+2,3))
  ALLOCATE(vcnt(nlon,nlat,nlev))
  ALLOCATE(mask(nlon,nlat,nlev,9))
  ALLOCATE(kmt3d(nlon,nlat,nlev))
  ALLOCATE(kmt2d(nlon,nlat))

  ! Create a 3D land sea map from kmt
  kmt3d=0
  kmt2d=0
  do k=1,nlev
    WHERE(k <= kmt) kmt3d(:,:,k) = 1
  enddo
  WHERE(kmt > 0) kmt2d = 1

  ! Do (3D)
  ! For each variable
  v3d=0
  WRITE(6,*) "Doing 3D part of minkowski flick."
  do n=1,nv3d ! Mostly needed for temperature and salinity (3,4)
    WRITE(6,*) "Doing variable: ", n
    mask=0
    vcnt=0
    ! Reset buffer
    buf3d=0

    ! flutter grid up, down left, right and diagonal.
    ! (Use v3d and v2d to store the boundary data during computation, then add it to v3d0 and v2d0 to get new grd)

    !STEVE: stick the data in the middle of the buffer
    buf3d(2:nlon+1,2:nlat+1,2:nlev+1) = v3d0(:,:,:,n)

    ! Follow keypad order for: 1..9
    do i=0,2
      do j=0,2
        do k=0,2
          where(kmt3d .lt. 1 .and. buf3d(1+i:nlon+i,1+j:nlat+j,1+k:nlev+k) > 0.0 )
            v3d(:,:,:,n) = v3d(:,:,:,n) + buf3d(1+i:nlon+i,1+j:nlat+j,1+k:nlev+k) ! sum all the values on this gridpoint that have water in them
            vcnt = vcnt + 1
          end where
        enddo
      enddo
    enddo

    ! If it intersects land (kmt<lev), then it's a boundary point (on the land side).
    ! Average the flutter values to get an approximate extrapolation value.
    where(vcnt > 0.0) ! (be careful not to divide by zero...)
      v3d(:,:,:,n) = v3d(:,:,:,n) / vcnt
    end where
  enddo

! Add back on to the pre-existing values
  v3d = v3d0 + v3d

  ! Do (2D)
  ! For each variable
  WRITE(6,*) "Doing 2D part of minkowski flick."
  v2d=0
  do n=1,nv2d
    mask=0
    vcnt=0
    ! Reset buffer
    buf2d=0

    ! flutter grid up, down left, right and diagonal.
    ! (Use v3d and v2d to store the boundary data during computation, then add it to v3d0 and v2d0 to get new grd)

    !STEVE: stick the data in the middle of the buffer
    buf2d(2:nlon+1,2:nlat+1,2) = v2d0(:,:,n)

    ! Follow keypad order for: 1..9
    do i=0,2
      do j=0,2
        do k=0,2
          where(kmt2d .lt. 1 .and. buf2d(1+i:nlon+i,1+j:nlat+j,1+k) > 0.0 )
            v2d(:,:,n) = v2d(:,:,n) + buf2d(1+i:nlon+i,1+j:nlat+j,1+k) ! sum all the values on this gridpoint that have water in them
            vcnt(:,:,ilev_sfc) = vcnt(:,:,ilev_sfc) + 1
          end where
        enddo
      enddo
    enddo

    ! If it intersects land (kmt<lev), then it's a boundary point (on the land side).
    ! Average the flutter values to get an approximate extrapolation value.
    where(vcnt(:,:,ilev_sfc) > 0.0) ! (be careful not to divide by zero...)
      v2d(:,:,n) = v2d(:,:,n) / vcnt(:,:,ilev_sfc)
    end where
  enddo

! Add back on to the pre-existing values
  v2d = v2d0 + v2d

  !Write to test
  if (dodebug) then
    CALL write_bingrd('test_mink_v3d.grd',v3d,v2d)
    CALL write_bingrd('test_mink_v3d0.grd',v3d0,v2d0)
  endif

  DEALLOCATE(buf)
  DEALLOCATE(buf3d)
  DEALLOCATE(buf2d)
  DEALLOCATE(vcnt)
  DEALLOCATE(mask)
  DEALLOCATE(kmt3d)
  DEALLOCATE(kmt2d)

  return
  endif
END SUBROUTINE minkowski_flick

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE grd_to_cor(v3d1,v2d1,v3d2,v2d2,days)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: days !STEVE: the number of days the model experiment will run for
  REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d), INTENT(IN) :: v3d1,v3d2 !analysis=1,guess=2
  REAL(r_size), DIMENSION(nlon,nlat,nv2d), INTENT(IN) :: v2d1,v2d2
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:,:) :: v3d, v3d3, blank3d1,blank3d2,bdypt3d
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: v2d, v2d3, blank2d1,blank2d2,bdypt2d
  INTEGER, ALLOCATABLE :: bdypt(:,:,:) !(nlon,nlat,nlev) !(used to identify ocean points that are adjacent to land)
  REAL(r_size) :: deltat !STEVE: divisor for model time steps
  REAL :: ainc_thresh
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  INTEGER :: nkmterr=0
! LOGICAL :: dodebug = .true.
  INTEGER, DIMENSION(3) :: errloc
  ! For analysis increment quality control:
  REAL(r_size) :: ainc_max_temp, ainc_min_temp, ainc_max_salt, ainc_min_salt


  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
  if (days .eq. 0) then
    WRITE(6,*) "Generating zero correctors file..."

    ! Write out the appropriate correctors for mom4
    v3d = 0.0
    v2d = 0.0
    CALL write_inc4('temp',v3d,v2d)
    CALL write_inc4('salt',v3d,v2d)
    CALL write_inc4('u',v3d,v2d)
    CALL write_inc4('v',v3d,v2d)
    CALL write_inc4('ssh',v3d,v2d)
    return
  endif

  if (kmt(1,1) < 0) then
    WRITE(6,*) "common_mom4.f90::grd_to_cor::"
    WRITE(6,*) "Error: the kmt.dta file has not been read into variable kmt"
    stop 1
  endif

  !STEVE: check background against kmt file
  if (.false. .and. dodebug) then
    nkmterr = 0
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
          !STEVE: match this with the update_kmt subroutine above so that the error is corrected within letkf.
          !if ( (kmt(i,j) >= k .and. v3d2(i,j,k,n) .eq. 0.0) .or. (kmt(i,j) < 1 .and. v3d2(i,j,k,n) > 0.0) ) then
          if ( kmt(i,j) < 1 .and. v3d2(i,j,k,n) > 0.000001 ) then
            nkmterr = nkmterr+1
            WRITE(6,*) "i,j,k,n = ", i,j,k,n
            WRITE(6,*) "kmt(i,j) = ", kmt(i,j)
            WRITE(6,*) "v3d2(i,j,k,:) = ", v3d2(i,j,k,:)
            WRITE(6,*) "MAXVAL(v3d2(i,j,k,:)) = ", MAXVAL(v3d2(i,j,k,:))
            WRITE(6,*) "-------------------------------------"
            !STEVE: add extra checks to see if it is really a coastal value
          endif
        enddo
      enddo
    enddo
    if (nkmterr > 0) then
      WRITE(6,*) "grd_to_cor:: nkmterr = ", nkmterr
      WRITE(6,*) "ERROR: STOP..."
      STOP(13)
    endif
  endif

  ! Get analysis increment:
  ! ADJUST OUTPUT TEMP AND SALT for model timestep:
  ! assuming model timestep of 1 hour:
  deltat = 1.0d0 !( 60.0d0 * 60.0d0 * 24.0d0 * REAL(days,r_size) )
  v3d = 0.0
  v2d = 0.0
  do n=1,nv3d
   !WRITE(6,*) "n = ",n
   do k=1,nlev
    WHERE (kmt >= k) v3d(:,:,k,n) = (v3d1(:,:,k,n) - v3d2(:,:,k,n)) / deltat
    WHERE (kmt < 1) v3d(:,:,k,n) = 0.0 !STEVE: shouldn't be necessary
    !WRITE(6,*) "k = ", k
    !WRITE(6,*) "maxval = ", maxval(v3d(:,:,k,n)) !*deltat
    !WRITE(6,*) "minval = ", minval(v3d(:,:,k,n)) !*deltat
   enddo
  enddo
  do n=1,nv2d
    WHERE (kmt > 0) v2d(:,:,n) = (v2d1(:,:,n) - v2d2(:,:,n)) / deltat
    WHERE (kmt < 1) v2d(:,:,n) = 0.0 !STEVE: shouldn't be necessary
  enddo

  !check for large corrections on the borders with land. This is an artifact of letkf with land constraints (investigate further)
  if (.true.) then ! Identify boundary points
   ALLOCATE(bdypt(nlon,nlat,nlev))
    bdypt = 0
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
          if (kmt(i,j) .ge. k) then
            if (i > 1 .and. j > 1 .and. kmt(i-1,j-1) .lt. k ) bdypt(i,j,k) = 1
            if (i > 1 .and. kmt(i-1,j) .lt. k ) bdypt(i,j,k) = 1
            if (j > 1 .and. kmt(i,j-1) .lt. k ) bdypt(i,j,k) = 1
            if (i < nlon .and. j < 1 .and. kmt(i+1,j+1) .lt. k ) bdypt(i,j,k) = 1
            if (i < nlon .and. kmt(i+1,j) .lt. k ) bdypt(i,j,k) = 1
            if (j < nlat .and. kmt(i,j+1) .lt. k ) bdypt(i,j,k) = 1
          endif
        enddo
      enddo
    enddo

    if (.false.) then !adjust boundary points to gues/background
!     ainc_thresh = 5.0 ! deg C
!     WRITE(6,*) "Correcting boundary points to 0.0 with threshold = ", ainc_thresh
!     WHERE(bdypt > 0 .and. abs(v3d1(:,:,:,3) - v3d2(:,:,:,3)) > ainc_thresh)
!       do n=1,nv3d
!         v3d(:,:,:,n) = 0.0
!       enddo
!     ENDWHERE
!     WHERE(bdypt(:,:,1) > 0 .and. abs(v3d1(:,:,1,3) - v3d2(:,:,1,3)) > ainc_thresh)
!       do n=1,nv2d
!         v2d(:,:,n) = 0.0
!       enddo
!     ENDWHERE
    endif

    ALLOCATE(bdypt3d(nlon,nlat,nlev,nv3d),bdypt2d(nlon,nlat,nv2d))

    bdypt3d(:,:,:,1) = bdypt
    bdypt3d(:,:,:,2) = bdypt
    bdypt3d(:,:,:,3) = bdypt
    bdypt3d(:,:,:,4) = bdypt
    bdypt2d(:,:,1) = bdypt(:,:,1)
    CALL write_bingrd4('bdypt.grd',bdypt3d,bdypt2d)
    DEALLOCATE(bdypt3d,bdypt2d)

    !check for blank analysis gridpoints
    if (.true.) then
    ALLOCATE(blank3d1(nlon,nlat,nlev,nv3d),blank2d1(nlon,nlat,nv2d))
    ALLOCATE(blank3d2(nlon,nlat,nlev,nv3d),blank2d2(nlon,nlat,nv2d))
    blank3d1 = 0
    blank3d2 = 0
    blank2d1 = 0
    blank2d2 = 0
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
          if ( kmt(i,j) > 0 .and. &
             & v3d1(i,j,k,3) == 0.0 .and. v3d2(i,j,k,3) /= 0.0 ) then
            WRITE(6,*) "MAJOR WARNING: analysis point is empty where gues point is nonzero."
            if (.true. .and. bdypt(i,j,k) > 0) then
              WRITE(6,*) "Boundary point,"
              WRITE(6,*) "Adjusting analysis to gues value..."
              v3d(i,j,k,:) = 0.0 !STEVE: updating with gues/background data since it's missing analysis here.
            endif
            blank3d1(i,j,k,:) = 1
            blank2d1(i,j,:) = 1
            WRITE(6,*) "i,j,k = ", i,j,k
            WRITE(6,*) "lon = ", lon(i)
            WRITE(6,*) "lat = ", lat(j)
            WRITE(6,*) "kmt(i,j) = ", kmt(i,j)
            if (i>1 .and. j>1) WRITE(6,*) "kmt(i-1,j-1) = ", kmt(i-1,j-1)
            if (i>1) WRITE(6,*) "kmt(i-1,j) = ", kmt(i-1,j)
            if (j>1) WRITE(6,*) "kmt(i,j-1) = ", kmt(i,j-1)
            if (i<nlon .and. j<nlat) WRITE(6,*) "kmt(i+1,j+1) = ", kmt(i+1,j+1)
            if (i<nlon) WRITE(6,*) "kmt(i+1,j) = ", kmt(i+1,j)
            if (j<nlat) WRITE(6,*) "kmt(i,j+1) = ", kmt(i,j+1)
            WRITE(6,*) "anal: v3d1(i,j,k,3) = ", v3d1(i,j,k,3)
            WRITE(6,*) "back: v3d2(i,j,k,3) = ", v3d2(i,j,k,3)
          elseif ( kmt(i,j) > 0 .and. &
                 & v3d2(i,j,k,3) == 0.0 .and. v3d1(i,j,k,3) /= 0.0 ) then
            WRITE(6,*) "MAJOR WARNING: gues point is empty where analysis point is nonzero."
            if (.true. .and. bdypt(i,j,k) > 0) then
              WRITE(6,*) "Boundary point,"
              WRITE(6,*) "Adjusting analysis to gues value..."
              v3d(i,j,k,:) = 0.0 !STEVE: updating with gues/background data since it's equal to 0.0 here (indicates land)
            endif
            blank3d2(i,j,k,:) = 1
            blank2d2(i,j,:) = 1
            WRITE(6,*) "i,j,k = ", i,j,k
            WRITE(6,*) "lon = ", lon(i)
            WRITE(6,*) "lat = ", lat(j)
            WRITE(6,*) "kmt(i,j) = ", kmt(i,j)
            if (i>1 .and. j>1) WRITE(6,*) "kmt(i-1,j-1) = ", kmt(i-1,j-1)
            if (i>1) WRITE(6,*) "kmt(i-1,j) = ", kmt(i-1,j)
            if (j>1) WRITE(6,*) "kmt(i,j-1) = ", kmt(i,j-1)
            if (i<nlon .and. j<nlat) WRITE(6,*) "kmt(i+1,j+1) = ", kmt(i+1,j+1)
            if (i<nlon) WRITE(6,*) "kmt(i+1,j) = ", kmt(i+1,j)
            if (j<nlat) WRITE(6,*) "kmt(i,j+1) = ", kmt(i,j+1)
            WRITE(6,*) "anal: v3d1(i,j,k,3) = ", v3d1(i,j,k,3)
            WRITE(6,*) "back: v3d2(i,k,k,3) = ", v3d2(i,j,k,3)
          endif
        enddo
      enddo
    enddo
    CALL write_bingrd4('blank_anal.grd',blank3d1,blank2d1)
    CALL write_bingrd4('blank_gues.grd',blank3d2,blank2d2)
    DEALLOCATE(blank3d1,blank2d1,blank3d2,blank2d2)
    endif

  ! APPLY AINC LIMITS
  ! STEVE: new feature, 12/31/10
  ! The analysis increments are occasionally very large due to the adaptive inflation,
  ! therefore a hard limit can be imposed here when creating the correction file.
  ainc_max_temp = days/deltat !deg C
  ainc_min_temp = -ainc_max_temp
  ainc_max_salt = (days/2)/deltat !psu
  ainc_min_salt = -ainc_max_salt
  !STEVE: may need to add limits to increments on U and V velocities

  ! Where temp ainc > ainc_max_temp, ainc = ainc_max_temp
  WHERE (v3d(:,:,:,3) > ainc_max_temp) v3d(:,:,:,3) = ainc_max_temp
  WHERE (v3d(:,:,:,3) < ainc_min_temp) v3d(:,:,:,3) = ainc_min_temp

  ! Where salt ainc > ainc_max_salt, ainc = ainc_max_salt
  WHERE (v3d(:,:,:,4) > ainc_max_salt) v3d(:,:,:,4) = ainc_max_salt
  WHERE (v3d(:,:,:,4) < ainc_min_salt) v3d(:,:,:,4) = ainc_min_salt

  !STEVE: ERROR CHECKS:
  if (dodebug) then
!!  WRITE(6,*) "Mean sfc background temp:    ", sum(v3d2(:,:,1,3))/(nlon*nlat)
!!  WRITE(6,*) "Mean sfc analysis_f temp:    ", sum(v3d3(:,:,1,3))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis temp:      ", sum(v3d1(:,:,1,3))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis increment: ", sum(v3d(:,:,1,3))/(nlon*nlat) *deltat
!  WRITE(6,*) " "
!  WRITE(6,*) "Mean sfc background salt:    ", sum(v3d2(:,:,1,4))/(nlon*nlat)
!!  WRITE(6,*) "Mean sfc analysis_f salt:    ", sum(v3d3(:,:,1,4))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis salt:      ", sum(v3d1(:,:,1,4))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis increment: ", sum(v3d(:,:,1,4))/(nlon*nlat) *deltat
!  WRITE(6,*) " "
!  WRITE(6,*) "Mean sfc background U:       ", sum(v3d2(:,:,1,1))/(nlon*nlat)
!!  WRITE(6,*) "Mean sfc analysis_f U:       ", sum(v3d3(:,:,1,1))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis U:         ", sum(v3d1(:,:,1,1))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis increment: ", sum(v3d(:,:,1,1))/(nlon*nlat) *deltat
!  WRITE(6,*) " "
!  WRITE(6,*) "Mean sfc background V:       ", sum(v3d2(:,:,1,2))/(nlon*nlat)
!!  WRITE(6,*) "Mean sfc analysis_f V:       ", sum(v3d3(:,:,1,2))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis V:         ", sum(v3d1(:,:,1,2))/(nlon*nlat)
!  WRITE(6,*) "Mean sfc analysis increment: ", sum(v3d(:,:,1,2))/(nlon*nlat) *deltat
  WRITE(6,*) " "
  WRITE(6,*) "At 220,65:"
! WRITE(6,*) "sfc background temp:    ", (v3d2(220,65,1,3))
! WRITE(6,*) "sfc analysis temp:      ", (v3d1(220,65,1,3))
! WRITE(6,*) "sfc analysis increment: ", (v3d(220,65,1,3))*deltat
  WRITE(6,*) " "
! WRITE(6,*) "sfc background salt:    ", (v3d2(220,65,1,4))
! WRITE(6,*) "sfc analysis salt:      ", (v3d1(220,65,1,4))
! WRITE(6,*) "sfc analysis increment: ", (v3d(220,65,1,4))*deltat
  WRITE(6,*) " "
! WRITE(6,*) "sfc background U:       ", (v3d2(220,65,1,1))
! WRITE(6,*) "sfc analysis U:         ", (v3d1(220,65,1,1))
! WRITE(6,*) "sfc analysis increment: ", (v3d(220,65,1,1))*deltat
  WRITE(6,*) " "
! WRITE(6,*) "sfc background V:       ", (v3d2(220,65,1,2))
! WRITE(6,*) "sfc analysis V:         ", (v3d1(220,65,1,2))
! WRITE(6,*) "sfc analysis increment: ", (v3d(220,65,1,2))*deltat
  WRITE(6,*) " "
  endif

  ainc_thresh = huge(1.0) !STEVE: matching this with 'gross_error' from common_obs_mom2 to test, should have same effect but don't...
  errloc = MAXLOC(ABS(v3d(:,:,:,3)))
  if (MAXVAL(ABS(v3d(:,:,:,3)))*deltat > ainc_thresh ) then
    WRITE(6,*) "common_mom4.f90::grd_to_cor:: "
    !STEVE: I want to keep running to see what happens...
    if (.true. .and. MAXVAL(ABS(v3d(:,:,:,3)))*deltat > ainc_thresh*2 ) then
      WRITE(6,*) "ERROR: max temp analysis increment = ", MAXVAL(ABS(v3d(:,:,:,3)))*deltat
      WRITE(6,*) "at location: MAXLOC(ABS(v3d(:,:,:,3))) = ", MAXLOC(ABS(v3d(:,:,:,3)))
      WRITE(6,*) "kmt(errloc) = ", kmt(errloc(1),errloc(2))
      WRITE(6,*) "bdypt(errloc) = ", bdypt(errloc(1),errloc(2),errloc(3))
      WRITE(6,*) "lon = ", lon(errloc(1))
      WRITE(6,*) "lat = ", lat(errloc(2))
      WRITE(6,*) "LETKF maximum analysis increment for temperature is over twice the warning threshold."
      WRITE(6,*) "ERROR:: STOP..."
      STOP(14)
    else
      !STEVE: I just wanted two levels of checks, 1 for a warning, and 1 for an error...
      WRITE(6,*) "WARNING: max temp analysis increment = ", MAXVAL(ABS(v3d(:,:,:,3)))*deltat
      WRITE(6,*) "at location: MAXLOC(ABS(v3d(:,:,:,3))) = ", MAXLOC(ABS(v3d(:,:,:,3)))
      errloc = MAXLOC(ABS(v3d(:,:,:,3)))
      WRITE(6,*) "kmt(errloc) = ", kmt(errloc(1),errloc(2))
      WRITE(6,*) "lon = ", lon(errloc(1))
      WRITE(6,*) "lat = ", lat(errloc(2))
      !STEVE: added for testing...
      WRITE(6,*) "ERROR:: STOP..."
      STOP(15)
    endif
  endif

    DEALLOCATE(bdypt)
  endif

  if (sum(abs(v3d(:,:,:,3))) < 0.0001/deltat ) then
    WRITE(6,*) "grd_to_cor :: No update produced for mom4 temp_increment.nc temperature..."
    WRITE(6,*) "WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!"
    !stop 1
  endif

  ! Write out the appropriate correctors for mom4p1
  CALL write_inc4('temp',v3d,v2d)
  CALL write_inc4('salt',v3d,v2d)
  CALL write_inc4('u',v3d,v2d)
  CALL write_inc4('v',v3d,v2d)
  CALL write_inc4('ssh',v3d,v2d)

END SUBROUTINE grd_to_cor

!----------------------------------------------------
subroutine write_correctors4(filename,temp,salt,ucur,vcur)
  character(*), INTENT(IN) :: filename
  REAL(r_sngl), dimension(nlon,nlat,nlev), INTENT(IN) :: temp, salt
  REAL(r_sngl), dimension(nlon,nlat,nlev), INTENT(IN) :: ucur, vcur
  integer :: numots ! number of time steps
  !integer, parameter :: ibm_rec=(nlon+2)*nlat*4 !STEVE: because mom4 requires the cyclic overlap longitudes
  integer, parameter :: ibm_rec = nlon*nlat*4 !STEVE: because mom4 requires the cyclic overlap longitudes
  integer :: fid = 79, irec, k
  WRITE(6,*) "write_correctors4: writing...", filename

  open(unit=fid,file=filename,status='new', &
       form='unformatted',access='direct',recl=ibm_rec)

  irec = 0
  do k=1,nlev
    irec = k
    ! mom2 correctors don't include the cyclic overlap longitudes
    write(fid,rec=irec) temp(:,:,k)
  enddo

  do k=1,nlev
    irec = k + nlev
    write(fid,rec=irec) salt(:,:,k)
  enddo

  do k=1,nlev
    irec = k + 2*nlev
    write(fid,rec=irec) ucur(:,:,k)
  enddo

  do k=1,nlev
    irec = k + 3*nlev
    write(fid,rec=irec) vcur(:,:,k)
  enddo

end subroutine write_correctors4
!----------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE set_kmt(file)
  CHARACTER(*), INTENT(IN) :: file
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)
  INTEGER :: i,j,k,n

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  CALL read_bingrd(file,v3d,v2d)

  kmt = 0
  do k=nlev,1,-1
    do j=1,nlat
      do i=1,nlon
        if ( kmt(i,j) < 1 .and. v3d(i,j,k,3) > 0.0 ) then !STEVE: wanted MAX(v3d(i,j,k,:)), but salinity is 35 on land due to unit conversion
          kmt(i,j) = k
        endif
      enddo
    enddo
  enddo

  DEALLOCATE(v3d,v2d)

END SUBROUTINE set_kmt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_kmt(v3d,v2d)
  REAL(r_size), INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size), INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: nkmterr = 0
  INTEGER :: i,j,k,n

  do k=nlev,1,-1
    do j=1,nlat
      do i=1,nlon
          if ( kmt(i,j) < 1 .and. v3d(i,j,k,3) > 0.0 ) then !STEVE: wanted MAX(v3d(i,j,k,:)), but salinity is 35 on land due to unit conversion
            nkmterr = nkmterr+1
            kmt(i,j) = k
          endif
      enddo
    enddo
  enddo
  WRITE(6,*) "common_mom4::update_kmt: KMT entries updated: ", nkmterr

  return
END SUBROUTINE update_kmt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_sfc_grid(infile,slon,slat,nslon,nslat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use netcdf
CHARACTER(*), INTENT(IN) :: infile
REAL(r_size), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: slon, slat
INTEGER, INTENT(OUT) :: nslon, nslat
INTEGER :: dimid, varid, ncid
INTEGER, DIMENSION(1) :: dimIDs
INTEGER :: ntimes
REAL(r_size), DIMENSION(:), ALLOCATABLE :: time
LOGICAL :: dodebug=.false.

! Open a surface flux data file
if (dodebug) WRITE(6,'(2A)') '  >> accessing file: ', infile
call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )

! Inquire the variable
if (dodebug) print *, "Calling NF90_INQ_VARID and nf90_inquire_variable..."
call check( NF90_INQ_VARID(ncid,'XAXIS_1',varid) )   ! Longitude for T-cell
call check( nf90_inquire_variable(ncid, varid, dimids = dimIDs) )
dimid=dimids(1)

! Inquire the dimension
if (dodebug) print *, "Calling nf90_inquire_dimension..."
call check( nf90_inquire_dimension(ncid, dimid, len = nslon) )
ALLOCATE(slon(nslon))

! Read lon arrays from file
if (dodebug) print *, "Calling NF90_GET_VAR..."
call check( NF90_GET_VAR(ncid,varid,slon) )
WRITE(6,*) "slon(1) = ", slon(1)
WRITE(6,*) "slon(nslon) = ", slon(nslon)

! Inquire the variable
if (dodebug) print *, "Calling NF90_INQ_VARID and nf90_inquire_variable..."
call check( NF90_INQ_VARID(ncid,'YAXIS_1',varid) )   ! Latitude for T-cell
call check( nf90_inquire_variable(ncid, varid, dimids = dimIDs) )
dimid=dimIDs(1)

! Inquire the dimension
if (dodebug) print *, "Calling nf90_inquire_dimension..."
call check( nf90_inquire_dimension(ncid, dimid, len = nslat) )
ALLOCATE(slat(nslat))

! Read lat arrays from file
if (dodebug) print *, "Calling NF90_GET_VAR..."
call check( NF90_GET_VAR(ncid,varid,slat) )
if (dodebug) WRITE(6,*) "slat(1) = ", slat(1)
if (dodebug) WRITE(6,*) "slat(nslat) = ", slat(nslat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read time arrays from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inquire the variable
if (dodebug) print *, "Calling NF90_INQ_VARID and nf90_inquire_variable..."
call check( NF90_INQ_VARID(ncid,'TIME',varid) )   ! Latitude for T-cell
call check( nf90_inquire_variable(ncid, varid, dimids = dimIDs) )
dimid=dimIDs(1)

! Inquire the dimension
if (dodebug) print *, "Calling nf90_inquire_dimension..."
call check( nf90_inquire_dimension(ncid, dimid, len = ntimes) )
ALLOCATE(time(ntimes))

if (dodebug) print *, "Calling NF90_GET_VAR..."
call check( NF90_GET_VAR(ncid,varid,time) )
if (dodebug) print *, "time(1) = ", time(1)
if (dodebug) print *, "time(nslon) = ", time(ntimes)

! Close netcdf file
call check( NF90_CLOSE(ncid) )

END SUBROUTINE read_sfc_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_sfc_data(infile,varname,ti,sfc_data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use netcdf
CHARACTER(*), INTENT(IN) :: infile, varname
INTEGER, INTENT(IN) :: ti ! timeslot
REAL(r_size), DIMENSION(:,:), INTENT(OUT) :: sfc_data
INTEGER :: dimid, varid, ncid
INTEGER, DIMENSION(3) :: dimIDs
REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
INTEGER :: ntimes
LOGICAL :: dodebug=.true.

! Open a surface flux data file
if (dodebug) WRITE(6,'(2A)') '  >> accessing file: ', infile
if (dodebug) print *, "read_sfc_data:: Opening: file, varname, ti = ", trim(infile), trim(varname), ti
call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read time arrays from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inquire the time variable
if (dodebug) print *, "Calling NF90_INQ_VARID and nf90_inquire_variable..."
call check( NF90_INQ_VARID(ncid,'TIME',varid) )   ! Time
! Inquire the time dimension
call check( nf90_inquire_variable(ncid, varid, dimids = dimIDs) )
dimid=dimIDs(1)

! Inquire the time dimension length
if (dodebug) print *, "Calling nf90_inquire_dimension..."
call check( nf90_inquire_dimension(ncid, dimid, len = ntimes) )
if (dodebug) print *, "The number of time steps are: ", ntimes

! Inquire the variable
if (dodebug) print *, "Calling NF90_INQ_VARID... for variable ",varname
call check( NF90_INQ_VARID(ncid,trim(varname),varid) )   ! Variable id

! Read data array from file
if (dodebug) print *, "Calling NF90_GET_VAR..."
if (dodebug) print *, "Allocate data4D..."
ALLOCATE(data4D(nslon,nslat,ntimes))
call check( NF90_GET_VAR(ncid,varid,data4D) ) ! Variable data
if (dodebug) print *, "data4D(1,1,ti) = ", data4D(1,1,ti)
if (dodebug) print *, "data4D(nslon,nslat,ti) = ", data4D(nslon,nslat,ti)
if (dodebug) print *, "Assigning data4D timeslot ",ti," to sfc_data"
sfc_data = data4D(:,:,ti)
if (dodebug) print *, "sfc_data(1,1) = ", sfc_data(1,1)
if (dodebug) print *, "sfc_data(nslon,nslat) = ", sfc_data(nslon,nslat)

! Close netcdf file
DEALLOCATE(data4D)
call check( NF90_CLOSE(ncid) )

END SUBROUTINE read_sfc_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_coeff(infile,xi,yi,xc,yc,lon,lat,NX,NY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use netcdf
IMPLICIT NONE
CHARACTER(*), INTENT(IN) :: infile
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: xi, yi
REAL(r_size), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: xc, yc
REAL(r_size), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: lon, lat
INTEGER, INTENT(OUT) :: NX, NY
! When we create netCDF files, variables and dimensions, we get back an ID for
! each one.
integer :: ncid, varid, dimids(2)
integer :: x_dimid, y_dimid, x_varid, y_varid
INTEGER :: varid_xi, varid_yi, varid_xc, varid_yc
LOGICAL :: dodebug=.false.

! Open the netCDF file.
if (dodebug) print *, "Opening file: ", infile
call check( nf90_open(infile, NF90_NOWRITE, ncid) )

! Inquire the variable
if (dodebug) print *, "Calling NF90_INQ_VARID and nf90_inquire_variable..."
call check( nf90_inq_dimid(ncid, 'X', x_dimid) )
call check( nf90_inq_dimid(ncid, 'Y', y_dimid) )

! Inquire the dimension
if (dodebug) print *, "Calling nf90_inquire_dimension..."
call check( nf90_inquire_dimension(ncid, x_dimid, len = NX) )
call check( nf90_inquire_dimension(ncid, y_dimid, len = NY) )
ALLOCATE(lon(NX))
ALLOCATE(lat(NY))
ALLOCATE(xi(NX))
ALLOCATE(yi(NY))
ALLOCATE(xc(NX))
ALLOCATE(yc(NY))

! Inquire the variable ids:
call check( nf90_inq_varid(ncid, 'xi', varid_xi) )
call check( nf90_inq_varid(ncid, 'yi', varid_yi) )
call check( nf90_inq_varid(ncid, 'xc', varid_xc) )
call check( nf90_inq_varid(ncid, 'yc', varid_yc) )
call check( nf90_inq_varid(ncid, 'X', x_varid) )
call check( nf90_inq_varid(ncid, 'Y', y_varid) )

! Write the coordinate data to the file.
if (dodebug) print *, "Calling nf90_get_var x and y..."
call check( nf90_get_var(ncid, x_varid, lon) )
call check( nf90_get_var(ncid, y_varid, lat) )

! Write the data to the file.
if (dodebug) print *, "Calling nf90_get_var..."
call check( nf90_get_var(ncid, varid_xi, xi))
call check( nf90_get_var(ncid, varid_yi, yi))
call check( nf90_get_var(ncid, varid_xc, xc))
call check( nf90_get_var(ncid, varid_yc, yc))

! Close the file. This frees up any internal netCDF resources
if (dodebug) print *, "Calling nf90_close..."
call check( nf90_close(ncid) )

END SUBROUTINE read_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE coeff_sfc2model(slon,slat,mlon,mlat,xi,yi,xc,yc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
REAL(r_size), DIMENSION(:), INTENT(IN) :: slon, slat, mlon, mlat
INTEGER, DIMENSION(:), INTENT(OUT) :: xi, yi
REAL(r_size), DIMENSION(:), INTENT(OUT) :: xc, yc
INTEGER :: i,j,ii,jj

! for each model longitude,
do i=1,nmlon
  print *, "i = ", i
  ! loop through each sfc longitude until we find the sfc cell
  xc(i) = 0
  do ii=1,nslon-1
    !STEVE: 'MODULO' converts to positive value
    if (MODULO(mlon(i),360.0) < slon(ii+1)) then
!     print *, "mlon(i) = ", mlon(i)
!     print *, "MODULO(mlon(i),360.0) = ", MODULO(mlon(i),360.0)
!     STOP(1)
      xc(i) = (MODULO(mlon(i),360.0) - slon(ii))/(slon(ii+1) - slon(ii))
      xi(i) = ii
      exit
    elseif (ii == nslon-1) then !slon: (0-359.375) !mlon: (-279.75 to 79.75)
      xc(i) = (MODULO(mlon(i),360.0) - slon(nslon))/(slon(1) - MOD(slon(nslon)-360.0,360.0))
      xi(i) = nslon
    endif
  enddo
  print *, "xc(i) = ", xc(i)
  !STEVE: debug
  if (xc(i) > 1) then ! .OR. xc(i) < 0.00001) then
    print *, "ii, i = ", ii, i
    print *, "xc(i) = ", xc(i)
    print *, "MODULO(mlon(i),360.0) < slon(ii+1)"
    print *, "mlon(i) = ", mlon(i)
    print *, "MODULO(mlon(i),360.0) = ", MODULO(mlon(i),360.0)
    print *, "slon(ii) = ", slon(ii)
    print *, "slon(ii+1) = ", slon(ii+1)
    print *, "slon(ii) = ", slon(ii)
    STOP(16)
  endif
  if ( .false. ) then !xi(i) .eq. nslon ) then ! .OR. xc(i) < 0.00001) then
    print *, "ii, i = ", ii, i
    print *, "xc(i) = ", xc(i)
    print *, "MODULO(mlon(i),360.0) < slon(ii+1)"
    print *, "mlon(i) = ", mlon(i)
    print *, "MODULO(mlon(i),360.0) = ", MODULO(mlon(i),360.0)
    print *, "slon(nslon) = ", slon(nslon)
    print *, "slon(1) = ", slon(1)
    print *, "slon(nslon) = ", slon(nslon)
    print *, "MOD(slon(nslon)-360.0,360.0) = ", MOD(slon(nslon)-360.0,360.0)
    STOP(17)
  endif
enddo

do j=1,nmlat !slon: (-89.522 to 89.522) !mlon: (-80.75 to 89.75)
  print *, "j = ", j
  ! loop through each sfc latitude until we find the sfc cell
  yc(j) = 0
  !STEVE:  If the model lat is within the max sfc lat, then:
  if (mlat(j) .le. slat(nslat)) then
    do jj=1,nslat-1
      if (mlat(j) < slat(jj+1)) then
        yc(j) = (slat(jj+1) - mlat(j))/(slat(jj+1) - slat(jj))
        yi(j) = jj
        exit
      endif
    enddo
  else !STEVE: otherwise, loop over the pole and interpolate
    yc(j) = (90.0 - slat(nslat) + 90.0 - mlat(j))/(90.0 - slat(nslat))
    yi(j) = nslat
  endif
  print *, "yc(j) = ", yc(j)
enddo

END SUBROUTINE coeff_sfc2model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE coeff_model2sfc(mlon,mlat,slon,slat,xi,yi,xc,yc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
REAL(r_size), DIMENSION(:), INTENT(IN) :: mlon, mlat, slon, slat
INTEGER, DIMENSION(:), INTENT(OUT) :: xi, yi
REAL(r_size), DIMENSION(:), INTENT(OUT) :: xc, yc
INTEGER :: i, j, ii, jj
REAL(r_size), PARAMETER :: offsetL = -279.75 , offsetR = 79.75
REAL(r_size) :: mlon_corr ! mlon corrected
INTEGER :: eidx, idx, nidx
INTEGER :: mzeroidx = 560 !563 !563-564
INTEGER :: s1steidx = 129
LOGICAL :: dodebug=.true.

!STEVE: if xx or yc > 1, then there is probably some kind of model error
!present.

!STEVE: do this in two pieces, because of the odd longitude grid (due to
!tripolar grid)
! for each sfc longitude,
do i=1,s1steidx !nslon
  ! loop through each sfc longitude until we find the sfc cell
  xc(i) = 0
  do ii=mzeroidx,nmlon-1
    !STEVE: 'MODULO' converts to positive value
    if ( slon(i) < MODULO(mlon(ii+1),360.0) ) then
      if (dodebug) then
        print *, "i = ", i
        print *, "MODULO(mlon(ii),360.0) = ", MODULO(mlon(ii),360.0)
        print *, "slon(i) = ", slon(i)
        print *, "MODULO(mlon(ii+1),360.0) = ", MODULO(mlon(ii+1),360.0)
      endif
      if (i > 1) then
        xc(i) = (MODULO(mlon(ii+1),360.0) - slon(i) )/(mlon(ii+1) - MODULO(mlon(ii),360.0))
      else
        xc(i) = (mlon(ii+1) - slon(i) )/(mlon(ii+1) - mlon(ii))
      endif
      xi(i) = ii
      eidx=i
      exit
    elseif (ii == nmlon-1 .AND. slon(i) < MODULO(mlon(1),360.0)) then
!     if (dodebug) then
        print *, "i = ", i
        print *, "MODULO(mlon(nmlon),360.0) = ", MODULO(mlon(nmlon),360.0)
        print *, "slon(i) = ", slon(i)
        print *, "MODULO(mlon(1),360.0) = ", MODULO(mlon(1),360.0)
      endif
      xc(i) = (MODULO(mlon(1),360.0) - slon(i))/(MODULO(mlon(1),360.0) - MODULO(mlon(nmlon),360.0))
      xi(i) = nmlon
      eidx=i
!   endif
  enddo
  if (dodebug) print *, "xc(i) = ", xc(i)
  if (xc(i) < 0.0001) then
    print *, "i, ii = ", i, ii
    STOP(18)
  endif
enddo

if (dodebug) print *, "eidx = ", eidx
! Next, do the part from offsetR to 360.0
do i=eidx+1,nslon
  ! loop through each sfc longitude until we find the sfc cell
  xc(i) = 0
  do ii=1,mzeroidx-1
    !STEVE: 'MODULO' converts to positive value
    if ( slon(i) < MODULO(mlon(ii+1),360.0) ) then
      if (dodebug) then
        print *, "i = ", i
        print *, "MODULO(mlon(ii),360.0) = ", MODULO(mlon(ii),360.0)
        print *, "slon(i) = ", slon(i)
        print *, "MODULO(mlon(ii+1),360.0) = ", MODULO(mlon(ii+1),360.0)
      endif
      xc(i) = (MODULO(mlon(ii+1),360.0) - slon(i))/(mlon(ii+1) - mlon(ii))
      xi(i) = ii
      exit
    endif
  enddo
  print *, "xc(i) = ", xc(i)
  if (xc(i) < 0.0001) then
    print *, "i, ii = ", i, ii
    STOP(19)
  endif
enddo

do j=1,nslat
! print *, "j = ", j
  ! loop through each sfc latitude until we find the sfc cell
  yc(j) = 0
  do jj=1,nmlat-1
    if (slat(j) < mlat(jj+1)) then
      if (dodebug) then
        print *, "j = ", j
        print *, "mlat(jj) = ", mlat(jj)
        print *, "slat(j) = ", slat(j)
        print *, "mlat(jj+1) = ", mlat(jj+1)
      endif
      yc(j) = (mlat(jj+1) - slat(j) )/(mlat(jj+1) - mlat(jj))
      yi(j) = jj
      exit
    endif
  enddo
  if (dodebug) print *, "yc(j) = ", yc(j)
enddo

END SUBROUTINE coeff_model2sfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE itpl_sfc2model(sfc_data,xi,yi,xc,yc,model_data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
REAL(r_size), DIMENSION(:,:), INTENT(IN) :: sfc_data
INTEGER, DIMENSION(:), INTENT(IN) :: xi, yi
REAL(r_size), DIMENSION(:), INTENT(IN) :: xc, yc
REAL(r_size), DIMENSION(:,:), INTENT(OUT) :: model_data
INTEGER :: i, j
REAL(r_size) :: f1, f2

! Apply a bilinear interpolation to get model_data:
print *, "Applying bilinear interpolation in itpl_sfc2model..."
model_data = 0
do j=1,nmlat-1
  do i=1,nmlon-1
!   print *, "i,j, sfc_data(i,j), xc, yc = ", i,j, sfc_data(i,j), xc(i), yc(j)
    f1 =  xc(i)*sfc_data(xi(i),yi(j))   + (1 - xc(i))*sfc_data(xi(i+1),yi(j))
    f2 =  xc(i)*sfc_data(xi(i),yi(j+1)) + (1 - xc(i))*sfc_data(xi(i+1),yi(j+1))
    model_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  enddo
enddo

!STEVE: handle boundary cases:

!STEVE: this case is complicated by the fact that the grids do not go to 90.00
!degrees. We have to loop over the pole and interpolate to the point that is
!180.0 degrees wrapped around the globe.
!Since the grids are evenly spaced in longitude, I'm just adding half of the
!grid points here. For a more general case, we'd need to find the exact
!longitudes.
!do j=nmlat-1,nmlat
j=nmlat
do i=1,nmlon-1
! f1 =  xc(i)*sfc_data(xi(i),yi(j))   + (1 - xc(i))*sfc_data(xi(i+1),yi(j))
! ii = MODULO(i+FLOOR(nslon/2.0)-1,nslon)+1 !STEVE: jump over the pole on slat
! grid
! f2 =  xc(ii)*sfc_data(xi(ii),yi(1)) + (1 - xc(ii))*sfc_data(xi(ii+1),yi(1))
! model_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  model_data(i,j) = model_data(i,j-1)
enddo
!enddo

!STEVE: This is the border where the mlon grid loops at 80/-280 (ish)
i=nmlon
do j=1,nmlat-1
  f1 =  xc(i)*sfc_data(xi(i),yi(j))   + (1 - xc(i))*sfc_data(xi(1),yi(j))
  f2 =  xc(i)*sfc_data(xi(i),yi(j+1)) + (1 - xc(i))*sfc_data(xi(1),yi(j+1))
  model_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
enddo

!STEVE: loop the corner of the grid in both indices
j=nmlat
i=nmlon
!f1 =  xc(i)*sfc_data(xi(i),yi(j))   + (1 - xc(i))*sfc_data(xi(1),yi(j))
!f2 =  xc(i)*sfc_data(xi(i),yi(1)) + (1 - xc(i))*sfc_data(xi(1),yi(1))
!model_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
model_data(i,j) = model_data(i,j-1)

!STEVE: the Southern pole boundary is not as much of an issue since there is no
!ocean there. Interpolating back from model to sfc may require more attention.

END SUBROUTINE itpl_sfc2model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE itpl_model2sfc(model_data,xi,yi,xc,yc,sfc_data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
REAL(r_size), DIMENSION(:,:), INTENT(IN) :: model_data
INTEGER, DIMENSION(:), INTENT(IN) :: xi, yi
REAL(r_size), DIMENSION(:), INTENT(IN) :: xc, yc
REAL(r_size), DIMENSION(:,:), INTENT(OUT) :: sfc_data
REAL(r_size) :: f1, f2
INTEGER :: i, j

! Apply a bilinear interpolation to get sfc_data:
sfc_data = 0
print *, "Applying bilinear interpolation in itpl_model2sfc..."
!print *, "nslon, nslat = ", nslon, nslat
do j=1,nslat-1
  do i=1,nslon-1
!   print *, "i,j, model_data(i,j), xc, yc = ", i,j, model_data(i,j), xc(i),
!   yc(j)
!   print *, "i,j, model_data(i,j), xi(i), yi(j) = ", i,j, model_data(i,j),
!   xi(i), yi(j)
!   print *, "i,j, model_data(i,j), xi(i+1), yi(j+1) = ", i,j, model_data(i,j),
!   xi(i+1), yi(j+1)
    f1 = xc(i)*model_data(xi(i),yi(j))   + (1 - xc(i))*model_data(xi(i+1),yi(j))
    f2 = xc(i)*model_data(xi(i),yi(j+1)) + (1 - xc(i))*model_data(xi(i+1),yi(j+1))
    sfc_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  enddo
enddo

!STEVE: handle boundary cases:

!STEVE: the model lat range goes higher N than the sfc lat range, so this should
!be ok.
j=nslat
do i=1,nslon-1
! f1 = xc(i)*model_data(xi(i),yi(j))   + (1 - xc(i))*model_data(xi(i+1),yi(j))
! f2 = xc(i)*model_data(xi(i),yi(1)) + (1 - xc(i))*model_data(xi(i+1),yi(1))
! sfc_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  sfc_data(i,j) =  sfc_data(i,j-1)
enddo

! This is covered similarly to the sfc2model case
i=nslon
do j=1,nslat-1
  f1 = xc(i)*model_data(xi(i),yi(j))   + (1 - xc(i))*model_data(xi(1),yi(j))
  f2 = xc(i)*model_data(xi(i),yi(j+1)) + (1 - xc(i))*model_data(xi(1),yi(j+1))
  sfc_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  if (.false.) then
    print *, "model2sfc:: i=nslon, j = ", j
    print *, "f1 = ", f1
    print *, "f2 = ", f2
    print *, "xi(i) = ", xi(i)
    print *, "xi(1) = ", xi(1)
    print *, "yi(j) = ", yi(j)
    print *, "yi(j+1) = ", yi(j+1)
    print *, "sfc_data(i,j) = ", sfc_data(i,j)
  endif
enddo

! Loop the corner of the grid in both indices
j=nslat
i=nslon
!f1 = xc(i)*model_data(xi(i),yi(j))   + (1 - xc(i))*model_data(xi(1),yi(j))
!f2 = xc(i)*model_data(xi(i),yi(1)) + (1 - xc(i))*model_data(xi(1),yi(1))
!sfc_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
sfc_data(i,j) = sfc_data(i,j-1)

! The S pole may require special attention since the model grid ends inside the
! sfc grid.
! Loop over the S pole to interpolate to point at +180.0 degrees.
j=1
do i=1,nslon-1
! f1 = xc(i)*model_data(xi(i),yi(j))   + (1 - xc(i))*model_data(xi(i+1),yi(j))
! ii = MODULO(i+FLOOR(nmlon/2.0)-1,nmlon)+1 !STEVE: jump over the pole on slat
! grid
! f2 = xc(ii)*model_data(xi(ii),yi(1)) + (1 - xc(ii))*model_data(xi(ii+1),yi(1))
! sfc_data(i,j) =  yc(j)*f1 + (1 - yc(j))*f2
  sfc_data(i,j) = sfc_data(i,j+1)
enddo

END SUBROUTINE itpl_model2sfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_sfc_data(outfile,varname,islot,data_out,lon,lat,NX,NY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use netcdf
IMPLICIT NONE
CHARACTER(*), INTENT(IN) :: outfile, varname
REAL(r_size), DIMENSION(:,:,:), INTENT(IN) :: data_out !time,lon,lat
REAL(r_size), DIMENSION(:) :: lon, lat
INTEGER, INTENT(IN) :: NX, NY, islot
! When we create netCDF files, variables and dimensions, we get back an ID for
! each one.
INTEGER, PARAMETER :: NDIMS=3
integer :: ncid, varid, dimids(NDIMS)
integer :: x_dimid, y_dimid, x_varid, y_varid, t_dimid, t_varid
INTEGER, DIMENSION(NDIMS) :: start, count

! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
if (islot .eq. 0) then
  print *, "Creating file: ", outfile
  call check( nf90_create(outfile, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  print *, "Calling nf90_def_dim..."
  call check( nf90_def_dim(ncid, "XAXIS_1", NX, x_dimid) )
  call check( nf90_def_dim(ncid, "YAXIS_1", NY, y_dimid) )
  call check( nf90_def_dim(ncid, "TIME", NF90_UNLIMITED, t_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ x_dimid, y_dimid, t_dimid /)

  ! Define the variable. 
  print *, "Calling nf90_def_var x and y..."
  call check( nf90_def_var(ncid, 'XAXIS_1', NF90_REAL, x_dimid, x_varid) )
  call check( nf90_def_var(ncid, 'YAXIS_1', NF90_REAL, y_dimid, y_varid) )
  call check( nf90_def_var(ncid, 'TIME', NF90_REAL, t_dimid, t_varid) )
  print *, "Calling nf90_def_var..."
  call check( nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid) )

  ! Define addtribute
  print *, "Calling nf90_put_att..."
  call check( nf90_put_att(ncid, x_varid, "axis", 'x') )
  call check( nf90_put_att(ncid, x_varid, "units", 'degrees_east') )
  call check( nf90_put_att(ncid, y_varid, "axis", 'y') )
  call check( nf90_put_att(ncid, y_varid, "units", 'degrees_north') )
  call check( nf90_put_att(ncid, t_varid, "axis", 't') )
  call check( nf90_put_att(ncid, t_varid, "units", 'days since 2010-12-02 12:00:00') ) !STEVE: insert correct date

  ! End define mode.
  print *, "Calling nf90_enddef..."
  call check( nf90_enddef(ncid) )

else
  call check( nf90_open(outfile, NF90_WRITE, ncid) )
endif

! Write the coordinate data to the file.
print *, "Calling nf90_put_var x and y..."
call check( nf90_put_var(ncid, x_varid, lon) )
call check( nf90_put_var(ncid, y_varid, lat) )
print *, "Calling nf90_put_var t..."
call check( nf90_put_var(ncid, t_varid, islot) )

! Write the data to the file.
print *, "Calling nf90_put_var data_out..."
!call check( nf90_put_var(ncid, varid, TRANSPOSE(data_out)) )
start = (/ 1, 1, 1 /)
count = (/ NX, NY, 1 /)
call check( nf90_put_var(ncid, varid, data_out, start=start,count=count) )

! Close the file. This frees up any internal netCDF resources
print *, "Calling nf90_close..."
call check( nf90_close(ncid) )

END SUBROUTINE write_sfc_data


END MODULE common_mom4
