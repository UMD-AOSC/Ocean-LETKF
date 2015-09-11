MODULE common_mom4
!===============================================================================
! MODULE: common_obs
! 
! PURPOSE: Common Information for MOM4
!
! USES:
!   use common
!   use params_model
!   use vars_model
!   use params_letkf
!
! !PUBLIC TYPES:
!                 implicit none
!                 [save]
!
!                 <type declaration>
!     
! !PUBLIC MEMBER FUNCTIONS:
!           <function>                     ! Description      
!
! !PUBLIC DATA MEMBERS:
!           <type> :: <variable>           ! Variable description
!
! DESCRIPTION: 
!   This module reads all observation data and stores in appropriate data structures
!
! !REVISION HISTORY:
!   04/26/2011 Steve PENNY converted to OCEAN for use with MOM4
!   10/15/2004 Takemasa Miyoshi  created
!-------------------------------------------------------------------------------
! $Author: Steve Penny $
!===============================================================================
  USE common
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_SFCFLUX

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------------
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
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
  REAL(r_size), DIMENSION(nlon,nlat) :: SSHclm_m
  REAL(r_sngl) :: buf4_2d(nlon,nlat)
! For AMOC computation
  REAL(r_size) :: zb(nlev)
  REAL(r_size) :: dz(nlev)

!STEVE: the following is for assimilating surface fluxes, it is not currently
!
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


SUBROUTINE set_common_mom4
!===============================================================================
! Initialize the module
!===============================================================================
  USE netcdf
  IMPLICIT NONE
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
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)

  ! Close the grid_spec.nc file:
  call check( NF90_CLOSE(ncid) )

  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

END SUBROUTINE set_common_mom4


!------------------------------------------------------------------------------
! File I/O
!------------------------------------------------------------------------------


SUBROUTINE check(status)
!===============================================================================
! Check the error status of the netcdf command
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
END SUBROUTINE check


SUBROUTINE read_etaclm(SSHclm_file,SSHclm_m)
!===============================================================================
! Read in the model mean climatology (e.g. 1991-1999), the real climatology 
! assumed to be subtracted already from the observation data.
!===============================================================================
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


SUBROUTINE read_diag(infile,v3d,v2d)
!===============================================================================
! Read in a set of netcdf-format mom4 restart files
!===============================================================================
  USE netcdf
  IMPLICIT NONE
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

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the T/S netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_diag:: just opened file ", tsfile

  !!! t
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'temp',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable temp"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_t) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable temp"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_diag:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ",MAXVAL(v3d(:,:,k,iv3d_t))
    enddo
  endif
! !STEVE: end

  !!! s
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'salt',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable salt"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_s) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable salt"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_diag:: tsfile = ", tsfile
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
    WRITE(6,'(A)') 'netCDF OPEN ERROR in read_diag for ',uvfile
    STOP(7)
  END IF
  WRITE(6,*) "read_diag:: just opened file ", uvfile

  !!! u
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable u"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_u) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable u"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_diag:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ",MAXVAL(v3d(:,:,k,iv3d_u))
    enddo
  endif
! !STEVE: end

  !!! v
  buf4=0.0
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  call check( NF90_GET_VAR(ncid,varid,buf4) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable v"
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        v3d(i,j,k,iv3d_v) = REAL(buf4(i,j,k),r_size)
      END DO
    END DO
  END DO
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable v"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_diag:: uvfile = ", uvfile
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
    WRITE(6,*) "read_diag:: just opened file ", sffile

    !!! SST
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'t_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable sfc temp"
    DO j=1,nlat
      DO i=1,nlon
        if (kmt(i,j) .ge. 1) v2d(i,j,iv2d_sst) = REAL(buf4(i,j,ilev_sfc),r_size) - t0c !kelvin
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable SST"

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SST"
      WRITE(6,*) "read_diag:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sst) = ",MAXVAL(v2d(:,:,iv2d_sst))
    endif
! !STEVE: end

    !!! SSS
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'s_surf',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable sfc salt"
    DO j=1,nlat
      DO i=1,nlon
        v2d(i,j,iv2d_sss) = REAL(buf4(i,j,ilev_sfc),r_size)
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable SSS"

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SSS"
      WRITE(6,*) "read_diag:: sffile = ", sffile
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
      if (dodebug) WRITE(6,*) "read_diag:: just got data for variable sea_lev"
      DO j=1,nlat
        DO i=1,nlon
          !STEVE: Hopefully reading in meters here... (data might be in cm)
          v2d(i,j,iv2d_ssh) = REAL(buf4(i,j,ilev_sfc),r_size)
        END DO
      END DO
      if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable SSH"

! !STEVE: debug
      if (dodebug) then
        WRITE(6,*) "POST-SSH"
        WRITE(6,*) "read_diag:: sffile = ", sffile
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
    WRITE(6,*) "read_diag:: just opened file ", bfile

    !!! SSH
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'eta_t',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4(:,:,ilev_sfc)) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable eta_t"
    DO j=1,nlat
      DO i=1,nlon
        !STEVE: Hopefully reading in meters here... (data might be in cm)
        v2d(i,j,iv2d_eta) = REAL(buf4(i,j,ilev_sfc),r_size)
      END DO
    END DO
    if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable SSH"

    ! Convert SSH eta stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,iv2d_eta) = v2d(:,:,iv2d_eta) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_diag:: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
      WRITE(6,*) "min val for level v2d(:,:,iv2d_eta) = ", MINVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  else
    WRITE(6,*) "read_diag:: Skipping SFC eta from: ", bfile
  endif altimetry

  ! For additional variables:
  ! E.g. surface fluxes, drifters

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  RETURN
END SUBROUTINE read_diag


SUBROUTINE read_grd4(infile,v3d,v2d)
!===============================================================================
! Read in a set of netcdf-format mom4 restart files in single precision
!===============================================================================
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
  REAL(r_size), PARAMETER :: vmax = 1.0e18

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

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

END SUBROUTINE read_grd4


SUBROUTINE write_grd4(outfile,v3d_in,v2d_in)
!===============================================================================
! Write out a set of netcdf-format restart files in single precision
!===============================================================================
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

END SUBROUTINE write_grd4


SUBROUTINE read_bingrd(filename,v3d,v2d)
!===============================================================================
! Read in an letkf grd-format binary file
!===============================================================================
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
END SUBROUTINE read_bingrd


SUBROUTINE read_bingrd4(filename,v3d,v2d)
!===============================================================================
! Read in an letkf grd-format binary file in single precision
!===============================================================================
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


SUBROUTINE write_bingrd4(filename,v3d,v2d)
!===============================================================================
! Write out an letkf grd-format binary file in single precision
!===============================================================================
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


SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
!===============================================================================
! Compute the ensemble mean
!===============================================================================
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
!STEVE: all of these are still direct copies from mom2 version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE minkowski_flick(v3d0,v2d0,v3d,v2d)
!===============================================================================
! Extend the water points onver the land boundaries to do computations near land more easily
!===============================================================================
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
! if (dodebug) then
!   CALL write_bingrd('test_mink_v3d.grd',v3d,v2d)
!   CALL write_bingrd('test_mink_v3d0.grd',v3d0,v2d0)
! endif

  DEALLOCATE(buf)
  DEALLOCATE(buf3d)
  DEALLOCATE(buf2d)
  DEALLOCATE(vcnt)
  DEALLOCATE(mask)
  DEALLOCATE(kmt3d)
  DEALLOCATE(kmt2d)

  endif
END SUBROUTINE minkowski_flick


SUBROUTINE set_kmt(file)
!===============================================================================
! Set the kmt (depth of water in grid cells)
!===============================================================================
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


SUBROUTINE update_kmt(v3d,v2d)
!===============================================================================
! Update the kmt data
!===============================================================================
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


END MODULE common_mom4
