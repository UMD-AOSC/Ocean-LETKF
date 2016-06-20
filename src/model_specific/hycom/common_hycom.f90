MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for HYCOM
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   06/08/2015 Steve Penny converted for use with HYCOM
!
!=======================================================================
  USE common
  USE params_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS
  USE vars_model,   ONLY: lon, lat, lev, lon2d, lat2d, dx0, dy0, lon0, lonf, lat0, latf, wrapgap, dx, dy, kmt, kmt0, phi0
  USE hycom_io,     ONLY: read_blkdat
  USE params_obs,   ONLY: SSHclm_file

  IMPLICIT NONE

  PUBLIC :: set_common_oceanmodel, read_diag, read_restart, write_restart
  PUBLIC :: read_bingrd, write_bingrd, read_bingrd4, write_bingrd4
  PUBLIC :: ensmean_grd

  PRIVATE
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
!
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: height(nlon,nlat,nlev)        !(OCEAN)(MOM6)
  REAL(r_size),SAVE :: sample_t(nlon,nlat,nlev)      !(OCEAN)(HYCOM)
  REAL(r_sngl),SAVE :: depth(nlon,nlat)              !(OCEAN)(MOM6) !STEVE: are floats in netcdf single precision? I think so.
  REAL(r_sngl),SAVE :: wet(nlon,nlat)                !(OCEAN)(MOM6)
  REAL(r_size),SAVE :: area_t(nlon,nlat)             !(OCEAN)
  CHARACTER(4),SAVE :: element(nv3d+nv2d+nv4d)
  REAL(r_sngl),ALLOCATABLE :: sigma(:)               !(HYCOM) reference density levels, will be assigned to lev()
  !STEVE: for Custom Localization
  INTEGER :: nobids(nij0*nlev)                       !(OCEAN)
  !STEVE: for verifying against input netcdf file
  INTEGER :: nlon0=0, nlat0=0, nlev0=0               !(OCEAN)
  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: ncundef = 1.267650600228229E+030 !1.0e18
  REAL(r_size), DIMENSION(nlon,nlat) :: SSHclm_m

  INTEGER :: islot
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: xc,yc
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: sfc_data
  INTEGER, DIMENSION(:), ALLOCATABLE :: xi,yi

  !For input/output of model binary files (converted form HYCOM ab-format):
! CHARACTER(slen) :: base  = 'fsd.bin'! (HYCOM)(All variables assembled in one binary file)
  !Otherwise, use netcdf files:
  CHARACTER(slen) :: tbase = '3zt.nc' ! Temperature
  CHARACTER(slen) :: sbase = '3zs.nc' ! Salinity
  CHARACTER(slen) :: ubase = '3zu.nc' ! Zonal Current
  CHARACTER(slen) :: vbase = '3zv.nc' ! Merional Current
  CHARACTER(slen) :: hbase = '3dh.nc' ! Layer Thickness
  CHARACTER(slen) :: ebase = 'fsd.nc' ! Sea Surface Height
  CHARACTER(slen) :: drbase           ! (DRIFTERS)

  LOGICAL :: dodebug = .true.

  ! For temporary dx and dy computations:
  REAL(r_size) :: dlon, dlat, d1, d2, d3

CONTAINS

SUBROUTINE set_common_oceanmodel
!===============================================================================
! Initialize the module
!===============================================================================
  USE netcdf

  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex

  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_s) = 'S   '             !(OCEAN)
  element(iv3d_h) = 'H   '             !(OCEAN)(MOM6)(HYCOM)
  element(nv3d+iv2d_ssh) = 'SSH '      !(OCEAN)
  element(nv3d+iv2d_ubt) = 'UBT '      !(OCEAN)(HYCOM)
  element(nv3d+iv2d_vbt) = 'VBT '      !(OCEAN)(HYCOM)
! element(nv3d+iv2d_sst) = 'SST '      !(OCEAN)
! element(nv3d+iv2d_sss) = 'SSS '      !(OCEAN)
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
      WRITE(6,*) "Exiting common_hycom.f90..."
      STOP(1)
    ENDIF
  else
    SSHclm_m = 0.0d0
  endif

  !STEVE: needed to index the observations, because HYCOM does not have a regular grid above 45ºN
  !       Instead, HYCOM has a 2D-field for the representation of the lon/lat coordinate.
  !       This requires an update to the lon and lat arrays to lon2d and lat2d, with subsequent
  !       updates in letkf_local.f90 and letkf_obs.f90. These are primarily used to sort and store
  !       observation locations for later lookup at the localization phase.
  dx0 = 0.239990234375 !360.0/nlon  !0.23999
! lon(1) = 74.23999
! lon(1) = 74.12024 !434.11975
  lon(1) = 74.2399902343750 
  do i=2,nlon
    lon(i) = lon(i-1) + dx0
  enddo
! dy0 = (180.0)/nlat  !0.096
! lat(1) = -78.608
! dy0 = (78.6080017089844 + 89.93317)/nlat
  dy0 = (78.6080017089844 + 90.0)/nlat
! lat(1) = -78.60800
  lat(1) = -78.6080017089844
  do j=2,nlat
    lat(j) = lat(j-1) + dy0
  enddo

  !
  ! Lon, Lat, f, orography
  !
!STEVE: this part adapted from ROMS, update for MOM4 (and later MOM6) netcdf files:
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile 
    WRITE(6,*) "Exiting common_hycom.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile
  CALL check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )

  CALL check( NF90_INQ_VARID(ncid,'Longitude',varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lon2d) )
  WRITE(6,*) "lon2d(1,1) = ", lon2d(1,1)
  WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)
  CALL check( NF90_INQ_VARID(ncid,'Latitude',varid) )   ! Latitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lat2d) )
  WRITE(6,*) "lat2d(1,1) = ", lat2d(1,1)
  WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)

  !-----------------------------------------------------------------------------
  ! Get the reference density levels:
  !-----------------------------------------------------------------------------
! CALL check( NF90_INQ_VARID(ncid,'Depth',varid) )      ! depth of T-cell
! CALL check( NF90_GET_VAR(ncid,varid,lev) )
  CALL read_blkdat(gridfile1,sigma)
  lev = sigma
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !-----------------------------------------------------------------------------
  !(HYCOM) Sample temperature to make kmt field (2d land-sea depth mask)
  !-----------------------------------------------------------------------------
  CALL check( NF90_INQ_VARID(ncid,'temperature',varid) )      
  CALL check( NF90_GET_VAR(ncid,varid,sample_t) )
  WRITE(6,*) "sample_t(1,1,1) = ", sample_t(1,1,1)
  WRITE(6,*) "sample_t(nlon,nlat,nlev) = ", sample_t(nlon,nlat,nlev)

  if (.false. .and. dodebug) then
  WRITE(6,*) "========================================="
  WRITE(6,*) "STEVE: debugging for HYCOM input grid..."
  WRITE(6,*) "lon = ", lon
  WRITE(6,*) "lat = ", lat
  WRITE(6,*) "lon2d(:,1) = ", lon2d(:,1)
  WRITE(6,*) "lon2d(1,:) = ", lon2d(1,:)
  WRITE(6,*) "lat2d(:,1) = ", lat2d(:,1)
  WRITE(6,*) "lat2d(1,:) = ", lat2d(1,:)
  WRITE(6,*) "MAXVAL(lat2d) = ", MAXVAL(lat2d)
  WRITE(6,*) "done lon/lat."
  WRITE(6,*) "========================================="
  endif

  !
  ! dx and dy
  !
!!STEVE:MOM6: open new gridfile (ocean_hgrid.nc, gridfile3)
  INQUIRE(FILE=trim(gridfile3),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile3
    WRITE(6,*) "Exiting common_hycom.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile3
  CALL check( NF90_OPEN(gridfile3,NF90_NOWRITE,ncid3) )

  if (.false.) then !STEVE: for MOM6: this is TEMPORARY, until I figure out how to use the ocean_hgrid.nc data
    dx=0.0
    dy=0.0
    area_t=0.0
    do j=2,nlat-1
      dlat = (lat(j+1) - lat(j-1))/2.0
      d2 = (sin(dlat/2.0))**2 
      d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) )
      dy(:,j) = re * d3
 
      do i=2,nlon-1
        dlon = (lon(i+1) - lon(i-1))/2.0
        d2 = cos(lat(i-1)) * cos(lat(i+1)) * (sin(dlon/2.0))**2
        d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) ) 
        dx(i,j) = re * d3
        area_t(i,j) = dx(i,j)*dy(i,j)
      enddo
    enddo
  endif

  if (.true.) then !STEVE: for HYCOM: this is TEMPORARY, better if it's read in
    WRITE(6,*) "Computing dx and dy for area_t..."
    dx=0.0
    dy=0.0
    area_t=0.0
    do i=1,nlon
      ! South pole: (j==1)
      dlat = abs( lat2d(i,2) - lat2d(i,1) )
      d2 = (sin(dlat/2.0))**2 
      d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) )
      dy(i,1) = re * d3
      ! All middle latitudes:
      do j=2,nlat-2
        dlat = abs( lat2d(i,j+1) - lat2d(i,j-1) )/2.0
        d2 = (sin(dlat/2.0))**2 
        d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) )
        dy(i,j) = re * d3
      enddo
      ! North pole: (j==nlat)
      dlat = abs( lat2d(i,nlat) - lat2d(i,nlat-1) )
      d2 = (sin(dlat/2.0))**2 
      d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) )
      dy(i,nlat) = re * d3
      j=nlat

      if (dodebug .and. i==nlon .and. j==nlat) then
        WRITE(6,*) "common_hycom.f90:: debug the area_t computation for dy:"
        WRITE(6,*) "i = ", i
        WRITE(6,*) "j = ", j
        WRITE(6,*) "lat2d(i,nlat)   = ", lat2d(i,nlat)
        WRITE(6,*) "lat2d(i,nlat-1) = ", lat2d(i,nlat-1)
        WRITE(6,*) "dlat = ", dlat
        WRITE(6,*) "d2 = ", d2
        WRITE(6,*) "d3 = ", d3
        WRITE(6,*) "dy(i,j) = ", dy(i,j)
      endif

      ! Compute dx and area_t for ALL latitudes:
      do j=1,nlat
        dlon = modulo( lon2d(modulo(i+1-1,nlon)+1,j) - lon2d(modulo(i-1-1,nlon)+1,j), 360.0 )/2.0
        d2 = cos(lat2d(modulo(i-1-1,nlon)+1,j)) * cos(lat2d(modulo(i+1-1,nlon)+1,j)) * (sin(dlon/2.0))**2
        d3 = 2 * atan2( sqrt(d2), sqrt(1-d2) ) 
        dx(i,j) = re * d3
        area_t(i,j) = dx(i,j)*dy(i,j)
        if (dodebug .and. i==1 .and. j==1) then
          WRITE(6,*) "common_hycom.f90:: debug the area_t computation for dx:"
          WRITE(6,*) "i = ", i
          WRITE(6,*) "j = ", j
          WRITE(6,*) "modulo(i+1-1,nlon)+1 = ", modulo(i+1-1,nlon)+1
          WRITE(6,*) "modulo(i-1-1,nlon)+1 = ", modulo(i-1-1,nlon)+1
          WRITE(6,*) "dlon = ", dlon
          WRITE(6,*) "d2 = ", d2
          WRITE(6,*) "d3 = ", d3
          WRITE(6,*) "dx(1,1) = ", dx(i,j)
          WRITE(6,*) "dy(1,1) = ", dy(i,j)
          WRITE(6,*) "area_t(1,1) = ", area_t(i,j)
        endif
      enddo

    enddo
  endif

! WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile3
  WRITE(6,*) "dx(1,1) = ", dx(1,1)
  WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
  WRITE(6,*) "dy(1,1) = ", dy(1,1)
  WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)
  WRITE(6,*) "area_t(1,1) = ", area_t(1,1)
  WRITE(6,*) "area_t(nlon,nlat) = ", area_t(nlon,nlat)

  !
  ! kmt data
  !
  !STEVE:MOM6: sum(h>0) (from MOM.res.nc) to find ocean layer depth, where depth > 0 (from ocean_togog.nc)
  !STEVE:HYCOM: this needs to correspond to the reference model layers,
  !             for now, it does not since the layers have been interpolated
  !             from 32 layer to 40 depths in the netcdf files.

  INQUIRE(FILE=trim(gridfile2),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile2
    WRITE(6,*) "Exiting common_hycom.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile2
  CALL check( NF90_OPEN(trim(gridfile2),NF90_NOWRITE,ncid2) )

  CALL check( NF90_INQ_VARID(ncid2,'bathymetry',varid) ) ! number of vertical T-cells
  CALL check( NF90_GET_VAR(ncid2,varid,phi0) )
  WRITE(6,*) "phi0(1,1) = ", phi0(1,1)
  WRITE(6,*) "phi0(nlon,nlat) = ", phi0(nlon,nlat)

  ! Find the 'kmt' value for the depth of the ocean grid cells (phi0)
  ! sum the layer thicknesses where depth > 0
  kmt=0
  do k=nlev,1,-1
    do j=1,nlat
      do i=1,nlon
!       if (height(i,j,k) > 0 .and. wet(i,j) > 0 .and. kmt(i,j) == 0 ) then
        if (kmt(i,j) == 0 .and. sample_t(i,j,k) < ncundef-2 ) then
!         phi0(i,j) = sum(height(i,j,1:k))
          kmt(i,j) = k
        endif
      enddo
    enddo
  enddo
  kmt0 = REAL(kmt,r_size)

  !
  ! Corioris parameter
  !
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)

  ! Close the grid files:
  CALL check( NF90_CLOSE(ncid) )
!!CALL check( NF90_CLOSE(ncid1) )
  CALL check( NF90_CLOSE(ncid2) )
  CALL check( NF90_CLOSE(ncid3) )

  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  !STEVE: for MOM4p1:
  wrapgap = 360.0d0 - abs(modulo(lon0,360.0)) - abs(modulo(lonf,360.0))

  if (dodebug) then
    WRITE(6,*) "common_hycom.f90::set_common_hycom:: lon0,lonf,lat0,latf = ", lon0,lonf,lat0,latf
  endif

  RETURN
END SUBROUTINE set_common_oceanmodel

!STEVE: ISSUE: need to update for HYCOM
SUBROUTINE read_etaclm(SSHclm_file,SSHclm_m)
!===============================================================================
! Read in (e.g.) 1993-1999 model climatology to subtract as done by AVISO
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
  CALL check( NF90_OPEN(SSHclm_file,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file

  buf4=0.0
  CALL check( NF90_INQ_VARID(ncid,'ssh',varid) )
  CALL check( NF90_GET_VAR(ncid,varid,buf4) )
  DO j=1,nlat
    DO i=1,nlon
      !STEVE: Hopefully reading in meters here... (data might be in cm)
      SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
    END DO
  END DO

  CALL check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm

!STEVE: add this:
SUBROUTINE read_diag(infile,v3d,v2d,prec)
!===============================================================================
! Read a archive/restart file in binary format extracted from hycom ab-format
!===============================================================================
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by HYCOM
  !       The t,s,u,v,(ssh) fields are read in from the ab archive files so that the o-f data can be computed
  !       by the obsop.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  USE netcdf
  USE hycom_io, ONLY: read_hycom, hycom_undef
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1=single, 2=double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4 !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8 !(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .true.
  CHARACTER(slen) :: varname, binfile

  if (prec == 1) then
    WRITE(6,*) "read_diag::  single precision requested, but option not available at this time."
    WRITE(6,*) "read_diag::  using double precision read..."
  elseif (prec == 2) then
    WRITE(6,*) "read_diag::  using double precision read..."
  else
    WRITE(6,*) "read_diag:: unsupported option, prec = ", prec
    STOP
  endif

  ! STEVE: this is provided externally at the moment
  binfile = trim(infile)//trim(base)
  if (dodebug) WRITE(6,*) "read_diag:: Pre-read_hycom: Reading file: ", trim(binfile)
  CALL read_hycom(binfile,v3d,v2d)
  if (dodebug) WRITE(6,*) "read_diag:: Post-read_hycom: Just read file: ", trim(binfile)

  if (dodebug) then
    WRITE(6,*) "read_diag:: Post-read_hycom"
    buf8 = v3d(:,:,:,iv3d_t)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_t)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_t)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_s)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_s)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_s)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_u)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_u)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_u)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_v)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_v)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_v)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_h)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_h)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_h)) = ", MINVAL(buf8(:,:,1))
    WRITE(6,*)
  endif

END SUBROUTINE read_diag

!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the HYCOM restart file.
  !       The t,s,u,v,(ssh) fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsope.f90 program prior
  !       to running letkf.f90, so the archive (diagnostic) files should not have to be touched during letkf runtime.
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1 = single, 2 = double
  REAL(r_size), ALLOCATABLE :: v3d8(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d8(:,:,:)

! !STEVE: for now, only use the HYCOM archive format. The output ab file will be converted to a restart file externally.
  ALLOCATE(v3d8(nlon,nlat,nlev,nv3d),v2d8(nlon,nlat,nv2d))
  CALL read_diag(infile,v3d8,v2d8,prec)
  v3d = REAL(v3d8,r_sngl)
  v2d = REAL(v2d8,r_sngl)
  DEALLOCATE(v3d8,v2d8)

END SUBROUTINE read_restart

!-----------------------------------------------------------------------
subroutine check(status)
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

!-----------------------------------------------------------------------
! Write a set of MOM6 restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d,v2d)
  !STEVE: This writes out the analysis to a pre-existing template netcdf file.
  !       IN THE FUTURE: output directly as a HYCOM restart file
  USE netcdf
  USE hycom_io, ONLY: write_hycom, hycom_undef
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  LOGICAL, PARAMETER :: do_physlimit=.true.
  CHARACTER(slen) :: binfile
  CHARACTER(3) :: MEM3

  ! STEVE: this is provided externally at the moment
  binfile = trim(outfile)//trim(base)
  WRITE(6,*) "common_hycom.f90::write_restart:: writing file: ", trim(binfile)

  WRITE(6,*) "common_hycom.f90::write_restart:: calling write_hycom..."
  CALL write_hycom(binfile,v3d,v2d)
  WRITE(6,*) "common_hycom.f90::write_restart:: Finished calling write_hycom."

END SUBROUTINE write_restart

!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_bingrd(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
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

!-----------------------------------------------------------------------
! Read a grid file with single precision
!-----------------------------------------------------------------------
SUBROUTINE read_bingrd4(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
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

!-----------------------------------------------------------------------
!-- Write a grid file -------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE write_bingrd(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:) !(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec
  LOGICAL, PARAMETER :: dodebug=.false.

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

!-----------------------------------------------------------------------
! Write a grid file in single precision
!-----------------------------------------------------------------------
SUBROUTINE write_bingrd4(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL, PARAMETER :: dodebug=.false.

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
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
  use hycom_io, ONLY: hycom_undef
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: cnt3d(nij,nlev,nv3d)
  INTEGER :: cnt2d(nij,nv2d)
  INTEGER :: i,k,m,n

! ALLOCATE(v3dm(nij,nlev,nv3d),v2dm(nij,nv2d))


  !STEVE: For HYCOM, we have to check whether the member point
  !       is undefined, because some members can differ in the
  !       grid definition.

  cnt3d = 0
  v3dm  = 0.0d0
  do n=1,nv3d
    do k=1,nlev
      do i=1,nij
        do m=1,member
          if (v3d(i,k,m,n) < hycom_undef) then
            v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
            cnt3d(i,k,n) = cnt3d(i,k,n) + 1
          endif
        enddo
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(cnt3d(i,k,n),r_size)
      enddo
    enddo
  enddo

  cnt2d = 0
  v2dm = 0.0d0
  do n=1,nv2d
    do i=1,nij
      do m=1,member
        if (v2d(i,m,n) < hycom_undef) then
          v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
          cnt2d(i,n) = cnt2d(i,n) + 1
        endif
      enddo
      v2dm(i,n) = v2dm(i,n) / REAL(cnt2d(i,n),r_size)
    enddo
  enddo

END SUBROUTINE ensmean_grd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!STEVE: additional subroutines for OCEAN
!STEVE: all of these are still direct copies from mom2 version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
! Adjustment available for near-coastline observations
!-----------------------------------------------------------------------
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


SUBROUTINE set_kmt(file)
!===============================================================================
! Set the kmt field that describes the water depth in grid cells
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
! Update the kmt field that describes the water depth in grid cells
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
  WRITE(6,*) "common_hycom.f90::update_kmt: KMT entries updated: ", nkmterr

END SUBROUTINE update_kmt

END MODULE common_oceanmodel

