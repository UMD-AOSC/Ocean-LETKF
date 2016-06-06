MODULE common_mom6
!=======================================================================
!
! [PURPOSE:] Common Information for MOM6
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!
!=======================================================================
  USE common
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
! GFDL MOM6:
  INTEGER,PARAMETER :: nlon=1440
  INTEGER,PARAMETER :: nlat=1080
  INTEGER,PARAMETER :: nlev=75
! MOM4 ncep2012 tripolar
! INTEGER,PARAMETER :: nlon=720 
! INTEGER,PARAMETER :: nlat=410 
! INTEGER,PARAMETER :: nlev=40 
!
  INTEGER,PARAMETER :: ilev_sfc=1
!
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s              !(OCEAN)
! INTEGER,PARAMETER :: nv3d=5 ! u,v,t,s,h            !(OCEAN)(MOM6)
  INTEGER,PARAMETER :: nv4d=3 ! x,y,z                !(OCEAN) STEVE: add t,x,y,z,id for DRIFTERS
! INTEGER,PARAMETER :: nv2d=3 ! ssh,sst,sss          !(OCEAN)
! INTEGER,PARAMETER :: nv2d=7 ! ssh/t/s, + sfc fluxes: taux,tauy,heat,freshwater
  INTEGER,PARAMETER :: nv2d=1 !3!4 ! ssh,sst,sss,eta      !(OCEAN) !(ALTIMETRY)
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4                      !(OCEAN)
! INTEGER,PARAMETER :: iv3d_h=5                      !(OCEAN)(MOM6)
                                                     !          From ocean_sbc.res.nc:
  INTEGER,PARAMETER :: iv2d_ssh=1                    !(OCEAN) ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0)
  INTEGER,PARAMETER :: iv2d_sst=2                    !(OCEAN) ! time averaged sst (Kelvin) passed to atmosphere/ice model
  INTEGER,PARAMETER :: iv2d_sss=3                    !(OCEAN) ! time averaged sss (psu) passed to atmosphere/ice models
  INTEGER,PARAMETER :: iv2d_eta=4                    !(OCEAN) ! eta sea surface perturbation from mom6's ocean_barotropic.res.nc restart file
  INTEGER,PARAMETER :: iv4d_x=1                      !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_y=2                      !(OCEAN) (DRIFTERS)
  INTEGER,PARAMETER :: iv4d_z=3                      !(OCEAN) (DRIFTERS)
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
  REAL(r_size),SAVE :: height(nlon,nlat,nlev)        !(OCEAN)(MOM6)
  REAL(r_size),SAVE :: kmt0(nlon,nlat)               !(OCEAN)
  REAL(r_sngl),SAVE :: depth(nlon,nlat)              !(OCEAN)(MOM6) !STEVE: are floats in netcdf single precision? I think so.
  REAL(r_sngl),SAVE :: wet(nlon,nlat)                !(OCEAN)(MOM6)
  REAL(r_size),SAVE :: area_t(nlon,nlat)             !(OCEAN)
  CHARACTER(4),SAVE :: element(nv3d+nv2d+nv4d)
  INTEGER, DIMENSION(nlon,nlat), SAVE     :: kmt=-1  !(OCEAN) STEVE: the bottom topography for mom6
  !STEVE: for Custom Localization
  INTEGER :: nobids(nij0*nlev)                       !(OCEAN)
  !STEVE: for verifying against input netcdf file
  INTEGER :: nlon0=0, nlat0=0, nlev0=0               !(OCEAN)
  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: vmax = 1.0e18
  !STEVE: For generalized grid
  REAL(r_size) :: lon0, lonf, lat0, latf
  REAL(r_size) :: wrapgap
  CHARACTER(14) :: SSHclm_file = 'aEtaCds9399.nc'
  REAL(r_size), DIMENSION(nlon,nlat) :: SSHclm_m
! REAL(r_sngl) :: buf4_2d(nlon,nlat)
  ! For grid_spec.nc data file:
  CHARACTER(slen) :: gridfile  = 'MOM.res.nc'
  CHARACTER(slen) :: gridfile1 = 'MOM.res_1.nc'
  CHARACTER(slen) :: gridfile2 = 'ocean_topog.nc'
  CHARACTER(slen) :: gridfile3 = 'ocean_hgrid.nc'
! CHARACTER(29) :: diagfile = 'ocean_hourly_YYYY_MM_DD_HH.nc'

! For AMOC computation
  REAL(r_size) :: zb(nlev)
  REAL(r_size) :: dz(nlev)

! INTEGER, PARAMETER :: nslot = 6
  INTEGER :: islot
  LOGICAL :: coeff_file_exists
  INTEGER, SAVE :: nslon=-1,nslat=-1
  INTEGER :: nmlon=nlon,nmlat=nlat
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: slon, slat, mlon, mlat
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: xc,yc
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: sfc_data
  INTEGER, DIMENSION(:), ALLOCATABLE :: xi,yi

  !STEVE: for debugging
  LOGICAL, PARAMETER :: dodebug = .false.

  !For input/output model files:
  CHARACTER(slen) :: tsbase = 'MOM.res.nc' !(and u, and h)
  CHARACTER(slen) :: uvbase = 'MOM.res_1.nc' !(v and ave_ssh/sfc)
  CHARACTER(slen) :: hbase
  CHARACTER(slen) :: drbase

  ! For temporary dx and dy computations:
  REAL(r_size) :: dlon, dlat, d1, d2, d3

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_mom6
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex

  WRITE(6,'(A)') 'Hello from set_common_mom6'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_s) = 'S   '             !(OCEAN)
! element(iv3d_h) = 'h   '             !(OCEAN)(MOM6)
  element(nv3d+iv2d_ssh) = 'SSH '      !(OCEAN)
  element(nv3d+iv2d_sst) = 'SST '      !(OCEAN)
  element(nv3d+iv2d_sss) = 'SSS '      !(OCEAN)
  if (DO_ALTIMETRY) then
    element(nv3d+iv2d_eta) = 'eta '      !(OCEAN)
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
      WRITE(6,*) "Exiting common_mom6.f90..."
      STOP(1)
    ENDIF
  endif

  !
  ! Lon, Lat, f, orography
  !
!STEVE: this part adapted from ROMS, update for MOM4 (and later MOM6) netcdf files:
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile 
    WRITE(6,*) "Exiting common_mom6.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile
  call check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )

  call check( NF90_INQ_VARID(ncid,'lonh',varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon) )
  WRITE(6,*) "lon(1) = ", lon(1)
  WRITE(6,*) "lon(nlon) = ", lon(nlon)
  call check( NF90_INQ_VARID(ncid,'lath',varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat) )
  WRITE(6,*) "lat(1) = ", lat(1)
  WRITE(6,*) "lat(nlat) = ", lat(nlat)
  call check( NF90_INQ_VARID(ncid,'Layer',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !
  ! dx and dy
  !
!!STEVE:MOM6: open new gridfile (ocean_hgrid.nc, gridfile3)
! INQUIRE(FILE=trim(gridfile3),EXIST=ex)
! IF(.not. ex) THEN
!   WRITE(6,*) "The file does not exist: ", gridfile3
!   WRITE(6,*) "Exiting common_mom6.f90..."
!   STOP(2)
! ENDIF
! WRITE(6,'(A)') '  >> accessing file: ', gridfile3
! call check( NF90_OPEN(gridfile3,NF90_NOWRITE,ncid3) )

! call check( NF90_INQ_VARID(ncid3,'dx',varid) )    ! width of T_cell (meters)
! call check( NF90_GET_VAR(ncid3,varid,dx) ) 
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MIN(dx) = ", MINVAL(dx)
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MAX(dx) = ", MAXVAL(dx)
! call check( NF90_INQ_VARID(ncid3,'dy',varid) )    ! height of T_cell (meters)
! call check( NF90_GET_VAR(ncid3,varid,dy) ) 
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MIN(dy) = ", MINVAL(dy)
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MAX(dy) = ", MAXVAL(dy)
! call check( NF90_INQ_VARID(ncid3,'area',varid) )        ! area of T_cell
! call check( NF90_GET_VAR(ncid3,varid,area_t) ) 
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MIN(area_t) = ", MINVAL(area_t)
! WRITE(6,*) "common_mom6:: ", trim(gridfile3), " MAX(area_t) = ", MAXVAL(area_t)

  if (.true.) then !STEVE: this is TEMPORARY, until I figure out how to use the ocean_hgrid.nc data
    i=0
    j=0
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

! WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile3
  WRITE(6,*) "dx(1,1) = ", dx(1,1)
  WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
  WRITE(6,*) "dy(1,1) = ", dy(1,1)
  WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)

  !
  ! kmt data
  !
  !STEVE:MOM6: sum(h>0) (from MOM.res.nc) to find ocean layer depth, where depth > 0 (from ocean_togog.nc)
  call check( NF90_INQ_VARID(ncid,'h',varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,height) )
  WRITE(6,*) "h(1,1,1) = ", height(1,1,1)
  WRITE(6,*) "h(nlon,nlat,nlev) = ", height(nlon,nlat,nlev)

  INQUIRE(FILE=trim(gridfile2),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile2
    WRITE(6,*) "Exiting common_mom6.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile2
  call check( NF90_OPEN(trim(gridfile2),NF90_NOWRITE,ncid2) )

  call check( NF90_INQ_VARID(ncid2,'depth',varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid2,varid,depth) )
  WRITE(6,*) "depth(1,1) = ", depth(1,1)
  WRITE(6,*) "depth(nlon,nlat) = ", depth(nlon,nlat)

! !STEVE:MOM6: from ocean_topog.nc:
  call check( NF90_INQ_VARID(ncid2,'wet',varid) )        ! land/sea flag (0=land) for T-cell
  call check( NF90_GET_VAR(ncid2,varid,wet) )
  WRITE(6,*) "wet(1,1) = ", wet(1,1)
  WRITE(6,*) "wet(nlon,nlat) = ", wet(nlon,nlat)

  ! Find the 'kmt' value for the depth of the ocean grid cells
  ! sum the layer thicknesses where depth > 0
! kmt0=SUM(h,MASK=depth>0,DIM=3)
  kmt=0
  do k=nlev,1,-1
    do j=1,nlat
      do i=1,nlon
        if (height(i,j,k) > 0 .and. wet(i,j) > 0 .and. kmt(i,j) == 0 ) then
          phi0(i,j) = sum(height(i,j,1:k))
          kmt(i,j) = k
        endif
      enddo
    enddo
  enddo
  kmt0 = REAL(kmt,r_size)

  !STEVE: needed for computing the AMOC based on the streamfunction calculation: (in MOM6, thickness from mom restart MOM.res.nc)
  !STEVE:(MOM6) This won't work in MOM6...
! call check( NF90_INQ_VARID(ncid,'h',varid) )      ! depth of T-cell
! call check( NF90_GET_VAR(ncid,varid,zb) )
! WRITE(6,*) "zb(1) = ", zb(1)
! WRITE(6,*) "zb(nlev) = ", zb(nlev)

  ! Compute dz:
! dz(1) = zb(1)
! do k=2,nlev
!   dz(k) = zb(k)-zb(k-1)
! enddo

  !
  ! Corioris parameter
  !
!$OMP PARALLEL WORKSHARE
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
!$OMP END PARALLEL WORKSHARE

  ! Close the grid files:
  call check( NF90_CLOSE(ncid) )
!!call check( NF90_CLOSE(ncid1) )
  call check( NF90_CLOSE(ncid2) )
! call check( NF90_CLOSE(ncid3) )

  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

  RETURN
END SUBROUTINE set_common_mom6

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
!-- Read a diagnostic file in mom6 netcdf format ---------------------------------------------------
SUBROUTINE read_diag(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by MOM6
  !       The t,s,u,v,(ssh) fields are read in from the files so that the o-f data can be computed
  !       by the obsope.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  USE netcdf
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
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  ! Assumed format of MOM6 diag file, specified in diag_table:
  !
  ! "ocean_hourly%4yr%2mo%2dy%2hr",1,"hours",1,"hours","time",1,"hours"
  !
  !'ocean_model','temp','temp','ocean_hourly%4yr%2mo%2dy%2hr','all','none','none',2
  !'ocean_model','salt','salt','ocean_hourly%4yr%2mo%2dy%2hr','all','none','none',2
  !'ocean_model','u','u','ocean_hourly%4yr%2mo%2dy%2hr','all','none','none',2
  !'ocean_model','v','v','ocean_hourly%4yr%2mo%2dy%2hr','all','none','none',2
  !'ocean_model','ssh','ssh','ocean_hourly%4yr%2mo%2dy%2hr','all','none','none',2
  !
  ! Example file:
  ! netcdf ocean_hourly_1900_01_02_10.nc {
  !  dimensions:
  !  	xh = 720 ;
  !  	yh = 540 ;
  !  	zl = 75 ;
  !  	time = UNLIMITED ; // (1 currently)
  !  	xq = 720 ;
  !  	yq = 540 ;
  !  variables:
  !  	double xh(xh) ;
  !  		xh:long_name = "h point nominal longitude" ;
  !  		xh:units = "degrees_E" ;
  !  		xh:cartesian_axis = "X" ;
  !  		xh:domain_decomposition = 1, 1440, 1, 720 ;
  !  	double yh(yh) ;
  !  		yh:long_name = "h point nominal latitude" ;
  !  		yh:units = "degrees_N" ;
  !  		yh:cartesian_axis = "Y" ;
  !  		yh:domain_decomposition = 1, 1080, 1, 540 ;
  !  	double zl(zl) ;
  !  		zl:long_name = "Layer pseaduo-depth, -z*" ;
  !  		zl:units = "meter" ;
  !  		zl:cartesian_axis = "Z" ;
  !  		zl:positive = "down" ;
  !  	double time(time) ;
  !  		time:long_name = "time" ;
  !  		time:units = "hours since 1900-01-01 00:00:00" ;
  !  		time:cartesian_axis = "T" ;
  !  		time:calendar_type = "NOLEAP" ;
  !  		time:calendar = "NOLEAP" ;
  !  	double xq(xq) ;
  !  		xq:long_name = "q point nominal longitude" ;
  !  		xq:units = "degrees_E" ;
  !  		xq:cartesian_axis = "X" ;
  !  		xq:domain_decomposition = 1, 1440, 1, 720 ;
  !  	double yq(yq) ;
  !  		yq:long_name = "q point nominal latitude" ;
  !  		yq:units = "degrees_N" ;
  !  		yq:cartesian_axis = "Y" ;
  !  		yq:domain_decomposition = 1, 1080, 1, 540 ;
  !  	float temp(time, zl, yh, xh) ;
  !  		temp:long_name = "Potential Temperature" ;
  !  		temp:units = "Celsius" ;
  !  		temp:missing_value = -1.e+34f ;
  !  		temp:_FillValue = -1.e+34f ;
  !  		temp:cell_methods = "time: point" ;
  !  	float salt(time, zl, yh, xh) ;
  !  		salt:long_name = "Salinity" ;
  !  		salt:units = "PSU" ;
  !  		salt:missing_value = -1.e+34f ;
  !  		salt:_FillValue = -1.e+34f ;
  !  		salt:cell_methods = "time: point" ;
  !  	float u(time, zl, yh, xq) ;
  !  		u:long_name = "Zonal velocity" ;
  !  		u:units = "meter  second-1" ;
  !  		u:missing_value = -1.e+34f ;
  !  		u:_FillValue = -1.e+34f ;
  !  		u:cell_methods = "time: point" ;
  !  	float v(time, zl, yq, xh) ;
  !  		v:long_name = "Meridional velocity" ;
  !  		v:units = "meter second-1" ;
  !  		v:missing_value = -1.e+34f ;
  !  		v:_FillValue = -1.e+34f ;
  !  		v:cell_methods = "time: point" ;
  !  	float ssh(time, yh, xh) ;
  !  		ssh:long_name = "Sea Surface Height" ;
  !  		ssh:units = "meter" ;
  !  		ssh:missing_value = -1.e+34f ;
  !  		ssh:_FillValue = -1.e+34f ;
  !  		ssh:cell_methods = "time: point" ;
  !  
  !  // global attributes:
  !  		:filename = "ocean_hourly_1900_01_02_10.nc.0000" ;
  !  		:NumFilesInSet = 4 ;
  !  		:title = "MOM_SIS_025_z" ;
  !  		:grid_type = "regular" ;
  !  		:grid_tile = "N/A" ;
  !  }

  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the T/S netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_diag :: just opened file ", infile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! t
  varname='temp'
  ivid=iv3d_t

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )

  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select

  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ",MAXVAL(v3d(:,:,k,iv3d_t))
    enddo
  endif
  ! !STEVE: end

  !!! s
  varname='salt'
  ivid=iv3d_s

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_s) = ",MAXVAL(v3d(:,:,k,iv3d_s))
    enddo
  endif
  ! !STEVE: end

  !!! u
  varname='u'
  ivid=iv3d_u

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ",MAXVAL(v3d(:,:,k,iv3d_u))
    enddo
  endif
  ! !STEVE: end

  !!! v
  varname='v'
  ivid=iv3d_v

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_v) = ",MAXVAL(v3d(:,:,k,iv3d_v))
    enddo
  endif
  ! !STEVE: end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! ssh
  varname='ssh'
  ivid=iv2d_ssh

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
      v2d(:,:,ivid) = REAL(buf4(:,:,1),r_size) - SSHclm_m(:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
      v2d(:,:,ivid) = buf8(:,:,1) - SSHclm_m(:,:)
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)

  call check( NF90_CLOSE(ncid) )
  
END SUBROUTINE read_diag

!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the MOM.res.nc, MOM.res_1.nc, etc. restart files produced by MOM6
  !       The t,s,u,v,(ssh) fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsope.f90 program prior
  !       to running letkf.f90, so the diagnostic files should not have to be touched during letkf runtime.
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1 = single, 2 = double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  CHARACTER(slen) :: tsfile,uvfile, sffile, bfile, drfile ! (TS) (UV) (SFC) (barotropic - eta) (DRIFTERS)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  !STEVE:
  REAL(r_size) :: meanSSH !STEVE: for temporary SSH estimate based on heat content
  REAL(r_size) :: videpth !STEVE: depth of vertically integrated heat content
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
  !LOGICAL, PARAMETER :: dodebug = .true.
  CHARACTER(3) :: MEM3
  CHARACTER(32) :: sfc_infile
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  tsfile = trim(infile)//'.'//trim(tsbase)
  uvfile = trim(infile)//'.'//trim(uvbase)
  bfile  = trim(infile)//'.'//trim(hbase)

  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

!! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the T/S netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart:: just opened file ", tsfile

  !!! t
  varname='Temp'
  ivid=iv3d_t

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable temp"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable temp"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_restart :: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ",MAXVAL(v3d(:,:,k,iv3d_t))
    enddo
  endif
! !STEVE: end

  !!! s
  varname='Salt'
  ivid=iv3d_s

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable salt"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable salt"

 !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_restart :: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_s) = ", MAXVAL(v3d(:,:,k,iv3d_s))
    enddo 
  endif
! !STEVE: end

  !!! h (MOM6)
! buf4=0.0
! call check( NF90_INQ_VARID(ncid,'h',varid) )
! call check( NF90_GET_VAR(ncid,varid,buf4) )
! if (dodebug) WRITE(6,*) "read_restart :: just got data for variable: layer thickness (h)"
! DO k=1,nlev
!   DO j=1,nlat
!     DO i=1,nlon
!       v3d(i,j,k,iv3d_h) = REAL(buf4(i,j,k),r_size)
!     END DO
!   END DO
! END DO
! if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable: layer thickness"

! !STEVE: debug
! if (dodebug) then
!   WRITE(6,*) "POST-h"
!   WRITE(6,*) "read_restart :: tsfile = ", tsfile
!   do k=1,nlev
!     WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_h) = ", MAXVAL(v3d(:,:,k,iv3d_h))
!   enddo 
! endif
! !STEVE: end

  !!! u
  varname='u'
  ivid=iv3d_u

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable u"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable u"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_restart :: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ",MAXVAL(v3d(:,:,k,iv3d_u))
    enddo
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the U/V netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )
  IF(istat /= NF90_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR in read_restart for ',uvfile
    STOP(7)
  END IF
  WRITE(6,*) "read_restart :: just opened file ", uvfile

  !!! v
  varname='v'
  ivid=iv3d_v

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

 if (dodebug) WRITE(6,*) "read_restart :: just got data for variable v"
 if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable v"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_restart :: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_v) = ", MAXVAL(v3d(:,:,k,iv3d_v))
    enddo 
  endif
! !STEVE: end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Use the forecast time-average SSH
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  altimetry : if(DO_ALTIMETRY) then

    !!! SSH
    varname='ave_ssh'
    ivid=iv2d_ssh

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    select case(prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
        v2d(:,:,ivid) = buf4(:,:,1) - REAL(SSHclm_m(:,:),r_sngl)
      case(2)
        buf8=0.0d0
        call check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
        v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl) - REAL(SSHclm_m(:,:),r_sngl)
    end select

    if (dodebug) WRITE(6,*) "read_restart :: just got data for variable ave_ssh"
    if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable SSH"

    ! Convert SSH eta stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,ivid) = v2d(:,:,iv2d_eta) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_restart :: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
      WRITE(6,*) "min val for level v2d(:,:,iv2d_eta) = ", MINVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

  else
    WRITE(6,*) "read_restart :: Skipping SFC eta from: ", bfile
  endif altimetry

  call check( NF90_CLOSE(ncid) )

  ! For additional variables:
  ! E.g. surface fluxes, drifters

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !STEVE: debug test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!
  RETURN
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
SUBROUTINE write_restart(outfile,v3d_in,v2d_in)
  USE netcdf
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:) !(nlon,nlat,nv2d)
  REAL(r_size), ALLOCATABLE :: t3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY)
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  LOGICAL, PARAMETER :: do_physlimit=.true.
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: mlon, mlat
  REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
  CHARACTER(slen) :: sfc_outfile, sfc_infile
  CHARACTER(3) :: MEM3

  !STEVE: this is the routine that writes out the individual analysis files for
  !       each esnsemble member in netcdf format.

  tsfile = trim(outfile)//'.'//trim(tsbase)
  uvfile = trim(outfile)//'.'//trim(uvbase)
  bfile  = trim(outfile)//'.'//trim(hbase)

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  v3d = REAL(v3d_in,r_size)
  v2d = REAL(v2d_in,r_size)

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

        if (v3d(i,j,k,iv3d_t) > 40.0) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = 40.0
        endif

        if (v3d(i,j,k,iv3d_s) < 0 ) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = 0.0
        endif

        if (v3d(i,j,k,iv3d_s) > 50.0) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = 50.0
        endif

        if (nv2d > 1 .and. k .eq. 1) then
          if (v2d(i,j,iv2d_sst) < -4) v2d(i,j,iv2d_sst) = -4.0
          if (v2d(i,j,iv2d_sst) > 40.0) v2d(i,j,iv2d_sst) = 40.0
          if (v2d(i,j,iv2d_sss) < 0) v2d(i,j,iv2d_sss) = 0.0
          if (v2d(i,j,iv2d_sss) > 50.0) v2d(i,j,iv2d_sss) = 50.0
        endif

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
!         WRITE(6,*) "common_mom6.f90::write_grd4:: ERROR: found NaN..."
!         WRITE(6,*) "v3d(i,j,k,m) contains NaN. i,j,k,m = ", i,j,k,m
!         STOP 1
!       endif
        enddo
      enddo
    enddo
  enddo
  endif

  call check( NF90_INQ_VARID(ncid,'Temp',varid) )

  !STEVE: debug DEBUG
  !STEVE: switch out the data to see if this writes properly
  !ALLOCATE(t3d(nlon,nlat,nlev,nv3d))
  !t3d(:,:,:,iv3d_t) = 1.0
  !call check( NF90_PUT_VAR(ncid,varid,t3d(:,:,:,iv3d_t)) )
  !DEALLOCATE(t3d)
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_t)) )

  !!! s
  call check( NF90_INQ_VARID(ncid,'Salt',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_s)) )

! !!! h (MOM6)
! call check( NF90_INQ_VARID(ncid,'h',varid) )
! call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_h)) )

  !!! u
  call check( NF90_INQ_VARID(ncid,'u',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_u)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! uv file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )

  !!! v
  call check( NF90_INQ_VARID(ncid,'v',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_v)) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the updated eta_t to analysis restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (DO_ALTIMETRY) then
    call check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    call check( NF90_INQ_VARID(ncid,'ave_ssh',varid) )

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,iv2d_eta) = v2d(:,:,iv2d_eta) + SSHclm_m(:,:)

    call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,iv2d_eta)) )
  endif

  call check( NF90_CLOSE(ncid) )

  RETURN
END SUBROUTINE write_restart

!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
! Read a grid file with single precision
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!-- Write a grid file -------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE write_bingrd(filename,v3d,v2d)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  WRITE(6,*) "common_mom6::update_kmt: KMT entries updated: ", nkmterr

  return
END SUBROUTINE update_kmt

END MODULE common_mom6
