MODULE common_oceanmodel
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

  USE params_letkf, ONLY: DO_ALTIMETRY, DO_SLA, DO_DRIFTERS

  IMPLICIT NONE

  PUBLIC

CONTAINS


!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_oceanmodel
  USE netcdf
  USE params_model, ONLY: SSHclm_file
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nlon2d, nlat2d
  USE vars_model,   ONLY: SSHclm_m, lon, lat, lon2d, lat2d, lev
  USE vars_model,   ONLY: dx, dy, area_t
  USE vars_model,   ONLY: height, kmt, depth, wet, phi0, kmt0
  USE vars_model,   ONLY: set_vars_model
  USE params_model, ONLY: initialize_params_model
  USE vars_model,   ONLY: initialize_vars_model
  USE params_model, ONLY: grid_lon_name, grid_lat_name, grid_lev_name
  USE params_model, ONLY: grid_lon2d_name, grid_lat2d_name
  USE params_model, ONLY: grid_wet_name, grid_depth_name, grid_height_name
  USE params_model, ONLY: grid_h_name
  USE params_model, ONLY: gridfile, gridfile2, gridfile3
  USE params_model, ONLY: grid_nlon_name, grid_nlat_name, grid_nlev_name

  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex
  REAL(r_size) :: dlat, dlon, d2, d3
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: work2d

  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'
  CALL initialize_params_model ! (checks to make sure it is initialized)

  !-----------------------------------------------------------------------------
  ! Open grid_spec.nc grid specification file for MOM6
  !-----------------------------------------------------------------------------
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", trim(gridfile)
    WRITE(6,*) "Exiting common_mom6.f90..."
    STOP (2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', trim(gridfile)
  call check( NF90_OPEN(trim(gridfile),NF90_NOWRITE,ncid) )

  !-----------------------------------------------------------------------------
  ! Read grid dimensions
  !-----------------------------------------------------------------------------
#ifdef DYNAMIC
  WRITE(6,*) "Reading x dimension (lon)..."
  call check( NF90_INQ_DIMID(ncid,trim(grid_nlon_name),dimid) )   ! Longitude for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlon) )
  WRITE(6,*) "nlon = ", nlon

  WRITE(6,*) "Reading y dimension (lat)..."
  call check( NF90_INQ_DIMID(ncid,trim(grid_nlat_name),dimid) )   ! Latitude for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlat) )
  WRITE(6,*) "nlat = ", nlat

  WRITE(6,*) "Reading z dimension (lev)..."
  call check( NF90_INQ_DIMID(ncid,trim(grid_nlev_name),dimid) )   ! Level for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlev) )
  WRITE(6,*) "nlev = ", nlev
#endif

  !-----------------------------------------------------------------------------
  ! Initialize the model variables (and allocate arrays if using #ifdef DYNAMIC)
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Initializing variables and allocating arrays..."
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  !-----------------------------------------------------------------------------
  ! Get grid variables
  !-----------------------------------------------------------------------------
  call check( NF90_INQ_VARID(ncid,trim(grid_lon_name),varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon) )
  WRITE(6,*) "lon(1) = ", lon(1)
  WRITE(6,*) "lon(nlon) = ", lon(nlon)
  call check( NF90_INQ_VARID(ncid,trim(grid_lat_name),varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat) )
  WRITE(6,*) "lat(1) = ", lat(1)
  WRITE(6,*) "lat(nlat) = ", lat(nlat)
  call check( NF90_INQ_VARID(ncid,trim(grid_lev_name),varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !-----------------------------------------------------------------------------
  ! lon2d, lat2, dx, and dy
  !-----------------------------------------------------------------------------
!!STEVE:MOM6: open new gridfile (ocean_hgrid.nc, gridfile3)
  INQUIRE(FILE=trim(gridfile3),EXIST=ex)
  if (.not. ex) then
    WRITE(6,*) "The file does not exist: ", trim(gridfile3)
    WRITE(6,*) "Exiting common_mom6.f90..."
    STOP (2)
  endif
  WRITE(6,'(A)') '  >> accessing file: ', trim(gridfile3)
  call check( NF90_OPEN(trim(gridfile3),NF90_NOWRITE,ncid3) )

  !STEVE: temporary until I figure out how to use ocean_hgrid 'supergrid' appropriately:
  ALLOCATE(work2d(nlon2d,nlat2d))

  WRITE(6,*) "Reading x coordinate (lon)..."
  call check( NF90_INQ_VARID(ncid3,trim(grid_lon2d_name),varid) )   ! Longitude for T-cell
! call check( NF90_GET_VAR(ncid3,varid,lon2d) )
  call check( NF90_GET_VAR(ncid3,varid,work2d) )    !STEVE: (ISSUE) ask GFDL for clarification
  !lon2d = work2d(1:nlon2d:2,1:nlat2d:2)             !STEVE: (ISSUE) ask GFDL for clarification
  lon2d = work2d(2:nlon2d:2,2:nlat2d:2)  !CDA: geolon (lon of T-Cell) calculated from supergrid x

  WRITE(6,*) "lon2d(1,1) = ", lon2d(1,1)
  WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)

  WRITE(6,*) "Reading y coordinate (lat)..."
  call check( NF90_INQ_VARID(ncid3,trim(grid_lat2d_name),varid) )   ! Latitude for T-cell
! call check( NF90_GET_VAR(ncid3,varid,lat2d) )
  call check( NF90_GET_VAR(ncid3,varid,work2d) )    !STEVE: (ISSUE) ask GFDL for clarification
  !lat2d = work2d(1:nlon2d:2,1:nlat2d:2)             !STEVE: (ISSUE) ask GFDL for clarification
  lat2d = work2d(2:nlon2d:2,2:nlat2d:2)    !CDA: geolat (lat of T-Cell) calculated from supergrid y

  WRITE(6,*) "lat2d(1,1) = ", lat2d(1,1)
  WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)

  !STEVE: I don't really need this anymore with the kdtree-based search
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
        d2 = cos(lat(j-1)) * cos(lat(j+1)) * (sin(dlon/2.0))**2
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

  !-----------------------------------------------------------------------------
  ! kmt data
  !-----------------------------------------------------------------------------
  !STEVE:MOM6: sum(h>0) (from MOM.res.nc) to find ocean layer <reference> depth, where depth > 0 (from ocean_togog.nc)
  call check( NF90_INQ_VARID(ncid,trim(grid_h_name),varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,height) )
  WRITE(6,*) "h(1,1,1) = ", height(1,1,1)
  WRITE(6,*) "h(nlon,nlat,nlev) = ", height(nlon,nlat,nlev)

  INQUIRE(FILE=trim(gridfile2),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", trim(gridfile2)
    WRITE(6,*) "Exiting common_mom6.f90..."
    STOP (2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', trim(gridfile2)
  call check( NF90_OPEN(trim(gridfile2),NF90_NOWRITE,ncid2) )

  call check( NF90_INQ_VARID(ncid2,trim(grid_depth_name),varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid2,varid,depth) )
  WRITE(6,*) "depth(1,1) = ", depth(1,1)
  WRITE(6,*) "depth(nlon,nlat) = ", depth(nlon,nlat)

  !-----------------------------------------------------------------------------
! !STEVE:MOM6: from ocean_topog.nc:
  !-----------------------------------------------------------------------------
  call check( NF90_INQ_VARID(ncid2,trim(grid_wet_name),varid) )        ! land/sea flag (0=land) for T-cell
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

  !-----------------------------------------------------------------------------
  ! Close the grid files:
  !-----------------------------------------------------------------------------
  call check( NF90_CLOSE(ncid) )
!!call check( NF90_CLOSE(ncid1) )
  call check( NF90_CLOSE(ncid2) )
! call check( NF90_CLOSE(ncid3) )

  ! Set model variables that depend on initialization and further processing.
  ! (e.g. lon0, lat0, lonf, latf, wrapgap, ...)
  CALL set_vars_model

  !-----------------------------------------------------------------------------
  ! Read the SSH model climatology if assimilating SLA's
  !-----------------------------------------------------------------------------
  if (DO_ALTIMETRY .and. DO_SLA) then
    INQUIRE(FILE=trim(SSHclm_file),EXIST=ex)
    if (ex) then
      ! Read in the model climatology
      CALL read_etaclm
    else
      WRITE(6,*) "The file does not exist: ", SSHclm_file
      WRITE(6,*) "Exiting common_mom6.f90..."
      STOP (1)
    endif
  endif

END SUBROUTINE set_common_oceanmodel


!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
SUBROUTINE read_etaclm

  USE netcdf
  USE params_model, ONLY: SSHclm_file, nlon, nlat
  USE vars_model,   ONLY: SSHclm_m
  
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
  do j=1,nlat
    do i=1,nlon
      !STEVE: Hopefully reading in meters here... (data might be in cm)
      SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
    enddo
  enddo

  call check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm


!STEVE: add this:
!-- Read a diagnostic file in mom6 netcdf format ---------------------------------------------------
SUBROUTINE read_diag(infile,v3d,v2d,prec_in)
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by MOM6
  !       The t,s,u,v,(ssh) fields are read in from the files so that the o-f data can be computed
  !       by the obsop_xxx.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  !STEVE: THUS, only the fields required for computing the OMFs need to be read in here.
  !       The full 3d model prognostic variables can be read in via read_restart
  USE netcdf
  USE vars_model,   ONLY: SSHclm_m
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_ssh
  USE params_model, ONLY: diag_lon_name, diag_lat_name, diag_lev_name
  USE params_model, ONLY: diag_temp_name, diag_salt_name
  USE params_model, ONLY: diag_u_name, diag_v_name, diag_h_name
  USE params_model, ONLY: diag_ssh_name, diag_height_name
  USE params_model, ONLY: diag_DO_temp, diag_DO_salt, diag_DO_u, diag_DO_v, diag_DO_ssh
  
  CHARACTER(*), INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  INTEGER :: prec ! precision, 1=single, 2=double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4 !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8 !(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,varid
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .false.
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  ! If prec is a provided argument, use indicated precision,
  if(PRESENT(prec_in))then
    prec = prec_in
  else
    ! otherwise default to single precision.
    prec = 1
  endif

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
  if (diag_DO_temp .or. diag_DO_salt) then
    WRITE(6,*) "read_diag:: opening file: ",infile
    call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
    WRITE(6,*) "read_diag :: just opened file ", infile
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!! t
  if (diag_DO_temp) then
    varname=diag_temp_name
    ivid=iv3d_t

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )

    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
    select case (prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4) )
        v3d(:,:,:,ivid) = buf4 !REAL(buf4,r_size)
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
  endif
 
  !!! s
  if (diag_DO_salt) then 
    varname=diag_salt_name
    ivid=iv3d_s

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
    select case (prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4) )
        v3d(:,:,:,ivid) = buf4 !REAL(buf4,r_size)
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
  endif

  !!! u
  if (diag_DO_u) then
    varname=diag_u_name
    ivid=iv3d_u

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
    select case (prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4) )
        v3d(:,:,:,ivid) = buf4 !REAL(buf4,r_size)
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
  endif

  !!! v
  if (diag_DO_v) then
    varname=diag_v_name
    ivid=iv3d_v

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
      select case (prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4) )
        v3d(:,:,:,ivid) = buf4 !REAL(buf4,r_size)
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
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! ssh
  if (diag_DO_ssh) then
    varname=diag_ssh_name
    ivid=iv2d_ssh

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
    select case (prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
        if (DO_SLA) then
!         v2d(:,:,ivid) = REAL(buf4(:,:,1),r_size) - SSHclm_m(:,:)
          v2d(:,:,ivid) = buf4(:,:,1) - SSHclm_m(:,:)
        else
          v2d(:,:,ivid) = buf4(:,:,1) !REAL(buf4(:,:,1),r_size)
        endif
      case(2)
        buf8=0.0d0
        call check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
        if (DO_SLA) then
          v2d(:,:,ivid) = buf8(:,:,1) - SSHclm_m(:,:)
        else
          v2d(:,:,ivid) = buf8(:,:,1)
        endif 
    end select
    if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)

  endif

  !STEVE: assuming all variables were in a single diagnostic file
  call check( NF90_CLOSE(ncid) )
  
END SUBROUTINE read_diag


!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the MOM.res.nc, MOM.res_1.nc, etc. restart files produced by MOM6
  !       The t,s,u,v,(ssh) fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsop_xxx.f90 program prior
  !       to running letkf.f90, so the diagnostic files should not have to be touched during letkf runtime.
  USE netcdf
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_ssh
  USE params_model, ONLY: iv3d_h, iv2d_eta
  USE params_letkf, ONLY: DO_UPDATE_H
  USE vars_model,   ONLY: SSHclm_m
  USE params_model, ONLY: rsrt_tsbase, rsrt_uvbase, rsrt_hbase
  USE params_model, ONLY: rsrt_lon_name, rsrt_lat_name, rsrt_lev_name
  USE params_model, ONLY: rsrt_temp_name, rsrt_salt_name
  USE params_model, ONLY: rsrt_u_name, rsrt_v_name
  USE params_model, ONLY: rsrt_h_name, rsrt_ssh_name
  
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1 = single, 2 = double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  CHARACTER(slen) :: tsfile,uvfile, sffile, bfile, drfile ! (TS) (UV) (SFC) (barotropic - eta) (DRIFTERS)
  INTEGER :: i,j,k
  INTEGER :: ncid,varid
  REAL(r_size), PARAMETER :: vmax = 1.0e18
  !STEVE:
  REAL(r_size) :: meanSSH !STEVE: for temporary SSH estimate based on heat content
  REAL(r_size) :: videpth !STEVE: depth of vertically integrated heat content
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .false.
  CHARACTER(3) :: MEM3
  CHARACTER(32) :: sfc_infile
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  tsfile = trim(infile)//'.'//trim(rsrt_tsbase)
  uvfile = trim(infile)//'.'//trim(rsrt_uvbase)
  bfile  = trim(infile)//'.'//trim(rsrt_hbase)

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

  !-----------------------------------------------------------------------------
  !!! t
  !-----------------------------------------------------------------------------
  varname=rsrt_temp_name
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

  !-----------------------------------------------------------------------------
  !!! s
  !-----------------------------------------------------------------------------
  varname=rsrt_salt_name
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

  if (DO_UPDATE_H) then
    !---------------------------------------------------------------------------
    !!! h
    !---------------------------------------------------------------------------
    varname=rsrt_h_name
    ivid=iv3d_h

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

    if (dodebug) WRITE(6,*) "read_restart :: just got data for variable thickness"
    if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable thickness"

   !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-H"
      WRITE(6,*) "read_restart :: tsfile = ", tsfile
      do k=1,nlev
        WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_h) = ", MAXVAL(v3d(:,:,k,iv3d_h))
      enddo
    endif
!   !STEVE: end
  endif

  !-----------------------------------------------------------------------------
  !!! u
  !-----------------------------------------------------------------------------
  varname=rsrt_u_name
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

  !-----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the U/V netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------
  call check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart :: just opened file ", uvfile

  !-----------------------------------------------------------------------------
  !!! v
  !-----------------------------------------------------------------------------
  varname=rsrt_v_name
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
    varname=rsrt_ssh_name
    ivid=iv2d_ssh

    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    select case(prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
        v2d(:,:,ivid) = buf4(:,:,1)
      case(2)
        buf8=0.0d0
        call check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
        v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl)
    end select

    if (dodebug) WRITE(6,*) "read_restart :: just got data for variable ssh"
    if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable SSH"

    ! Convert SSH eta stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    if (DO_SLA) v2d(:,:,ivid) = v2d(:,:,iv2d_eta) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_restart :: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
      WRITE(6,*) "min val for level v2d(:,:,iv2d_eta) = ", MINVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

  else
    WRITE(6,*) "read_restart :: Skipping reading eta. Enabled only when using DO_ALTIMETRY=.true."
  endif altimetry

  call check( NF90_CLOSE(ncid) )

  ! For additional variables:
  ! E.g. surface fluxes, drifters

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

END SUBROUTINE read_restart


!-----------------------------------------------------------------------
! Write a set of MOM6 restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d_in,v2d_in,prec)
  USE netcdf
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_ssh
  USE params_model, ONLY: iv3d_h, iv2d_eta
  USE params_letkf, ONLY: DO_UPDATE_H, DO_SLA
  USE vars_model,   ONLY: SSHclm_m
  USE params_model, ONLY: rsrt_tsbase, rsrt_uvbase, rsrt_hbase
  USE params_model, ONLY: rsrt_lon_name, rsrt_lat_name, rsrt_lev_name
  USE params_model, ONLY: rsrt_temp_name, rsrt_salt_name
  USE params_model, ONLY: rsrt_u_name, rsrt_v_name
  USE params_model, ONLY: rsrt_h_name, rsrt_ssh_name
  USE params_model, ONLY: do_physlimit
  USE params_model, ONLY: min_t, max_t, min_s, max_s
  USE params_model, ONLY: min_u, max_u, min_v, max_v
  USE params_model, ONLY: min_eta, max_eta
  USE params_model, ONLY: min_h, max_h
  
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1 = single, 2 = double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY)
  INTEGER :: ncid,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: mlon, mlat
  REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
  CHARACTER(slen) :: sfc_outfile, sfc_infile
  CHARACTER(3) :: MEM3
  LOGICAL :: dodebug = .true.

  !STEVE: this is the routine that writes out the individual analysis files for
  !       each esnsemble member in netcdf format.

  tsfile = trim(outfile)//'.'//trim(rsrt_tsbase)
  uvfile = trim(outfile)//'.'//trim(rsrt_uvbase)
  bfile  = trim(outfile)//'.'//trim(rsrt_hbase)

  if (prec == 1) then
    WRITE(6,*) "common_mom6.f90::write_restart:: input argument prec=1 (single precision) is not yet supported. Update write_restart subroutine. EXITING..."
    STOP (24)
  endif

  select case(prec)
    case(1)
      if (dodebug) WRITE(6,*) "write_restart::ALLOCATE buf4..."
      ALLOCATE(buf4(nlon,nlat,nlev))
      if (dodebug) WRITE(6,*) "Done ALLOCATE buf4."
    case(2)
      if (dodebug) WRITE(6,*) "write_restart::ALLOCATE buf8..."
      ALLOCATE(buf8(nlon,nlat,nlev))
      if (dodebug) WRITE(6,*) "Done ALLOCATE buf8."
  end select

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !STEVE: open temp/salt file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (dodebug) WRITE(6,*) "Opening file for temperature and salinity..."
  call check( NF90_OPEN(tsfile,NF90_WRITE,ncid) )
  if (dodebug) WRITE(6,*) "Opened file for temperature and salinity."

  if (dodebug) WRITE(6,*) "Get id for temperature..."
  call check( NF90_INQ_VARID(ncid,rsrt_temp_name,varid) )

  buf8 = REAL(v3d_in(:,:,:,iv3d_t),r_size)

  ! Apply physical limits to the output analysis:
  if (do_physlimit) then
    WHERE (buf8 < min_t) buf8 = min_t
    WHERE (buf8 > max_t) buf8 = max_t
  endif

  if (dodebug) WRITE(6,*) "Write data for temperature..."
  call check( NF90_PUT_VAR(ncid,varid,buf8) )

  !-----------------------------------------------------------------------------
  !!! s
  !-----------------------------------------------------------------------------
  if (dodebug) WRITE(6,*) "Get id for salinity..."
  call check( NF90_INQ_VARID(ncid,rsrt_salt_name,varid) )

  buf8 = REAL(v3d_in(:,:,:,iv3d_s),r_size)

  ! Apply physical limits to the output analysis:
  if (do_physlimit) then
    WHERE (buf8 < min_s) buf8 = min_s
    WHERE (buf8 > max_s) buf8 = max_s
  endif

  if (dodebug) WRITE(6,*) "Write data for salinity..."
  call check( NF90_PUT_VAR(ncid,varid,buf8) )

  !-----------------------------------------------------------------------------
  !!! h (MOM6)
  !-----------------------------------------------------------------------------
  if (DO_UPDATE_H) then
    if (dodebug) WRITE(6,*) "Get id for layer thickness (h)..."
    call check( NF90_INQ_VARID(ncid,rsrt_h_name,varid) )
    if (dodebug) WRITE(6,*) "Write data for layer thickness (h)..."

    buf8 = REAL(v3d_in(:,:,:,iv3d_h),r_size)

    ! Apply physical limits to the output analysis:
    if (do_physlimit) then
      WHERE (buf8 < min_h) buf8 = min_h
      WHERE (buf8 > max_h) buf8 = max_h
    endif

    call check( NF90_PUT_VAR(ncid,varid,buf8) )

  endif

  !-----------------------------------------------------------------------------
  !!! u
  !-----------------------------------------------------------------------------
  if (dodebug) WRITE(6,*) "Get id for u component of velocity..."
  call check( NF90_INQ_VARID(ncid,rsrt_u_name,varid) )

  buf8 = REAL(v3d_in(:,:,:,iv3d_u),r_size)

  ! Apply physical limits to the output analysis:
  if (do_physlimit) then
    WHERE (buf8 < min_u) buf8 = min_u
    WHERE (buf8 > max_u) buf8 = max_u
  endif

  if (dodebug) WRITE(6,*) "Write data for u component of velocity..."
  call check( NF90_PUT_VAR(ncid,varid,buf8) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! uv file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (dodebug) WRITE(6,*) "Opening file for velocity and surface height..."
  call check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )

  !-----------------------------------------------------------------------------
  !!! v
  !-----------------------------------------------------------------------------
  if (dodebug) WRITE(6,*) "Get id for v component of velocity..."
  call check( NF90_INQ_VARID(ncid,rsrt_v_name,varid) )

  buf8 = REAL(v3d_in(:,:,:,iv3d_v),r_size)

  ! Apply physical limits to the output analysis:
  if (do_physlimit) then
    WHERE (buf8 < min_v) buf8 = min_v
    WHERE (buf8 > max_v) buf8 = max_v
  endif

  if (dodebug) WRITE(6,*) "Write data for v component of velocity..."
  call check( NF90_PUT_VAR(ncid,varid,buf8) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the updated eta_t to analysis restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (DO_ALTIMETRY) then
! if (dodebug) WRITE(6,*) "Opening file for surface height..."
!   call check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    if (dodebug) WRITE(6,*) "Get id for ssh..."
    call check( NF90_INQ_VARID(ncid,rsrt_ssh_name,varid) )

    buf8(:,:,1) = REAL(v2d_in(:,:,iv2d_eta),r_size)

    ! Apply physical limits to the output analysis:
    if (do_physlimit) then
      WHERE (buf8(:,:,1) < min_eta) buf8(:,:,1) = min_eta
      WHERE (buf8(:,:,1) > max_eta) buf8(:,:,1) = max_eta
    endif

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    if (DO_SLA) then
      buf8(:,:,1) = buf8(:,:,1) + SSHclm_m(:,:)
    endif

    if (dodebug) WRITE(6,*) "Write data for ssh..."
    call check( NF90_PUT_VAR(ncid,varid,buf8(:,:,1)) )
  endif

  call check( NF90_CLOSE(ncid) )

END SUBROUTINE write_restart


!-----------------------------------------------------------------------
! Read a grid file with single precision
!-----------------------------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: nij0
  
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
  do n=1,nv3d
    do k=1,nlev
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    enddo
  enddo

  do n=1,nv2d
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  enddo

  CLOSE(iunit)

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

END SUBROUTINE read_grd


!-----------------------------------------------------------------------
! Write a grid file in single precision
!-----------------------------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: nij0

  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL, PARAMETER :: dodebug=.false.

  if (dodebug) print *, "write_grd4:: open filename = ",filename
  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_grd4:: nij0,iolength = ", nij0,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  do n=1,nv3d
    do k=1,nlev
      if (dodebug) print *, "write_grd4:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    enddo
  enddo

  do n=1,nv2d
    if (dodebug) print *, "write_grd4:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  enddo

  CLOSE(iunit)

END SUBROUTINE write_grd


!-----------------------------------------------------------------------
! Monitor  (ISSUE) Depricated. Subroutine slated for removal.
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
 
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: nv3d, nv2d
  USE params_model, ONLY: element

  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  do k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    do n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    enddo
  enddo

  do n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  enddo

END SUBROUTINE monit_grd


!STEVE: put this somewhere more general (non model-specific) (ISSUE)
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  
  USE params_model, ONLY: nlev
  USE params_model, ONLY: nv3d, nv2d

  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

! ALLOCATE(v3dm(nij,nlev,nv3d),v2dm(nij,nv2d))

  do n=1,nv3d
    do k=1,nlev
      do i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        do m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        enddo
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      enddo
    enddo
  enddo

  do n=1,nv2d
    do i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      do m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      enddo
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    enddo
  enddo

END SUBROUTINE ensmean_grd


!-----------------------------------------------------------------------
subroutine check(status)
  USE netcdf
  
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

END MODULE common_oceanmodel
