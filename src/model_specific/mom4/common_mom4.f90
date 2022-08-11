MODULE common_oceanmodel
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

  IMPLICIT NONE

  PUBLIC

  REAL(r_size), DIMENSION(:), ALLOCATABLE :: xc,yc
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: sfc_data
  INTEGER, DIMENSION(:), ALLOCATABLE :: xi,yi

  !STEVE: for debugging
  LOGICAL, PRIVATE :: doverbose = .true.

CONTAINS


SUBROUTINE set_common_oceanmodel
!===============================================================================
! Initialize the module
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER :: ncid,varid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex
  LOGICAL :: dodebug = .false.
  !STEVE: for verifying against input netcdf file

  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'
  CALL initialize_params_model ! (checks to make sure it is initialized)
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  if (DO_SLA) then
    INQUIRE(FILE=trim(SSHclm_file),EXIST=ex)
    if (ex) then
      ! Read in the model climatology
      CALL read_etaclm()
    ELSE
      WRITE(6,*) "The model SSH climatology file does not exist: ", SSHclm_file
      WRITE(6,*) "Exiting set_common_mom4.f90..."
      STOP (1)
    ENDIF
  endif

  !
  ! Lon, Lat, f, orography
  !
!STEVE: this part adapted from ROMS, update from MOM4 netcdf files:
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  if (.not. ex) then
    WRITE(6,*) "The file does not exist: ", gridfile 
    WRITE(6,*) "Exiting set_common_mom4.f90..."
    STOP (2)
  ENDIF

  if (doverbose) then
    WRITE(6,'(A)') '  >> accessing file: ', gridfile
  endif

!(UPDATE) with lon2d/lat2d search for observation (LON2D/LAT2D)
  !-----------------------------------------------------------------------------
  ! Get 1d and 2d longitude fields:
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )
  CALL check( NF90_INQ_VARID(ncid,'grid_x_T',varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lon) )

  if (doverbose) then
    WRITE(6,*) "lon(1) = ", lon(1)
    WRITE(6,*) "lon(nlon) = ", lon(nlon)
  endif

  CALL check( NF90_INQ_VARID(ncid,'x_T',varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lon2d) )

  if (doverbose) then
    WRITE(6,*) "lon2d(1,1)       = ", lon2d(1,1)
    WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)
  endif

  !-----------------------------------------------------------------------------
  ! Get 1d and 2d latitude fields:
  !-----------------------------------------------------------------------------
  CALL check( NF90_INQ_VARID(ncid,'grid_y_T',varid) )   ! Latitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lat) )

  if (doverbose) then
    WRITE(6,*) "lat(1) = ", lat(1)
    WRITE(6,*) "lat(nlat) = ", lat(nlat)
  endif

  CALL check( NF90_INQ_VARID(ncid,'y_T',varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lat2d) )

  if (doverbose) then
    WRITE(6,*) "lat2d(1,1)       = ", lat2d(1,1)
    WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)
  endif

  !-----------------------------------------------------------------------------
  ! Get 1d and 2d depths of cell:
  !-----------------------------------------------------------------------------
  CALL check( NF90_INQ_VARID(ncid,'zt',varid) )      ! depth of T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lev) )

  if (doverbose) then
    WRITE(6,*) "lev(1) = ", lev(1)
    WRITE(6,*) "lev(nlev) = ", lev(nlev)
  endif

  CALL check( NF90_INQ_VARID(ncid,'depth_t',varid) )      ! depth of T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lev2d) )

  if (doverbose) then
    WRITE(6,*) "lev2d(1,1)       = ", lev2d(1,1)
    WRITE(6,*) "lev2d(nlon,nlat) = ", lev2d(nlon,nlat)
  endif

  !-----------------------------------------------------------------------------
  ! dx and dy
  !-----------------------------------------------------------------------------
  CALL check( NF90_INQ_VARID(ncid,'ds_01_21_T',varid) )    ! width of T_cell (meters)
  CALL check( NF90_GET_VAR(ncid,varid,dx) ) 
  CALL check( NF90_INQ_VARID(ncid,'ds_10_12_T',varid) )    ! height of T_cell (meters)
  CALL check( NF90_GET_VAR(ncid,varid,dy) ) 

  if (doverbose) then
    WRITE(6,*) "set_common_oceanmodel:: grid_spec.nc MIN(dx) = ", MINVAL(dx)
    WRITE(6,*) "set_common_oceanmodel:: grid_spec.nc MAX(dx) = ", MAXVAL(dx)
    WRITE(6,*) "set_common_oceanmodel:: grid_spec.nc MIN(dy) = ", MINVAL(dy)
    WRITE(6,*) "set_common_oceanmodel:: grid_spec.nc MAX(dy) = ", MAXVAL(dy)
  endif

  !
  ! kmt data
  !
  CALL check( NF90_INQ_VARID(ncid,'num_levels',varid) ) ! number of vertical T-cells
  CALL check( NF90_GET_VAR(ncid,varid,kmt0) )
  if (doverbose) then
    WRITE(6,*) "kmt0(1,1) = ", kmt0(1,1)
    WRITE(6,*) "kmt0(nlon,nlat) = ", kmt0(nlon,nlat)
  endif
  kmt = NINT(kmt0)

  if (doverbose) then
    WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile
    WRITE(6,*) "dx(1,1) = ", dx(1,1)
    WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
    WRITE(6,*) "dy(1,1) = ", dy(1,1)
    WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)
  endif

  !STEVE: needed for computing the AMOC based on the streamfunction calculation:
  CALL check( NF90_INQ_VARID(ncid,'zb',varid) )      ! depth of T-cell
  CALL check( NF90_GET_VAR(ncid,varid,zb) )

  if (doverbose) then
    WRITE(6,*) "zb(1) = ", zb(1)
    WRITE(6,*) "zb(nlev) = ", zb(nlev)
  endif

  ! Compute dz:
  dz(1) = zb(1)
  do k=2,nlev
    dz(k) = zb(k)-zb(k-1)
  enddo

  ! Close the grid_spec.nc file:
  CALL check( NF90_CLOSE(ncid) )

  ! Set model variables that depend on initialization and further processing.
  ! (e.g. lon0, lat0, lonf, latf, wrapgap, ...)
  CALL set_vars_model 

END SUBROUTINE set_common_oceanmodel


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
    stop "NetCDF read error! EXITING..."
  end if
END SUBROUTINE check


SUBROUTINE read_etaclm()
!===============================================================================
! Read in the model mean climatology (e.g. 1991-1999), the real climatology 
! assumed to be subtracted already from the observation data.
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
  IMPLICIT NONE

  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: i,j
  INTEGER :: ncid, varid

  ! read the model SSH climatology netcdf file
  ! read into: SSHclm_m
  CALL check( NF90_OPEN(SSHclm_file,NF90_NOWRITE,ncid) )

  if (doverbose) then
    WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file
  endif

  buf4=0.0
  CALL check( NF90_INQ_VARID(ncid,'ssh',varid) )
  CALL check( NF90_GET_VAR(ncid,varid,buf4) )
  do j=1,nlat
    do i=1,nlon
      !STEVE: Hopefully reading in meters here... (data might be in cm)
      SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
    enddo
  enddo

  CALL check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm


SUBROUTINE read_history(infile,v3d,v2d)
!===============================================================================
! Read in a set of netcdf-format mom4 history files
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(INOUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(INOUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: ncid,varid
  INTEGER :: i,j,k
  LOGICAL :: dodebug = .false.

  ! For now, just read the mixed layer depth if it is an option, otherwise return
  if (DO_MLD) then
    WRITE(6,*) "Reading mixed layer depth from history file..."

    !---------------------------------------------------------------------------
    ! Open the T/S netcdf restart file
    !---------------------------------------------------------------------------
    CALL check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
 
    if (doverbose) then
      WRITE(6,*) "read_history:: just opened file ", infile
    endif  

    !---------------------------------------------------------------------------
    ! mixed layer depth (mld)
    !---------------------------------------------------------------------------
    buf4=0.0
    CALL check( NF90_INQ_VARID(ncid,'mld',varid) )
    CALL check( NF90_GET_VAR(ncid,varid,buf4) )
    if (dodebug) WRITE(6,*) "read_history:: just got data for variable mld"
    v2d(:,:,iv2d_mld) = buf4
    if (dodebug) WRITE(6,*) "read_history:: finished processing data for variable: mld"
  
    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-MLD:"
      WRITE(6,*) "read_history:: infile = ", infile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_mld) = ",MAXVAL(v2d(:,:,iv2d_mld))
    endif
    ! !STEVE: end

  endif 

END SUBROUTINE read_history


SUBROUTINE read_diag(infile,v3d,v2d,prec_in)
!===============================================================================
! Read in a set of netcdf-format mom4 diagnostic files (for now, we're using the restarts)
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  INTEGER :: prec ! precision, 1=single, 2=double
  REAL(r_sngl), ALLOCATABLE :: v3d0(:,:,:,:)
  REAL(r_sngl), ALLOCATABLE :: v2d0(:,:,:)

  ! If prec is a provided argument, use indicated precision,
  if(PRESENT(prec_in))then
    prec = prec_in
  else
    ! otherwise default to single precision.
    prec = 1
  endif

  ALLOCATE(v3d0(nlon,nlat,nlev,nv3d), v2d0(nlon,nlat,nv2d))
  CALL read_restart(infile,v3d0,v2d0,prec)
  v3d = REAL(v3d0,r_size)
  v2d = REAL(v2d0,r_size)

END SUBROUTINE read_diag


!SUBROUTINE read_grd4(infile,v3d,v2d)
!!===============================================================================
!! Read subroutine for backwards compatibility
!!===============================================================================
!  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: infile
!  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
!
!  CALL read_restart(infile,v3d,v2d,1)
!
!END SUBROUTINE read_grd4


SUBROUTINE read_restart(infile,v3d,v2d,prec_in)
!===============================================================================
! Read in a set of netcdf-format mom4 restart files in single precision
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, OPTIONAL, INTENT(IN) :: prec_in ! precision, 1 = single, 2 = double

  INTEGER :: prec
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile,hsfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY: ocean_barotropic.res.nc, contains eta) 
                                                                ! (history, contains mld)
  INTEGER :: ncid,varid,ivid
  CHARACTER(slen) :: varname
  INTEGER :: i,j,k
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
  LOGICAL :: dodebug = .false.
  REAL(r_size), PARAMETER :: vmax = 1.0e18

  
  !-----------------------------------------------------------------------------
  ! Set default precision:
  !-----------------------------------------------------------------------------
  if (.not. present(prec_in)) then  
    prec = 1 
  else
    prec = prec_in
  endif
 
  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

  tsfile = trim(infile)//'.'//trim(ts_basefile) !ocean_temp_salt.res.nc'
  uvfile = trim(infile)//'.'//trim(uv_basefile) !ocean_velocity.res.nc'
  sffile = trim(infile)//'.'//trim(sf_basefile) !ocean_sbc.res.nc'
  bfile  = trim(infile)//'.'//trim(sh_basefile) !ocean_barotropic.res.nc'
  hsfile  = trim(infile)//'.'//trim(hs_basefile) !ocean_TS.nc'

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  CALL check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )

  !-----------------------------------------------------------------------------
  ! t
  !-----------------------------------------------------------------------------
  varname='temp'
  ivid=iv3d_t 

  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      CALL check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      CALL check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-TEMP"
    WRITE(6,*) "read_restart:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t) = ", MAXVAL(v3d(:,:,k,iv3d_t))
    enddo 
  endif
! !STEVE: end

  !-----------------------------------------------------------------------------
  ! s
  !-----------------------------------------------------------------------------
  varname='salt'
  ivid=iv3d_s

  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      CALL check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      CALL check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_restart:: tsfile = ", tsfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_s) = ", MAXVAL(v3d(:,:,k,iv3d_s))
    enddo 
  endif
! !STEVE: end

  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Open file with UV data:
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )

  !-----------------------------------------------------------------------------
  ! u
  !-----------------------------------------------------------------------------
  varname='u'
  ivid=iv3d_u

  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      CALL check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      CALL check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-U"
    WRITE(6,*) "read_restart:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_u) = ", MAXVAL(v3d(:,:,k,iv3d_u))
    enddo 
  endif
! !STEVE: end

  !-----------------------------------------------------------------------------
  ! v
  !-----------------------------------------------------------------------------
  varname='v'
  ivid=iv3d_v

  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      CALL check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      CALL check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-V"
    WRITE(6,*) "read_restart:: uvfile = ", uvfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_v) = ", MAXVAL(v3d(:,:,k,iv3d_v))
    enddo 
  endif
! !STEVE: end

  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Set the SST and SSS data for the SFC
  !-----------------------------------------------------------------------------
  if (.false.) then
    v2d(:,:,iv2d_sst) = v3d(:,:,ilev_sfc,iv3d_t)
    v2d(:,:,iv2d_sss) = v3d(:,:,ilev_sfc,iv3d_s)
  else
    CALL check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
    WRITE(6,*) "read_restart:: just opened file ", sffile

    !---------------------------------------------------------------------------
    ! SST
    !---------------------------------------------------------------------------
    varname='t_surf'
    ivid=iv2d_sst
    CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) ) 
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)
    
    select case(prec)
      case(1)
        buf4=0.0
        CALL check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) ) 
        WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = buf4(:,:,1) - t0c ! kelvin to dec C
      case(2)
        buf8=0.0d0
        CALL check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) ) 
        WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl) - t0c ! kelvin to deg C
    end select
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid


! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SST"
      WRITE(6,*) "read_restart:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sst) = ", MAXVAL(v2d(:,:,iv2d_sst))
    endif
! !STEVE: end

    !---------------------------------------------------------------------------
    ! SSS
    !---------------------------------------------------------------------------
    varname='s_surf'
    ivid=iv2d_sss
    CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)

    select case(prec)
      case(1)
        buf4=0.0
        CALL check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
        v2d(:,:,ivid) = buf4(:,:,1)
      case(2)
        buf8=0.0d0
        CALL check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
        v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl)
    end select
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid

! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SSS"
      WRITE(6,*) "read_restart:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sss) = ", MAXVAL(v2d(:,:,iv2d_sss))
    endif
! !STEVE: end

    CALL check( NF90_CLOSE(ncid) )
  endif
  
  !-----------------------------------------------------------------------------
  ! Open the ALTIMETRY netcdf diagnostic file (SSH)
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart:: just read file ", sffile

  !-----------------------------------------------------------------------------
  ! SSH
  !-----------------------------------------------------------------------------
  varname='sea_lev'
  ivid=iv2d_ssh
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)

  select case(prec)
    case(1)
      buf4=0.0
      CALL check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
      v2d(:,:,ivid) = buf4(:,:,1)
    case(2)
      buf8=0.0d0
      CALL check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
      v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl)
  end select
  if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SSH"
    WRITE(6,*) "read_restart:: sffile = ", sffile
    WRITE(6,*) "max val for level v3d(:,:,iv2d_ssh) = ", MAXVAL(v2d(:,:,ivid))
  endif
! !STEVE: end

  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Do altimetry
  !-----------------------------------------------------------------------------
  altimetry : if(DO_ALTIMETRY) then
    !STEVE: use the sea level perturbation from ocean_barotropic.res.nc
    CALL check( NF90_OPEN(bfile,NF90_NOWRITE,ncid) )
    if (doverbose) WRITE(6,*) "read_restart:: just opened file ", bfile

    varname='eta_t'
    ivid=iv2d_eta
    CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable eta_t"
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable SSH"

    !---------------------------------------------------------------------------
    ! SSH
    !---------------------------------------------------------------------------
    select case(prec)
      case(1)
        buf4=0.0
        CALL check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) )
        v2d(:,:,ivid) = buf4(:,:,1)
      case(2)
        buf8=0.0d0
        CALL check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) )
        v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl)
    end select

    ! Convert SSH eta stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    if (DO_SLA) v2d(:,:,ivid) = v2d(:,:,ivid) - SSHclm_m(:,:)

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_restart:: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

    CALL check( NF90_CLOSE(ncid) )
  endif altimetry

  mld : if (DO_MLD .and. .not. DO_MLD_MAXSPRD) then
    ! Add the ensemble of mixed-layer depths to the model state for analysis
    CALL read_history(hsfile,v3d,v2d)
  endif mld

  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

! !STEVE: debug test
  if (.false.) then
    testfile = "test_read4.grd"
!   CALL write_grd4(trim(testfile),v3d,v2d)

    iunit=55
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=testfile,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

    WRITE(6,*) "Writing to", testfile
    irec=1
    do n=1,nv3d
      do k=1,nlev
        WRITE(6,*) "n, k, irec = ", n, k, irec
        WRITE(6,*) "max v3d(n) = ", MAXVAL(v3d(:,:,k,n))
        WRITE(iunit,REC=irec) v3d(:,:,k,n)
        irec = irec + 1
      enddo
    enddo

    do n=1,nv2d
      WRITE(iunit,REC=irec) v2d(:,:,n)
      irec = irec + 1
    enddo
    CLOSE(iunit)

    WRITE(6,*) "Initially read from file: ", infile
    WRITE(6,*) "STOP 10"
    STOP (10)
  endif
! !STEVE: debug end

END SUBROUTINE read_restart


!SUBROUTINE write_grd4(outfile,v3d_in,v2d_in)
!!===============================================================================
!! Write subroutine for backwards compatibility
!!===============================================================================
!  USE params_model
!  USE vars_model
!  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: outfile
!  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
!
!  CALL write_restart(outfile,v3d_in,v2d_in) !,1)
!
!END SUBROUTINE write_grd4


SUBROUTINE write_restart(outfile,v3d_in,v2d_in) !,prec_in)
!===============================================================================
! Write out a set of netcdf-format restart files in single precision
!===============================================================================
  USE netcdf
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
! INTEGER, OPTIONAL, INTENT(IN) :: prec_in

! INTEGER :: prec
! REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
! REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  REAL(r_sngl), ALLOCATABLE :: v3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2d(:,:,:) !(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: t3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  CHARACTER(slen) :: tsfile,uvfile, sffile,drfile, bfile ! (TS) (UV) (SFC) (DRIFTERS) (ALTIMETRY)
  INTEGER :: ncid,varid,ivid
  CHARACTER(slen) :: varname
  INTEGER :: m,k,j,i !STEVE: for debugging

  ! Set default precision:
! if (.not. present(prec_in)) then
!   prec = 1
! else
!   prec = prec_in
! endif

  !STEVE: this is the routine that writes out the individual analysis files for
  !       each esnsemble member in netcdf format.

  tsfile = trim(outfile)//'.'//trim(ts_basefile) !ocean_temp_salt.res.nc'
  uvfile = trim(outfile)//'.'//trim(uv_basefile) !ocean_velocity.res.nc'
  sffile = trim(outfile)//'.'//trim(sf_basefile) !ocean_sbc.res.nc'
  bfile  = trim(outfile)//'.'//trim(sh_basefile) !ocean_barotropic.res.nc'

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

        if (v3d(i,j,k,iv3d_t) < min_t) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = min_t
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sst) < min_t) v2d(i,j,iv2d_sst) = min_t

        if (v3d(i,j,k,iv3d_t) > max_t) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = max_t
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sst) > max_t) v2d(i,j,iv2d_sst) = max_t

        if (v3d(i,j,k,iv3d_s) < min_s ) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = min_s
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sss) < min_s) v2d(i,j,iv2d_sss) = min_s

        if (v3d(i,j,k,iv3d_s) > max_s) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = max_s
        endif

        if (k .eq. 1 .and. v2d(i,j,iv2d_sss) > max_s) v2d(i,j,iv2d_sss) = max_s

      enddo
    enddo
  enddo
  endif
  
  !-----------------------------------------------------------------------------
  ! Open temp/salt file
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(tsfile,NF90_WRITE,ncid) )

  !!! t
  varname = 'temp'
  ivid = iv3d_t
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  !!! s
  varname = 'salt'
  ivid = iv3d_s
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Open uv file
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )

  !!! u
  varname = 'u'
  ivid = iv3d_u
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  !!! v
  varname = 'v'
  ivid = iv3d_v
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Open sfc file
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(sffile,NF90_WRITE,ncid) )
  
  !!! SST
  varname = 't_surf'
  ivid = iv2d_sst
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = v2d(:,:,ivid) + t0c ! kelvin
  CALL check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )

  !!! SSS
  varname = 's_surf'
  ivid = iv2d_sss
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )

  !!! SSH
  varname = 'sea_lev'
  ivid = iv2d_ssh
  CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  CALL check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )
  CALL check( NF90_CLOSE(ncid) )

  !-----------------------------------------------------------------------------
  ! Write the updated eta_t to analysis restart file
  !-----------------------------------------------------------------------------
  if (DO_ALTIMETRY) then
    CALL check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    varname = 'eta_t'
    ivid = iv2d_eta
    CALL check( NF90_INQ_VARID(ncid,trim(varname),varid) )

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,ivid) = v2d(:,:,ivid) + SSHclm_m(:,:)

    CALL check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )
    CALL check( NF90_CLOSE(ncid) )
  endif

END SUBROUTINE write_restart


!SUBROUTINE read_grd(filename,v3d,v2d)
!!===============================================================================
!! Read in an letkf grd-format binary file
!!===============================================================================
!  USE params_model
!  USE vars_model
!  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
!  REAL(r_sngl), ALLOCATABLE :: buf4(:,:) !(nlon,nlat)
!  INTEGER :: iunit,iolen
!  INTEGER :: k,n,irec
!
!  ALLOCATE(buf4(nlon,nlat))
!
!  iunit=11
!  INQUIRE(IOLENGTH=iolen) iolen
!  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
!
!  irec=1
!  do n=1,nv3d
!    do k=1,nlev
!      READ(iunit,REC=irec) buf4
!      irec = irec + 1
!      v3d(:,:,k,n) = REAL(buf4,r_size)
!    enddo
!  enddo
!
!  do n=1,nv2d
!    READ(iunit,REC=irec) buf4
!    irec = irec + 1
!    v2d(:,:,n) = REAL(buf4,r_size)
!  enddo
!
!  CLOSE(iunit)
!
!END SUBROUTINE read_grd


SUBROUTINE read_grd(filename,v3d,v2d)
!===============================================================================
! Read in an letkf grd-format binary file in single precision
!===============================================================================
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
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

END SUBROUTINE read_grd


SUBROUTINE write_grd(filename,v3d,v2d)
!===============================================================================
! Write out an letkf grd-format binary file in single precision
!===============================================================================
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL :: dodebug=.false.

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

!STEVE: move this somewhere more general (ISSUE)
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
!===============================================================================
! Compute the ensemble mean
!===============================================================================
  USE params_model
  USE vars_model
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA, DO_MLD_MAXSPRD  
  IMPLICIT NONE
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


END MODULE common_oceanmodel
