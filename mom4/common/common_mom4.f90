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
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------------
  REAL(r_size),SAVE :: dy2(nlat)
  REAL(r_size),SAVE :: fcori(nlat)
  REAL(r_size),SAVE :: wet(nlon,nlat)                !(OCEAN)
  REAL(r_size),SAVE :: area_t(nlon,nlat)             !(OCEAN)
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

  INTEGER :: islot
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: xc,yc
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: sfc_data
  INTEGER, DIMENSION(:), ALLOCATABLE :: xi,yi

  !STEVE: for debugging
  LOGICAL :: dodebug = .false.
  LOGICAL :: doverbose = .true.

CONTAINS


SUBROUTINE set_common_mom4
!===============================================================================
! Initialize the module
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER :: ncid,varid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex

  WRITE(6,'(A)') 'Hello from set_common_mom4'

  if (DO_ALTIMETRY) then
    INQUIRE(FILE=trim(SSHclm_file),EXIST=ex)
    if (ex) then
      ! Read in the model climatology
      CALL read_etaclm(SSHclm_file,SSHclm_m)
    ELSE
      WRITE(6,*) "The model SSH climatology file does not exist: ", SSHclm_file
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
  if (.not. ex) then
    WRITE(6,*) "The file does not exist: ", gridfile 
    WRITE(6,*) "Exiting common_mom4.f90..."
    STOP(2)
  ENDIF

  if (doverbose) then
    WRITE(6,'(A)') '  >> accessing file: ', gridfile
  endif

  call check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )
  call check( NF90_INQ_VARID(ncid,'grid_x_T',varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon) )

  if (doverbose) then
    WRITE(6,*) "lon(1) = ", lon(1)
    WRITE(6,*) "lon(nlon) = ", lon(nlon)
  endif

  call check( NF90_INQ_VARID(ncid,'grid_y_T',varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat) )

  if (doverbose) then
    WRITE(6,*) "lat(1) = ", lat(1)
    WRITE(6,*) "lat(nlat) = ", lat(nlat)
  endif

  call check( NF90_INQ_VARID(ncid,'zt',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )

  if (doverbose) then
    WRITE(6,*) "lev(1) = ", lev(1)
    WRITE(6,*) "lev(nlev) = ", lev(nlev)
  endif

  !-----------------------------------------------------------------------------
  ! dx and dy
  !-----------------------------------------------------------------------------
  call check( NF90_INQ_VARID(ncid,'ds_01_21_T',varid) )    ! width of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dx) ) 
  call check( NF90_INQ_VARID(ncid,'ds_10_12_T',varid) )    ! height of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dy) ) 
  call check( NF90_INQ_VARID(ncid,'area_T',varid) )        ! area of T_cell
  call check( NF90_GET_VAR(ncid,varid,area_t) ) 

  if (doverbose) then
    WRITE(6,*) "common_mom4:: grid_spec.nc MIN(dx) = ", MINVAL(dx)
    WRITE(6,*) "common_mom4:: grid_spec.nc MAX(dx) = ", MAXVAL(dx)
    WRITE(6,*) "common_mom4:: grid_spec.nc MIN(dy) = ", MINVAL(dy)
    WRITE(6,*) "common_mom4:: grid_spec.nc MAX(dy) = ", MAXVAL(dy)
    WRITE(6,*) "common_mom4:: grid_spec.nc MIN(area_t) = ", MINVAL(area_t)
    WRITE(6,*) "common_mom4:: grid_spec.nc MAX(area_t) = ", MAXVAL(area_t)
  endif

  !
  ! kmt data
  !
  call check( NF90_INQ_VARID(ncid,'num_levels',varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,kmt0) )
  if (doverbose) then
    WRITE(6,*) "kmt0(1,1) = ", kmt0(1,1)
    WRITE(6,*) "kmt0(nlon,nlat) = ", kmt0(nlon,nlat)
  endif
  kmt = NINT(kmt0)
  call check( NF90_INQ_VARID(ncid,'wet',varid) )        ! land/sea flag (0=land) for T-cell
  call check( NF90_GET_VAR(ncid,varid,wet) )

  if (doverbose) then
    WRITE(6,*) "wet(1,1) = ", wet(1,1)
    WRITE(6,*) "wet(nlon,nlat) = ", wet(nlon,nlat)
  endif

  if (doverbose) then
    WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile
    WRITE(6,*) "dx(1,1) = ", dx(1,1)
    WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
    WRITE(6,*) "dy(1,1) = ", dy(1,1)
    WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)
  endif

  !STEVE: needed for computing the AMOC based on the streamfunction calculation:
  call check( NF90_INQ_VARID(ncid,'zb',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,zb) )

  if (doverbose) then
    WRITE(6,*) "zb(1) = ", zb(1)
    WRITE(6,*) "zb(nlev) = ", zb(nlev)
  endif

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
    stop "NetCDF read error! EXITING..."
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

  if (doverbose) then
    WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file
  endif

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


SUBROUTINE read_history(infile,v3d,v2d)
!===============================================================================
! Read in a set of netcdf-format mom4 history files
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(INOUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(INOUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: ncid,varid
  INTEGER :: i,j,k

  ! For now, just read the mixed layer depth if it is an option, otherwise return
  if (DO_MLD) then
    WRITE(6,*) "Reading mixed layer depth from history file..."

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Open the T/S netcdf restart file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
 
    if (doverbose) then
      WRITE(6,*) "read_diag:: just opened file ", infile
    endif  

    !!! mixed layer depth (mld)
    buf4=0.0
    call check( NF90_INQ_VARID(ncid,'mld',varid) )
    call check( NF90_GET_VAR(ncid,varid,buf4) )
    if (dodebug) WRITE(6,*) "read_diag:: just got data for variable mld"
    v2d(:,:,iv2d_mld) = buf4
    if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable: mld"
  
    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-MLD:"
      WRITE(6,*) "read_diag:: infile = ", infile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_mld) = ",MAXVAL(v2d(:,:,iv2d_mld))
    endif
    ! !STEVE: end

  endif 

END SUBROUTINE read_history


SUBROUTINE read_diag(infile,v3d,v2d)
!===============================================================================
! Read in a set of netcdf-format mom4 diagnostic files (for now, we're using the restarts)
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3d0(:,:,:,:)
  REAL(r_sngl), ALLOCATABLE :: v2d0(:,:,:)

  ALLOCATE(v3d0(nlon,nlat,nlev,nv3d), v2d0(nlon,nlat,nv2d))
  CALL read_restart(infile,v3d0,v2d0,1)
  v3d = REAL(v3d0,r_size)
  v2d = REAL(v2d0,r_size)

END SUBROUTINE read_diag


SUBROUTINE read_grd4(infile,v3d,v2d)
!===============================================================================
! Read subroutine for backwards compatibility
!===============================================================================
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)

  CALL read_restart(infile,v3d,v2d,1)

END SUBROUTINE read_grd4


SUBROUTINE read_restart(infile,v3d,v2d,prec_in)
!===============================================================================
! Read in a set of netcdf-format mom4 restart files in single precision
!===============================================================================
  USE netcdf
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

  
  ! Set default precision:
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

  call check( NF90_OPEN(tsfile,NF90_NOWRITE,ncid) )

  !!! t
  varname='temp'
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
  varname='salt'
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

  !-----------------------------------------------------------------------------
  ! Open file with UV data:
  !-----------------------------------------------------------------------------
  call check( NF90_OPEN(uvfile,NF90_NOWRITE,ncid) )

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
    WRITE(6,*) "read_restart:: just opened file ", sffile

    !!! SST
    varname='t_surf'
    ivid=iv2d_ssh
    call check( NF90_INQ_VARID(ncid,trim(varname),varid) ) 
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)
    
    select case(prec)
      case(1)
        buf4=0.0
        call check( NF90_GET_VAR(ncid,varid,buf4(:,:,1)) ) 
        WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = buf4(:,:,1) - t0c ! kelvin to dec C
      case(2)
        buf8=0.0d0
        call check( NF90_GET_VAR(ncid,varid,buf8(:,:,1)) ) 
        WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = REAL(buf8(:,:,1),r_sngl) - t0c ! kelvin to deg C
    end select
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid


! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-SST"
      WRITE(6,*) "read_grd4:: sffile = ", sffile
      WRITE(6,*) "max val for level v3d(:,:,iv2d_sst) = ", MAXVAL(v2d(:,:,iv2d_sst))
    endif
! !STEVE: end

    !!! SSS
    varname='s_surf'
    ivid=iv2d_ssh
    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)

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
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid

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
  call check( NF90_OPEN(sffile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart:: just read file ", sffile

  !!! SSH
  varname='sea_lev'
  ivid=iv2d_ssh
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_restart:: just got data for variable ", trim(varname)

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
  if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable ", ivid

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SSH"
    WRITE(6,*) "read_restart:: sffile = ", sffile
    WRITE(6,*) "max val for level v3d(:,:,iv2d_ssh) = ", MAXVAL(v2d(:,:,ivid))
  endif
! !STEVE: end

  call check( NF90_CLOSE(ncid) )

  altimetry : if(DO_ALTIMETRY) then
    !STEVE: use the sea level perturbation from ocean_barotropic.res.nc
    call check( NF90_OPEN(bfile,NF90_NOWRITE,ncid) )
    if (doverbose) WRITE(6,*) "read_grd4:: just opened file ", bfile

    varname='eta_t'
    ivid=iv2d_eta
    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
    if (dodebug) WRITE(6,*) "read_restart:: just got data for variable eta_t"
    if (dodebug) WRITE(6,*) "read_restart:: finished processing data for variable SSH"

    !!! SSH
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

    ! !STEVE: debug
    if (dodebug) then
      WRITE(6,*) "POST-eta"
      WRITE(6,*) "read_grd4:: bfile = ", bfile
      WRITE(6,*) "max val for level v2d(:,:,iv2d_eta) = ", MAXVAL(v2d(:,:,iv2d_eta))
    endif
    ! !STEVE: end

    call check( NF90_CLOSE(ncid) )
  endif altimetry

  mld : if (DO_MLD) then
    ! Add the ensemble of mixed-layer depths to the model state for analysis
    CALL read_history(hsfile,v3d,v2d)
  endif mld

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
    STOP(10)
  endif
! !STEVE: debug end

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

END SUBROUTINE read_restart


SUBROUTINE write_grd4(outfile,v3d_in,v2d_in)
!===============================================================================
! Write subroutine for backwards compatibility
!===============================================================================
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)

  CALL write_restart(outfile,v3d_in,v2d_in) !,1)

END SUBROUTINE write_grd4


SUBROUTINE write_restart(outfile,v3d_in,v2d_in) !,prec_in)
!===============================================================================
! Write out a set of netcdf-format restart files in single precision
!===============================================================================
  USE netcdf
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !STEVE: open temp/salt file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(tsfile,NF90_WRITE,ncid) )

  !!! t
  varname = 'temp'
  ivid = iv3d_t
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  !!! s
  varname = 'salt'
  ivid = iv3d_s
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! uv file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(uvfile,NF90_WRITE,ncid) )

  !!! u
  varname = 'u'
  ivid = iv3d_u
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  !!! v
  varname = 'v'
  ivid = iv3d_v
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,ivid)) )

  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! sfc file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(sffile,NF90_WRITE,ncid) )
  
  !!! SST
  varname = 't_surf'
  ivid = iv2d_sst
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  WHERE (kmt(:,:) .ge. 1) v2d(:,:,ivid) = v2d(:,:,ivid) + t0c ! kelvin
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )

  !!! SSS
  varname = 's_surf'
  ivid = iv2d_sss
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )

  !!! SSH
  varname = 'sea_lev'
  ivid = iv2d_ssh
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )
  call check( NF90_CLOSE(ncid) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the updated eta_t to analysis restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (DO_ALTIMETRY) then
    call check( NF90_OPEN(bfile,NF90_WRITE,ncid) )
    varname = 'eta_t'
    ivid = iv2d_eta
    call check( NF90_INQ_VARID(ncid,trim(varname),varid) )

    ! Convert SSH stored in v2d to climatological Sea Level Anomaly (SLA) by subtracting pre-computed model climatology
    v2d(:,:,ivid) = v2d(:,:,ivid) + SSHclm_m(:,:)

    call check( NF90_PUT_VAR(ncid,varid,v2d(:,:,ivid)) )
    call check( NF90_CLOSE(ncid) )
  endif

END SUBROUTINE write_restart


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
  do n=1,nv3d
    do k=1,nlev
      READ(iunit,REC=irec) buf4
      irec = irec + 1
      v3d(:,:,k,n) = REAL(buf4,r_size)
    enddo
  enddo

  do n=1,nv2d
    READ(iunit,REC=irec) buf4
    irec = irec + 1
    v2d(:,:,n) = REAL(buf4,r_size)
  enddo

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
  do n=1,nv3d
    do k=1,nlev
      if (dodebug) print *, "write_bingrd4:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    enddo
  enddo

  do n=1,nv2d
    if (dodebug) print *, "write_bingrd4:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  enddo

  CLOSE(iunit)

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
