MODULE read_ecmwf_fdbk
!===============================================================================
! This program reads netcdf to convert ECMWF NEMOVAR-feedback files containing T/S
! profiles of temperature and salinity to a format readable by letkf. 
! The obsop.f90 program uses this module to read these in directly.
!
! Observation errors should be computed externally using Dave Behringer's technique based
! on the temperature gradients, or read from an observation error file
! (e.g. when doing adaptive obs error)
!
! Either the in situ temperatures will have to be converted to potential temperature,
! or the obs operator will have to transform the model potential temperature to
! in situ equivalent. The latter will be easier because there is no guarantee
! that temperature and salinity are observed simultaneously.
!
! For the ECMWF feedback files, I'm assuming the observation is converted
! to potential temperature.
!
!===============================================================================

USE common,                     ONLY: r_sngl, r_size, slen
USE params_obs,                 ONLY: id_t_obs, id_s_obs
USE compute_profile_error,      ONLY: cmpTz

IMPLICIT NONE

PUBLIC :: read_fdbk_nc, argo_data

PRIVATE

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n
REAL(r_size) :: se0, seF

TYPE argo_data
  REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: hxb       ! Model equivalent ! (NEMO) added
  REAL(r_size) :: lev       ! grid z level
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  CHARACTER(9) :: plat      ! Platform
  CHARACTER(3) :: ptyp      ! Profile type
  CHARACTER(3) :: sid       ! Source id
  CHARACTER(1) :: qkey      ! Quality key
  INTEGER :: qc             ! Quality control flag  ! (NEMO) added
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
  INTEGER :: id     ! id number used in observation files to identify the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE argo_data

LOGICAL :: DO_READ_OBE      ! Read the observation error estimate from the file
LOGICAL :: DO_COMPUTE_OBE   ! Compute the observation error estimate using the NCEP approach w/ vertical gradient

INTEGER :: undef = 99999

!! Write letkf file
!do i=1,nobs
!!STEVE: the following are required for miyoshi's letkf observation input format:
!!1 = obelm
!!2 = lon
!!3 = lat
!!4 = lev
!!5 = value
!!6 = oberr
! wk(1) = obs_data(i)%typ
! wk(2) = obs_data(i)%x_grd(1)
! wk(3) = obs_data(i)%x_grd(2)
! wk(4) = obs_data(i)%x_grd(3)
! wk(5) = obs_data(i)%value
! wk(6) = obs_data(i)%oerr
! WRITE(fid) wk
!enddo

CONTAINS

SUBROUTINE read_fdbk_nc(infile,typ,obs_data,nobs,obs_name,mod_name,opt_in)
!===============================================================================
! Read the argo profile data
! contained in the NEMOVAR netcdf feedback (fdbk) file used by ECMWF
! Option: read in observation (0), or model equivalent (1)
!===============================================================================

USE netcdf

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
TYPE(argo_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs
INTEGER, INTENT(IN), OPTIONAL :: opt_in
CHARACTER(*), INTENT(IN) :: obs_name, mod_name !STEVE: observation name and model equivalent name in netcdf file

! Other variables:
INTEGER :: opt
REAL :: err
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
CHARACTER(NF90_MAX_NAME) :: dimname
INTEGER :: grid_z_dim
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour
CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: ptyp, sid
CHARACTER(9), ALLOCATABLE, DIMENSION(:) :: plat
CHARACTER(1), ALLOCATABLE, DIMENSION(:) :: qkey
INTEGER, ALLOCATABLE, DIMENSION(:) :: qc
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: ovals, mvals, stde
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: depth ! profile depths  (dimensioned: nobs, depth)
REAL(r_size) :: val, hxb
INTEGER :: cnt, nlv
LOGICAL :: dodebug=.true.
REAL(r_size) :: missing_value=99999 ! Missing value flag in the data
REAL(r_sngl) :: max_value=999       ! Maximum value condition
REAL(r_sngl) :: max_depth = 9999    ! Maximum depth condition

CHARACTER(slen) :: tobs_name, tmod_name
CHARACTER(slen) :: sobs_name, smod_name

! Set option to identify whether we are reading the observation or the model equivalent
! from this NEMOVAR feedback file.
if (present(opt_in)) then
  opt = opt_in
else
  opt = 1
endif

if (typ .eq. id_t_obs) then
! Example possible options:
! tobs_name = 'POTM_OBS'
! tmod_name = 'POTM_Hx0'
! tmod_name = 'POTM_Hx0qc'
! tmod_name = 'POTM_Hx'
  tobs_name = trim(obs_name)
  tmod_name = trim(mod_name) 
elseif (typ .eq. id_s_obs) then
! Example possible options:
! sobs_name = 'PSAL_OBS'
! smod_name = 'PSAL_Hx0'
! smod_name = 'PSAL_Hx0qc'
! smod_name = 'PSAL_Hx'
  sobs_name = trim(obs_name)
  smod_name = trim(mod_name) 
endif

if (opt==0) then
  DO_READ_OBE = .false.
  DO_COMPUTE_OBE = .false.
elseif (opt==1) then
  ! To read the fdbk file
  DO_READ_OBE = .false.
  DO_COMPUTE_OBE = .true.
elseif (opt==2) then
  ! To read the fdbkqc file
  DO_READ_OBE = .true.
  DO_COMPUTE_OBE = .false.
endif

!-------------------------------------------------------------------------------
! Open netcdf file
!-------------------------------------------------------------------------------
istat = NF90_OPEN(infile,NF90_NOWRITE,ncid)
if (istat /= NF90_NOERR) then
  WRITE(6,'(A)') 'netCDF OPEN ERROR on ', infile
  STOP
endif

!-------------------------------------------------------------------------------
! Read the number of records
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'N_OBS',dimid)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID N_OBS failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid,dimname,cnt)
if (istat /= NF90_NOERR) then
  print *, "count NF90_INQUIRE_DIMENSION failed"
  STOP
endif
print *, "**************************"
print *, "Record Count is = ", cnt
print *, "**************************"

!-------------------------------------------------------------------------------
! Read the number of stored depths
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'N_LEVELS',dimid)
if (istat /= NF90_NOERR) then 
  print *, "NF90_INQ_DIMID N_LEVELS failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid,dimname,nlv)
if (istat /= NF90_NOERR) then
  print *, "nlv NF90_INQUIRE_DIMENSION failed"
  STOP
endif
print *, "**************************"
print *, "Number of depths is = ", nlv
print *, "**************************"

!-------------------------------------------------------------------------------
! Read the depths
!-------------------------------------------------------------------------------
ALLOCATE(depth(nlv,cnt))
istat = NF90_INQ_VARID(ncid,'DEPTH',varid)   

if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID DEPTH failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,depth)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR DEPTH failed"
  STOP
ELSE
  print *, "depth(:,:) = ", depth(:,:)
!  STOP 0
endif

!-------------------------------------------------------------------------------
! Read the longitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(xlon(cnt))
istat = NF90_INQ_VARID(ncid,'LONGITUDE',varid)  
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID LONGITUDE failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,xlon)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR LONGITUDE failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the latitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(ylat(cnt))
istat = NF90_INQ_VARID(ncid,'LATITUDE',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID LATITUDE failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,ylat)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR LATITUDE failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the time as 'hour'
!-------------------------------------------------------------------------------
ALLOCATE(hour(cnt))
istat = NF90_INQ_VARID(ncid,'JULD',varid)   
! relative julian days with decimal part (as parts of day)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID JULD failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,hour)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR JULD failed"
  STOP
endif
! Assign as an hour:
!STEVE (ISSUE) - Need the datetime also to give a unique timestamp
!                if the files contain multiple days of data
!hour = (hour - floor(hour))*24
!STEVE: for now, using JULD in place of 'hour' this can be used for temporal smoothing

!-------------------------------------------------------------------------------
! Read the platform
!-------------------------------------------------------------------------------
!ALLOCATE(plat(cnt))
!istat = NF90_INQ_VARID(ncid,'STATION_IDENTIFIER',varid)   
!if (istat /= NF90_NOERR) then
!  print *, "NF90_INQ_VARID STATION_IDENTIFIER failed"
!  STOP
!endif
!istat = NF90_GET_VAR(ncid,varid,plat)
!if (istat /= NF90_NOERR) then
!  print *, "NF90_GET_VAR plat failed"
!  STOP
!else
!  print *, "plat(1) = ", plat(1)
!endif

!-------------------------------------------------------------------------------
! Read the profile type
!-------------------------------------------------------------------------------
!ALLOCATE(ptyp(cnt))
!istat = NF90_INQ_VARID(ncid,'STATION_TYPE',varid)   
!if (istat /= NF90_NOERR) then
!  print *, "NF90_INQ_VARID STATION_TYPE failed"
!  STOP
!endif
!istat = NF90_GET_VAR(ncid,varid,ptyp)
!if (istat /= NF90_NOERR) then
!  print *, "NF90_GET_VAR ptyp failed"
!  STOP
!else
!  print *, "ptyp(1) = ", ptyp(1)
!endif

!-------------------------------------------------------------------------------
! Read the Source ID
!-------------------------------------------------------------------------------
!ALLOCATE(sid(cnt))
!istat = NF90_INQ_VARID(ncid,'ORIGINAL_FILE_INDEX',varid)   
!if (istat /= NF90_NOERR) then
!  print *, "NF90_INQ_VARID ORIGINAL_FILE_INDEX failed"
!  STOP
!endif
!istat = NF90_GET_VAR(ncid,varid,sid)
!if (istat /= NF90_NOERR) then
!  print *, "NF90_GET_VAR sid failed"
!  STOP
!else
!  print *, "sid(1) = ", sid(1)
!endif

!-------------------------------------------------------------------------------
! Read the Quality Key
!-------------------------------------------------------------------------------
ALLOCATE(qc(cnt))
if (typ .eq. id_t_obs) then
  istat = NF90_INQ_VARID(ncid,'POTM_QC',varid)   
elseif (typ .eq. id_s_obs) then
  istat = NF90_INQ_VARID(ncid,'PSAL_QC',varid)
endif
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID PXXX_QC failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,qc)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR qc failed"
  STOP
else
  print *, "qc(1) = ", qc(1)
endif

!-------------------------------------------------------------------------------
! Read the observed profile data (temperature or salinity)
!-------------------------------------------------------------------------------
ALLOCATE(ovals(nlv,cnt))
if (typ .eq. id_t_obs) then
  istat = NF90_INQ_VARID(ncid,trim(tobs_name),varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID failed for ", trim(tobs_name)
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,ovals)
  if (dodebug) print *, "temp(:,1) = ", ovals(:,1)
elseif (typ .eq. id_s_obs) then
  istat = NF90_INQ_VARID(ncid,trim(sobs_name),varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID failed for ", trim(sobs_name)
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,ovals)
  if (dodebug) print *, "salt(:,1) = ", ovals(:,1)
endif
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! Read the model equivalanet profile data (temperature or salinity)
!-------------------------------------------------------------------------------
ALLOCATE(mvals(nlv,cnt))
mvals = 0.0d0
if (opt > 0) then
  if (typ .eq. id_t_obs) then
    istat = NF90_INQ_VARID(ncid,trim(tmod_name),varid)   
    if (istat /= NF90_NOERR) then
      print *, "NF90_INQ_VARID failed for ", trim(tmod_name)
      STOP
    endif
    istat = NF90_GET_VAR(ncid,varid,mvals)
    if (dodebug) print *, "temp(:,1) = ", mvals(:,1)
  elseif (typ .eq. id_s_obs) then
    istat = NF90_INQ_VARID(ncid,trim(smod_name),varid)   
    if (istat /= NF90_NOERR) then
      print *, "NF90_INQ_VARID failed for ", trim(smod_name)
      STOP
    endif
    istat = NF90_GET_VAR(ncid,varid,mvals)
    if (dodebug) print *, "salt(:,1) = ", mvals(:,1)
  endif
  if (istat /= NF90_NOERR) STOP
endif

!-------------------------------------------------------------------------------
! Read the observed profile standard errors (temperature or salinity)
!-------------------------------------------------------------------------------
if (DO_READ_OBE) then
  ALLOCATE(stde(nlv,cnt))
  if (typ .eq. id_t_obs) then
    istat = NF90_INQ_VARID(ncid,'POTM_OBEa1qc',varid)   
    if (istat /= NF90_NOERR) then
      print *, "NF90_INQ_VARID failed for ", 'POTM_OBEa1qc'
      STOP
    endif
    istat = NF90_GET_VAR(ncid,varid,stde)
    if (dodebug) print *, "temp error(:,1) = ", stde(:,1)
  elseif (typ .eq. id_s_obs) then
    istat = NF90_INQ_VARID(ncid,'PSAL_OBEa1qc',varid)   
    if (istat /= NF90_NOERR) then
      print *, "NF90_INQ_VARID failed for ", 'PSAL_OBEa1qc'
      STOP
    endif
    istat = NF90_GET_VAR(ncid,varid,stde)
    if (dodebug) print *, "salt error(:,1) = ", stde(:,1)
  endif
  if (istat /= NF90_NOERR) STOP
endif

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
print *, "Finished reading netCDF file, closing..."
istat = NF90_CLOSE(ncid)
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! Either compute the standard error as Behringer did at NCEP based on the
! vertical temperature gradient (DO_COMPUTE_OBE), OR use the precomputed 
! standard error from the NEMOVAR observation operator (DO_READ_OBE).
!-------------------------------------------------------------------------------
if (DO_COMPUTE_OBE) then
  ALLOCATE(stde(nlv,cnt))
  stde=0.0
  if (typ .eq. id_t_obs) then
    se0 = 1.0
    seF = 1.5
  elseif (typ .eq. id_s_obs) then
    se0 = 0.05
    seF = 0.15
  endif
endif

! CONVERT netcdf data to argo_data format:
print *, "Finished reading netCDF file, formatting data..."

! Loop through all profiles and create a new observation for each unique data point
nobs = cnt * nlv
print *, "ALLOCATING obs_data with nobs = ", nobs
ALLOCATE(obs_data(nobs))
n = 0
do i=1,cnt
  if (dodebug) print *, "i = ", i
  if (DO_COMPUTE_OBE) then
    CALL cmpTz(stde(:,i),se0,seF,ovals(:,i),depth(:,i),nlv,missing_value)  
      !STEVE: I would prefer to call this in the calling function, but it's
      !       easier here where the data is still organized in profiles
  endif
  if (DO_COMPUTE_OBE .or. DO_READ_OBE) then
    if (dodebug) print *, "read_ecmwf_fdbk.f90:: stde(:,i) = ", stde(:,i)
  endif
! if (dodebug) STOP(1)

  do k=1,nlv

    ! Assign observation value:
    val = ovals(k,i)
    if (opt > 0) then
      hxb = mvals(k,i)
    else
      hxb = missing_value
    endif

    ! Assign observation error:
    if (DO_COMPUTE_OBE .or. DO_READ_OBE) then
      err = stde(k,i)
      ! NaN check:
!     if (isnan(err)) then
!       err = max(se0,seF)
!       print *, "WARNING: NaN computed for stde, using instead err = ", err
!       print *, "This is a sign of something wrong in the read_ecmwf_fdbk.f90 routine."
!       print *, "Exiting..."
!       STOP('NaN stde computed.')
!     endif
    else
      err = 0.0
    endif

    if (abs(val) < abs(max_value) .and. abs(val) < abs(missing_value-1) .and. depth(k,i) < max_depth) then
      n = n+1
      if (dodebug) print *, "n,lon,lat,depth,hour,val,err,hxb,qc = ", n,xlon(i),ylat(i),depth(k,i),hour(i),val,err,hxb,qc(i)

      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = xlon(i)
      obs_data(n)%x_grd(2) = ylat(i)
      obs_data(n)%x_grd(3) = depth(k,i)
      obs_data(n)%hour = hour(i)
      obs_data(n)%value = val
      obs_data(n)%hxb   = hxb
      obs_data(n)%oerr = err
      obs_data(n)%rid = i         ! record id
      obs_data(n)%lid = k         ! level id
!     obs_data(n)%plat = plat(i)
!     obs_data(n)%ptyp = ptyp(i)
!     obs_data(n)%sid  = sid(i)
!     obs_data(n)%qkey = qkey(i)
      obs_data(n)%qc = qc(i)
    endif
  enddo
enddo
nobs = n
if (dodebug) print *, "nobs = ", nobs

! Explicitly deallocate temporary arrays
if (ALLOCATED(depth)) DEALLOCATE(depth)
if (ALLOCATED(xlon))  DEALLOCATE(xlon)
if (ALLOCATED(ylat))  DEALLOCATE(ylat)
if (ALLOCATED(hour))  DEALLOCATE(hour)
!DEALLOCATE(plat)
!DEALLOCATE(ptyp)
!DEALLOCATE(sid)
!DEALLOCATE(qkey)
if (ALLOCATED(qc))    DEALLOCATE(qc)
if (ALLOCATED(ovals))  DEALLOCATE(ovals)
if (ALLOCATED(mvals))  DEALLOCATE(mvals)
if (ALLOCATED(stde))  DEALLOCATE(stde)

if (dodebug) print *, "Temporary arrays deallocated."
if (dodebug) print *, "Returning from read_fdbk_nc..."

END SUBROUTINE read_fdbk_nc

END MODULE read_ecmwf_fdbk
