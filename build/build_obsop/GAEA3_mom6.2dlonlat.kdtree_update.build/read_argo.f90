MODULE read_argo
!===============================================================================
! This program reads netcdf to convert profiles of temperature and salinity to 
! a format readable by letkf. The obsop.f90 program uses this module to read these
! in directly.
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
!===============================================================================

USE common,                     ONLY: r_sngl, r_size, slen
USE params_model,               ONLY: nlev
USE params_obs,                 ONLY: id_t_obs, id_s_obs
USE compute_profile_error,      ONLY: cmpTz

IMPLICIT NONE

PUBLIC :: read_argo_nc, argo_data

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n
REAL(r_size) :: se0, seF

TYPE argo_data
  REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: lev       ! grid z level
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  CHARACTER(9) :: plat      ! Platform
  CHARACTER(3) :: ptyp      ! Profile type
  CHARACTER(3) :: sid       ! Source id
  CHARACTER(1) :: qkey      ! Quality key
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
  INTEGER :: id     ! id number used in observation files to identify the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE argo_data

TYPE(argo_data), ALLOCATABLE, DIMENSION(:) :: obs_data

!! Read temp data
!infile = trim(indir1)//'/'//YYYY//MM//DD//'tmp.nc' !STEVE: sample
!typ = id_t_obs
!CALL read_argo_nc(infile,typ,obs_data,nobs)
!print *, "temp nobs = ", nobs

!! Write letkf file
!do i=1,nobs
!!STEVE: the following are required for miyoshi's letkf observation input:
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

SUBROUTINE read_argo_nc(infile,typ,obs_data,nobs)
!===============================================================================
! Read the argo profile data
!===============================================================================

USE netcdf

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
TYPE(argo_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs

! Other variables:
REAL :: err
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
CHARACTER(NF90_MAX_NAME) :: dimname
INTEGER :: grid_z_dim, len1, len2
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour
CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: ptyp, sid
CHARACTER(9), ALLOCATABLE, DIMENSION(:) :: plat
CHARACTER(1), ALLOCATABLE, DIMENSION(:) :: qkey
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: vals, stde
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: depth ! profile depths
REAL(r_size) :: val
INTEGER :: cnt, nlv
LOGICAL :: dodebug=.true.
REAL(r_size) :: missing_value=-99.0
REAL(r_sngl) :: mvc=999

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
istat = NF90_INQ_DIMID(ncid,'count',dimid)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID count failed"
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
istat = NF90_INQ_DIMID(ncid,'grid_z',dimid)
if (istat /= NF90_NOERR) then 
  print *, "NF90_INQ_DIMID grid_z failed"
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
! Read the length of string type 1
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'len1',dimid)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID dimid failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid,dimname,len1)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQUIRE_DIMENSION len1 failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the length of string type 2
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'len2',dimid)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID len2 failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid,dimname,len2)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQUIRE_DIMENSION len2 failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the depths
!-------------------------------------------------------------------------------
ALLOCATE(depth(nlv))
istat = NF90_INQ_VARID(ncid,'grid_z',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID grid_z failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,depth)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR depth failed"
  STOP
ELSE
  print *, "depth(:) = ", depth(:)
!  STOP 0
endif

!-------------------------------------------------------------------------------
! Read the longitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(xlon(cnt))
istat = NF90_INQ_VARID(ncid,'xlon',varid)  
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID xlon failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,xlon)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR xlon failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the latitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(ylat(cnt))
istat = NF90_INQ_VARID(ncid,'ylat',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID ylat failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,ylat)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR ylat failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the time as 'hour'
!-------------------------------------------------------------------------------
ALLOCATE(hour(cnt))
istat = NF90_INQ_VARID(ncid,'hour',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID hour failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,hour)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR hour failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the platform
!-------------------------------------------------------------------------------
ALLOCATE(plat(cnt))
istat = NF90_INQ_VARID(ncid,'plat',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID plat failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,plat)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR plat failed"
  STOP
else
  print *, "plat(1) = ", plat(1)
endif

!-------------------------------------------------------------------------------
! Read the profile type
!-------------------------------------------------------------------------------
ALLOCATE(ptyp(cnt))
istat = NF90_INQ_VARID(ncid,'ptyp',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID ptyp failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,ptyp)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR ptyp failed"
  STOP
else
  print *, "ptyp(1) = ", ptyp(1)
endif

!-------------------------------------------------------------------------------
! Read the Source ID
!-------------------------------------------------------------------------------
ALLOCATE(sid(cnt))
istat = NF90_INQ_VARID(ncid,'sid',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID sid failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,sid)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR sid failed"
  STOP
else
  print *, "sid(1) = ", sid(1)
endif

!-------------------------------------------------------------------------------
! Read the Quality Key
!-------------------------------------------------------------------------------
ALLOCATE(qkey(cnt))
istat = NF90_INQ_VARID(ncid,'qkey',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID qkey failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,qkey)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR qkey failed"
  STOP
else
  print *, "qkey(1) = ", qkey(1)
endif

!-------------------------------------------------------------------------------
! Read the observed profile data (temperature or salinity)
!-------------------------------------------------------------------------------
ALLOCATE(vals(nlv,cnt))
if (typ .eq. id_t_obs) then
  istat = NF90_INQ_VARID(ncid,'temp',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID temp failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,vals)
  if (dodebug) print *, "temp(:,1) = ", vals(:,1)
elseif (typ .eq. id_s_obs) then
  istat = NF90_INQ_VARID(ncid,'salt',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID salt failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,vals)
  if (dodebug) print *, "salt(:,1) = ", vals(:,1)
endif
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
print *, "Finished reading netCDF file, closing..."
istat = NF90_CLOSE(ncid)
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! STEVE: need to compute the standard error as Dave did based on vertical temperature gradient
!-------------------------------------------------------------------------------
ALLOCATE(stde(nlv,cnt))
stde=0.0 !STEVE: placeholder for now. Perhaps better to do this in the calling function (e.g. obsop_tprof.f90)
if (typ .eq. id_t_obs) then
  se0 = 1.0
  seF = 1.5
elseif (typ .eq. id_s_obs) then
  se0 = 0.05
  seF = 0.15
endif

! CONVERT netcdf data to argo_data format:
print *, "Finished reading netCDF file, formatting data..."

nobs = cnt * nlv
print *, "ALLOCATING obs_data with nobs = ", nobs
ALLOCATE(obs_data(nobs))
n = 0
do i=1,cnt
  if (dodebug) print *, "i = ", i
  CALL cmpTz(stde(:,i),se0,seF,vals(:,i),depth(:),nlv,missing_value)  !STEVE: I would prefer to call this in the calling function, but it's
                                                                        !       easier here where the data is still organized in profiles
  if (dodebug) print *, "read_argo.f90:: stde(:,i) = ", stde(:,i)
! if (dodebug) STOP(1)

  do k=1,nlev
    val = vals(k,i)
    err = stde(k,i)
    if (val < mvc .and. depth(k) < mvc .and. val > missing_value) then
      n = n+1
      if (dodebug) print *, "n,lon,lat,depth,hour,val,err,plat,ptyp,sid,qkey = ", n,xlon(i),ylat(i),depth(k),hour(i),val,err,plat(i), ptyp(i), sid(i), qkey(i)
      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = xlon(i)
      obs_data(n)%x_grd(2) = ylat(i)
      obs_data(n)%x_grd(3) = depth(k)
      obs_data(n)%hour = hour(i)
      obs_data(n)%value = val
      obs_data(n)%oerr = err
      obs_data(n)%rid = i         ! record id
      obs_data(n)%lid = k         ! level id
      obs_data(n)%plat = plat(i)
      obs_data(n)%ptyp = ptyp(i)
      obs_data(n)%sid  = sid(i)
      obs_data(n)%qkey = qkey(i)
    endif
  enddo
enddo
nobs = n
if (dodebug) print *, "nobs = ", nobs

END SUBROUTINE read_argo_nc

END MODULE read_argo
