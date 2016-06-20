PROGRAM prof2letkf_rawnc
!===============================================================================
! This program converts netcdf profiles of temperature and salinity to 
! a format readable by letkf. Eventually, the obsop.f90 should read these
! in directly.
!
! Observation errors are computed here using Dave Behringer's technique based
! on the temperature gradients, or read from an observation error file
! (e.g. when doing adaptive obs error)
!
! Either the in situ temperatures will have to be converted to potential temperature,
! or the obs operator will have to transform the model potential temperature to
! in situ equivalent. The latter will be easier because there is no guarantee
! that temperature and salinity are observed simultaneously.
!
!===============================================================================

USE common,       ONLY: r_sngl, r_size, slen
USE params_model, ONLY: nlev
USE params_obs,   ONLY: id_t_obs, id_s_obs

IMPLICIT NONE

INTEGER :: nobs, nobs0
INTEGER :: year,month,day,hour
INTEGER :: start_year, end_year
INTEGER :: start_day, end_day
INTEGER :: start_month, end_month
INTEGER :: start_hour, end_hour
CHARACTER(4) :: YYYY
CHARACTER(2) :: MM, DD, HH
CHARACTER(slen) :: indir, outdir, indir1, indir2
CHARACTER(3) :: obtype = 'tmp' !tmp, sal, alt

INTEGER :: i,j,k,n
CHARACTER(slen) :: infile, outfile
INTEGER :: typ
!REAL(r_size) :: wk(6)
REAL(r_sngl) :: wk(6)
INTEGER, PARAMETER :: fid = 90

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

INTEGER, PARAMETER :: maxcnt = 2000
TYPE(argo_data), DIMENSION(maxcnt*nlev) :: obs_data

! NCEP mom4p1 40 levels:
REAL(r_size), DIMENSION(nlev) :: grid_z
grid_z = (/ 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 125.0, 135.0, 145.0, 155.0, &
    165.0, 175.0, 185.0, 195.0, 205.0, 215.0, 225.0, 238.4779, 262.2945, 303.0287, &
    366.7978, 459.091, 584.6193, 747.187, 949.5881, 1193.53, 1479.588, &
    1807.187, 2174.619, 2579.091, 3016.798, 3483.029, 3972.294, 4478.478 /)

! Process the command line
CALL process_command_line

! Read temp data
!infile = '20110101.tmpa.mom' !STEVE: sample
infile = trim(indir1)//'/'//YYYY//MM//DD//'tmp.nc' !STEVE: sample
nobs = 0
typ = id_t_obs
CALL read_argo_nc(infile,obs_data,nobs,typ)
print *, "temp nobs = ", nobs
nobs0 = nobs

! Read salt data
!infile = '20110101.sala.mom' !STEVE: sample
infile = trim(indir2)//'/'//YYYY//MM//DD//'sal.nc' !STEVE: sample
typ = id_s_obs
CALL read_argo_nc(infile,obs_data,nobs,typ)
print *, "salt nobs = ", nobs - nobs0

! Open letkf file
!outfile = '20110101.dat'
outfile = trim(outdir)//'/'//YYYY//MM//DD//'.dat' !STEVE: sample
OPEN(fid,FILE=outfile,FORM='unformatted',ACCESS='sequential')

! Write letkf file
do i=1,nobs
!STEVE: the following are required for miyoshi's letkf observation input:
!1 = obelm
!2 = lon
!3 = lat
!4 = lev
!5 = value
!6 = oberr
 wk(1) = obs_data(i)%typ
 wk(2) = obs_data(i)%x_grd(1)
 wk(3) = obs_data(i)%x_grd(2)
 wk(4) = obs_data(i)%x_grd(3)
 wk(5) = obs_data(i)%value
 wk(6) = obs_data(i)%oerr

 WRITE(fid) wk
enddo

! Close letkf file
CLOSE(fid)

CONTAINS


SUBROUTINE read_argo_nc(infile,obs_data,nobs,typ)
!===============================================================================
! Read the argo profile data
!===============================================================================

USE netcdf

IMPLICIT NONE

TYPE(argo_data), INTENT(INOUT), DIMENSION(:) :: obs_data
INTEGER, INTENT(INOUT) :: nobs
CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ

! Other variables:
REAL :: err
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
CHARACTER(NF90_MAX_NAME) :: dimname
INTEGER :: grid_z_dim, len1, len2
!! STEVE: later, make these allocatable so 'maxcnt' isn't necessary
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour
CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: ptyp, sid
CHARACTER(9), ALLOCATABLE, DIMENSION(:) :: plat
CHARACTER(1), ALLOCATABLE, DIMENSION(:) :: qkey
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: temp, salt, stde
!REAL(r_size), DIMENSION(nlev,maxcnt) :: temp, salt, stde
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: depth ! profile depths
REAL(r_size) :: val
INTEGER :: cnt, nlv
LOGICAL :: dodebug=.true.
REAL(r_size) :: missing_value=9.969209968386869E+036
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
istat = NF90_INQ_DIMID(ncid,'depth',dimid)
if (istat /= NF90_NOERR) then 
  print *, "NF90_INQ_DIMID depth failed"
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
ALLOCATE(depth(nlv,cnt))
istat = NF90_INQ_VARID(ncid,'depth',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID depth failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,depth)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR depth failed"
  STOP
ELSE
  print *, "depth(:,1) = ", depth(:,1)
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
if (typ .eq. id_t_obs) then
  ALLOCATE(temp(nlv,cnt))
  istat = NF90_INQ_VARID(ncid,'temp',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID temp failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,temp)
  if (dodebug) print *, "temp(:,1) = ", temp(:,1)
elseif (typ .eq. id_s_obs) then
  ALLOCATE(salt(nlv,cnt))
  istat = NF90_INQ_VARID(ncid,'salt',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID salt failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,salt)
  if (dodebug) print *, "salt(:,1) = ", salt(:,1)
endif
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
istat = NF90_CLOSE(ncid)
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! STEVE: need to compute the standard error as Dave did based on vertical temperature gradient
!-------------------------------------------------------------------------------
ALLOCATE(stde(nlv,cnt))
stde=0.0 !STEVE: placeholder for now.


! CONVERT netcdf data to argo_data format:
print *, "Finished reading netCDF file, formatting data..."
n=nobs
!print *, "temp = ", temp
!do i=1,cnt
!  print *, "temp(:,",i,") = ", temp(:,i)
!enddo
!STOP
do i=1,cnt
  do k=1,nlev
    if (typ .eq. id_t_obs) then
      val = temp(k,i)
      err = stde(k,i) !t_eprof(k)
    elseif (typ .eq. id_s_obs) then
      val = salt(k,i)
      err = stde(k,i) !s_eprof(k)
    endif
    if (val < mvc .and. depth(k,i) < mvc) then
      n = n+1
      !print *, "n,i,k,depth, val = ", n,i,k,grid_z(k), val
      print *, "n,lon,lat,depth,hour,val,err,plat,ptyp,sid,qkey = ", n,xlon(i),ylat(i),depth(k,i),hour(i),val,err,plat(i), ptyp(i), sid(i), qkey(i)
      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = xlon(i)
      obs_data(n)%x_grd(2) = ylat(i)
      obs_data(n)%x_grd(3) = depth(k,i)
      obs_data(n)%hour = hour(i)
      obs_data(n)%value = val
      obs_data(n)%oerr = err
      obs_data(n)%rid = i
      obs_data(n)%lid = k
      obs_data(n)%plat = plat(i)
      obs_data(n)%ptyp = ptyp(i)
      obs_data(n)%sid  = sid(i)
      obs_data(n)%qkey = qkey(i)
    endif
  enddo
enddo
nobs = n

END SUBROUTINE read_argo_nc


SUBROUTINE process_command_line
!===============================================================================
! Proces the command line arguments
!===============================================================================

IMPLICIT NONE

CHARACTER(slen) :: arg1,arg2
INTEGER :: ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!^M
! inputs are in the format "-x xxx"^M
DO i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "Argument ", i, " = ",TRIM(arg1)
 
  select case (arg1)
    case ('-ob')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obtype = arg2
      print *, "obtype = ", obtype
    case ('-y')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) year
      YYYY = arg2
      print *, "year = ", arg2
    case('-m')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) month
      MM = arg2
      print *, "month = ", arg2
    case('-d')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) day
      DD = arg2
      print *, "day = ", arg2
    case('-h')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) hour
      HH = arg2
      print *, "hour = ", arg2
    case ('-sy')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) start_year
      print *, "start_year = ", start_year
    case('-sm')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) start_month
      print *, "start_month = ", start_month
    case('-sd')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) start_day
      print *, "start_day = ", start_day
    case('-sh')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) start_hour
      print *, "start_hour = ", start_hour
    case('-ey')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) end_year
      print *, "end_year = ", end_year
    case('-em')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) end_month
      print *, "end_month = ", end_month
    case('-ed')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) end_day
      print *, "end_day = ", end_day
    case('-eh')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) end_hour
      print *, "end_hour = ", end_hour
    case('-indir')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      indir1 = trim(arg2)
    case('-indir1')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      indir1 = trim(arg2)
    case('-indir2')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      indir2 = trim(arg2)
    case('-outdir')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      outdir = trim(arg2)
  end select
END DO

END SUBROUTINE process_command_line

END PROGRAM prof2letkf_rawnc
