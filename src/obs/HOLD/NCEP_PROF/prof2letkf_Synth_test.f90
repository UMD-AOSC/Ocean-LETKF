PROGRAM prof2letkf
! This program converts netcdf profiles of temperature and salinity to 
! a format readable by letkf. Eventually, letkf may read these netcdf profiles directly.
!
! Observation errors are read from the observation file
!
!USE common
!USE common_obs_mom4
!USE common_mom4

IMPLICIT NONE

INTEGER, PARAMETER :: nlev = 40
INTEGER :: nobs, nobs0
!STEVE: from common_obs_mom4.f90...
INTEGER, PARAMETER :: r_sngl = kind(0.0e0)
INTEGER, PARAMETER :: r_size = kind(0.0d0)
INTEGER,PARAMETER :: id_u_obs=2819
INTEGER,PARAMETER :: id_v_obs=2820
INTEGER,PARAMETER :: id_t_obs=3073
INTEGER,PARAMETER :: id_q_obs=3330
INTEGER,PARAMETER :: id_rh_obs=3331
INTEGER,PARAMETER :: id_ps_obs=14593
INTEGER,PARAMETER :: id_rain_obs=9999
INTEGER,PARAMETER :: id_s_obs=5521      !(OCEAN)
INTEGER,PARAMETER :: id_ssh_obs=5526    !(OCEAN)
INTEGER,PARAMETER :: id_sst_obs=5525    !(OCEAN)
INTEGER,PARAMETER :: id_sss_obs=5522    !(OCEAN)
INTEGER,PARAMETER :: id_x_obs=6622      !(OCEAN) (DRIFTERS)
INTEGER,PARAMETER :: id_y_obs=6666      !(OCEAN) (DRIFTERS)
INTEGER,PARAMETER :: id_z_obs=6611      !(OCEAN) (DRIFTERS)

INTEGER, PARAMETER :: slen = 512
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

TYPE xbt_data
  REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: lev       ! grid z level
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
  INTEGER :: id     ! id number used in observation files to identify the date of the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
  LOGICAL :: kept    ! tells letkf whether this obs is kept for assimilation
END TYPE xbt_data

INTEGER, PARAMETER :: maxcnt = 2000
TYPE(xbt_data), DIMENSION(maxcnt*nlev) :: obs_data

! Process the command line
CALL process_command_line

! Read temp data
!infile = '20110101.tmpa.mom' !STEVE: sample
infile = trim(indir1)//'/'//YYYY//MM//DD//'tmp.nc' !STEVE: sample
nobs = 0
typ = id_t_obs
CALL read_obs_nc(infile,obs_data,nobs,typ)
print *, "temp nobs = ", nobs
nobs0 = nobs

! Read salt data
!infile = '20110101.sala.mom' !STEVE: sample
infile = trim(indir2)//'/'//YYYY//MM//DD//'sal.nc' !STEVE: sample
typ = id_s_obs
CALL read_obs_nc(infile,obs_data,nobs,typ)
print *, "salt nobs = ", nobs - nobs0

! Open letkf file
!outfile = '20110101.dat'
outfile = trim(outdir)//'/'//YYYY//MM//DD//'.dat' !STEVE: sample
OPEN(90,FILE=outfile,FORM='unformatted',ACCESS='sequential')

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

 WRITE(90) wk
enddo

! Close letkf file
CLOSE(90)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_obs_godasErr(infile,obs_data,nobs,typ)
TYPE(xbt_data), INTENT(INOUT), DIMENSION(:) :: obs_data
INTEGER, INTENT(INOUT) :: nobs
CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
! Other variables:
INTEGER :: i,j,k,n
REAL(r_size), DIMENSION(nlev) :: grid_z
REAL(r_size) :: val
REAL :: err
INTEGER :: count
! For reading data:
INTEGER :: fw = 33
INTEGER :: n4 = 4, n8 = 8, n24 = 24, nb, nprf, kd
INTEGER, DIMENSION(6) :: ibuf
REAL, DIMENSION(2) :: buf
REAL, DIMENSION(2*40) :: aa

! NCEP mom4p1 40 levels:
grid_z = (/ 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 125.0, 135.0, 145.0, 155.0, &
    165.0, 175.0, 185.0, 195.0, 205.0, 215.0, 225.0, 238.4779, 262.2945, 303.0287, &
    366.7978, 459.091, 584.6193, 747.187, 949.5881, 1193.53, 1479.588, &
    1807.187, 2174.619, 2579.091, 3016.798, 3483.029, 3972.294, 4478.478 /)

! Open observation file:
OPEN(fw,FILE=infile,form='unformatted',access='STREAM')
READ(fw) n4
print *, "n4 = ", n4
READ(fw) nprf
print *, "nprf = ", nprf
READ(fw) n4
print *, "n4 = ", n4

! Read data:
! CONVERT godas M4-formatted data to xbt_data format:
n=nobs
do i=1,nprf
  ! Read data:
  READ(fw) n24
  print *, "n24 = ", n24
  READ(fw) ibuf
  print *, "ibuf = ", ibuf
  READ(fw) n24
  print *, "n24 = ", n24

  READ(fw) n8
  print *, "n8 = ", n8
  READ(fw) buf
  print *, "buf = ", buf
  READ(fw) n8
  print *, "n8 = ", n8

  READ(fw) nb
  print *, "nb = ", nb
  READ(fw) aa(1:2*ibuf(6))
  print *, "aa = ", aa(1:2*ibuf(6))
  READ(fw) nb
  print *, "nb = ", nb

  !STEVE: Loop through levels
  do k=1,2*ibuf(6),2
    if (typ .eq. id_t_obs) then
      val = aa(k)
      err = aa(k+1)**(-0.5)
    elseif (typ .eq. id_s_obs) then
      val = aa(k)
      err = aa(k+1)**(-0.5)
    endif
    if (val > -99) then
      n = n+1
      !print *, "n,i,k,depth, val = ", n,i,k,grid_z(k), val
      print *, "n,lon,lat,depth,val,err = ",n,buf(1),buf(2),grid_z((k+1)/2),val,err
      print *, "year,month,day,hour,minute,kd = ",ibuf(1:6)
      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = buf(1)
      obs_data(n)%x_grd(2) = buf(2)
      obs_data(n)%x_grd(3) = grid_z((k+1)/2)
      obs_data(n)%hour = ibuf(4)
      obs_data(n)%value = val
      obs_data(n)%oerr = err
      ! Dave Behringer: If  se is an estimate of the error, the field that's
      ! written out is 1.0/se^2
      obs_data(n)%rid = i
      obs_data(n)%lid = (k+1)/2
    endif
  enddo
enddo
nobs = n

END SUBROUTINE read_obs_godasErr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_obs_nc(infile,obs_data,nobs,typ)
IMPLICIT NONE
INCLUDE 'netcdf.inc'
TYPE(xbt_data), INTENT(INOUT), DIMENSION(:) :: obs_data
INTEGER, INTENT(INOUT) :: nobs
CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ

! Other variables:
REAL :: err
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
CHARACTER(NF_MAX_NAME) :: dimname
INTEGER, PARAMETER :: maxcnt = 1000
INTEGER, PARAMETER :: maxit = 1000
INTEGER :: grid_z_dim, len1, len2
REAL(r_size), DIMENSION(maxcnt) :: xlon, ylat, hour
CHARACTER(3), DIMENSION(maxcnt) :: ptyp, sid
CHARACTER(9), DIMENSION(maxcnt) :: plat
CHARACTER(1), DIMENSION(maxcnt) :: qkey
!REAL(r_size), DIMENSION(maxcnt,nlev) :: temp, salt
REAL(r_size), DIMENSION(nlev,maxcnt) :: temp, salt, stde
REAL(r_size), DIMENSION(nlev) :: grid_z
REAL(r_size) :: val
INTEGER :: count

!! STEVE: later, make these allocatable so 'maxcnt' isn't necessary
!REAL(r_size), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour
!CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: ptyp, sid
!CHARACTER(9), ALLOCATABLE, DIMENSION(:) :: plat
!CHARACTER(1), ALLOCATABLE, DIMENSION(:) :: qkey
!REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: temp, salt, stde

! Open netcdf file
istat = NF_OPEN(infile,NF_NOWRITE,ncid)
IF(istat /= NF_NOERR) THEN
  WRITE(6,'(A)') 'netCDF OPEN ERROR on ', infile
  STOP
END IF

!grid_z = 40 ;
!len1 = 3 ;
!len2 = 9 ;
!count = UNLIMITED ; // (429 currently)

! verify grid dimensions from grid_spec.nc file
istat = NF_INQ_DIMID(ncid,'grid_z',dimid)
IF(istat /= NF_NOERR) THEN 
  print *, "NF_INQ_DIMID grid_z failed"
  STOP
ENDIF
istat = NF_INQ_DIM(ncid,dimid,dimname,grid_z_dim)
IF(istat /= NF_NOERR) THEN
  print *, "grid_z_dim NF_INQ_DIM failed"
  STOP
ENDIF
if (grid_z_dim .ne. nlev) then
  print *, "ERROR: nlev must equal grid_z_dim => ", nlev, grid_z_dim
endif
istat = NF_INQ_DIMID(ncid,'len1',dimid)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_DIMID dimid failed"
  STOP
ENDIF
istat = NF_INQ_DIM(ncid,dimid,dimname,len1)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_DIM len1 failed"
  STOP
ENDIF
istat = NF_INQ_DIMID(ncid,'len2',dimid)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_DIMID len2 failed"
  STOP
ENDIF
istat = NF_INQ_DIM(ncid,dimid,dimname,len2)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_DIM len2 failed"
  STOP
ENDIF

! Read data
istat = NF_INQ_VARID(ncid,'grid_z',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID grid_z failed"
  STOP
ENDIF
istat = NF_GET_VAR_DOUBLE(ncid,varid,grid_z)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE grid_z failed"
  STOP
!ELSE
!  print *, "grid_z = ", grid_z
!  STOP 0
ENDIF
istat = NF_INQ_VARID(ncid,'xlon',varid)  
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID xlon failed"
  STOP
ENDIF
istat = NF_GET_VAR_DOUBLE(ncid,varid,xlon)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE xlon failed"
  STOP
ENDIF
istat = NF_INQ_VARID(ncid,'ylat',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID ylat failed"
  STOP
ENDIF
istat = NF_GET_VAR_DOUBLE(ncid,varid,ylat)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE ylat failed"
  STOP
ENDIF
istat = NF_INQ_VARID(ncid,'hour',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID hour failed"
  STOP
ENDIF
istat = NF_GET_VAR_DOUBLE(ncid,varid,hour)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE hour failed"
  STOP
ENDIF
istat = NF_INQ_VARID(ncid,'plat',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID plat failed"
  STOP
ENDIF
istat = NF_GET_VAR_TEXT(ncid,varid,plat)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE plat failed"
  STOP
ENDIF

istat = NF_INQ_VARID(ncid,'ptyp',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID ptyp failed"
  STOP
ENDIF
istat = NF_GET_VAR_TEXT(ncid,varid,ptyp)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE ptyp failed"
  STOP
ENDIF
istat = NF_INQ_VARID(ncid,'sid',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID sid failed"
  STOP
ENDIF
istat = NF_GET_VAR_TEXT(ncid,varid,sid)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE sid failed"
  STOP
ENDIF
istat = NF_INQ_VARID(ncid,'qkey',varid)   
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID qkey failed"
  STOP
ENDIF
istat = NF_GET_VAR_TEXT(ncid,varid,qkey)
IF(istat /= NF_NOERR) THEN
  print *, "NF_GET_VAR_DOUBLE qkey failed"
  STOP
ENDIF

if (typ .eq. id_t_obs) then
  istat = NF_INQ_VARID(ncid,'temp',varid)   
  IF(istat /= NF_NOERR) THEN
    print *, "NF_INQ_VARID temp failed"
    STOP
  ENDIF
  istat = NF_GET_VAR_DOUBLE(ncid,varid,temp)
elseif (typ .eq. id_s_obs) then
  istat = NF_INQ_VARID(ncid,'salt',varid)   
  IF(istat /= NF_NOERR) THEN
    print *, "NF_INQ_VARID salt failed"
    STOP
  ENDIF
  istat = NF_GET_VAR_DOUBLE(ncid,varid,salt)
endif
IF(istat /= NF_NOERR) STOP

istat = NF_INQ_VARID(ncid,'stde',varid)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_VARID stde failed"
  STOP
ENDIF
istat = NF_GET_VAR_DOUBLE(ncid,varid,stde)
IF(istat /= NF_NOERR) STOP

!STEVE: need to read the record length:
!* get ID of record dimension (among other things)
!      CALL NCINQ(NCID, NDIMS, NVARS, NGATTS, RECID, RCODE)
!* get record dimension name and current size
!      CALL NCDINQ(NCID, RECID, RECNAME, NRECS, RCODE)
istat = NF_INQ_DIMID(ncid,'count',dimid)
IF(istat /= NF_NOERR) THEN
  print *, "NF_INQ_DIMID count failed"
  STOP
ENDIF
istat = NF_INQ_DIM(ncid,dimid,dimname,count)
IF(istat /= NF_NOERR) THEN
  print *, "count NF_INQ_DIM failed"
  STOP
ENDIF
print *, "**************************"
print *, "Record Count is = ", count
print *, "**************************"
if (maxcnt < count) then
  print *, "ERROR: maxcnt must be increased to > ", count
  stop 1
endif


! Close netcdf file
istat = NF_CLOSE(ncid)
IF(istat /= NF_NOERR) STOP

! CONVERT netcdf data to xbt_data format:
print *, "Finished reading netCDF file, formatting data..."
n=nobs
!print *, "temp = ", temp
!do i=1,count
!  print *, "temp(:,",i,") = ", temp(:,i)
!enddo
!STOP
do i=1,count
  !STEVE: need to fix this line to match it with the <temporary> time dimension (which is UNLIMITED)
  !if (i > 429) EXIT
  !if (i > count) EXIT
  do k=1,nlev
    if (typ .eq. id_t_obs) then
      val = temp(k,i)
      err = stde(k,i) !t_eprof(k)
    elseif (typ .eq. id_s_obs) then
      val = salt(k,i)
      err = stde(k,i) !s_eprof(k)
    endif
    if (val > -99) then
      n = n+1
      !print *, "n,i,k,depth, val = ", n,i,k,grid_z(k), val
      print *, "n,lon,lat,depth,hour,val,err = ", n,xlon(i),ylat(i),grid_z(k),hour(i),val,err
      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = xlon(i)
      obs_data(n)%x_grd(2) = ylat(i)
      obs_data(n)%x_grd(3) = grid_z(k)
      obs_data(n)%hour = hour(i)
      obs_data(n)%value = val
      obs_data(n)%oerr = err
      obs_data(n)%rid = i
      obs_data(n)%lid = k
    endif
  enddo
enddo
nobs = n

END SUBROUTINE read_obs_nc

!################################################################################################
SUBROUTINE process_command_line

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

END PROGRAM prof2letkf
