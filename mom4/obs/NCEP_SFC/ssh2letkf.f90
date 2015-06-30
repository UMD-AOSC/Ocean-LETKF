PROGRAM ssh2letkf

USE netcdf

IMPLICIT NONE

 INTEGER, PARAMETER :: nlev = 40
 INTEGER, PARAMETER :: nlon = 720, nlat = 410
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
 INTEGER,PARAMETER :: id_eta_obs=5351    !(OCEAN)
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
CHARACTER(slen) :: indir='', outdir=''
REAL :: aid,a0,a1 !date id for altimetry data
CHARACTER(slen) :: afile !file with altimetry data

INTEGER :: i,j,k,n
CHARACTER(slen) :: infile, outfile
REAL(r_size) :: err
INTEGER :: typ
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

!REAL(r_size) :: ssh_model_clim(nlon,nlat)     !(OCEAN)
INTEGER :: ncid, varid !netcdf
LOGICAL :: dodebug = .true.


! Process the command line
CALL process_command_line
infile = trim(indir)//'/'//trim(afile)
!
! Read the SSH climatology for 1993-1999 to add to the observation data, as an
! alternative to removing from model SSH during
! assimilation. This is a GODAS convention we are temporarily adopting here
!
!WRITE(6,'(A)') '  >> accessing file: aEtaCds9399.nc'
!call check( NF90_OPEN('aEtaCds9399.nc',NF90_NOWRITE,ncid) )
!call check( NF90_INQ_VARID(ncid,'ssh',varid) )      ! depth of T-cell
!call check( NF90_GET_VAR(ncid,varid,ssh_model_clim) )
!WRITE(6,*) "ssh_model_clim(1,1) = ", ssh_model_clim(1,1)
!WRITE(6,*) "ssh_model_clim(nlon,nlat) = ", ssh_model_clim(nlon,nlat)

! Read alt data
!infile = 'j1_c022.txt'
!infile = trim(indir)//'/'//YYYY//MM//DD//'tmp.nc' !STEVE: sample
nobs = 0
err = 0.04d0  ! 4 cm !STEVE: WARNING! Change this to GODAS values!
print *, "WARNING! must change this err=0.04 to appropriate GODAS values!"
typ = id_eta_obs !id_ssh_obs
CALL read_obs_ascii(infile,obs_data,nobs,err,typ)
print *, "ssh nobs = ", nobs
nobs0 = nobs

outfile = trim(outdir)//'/'//YYYY//MM//DD//'.dat' !STEVE: sample
print *, "outdir = ", outdir
print *, "YYYY = ", YYYY
print *, "MM = ", MM
print *, "DD = ", DD
print *, "Reading file for output in letkf format: ", outfile
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
!wk(7) = obs_data(i)%hour
!wk(8) = !obs id
 if (dodebug) print *, wk(1), wk(2), wk(3), wk(4), wk(5), wk(6) !, wk(7), wk(8)

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

SUBROUTINE read_obs_ascii(infile,obs_data,nobs,err,typ)
IMPLICIT NONE
TYPE(xbt_data), INTENT(INOUT), DIMENSION(:) :: obs_data
INTEGER, INTENT(INOUT) :: nobs
CHARACTER(*), INTENT(IN) :: infile
REAL(r_size), INTENT(IN) :: err
INTEGER, INTENT(IN) :: typ

! Other variables:
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
INTEGER, PARAMETER :: maxcnt = 50000
INTEGER, PARAMETER :: maxit = 1000
INTEGER :: grid_z_dim, len1, len2
REAL(r_size), DIMENSION(maxcnt) :: xlon, ylat, ascdsc, dyssnc, ssh, samples
INTEGER, DIMENSION(maxcnt) :: track
!REAL(r_size), DIMENSION(nlev) :: grid_z
REAL(r_size) :: val
CHARACTER(64) :: dat
INTEGER :: count
REAL(r_size) :: chkaid
LOGICAL :: dodebug=.true.

! Open file
print *, "Opening altimetry infile: ", infile
open (80,file=infile)

! Fields are:
! These are cycle files in ascii.  Each line has one data point (a 1-degree average along the track): <lat> <long> <asc/desc> <day>
n = 0
do
    read (80,'(a)',end=101) dat
    read(dat(27:34),'(f7.2)') chkaid
    if (chkaid < a0 .or. a1 < chkaid) CYCLE  !Only get the correct date
    if (dodebug) print *, "a0, a1, chkaid = ", a0, a1, chkaid
    n = n+1
!   print *, "line n = ", n
!   print *, "dat = ", dat
    ! lat = dat(1:9)
    ! lon = dat(10:17)
    ! asc/desc = dat(18:26)
    ! days since 1/1/1985 = dat(27:34)
    ! ssh rel to 93-99 ave = dat(35:42)
    ! #samples = dat(43:50)
    ! track# = dat(51:54)
    ! ??? = dat(55:59)
    ! ??? = dat(60:63)
    read(dat(1:9),'(f7.2)') ylat(n)
    read(dat(10:17),'(f7.2)') xlon(n)
    read(dat(18:26),'(f7.2)') ascdsc(n)
    read(dat(27:34),'(f7.2)') dyssnc(n)
    read(dat(35:42),'(f7.2)') ssh(n)
    read(dat(43:50),'(f7.2)') samples(n)
    read(dat(51:54),'(i4)') track(n)
    
!STEVE: these lines are from another sample file:
!   if (dat(1:4) == 'time') write (2,'(a)') trim(dat)
!   i = index(dat,' = ')
!   if (i == 0) cycle
!   if (n == NREC) exit
end do
101 continue
count = n

print *, "**************************"
print *, "Record Count is = ", count
print *, "**************************"
if (maxcnt < count) then
  print *, "ERROR: maxcnt must be increased to > ", count
  stop 1
endif

! Close file
close(80)

! CONVERT ascii data to xbt_data format:
print *, "Finished reading file, formatting data..."
n=nobs
do i=1,count
  n = n+1
  !print *, "n,i,k,depth, val = ", n,i,k,grid_z(k), val
  hour = (dyssnc(i) - FLOOR(dyssnc(i))) * 24
! print *, "n,lon,lat,depth,hour,val = ", n,xlon(i),ylat(i),hour,val
  obs_data(n)%typ = typ
  obs_data(n)%x_grd(1) = xlon(i)
  obs_data(n)%x_grd(2) = ylat(i)
  obs_data(n)%x_grd(3) = 0 !SFC data
  obs_data(n)%hour = hour
  obs_data(n)%value = ssh(i)/100 !STEVE: converting from cm to m
  obs_data(n)%oerr = err
  obs_data(n)%rid = i
  obs_data(n)%lid = k
enddo
nobs = n

END SUBROUTINE read_obs_ascii

!################################################################################################
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

  RETURN
END SUBROUTINE itpl_2d

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
    case ('-aid')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) aid
      a0 = FLOOR(aid)
      a1 = FLOOR(aid+1)-0.001
      print *, "aid = ", aid
    case ('-a0')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) a0
      print *, "a0 = ", a0
    case ('-a1')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) a1
      print *, "a1 = ", a1
    case ('-afile')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      afile = arg2
      print *, "afile = ", afile
    case ('-y')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) year
      WRITE(YYYY,'(i4.4)') year
      print *, "year = ", arg2
      print *, "YYYY = ", YYYY
    case('-m')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) month
      WRITE(MM,'(i2.2)') month
      print *, "month = ", arg2
      print *, "MM = ", MM
    case('-d')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) day
      WRITE(DD,'(i2.2)') day
      print *, "day = ", arg2
      print *, "DD = ", DD
    case('-h')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) hour
      WRITE(HH,'(i2.2)') hour
      print *, "hour = ", arg2
      print *, "HH = ", HH
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
      indir = trim(arg2)
    case('-outdir')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      outdir = trim(arg2)
  end select
END DO

END SUBROUTINE process_command_line

subroutine check(status)
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

END PROGRAM ssh2letkf
