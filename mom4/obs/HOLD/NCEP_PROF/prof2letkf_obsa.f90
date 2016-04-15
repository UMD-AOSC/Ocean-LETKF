PROGRAM prof2letkf
!===============================================================================
! PROGRAM: prof2letkf
! 
! USES: none
!
! PUBLIC TYPES:
!                 implicit none
!                 [save]
!
!                 <type declaration>
!     
! PUBLIC MEMBER FUNCTIONS:
!           <function>                     ! Description      
!
! PUBLIC DATA MEMBERS:
!           <type> :: <variable>           ! Variable description
!
! DESCRIPTION: 
!   This program converts netcdf profiles of temperature and salinity to 
!   a format readable by letkf. Eventually, letkf may read these netcdf profiles directly.
! 
! !REVISION HISTORY:
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny $
!===============================================================================
!
! Observation errors are read from the observation file
!

  IMPLICIT NONE

  INTEGER, PARAMETER :: nlev = 40
  INTEGER :: nobs, nobs0
  !STEVE: from common_obs_mom4.f90...
  INTEGER, PARAMETER :: r_sngl = kind(0.0e0)
  INTEGER, PARAMETER :: r_size = kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
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

  INTEGER :: i,j,k,n
  CHARACTER(slen) :: infile, outfile
  INTEGER :: typ
  !REAL(r_size) :: wk(6)
  REAL(r_sngl) :: wk(6)
  LOGICAL :: ex ! For checking if files exist
  INTEGER, PARAMETER :: fid=90

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
  infile = trim(indir1)//'/'//YYYY//MM//DD//'.tmpa.mom' !STEVE: sample
  nobs = 0
  typ = id_t_obs
  print *, "Calling: read_obsa for, ", typ
  CALL read_obsa(infile,obs_data,nobs,typ)
  print *, "temp nobs = ", nobs
  nobs0 = nobs

  ! Read salt data
  !infile = '20110101.sala.mom' !STEVE: sample
  infile = trim(indir2)//'/'//YYYY//MM//DD//'.sala.mom' !STEVE: sample
  INQUIRE(FILE=infile,EXIST=ex)
  if (ex) then
    typ = id_s_obs
    print *, "Calling: read_obsa for, ", typ
    CALL read_obsa(infile,obs_data,nobs,typ)
    print *, "salt nobs = ", nobs - nobs0
  else
    print *, "WARNING !!!! Salinity data file does not exist!"
    print *, "WARNING !!!! Check to make sure there is not a data file for this date!"
  endif

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

SUBROUTINE read_obsa(infile,obs_data,nobs,typ)
!===============================================================================
! Read in the godas obs data format (including obs error estimate)
!===============================================================================
  IMPLICIT NONE
  TYPE(xbt_data), INTENT(INOUT), DIMENSION(:) :: obs_data
  INTEGER, INTENT(INOUT) :: nobs
  CHARACTER(*), INTENT(IN) :: infile
  INTEGER, INTENT(IN) :: typ

  INTEGER :: nu=53
  REAL, DIMENSION(nlev) :: val, err
  INTEGER :: plti, pltf
  INTEGER :: year_a, month_a, day_a, hour_a, minute_a, kd
  INTEGER :: n,i,j,k,icnt,oidx
  REAL :: lon, lat
  ! NCEP mom4p1 40 levels:
  REAL(r_size), DIMENSION(nlev) :: grid_z = &
    (/ 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 125.0, 135.0, 145.0, 155.0, &
    165.0, 175.0, 185.0, 195.0, 205.0, 215.0, 225.0, 238.4779, 262.2945, 303.0287, &
    366.7978, 459.091, 584.6193, 747.187, 949.5881, 1193.53, 1479.588, &
    1807.187, 2174.619, 2579.091, 3016.798, 3483.029, 3972.294, 4478.478 /)

  print *, "Reading from file: ", infile

  !Read the input file
  print *, "Opening file..."
  open(nu,FILE=trim(infile),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  print *, "Reading icnt..."
  read (nu) icnt
  print *, "icnt = ", icnt
  oidx=nobs
  do n=1,icnt
    print *, "n = ", n
    val = 0.0
    err = 0.0
!   STEVE: the old version needs these:
    print *, "Reading plti, pltf..."
    read (nu) plti, pltf
    print *, "plti, pltf = ", plti, pltf
    read (nu) year_a, month_a, day_a, hour_a, minute_a, kd
!   print *, "Reading year, month, day, hour, minute, kd ..."
    print *, "year, month, day, hour, minute, kd = ", year_a, month_a, day_a, hour_a, minute_a, kd
!   print *, "Reading lon, lat"
    read (nu) lon, lat
    print *, "lon, lat = ", lon, lat
!   print *, "Reading val, err..."
    read (nu) (val(k), err(k), k=1,kd)
    print *, "val = ", val
    print *, "1/sqrt(err) = ", 1/(sqrt(err))
    print *, "val(1), err(1) = ", val(1), err(1)
    print *, "1/sqrt(err(1)) = ", 1/sqrt(err(1))

    do k=1,kd
      if (val(k) < 99) then !STEVE: if not NaN or empty value
        oidx=oidx+1
        obs_data(oidx)%typ = typ
        obs_data(oidx)%x_grd(1) = lon
        obs_data(oidx)%x_grd(2) = lat
        obs_data(oidx)%x_grd(3) = grid_z(k) !(k+1)/2)
        obs_data(oidx)%hour = hour_a*1.0d0 + (minute_a/60.0d0)
        obs_data(oidx)%value = val(k)
        ! Dave Behringer: If se is an estimate of the error, the field that's
        ! written out for the godas input is 1.0/se^2
        obs_data(oidx)%oerr = 1/sqrt(err(k))
        obs_data(oidx)%rid = n
        obs_data(oidx)%lid = (k+1)/2
      endif
    enddo
  enddo
  ! Output number of observations:
  nobs = oidx
END SUBROUTINE read_obsa

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
