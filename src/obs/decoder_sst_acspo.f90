PROGRAM decoder_sst_acspo
  USE common
  USE params_model
  USE vars_model
  USE common_oceanmodel
  !USE vars_obs
  USE common_obs_oceanmodel
  USE read_sst_acspo,     ONLY: read_acspo_nc, sst_acspo_data

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Command line inputs:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)

  !-----------------------------------------------------------------------------
  ! Obs data arrays
  !-----------------------------------------------------------------------------
  TYPE(sst_acspo_data), ALLOCATABLE :: obs_data(:)
  INTEGER :: nobs
  CHARACTER(10) :: Syyyymmddhh = "YYYYMMDDHH"
  !                 1234567890
  INTEGER :: delta_seconds = -999 ! max difference from Syyymmddhh
  CHARACTER(256) :: valid_id = "NONE"

  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  REAL(r_size), ALLOCATABLE :: obhr(:)
  INTEGER     , ALLOCATABLE :: oqc(:)

  !-----------------------------------------------------------------------------
  ! Miscellaneous
  !-----------------------------------------------------------------------------
  INTEGER :: n, i

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  !STEVE: to adjust writing to output file
  LOGICAL :: print1st = .true.

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  INTEGER :: min_quality_level=5  ! CDA

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)


  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  if (delta_seconds>0) then
      CALL read_acspo_nc(trim(obsinfile), min_quality_level, obs_data, nobs, valid_id, &
                         Syyyymmddhh, delta_seconds)
  else
      CALL read_acspo_nc(trim(obsinfile), min_quality_level, obs_data, nobs, valid_id)
  end if

  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  ALLOCATE( obhr(nobs) )

  print *, "obsop_sst.f90:: starting nobs = ", nobs
  do i=1,nobs
    elem(i) = obs_data(i)%typ
    rlon(i) = obs_data(i)%x_grd(1)
    rlat(i) = obs_data(i)%x_grd(2)
    rlev(i) = 0.d0 
    odat(i) = obs_data(i)%value
    oerr(i) = obs_data(i)%oerr
    ohx(i)  = 0
    oqc(i)  = 1
    obhr(i) = obs_data(i)%hour
  enddo
  DEALLOCATE(obs_data)

  if (print1st) then  
    i=1
    print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
    i=nobs
    print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
  endif

  call write_obs3(trim(obsoutfile),nobs,elem,rlon,rlat,rlev, &
                  odat,oerr,obhr,oqc,qcflag_in=.true.)
 
  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr )

CONTAINS

SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In obsop_sst_podaac.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-minqc')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) min_quality_level
    case('-qcyyyymmddhh')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      !read (arg2,*) min_quality_level
      Syyyymmddhh = arg2
    case('-validid')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      !read (arg2,*) min_quality_level
      valid_id = trim(arg2)
    case('-maxdt')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) delta_seconds
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM decoder_sst_acspo
