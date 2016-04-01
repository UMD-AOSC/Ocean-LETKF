PROGRAM letkf
!==============================================================================
!
! [PURPOSE:] Main program of LETKF
!
!===============================================================================
! MODULE: letkf_local
! 
! USES:
!  use common
!  use common_mpi
!  use common_mom4
!  use common_mpi_mom4
!  use common_letkf
!  use letkf_obs
!  use letkf_tools
!  use params_letkf
!  use params_model
!  use params_obs
!
! DESCRIPTION: 
!   This is the main program for the letkf data assimilation.
!
! REVISION HISTORY:
!   01/16/2009 Takemasa Miyoshi created for atmospheric analysis
!   04/26/2011 Steve Penny converted to OCEAN for use with mom4
!   03/18/2014 Steve Penny adapted to use on Gaea at NCEP/GFDL
!   07/08/2015 uncomment all the drifters subroutine
!
!-------------------------------------------------------------------------------
! $Authors: Steve Penny, Takemasa Miyoshi $
!===============================================================================

  USE common
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE letkf_obs
  USE letkf_tools
  USE params_letkf
  USE params_model
  USE params_obs

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: gues4d(:,:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal4d(:,:,:,:,:)
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(9) :: stdoutf='NOUT-0000'
  CHARACTER(4) :: guesf='gs00'
  
  INTEGER :: ij, m !STEVE: for debugging
  LOGICAL :: ex
  INTEGER :: fid=21
  LOGICAL :: dortout=.true.    ! Force 'realtime' output (helps with parallel debugging)
  LOGICAL :: dodebug0=.false.  ! Debug flag for various routines

  NAMELIST /params_model_nml/ gridfile, &  ! MOM4 grid_spec.nc file 
                              SSHclm_file  ! model ssh climatology for altimetry assimilation
  NAMELIST /params_obs_nml/   obs1nrec, &  ! number of records in obs.dat type file
                              obs2nrec     ! number of records in obs2.dat type file
  NAMELIST /params_letkf_nml/ nslots, &              ! Number of time slots for 4D assimilation
                              nbslot, &              ! Index of base timeslot (time at which to form analysis)
                              sigma_obs, &           ! Sigma-radius (half-width) for horizontal localization at the equator (m)
                              sigma_obs0, &          ! Sigma-radius for horizontal localization at the poles (m)
                              sigma_obsv, &          ! Sigma-radius for vertical localization (m)
                              sigma_obst, &          ! Sigma-radius for temporal localization (not activated)
                              gross_error, &         ! number of standard deviations for quality control (all outside removed)
                              DO_DRIFTERS, &         ! logical flag to do lagrangian drifters assimilation
                              DO_ALTIMETRY, &        ! logical flag to do altimetry data assimilation
                              DO_SLA, &              ! logical flag to use SLA for altimetry
                              DO_ADT, &              ! logical flag to use Absolute Dynamic Topography (ADT) for altimetry
                              DO_NO_VERT_LOC, &      ! logical flag to skip all vertical localization and project weights (default)
                              DO_MLD, &              ! logical flag to use 2-layer vertical localization: (1) SFC and mixed layer, (2) below mixed layer
                              DO_REMOVE_65N, &       ! option to remove points above 65ÂN in letkf instead of in observation operators
                              localization_method, & ! localization method to be used in letkf_local.f90
                              cov_infl_mul, &        ! multiplicative inflation factor (default=1.0, i.e. none)
                              sp_infl_add            ! additive inflation factor (default none)

!------------------------------------------------------------------------------
! Initial settings
!------------------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  WRITE(stdoutf(6:9), '(I4.4)') myrank
  WRITE(6,'(3A,I4.4)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  !STEVE: if it halts here, it probably means the nlon, nlat and nlev in common_mom4
  !       have not been set properly for this model's grid
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I4.4,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

  !----------------------------------------------------------------------------
  ! Read in namelist parameters
  !----------------------------------------------------------------------------
  INQUIRE(FILE="input.nml", EXIST=ex)
  if (ex) then
    open(fid,file="input.nml", status='OLD') !, delim='APOSTROPHE')
    read(fid,nml=params_model_nml)
    read(fid,nml=params_obs_nml)
    read(fid,nml=params_letkf_nml)
  endif
  if (DO_ADT .or. DO_SLA) DO_ALTIMETRY=.true.
  WRITE(6,*) "================================================================="
  WRITE(6,*) "Namelist inputs:"
  WRITE(6,*) "================================================================="
  write(6,params_model_nml)
  write(6,params_obs_nml)
  write(6,params_letkf_nml)
  WRITE(6,*) "================================================================="

  !----------------------------------------------------------------------------
  ! Print header
  !----------------------------------------------------------------------------
  WRITE(6,'(A)') '================================================='
  WRITE(6,'(A)') '  THE LOCAL ENSEMBLE TRANSFORM KALMAN FILTER     '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF        '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF            '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF         '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF            '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF            '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') '  Developed for NCEP use by Steve Penny (2014)   '
  WRITE(6,'(A)') '  Developed for the OCEAN by Steve Penny (2011)  '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') '  Based on original code by T. Miyoshi           '
  WRITE(6,'(A)') '  for the SPEEDY Atmospheric GCM                 '
  WRITE(6,'(A)') '  and algorithms by Ott (2004) and Hunt (2007)   '
  WRITE(6,'(A)') '================================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS                   '
  WRITE(6,'(A)') ' ----------------------------------------------- '
  WRITE(6,'(A,I15)')   '   nbv        :',nbv
  WRITE(6,'(A,I15)')   '   nslots     :',nslots
  WRITE(6,'(A,I15)')   '   nbslot     :',nbslot
  WRITE(6,'(A,F15.2)') '   sigma_obs  :',sigma_obs
  WRITE(6,'(A,F15.2)') '   sigma_obs0 :',sigma_obs0
  WRITE(6,'(A,F15.2)') '   sigma_obsv :',sigma_obsv
  WRITE(6,'(A,F15.2)') '   sigma_obst :',sigma_obst
  WRITE(6,'(A)') '================================================='

  !-----------------------------------------------------------------------------
  ! Initialize modules
  !-----------------------------------------------------------------------------
  CALL set_common_mom4
  CALL set_common_mpi_mom4

  !-----------------------------------------------------------------------------
  ! Allocate dynamic arrays
  !-----------------------------------------------------------------------------
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(anal2d(nij1,nbv,nv2d))

  !-----------------------------------------------------------------------------
  ! Check timer for initialization
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !-----------------------------------------------------------------------------
  ! Observations
  !-----------------------------------------------------------------------------
  CALL set_letkf_obs 
  !STEVE: calls read_grd, then read_ens_mpi calls read_grd4 below. This may be an inefficiency
  !STEVE: Perhaps call once here then output v3d and v2d for use as gues3d and gues2d below.

  !-----------------------------------------------------------------------------
  ! Check timer for initializing observations
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !-----------------------------------------------------------------------------
  ! Read forecast ensemble
  !-----------------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(guesf(3:4),'(I2.2)') nbslot
  WRITE(6,*) "From letkf.f90, calling read_ens_mpi..."
  CALL read_ens_mpi(guesf,nbv,gues3d,gues2d)
  WRITE(6,*) "From letkf.f90, finished calling read_ens_mpi..."

  if (dodebug0) then ! Test processing of forecast ensemble, write, and quit
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    WRITE(6,*) "From letkf.f90, calling write_ens_mpi..."
    CALL write_ens_mpi('anal',nbv,gues3d,gues2d)
    WRITE(6,*) "From letkf.f90, finished calling write_ens_mpi..."
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    STOP 4
  endif

  !-----------------------------------------------------------------------------
  ! Check timer for reading forecast ensemble
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_ENS_MPI):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !-----------------------------------------------------------------------------
  ! Write ensemble mean and spread
  !-----------------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('gues',nbv,gues3d,gues2d)

  !STEVE: debug
  if (dodebug0) CALL write_ens_mpi_grd('test',1,gues3d,gues2d)
  
  !-----------------------------------------------------------------------------
  ! Check timer for writing forecast ensemble
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_GUES):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !------------------------------------------------------------------------------
  ! Data Assimilation (MAIN)
  !------------------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL das_letkf(gues3d,gues2d,anal3d,anal2d)

  !-----------------------------------------------------------------------------
  ! Check timer for computing letkf analysis
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_LETKF):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !--
  ! Could call 3DVar here for hybrid...
  !--

  !----------------------------------------------------------------------------
  ! Write analysis ensemble
  !----------------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ens_mpi('anal',nbv,anal3d,anal2d)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('anal',nbv,anal3d,anal2d)

  !-----------------------------------------------------------------------------
  ! Check timer for writing analysis ensemble, mean, and spread
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

  !----------------------------------------------------------------------------
  ! Finalize the MPI
  !----------------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

  !-----------------------------------------------------------------------------
  ! Check timer for total runtime
  !-----------------------------------------------------------------------------
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(FINAL):',rtimer,rtimer-rtimer00
  if (dortout) then !STEVE: force output to file
    CLOSE(6)
    OPEN(6,FILE=stdoutf,POSITION='APPEND',STATUS = 'OLD')
  endif
  rtimer00=rtimer

END PROGRAM letkf
