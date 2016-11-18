PROGRAM 3dvar_main
!==============================================================================
!
! [PURPOSE:] Main program for executing 3dvar
!
!===============================================================================
!
! DESCRIPTION: 
!   This is the main program for the 3dvar data assimilation.
!
! REVISION HISTORY:
!   09/03/2016 Steve Penny created
!
!-------------------------------------------------------------------------------
! $Authors: Steve Penny
!===============================================================================

  USE common
  USE common_mpi
  USE common_oceanmodel
  USE common_mpi_oceanmodel
  USE params_3dvar
  USE params_model
  USE params_obs
  USE input_nml_oceanmodel, ONLY: read_input_namelist

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)     :: gues2d, anal2d
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:)   :: gues3d, anal3d
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: gues4d, anal4d
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(9) :: stdoutf='NOUT-0000'
  CHARACTER(4) :: guesf='gs00'
  
  INTEGER :: ij, m !STEVE: for debugging
  LOGICAL :: dortout=.true.    ! Force 'realtime' output (helps with parallel debugging)
  LOGICAL :: dodebug0=.false.  ! Debug flag for various routines

!------------------------------------------------------------------------------
! Initial settings
!------------------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  WRITE(stdoutf(6:9), '(I4.4)') myrank
  WRITE(6,'(3A,I4.4)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  !STEVE: if it halts here, it probably means the nlon, nlat and nlev in common_oceanmodel
  !       have not been set properly for this model's grid
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I4.4,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

  !----------------------------------------------------------------------------
  ! Read in namelist parameters
  !----------------------------------------------------------------------------
  CALL read_input_namelist !STEVE: this has been moved to input_nml_{oceanmodel}.f90 since it needed to be slightly different for each model

  !----------------------------------------------------------------------------
  ! Print header
  !----------------------------------------------------------------------------
  WRITE(6,'(A)') '================================================='
  WRITE(6,'(A)') ' 3DVar DA using preconditioned conjugate gradient'
  WRITE(6,'(A)') ' 3DVar DA using domain decomposition             '
  WRITE(6,'(A)') ' 3DVar DA using multigrid methods                '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') ' 33333   DDDD    VV     VV    A     RRRR         '
  WRITE(6,'(A)') '      3  D   DD   VV   VV    A A    R   R        '
  WRITE(6,'(A)') '   333   D    D    VV VV     AAA    RRRR         '
  WRITE(6,'(A)') '      3  D   DD     VVV     A   A   R  RR        '
  WRITE(6,'(A)') ' 33333   DDDD        V     A     A  R   R        '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') '  Developed for NCEP use by Steve Penny (2016)   '
  WRITE(6,'(A)') '                                                 '
  WRITE(6,'(A)') '================================================='
  WRITE(6,'(A)') '              3DVar PARAMETERS                   '
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
  CALL set_common_oceanmodel
  CALL set_common_mpi_oceanmodel

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
  CALL set_3dvar_obs 

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
  ! Read background error covariance information
  !-----------------------------------------------------------------------------

  ! Read in EOFs 

  ! Read in variance magnitudes for all variables

  ! read in error covariance across variables

  !-----------------------------------------------------------------------------
  ! Read climatologcal forecast ensemble
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
  CALL das_3dvar(gues3d,gues2d,anal3d,anal2d)

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

END PROGRAM 3dvar_main
