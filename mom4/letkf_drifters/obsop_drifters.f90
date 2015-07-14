PROGRAM obsop_drifters

  use letkf_drifters 
  use common_obs_mom4.f90 ! use get_nobs 

  IMPLICIT NONE
  CHARACTER(slen) :: obsinfile = 'obsin_drifters.dat'	!IN (default) observation data of drifters
  CHARACTER(slen) :: guesfile = 'drifters_out.txt'	!IN data from model space in folder "DRIFTERS"
  CHARACTER(slen) :: obsoutfile = 'obsout_drifters.dat'	!OUT (default) datafile to be passed to letkf
  
  REAL(r_size), ALLOCATABLE :: elem(:)	!elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: obid(:)	! odat
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: dxyz(:)	! x',y',or z',ohx, necessary for drifters?
  INTEGER     , ALLOCATABLE :: oqc(:)	! quality control, necessary for drifters?
  REAL(r_size), ALLOCATABLE :: otime(:) ! time
  REAL(r_size), ALLOCATABLE :: v4d_all(:,;,:)
  INTEGER :: nobs  ! total number of drifters observation
  REAL(r_size) :: ri, rj, rk
  INTEGER :: n

  !-----------------------------------------------------------------------------
  ! Initialize the common_mom4 module, and process command line options
  !-----------------------------------------------------------------------------
  CALL set_common_drifters 
  CALL process_command_line ! ??, what does this for should I write the subroutine

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  CALL get_obs(obsinfile,6,nobs) ! ??, what does nrec mean. In this case it is 6 or 8?
  ALLOCATE(  elem(nobs)  )
  ALLOCATE(  rlon(nobs)  )
  ALLOCATE(  rlat(nobs)  )
  ALLOCATE(  rlev(nobs)  )
  ALLOCATE(  obid(nobs)  ) ! odat
  ALLOCATE(  oerr(nobs)  ) ! ??, read from observation?
  ALLOCATE(  dxyz(nobs)  ) ! ohx
  ALLOCATE(   oqc(nobs)  )  
  ALLOCATE(  otime(nobs) ) ! otime
  CALL read_obs2(trim(obsinfile),nobs,elem,rlon,rlat,rlev,obid,oerr,otime)

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v4d_all(num_drifters,num_times,nv4d)  )
  CALL read_drifters(trim(guesfile),v4d_all) ! read_drifters is in letkf_drifters.f90, read data from model space.

  !-----------------------------------------------------------------------------
  ! Cycle through all observations
  !-----------------------------------------------------------------------------
  dxyz=0.0d0
  DO n=1,nobs
    !---------------------------------------------------------------------------
    ! Convert the physical coordinate to model grid coordinate (note: real, not integer)
    !---------------------------------------------------------------------------
    CALL drift2ijk(elem(n),obid(n),otime(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) 
    ! ri: drifter id index
    ! rj: time idex for interpolation
    ! rk: var (x,y,x) index if elem == iv4d_x,y,z, then rk = 1,2,3 so on

    !---------------------------------------------------------------------------
    ! Filter in the tripolar region until localization is examined in the arctic !(ISSUE)
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Filter out observations that are out of range for the grid
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Check the observation against boundaries
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! observation operator (computes H(x)) for specified member
    !---------------------------------------------------------------------------
    CALL Trans_DtoY_drift(elem(n),ri,rj,rk,v4d,dxyz(n))

  enddo ! end do n=1,nobs

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated innovations to file
  !-----------------------------------------------------------------------------
  CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,obid,oerr,dxyz,oqc,time)

SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
END SUBROUTINE process_command_line
END PROGRAM obsop_drifters
