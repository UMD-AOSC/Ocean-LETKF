PROGRAM obsop_adt
!===============================================================================
! PROGRAM: obsop_adt
! 
! USES:
!  use common
!  use params_model
!  use vars_model
!  use common_mom4
!  use params_obs
!  use vars_obs
!  use common_obs_mom4
!  use params_letkf,     ONLY: DO_ALTIMETRY, DO_ADT
!
!
! DESCRIPTION: 
!   This program acts as the observation operator for the Altimetry 
!   Absolute Dynamic Topography (ADT). It inputs observations of ADT (yo)
!   and a single forecast (xf) and computes the H(xf) used to form the 
!   innovations (yo-H(xf)) associated with that member.
!   NOTE: The innovations are computed during the execution of letkf.f90
!
! USAGE:
!   A separate instance is run independently for each member and timeslot
!   All observations are typically preprocessed to the letkf obs format, then
!   here they are read in and converted to the letkf obs2 format w/ H(x) data 
!   for each member added in a new column.
! 
! !REVISION HISTORY:
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
!   04/03/2013 Takemasa Miyoshi created for SPEEDY atmospheric model.
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny, Takemasa Miyoshi $
!===============================================================================

  USE common
  USE params_model
  USE vars_model
  USE common_mom4
  USE params_obs,                ONLY: nobs, id_t_obs, id_s_obs, id_u_obs, id_v_obs, id_eta_obs, DO_REMOVE_65N
  USE vars_obs
  USE common_obs_mom4
  USE params_letkf,              ONLY: DO_ALTIMETRY, DO_ADT

  IMPLICIT NONE

  CHARACTER(slen) :: obsinfile='obsin.dat'    !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  CHARACTER(slen) :: aoerinfile='oer_inp.grd' !IN (default)
  CHARACTER(slen) :: aoeroutfile='oer_inp.grd'!OUT(default)
  
  
  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  REAL(r_size), ALLOCATABLE :: obhr(:)
  INTEGER     , ALLOCATABLE :: oqc(:)
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size) :: dk,tg,qg
  REAL(r_size) :: ri,rj,rk
  INTEGER :: n

  REAL(r_size) :: obserr_scaling=1.0d0 !STEVE: use this to scale the input observations
  REAL(r_size) :: obs_randselect=1.0d0 !STEVE: use this to select random input observations
  REAL(r_size), DIMENSION(1) :: rand

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  INTEGER :: bdyobs=2                 !STEVE: use of boundary obs.
                                      !       1 := less restrictive, remove obs inside boundary
                                      !       2 := remove all observations touching a boundary
  LOGICAL :: debug_obsfilter = .false.
  LOGICAL :: debug_hdxf_0 = .false.   ! This error occured because there was not a model representation 
                                      ! of the observed value (i.e. SST obs with no SST model field).
                                      ! Solution was to populate a SST model field (v2d) with 
                                      ! surface temp data from the model (v3d(:,:,1)).
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
! LOGICAL :: dodebug = .false.

  INTEGER :: cnt_obs_u=0, cnt_obs_v=0, cnt_obs_t=0, cnt_obs_s=0
  INTEGER :: cnt_obs_ssh=0, cnt_obs_sst=0, cnt_obs_sss=0, cnt_obs_eta=0
  INTEGER :: cnt_obs_x=0, cnt_obs_y=0, cnt_obs_z=0
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0, cnt_oerlt0=0, cnt_altlatrng=0
  !STEVE: for adaptive obs:
  LOGICAL :: oerfile_exists
  REAL(r_size) :: oberr

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  ! For ADT assimilation:
  REAL(r_size) :: offset  = 0.0d0 !GUILLAUME:ADT
  REAL(r_size) :: cnt_adt = 0.0d0 !GUILLAUME:ADT
                                  !NOTE: while the 'cnt_' naming convention used above
                                  !      is applied for counting removed observations,
                                  !      the cnt_adt variable counts kept observations.

  !-----------------------------------------------------------------------------
  ! Initialize the common_mom4 module, and process command line options
  !-----------------------------------------------------------------------------
  CALL set_common_mom4
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  CALL get_nobs(obsinfile,6,nobs)
  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  ALLOCATE( obhr(nobs) )
  if (obs2nrec==8) then
    CALL read_obs(trim(obsinfile),nobs,elem,rlon,rlat,rlev,odat,oerr)
  elseif (obs2nrec==9) then
    CALL read_obs(trim(obsinfile),nobs,elem,rlon,rlat,rlev,odat,oerr,obhr)
  else
    WRITE(6,*) "obsop.f90:: no read_obs option for obs2nrec = ", obs2nrec
    STOP 95
  endif

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  CALL read_diag(guesfile,v3d,v2d)
  WRITE(6,*) '****************'

  !-----------------------------------------------------------------------------
  ! Cycle through all observations
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Cycling through observations..."

  !GUILLAUME:ADT: initialize adt counter
  if (DO_ADT) then 
    offset=0.0d0  
    cnt_adt=0.0d0
  endif

  ohx=0.0d0
  oqc=0
  DO n=1,nobs
    !---------------------------------------------------------------------------
    ! Count bad obs errors associated with observations, and skip the ob
    !---------------------------------------------------------------------------
    if (oerr(n) <= 0) then
      !STEVE: this occurred for a few synthetic obs, perhaps due to Dave's code generating obs errors
      cnt_oerlt0 = cnt_oerlt0 + 1
      CYCLE
    endif

    !---------------------------------------------------------------------------
    ! Convert the physical coordinate to model grid coordinate (note: real, not integer)
    !---------------------------------------------------------------------------
    CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
   
    !---------------------------------------------------------------------------
    ! Filter in the tripolar region until localization is examined in the arctic !(ISSUE)
    !---------------------------------------------------------------------------
    if (DO_REMOVE_65N .and. rlat(n) > 65) then
      if (verbose) WRITE(6,'(A)') "Latitude above 65.0N, in tripolar region. Removing observation..."
      cnt_triout = cnt_triout + 1
      CYCLE
    endif

    !---------------------------------------------------------------------------
    ! Filter out observations that are out of range for the grid
    !---------------------------------------------------------------------------
    if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) then
      if (verbose) WRITE(6,'(A)') '* X-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', olon=', rlon(n)
      cnt_xout = cnt_xout + 1
      CYCLE
    endif
    if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) then
      if (verbose) WRITE(6,'(A)') '* Y-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', olat=',rlat(n)
      cnt_yout = cnt_yout + 1
      CYCLE
    endif
    !STEVE: check against kmt, not nlev (OCEAN)
    if (CEILING(rk) > nlev) then
      CALL itpl_2d(kmt0,ri,rj,dk)
      WRITE(6,'(A)') '* Z-coordinate out of range'
      WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
           & '*   rk=',rk,', olev=',rlev(n),&
           & ', (lon,lat)=(',rlon(n),',',rlat(n),'), kmt0=',dk
      cnt_zout = cnt_zout + 1
      CYCLE
    endif
    if (CEILING(rk) < 2 .AND. rk < 1.00001d0) then   !(OCEAN)
      rk = 1.00001d0                                 !(OCEAN)
    endif                                            !(OCEAN)

    !---------------------------------------------------------------------------
    ! Check the observation against boundaries
    !---------------------------------------------------------------------------
    !STEVE: Check to make sure it's in the ocean, as determined       (OCEAN)
    !       by mom4's topography map.
    ! (note: must do it after coordinate checks, or the coordinate
    !        could be outside of the range of the kmt array)
    boundary_points : if (ri > nlon) then
      !STEVE: I have to check what it does for this case...
      !       but it causes an error in the next line if ri > nlon
      if (verbose) WRITE(6,'(A)') '* coordinate is not on mom4 model grid: ri > nlon'
      cnt_rigtnlon = cnt_rigtnlon + 1
      if (cnt_rigtnlon <= 1) then
        WRITE(6,*) "STEVE: ri > nlon (cnt_rigtnlon)"
        WRITE(6,*) "ri = ", ri
        WRITE(6,*) "nlon = ", nlon
        WRITE(6,*) "rj = ", rj
        WRITE(6,*) "elem(n) = ", elem(n)
        WRITE(6,*) "rlon(n) = ", rlon(n)
        WRITE(6,*) "rlat(n) = ", rlat(n)
        WRITE(6,*) "rlev(n) = ", rlev(n)
        WRITE(6,*) "rk = ", rk
      endif
      CYCLE
    else
      !STEVE: check this, case 1 allows more observations, case 2 is more restrictive
      select case (bdyobs)
      case(1)
        if (kmt(NINT(ri),NINT(rj)) .lt. 1) then
          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is on or too close to land, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
            WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          CYCLE
        elseif (kmt(NINT(ri),NINT(rj)) .lt. rk) then
          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is on or too close to underwater topography, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
            WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          CYCLE
        endif
      case(2)
        if(kmt(CEILING(ri),CEILING(rj)) .lt. 1 .or. &
             kmt(CEILING(ri),FLOOR(rj)) .lt. 1 .or. &
             kmt(FLOOR(ri),CEILING(rj)) .lt. 1 .or. &
             kmt(FLOOR(ri),FLOOR(rj)) .lt. 1) THEN

          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is too close to land, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          CYCLE
        elseif(kmt(CEILING(ri),CEILING(rj)) .lt. rk .or. &
                  kmt(CEILING(ri),FLOOR(rj)) .lt. rk .or. &
                  kmt(FLOOR(ri),CEILING(rj)) .lt. rk .or. &
                  kmt(FLOOR(ri),FLOOR(rj)) .lt. rk) THEN

          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is too close to underwater topography, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          CYCLE
        endif
      end select
    endif boundary_points
        
    !---------------------------------------------------------------------------
    ! (OPTIONAL) Scale the observation errors
    !---------------------------------------------------------------------------
    scale_obs : if (abs(obserr_scaling - 1.0d0) .gt. TINY(1.0d0) ) then
      oerr(n) = oerr(n) * obserr_scaling
    endif scale_obs

    !---------------------------------------------------------------------------
    ! (OPTIONAL) Thin the observations
    !---------------------------------------------------------------------------
    ! Select random subset of observations to speed up processing
    thin_obs : if (abs(obs_randselect - 1) > TINY(1.0d0)) then
      CALL com_rand(1,rand)
      if (rand(1) > obs_randselect) CYCLE 
    endif thin_obs

    !---------------------------------------------------------------------------
    ! observation operator (computes H(x)) for specified member
    !---------------------------------------------------------------------------
    CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))

    !GUILLAUME:ADT
    if (DO_ADT .AND. elem(n).eq.id_eta_obs) then
      offset=offset+odat(n)-ohx(n)
      cnt_adt=cnt_adt+1
    end if
    
    if (dodebug .and. DO_ALTIMETRY .and. elem(n) .eq. id_eta_obs) then
      WRITE(6,*) "post-Trans_XtoY:: id_eta_obs, ohx(n) = ", ohx(n)
    endif
    oqc(n) = 1
  enddo !1:nobs

  !GUILLAUME:ADT
  if (DO_ADT) then
    offset=offset/cnt_adt
    WRITE(6,*) 'ADT OFFSET = ',offset
    !GUILLAUME:ADT: Remove offset from adt innovations
    where (elem == id_eta_obs)
       odat=odat-offset
       !ohx=ohx-offset
    end where
  endif

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_adt.f90:: observations removed for:"
  WRITE(6,*) "cnt_oerlt0 = ", cnt_oerlt0
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_zout = ", cnt_zout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_altlatrng = ", cnt_altlatrng

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated innovations to file
  !-----------------------------------------------------------------------------
  CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)

  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,v3d,v2d )

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
  PRINT *, "In grd2cor.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-gues')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      guesfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-aoerin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      aoerinfile = arg2
    case('-aoerout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      aoeroutfile = arg2
    case('-alt')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_ALTIMETRY
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM obsop_adt
