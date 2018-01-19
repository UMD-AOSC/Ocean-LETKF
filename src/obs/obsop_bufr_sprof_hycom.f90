PROGRAM obsop_bufr_sprof
!===============================================================================
! PROGRAM: obsop_bufr_sprof
! 
! USES:
!  use common
!  use params_model
!  use vars_model
!  use common_mom4
!  use params_obs
!  use vars_obs
!  use common_obs_mom4
!
!
! DESCRIPTION: 
!   This program acts as the observation operator for temperature profiles. 
!   It inputs observations of temperature profiles (yo)
!   and a single forecast (xf) and computes the H(xf) used to form the 
!   innovations (yo-H(xf)) associated with that member.
!
! USAGE:
!   A separate instance is run independently for each member and timeslot
!   All observations are typically preprocessed to the letkf obs format, then
!   here they are read in and converted to the letkf obs2 format w/ H(x) data 
!   for each member added in a new column.
! 
! !REVISION HISTORY:
!   06/06/2017 Jili Dong modified for HYCOM 
!   03/17/2016 Steve Penny separated obsop.f90 into separate files
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
!   04/03/2013 Takemasa Miyoshi created for SPEEDY atmospheric model.
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny, Takemasa Miyoshi $
!===============================================================================

  USE common
  USE params_model
  USE vars_model
!  USE common_mom4
  USE common_oceanmodel
  USE params_obs,                ONLY: nobs, id_t_obs, id_s_obs, id_u_obs, id_v_obs, id_eta_obs
  USE vars_obs
!  USE common_obs_mom4
  USE common_obs_oceanmodel
  USE read_bufr,                 ONLY: bufr_data, read_bufr_dumpjb, maxobs


  IMPLICIT NONE

  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  
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
  INTEGER :: i,j,k,n



  REAL(r_size) :: obserr_scaling=1.0d0 !STEVE: use this to scale the input observations
  REAL(r_size) :: obs_randselect=1.0d0 !STEVE: use this to select random input observations
  REAL(r_size), DIMENSION(1) :: rand

  LOGICAL :: remap_obs_coords = .true.

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
  LOGICAL :: print1st = .true.
  LOGICAL :: dodebug1 = .false.
  LOGICAL :: DO_REMOVE_65N = .true.


  INTEGER :: cnt_obs_u=0, cnt_obs_v=0, cnt_obs_t=0, cnt_obs_s=0
  INTEGER :: cnt_obs_ssh=0, cnt_obs_sst=0, cnt_obs_sss=0, cnt_obs_eta=0
  INTEGER :: cnt_obs_x=0, cnt_obs_y=0, cnt_obs_z=0                       
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0, cnt_oerlt0=0, cnt_altlatrng=0

  ! JILI 
  INTEGER ::rim1,rjm1 

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  INTEGER :: lunit=20  !LUBFR for bufr library (i.e. indicates that file 'fort.20' stores the data)
  TYPE(bufr_data), ALLOCATABLE, DIMENSION(:) :: obs_data

  !-----------------------------------------------------------------------------
  ! Initialize the common_mom4 module, and process command line options
  !-----------------------------------------------------------------------------
! JILI For HYCOM use common ocean model
!  CALL set_common_mom4
  CALL set_common_oceanmodel

  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  !STEVE: for lack of a better approach, for now I'm just guessing at the
  !possible number of observations in one given day:
  nobs = maxobs

  print *, "ALLOCATING obs_data with nobs = ", nobs
  ALLOCATE(obs_data(nobs))

  ! Read bufr observations:
  CALL read_bufr_dumpjb(lunit,id_s_obs,obs_data,nobs)
  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  ALLOCATE( obhr(nobs) )

  print *, "obsop_bufr_sprof.f90:: starting nobs = ", nobs
  do i=1,nobs
    elem(i) = obs_data(i)%typ
    rlon(i) = obs_data(i)%x_grd(1)
    rlat(i) = obs_data(i)%x_grd(2)
    rlev(i) = obs_data(i)%x_grd(3)
    odat(i) = obs_data(i)%value
    oerr(i) = obs_data(i)%oerr
    ohx(i)  = 0
    oqc(i)  = 0
    obhr(i) = obs_data(i)%hour
  enddo

  if (print1st .and. nobs > 0) then 
    print *, "elem(1),rlon(1),rlat(1),obhr(1) = ", elem(1),rlon(1),rlat(1),obhr(1)
    print *, "odat(1:40) = ", odat(1:40)
    print *, "oerr(1:40) = ", oerr(1:40)
  endif

  if (dodebug1) print *, "deallocate obs_data..."
  DEALLOCATE(obs_data)

  !-----------------------------------------------------------------------------
  ! Update the coordinate to match the model grid
  ! (Stolen from common_obs_mom4.f90::read_obs
  !-----------------------------------------------------------------------------
  if (remap_obs_coords) then
    if (dodebug1) print *, "Update lon coordinate to match model grid..."
    !CALL center_obs_coords(rlon,oerr,nobs)
  endif

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  WRITE(6,*) "Reading background/forecast data on model grid..."
  CALL read_diag(guesfile,v3d,v2d)
  WRITE(6,*) '****************'

  !-----------------------------------------------------------------------------
  ! Cycle through all observations
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Cycling through observations..."
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
!JILI For HYCOM just get i and j, vertical interpolation will be done later
!wihout rk 
!    CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
    CALL phys2ij(elem(n),rlon(n),rlat(n),rlev(n),ri,rj) !(OCEAN)

   
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
    ! JILI for HYCOM, use observations further from the boundary to be safe             
!    if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) then
    if (CEILING(ri) < 2 .OR. nlon-1 < CEILING(ri)) then
      if (verbose) WRITE(6,'(A)') '* X-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', olon=', rlon(n)
      cnt_xout = cnt_xout + 1
      CYCLE
    endif
!    if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) then
    if (CEILING(rj) < 2 .OR. nlat-1 < CEILING(rj)) then
      if (verbose) WRITE(6,'(A)') '* Y-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', olat=',rlat(n)
      cnt_yout = cnt_yout + 1
      CYCLE
    endif
    !STEVE: check against kmt, not nlev (OCEAN)
    !JILI: For HYCOM, use four neighbouring points bathymetry
    !to check if obs. is above the shallowest bathymetry
!    if (CEILING(rk) > nlev) then
     rim1=FLOOR(ri) 
     rjm1=FLOOR(rj)
!      CALL itpl_2d(kmt0,ri,rj,dk)
    if (minval(phi0(rim1:rim1+1,rjm1:rjm1+1)) < rlev(n)) THEN                        
      WRITE(6,'(A)') '* Z-coordinate out of range'
!      WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
!           & '*   rk=',rk,', olev=',rlev(n),&
!           & ', (lon,lat)=(',rlon(n),',',rlat(n),'), kmt0=',dk

      WRITE(6,'(A,A,F10.2,A,F6.2,A,F6.2,A)') &
           & '* ',', olev=',rlev(n),&
           & ', (lon,lat)=(',rlon(n),',',rlat(n),')'


      cnt_zout = cnt_zout + 1
      CYCLE
    endif
! JILI remove the rk part
!    if (CEILING(rk) < 2 .AND. rk < 1.00001d0) then   !(OCEAN)
!      rk = 1.00001d0                                 !(OCEAN)
!    endif                                            !(OCEAN)

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
      if (verbose) WRITE(6,'(A)') '* coordinate is not on ocean model grid: ri > nlon'
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
!        WRITE(6,*) "rk = ", rk
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
        ! JILI THis is already in the previous code
        !elseif(kmt(CEILING(ri),CEILING(rj)) .lt. rk .or. &
        !          kmt(CEILING(ri),FLOOR(rj)) .lt. rk .or. &
        !          kmt(FLOOR(ri),CEILING(rj)) .lt. rk .or. &
        !          kmt(FLOOR(ri),FLOOR(rj)) .lt. rk) THEN

        !  if (debug_obsfilter) then
        !    WRITE(6,'(A)') '* coordinate is too close to underwater topography, according to kmt'
        !    WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
        !    WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
        !  endif
        !  cnt_nearland = cnt_nearland + 1
        !  CYCLE
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
    ! JILI change rk to rlev for HYCOM
    !CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
    CALL Trans_XtoY(elem(n),ri,rj,rlev(n),v3d,v2d,ohx(n))

    oqc(n) = 1
  enddo !1:nobs

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_bufr_sprof.f90:: observations removed for:"
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
!JILI  skip all obs with bad qc flag by setting qcflag_in as true in write_obs2
  CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,.true.)

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
  PRINT *, "In obsop_bufr_sprof.f90::"
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
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case('-u')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) lunit
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug1
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM obsop_bufr_sprof
