PROGRAM obsop_sss
  USE common
  USE params_model
  USE vars_model
  USE common_oceanmodel
  USE params_obs,                ONLY: nobs, id_sss_obs
  USE params_obs,                ONLY: DO_REMOVE_65N
  USE vars_obs
  USE common_obs_oceanmodel
  USE read_smap,     ONLY: sss_data, read_jpl_smap_l2_sss_h5, &
                                     read_jpl_smap_l3_sss_nc
#ifdef DYNAMIC
  USE input_nml_oceanmodel,      ONLY: read_input_namelist
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Command line inputs:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  REAL(r_size) :: obserr_scaling=1.0d0 !STEVE: use this to scale the input observations
  REAL(r_size) :: obs_randselect=1.0d0 !STEVE: use this to select random input observations
  INTEGER      :: obs_level = -1       !CDA: read level-2 smap if obs_level=2,
                                       !CDA:      level-3 smap if obs_level=3

  !-----------------------------------------------------------------------------
  ! Obs data arrays
  !-----------------------------------------------------------------------------
  LOGICAL :: BINARY_INPUT = .FALSE.
  CHARACTER(10) :: Syyyymmddhh = "YYYYMMDDHH"
  INTEGER :: delta_seconds = -999
  TYPE(sss_data), ALLOCATABLE :: obs_data(:)

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
  ! Model background data arrays
  !-----------------------------------------------------------------------------
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)

  !-----------------------------------------------------------------------------
  ! Miscellaneous
  !-----------------------------------------------------------------------------
  REAL(r_size) :: dk,tg,qg
  REAL(r_size) :: ri,rj,rk
  INTEGER :: n, i,j,k
  REAL(r_size), DIMENSION(1) :: rand
  
  LOGICAL :: remap_obs_coords = .true.

  INTEGER :: bdyobs=2                 !STEVE: use of boundary obs.
                                      !       1 := less restrictive, remove obs inside boundary
                                      !       2 := remove all observations touching a boundary

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  LOGICAL :: debug_obsfilter = .false.
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
  LOGICAL :: dodebug1 = .false.
  LOGICAL :: print1st = .true.

  INTEGER :: cnt_obs_u=0, cnt_obs_v=0, cnt_obs_t=0, cnt_obs_s=0
  INTEGER :: cnt_obs_ssh=0, cnt_obs_sst=0, cnt_obs_sss=0, cnt_obs_eta=0
  INTEGER :: cnt_obs_x=0, cnt_obs_y=0, cnt_obs_z=0
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0, cnt_oerlt0=0, cnt_altlatrng=0

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  INTEGER :: typ = id_sss_obs
  LOGICAL :: DO_SUPEROBS = .false.
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: superobs, delta, M2 ! for online computation of the mean and variance
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: supercnt
  INTEGER :: idx
  INTEGER :: cnt_obs_thinning = 0
  REAL(r_size) :: min_oerr = 1.0 !�C !�C

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
#ifdef DYNAMIC
  CALL read_input_namelist
#endif
  CALL set_common_oceanmodel
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  call create_empty_obsfile(trim(obsoutfile))

  if (BINARY_INPUT) then ![processed binary input]
     print *, "obsop_sss.f90:: reading binary file=",trim(obsinfile)
     CALL get_nobs(trim(obsinfile),8,nobs,errIfNoObs=.false.)

     if (nobs<=0) then
        print*, "obsop_sss.f90: nobs=",nobs,"<=0. Exit now..."
        STOP (0)
     endif

     ALLOCATE( elem(nobs) )
     ALLOCATE( rlon(nobs) )
     ALLOCATE( rlat(nobs) )
     ALLOCATE( rlev(nobs) )
     ALLOCATE( odat(nobs) )
     ALLOCATE( oerr(nobs) )
     ALLOCATE( ohx(nobs) )
     ALLOCATE( oqc(nobs) )
     ALLOCATE( obhr(nobs) )

     CALL read_obs3(trim(obsinfile), nobs, elem, rlon, rlat, rlev, odat, oerr, obhr)
     oqc(:) = 0    ! now set all obs as bad. 
     ohx(:) = 0.d0
  else ! [nc input]
      SELECT CASE(obs_level)
        CASE(2)
          if (delta_seconds>0) then
             CALL read_jpl_smap_l2_sss_h5(trim(obsinfile), obs_data, nobs, &
                                           Syyyymmddhh, delta_seconds)
          else
             CALL read_jpl_smap_l2_sss_h5(trim(obsinfile), obs_data, nobs)
          end if
        CASE(3)
          CALL read_jpl_smap_l3_sss_nc(trim(obsinfile), obs_data, nobs)
        CASE DEFAULT
          WRITE(6,*) "[err] obsop_sss.f90::Unsupported obs_level=", obs_level, &
                     ". obs_level supported should either be 2 or 3"
          STOP (10)
      END SELECT

      ALLOCATE( elem(nobs) )
      ALLOCATE( rlon(nobs) )
      ALLOCATE( rlat(nobs) )
      ALLOCATE( rlev(nobs) )
      ALLOCATE( odat(nobs) )
      ALLOCATE( oerr(nobs) )
      ALLOCATE( ohx(nobs) )
      ALLOCATE( oqc(nobs) )
      ALLOCATE( obhr(nobs) )

      print *, "obsop_sss.f90:: starting nobs = ", nobs
      do i=1,nobs
        elem(i) = obs_data(i)%typ
        rlon(i) = obs_data(i)%x_grd(1)
        rlat(i) = obs_data(i)%x_grd(2)
        rlev(i) = 0.d0
        odat(i) = obs_data(i)%value
        oerr(i) = obs_data(i)%oerr
        ohx(i)  = 0
        oqc(i)  = 0
        obhr(i) = obs_data(i)%hour
      enddo
      DEALLOCATE(obs_data)
  end if

  if (print1st) then  
    print *, "elem(1),rlon(1),rlat(1),obhr(1) = ", elem(1),rlon(1),rlat(1),obhr(1)
    print *, "odat(1:40) = ", odat(1:40)
    print *, "oerr(1:40) = ", oerr(1:40)
  endif

  !-----------------------------------------------------------------------------
  ! Update the coordinate to match the model grid
  !-----------------------------------------------------------------------------
  if (remap_obs_coords) then
    CALL center_obs_coords(rlon,oerr,nobs)
  endif

  !-----------------------------------------------------------------------------
  ! Bin the obs and estimate the obs error based on bin standard deviations
  !-----------------------------------------------------------------------------
  if (DO_SUPEROBS) then
    ALLOCATE(superobs(nlon,nlat),delta(nlon,nlat),M2(nlon,nlat),supercnt(nlon,nlat))
    superobs=0.0; delta = 0.0; M2 = 0.0; supercnt = 0

    print *, "Computing superobs..."

    !STEVE: (ISSUE) this may have problems in the arctic for the tripolar grid,
    !               or for any irregular grid in general.
    do n=1,nobs ! for each ob,
!     if (dodebug1) print *, "n = ", n

      if (DO_REMOVE_65N .and. rlat(n) > 65) CYCLE
      !! Scan the longitudes
      !do i=1,nlon-1
      !  if (lon(i+1) > rlon(n)) exit
      !enddo
      !! Scan the latitudes
      !do j=1,nlat-1
      !  if (lat(j+1) > rlat(n)) exit
      !enddo
      !CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
      CALL phys2ij_nearest(rlon(n),rlat(n),ri,rj) ! OCEAN

      if (CEILING(ri) < 1 .OR. nlon < CEILING(ri)) CYCLE
      if (CEILING(rj) < 1 .OR. nlat < CEILING(rj)) CYCLE

      i = NINT(ri); j = NINT(rj)
      supercnt(i,j) = supercnt(i,j) + 1
      delta(i,j)    = odat(n) - superobs(i,j)
      superobs(i,j) = superobs(i,j) + delta(i,j)/supercnt(i,j)
      M2(i,j)       = M2(i,j) + delta(i,j)*(odat(n) - superobs(i,j))
!     if (dodebug1) print *, "supercnt(",i,",",j,") = ", supercnt(i,j)
    enddo
    !"superobs" contains the mean
    !M2 contains the variance:
    WHERE(supercnt > 1) M2 = M2 / (supercnt - 1)

    idx=0
    do j=1,nlat-1
      do i=1,nlon-1
        if (supercnt(i,j) > 1 .and. wet(i,j)>0.5) then ! only retain ocn pts defined by MOM
          idx = idx+1
          if (dodebug1) print *, "idx = ", idx
          odat(idx) = superobs(i,j)
          oerr(idx) = obserr_scaling*(min_oerr + SQRT(M2(i,j)))
          rlon(idx) = lon(i) !(lon(i+1)+lon(i))/2.0d0
          rlat(idx) = lat(j) !(lat(j+1)+lat(j))/2.0d0
          rlev(idx) = 0.d0
          elem(idx) = id_sss_obs
          if (dodebug1) print *, "odat(idx) = ", odat(idx)
          if (dodebug1) print *, "oerr(idx) = ", oerr(idx)
          if (dodebug1) print *, "rlon(idx) = ", rlon(idx)
          if (dodebug1) print *, "rlat(idx) = ", rlat(idx)
          if (dodebug1) print *, "ocnt(idx) = ", supercnt(i,j)
        endif
      enddo
    enddo
    print *, "DO_SUPEROBS:: superobs reducing from ",nobs," to ",idx, " observations."
    print *, "with min obs error = ", MINVAL(oerr(1:idx))
    print *, "with max obs error = ", MAXVAL(oerr(1:idx))
    nobs = idx
  endif

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  !STEVE: (ISSUE) this can be improved by only reading in the file for needed for this data.
  CALL read_diag(guesfile,v3d,v2d)
  WRITE(6,*) '****************'

  !-----------------------------------------------------------------------------
  ! Cycle through all observations
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Cycling through observations..."

  ohx=0.0d0
  oqc=0
  idx=0
  DO n=1,nobs
    ! Thin observations
    ! Select random subset of observations to speed up processing
    !STEVE: I know this would be more efficient up above when obs are first introduced.
    thin : if (obs_randselect < 1.0d0) then
      CALL com_rand(1,rand)
      if (rand(1) > obs_randselect) then
        cnt_obs_thinning = cnt_obs_thinning + 1
        CYCLE
      endif
    endif thin

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
    ! observation operator (computes H(x)) for specified member
    !---------------------------------------------------------------------------
    CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
    idx=idx+1
    oqc(n) = 1
  enddo !1:nobs

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_sss.f90::"
  WRITE(6,*) "observations at start = ", nobs
  WRITE(6,*) "== observations removed for: =="
  WRITE(6,*) "cnt_obs_thinning = ", cnt_obs_thinning
  WRITE(6,*) "cnt_oerlt0 = ", cnt_oerlt0
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_zout = ", cnt_zout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_altlatrng = ", cnt_altlatrng
  WRITE(6,*) "==============================="
  WRITE(6,*) "observations kept = ", idx
  WRITE(6,*) "==============================="

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated innovations to file
  !-----------------------------------------------------------------------------
  CALL write_obs2(obsoutfile,nobs,elem(1:nobs), &
                                  rlon(1:nobs), &
                                  rlat(1:nobs), &
                                  rlev(1:nobs), &
                                  odat(1:nobs), &
                                  oerr(1:nobs), &
                                   ohx(1:nobs), &
                                   oqc(1:nobs), &
                                  obhr(1:nobs), &
                                  lappend_in = .true.)

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
  PRINT *, "In obsop_sss.f90::"
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
    case('-superob')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_SUPEROBS
    case('-binary')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) BINARY_INPUT
    case('-thin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) obs_randselect
    case('-scale')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) obserr_scaling
    case('-minerr')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) min_oerr
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug1
    case('-olevel')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) obs_level
    case('-qcyyyymmddhh')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      !read (arg2,*) min_quality_level
      Syyyymmddhh = arg2
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

END PROGRAM obsop_sss
