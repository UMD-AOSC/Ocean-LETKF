PROGRAM obsope
!=======================================================================
!
! [PURPOSE:] Main program of observation operator
!
! A separate instance is run independently for each member and timeslot
!
! [HISTORY:]
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
!   04/03/2013 Takemasa Miyoshi created
!
!=======================================================================
  USE common
  USE common_mom4
  USE common_obs_mom4

  IMPLICIT NONE
  CHARACTER(slen) :: obsinfile='obsin.dat'    !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  CHARACTER(slen) :: aoerinfile='oer_inp.grd' !IN (default)
  CHARACTER(slen) :: aoeroutfile='oer_inp.grd'!OUT(default)
  REAL(r_size),ALLOCATABLE :: elem(:)
  REAL(r_size),ALLOCATABLE :: rlon(:)
  REAL(r_size),ALLOCATABLE :: rlat(:)
  REAL(r_size),ALLOCATABLE :: rlev(:)
  REAL(r_size),ALLOCATABLE :: odat(:)
  REAL(r_size),ALLOCATABLE :: oerr(:)
  REAL(r_size),ALLOCATABLE :: ohx(:)
  INTEGER,ALLOCATABLE :: oqc(:)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: o3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: o2d(:,:,:)
  REAL(r_size) :: dk,tg,qg
  INTEGER :: nobs
  REAL(r_size) :: ri,rj,rk
  INTEGER :: n
  !STEVE:
  !STEVE: for (DRIFTERS)
  REAL(r_size),ALLOCATABLE,SAVE :: obsid(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obstime(:)
  !STEVE:
  REAL(r_size),DIMENSION(nid_obs),PARAMETER :: & !STEVE: use this to scale the input observations
              obserr_scaling=(/ 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0 /)
!             obserr_scaling=(/ 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.0d0, 0.5d0, 1.0d0 /)
                              ! U       V       Temp    Salt    SSH    SST    SSS 
  REAL(r_size),DIMENSION(nid_obs),PARAMETER :: & !STEVE: use this to select random input observations
              obs_randselect=(/ 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.00d0 /)
!             obs_randselect=(/ 1.00d0, 1.00d0, 1.00d0, 1.00d0, 1.00d0, 0.10d0, 1.00d0 /)
                              ! U       V       Temp    Salt    SSH     SST     SSS 
  REAL(r_size), DIMENSION(1) :: rand

  INTEGER :: bdyobs=2 !STEVE: use of boundary obs. 1 = less restrictive, remove obs inside boundary; 2 = remove all observations touching a boundary
  LOGICAL :: debug_obsfilter = .false.
  LOGICAL :: debug_hdxf_0 = .false.   !This error occured because there was not a model representation of the observed value (i.e. SST obs with no SST model field)
                                      ! Solution was to populate a SST model field (v2d) with surface temp data from the model (v3d(:,:,1))
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
  !STEVE: for using satellite altimeter data of SSH
! LOGICAL :: DO_ALTIMETRY !STEVE: now in common_mom4.f90
  REAL(r_size), DIMENSION(nlon,nlat) :: Lxbar
  REAL(r_size) :: HLxbar
  REAL(r_sngl), DIMENSION(nlon,nlat) :: buf4
  INTEGER :: iolen, fid=13
  CHARACTER(slen) :: Lxfile='Lxbar.grd'
  LOGICAL :: ex

  CALL set_common_mom4
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  CALL get_nobs(obsinfile,6,nobs)
  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  CALL read_obs(trim(obsinfile),nobs,elem,rlon,rlat,rlev,odat,oerr)

  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  CALL read_grd(trim(guesfile),v3d,v2d)

  !!STEVE: for adaptive observation error:
  INQUIRE(FILE=aoerinfile, EXIST=oerfile_exists)
  if (oerfile_exists) then
    ALLOCATE(o3d(nlon,nlat,nlev,nv3d),o2d(nlon,nlat,nv2d))
    CALL read_bingrd(trim(aoerinfile),o3d,o2d)
  endif

! altimetry_file : if (DO_ALTIMETRY) then
!   INQUIRE(FILE=trim(Lxfile),EXIST=ex)
!   IF(ex) THEN
!     ! Read Lxbar from file (precomputed)
!     INQUIRE(IOLENGTH=iolen) iolen
!     WRITE(6,*) "iolen = ", iolen
!     WRITE(6,*) "Reading file: ", trim(Lxfile)
!     open(fid,FILE=trim(Lxfile),FORM='unformatted') !,ACCESS='direct',RECL=nlon*nlat*iolen)
!!    read(fid,REC=1) buf4
!     read(fid) buf4
!     WRITE(6,*) "MAXVAL(buf4) = ", MAXVAL(buf4)
!     WRITE(6,*) "MINVAL(buf4) = ", MINVAL(buf4)
!     Lxbar = REAL(buf4,r_size)
!     close(fid)
!     WRITE(6,*) "MAXVAL(Lxbar) = ", MAXVAL(Lxbar)
!     WRITE(6,*) "MINVAL(Lxbar) = ", MINVAL(Lxbar)
!     if (MAXVAL(Lxbar) < epsilon(1.0d0) .and. MINVAL(Lxbar) < epsilon(1.0d0) ) then
!       print *, "Lxbar from file is zero (0.0d0) in all values. Exiting..."
!       stop(1)
!     endif
!   ELSE
!     WRITE(6,*) "File does not exist: ", trim(Lxfile)
!     WRITE(6,*) "Exiting..."
!     stop(2)
!   ENDIF
! else
!   WRITE(6,*) "Not using ALTIMETRY, thus not reading: ", trim(Lxfile)
! endif altimetry_file

  ohx=0.0d0
  oqc=0
  DO n=1,nobs
    if (oerr(n) <= 0) then
      !STEVE: this occurred for a few synthetic obs, perhaps due to Dave's code generating obs errors
      cnt_oerlt0 = cnt_oerlt0 + 1
      CYCLE
    endif

    CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
   
    !STEVE: 5/24/2013, removing obs in the tripolar region due to
        !'saturation vapor pressure table overflow' errors that may be due to
        !problems in the ice model, or possible due to how the tripolar grid is
        !handled in letkf.
        !STEVE: 3/17/2014: It was the ice model, and probably due to perturbations in the ice model.
        !                  For now, the ice model has been removed.
        if (rlat(n) > 65) then
            if (verbose) WRITE(6,'(A)') "Latitude above 65.0N, in tripolar region. Removing observation..."
            cnt_triout = cnt_triout + 1
            CYCLE
            !STEVE: 9/5/2013, trying to keep obs, but increasing error to keep
            !the increment as smooth as possible.
            if (verbose) WRITE(6,'(A)') "Latitude above 60, in tripolar region. Increasing obs error..."
            oerr(n) = oerr(n)*3.0d0 !*gross_error
!       elseif (rlat(n) < -65) then
!           !STEVE: 9/6/2013, trying to keep obs, but increasing error to keep
!           !the increment as smooth as possible.
!           if (verbose) WRITE(6,'(A)') "Latitude below 65.0S. Increasing obs error..."
!           oerr(n) = oerr(n)*3.0d0 !*gross_error
        endif

    IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) THEN
      if (verbose) WRITE(6,'(A)') '* X-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', olon=', rlon(n)
      cnt_xout = cnt_xout + 1
      CYCLE
    END IF
    IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) THEN
      if (verbose) WRITE(6,'(A)') '* Y-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', olat=',rlat(n)
      cnt_yout = cnt_yout + 1
      CYCLE
    END IF
    !STEVE: check against kmt, not nlev (OCEAN)
    IF(CEILING(rk) > nlev) THEN
      CALL itpl_2d(kmt0,ri,rj,dk)
      WRITE(6,'(A)') '* Z-coordinate out of range'
      WRITE(6,'(A,F6.2,A,F10.2,A,F6.2,A,F6.2,A,F10.2)') &
           & '*   rk=',rk,', olev=',rlev(n),&
           & ', (lon,lat)=(',rlon(n),',',rlat(n),'), kmt0=',dk
      cnt_zout = cnt_zout + 1
      CYCLE
    END IF
    IF(CEILING(rk) < 2 .AND. rk < 1.00001d0) THEN    !(OCEAN)
      rk = 1.00001d0                                 !(OCEAN)
    END IF                                           !(OCEAN)

    !STEVE: Check to make sure it's in the ocean, as determined       (OCEAN)
    !       by mom4's topography map.
    ! (note: must do it after coordinate checks, or the coordinate
    !        could be outside of the range of the kmt array)
    boundary_points : IF (ri > nlon) THEN
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
        !STEVE: check this, case 1 allows more, case 2 is more restrictive
        select case (bdyobs)
        case(1)
          IF(kmt(NINT(ri),NINT(rj)) .lt. 1) THEN
            if (debug_obsfilter) then
              WRITE(6,'(A)') '* coordinate is on or too close to land, according to kmt'
              WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
              WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
            endif
            cnt_nearland = cnt_nearland + 1
            CYCLE
          ELSEIF (kmt(NINT(ri),NINT(rj)) .lt. rk) THEN
            if (debug_obsfilter) then
              WRITE(6,'(A)') '* coordinate is on or too close to underwater topography, according to kmt'
              WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
              WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
            endif
            cnt_nearland = cnt_nearland + 1
            CYCLE
          ENDIF
        case(2)
          IF(kmt(CEILING(ri),CEILING(rj)) .lt. 1 .or. &
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
           ELSEIF(kmt(CEILING(ri),CEILING(rj)) .lt. rk .or. &
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
           ENDIF
         end select
       ENDIF boundary_points
        
       scale_obs : if (abs(maxval(obserr_scaling) - 1.0d0) .gt. (TINY(1.0d0)) .or. &
                       abs(minval(obserr_scaling) - 1.0d0) .gt. (TINY(1.0d0)) .or. &
                       abs(maxval(obs_randselect) - 1.0d0) .gt. (TINY(1.0d0)) .or. &
                       abs(minval(obs_randselect) - 1.0d0) .gt. (TINY(1.0d0)) ) then

          !STEVE: Process and correct observational error: 
          !STEVE: Scale observation error if requested.
          !STEVE: Select random subset of observations if requested.
          !STEVE: if the observation value is far off from the mean background, then drop it.
          !STEVE: NOTE: this might be more efficient if added above in more parallel
          !section
          SELECT CASE(NINT(elem(n)))
            CASE (id_u_obs) !u
              if (abs(obserr_scaling(iv3d_u) - 1.0) > TINY(1.0d0)) &
                oerr(n) = oerr(n) * obserr_scaling(iv3d_u)
            CASE (id_v_obs) !v
              if (abs(obserr_scaling(iv3d_v) - 1.0) > TINY(1.0d0)) &
                oerr(n) = oerr(n) * obserr_scaling(iv3d_v)
            CASE (id_t_obs) !temp
              if (abs(obserr_scaling(iv3d_t) - 1.0) > TINY(1.0d0)) then
!               WRITE(6,*) "OLD temp oerr = ", oerr(n)
                oerr(n) = oerr(n) * obserr_scaling(iv3d_t)
!               WRITE(6,*) "NEW temp oerr = ", oerr(n)
              endif
            CASE (id_s_obs) !salt
              if (abs(obserr_scaling(iv3d_s) - 1.0) > TINY(1.0d0)) then
!               WRITE(6,*) "OLD salt oerr = ", oerr(n)
                oerr(n) = oerr(n) * obserr_scaling(iv3d_s)
!               WRITE(6,*) "NEW salt oerr = ", oerr(n)
              endif
            CASE (id_ssh_obs) !ssh
              if (abs(obserr_scaling(iv2d_ssh) - 1.0) > TINY(1.0d0)) &
                oerr(n) = oerr(n) * obserr_scaling(iv2d_ssh)
            CASE (id_sst_obs) !surface temp
              ! Select random subset of observations to speed up processing
              if (abs(obs_randselect(iv2d_sst) - 1) > TINY(1.0d0)) then
                CALL com_rand(1,rand)
                if (rand(1) > obs_randselect(iv2d_sst)) CYCLE
              endif
              if (abs(obserr_scaling(iv2d_sst) - 1.0) > TINY(1.0d0)) then
!               WRITE(6,*) "OLD sst  oerr = ", oerr(n)
                oerr(n) = oerr(n) * obserr_scaling(iv2d_sst)
!               WRITE(6,*) "NEW sst  oerr = ", oerr(n)
              endif
            CASE (id_sss_obs) !surface salt
              if (abs(obserr_scaling(iv2d_sss) - 1.0) > TINY(1.0d0)) &
                oerr(n) = oerr(n) * obserr_scaling(iv2d_sss)
            CASE DEFAULT
              WRITE(6,*) "for n = ", n
                WRITE(6,*) "WARNING: for adaptive observation error, no support for observation type: ", elem(n)
          END SELECT
       endif scale_obs

       !STEVE: apply adaptive observation error (OCEAN)
       !       I'm assuming that any type of Kalman filter
       !       combining the old and the new observation errors
       !       that should be applied will have been done externally,
       !       prior to reading it in here.
       ! adapt obs adapt_obs 
       adapt_obs_error : if (oerfile_exists) then
          SELECT CASE(NINT(elem(n)))
          CASE (id_u_obs) !u
            oberr = o3d(NINT(ri),NINT(rj),NINT(rk),iv3d_u) 
          CASE (id_v_obs) !v
            oberr = o3d(NINT(ri),NINT(rj),NINT(rk),iv3d_v) 
          CASE (id_t_obs) !temp
            oberr = o3d(NINT(ri),NINT(rj),NINT(rk),iv3d_t) 
          CASE (id_s_obs) !salt
            oberr = o3d(NINT(ri),NINT(rj),NINT(rk),iv3d_s) 
          CASE (id_ssh_obs) !ssh
            oberr = o2d(NINT(ri),NINT(rj),iv2d_ssh) 
          CASE (id_sst_obs) !surface temp
            oberr = o2d(NINT(ri),NINT(rj),iv2d_sst) 
          CASE (id_sss_obs) !surface salt
            oberr = o2d(NINT(ri),NINT(rj),iv2d_sss) 

!!(DRIFTERS)
!!STEVE: not sure how to do this yet for drifters...
!!        CASE (id_x_obs) !x
!!          oberr = o4d(NINT(ri),NINT(rj),NINT(rk),iv4d_x) 
!!        CASE (id_y_obs) !y
!!          oberr = o4d(NINT(ri),NINT(rj),NINT(rk),iv4d_y) 
!!        CASE (id_z_obs) !z
!!          oberr = o4d(NINT(ri),NINT(rj),NINT(rk),iv4d_z) 
          CASE DEFAULT
            WRITE(6,*) "for n = ", n
            WRITE(6,*) "WARNING: for adaptive observation error, no support for observation type: ", elem(n)
          END SELECT
          if (oberr > 0) oerr(n) = oberr
        endif adapt_obs_error 

    !STEVE: since I'm not relying on the linear steric operator, I'm going to keep all the altimetry data for now...
!   altimetry : if (DO_ALTIMETRY .and. elem(n) .eq. id_eta_obs ) then
!     !STEVE: outside of this range, the linear steric operator is questionable
!     if (rlat(n) > 40 .or. rlat(n) < -30) then
!       cnt_altlatrng = cnt_altlatrng + 1
!       CYCLE
!     endif
!
!     !STEVE: since it's easier to code, I'm doing:
!     ! (SSH_o - CLM_o) - (SSH_m - CLM_m) + HLxb - HLx
!     ! (SLA_o - SLA_m + HLxb) - HLx
!     ! SLA_o and SLA_m are read in from file (climatology is subtracted from model upon reading in)
!     ! Here I'm just computing the first part:
!     ! yo^ = (SLA_o - SLA_m + HLxb)
!
!     ! obsdep (yo^ - HLx_mean) will then be computed later in letkf_obs.f90
!     ! when we finally get to the H operator (Trans_XY), we're just computing (HLx)
!
!     ! Create adjusted observation  to compare with H(x_SLA + HLx)
!     CALL itpl_2d(Lxbar(:,:),ri,rj,HLxbar)
!     if (dodebug) then
!       WRITE(6,*) "Lxbar(floor(ri),floor(rj)) = ", Lxbar(floor(ri),floor(rj))
!       WRITE(6,*) "Lxbar(ceiling(ri),floor(rj)) = ", Lxbar(ceiling(ri),floor(rj))
!       WRITE(6,*) "Lxbar(ceiling(ri),ceiling(rj)) = ", Lxbar(ceiling(ri),ceiling(rj))
!       WRITE(6,*) "Lxbar(floor(ri),ceiling(rj)) = ", Lxbar(floor(ri),ceiling(rj))
!       WRITE(6,*) "HLxbar = ", HLxbar
!       WRITE(6,*) "odat(n) = ", odat(n)
!       WRITE(6,*) "adjusted ob = ", odat(n) + HLxbar
!     endif
!     odat(n) = odat(n) + HLxbar !i.e.: yo_SLA + H(Lxb)
!                                !STEVE: this essentially puts the observed climatological anomaly on par with the modeled clim. anomaly
!     ! We will output the modified observation (odat) and the H(xb) for this member as: ohx
!     
!   endif altimetry
 
    !
    ! observational operator (computes H(x), or for altimetry, H(x_SLA + Lx))
    !
    CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
    if (dodebug .and. DO_ALTIMETRY .and. elem(n) .eq. id_eta_obs) then
      WRITE(6,*) "post-Trans_XtoY:: id_eta_obs, ohx(n) = ", ohx(n)
    endif
    oqc(n) = 1
  END DO  !1:nobs

  WRITE(6,*) "In letkf_obs.f90:: observations removed for:"
  WRITE(6,*) "cnt_oerlt0 = ", cnt_oerlt0
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_zout = ", cnt_zout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_altlatrng = ", cnt_altlatrng

  CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc)

  if (ALLOCATED(o3d)) DEALLOCATE(o3d)
  if (ALLOCATED(o2d)) DEALLOCATE(o2d)
  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,v3d,v2d )

CONTAINS

SUBROUTINE process_command_line

IMPLICIT NONE

INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
DO i=1,COMMAND_ARGUMENT_COUNT(),2
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
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
ENDDO

END SUBROUTINE process_command_line

END PROGRAM obsope
