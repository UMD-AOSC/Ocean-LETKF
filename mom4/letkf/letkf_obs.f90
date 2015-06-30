MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!            Reads in all observational data
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   04/26/2011 Steve PENNY converted to OCEAN for use with MOM4
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_mom4
  USE common_obs_mom4
  USE common_mpi_mom4
  USE common_letkf
  !(DRIFTERS)
! USE letkf_drifters

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs
  !STEVE: making these namelist accessible:
  INTEGER :: nslots=5 ! number of time slots for 4D-LETKF
  INTEGER :: nbslot=5 !1 !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
  REAL(r_size) :: sigma_obs=720.0d3 !3x Rossby Radius of Deformation at Equ. according to Chelton
  REAL(r_size) :: sigma_obs0=200.0d3 !20x Rossby Radius of Deformation at Pole, according to Chelton
  REAL(r_size) :: sigma_obsv=1000.0d0  !STEVE: doesn't matter if using option "DO_NO_VERT_LOC"
  REAL(r_size) :: sigma_obst=5.0d0
! INTEGER,PARAMETER :: nslots=5 ! number of time slots for 4D-LETKF
! INTEGER,PARAMETER :: nbslot=5 !1 !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
! REAL(r_size),PARAMETER :: sigma_obs=720.0d3 !3x Rossby Radius of Deformation at Equ. according to Chelton
! REAL(r_size),PARAMETER :: sigma_obs0=200.0d3 !20x Rossby Radius of Deformation at Pole, according to Chelton
! REAL(r_size),PARAMETER :: sigma_obsv=1000.0d0  !STEVE: doesn't matter if using option "DO_NO_VERT_LOC"
! REAL(r_size),PARAMETER :: sigma_obst=5.0d0
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),ALLOCATABLE,SAVE :: dlon_zero(:)
  REAL(r_size),SAVE :: dlat_zero
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
  INTEGER, ALLOCATABLE, SAVE :: obs_useidx(:) !STEVE: general version of nobs_use()
  INTEGER,SAVE :: nobsgrd(nlon,nlat)
  !STEVE: for (DRIFTERS)
  REAL(r_size),ALLOCATABLE,SAVE :: obsid(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obstime(:)
  !STEVE:
  LOGICAL :: debug_hdxf_0 = .true.   !This error occured because there was not a model representation of the observed value (i.e. SST obs with no SST model field)
                                     ! Solution was to populate a SST model field (v2d) with surface temp data from the model (v3d(:,:,1))
  INTEGER :: cnt_obs_u, cnt_obs_v, cnt_obs_t, cnt_obs_s, cnt_obs_x, cnt_obs_y, cnt_obs_z, cnt_obs_ssh, cnt_obs_eta, cnt_obs_sst, cnt_obs_sss
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0
  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0
  !STEVE: for adaptive obs:
  LOGICAL :: oerfile_exists
  !STEVE: for using satellite altimeter data of SSH
! LOGICAL :: DO_ALTIMETRY !STEVE: now in common_mom4.f90
! REAL(r_size),PARAMETER :: gross_error=10.0d0 !3.0d0 ! number of standard deviations   (Use for OSSEs)
                                                      ! used to filter out observations
  REAL(r_size) :: gross_error=3.0d0 ! number of standard deviations

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  IMPLICIT NONE
  REAL(r_size) :: dz,tg,qg
  REAL(r_size) :: ri,rj,rk
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmpi(:)
  REAL(r_size),ALLOCATABLE :: tmpj(:)
  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  REAL(r_size),ALLOCATABLE :: tmpid(:)   !(DRIFTERS)
  REAL(r_size),ALLOCATABLE :: tmptime(:) !(DRIFTERS)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: tmp2id(:)   !(DRIFTERS)
  REAL(r_size),ALLOCATABLE :: tmp2time(:) !(DRIFTERS)
  INTEGER :: nobslots(nslots)
  INTEGER :: n,i,j,ierr,islot,nn,l,im
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(12) :: obsfile='obsTTNNN.dat'

  !STEVE: for adaptive observation error:
  REAL(r_size) :: tmpoerr 
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
  LOGICAL :: dodebug = .true.
  !STEVE: for obs qc:
  REAL(r_size) :: hdx2,mstd
  INTEGER :: gross_cnt,gross_2x_cnt
  !STEVE: for DO_ALTIMETRY
  REAL(r_size) :: SSH_CLM_m !SSHclm_ij

  WRITE(6,'(A)') 'Hello from set_letkf_obs'

  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dlat_zero = dist_zero / pi / re * 180.0d0
  ALLOCATE(dlon_zero(nij1))
  DO i=1,nij1
    dlon_zero(i) = dlat_zero / COS(pi*lat1(i)/180.0d0)
  END DO

  IF(myrank == 0) THEN !Assuming all members have the identical obs records
    DO islot=1,nslots
      im = myrank+1
      WRITE(obsfile(4:8),'(I2.2,I3.3)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading an obs2-formatted file ',obsfile
      CALL get_nobs(obsfile,8,nobslots(islot))
    ENDDO
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nobslots,nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  nobs = SUM(nobslots)
  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'

  IF(nobs == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    RETURN
  END IF

!
! INITIALIZE GLOBAL VARIABLES
!
  ALLOCATE( tmpelm(nobs) )
  ALLOCATE( tmplon(nobs) )
  ALLOCATE( tmplat(nobs) )
  ALLOCATE( tmplev(nobs) )
  ALLOCATE( tmpdat(nobs) )
  ALLOCATE( tmperr(nobs) )
  ALLOCATE( tmpk(nobs) )
  ALLOCATE( tmpdep(nobs) )
  ALLOCATE( tmphdxf(nobs,nbv) )
  ALLOCATE( tmpqc0(nobs,nbv) )
  ALLOCATE( tmpqc(nobs) )
  ALLOCATE( tmpid(nobs) )   !(DRIFTERS)
  ALLOCATE( tmptime(nobs) ) !(DRIFTERS)
  tmpqc0 = 0
  tmphdxf = 0.0d0
  tmperr = 0.0d0

!
! LOOP of timeslots
!
  nn=0
  timeslots0: DO islot=1,nslots
    IF(nobslots(islot) == 0) CYCLE
    l=0
    DO
      im = myrank+1 + nprocs * l
      IF(im > nbv) EXIT
      WRITE(obsfile(4:8),'(I2.2,I3.3)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile
      CALL read_obs2(obsfile,nobslots(islot),&
       & tmpelm(nn+1:nn+nobslots(islot)),tmplon(nn+1:nn+nobslots(islot)),&
       & tmplat(nn+1:nn+nobslots(islot)),tmplev(nn+1:nn+nobslots(islot)),&
       & tmpdat(nn+1:nn+nobslots(islot)),tmperr(nn+1:nn+nobslots(islot)),&
       & tmphdxf(nn+1:nn+nobslots(islot),im),tmpqc0(nn+1:nn+nobslots(islot),im))
      l = l+1
    ENDDO
    nn = nn + nobslots(islot)
  ENDDO timeslots0

  WRITE(6,*) "Commencing collecting obs on all procs..."
  !STEVE: broadcast the 1d arrays from root onto all procs
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_BCAST's..."
  CALL MPI_BCAST( tmpelm, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  !STEVE: just to be safe, calling MPI_BARRIER after each MPI_BCAST
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplon, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplat, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmplev, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmpdat, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST( tmperr, nobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !STEVE: compile the tmphdxf array on all procs
  ALLOCATE(wk2d(nobs,nbv))
  wk2d = tmphdxf
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_ALLREDUCE..."
  CALL MPI_ALLREDUCE(wk2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(wk2d)

  !STEVE: compile the tmpqc0 array on all procs
  ALLOCATE(iwk2d(nobs,nbv))
  iwk2d = tmpqc0
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(6,*) "Calling MPI_ALLREDUCE..."
  CALL MPI_ALLREDUCE(iwk2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(iwk2d)
  WRITE(6,*) "Finished collecting obs on all procs."

  WRITE(6,*) "STEVE: DEBUGGING..."
  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS, MYRANK = ',myrank
  WRITE(6,*) "tmphdxf(1,:) = ", tmphdxf(1,:)
  WRITE(6,*) "tmphdxf(2,:) = ", tmphdxf(2,:)
  WRITE(6,*) "tmphdxf(3,:) = ", tmphdxf(3,:)
  WRITE(6,*) "..."
  WRITE(6,*) "tmphdxf(nobs,:) = ", tmphdxf(nobs,:)
  WRITE(6,*)
  n=1
  WRITE(6,*) "For n=1,"
  WRITE(6,*) "MINVAL(tmpqc0(n,:)) = ",MINVAL(tmpqc0(n,:))
  WRITE(6,*) "tmpelm(n) = ", tmpelm(n) 
  WRITE(6,*) "tmplon(n) = ", tmplon(n) 
  WRITE(6,*) "tmplat(n) = ", tmplat(n) 
  WRITE(6,*) "tmplev(n) = ", tmplev(n) 
  WRITE(6,*) "tmpdat(n) = ", tmpdat(n) 
  WRITE(6,*) "tmperr(n) = ", tmperr(n) 
  WRITE(6,*)
  n=nobs
  WRITE(6,*) "For n=nobs=",nobs
  WRITE(6,*) "MINVAL(tmpqc0(n,:)) = ",MINVAL(tmpqc0(n,:))
  WRITE(6,*) "tmpelm(n) = ", tmpelm(n) 
  WRITE(6,*) "tmplon(n) = ", tmplon(n) 
  WRITE(6,*) "tmplat(n) = ", tmplat(n) 
  WRITE(6,*) "tmplev(n) = ", tmplev(n) 
  WRITE(6,*) "tmpdat(n) = ", tmpdat(n) 
  WRITE(6,*) "tmperr(n) = ", tmperr(n) 
  WRITE(6,*) "STEVE: END DEBUGGING."
  WRITE(6,*)

  !STEVE: After processing ensemble members, apply some actions based on
  !forecast mean, to all observations
  cnt_obs_u = 0
  cnt_obs_v = 0
  cnt_obs_t = 0
  cnt_obs_s = 0
  cnt_obs_x = 0
  cnt_obs_y = 0
  cnt_obs_z = 0
  cnt_obs_ssh = 0
  cnt_obs_eta = 0
  cnt_obs_sst = 0
  cnt_obs_sss = 0
  gross_cnt = 0
  gross_2x_cnt = 0
  
  WRITE(6,*) "Processing tmphdxf for n=1 to n=nobs=",nobs 
  WRITE(6,*) "and filtering bad observations..."

if (.true.) then
  !STEVE: this is the original version

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
  DO n=1,nobs
    tmpqc(n) = MINVAL(tmpqc0(n,:))
    IF(tmpqc(n) /= 1) CYCLE
    tmpdep(n) = tmphdxf(n,1) !note: tmpdep is just used as a dummy variable to compute the mean over the next few lines
    DO i=2,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdxf (perturbations from mean)
    END DO
    ! Now, tmpdep is defined appropriately as the obs departure from mean background
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
      tmpqc(n) = 0
      gross_cnt = gross_cnt + 1
    END IF

    !STEVE: as a check, count the number of each type of observation
    if (tmpelm(n) .eq. id_u_obs) cnt_obs_u = cnt_obs_u + 1
    if (tmpelm(n) .eq. id_v_obs) cnt_obs_v = cnt_obs_v + 1
    if (tmpelm(n) .eq. id_t_obs) cnt_obs_t = cnt_obs_t + 1
    if (tmpelm(n) .eq. id_s_obs) cnt_obs_s = cnt_obs_s + 1
    if (tmpelm(n) .eq. id_ssh_obs) cnt_obs_ssh = cnt_obs_ssh + 1
    if (tmpelm(n) .eq. id_eta_obs) cnt_obs_eta = cnt_obs_eta + 1
    if (tmpelm(n) .eq. id_sst_obs) cnt_obs_sst = cnt_obs_sst + 1
    if (tmpelm(n) .eq. id_sss_obs) cnt_obs_sss = cnt_obs_sss + 1
    !(DRIFTERS)
    if (tmpelm(n) .eq. id_x_obs) cnt_obs_x = cnt_obs_x + 1
    if (tmpelm(n) .eq. id_y_obs) cnt_obs_y = cnt_obs_y + 1
    if (tmpelm(n) .eq. id_z_obs) cnt_obs_z = cnt_obs_z + 1

  END DO
!$OMP END PARALLEL DO
  DEALLOCATE(tmpqc0)

else

  !STEVE: this is the augmented version I created to do more
  !       based on ensemble spread

  DO n=1,nobs
    !WRITE(6,*) "n = ", n

    tmpqc(n) = MINVAL(tmpqc0(n,:))
    IF(tmpqc(n) /= 1) CYCLE
    tmpdep(n) = tmphdxf(n,1)
    DO i=2,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)

    hdx2=0
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
      hdx2 = hdx2 + tmphdxf(n,i)**2 !STEVE: for std dev
      !STEVE: SUBTRACT THE MEAN
      !STEVE: make sure none are zero
      if ( debug_hdxf_0 .AND. tmphdxf(n,i) == 0 ) then
        WRITE(6,*) "================================================================="
        WRITE(6,*) "letkf_obs.f90:: WARNING: tmphdxf(n,i) == 0"
        WRITE(6,*) "n,nobs = ", n,nobs
        WRITE(6,*) "i,nbv = ", i,nbv
        WRITE(6,*) "tmphdxf(n,i) = ", tmphdxf(n,i)
        WRITE(6,*) "tmpdep(n) = ", tmpdep(n)
        WRITE(6,*) "hdx2 = ", hdx2
        WRITE(6,*) "This is later used as a divisor for adaptive inflation."
!       WRITE(6,*) "EXITING... (STOP 10)"
        WRITE(6,*) "tmphdxf(n,:) = ", tmphdxf(n,:)
!       STOP(10)
        WRITE(6,*) "================================================================="
      endif
    END DO
    !STEVE: FORM THE DEPARTURE (OBS INNOVATION)
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx

    ! Use the max of obs and background error to test for compliance of ob
    hdx2 = hdx2 / REAL(nbv-1,r_size) !STEVE: for std dev
    hdx2 = SQRT(hdx2) !STEVE: std dev
    mstd = MAX(tmperr(n),hdx2)
           !STEVE: This is so that if the ensemble spread collapses, then it
           !will still not throw out obs that are outside the spread. Otherwise,
           !the model should grow the spread greater than the obs error.

    qc_obs : IF(ABS(tmpdep(n)) > gross_error*mstd) THEN !gross error
      !tmpqc(n) = 0
      !STEVE: changing this to gradual adjustment method
      ! Rather than removing the observation, increase the obs error so it 
      ! satisfies the QC condition...
      ! STEVE: new feature, 12/31/10
      if (ABS(tmpdep(n)) > 2.0*gross_error*mstd) then
        tmpqc(n) = 0
        gross_2x_cnt = gross_2x_cnt + 1
      else !STEVE: scale up the observational error
        !tmperr(n) = ABS(tmpdep(n))-ABS(hdx2)
        !tmperr(n) = ABS(tmpdep(n))/gross_error
        !!STEVE: more agressive scaling
        !tmperr(n) = (1 + ABS(tmpdep(n))/gross_error)**2 - 1.0 
        ! 9/29/14: Scale up the obs error
        tmperr(n) = tmperr(n) * ABS(tmpdep(n)) / (gross_error*mstd)
        !!STEVE:'just barely satisfying qc' scaling
        gross_cnt = gross_cnt + 1
      endif
        !STEVE: How does adaptive inflation interact with these points?
    ENDIF qc_obs

    !STEVE: as a check, count the number of each type of observation
    if (tmpelm(n) .eq. id_u_obs) cnt_obs_u = cnt_obs_u + 1
    if (tmpelm(n) .eq. id_v_obs) cnt_obs_v = cnt_obs_v + 1
    if (tmpelm(n) .eq. id_t_obs) cnt_obs_t = cnt_obs_t + 1
    if (tmpelm(n) .eq. id_s_obs) cnt_obs_s = cnt_obs_s + 1
    if (tmpelm(n) .eq. id_ssh_obs) cnt_obs_ssh = cnt_obs_ssh + 1
    if (tmpelm(n) .eq. id_eta_obs) cnt_obs_eta = cnt_obs_eta + 1
    if (tmpelm(n) .eq. id_sst_obs) cnt_obs_sst = cnt_obs_sst + 1
    if (tmpelm(n) .eq. id_sss_obs) cnt_obs_sss = cnt_obs_sss + 1
    !(DRIFTERS)
    if (tmpelm(n) .eq. id_x_obs) cnt_obs_x = cnt_obs_x + 1
    if (tmpelm(n) .eq. id_y_obs) cnt_obs_y = cnt_obs_y + 1
    if (tmpelm(n) .eq. id_z_obs) cnt_obs_z = cnt_obs_z + 1

  END DO
  DEALLOCATE(tmpqc0)

endif

  WRITE(6,'(I10,A)') SUM(tmpqc),' OBSERVATIONS TO BE ASSIMILATED'
  !STEVE:
  WRITE(6,*) "cnt_obs_u = ", cnt_obs_u
  WRITE(6,*) "cnt_obs_v = ", cnt_obs_v
  WRITE(6,*) "cnt_obs_t = ", cnt_obs_t
  WRITE(6,*) "cnt_obs_s = ", cnt_obs_s
  WRITE(6,*) "cnt_obs_x = ", cnt_obs_x
  WRITE(6,*) "cnt_obs_y = ", cnt_obs_y
  WRITE(6,*) "cnt_obs_z = ", cnt_obs_z
  WRITE(6,*) "cnt_obs_ssh = ", cnt_obs_ssh
  WRITE(6,*) "cnt_obs_eta = ", cnt_obs_eta
  WRITE(6,*) "cnt_obs_sst = ", cnt_obs_sst
  WRITE(6,*) "cnt_obs_sss = ", cnt_obs_sss
  WRITE(6,*) "gross_cnt = ", gross_cnt
  WRITE(6,*) "gross_2x_cnt = ", gross_2x_cnt

  cnt_obs(iv3d_u) = cnt_obs_u
  cnt_obs(iv3d_v) = cnt_obs_v
  cnt_obs(iv3d_t) = cnt_obs_t
  cnt_obs(iv3d_s) = cnt_obs_s
  cnt_obs(nv3d+iv2d_ssh) = cnt_obs_ssh
  cnt_obs(nv3d+iv2d_eta) = cnt_obs_eta
  cnt_obs(nv3d+iv2d_sst) = cnt_obs_sst
  cnt_obs(nv3d+iv2d_sss) = cnt_obs_sss

  CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)

!STEVE: maybe use this if there are enough observations...
!
! temporal observation localization
!
! If nbslot == 5, multiplier is approximately
! At islot = 1, ~ 1.04
!          = 2, ~ 1.02
!          = 3, ~ 1.01
!          = 4, ~ 1.0025
!          = 5, ~ 1.0
!
! STEVE: watch out if using this for adaptive observation error, this scaling
! may cause problems if applied here...
! STEVE: commenting out because the observations are not having enough of an
! impact. It may be due to the increased error on the observations.
! PLUS, the temporal correlation scales are much longer than 5 days, so we can just ignore this

!  nn = 0
!  DO islot=1,nslots
!    if ( islot .ne. nbslot ) then
!      tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
!                                      & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
!    endif
!    nn = nn + nobslots(islot)
!  END DO

!
! SELECT OBS IN THE NODE
!
  nn = 0
  !STEVE: first, remove all of the Quality-Controlled data
  DO n=1,nobs
    IF(tmpqc(n) /= 1) CYCLE
!    IF(tmplat(n) < MINVAL(lat1) .OR. MAXVAL(lat1) < tmplat(n)) THEN
!      dlat = MIN( ABS(MINVAL(lat1)-tmplat(n)),ABS(MAXVAL(lat1)-tmplat(n)) )
!      IF(dlat > dlat_zero) CYCLE
!    END IF
!    IF(tmplon(n) < MINVAL(lon1) .OR. MAXVAL(lon1) < tmplon(n)) THEN
!      dlon1 = ABS(MINVAL(lon1) - tmplon(n))
!      dlon1 = MIN(dlon1,360.0d0-dlon1)
!      dlon2 = ABS(MAXVAL(lon1) - tmplon(n))
!      dlon2 = MIN(dlon2,360.0d0-dlon2)
!      dlon =  MIN(dlon1,dlon2) &
!         & * pi*re*COS(tmplat(n)*pi/180.d0)/180.0d0
!      IF(dlon > dist_zero) CYCLE
!    END IF
    nn = nn+1
    tmpelm(nn) = tmpelm(n)
    tmplon(nn) = tmplon(n)
    tmplat(nn) = tmplat(n)
    tmplev(nn) = tmplev(n)
    tmpdat(nn) = tmpdat(n)
    tmperr(nn) = tmperr(n)
    tmpk(nn) = tmpk(n)
    tmpdep(nn) = tmpdep(n)
    tmphdxf(nn,:) = tmphdxf(n,:)
    tmpqc(nn) = tmpqc(n)
    tmpid(nn) = tmpid(n)     !(DRIFTERS)
    tmptime(nn) = tmptime(n) !(DRIFTERS)
  END DO
  nobs = nn
  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank

!
! SORT
!
  ALLOCATE( tmp2elm(nobs) )
  ALLOCATE( tmp2lon(nobs) )
  ALLOCATE( tmp2lat(nobs) )
  ALLOCATE( tmp2lev(nobs) )
  ALLOCATE( tmp2dat(nobs) )
  ALLOCATE( tmp2err(nobs) )
  ALLOCATE( tmp2dep(nobs) )
  ALLOCATE( tmp2hdxf(nobs,nbv) )
  ALLOCATE( tmp2id(nobs) )    !(DRIFTERS)
  ALLOCATE( tmp2time(nobs) )  !(DRIFTERS)
  ALLOCATE( obselm(nobs) )
  ALLOCATE( obslon(nobs) )
  ALLOCATE( obslat(nobs) )
  ALLOCATE( obslev(nobs) )
  ALLOCATE( obsdat(nobs) )
  ALLOCATE( obserr(nobs) )
  ALLOCATE( obsdep(nobs) )
  ALLOCATE( obshdxf(nobs,nbv) )
  ALLOCATE( obsid(nobs) )    !(DRIFTERS)
  ALLOCATE( obstime(nobs) )  !(DRIFTERS)

  nobsgrd = 0
  nj = 0
  ! Count the number of observations within each latitude range
  DO j=1,nlat-1
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nj(j) = nj(j) + 1
    END DO
  END DO
  ! Record cumulative sum of observations up to this latitude
  ! Creates the basis for an indexing of observations from lat to lat
  DO j=1,nlat-1
    njs(j) = SUM(nj(0:j-1))
  END DO

  ! Rearrange observations by latitude
  DO j=1,nlat-1
    nn = 0
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
!     IF(tmplon(n) >= lon(nlon)-EPSILON(1.0d0)) CYCLE   !STEVE: I added this to align with the same condition in the code above
                                                        !       Otherwise, sometimes nn /= nj(j)
      nn = nn + 1
      tmp2elm(njs(j)+nn) = tmpelm(n)
      tmp2lon(njs(j)+nn) = tmplon(n)
      tmp2lat(njs(j)+nn) = tmplat(n)
      tmp2lev(njs(j)+nn) = tmplev(n)
      tmp2dat(njs(j)+nn) = tmpdat(n)
      tmp2err(njs(j)+nn) = tmperr(n)
!      tmp2k(njs(j)+nn) = tmpk(n)
      tmp2dep(njs(j)+nn) = tmpdep(n)
      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
      tmp2id(njs(j)+nn) = tmpid(n)     !(DRIFTERS)
      tmp2time(njs(j)+nn) = tmptime(n) !(DRIFTERS)
    END DO
  END DO

  ! For each latitude, identify the number of obs per longitude.
  ! Then, rearrange observations by longitude within each latitude step
  DO j=1,nlat-1
    IF(nj(j) == 0) THEN
      nobsgrd(:,j) = njs(j)
      CYCLE
    END IF
    nn = 0
    DO i=1,nlon
      DO n=njs(j)+1,njs(j)+nj(j)

        ! Find the correct longitude bin for this observation...
        IF(i < nlon) THEN
          IF(tmp2lon(n) < lon(i) .OR. lon(i+1) <= tmp2lon(n)) CYCLE
        ELSE
! STEVE: this is causing nn /= nj(j), the error thrown below.
!        We need these points that are skipped, otherwise there are
!        blank entries in the obselm etc. arrays, and this will
!        lead to problems during the main letkf algorithm.
!        Another solution may be to cut out all the empty entries
!        by changing the obsxxx indicies.
!
          IF(tmp2lon(n) < lon(nlon)) CYCLE

          !STEVE: debugging
          IF(.false.) THEN
            WRITE(6,*) "n, nn, njs(j), nj(j) = ", n, nn, njs(j), nj(j)
            WRITE(6,*) "KEEPING, i == nlon == ", i, nlon
            WRITE(6,*) "tmp2lon(n) = ", tmp2lon(n)
            WRITE(6,*) "lon(nlon) = ", lon(nlon)
            !WRITE(6,*) "either tmp2lon(n) >= lon(nlon) .OR. 360.0d0 > tmp2lon(n)"
            WRITE(6,*) "tmp2lon(n) >= lon(nlon)"
            WRITE(6,*) "========================================================"
          ENDIF
        END IF
        nn = nn + 1
        obselm(njs(j)+nn) = tmp2elm(n)
        obslon(njs(j)+nn) = tmp2lon(n)
        obslat(njs(j)+nn) = tmp2lat(n)
        obslev(njs(j)+nn) = tmp2lev(n)
        obsdat(njs(j)+nn) = tmp2dat(n)
        obserr(njs(j)+nn) = tmp2err(n)
        obsdep(njs(j)+nn) = tmp2dep(n)
        obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
        obsid(njs(j)+nn) = tmp2id(n)     !(DRIFTERS)
        obstime(njs(j)+nn) = tmp2time(n) !(DRIFTERS)
      END DO
      
      ! This now contains the accumulated count of obs up to this lat, up to this lon
      nobsgrd(i,j) = njs(j) + nn
    END DO

    IF(nn /= nj(j)) THEN
      WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(6,'(F6.2,A,F6.2)') lat(j),'<= LAT <',lat(j+1)
      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j))),'<= OBSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j)))
      WRITE(6,*) "j = ", j
      WRITE(6,*) "njs(j) = ", njs(j)
      WRITE(6,*) "nj(j) = ", nj(j)
      !STEVE: this is bad, something is wrong
      WRITE(6,*) "STEVE: this error will cause matrix eigenvalue < 0 error."
      stop 3
    END IF

  END DO

  DEALLOCATE( tmp2elm )
  DEALLOCATE( tmp2lon )
  DEALLOCATE( tmp2lat )
  DEALLOCATE( tmp2lev )
  DEALLOCATE( tmp2dat )
  DEALLOCATE( tmp2err )
  DEALLOCATE( tmp2dep )
  DEALLOCATE( tmp2hdxf )
  DEALLOCATE( tmp2id )    !(DRIFTERS)
  DEALLOCATE( tmp2time )  !(DRIFTERS)
  DEALLOCATE( tmpelm )
  DEALLOCATE( tmplon )
  DEALLOCATE( tmplat )
  DEALLOCATE( tmplev )
  DEALLOCATE( tmpdat )
  DEALLOCATE( tmperr )
  DEALLOCATE( tmpk )
  DEALLOCATE( tmpdep )
  DEALLOCATE( tmphdxf )
  DEALLOCATE( tmpqc )
  DEALLOCATE( tmpid )     !(DRIFTERS)
  DEALLOCATE( tmptime )   !(DRIFTERS)

  RETURN
END SUBROUTINE set_letkf_obs

!-----------------------------------------------------------------------
! Monitor departure from gues/anal mean
!-----------------------------------------------------------------------
SUBROUTINE monit_mean(file)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:) !(nlon,nlat,nv2d)
  REAL(r_size) :: elem
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_ssh,bias_sst,bias_sss !(OCEAN)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_s,rmse_ssh,rmse_sst,rmse_sss !(OCEAN)
  REAL(r_size) :: hdxf,dep,ri,rj,rk
  INTEGER :: n,iu,iv,it,is,issh,isst,isss !(OCEAN)
  CHARACTER(11) :: filename='filexxx.grd'


  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  rmse_u  = 0.0d0
  rmse_v  = 0.0d0
  rmse_t  = 0.0d0
  rmse_s  = 0.0d0  !(OCEAN)
  rmse_ssh = 0.0d0 !(OCEAN)
  rmse_sst = 0.0d0 !(OCEAN)
  rmse_sss = 0.0d0 !(OCEAN)
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_s = 0.0d0   !(OCEAN)
  bias_ssh = 0.0d0 !(OCEAN)
  bias_sst = 0.0d0 !(OCEAN)
  bias_sss = 0.0d0 !(OCEAN)
  iu  = 0
  iv  = 0
  it  = 0
  is  = 0
  issh = 0
  isst = 0
  isss = 0

  WRITE(filename(1:7),'(A4,A3)') file,'_me'
  CALL read_bingrd(filename,v3d,v2d)

  DO n=1,nobs
   !HOLDIT 
    CALL phys2ijk(obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)    !(OCEAN)
    IF(CEILING(rk) > nlev) CYCLE
!   IF(FLOOR(rk) < 1) CYCLE
!   IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ssh_obs) THEN              !(OCEAN)
!     IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
!       rk = 1.00001d0
!     ELSE
!       CYCLE
!     END IF
!   END IF
    IF(CEILING(rk) < 2 .AND. rk < 1.00001d0) THEN
      rk = 1.00001d0
    END IF
    CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,hdxf)                   !(OCEAN)
    dep = obsdat(n) - hdxf
    SELECT CASE(NINT(obselm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep**2
      bias_u = bias_u + dep
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep**2
      bias_v = bias_v + dep
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep**2
      bias_t = bias_t + dep
      it = it + 1
    CASE(id_s_obs)
      rmse_s = rmse_s + dep**2        !(OCEAN)
      bias_s = bias_s + dep           !(OCEAN)
      is = is + 1                     !(OCEAN)
    CASE(id_ssh_obs)                  !(OCEAN)
      rmse_ssh = rmse_ssh + dep**2    !(OCEAN)
      bias_ssh = bias_ssh + dep       !(OCEAN)
      issh = issh + 1                 !(OCEAN)
    CASE(id_sst_obs)                  !(OCEAN)
      rmse_sst = rmse_sst + dep**2    !(OCEAN)
      bias_sst = bias_sst + dep       !(OCEAN)
      isst = isst + 1                 !(OCEAN)
    CASE(id_sss_obs)                  !(OCEAN)
      rmse_sss = rmse_sss + dep**2    !(OCEAN)
      bias_sss = bias_sss + dep       !(OCEAN)
      isss = isss + 1                 !(OCEAN)
    END SELECT
  END DO

  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(is == 0) THEN                            !(OCEAN)
    rmse_s = undef                            !(OCEAN)
    bias_s = undef                            !(OCEAN)
  ELSE
    rmse_s = SQRT(rmse_s / REAL(is,r_size))   !(OCEAN)
    bias_s = bias_s / REAL(is,r_size)         !(OCEAN)
  END IF

  IF(issh == 0) THEN                               !(OCEAN)
    rmse_ssh = undef                               !(OCEAN)
    bias_ssh = undef                               !(OCEAN)
  ELSE
    rmse_ssh = SQRT(rmse_ssh / REAL(issh,r_size))  !(OCEAN)
    bias_ssh = bias_ssh / REAL(issh,r_size)        !(OCEAN)
  END IF
  IF(isst == 0) THEN                               !(OCEAN)
    rmse_sst = undef                               !(OCEAN)
    bias_sst = undef                               !(OCEAN)
  ELSE
    rmse_sst = SQRT(rmse_sst / REAL(isst,r_size))  !(OCEAN)
    bias_sst = bias_sst / REAL(isst,r_size)        !(OCEAN)
  END IF
  IF(isss == 0) THEN                               !(OCEAN)
    rmse_sss = undef                               !(OCEAN)
    bias_sss = undef                               !(OCEAN)
  ELSE
    rmse_sss = SQRT(rmse_sss / REAL(isss,r_size))  !(OCEAN)
    bias_sss = bias_sss / REAL(isss,r_size)        !(OCEAN)
  END IF

  WRITE(6,'(3A)') '== PARTIAL OBSERVATIONAL DEPARTURE (',file,') ========================================'
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','SST','SSS'                          !(OCEAN)
  WRITE(6,'(7ES12.3)') bias_u,bias_v,bias_t,bias_s,bias_ssh,bias_sst,bias_sss  !(OCEAN)
  WRITE(6,'(7ES12.3)') rmse_u,rmse_v,rmse_t,rmse_s,rmse_ssh,rmse_sst,rmse_sss  !(OCEAN)
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS ========================================================'
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','SST','SSS'                          !(OCEAN)
  WRITE(6,'(7I12)') iu,iv,it,is,issh,isst,isss                                 !(OCEAN)
  WRITE(6,'(A)') '=================================================================================='

  DEALLOCATE(v3d,v2d)

  RETURN

END SUBROUTINE monit_mean

SUBROUTINE get_hdxa(anal3dg,anal2dg,hdxa) !,depa)
REAL(r_size), INTENT(IN) :: anal3dg(nlon,nlat,nlev,nv3d)
REAL(r_size), INTENT(IN) :: anal2dg(nlon,nlat,nv2d)
REAL(r_size),INTENT(OUT) :: hdxa(nobs)
!REAL(r_size),INTENT(OUT) :: depa(nobs)
REAL(r_size),ALLOCATABLE :: wk2d(:,:)
REAL(r_size),ALLOCATABLE :: wk1d(:)
REAL(r_size) :: ri,rj,rk
INTEGER :: n,i,ierr,obsper,start,finish
INTEGER :: prntmod=HUGE(1)
LOGICAL :: dodebug = .true.

hdxa = 0.0d0

! l=0
! DO
!   im = myrank+1 + nprocs * l
!   IF(im > nbv) EXIT

!   WRITE(analfile(5:7),'(I3.3)') im
!   WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',analfile
!   CALL read_grd(analfile,v3d,v2d)

!   WRITE(6,*) "letkf_tools.f90::adapt_obserr: calling get_hdxa..."
!   CALL get_hdxa(v3d,v2d,hdxa(:,im),obsdep_a)
!     
!   l = l+1
! END DO

WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is computing Trans_XtoY for get_hdxa'
WRITE(6,*) "myrank+1,nobs,nprocs = ", myrank+1,nobs,nprocs
obsper = CEILING(REAL(nobs)/REAL(nprocs))
if (myrank + 1 .lt. nprocs)  then
  start = myrank*obsper+1
  finish = (myrank+1)*obsper 
else
  start = myrank*obsper+1
  finish = nobs
endif
WRITE(6,*) "obsper = ", obsper
WRITE(6,*) "start,finish = ", start, finish

if (dodebug) prntmod = NINT(obsper/2.0)

! Process all obs allocated to this processor
DO n=start,finish
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "n = ", n 
  ! interpolation
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "Calling phys2ijk..."
  CALL phys2ijk(obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "ri,rj,rk = ", ri,rj,rk
  ! observational operator
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "Calling Trans_XtoY..."
  CALL Trans_XtoY(obselm(n),ri,rj,rk,anal3dg,anal2dg,hdxa(n))
  if (MOD(n,prntmod) .eq. 0) then
    WRITE(6,*) "anal3d(FLOOR(ri),FLOOR(rj),FLOOR(rk),:) = ", anal3dg(FLOOR(ri),FLOOR(rj),FLOOR(rk),:)
    WRITE(6,*) "anal3d(CEILING(ri),FLOOR(rj),FLOOR(rk),:) = ", anal3dg(CEILING(ri),FLOOR(rj),FLOOR(rk),:)
    WRITE(6,*) "anal3d(FLOOR(ri),CEILING(rj),FLOOR(rk),:) = ", anal3dg(FLOOR(ri),CEILING(rj),FLOOR(rk),:)
    WRITE(6,*) "anal3d(FLOOR(ri),FLOOR(rj),CEILING(rk),:) = ", anal3dg(FLOOR(ri),FLOOR(rj),CEILING(rk),:)
    WRITE(6,*) "anal3d(CEILING(ri),CEILING(rj),FLOOR(rk),:) = ", anal3dg(CEILING(ri),CEILING(rj),FLOOR(rk),:)
    WRITE(6,*) "anal3d(CEILING(ri),FLOOR(rj),CEILING(rk),:) = ", anal3dg(CEILING(ri),FLOOR(rj),CEILING(rk),:)
    WRITE(6,*) "anal3d(FLOOR(ri),CEILING(rj),CEILING(rk),:) = ", anal3dg(FLOOR(ri),CEILING(rj),CEILING(rk),:)
    WRITE(6,*) "anal3d(CEILING(ri),CEILING(rj),CEILING(rk),:) = ", anal3dg(CEILING(ri),CEILING(rj),CEILING(rk),:)
  endif
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "hdxa(n) = ", hdxa(n)
  if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "obsdat(n) = ", obsdat(n)
! if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "pre-calc: depa(n) = ", depa(n)

! if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "Calculating dep..."
! depa(n) = obsdat(n) - hdxa(n) ! y-Hx
! if (MOD(n,prntmod) .eq. 0) WRITE(6,*) "depa(n) = ", depa(n)
ENDDO

WRITE(6,*) "MPI_BARRIER"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
WRITE(6,*) "ALLOCATE(wk1d(nobs))..."
ALLOCATE(wk1d(nobs))
wk1d=hdxa
WRITE(6,*) "MPI_BARRIER"
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
WRITE(6,*) "MPI_ALLREDUCE"
CALL MPI_ALLREDUCE(wk1d,hdxa,nobs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
WRITE(6,*) "Calling MPI_BARRIER..."
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
WRITE(6,*) "finished MPI_ALLREDUCE"
DEALLOCATE(wk1d)

!WRITE(6,*) "MPI_BARRIER"
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!WRITE(6,*) "ALLOCATE(wk1d(nobs))..."
!ALLOCATE(wk1d(nobs))
!wk1d = depa
!WRITE(6,*) "MPI_BARRIER"
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!WRITE(6,*) "MPI_ALLREDUCE"
!CALL MPI_ALLREDUCE(wk1d,depa,nobs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!WRITE(6,*) "MPI_BARRIER"
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!DEALLOCATE(wk1d)

!STEVE: may need this if I want to do for full ensemble:
!DO n=1,nobs
!  depa(n) = hdxa(n,1)
!  DO i=2,nbv
!    depa(n) = depa(n) + hdxa(n,i)
!  END DO
!  depa(n) = depa(n) / REAL(nbv,r_size)
!  DO i=1,nbv
!    hdxa(n,i) = hdxa(n,i) - depa(n) ! Hdx
!    !STEVE: make sure none are zero
!    if ( debug_hdxf_0 .AND. hdxa(n,i) == 0 ) then
!      WRITE(6,*) "get_hdxa:: WARNING: hdxa(n,i) == 0"
!      WRITE(6,*) "This is later used as a divisor for adaptive inflation."
!    endif
!  END DO
!  depa(n) = obsdat(n) - depa(n) ! y-Hx
!ENDDO

END SUBROUTINE get_hdxa

END MODULE letkf_obs
