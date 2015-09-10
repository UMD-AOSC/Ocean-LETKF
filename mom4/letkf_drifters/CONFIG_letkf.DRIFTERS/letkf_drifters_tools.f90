MODULE letkf_drifters_tools
! Original Author: Prof. Stephen G. Penny
!                  Department of Atmospheric and Oceanic Science
!                  University of Maryland College Park
! Edits:
!        06/29/2015: moved code to UMD DT2 (STEVE), compiled successfully
!        07/02/2015: moved code back to Gaea to synch with GitHub repository
!
! Additional Edits:
!
!
!
! Function:
! Normally, when analyzing the model grid, the gridpoints themselves are the
! model space while the data values are the temperature, salinity, u and v
! velocities, etc. For the drifters, it is the drifter_ids that are the model
! space, and the data values are the positions themselves, as specified by the
! x, y and z-coordinates.
!
!
USE netcdf
USE common
USE common_mpi  
USE common_letkf    !nbv := number of ensemble members
USE letkf_obs, ONLY : nobs
USE letkf_drifters_local, ONLY : obs_local  !(DRIFTERS) NOTE: this is a slightly modified version of letkf_local.f90
USE params_letkf
USE params_model
USE params_obs

!STEVE: for (DRIFTERS)
REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: v4d          !(DRIFTERS) num_drifters x num_times x nbv x nv4d 
                                                                    !(DRIFTERS) value stored is coordinate
REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_idx        !(DRIFTERS) dim=nid1, this is read in and stored during scatter
REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_ids        !(DRIFTERS) dim=num_drifters, this is read in from model output
REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_times      !(DRIFTERS) dim=num_drifters, read in from model output
REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_ti         !(DRIFTERS) dim=num_drifters, time index (REAL)
INTEGER, SAVE :: num_drifters !read in from model output, updated by observation input, used to dimension prev. arrays
INTEGER, SAVE :: num_times    !read in from model output, updated by observation input, used to dimension prev. arrays
!LOGICAL, PARAMETER :: DO_DRIFTERS = .false. !STEVE: this is in common_mom4.f90

! Inflation parameter for (DRIFTERS)
!REAL(r_size),PARAMETER :: cov_infl_mul = 1.01d0 !-1.01d0 !multiplicative inflation

! For parallel mpi operations:
INTEGER,SAVE :: nid1                             !STEVE: this is the number of gridpoints to run on this (myrank) processor
INTEGER,SAVE :: nid1max                          !STEVE: the largest number of gridpoints on any 1 processor
INTEGER,ALLOCATABLE,SAVE :: nid1node(:)

! Variable localization:
! INTEGER, PARAMETER :: nid_dobs = nv4d
REAL(r_size),PARAMETER :: var_local(nv4d,nid_dobs) = RESHAPE( &
!           X      Y      Z
   & (/ 1.0d0, 1.0d0, 1.0d0,  & ! X
   &    1.0d0, 1.0d0, 1.0d0,  & ! Y
   &    1.0d0, 1.0d0, 1.0d0 /)& ! Z
   & ,(/nv4d,nid_dobs/))

CONTAINS

! The first routine called by letkf.f90:
SUBROUTINE set_common_drifters
  IMPLICIT NONE
  INTEGER :: i,j,n !,ibv 
  CHARACTER(32) :: drifile

  WRITE(6,'(A)') 'Hello from set_common_drifters'

  !test_drifters:
  num_drifters=5
  num_times=45
  !ndr0=num_drifters*num_times

  ALLOCATE(drifter_times(num_times))
  ALLOCATE(drifter_ids(num_drifters))

  !alternative: read this from a file

  ! Next, split up the drifters to the different processors:
  i = MOD(num_drifters,nprocs)
  nid1max = (num_drifters - i)/nprocs + 1
  WRITE(6,*) "nid1max = ", nid1max
  IF(myrank < i) THEN
    nid1 = nid1max
  ELSE
    nid1 = nid1max - 1
  END IF
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of drifter ids: nid1=',nid1
  ALLOCATE(nid1node(nprocs))
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nid1node(n) = nid1max
    ELSE
      nid1node(n) = nid1max - 1
    END IF
  END DO
  WRITE(6,*) "nid1node...."
  WRITE(6,*) nid1node

  ! For this processor, allocate enough space for its assigned drifters
  !ALLOCATE(v4d(nid1max,num_times,nbv,nv4d))

  ! Read the ensemble of drifter files into this data structure
  !CALL read_ens_drifters('gues',v4d)
  RETURN

END SUBROUTINE set_common_drifters

SUBROUTINE read_obs2_drifters(cfile,nn,elem,rlon,rlat,rlev,obid,oerr,otime)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn  ! LUYU: number of drifters
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: otime(nn)
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number: id_x_obs, id_y_obs, id_z_obs
  INTEGER,INTENT(OUT) :: obid(nn)
  REAL(r_sngl) :: wk(7)
  INTEGER :: n,iunit  

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  !read(iunit,*)  ! LUYU: necessary, if the data format is the same as the one of background data.
  !read(iunit,*) 
  DO n=1,nn
    READ(iunit) wk
!   SELECT CASE(NINT(wk(1)))
!   CASE(id_u_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_v_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_t_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_q_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!  END SELECT
    elem(n) = INT(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    obid(n) = INT(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    otime(n) = REAL(wk(7),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs2_drifters

SUBROUTINE write_obs2_drifters(cfile,nn,elem,rlon,rlat,rlev,obid,oerr,ohx,oqc,otime)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: ohx(nn)
  INTEGER,INTENT(IN) :: obid(nn)
  INTEGER,INTENT(IN) :: oqc(nn)
  REAL(r_size),INTENT(IN) :: otime(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit
  LOGICAL :: dodebug=.false.

  iunit=92
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    wk(1) = REAL(elem(n),r_sngl)  ! Type of observation x,y,z
    wk(2) = REAL(rlon(n),r_sngl)  ! Ob lon
    wk(3) = REAL(rlat(n),r_sngl)  ! Ob lat
    wk(4) = REAL(rlev(n),r_sngl)  ! Ob level
    wk(5) = INT(obid(n))  ! ID for the drifters
    wk(6) = REAL(oerr(n),r_sngl)  ! Estimated observation error
    wk(7) = REAL(ohx(n),r_sngl)   ! Model forecast transformed to observation space: H(xb)
    wk(8) = INT(oqc(n))   ! Quality control ID (1==keep, 0==discard) for use in assimilation
    wk(9) = REAL(otime(n),r_sngl)  ! Ob time
    if (dodebug) PRINT '(I6,2F7.2,F10.2,4F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN

END SUBROUTINE write_obs2_drifters

! The second routine called by letkf.f90:
SUBROUTINE das_drifters(gues4d,anal4d)
  ! das_letkf went through all model grid points. Now we're going to go through
  ! each drifter_id, as if we had appended them to the model state vector.
  
  USE common_letkf, ONLY: letkf_core
  USE params_model

  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues4d(nid1,num_times,nbv,nv4d) ! background ensemble
  REAL(r_size),INTENT(OUT) :: anal4d(nid1,num_times,nbv,nv4d) ! analysis ensemble
  REAL(r_size),ALLOCATABLE :: mean4d(:,:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv4d)
  LOGICAL :: ex
  INTEGER :: id,it,n,m,i,j,k,nobsl,ierr
  !STEVE: for debugging
  LOGICAL :: dodebug = .false.
  INTEGER :: nn
  !STEVE: added for drifters:
  INTEGER :: nobstotal,itim

  WRITE(6,'(A)') 'Hello from das_drifters'
  nobstotal = nobs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs

  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    anal4d = gues4d
    RETURN
  ELSE                   !(OCEAN)
    anal4d = 0.0d0       !(OCEAN)
  END IF

  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean4d(nid1,num_times,nv4d))
  CALL ensmean_drifters(nbv,nid1,gues4d,mean4d)
  !STEVE: nid1 is each one of the grid points

  DO n=1,nv4d
    DO m=1,nbv
      DO k=1,num_times
        DO i=1,nid1
          gues4d(i,k,m,n) = gues4d(i,k,m,n) - mean4d(i,k,n)
        END DO
      END DO
    END DO
  END DO

  !
  ! multiplicative inflation
  !
  !STEVE: (using simple constant inflation for now.)

  !
  ! MAIN ASSIMILATION LOOP
  !
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE(hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  WRITE(6,*) "... done."

  DO it=1,num_times !STEVE: go through every possible coordinate of the grid in list form...
    if (dodebug) WRITE(6,*) "it = ", it

    DO id=1,nid1
      if (dodebug) WRITE(6,*) "id= ", id
      itim=NINT(drifter_ti(it)) ! The index of the exact time of the observation

    ! For each coordinate, x,y,and z:
    DO n=1,nv4d
      ! Find the observations around this point.

      CALL obs_local(id,it,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)
      parm = cov_infl_mul !STEVE: keeping it simple

      CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))

      ! The ":" in place of "itim" implies that all times are affected equally by
      ! the observed drifter location.
      DO m=1,nbv
        anal4d(id,it,m,n) = mean4d(id,it,n)
        DO k=1,nbv
          anal4d(id,it,m,n) = anal4d(id,it,m,n) + gues4d(id,it,k,n) * trans(k,m,n)
        END DO
      END DO

    ENDDO

    ENDDO

 ENDDO

  DEALLOCATE(hdxf,rdiag,rloc,dep)

! If there are observations that weren't on the model grid, add them to the
! output list.

END SUBROUTINE das_drifters

!-----------------------------------------------------------------------
! Transformation from model (drifter_id) space to observation space (i.e. H-operator)
!-----------------------------------------------------------------------
SUBROUTINE Trans_DtoY(elm,ri,rj,rk,v4d_in,yobs)        !(OCEAN)
  USE common_obs_mom4
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_sngl),INTENT(IN) :: v4d_in(num_drifters,num_times,nv4d)
  REAL(r_size),INTENT(OUT) :: yobs  
  INTEGER :: i,j,k
  INTEGER :: intelm

  intelm = NINT(elm)
  SELECT CASE (intelm)
  CASE(id_x_obs)   ! (OCEAN) (DRIFTERS)
    !STEVE: interpolate the time-based drifter information to the observation time.
    i = NINT(ri)
    k = iv4d_x
    CALL itpl_4d(v4d_in(i,:,k),rj,yobs) !STEVE: here, rj == rt (REAL time index), rk == x,y, or z (1,2,3), and ri == drifter_id
  CASE(id_y_obs)   ! (OCEAN) (DRIFTERS)
    i = NINT(ri)
    k = iv4d_y
    CALL itpl_4d(v4d_in(i,:,k),rj,yobs)
  CASE(id_z_obs)   ! (OCEAN) (DRIFTERS)
    i = NINT(ri)
    k = iv4d_z
    CALL itpl_4d(v4d_in(i,:,k),rj,yobs)
  CASE DEFAULT
    print *, "ERROR::Trans_DtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_x_obs, id_y_obs, id_z_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

  RETURN
END SUBROUTINE Trans_DtoY

SUBROUTINE drift2ijk(elem,oid,otime,rlon,rlat,rlev,ri,rj,rk)     !(OCEAN)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  INTEGER     ,INTENT(IN) :: oid
  REAL(r_size),INTENT(IN) :: otime 
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! depth
  REAL(r_size),INTENT(OUT) :: ri  ! drifter id index
  REAL(r_size),INTENT(OUT) :: rj  ! time index
  REAL(r_size),INTENT(OUT) :: rk  ! var index (x,y,z)
  REAL(r_size) :: aj
  INTEGER :: i,j,k
  LOGICAL :: dodebug = .false.

  !(DRIFTERS) num_times x nv4d x num_drifters
  ! If this is a drifter observation, then what we need is a time interpolation:
  if (elem .eq. id_x_obs) then   ! ??, id_x_obs or iv4d_x
    rk = 1
  elseif (elem .eq. id_y_obs) then  ! ??, id_y_obs
    rk = 2
  elseif (elem .eq. id_z_obs) then  ! ??, id_z_obs
    rk = 3
  else
    WRITE(6,*) "This observation is not a drifter! Stopping..."
    STOP(3)
  endif

  ! Identify drifter and time indices:
  ! Loop through drifter id's until we find the right drifter:
  do k=1,num_drifters
    if (drifter_ids(i) .eq. oid) exit
    ! drifter_ids is an array of ids corrresponding the the ids in the
    ! observations. This must be input via the model input file, and the
    ! drifter ids in the model must match the observation ids.
  enddo
  ri = i

  ! Loop through time indexes until we find the right time index and then get
  ! the fractional time index:
  do j=1,num_drifters
    if (drifter_times(j) > otime) exit
  enddo
  !STEVE: now apply the interpolation at the identified drifter id
  aj = (otime - drifter_times(j-1)) / (drifter_times(j) - drifter_times(j-1))
  rj = REAL(j-1,r_size) + aj

  return
END SUBROUTINE drift2ijk

!STEVE: for assimilating (DRIFTERS)
! NOTE: input data for one drifter id, at all times
SUBROUTINE itpl_4d(var,rt,var5)
  IMPLICIT NONE
  REAL(r_sngl),INTENT(IN) :: var(:) ! dim = num_times
! REAL(r_size),INTENT(IN) :: var(:,:) ! num_times,3 = (nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: rt ! model time matched to observation timestamp
  REAL(r_size),INTENT(OUT) :: var5
! REAL(r_size),DIMENSION(3),INTENT(OUT) :: var5
  REAL(r_size) :: an
  INTEGER :: n

  ! time index: 
  n = CEILING(rt) 
  an = rt - REAL(n-1,r_size)

  ! Linearly interpolate drifter position information between two drifter times
  var5 = (1-an) * var(n-1) + (an) * var(n) ! x, y, or z-coordinate, depending on which is input
  
  RETURN
END SUBROUTINE itpl_4d

SUBROUTINE read_ens_drifters(file,v4d_sub)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file !(gues)
  REAL(r_size),INTENT(OUT) :: v4d_sub(nid1,num_times,nbv,nv4d)
  REAL(r_sngl) :: v4d_all(num_drifters,num_times,nv4d)
  INTEGER :: l,n,ll,im,imm
  CHARACTER(11) :: filename='file000'
  !STEVE: debug
  LOGICAL :: dodebug = .true.

  !STEVE: for my purposes, nbv << nprocs, so ll==1 and this loop just reads the
  !drifters for one ensemble member. Then it must distribute to all procs.
  ll = CEILING(REAL(nbv)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs  
    IF(im <= nbv) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'In common_mpi_mom4.f90::read_ens_drifters, MYRANK',myrank,' is reading a file ',filename, 'nprocs=', nprocs
      CALL read_drifters(filename,v4d_all)
    END IF

    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= nbv) THEN
        WRITE(6,*) "In common_mpi_mom4.f90::read_ens_drifters, calling scatter_drifters_mpi..."
        CALL scatter_drifters_mpi(n,v4d_all,v4d_sub(:,:,im,:))
        WRITE(6,*) "In common_mpi_mom4.f90::read_ens_drifters, finished calling scatter_drifters_mpi."
      END IF
    END DO
  END DO

  RETURN

END SUBROUTINE read_ens_drifters

SUBROUTINE scatter_drifters_mpi(nrank,v4d_all,v4d_sub)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN)  :: v4d_all(num_drifters,num_times,nv4d)
  REAL(r_size),INTENT(OUT) :: v4d_sub(nid1,num_times,nv4d)
  REAL(r_sngl) :: bufs(nid1max,num_times*nv4d,nprocs)
  REAL(r_sngl) :: bufr(nid1max,num_times*nv4d)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ns = nid1max * (num_times * nv4d)
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv4d
      DO k=1,num_times
        j = j+1
        CALL drft_to_buf(v4d_all(:,k,n),bufs(:,j,:))
      END DO
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  j=0
  DO n=1,nv4d
    DO k=1,num_times
      j = j+1
      v4d_sub(:,k,n) = REAL(bufr(1:nid1,j),r_size)
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN

END SUBROUTINE scatter_drifters_mpi

!-----------------------------------------------------------------------
! drifter id data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE drft_to_buf(drft,buf)
  REAL(r_sngl),INTENT(IN) :: drft(num_drifters)
  REAL(r_sngl),INTENT(OUT) :: buf(nid1max,nprocs)
  INTEGER :: i,j,m,nid

  DO m=1,nprocs
    DO i=1,nid1node(m)
      j = m-1 + nprocs * (i-1)
      nid = MOD(j,num_drifters) + 1
      buf(i,m) = drft(nid)
    END DO
  END DO

  RETURN
END SUBROUTINE drft_to_buf

SUBROUTINE read_dimension(file,num_drifters,num_times)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: file 
  INTEGER, INTENT(OUT) :: num_drifters,num_times
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: fid = 33
  CHARACTER(24) :: drfile

  drfile = trim(file)//'.drifters_inp.txt' ! (DRIFTERS)

  print *, 'Hello from read_dimension.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  close(fid)

  RETURN
END SUBROUTINE read_dimension

SUBROUTINE read_drifters(file,v4d_all)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: file 
  REAL(r_sngl),INTENT(OUT) :: v4d_all(num_drifters,num_times,nv4d)
  REAL :: dlon, dlat, ddepth, dtemp, dsalt, dtime
  INTEGER :: ditime, dids
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: di, ti
  INTEGER :: fid = 33, dodebug = 0
  CHARACTER(24) :: drfile

  drfile = trim(file)//'.drifters_inp.txt' ! (DRIFTERS)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *, 'Hello from read_drifters.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  read(fid,*) header_line

  ! Read all positions (and possibly temp and salt)  of each drifter:
  DO di=1,num_drifters
    DO ti=1,num_times
      read(fid,'(I12,6F12.4,I12)') dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      drifter_ids(di) = dids
      IF (dodebug .eq. 1) THEN
        print *, dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      END IF
      IF (di .eq. 1) drifter_times(ti) = dtime  !
      v4d_all(di,ti,1) = dlon
      v4d_all(di,ti,2) = dlat
      v4d_all(di,ti,3) = ddepth
      IF (nv4d .ge. 5) THEN
        ! If we have temperature and slinity observations at each position,
        ! we can assimilate this data too
        v4d_all(di,ti,4) = dtemp
        v4d_all(di,ti,5) = dsalt
      END IF
    END DO
  END DO
  close(fid)

END SUBROUTINE read_drifters

SUBROUTINE write_drifters(file,v4d_all)
  USE netcdf
  IMPLICIT NONE
! (DRIFTERS)
  CHARACTER(*),INTENT(IN) :: file
  REAL(r_sngl), INTENT(IN) :: v4d_all(num_drifters,num_times,nv4d)
  INTEGER :: nd, np
  INTEGER :: nd_dimid, np_dimid, dimids(2)
  INTEGER :: pos_varid, ids_varid
  REAL, DIMENSION(:,:), ALLOCATABLE :: positions
  INTEGER, DIMENSION(:), ALLOCATABLE :: ids
  CHARACTER(24) :: drfile
  INTEGER :: ncid

  drfile = trim(file)//'.drifters_inp.nc' ! (DRIFTERS)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write drifter positions in input file readable by mom4p1

  !Collect positions and times of drifters
  if (.false.) then ! FOR TESTING
    print *, "letkf_drifters.f90::write_drifters:: WARNING! OUTPUTTING only TEST SAMPLE DATA."
    nd = 3
    np = 5
    ALLOCATE(positions(nd,np),ids(np))
    positions = RESHAPE((/200, 10, 1000, &
                200,  0, 1000, &
                200,-10, 1000, &
                180,  5, 1000/), (/nd,np/) )
    ids = (/12, 23, 34, 45/)
  else
    nd = nv4d
    np = num_drifters
    ALLOCATE(positions(nd,np),ids(np))
    WRITE(6,*) "#### BEFORE####"
    WRITE(6,*) v4d_all(1:np,num_times,1:nd)
    positions = RESHAPE(TRANSPOSE(v4d_all(1:np,num_times,1:nd)),(/nd,np/)) 
    WRITE(6,*) "#### AFTER####"
    WRITE(6,*) positions
    ids = drifter_ids
  endif

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  call check( nf90_create(drfile, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "nd", nd, nd_dimid) )
  call check( nf90_def_dim(ncid, "np", np, np_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ nd_dimid, np_dimid /)

  ! Define the variable.
  call check( nf90_def_var(ncid, "positions", NF90_DOUBLE, dimids, pos_varid))

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, pos_varid, "names", "lon lat depth") )
  call check( nf90_put_att(ncid, pos_varid, "units", "deg_E deg_N meters") )

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).
  call check( nf90_def_var(ncid, "ids", NF90_INT, np_dimid, ids_varid) )

  ! Add global attributes with NF90_GLOBAL
  call check( nf90_put_att(ncid, NF90_GLOBAL, "velocity_names", "u v w") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_names", "lon lat depth temp salt") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_units", "deg_E deg_N meters Celsius PSU") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "time_units", "seconds") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "title", "LETKF analyzed positions for drifters, for input into MOM4p1") )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the data to the file.
  call check( nf90_put_var(ncid, pos_varid, positions) )
  call check( nf90_put_var(ncid, ids_varid, ids) )

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  ! DEALLOCATE
  DEALLOCATE(positions,ids)

END SUBROUTINE write_drifters

!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_drifters(file,v4d_sub)
! INCLUDE 'netcdf.inc' !STEVE: for NaN correction (OCEAN)
  use netcdf
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v4d_sub(nid1,num_times,nbv,nv4d)
  REAL(r_sngl) :: v4d_all(num_drifters,num_times,nv4d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000'
  INTEGER :: i,j,k,m !STEVE: for debugging
  LOGICAL :: verbose = .true.
  INTEGER :: convcnt = 0

  ll = CEILING(REAL(nbv)/REAL(nprocs))
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= nbv) THEN
        WRITE(6,*) "In letkf_drifters.f90::write_ens_drifters, calling gather_drifters_mpi..."
        CALL gather_drifters_mpi(n,v4d_sub(:,:,im,:),v4d_all)
        WRITE(6,*) "In letkf_drifters.f90::write_ens_drifters, finished calling gather_drifters_mpi."
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= nbv) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing file: ',filename

      !STEVE: debug
      WRITE(6,*) "letkf_drifters.f90::write_ens_drifters:: MAXVAL(ABS(v4d(:,:,:,iv4d_x))) = ", MAXVAL(ABS(v4d_all(:,:,iv4d_x)))
      WRITE(6,*) "letkf_drifters.f90::write_ens_drifters:: MAXVAL(ABS(v4d(:,:,:,iv4d_y))) = ", MAXVAL(ABS(v4d_all(:,:,iv4d_y)))
      WRITE(6,*) "letkf_drifters.f90::write_ens_drifters:: MAXVAL(ABS(v4d(:,:,:,iv4d_z))) = ", MAXVAL(ABS(v4d_all(:,:,iv4d_z)))

      CALL write_drifters(filename,v4d_all)
    END IF
  END DO

  RETURN
END SUBROUTINE write_ens_drifters

SUBROUTINE gather_drifters_mpi(nrank,v4d,v4d_all)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v4d(nid1,num_times,nv4d)
  REAL(r_sngl),INTENT(OUT) :: v4d_all(num_drifters,num_times,nv4d)
  REAL(r_sngl) :: bufs(nid1max,num_times*nv4d)
  REAL(r_sngl) :: bufr(nid1max,num_times*nv4d,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nid1max * (num_times*nv4d)
  nr = ns
  j=0
  DO n=1,nv4d
    DO k=1,num_times
      j = j+1
      bufs(1:nid1,j) = REAL(v4d(:,k,n),r_sngl)
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv4d
      DO k=1,num_times
        j = j+1
        CALL buf_to_drft(bufr(:,j,:),v4d_all(:,k,n))
      END DO
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_drifters_mpi

!-----------------------------------------------------------------------
! buffer -> drifters list
!-----------------------------------------------------------------------
SUBROUTINE buf_to_drft(buf,drft)
  REAL(r_sngl),INTENT(IN) :: buf(nid1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: drft(num_drifters)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nid1node(m)
      j = m-1 + nprocs * (i-1)
      nid = MOD(j,num_drifters) + 1
      drft(nid) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_drft

!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_drifters(file,v4d_sub)
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(IN) :: v4d_sub(nid1,num_times,nbv,nv4d)
  REAL(r_size) :: v4dm(nid1,num_times,nv4d)
  REAL(r_size) :: v4ds(nid1,num_times,nv4d)
  REAL(r_sngl) :: v4d_all(num_drifters,num_times,nv4d)
  INTEGER :: i,k,m,n,j
  CHARACTER(11) :: filename='file000.grd'
  
  WRITE(6,*) "Hello from write_ensmspr_drifters..."
  WRITE(6,*) "Calling ensmean_drifters....."
  CALL ensmean_drifters(nbv,nid1,v4d_sub,v4dm)
  WRITE(6,*) "Finished calling ensmean_drifters and start calling gather_drifters_mpi..."

  CALL gather_drifters_mpi(0,v4dm,v4d_all)
  WRITE(6,*) "Finished calling gather_drifters_mpi..."
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_bindrf4(filename,v4d_all)
  END IF

  DO n=1,nv4d
    DO k=1,num_times
      DO i=1,nid1
        v4ds(i,k,n) = (v4d_sub(i,k,1,n)-v4dm(i,k,n))**2
        DO m=2,nbv
          v4ds(i,k,n) = v4ds(i,k,n) + (v4d_sub(i,k,m,n)-v4dm(i,k,n))**2
        END DO
        v4ds(i,k,n) = SQRT(v4ds(i,k,n) / REAL(nbv-1,r_size))
      END DO
    END DO
  END DO

  CALL gather_drifters_mpi(0,v4ds,v4d_all)
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_bindrf4(filename,v4d_all)
  END IF

  RETURN
END SUBROUTINE write_ensmspr_drifters

SUBROUTINE write_bindrf4(filename,v4d)
!===============================================================================
! Write out an letkf grd-format binary file in single precision
!===============================================================================
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v4d(num_drifters,num_times,nv4d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL :: dodebug=.false.

  if (dodebug) print *, "write_bindrf4:: open filename = ",filename
  iunit=56
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_bindrf4:: nij,iolength = ", num_drifters*num_times,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=num_drifters*num_times*iolen)

  irec=1

  DO n=1,nv4d
    if (dodebug) print *, "write_bindrf4:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v4d(i,j,n),i=1,num_drifters),j=1,num_times)
    irec = irec + 1
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_bindrf4

!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_drifters(mem,nid,v4d,v4dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mem ! number of ensemble members
  INTEGER,INTENT(IN) :: nid
  REAL(r_size),INTENT(IN) :: v4d(nid,num_times,mem,nv4d)
  REAL(r_size),INTENT(OUT) :: v4dm(nid,num_times,nv4d)
  INTEGER :: i,k,m,n

  DO n=1,nv4d
    DO k=1,num_times
      DO i=1,nid
        v4dm(i,k,n) = v4d(i,k,1,n)
        DO m=2,mem
          v4dm(i,k,n) = v4dm(i,k,n) + v4d(i,k,m,n)
        END DO
        v4dm(i,k,n) = v4dm(i,k,n) / REAL(nbv,r_size)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_drifters

SUBROUTINE check(status)
!===============================================================================
! Check the error status of the netcdf command
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
END SUBROUTINE check


END MODULE letkf_drifters_tools
