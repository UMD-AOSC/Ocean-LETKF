MODULE common_mpi_oceanmodel
!=======================================================================
!
! [PURPOSE:] MPI procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!   04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   09/11/2015 Steve Penny converted back for use with MOM4 due to MPI upgrades from MOM6 version
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_oceanmodel

  IMPLICIT NONE

  PUBLIC

  INTEGER,SAVE :: nij1                  !STEVE: this is the number of gridpoints to run on this (myrank) processor
  INTEGER,SAVE :: nij1max               !STEVE: the largest number of gridpoints on any 1 processor
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: kmt1(:)         !(OCEAN)
  REAL(r_size),ALLOCATABLE,SAVE :: dx1(:),dy1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: i1(:),j1(:)     !(OCEAN) splits grid coordinates out into list like ijs

  INTEGER,PARAMETER,PRIVATE :: mpibufsize=1000 !600 !(this worked as 'safe' with 480 procs on Gaea) !200 !1000  
                                        !STEVE: this fixes the problem of bad output when using many nodes, due to default MPI buffer size
CONTAINS


!-----------------------------------------------------------------------
! Initialize the mpi processing
!-----------------------------------------------------------------------
SUBROUTINE set_common_mpi_oceanmodel

  USE params_model, ONLY: nlon, nlat
  USE vars_model,   ONLY: dx, dy
  USE vars_model,   ONLY: lon, lat, kmt0, phi0
  USE vars_model,   ONLY: lon2d, lat2d !, lev2d !STEVE: I'm assuming I'll need this one day

  IMPLICIT NONE

  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d) != 0 !STEVE: initializing
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d) != 0      !STEVE: initializing
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j,n                                 !(OCEAN)
  INTEGER :: l, ll, mstart, mend, nv0 !(OCEAN)
  LOGICAL :: dodebug = .true.

  WRITE(6,'(A)') 'Hello from set_common_mpi_oceanmodel'
  i = MOD(nlon*nlat,nprocs)
  nij1max = (nlon*nlat - i)/nprocs + 1
  WRITE(6,*) "nij1max = ", nij1max
  WRITE(6,*) "mpibufsize = ", mpibufsize

  if (mpibufsize > nij1max) then
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_fast."
  else
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_safe."
  endif

  ! First, identify the number of gridpoints to use on this processor
  if (myrank < i) then
    nij1 = nij1max
  else
    nij1 = nij1max - 1
  endif
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  do n=1,nprocs
    if (n-1 < i) then
      nij1node(n) = nij1max
    else
      nij1node(n) = nij1max - 1
    endif
  enddo

  ! Next, allocate arrays to hold grid, land/sea mask, etc. information
  if (dodebug) WRITE(6,*) "ALLOCATING fields to convert to vectorized form..."
  ALLOCATE(phi1(nij1))
  ALLOCATE(kmt1(nij1))               !(OCEAN)
  ALLOCATE(dx1(nij1))
  ALLOCATE(dy1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))
  ALLOCATE(i1(nij1))                 !(OCEAN)
  ALLOCATE(j1(nij1))                 !(OCEAN)

  nv0=8
  ALLOCATE(v2d(nij1,nv0))
  ALLOCATE(v2dg(nlon,nlat,nv0))

!(UPDATE) with lon2d/lat2d search for observation (LON2D/LAT2D)
  ! Distribute that data to the appropriate processors
  if (dodebug) WRITE(6,*) "Converting dx, dy, lon, lat, i, j, phi0, and kmt0 to vectorized form..."
  v2dg=0.0
  do j=1,nlat
    ! 2D and 1D Data stored in first layer:
    v2dg(:,j,1) = SNGL(dx(:,j))
    v2dg(:,j,2) = SNGL(dy(:,j))
!   v2dg(:,j,3) = SNGL(lon(:))
!   v2dg(:,j,4) = SNGL(lat(j)) !(Single value promoted to array)
    ! 2D Data stored in second layer:
    v2dg(:,j,3) = SNGL(lon2d(:,j))
    v2dg(:,j,4) = SNGL(lat2d(:,j))
    v2dg(:,j,5) = SNGL(kmt0(:,j))
    v2dg(:,j,6) = SNGL(phi0(:,j))  !STEVE: WARNING - this needs to be generalized to the specified dimensions
    !STEVE: For custom localization: (need to know how the grid points are distributed per node)
    do i=1,nlon                             !(OCEAN)
      v2dg(i,j,7) = REAL(i,r_sngl)          !(OCEAN)
    enddo                                   !(OCEAN)
    v2dg(:,j,8) = REAL(j,r_sngl)            !(OCEAN)
    !v2dg(:,j,9) = SNGL(lev2d(:,j))
  enddo
  if (dodebug) WRITE(6,*) "Calling scatter_grd_mpi_small..."
  CALL scatter_grd_mpi_small(0,v2dg,v2d,nlon,nlat,nv0)

  dx1(:)  = v2d(:,1)
  dy1(:)  = v2d(:,2)
  lon1(:) = v2d(:,3)
  lat1(:) = v2d(:,4)
  kmt1(:) = v2d(:,5)                 !(OCEAN)
  phi1(:) = v2d(:,6)                 !(OCEAN)
  i1(:)   = v2d(:,7)                 !(OCEAN)
  j1(:)   = v2d(:,8)                 !(OCEAN)
  !lev1(:) = v2d(:,9)

  DEALLOCATE(v2d,v2dg)
  if (dodebug) WRITE(6,*) "Finished set_common_mpi_oceanmodel..."

END SUBROUTINE set_common_mpi_oceanmodel


!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  LOGICAL :: dodebug = .true.

  if (mpibufsize > nij1max) then
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_fast. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  else
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_safe. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  endif
  if (dodebug) WRITE(6,*) "Finished scatter_grd_mpi..."

END SUBROUTINE scatter_grd_mpi


!-----------------------------------------------------------------------
! 'Safe' version of routine to scatter grid points to all processors
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: tmp(:,:) !(nij1max,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:) !(mpibufsize,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:) !(mpibufsize)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ALLOCATE(tmp(nij1max,nprocs),bufs(mpibufsize,nprocs),bufr(mpibufsize))

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  do n=1,nv3d
    do k=1,nlev
      if (myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
      do iter=1,niter
        if (myrank == nrank) then
          i = mpibufsize * (iter-1)
          do j=1,mpibufsize
            i=i+1
            if (i > nij1) EXIT
            bufs(j,:) = tmp(i,:)
          enddo
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                       & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        i = mpibufsize * (iter-1)
        do j=1,mpibufsize
          i=i+1
          if (i > nij1) EXIT
          v3d(i,k,n) = REAL(bufr(j),r_size)
        enddo
      enddo
    enddo
  enddo

  do n=1,nv2d
    if (myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
    do iter=1,niter
      if (myrank == nrank) then
        i = mpibufsize * (iter-1)
        do j=1,mpibufsize
          i=i+1
          if (i > nij1) EXIT
          bufs(j,:) = tmp(i,:)
        enddo
      endif
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      i = mpibufsize * (iter-1)
      do j=1,mpibufsize
        i=i+1
        if (i > nij1) EXIT
        v2d(i,n) = REAL(bufr(j),r_size)
      enddo
    enddo
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

END SUBROUTINE scatter_grd_mpi_safe


!-----------------------------------------------------------------------
! 'Fast' version of routine to scatter grid points to all processors
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: nlevall

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:,:) !(nij1max,nlevall,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nlevall,nprocs), bufr(nij1max,nlevall))

  ns = nij1max * nlevall
  nr = ns
  if (myrank == nrank) then
    j=0
    do n=1,nv3d
      do k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      enddo
    enddo

    do n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  j=0
  do n=1,nv3d
    do k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    enddo
  enddo

  do n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

END SUBROUTINE scatter_grd_mpi_fast


!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  LOGICAL :: dodebug = .true.

  if (mpibufsize > nij1max) then
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_fast. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  else
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_safe. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  endif

END SUBROUTINE gather_grd_mpi


!-----------------------------------------------------------------------
! 'Safe' version of routine to gather grid points from all processors
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: tmp(:,:) !(nij1max,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufs(:) !(mpibufsize)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(mpibufsize,nprocs)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ALLOCATE(tmp(nij1max,nprocs),bufs(mpibufsize),bufr(mpibufsize,nprocs))

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  do n=1,nv3d
    do k=1,nlev
      do iter=1,niter
        i = mpibufsize * (iter-1)
        do j=1,mpibufsize
          i=i+1
          if (i > nij1) EXIT
          bufs(j) = REAL(v3d(i,k,n),r_sngl)
        enddo
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                      & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        if (myrank == nrank) then
          i = mpibufsize * (iter-1)
          do j=1,mpibufsize
            i=i+1
            if (i > nij1) EXIT
            tmp(i,:) = bufr(j,:)
          enddo
        endif
      enddo
      if (myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
    enddo
  enddo

  do n=1,nv2d
    do iter=1,niter
      i = mpibufsize * (iter-1)
      do j=1,mpibufsize
        i=i+1
        if (i > nij1) EXIT
        bufs(j) = REAL(v2d(i,n),r_sngl)
      enddo
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      if (myrank == nrank) then
        i = mpibufsize * (iter-1)
        do j=1,mpibufsize
          i=i+1
          if (i > nij1) EXIT
          tmp(i,:) = bufr(j,:)
        enddo
      endif
    enddo
    if (myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

END SUBROUTINE gather_grd_mpi_safe


!-----------------------------------------------------------------------
! 'Fast' version of routine to gather grid points from all processors
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: nlevall

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:) !(nij1max,nlevall)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:,:) !(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nlevall),bufr(nij1max,nlevall,nprocs))

  ns = nij1max * nlevall
  nr = ns
  j=0
  do n=1,nv3d
    do k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    enddo
  enddo

  do n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  if (myrank == nrank) then
    j=0
    do n=1,nv3d
      do k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      enddo
    enddo

    do n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

END SUBROUTINE gather_grd_mpi_fast

!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ens_mpi(file,member,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im,mstart,mend
  CHARACTER(11) :: filename='file000'
  !STEVE: debug
  LOGICAL :: dodebug = .true.

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  do l=1,ll
    im = myrank+1 + (l-1)*nprocs
    if (im <= member) then
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'In common_mpi_mom4.f90::read_ens_mpi, MYRANK ',myrank,' is reading a file ',filename
      CALL read_restart(filename,v3dg,v2dg,1) !STEVE: 20150317, trying this out...
!     CALL read_grd4(filename,v3dg,v2dg) !STEVE: 20130709, trying this out...
!     CALL read_grd(filename,v3dg,v2dg)  !STEVE: causes type problem in scatter_grd_mpi
!     if (.false.) CALL write_grd4('test.'//filename,v3dg,v2dg) 
    endif

    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    if (dodebug) WRITE(6,*) "In common_mpi_mom4.f90::read_ens_mpi, calling scatter_grd_mpi_alltoall..."
    CALL scatter_grd_mpi_alltoall(mstart,mend,member,v3dg,v2dg,v3d,v2d)

  enddo

  DEALLOCATE(v3dg,v2dg)

END SUBROUTINE read_ens_mpi

!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,member,v3d,v2d)

  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d

  IMPLICIT NONE

! INCLUDE 'netcdf.inc' !STEVE: for NaN correction (OCEAN)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000'
  INTEGER :: i,j,k,m !STEVE: for debugging
  LOGICAL :: verbose = .true.
  INTEGER :: convcnt = 0
  INTEGER :: mstart,mend

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  do l=1,ll
    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)

    im = myrank+1 + (l-1)*nprocs
    if (im <= member) then
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing file: ',filename

      !STEVE: debug
!     WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))
!     WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_s))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_s)))

      CALL write_restart(filename,v3dg,v2dg)
!     CALL write_grd4(filename,v3dg,v2dg)
    endif
  enddo

  DEALLOCATE(v3dg,v2dg)

END SUBROUTINE write_ens_mpi

!!!!!!!!!
!STEVE: for debugging grid:
!!!!!!!!!
SUBROUTINE write_ens_mpi_grd(file,member,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000.grd'
  INTEGER :: i,j,k,m !STEVE: for debugging
  LOGICAL :: verbose = .false.
  INTEGER :: convcnt = 0

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  do l=1,ll
    do n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      if (im <= member) then
        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
      endif
    enddo

    im = myrank+1 + (l-1)*nprocs
    if (im <= member) then
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename

      !STEVE: debug
!     print *, "common_mpi_mom4.f90::write_ens_mpi_grd:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))

!     CALL write_grd4(filename,v3dg,v2dg)
      CALL write_grd(filename,v3dg,v2dg)
    endif

  enddo

  DEALLOCATE(v3dg,v2dg)

END SUBROUTINE write_ens_mpi_grd

!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  USE params_model, ONLY: nlon, nlat
  IMPLICIT NONE
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  do m=1,nprocs
    do i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    enddo
  enddo

END SUBROUTINE grd_to_buf


!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  USE params_model, ONLY: nlon, nlat
  IMPLICIT NONE
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  do m=1,nprocs
    do i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    enddo
  enddo

END SUBROUTINE buf_to_grd


!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size), ALLOCATABLE :: v3dm(:,:,:) !(nij1,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2dm(:,:) !(nij1,nv2d)
  REAL(r_size), ALLOCATABLE :: v3ds(:,:,:) !(nij1,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2ds(:,:) !(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: i,k,m,n,j,l,ll,im,mstart,mend
  CHARACTER(11) :: filename='file000.grd'

  ALLOCATE(v3dm(nij1,nlev,nv3d),v2dm(nij1,nv2d))
  ALLOCATE(v3ds(nij1,nlev,nv3d),v2ds(nij1,nv2d))
  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)

! ll = CEILING(REAL(member)/REAL(nprocs))
! do l=1,ll
!   mstart = 1 + (l-1)*nprocs
!   mend = MIN(l*nprocs, member)
!   CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
! enddo
  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)

  if (myrank == 0) then
    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    WRITE(6,'(A,I3.3,2A)') 'write_ensmspr_mpi::MYRANK ',myrank,' is writing a file ',filename
!   CALL write_grd4(filename,v3dg,v2dg)
    CALL write_grd(filename,v3dg,v2dg)
  endif

  do n=1,nv3d
    do k=1,nlev
      do i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        do m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        enddo
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      enddo
    enddo
  enddo

  do n=1,nv2d
    do i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      do m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      enddo
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    enddo
  enddo

! do l=1,ll
!   mstart = 1 + (l-1)*nprocs
!   mend = MIN(l*nprocs, member)
!   CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
! enddo
  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)

  if (myrank == 0) then
    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
!   CALL write_grd4(filename,v3dg,v2dg)
    CALL write_grd(filename,v3dg,v2dg)
  endif

  DEALLOCATE(v3dm,v2dm,v3ds,v2ds,v3dg,v2dg)

END SUBROUTINE write_ensmspr_mpi


!-----------------------------------------------------------------------
! Scatter gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_alltoall(mstart,mend,member,v3dg,v2dg,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)
  LOGICAL :: dodebug = .false.

  mcount = mend - mstart + 1
  if (mcount > nprocs .OR. mcount <= 0) STOP

  if (dodebug) WRITE(6,*) "ALLOCATE bufs3,bufr3,bufs2,bufr2..."
  ALLOCATE( bufs3(nij1max,nlev,nprocs) , bufr3(nij1max,nlev,mcount) )
  ALLOCATE( bufs2(nij1max,nprocs)    , bufr2(nij1max,mcount) )

  if (dodebug) WRITE(6,*) "Cycling n=1,nv3d..."
  do n=1,nv3d
    if (dodebug) WRITE(6,*) "n = ",n
    if (myrank < mcount) then
      do k=1,nlev
        if (dodebug) WRITE(6,*) "grd_to_buf :: k = ", k
        CALL grd_to_buf(v3dg(:,:,k,n),bufs3(:,k,:)) !,nlon,nlat,nij1max,nij1node)
      enddo
    endif

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (mcount == nprocs) then
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs3, nij1max*nlev, MPI_REAL, &
                        bufr3, nij1max*nlev, MPI_REAL, MPI_COMM_WORLD, ierr)
    else
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max*nlev,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs3, ns, nst, MPI_REAL, &
                         bufr3, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    endif

    do m = mstart,mend
      do k=1,nlev
         if (dodebug) WRITE(6,*) "m,k,n = ",m,k,n
         if (dodebug) WRITE(6,*) "Assigning v3d(:,k,m,n) = REAL(bufr3(1:nij1,k,m-mstart+1),r_size)"
         if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
         v3d(:,k,m,n) = REAL(bufr3(1:nij1,k,m-mstart+1),r_size)
      END  DO
    enddo
  enddo

  if (dodebug) WRITE(6,*) "Cycling n=1,nv2d..."
  do n=1,nv2d
    if (myrank < mcount) then
      if (dodebug) WRITE(6,*) "grd_to_buf :: 2D"
      CALL grd_to_buf(v2dg(:,:,n),bufs2(:,:)) !,nlon,nlat,nij1max,nij1node)
    endif

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (mcount == nprocs) then
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    else
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    endif

    do m = mstart,mend
      if (dodebug) WRITE(6,*) "m,n = ", m,n
      if (dodebug) WRITE(6,*) "Assigning v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)"
      if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
      v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)
    enddo
  enddo

  if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr3,bufs3,bufr2,bufs2)

END SUBROUTINE scatter_grd_mpi_alltoall


!-----------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  mcount = mend - mstart + 1
  if (mcount > nprocs .OR. mcount <= 0) STOP

  ALLOCATE( bufs3(nij1max,nlev,mcount) , bufr3(nij1max,nlev,nprocs) )
  ALLOCATE( bufs2(nij1max,mcount)    , bufr2(nij1max,nprocs) )

  do n=1,nv3d
    do m = mstart,mend
      do k=1,nlev
        bufs3(1:nij1,k,m-mstart+1) = REAL(v3d(:,k,m,n),r_sngl)
      enddo
    enddo

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (mcount == nprocs) then
      CALL MPI_ALLTOALL(bufs3, nij1max*nlev, MPI_REAL, &
                        bufr3, nij1max*nlev, MPI_REAL, MPI_COMM_WORLD, ierr)
    else
      CALL set_alltoallv_counts(mcount,nij1max*nlev,ns,nst,nr,nrt)
      CALL MPI_ALLTOALLV(bufs3, ns, nst, MPI_REAL, &
                         bufr3, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    endif

    if (myrank < mcount) then
      do k=1,nlev
        CALL buf_to_grd(bufr3(:,k,:),v3dg(:,:,k,n))
      enddo
    endif
  enddo

  do n=1,nv2d
    do m = mstart,mend
      do k=1,nlev
        bufs2(1:nij1,m-mstart+1) = REAL(v2d(:,m,n),r_sngl)
      enddo
    enddo

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (mcount == nprocs) then
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    else
      CALL set_alltoallv_counts(mcount,nij1max,ns,nst,nr,nrt)
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    endif

    if (myrank < mcount) then
      CALL buf_to_grd(bufr2(:,:),v2dg(:,:,n))
    endif
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr3,bufs3,bufr2,bufs2)

END SUBROUTINE gather_grd_mpi_alltoall


!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
SUBROUTINE set_alltoallv_counts(mcount,ngpblock,n_ens,nt_ens,n_mem,nt_mem)
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: mcount,ngpblock
  INTEGER,INTENT(OUT) :: n_ens(nprocs),nt_ens(nprocs),n_mem(nprocs),nt_mem(nprocs)
  INTEGER :: p

  n_ens = 0
  nt_ens = 0
  n_mem = 0
  nt_mem = 0
  do p=1,mcount
    n_ens(p) = ngpblock
    if (myrank+1 == p) then
      n_mem(:) = ngpblock
    endif
  enddo
  do p=2,nprocs
    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
  enddo

END SUBROUTINE set_alltoallv_counts


!-----------------------------------------------------------------------
! Scatter a smaller sized gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
!STEVE: not using at the moment, didn't lead to better performance.
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_smalltoall(mstart,mend,member,v2dg,v2d,nx,ny,nv)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: mstart,mend,member
  INTEGER,INTENT(IN) :: nx,ny,nv
  REAL(r_sngl),INTENT(IN) :: v2dg(nx,ny,nv)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)
  LOGICAL :: dodebug = .true.

  mcount = mend - mstart + 1
  if (mcount > nprocs .OR. mcount <= 0) STOP

  if (dodebug) WRITE(6,*) "ALLOCATE bufs2,bufr2..."
  ALLOCATE( bufs2(nij1max,nprocs)    , bufr2(nij1max,mcount) )

  if (dodebug) WRITE(6,*) "Cycling n=1,nv..."
  do n=1,nv
    if (myrank < mcount) then
      if (dodebug) WRITE(6,*) "grd_to_buf :: 2D"
      CALL grd_to_buf(v2dg(:,:,n),bufs2(:,:)) !,nlon,nlat,nij1max,nij1node)
    endif

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (mcount == nprocs) then
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    else
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    endif

    do m = mstart,mend
      if (dodebug) WRITE(6,*) "m,n = ", m,n
      if (dodebug) WRITE(6,*) "Assigning v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)"
      if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
      v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)
    enddo
  enddo

  if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr2,bufs2)
  
END SUBROUTINE scatter_grd_mpi_smalltoall


!-----------------------------------------------------------------------
! Scatter grid but use only the 2D arrays for improved performance
! (since in the initialization, that is all we need)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi_small(nrank,v2dg,v2d,nx,ny,nv)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nrank
  INTEGER,INTENT(IN) :: nx,ny,nv
  REAL(r_sngl),INTENT(IN) :: v2dg(nx,ny,nv)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:,:) !(nij1max,nlevall,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nv,nprocs), bufr(nij1max,nv))

  ns = nij1max * nv
  nr = ns
  if (myrank == nrank) then
    j=0
    do n=1,nv
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  j=0
  do n=1,nv
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  enddo

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

END SUBROUTINE scatter_grd_mpi_small

END MODULE common_mpi_oceanmodel
