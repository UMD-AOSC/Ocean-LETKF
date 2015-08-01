MODULE common_mpi_mom4
!===============================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!   04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!
!===============================================================================
  USE common
  USE common_mpi
  USE common_mom4

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: mpibufsize=1000 !600 !(this worked as 'safe' with 480 procs on Gaea) !200 !1000  !STEVE: this fixes the problem of bad output when using over 6 nodes default=1000,mom2(mpich2)=200
  INTEGER,SAVE :: nij1                  !STEVE: this is the number of gridpoints to run on this (myrank) processor
  INTEGER,SAVE :: nij1max               !STEVE: the largest number of gridpoints on any 1 processor
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: kmt1(:)         !(OCEAN)
  REAL(r_size),ALLOCATABLE,SAVE :: dx1(:),dy1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: i1(:),j1(:)     !(OCEAN) splits grid coordinates out into list like ijs

CONTAINS

SUBROUTINE set_common_mpi_mom4
!===============================================================================
! Initialize this module
!===============================================================================
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d) != 0 !STEVE: initializing
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d) != 0      !STEVE: initializing
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j,n                                 !(OCEAN)

  WRITE(6,'(A)') 'Hello from set_common_mpi_mom4'
  i = MOD(nlon*nlat,nprocs)
  nij1max = (nlon*nlat - i)/nprocs + 1
  WRITE(6,*) "nij1max = ", nij1max
  WRITE(6,*) "mpibufsize = ", mpibufsize

  IF(mpibufsize > nij1max) THEN
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_fast."
  ELSE
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_safe."
  END IF

  IF(myrank < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO

  ALLOCATE(phi1(nij1))
  ALLOCATE(kmt1(nij1))               !(OCEAN)
  ALLOCATE(dx1(nij1))
  ALLOCATE(dy1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))
  ALLOCATE(i1(nij1))                 !(OCEAN)
  ALLOCATE(j1(nij1))                 !(OCEAN)

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))
  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2dg(nlon,nlat,nv2d))
  DO j=1,nlat
    v3dg(:,j,1,1) = SNGL(dx(:,j))
    v3dg(:,j,1,2) = SNGL(dy(:,j))
    v3dg(:,j,1,3) = SNGL(lon(:))
    v3dg(:,j,1,4) = SNGL(lat(j)) !(Single value promoted to array)
    !STEVE: For custom localization: (need to know how the grid points are distributed per node)
    DO i=1,nlon                      !(OCEAN)
      v3dg(i,j,2,3) = REAL(i,r_sngl)        !(OCEAN)
    ENDDO                            !(OCEAN)
    v3dg(:,j,2,4) = REAL(j,r_sngl)          !(OCEAN)
  END DO
  v2dg(:,:,1) = SNGL(phi0(:,:))
  v2dg(:,:,2) = SNGL(kmt0(:,:))
  CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)
  dx1(:) = v3d(:,1,1)
  dy1(:) = v3d(:,1,2)
  lon1(:) = v3d(:,1,3)
  lat1(:) = v3d(:,1,4)
  i1(:) = v3d(:,2,3)                 !(OCEAN)
  j1(:) = v3d(:,2,4)                 !(OCEAN)
  phi1 = v2d(:,1)
  kmt1 = v2d(:,2)                    !(OCEAN)

  DEALLOCATE(v3d,v2d,v3dg,v2dg)

  RETURN
END SUBROUTINE set_common_mpi_mom4

SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
!===============================================================================
! Scatter gridded data to processes (nrank -> all)
!===============================================================================
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  LOGICAL :: dodebug = .true.

  IF(mpibufsize > nij1max) THEN
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_fast. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  ELSE
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_safe. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  END IF

  RETURN
END SUBROUTINE scatter_grd_mpi

SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!===============================================================================
! Scatter using a smaller buffer
!===============================================================================
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

  DO n=1,nv3d
    DO k=1,nlev
      IF(myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
      DO iter=1,niter
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1) EXIT
            bufs(j,:) = tmp(i,:)
          END DO
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                       & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          v3d(i,k,n) = REAL(bufr(j),r_size)
        END DO
      END DO
    END DO
  END DO

  DO n=1,nv2d
    IF(myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
    DO iter=1,niter
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          bufs(j,:) = tmp(i,:)
        END DO
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        v2d(i,n) = REAL(bufr(j),r_size)
      END DO
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_safe

SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!===============================================================================
! Scatter using a larger buffer
!===============================================================================
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
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_fast

SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
!===============================================================================
! Gather gridded data (all -> nrank)
!===============================================================================
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  LOGICAL :: dodebug = .true.

  IF(mpibufsize > nij1max) THEN
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_fast. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  ELSE
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_safe. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  END IF

  RETURN
END SUBROUTINE gather_grd_mpi

SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
!===============================================================================
! Gather the gridded data using a small buffer
!===============================================================================
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

  DO n=1,nv3d
    DO k=1,nlev
      DO iter=1,niter
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          bufs(j) = REAL(v3d(i,k,n),r_sngl)
        END DO
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                      & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1) EXIT
            tmp(i,:) = bufr(j,:)
          END DO
        END IF
      END DO
      IF(myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    DO iter=1,niter
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        bufs(j) = REAL(v2d(i,n),r_sngl)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          tmp(i,:) = bufr(j,:)
        END DO
      END IF
    END DO
    IF(myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

  RETURN
END SUBROUTINE gather_grd_mpi_safe

SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
!===============================================================================
! Gather the gridded data using a large buffer
!===============================================================================
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
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

  RETURN
END SUBROUTINE gather_grd_mpi_fast

SUBROUTINE read_ens_mpi(file,member,v3d,v2d)
!===============================================================================
! Read ensemble data and distribute to processes
!===============================================================================
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000'
  !STEVE: debug
  LOGICAL :: dodebug = .true.
! REAL(r_sngl) :: v3dgT(nlon,nlat,nlev,nv3d)
! REAL(r_sngl) :: v2dgT(nlon,nlat,nv2d)

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'In common_mpi_mom4.f90::read_ens_mpi, MYRANK ',myrank,' is reading a file ',filename
      CALL read_grd4(filename,v3dg,v2dg) !STEVE: 20130709, trying this out...
!     CALL read_grd(filename,v3dg,v2dg)  !STEVE: causes type problem in scatter_grd_mpi
!     if (.false.) CALL write_grd4('test.'//filename,v3dg,v2dg) 
    END IF

    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        WRITE(6,*) "In common_mpi_mom4.f90::read_ens_mpi, calling scatter_grd_mpi..."
        CALL scatter_grd_mpi(n,v3dg,v2dg,v3d(:,:,im,:),v2d(:,im,:))
        WRITE(6,*) "In common_mpi_mom4.f90::read_ens_mpi, finished calling scatter_grd_mpi."
!       if (dodebug) then
!         v3dgT=0
!         v2dgT=0
!         CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dgT,v2dgT)
!         WRITE(filename(1:7),'(A4,I3.3)') file,im
!         WRITE(6,*) "Writing to file: test."//filename 
!         CALL write_grd4('test.'//filename,v3dgT,v2dgT) 
!       endif
      END IF
    END DO
  END DO

  DEALLOCATE(v3dg,v2dg)

  RETURN
END SUBROUTINE read_ens_mpi

SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
!===============================================================================
! Write ensemble data after collecting data from processes
!===============================================================================
  INCLUDE 'netcdf.inc' !STEVE: for NaN correction (OCEAN)
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

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        WRITE(6,*) "In common_mpi_mom4.f90::write_ens_mpi, calling gather_grd_mpi..."
        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
        WRITE(6,*) "In common_mpi_mom4.f90::write_ens_mpi, finished calling gather_grd_mpi."
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing file: ',filename

      !STEVE: debug
      WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))
      WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_s))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_s)))

      CALL write_grd4(filename,v3dg,v2dg)
    END IF
  END DO

  DEALLOCATE(v3dg,v2dg)

  RETURN
END SUBROUTINE write_ens_mpi

SUBROUTINE write_ens_mpi_grd(file,member,v3d,v2d)
!===============================================================================
! For debugging the grid
!===============================================================================
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
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename

      !STEVE: debug
      print *, "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))

      CALL write_bingrd4(filename,v3dg,v2dg)
    END IF

  END DO

  DEALLOCATE(v3dg,v2dg)

  RETURN
END SUBROUTINE write_ens_mpi_grd

SUBROUTINE grd_to_buf(grd,buf)
!===============================================================================
! gridded data -> buffer
!===============================================================================
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  RETURN
END SUBROUTINE grd_to_buf

SUBROUTINE buf_to_grd(buf,grd)
!===============================================================================
! buffer -> gridded data
!===============================================================================
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd

SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
!===============================================================================
! STORING DATA (ensemble mean and spread)
!===============================================================================
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
  INTEGER :: i,k,m,n,j
  CHARACTER(11) :: filename='file000.grd'

  ALLOCATE(v3dm(nij1,nlev,nv3d),v2dm(nij1,nv2d))
  ALLOCATE(v3ds(nij1,nlev,nv3d),v2ds(nij1,nv2d))
  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)

  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    WRITE(6,'(A,I3.3,2A)') 'write_ensmspr_mpi::MYRANK ',myrank,' is writing a file ',filename
    CALL write_bingrd4(filename,v3dg,v2dg)
  END IF

  DO n=1,nv3d
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      END DO
    END DO
  END DO

  DO n=1,nv2d
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    END DO
  END DO

  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_bingrd4(filename,v3dg,v2dg)
  END IF

  DEALLOCATE(v3dm,v2dm,v3ds,v2ds,v3dg,v2dg)

  RETURN
END SUBROUTINE write_ensmspr_mpi

END MODULE common_mpi_mom4
