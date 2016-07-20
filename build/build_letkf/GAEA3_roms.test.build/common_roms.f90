MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for ROMS
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   02/02/2009 Takemasa Miyoshi  modified for ROMS
!   06/09/2016 Steve Penny modified for merge with Ocean-LETKF package
!
!=======================================================================
  USE common
  IMPLICIT NONE
  PUBLIC

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_oceanmodel
  USE params_model, ONLY: gridfile, initialize_params_model
  USE params_model, ONLY: nlev
  USE vars_model,   ONLY: phi0, initialize_vars_model
  USE vars_model,   ONLY: lon2d, lat2d, fcori2d
  USE vars_model,   ONLY: mskrho, msku, mskv, mskpsi
  USE vars_model,   ONLY: kmt, kmt0

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid
  INTEGER :: i,j,k

  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'
  CALL initialize_params_model ! (checks to make sure it is initialized)
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  !
  ! Lon, Lat, f, orography
  !
  WRITE(6,'(A)') '  >> accessing file: ', gridfile
  istat = NF_OPEN(gridfile,NF_NOWRITE,ncid)
  if (istat /= NF_NOERR) then
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  endif

  !STEVE: mostly only the lon2d and lat2d fields are used. 
  !       The mask_rho is used to create the land/sea mask
  istat = NF_INQ_VARID(ncid,'lon_rho',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lon2d)
  istat = NF_INQ_VARID(ncid,'lat_rho',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lat2d)
  istat = NF_INQ_VARID(ncid,'f',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,fcori2d) !STEVE: at the moment, this isn't used anywhere
  istat = NF_INQ_VARID(ncid,'h',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,phi0)
  istat = NF_INQ_VARID(ncid,'mask_rho',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,mskrho)
  istat = NF_INQ_VARID(ncid,'mask_u',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,msku)
  istat = NF_INQ_VARID(ncid,'mask_v',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,mskv)
  istat = NF_INQ_VARID(ncid,'mask_psi',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,mskpsi)
  istat = NF_CLOSE(ncid)

  !STEVE: need 'kmt' field giving 'wet' versus 'dry' points
  kmt = 0
  WHERE(mskrho .gt. 0) kmt = nlev
  kmt0 = REAL(kmt,r_size) ! scattered to multiple processes in common_mpi_roms

END SUBROUTINE set_common_oceanmodel


!-----------------------------------------------------------------------
! File I/O (netCDF)
!-----------------------------------------------------------------------
SUBROUTINE read_diag(infile,v3d,v2d)
!===============================================================================
! Read in a set of (hourly or daily) diagnostic files (for now, we're using the restarts)
!===============================================================================
  USE netcdf
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3d0(:,:,:,:)
  REAL(r_sngl), ALLOCATABLE :: v2d0(:,:,:)
! REAL(r_size), ALLOCATABLE :: v3d0(:,:,:,:)
! REAL(r_size), ALLOCATABLE :: v2d0(:,:,:)

  ALLOCATE(v3d0(nlon,nlat,nlev,nv3d), v2d0(nlon,nlat,nv2d))
  CALL read_restart(infile,v3d0,v2d0) !,1)
  v3d = REAL(v3d0,r_size)
  v2d = REAL(v2d0,r_size)

END SUBROUTINE read_diag


!-----------------------------------------------------------------------
! Read a model restart file
!-----------------------------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: iv3d_t, iv3d_s, iv3d_u, iv3d_v
  USE params_model, ONLY: iv2d_z, iv2d_ubar, iv2d_vbar, iv2d_hbl
  USE params_model, ONLY: nv3d, nv2d

  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  CALL read_restart(filename,v3d,v2d)

END SUBROUTINE read_grd


SUBROUTINE read_restart(filename,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: iv3d_t, iv3d_s, iv3d_u, iv3d_v
  USE params_model, ONLY: iv2d_z, iv2d_ubar, iv2d_vbar, iv2d_hbl
  USE params_model, ONLY: nv3d, nv2d
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,:,iv2d_z))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(1:nlon-1,:,iv2d_ubar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (ubar)'
    STOP
  END IF
  v2d(nlon,:,iv2d_ubar) = 0.0
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,1:nlat-1,iv2d_vbar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (vbar)'
    STOP
  END IF
  v2d(:,nlat,iv2d_vbar) = 0.0
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(1:nlon-1,:,:,iv3d_u))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (u)'
    STOP
  END IF
  v3d(nlon,:,:,iv3d_u) = 0.0
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,1:nlat-1,:,iv3d_v))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (v)'
    STOP
  END IF
  v3d(:,nlat,:,iv3d_v) = 0.0
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (salt)'
    STOP
  END IF
  !!! Hsbl
  istat = NF_INQ_VARID(ncid,'Hsbl',varid)
  istat = NF_GET_VAR_REAL(ncid,varid,v2d(:,:,iv2d_hbl))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (Hsbl)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

END SUBROUTINE read_restart


!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: iv3d_t, iv3d_s, iv3d_u, iv3d_v
  USE params_model, ONLY: iv2d_z, iv2d_ubar, iv2d_vbar, iv2d_hbl
  USE params_model, ONLY: nv3d, nv2d
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  CALL write_restart(filename,v3d,v2d)

END SUBROUTINE write_grd



SUBROUTINE write_restart(filename,v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: iv3d_t, iv3d_s, iv3d_u, iv3d_v
  USE params_model, ONLY: iv2d_z, iv2d_ubar, iv2d_vbar, iv2d_hbl
  USE params_model, ONLY: nv3d, nv2d
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: ncid,istat,varid
  istat = NF_OPEN(filename,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR -LPP'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'zeta',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,:,iv2d_z))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'ubar',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(1:nlon-1,:,iv2d_ubar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (ubar)'
    STOP
  END IF
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vbar',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,1:nlat-1,iv2d_vbar))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (vbar)'
    STOP
  END IF
  !!! u
  istat = NF_INQ_VARID(ncid,'u',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(1:nlon-1,:,:,iv3d_u))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (u)'
    STOP
  END IF
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,1:nlat-1,:,iv3d_v))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (v)'
    STOP
  END IF
  !!! t
  istat = NF_INQ_VARID(ncid,'temp',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'salt',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (salt)'
    STOP
  END IF
  !!! Hsbl
  istat = NF_INQ_VARID(ncid,'Hsbl',varid)
  istat = NF_PUT_VAR_REAL(ncid,varid,v2d(:,:,iv2d_hbl))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (Hsbl)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

END SUBROUTINE write_restart

!-----------------------------------------------------------------------
! Depth
!-----------------------------------------------------------------------
SUBROUTINE calc_depth(zeta,bottom,depth)
  USE params_model, ONLY: nlev
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: zeta
  REAL(r_size),INTENT(IN) :: bottom
  REAL(r_size),INTENT(OUT) :: depth(nlev)
  REAL(r_size) :: cs(nlev)
  REAL(r_size) :: surfd
  INTEGER      :: k
  INTEGER      :: KI
  REAL         :: K1
  REAL(r_size),SAVE :: SC_R(nlev)
  REAL(r_size),SAVE :: SC_W(nlev)
  REAL(r_size),SAVE :: CFF1,CFF2
!  REAL(r_size),PARAMETER :: THETA_B=5.    ! LPP - Same parameters used to generate original ROMS Sigma Coord
!  REAL(r_size),PARAMETER :: THETA_S=0.6   ! LPP - Same parameters used to generate original ROMS Sigma Coord
   REAL(r_size),PARAMETER :: THETA_B=0.6   ! Leo
   REAL(r_size),PARAMETER :: THETA_S=5.    ! Leo

      IF(THETA_S.NE.0.D0) THEN
       CFF1=1.D0/DSINH(THETA_S)
       CFF2=0.5D0/DTANH(0.5D0*THETA_S)
      ELSE
       CFF1=0.D0
       CFF2=0.D0
      ENDIF

      DO K=nlev,1,-1
      KI=(nlev-K)+1
      K1=REAL(-K+1.D0)
      SC_R(KI)=(K1-0.5D0)/REAL(nlev)
      ENDDO

      DO K=1,nlev
      cs(K)=(1.D0-THETA_B)*CFF1*DSINH(THETA_S*SC_R(K))+ &
     & THETA_B*(CFF2*DTANH(THETA_S*(SC_R(K)+0.5D0))-0.5D0)
      IF(THETA_S.EQ.0.D0) cs(K)=SC_R(K)
      ENDDO
      
      DO k=1,nlev
      depth(k) = zeta+(zeta+bottom)* &
     & (10.0d0*(REAL(k,r_size)-REAL(nlev,r_size)-0.5d0)/ &
     & REAL(nlev,r_size)+cs(k)*bottom)/(bottom+10.0d0)
      END DO
      surfd = depth(nlev)
      DO k=1,nlev
        depth(k) = depth(k) - surfd
      END DO

END SUBROUTINE calc_depth

!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

! DO k=1,nlev
!   WRITE(6,'(I2,A)') k,'th level'
!   DO n=1,nv3d
!     WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
!   END DO
! END DO

! DO n=1,nv2d
!   WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
! END DO

END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  USE params_model, ONLY: nlev, nv3d, nv2d
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

END SUBROUTINE ensmean_grd

END MODULE common_oceanmodel
