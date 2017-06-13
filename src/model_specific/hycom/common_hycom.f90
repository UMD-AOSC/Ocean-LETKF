MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for HYCOM
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   06/08/2015 Steve Penny converted for use with HYCOM
!
!=======================================================================
  USE common
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_SLA
  USE vars_model,   ONLY: lon, lat, lev, lon2d, lat2d, lon0, lonf, lat0, latf, wrapgap, dx, dy, kmt, kmt0, phi0
  USE vars_model,   ONLY: pmsk,umsk,vmsk 
  USE hycom_io,     ONLY: read_blkdat
  USE hycom_io,     ONLY: get_hycom_depth,get_hycom_sample 
  USE params_model, ONLY: SSHclm_file
  USE vars_model,   ONLY: SSHclm_m

  IMPLICIT NONE

  PUBLIC :: set_common_oceanmodel, read_diag, read_restart, write_restart
  PUBLIC :: read_grd, write_grd
  PUBLIC :: ensmean_grd

  PRIVATE

CONTAINS

SUBROUTINE set_common_oceanmodel
!===============================================================================
! Initialize the module
!===============================================================================
  USE netcdf
  USE params_model, ONLY: gridfile, gridfile1, gridfile2, gridfile3, ncundef
  USE vars_model,   ONLY: area_t
  USE vars_model,   ONLY: fcori, fcori2d
  !USE vars_model,   ONLY: pmsk, umsk, vmsk 
  USE params_model, ONLY: initialize_params_model
  USE vars_model,   ONLY: initialize_vars_model

  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  LOGICAL :: ex

  ! For temporary dx and dy computations:
  REAL(r_size) :: dlon, dlat, d1, d2, d3

  REAL(r_size),SAVE :: sample_t(nlon,nlat,nlev)      !(OCEAN)(HYCOM)
  REAL(r_sngl),ALLOCATABLE :: sigma(:)               !(HYCOM) reference density levels, will be assigned to lev()

  ! For debugging:
  LOGICAL :: dodebug = .true.

  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'
  CALL initialize_params_model ! (checks to make sure it is initialized)
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  !STEVE: needed to index the observations, because HYCOM does not have a regular grid above 45ºN
  !       Instead, HYCOM has a 2D-field for the representation of the lon/lat coordinate.
  !       This requires an update to the lon and lat arrays to lon2d and lat2d, with subsequent
  !       updates in letkf_local.f90 and letkf_obs.f90. These are primarily used to sort and store
  !       observation locations for later lookup at the localization phase.

  !
  ! Lon, Lat, f, orography
  !
!STEVE: this part adapted from ROMS, update for MOM4 (and later MOM6) netcdf files:
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim("regional.grid.a"),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", "regional.grid.a" 
    WRITE(6,*) "Exiting common_hycom.f90..."
    STOP(2)
  ENDIF
  call get_hycom_depth(phi0,lat2d,lon2d,dx,dy) 
  WRITE(6,'(A)') '  >> accessing file: ', "regional.grid.a" 
  WRITE(6,*) "lon2d(1,1) = ", lon2d(1,1)
  WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)
  WRITE(6,*) "lat2d(1,1) = ", lat2d(1,1)
  WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)

  lat=lat2d(2,:)
  lon=lon2d(:,2)

  !-----------------------------------------------------------------------------
  ! Get the reference density levels:
  !-----------------------------------------------------------------------------
  CALL read_blkdat(gridfile1,sigma)
  lev = sigma
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !-----------------------------------------------------------------------------
  !(HYCOM) Sample temperature to make kmt field (2d land-sea depth mask)
  !-----------------------------------------------------------------------------

  ! JILI read HYCOM sample to get mask, which is needed by HYCOM output 
  CALL get_hycom_sample(sample_t,pmsk,umsk,vmsk)

  WRITE(6,*) "sample_t(1,1,1) = ", sample_t(1,1,1)
  WRITE(6,*) "sample_t(nlon,nlat,nlev) = ", sample_t(nlon,nlat,nlev)

  if (.false. .and. dodebug) then
  WRITE(6,*) "========================================="
  WRITE(6,*) "STEVE: debugging for HYCOM input grid..."
  WRITE(6,*) "lon = ", lon
  WRITE(6,*) "lat = ", lat
  WRITE(6,*) "lon2d(:,1) = ", lon2d(:,1)
  WRITE(6,*) "lon2d(1,:) = ", lon2d(1,:)
  WRITE(6,*) "lat2d(:,1) = ", lat2d(:,1)
  WRITE(6,*) "lat2d(1,:) = ", lat2d(1,:)
  WRITE(6,*) "MAXVAL(lat2d) = ", MAXVAL(lat2d)
  WRITE(6,*) "done lon/lat."
  WRITE(6,*) "========================================="
  endif

  !
  ! dx and dy
  !

  if (.true.) then !STEVE: for HYCOM: this is TEMPORARY, better if it's read in
    WRITE(6,*) "Computing area_t..."
!    dx=0.0
!    dy=0.0
    area_t=dx*dy
  endif

! WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile3
  WRITE(6,*) "dx(1,1) = ", dx(1,1)
  WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
  WRITE(6,*) "dy(1,1) = ", dy(1,1)
  WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)
  WRITE(6,*) "area_t(1,1) = ", area_t(1,1)
  WRITE(6,*) "area_t(nlon,nlat) = ", area_t(nlon,nlat)

  !
  ! kmt data
  !
  !STEVE:MOM6: sum(h>0) (from MOM.res.nc) to find ocean layer depth, where depth > 0 (from ocean_togog.nc)
  !STEVE:HYCOM: this needs to correspond to the reference model layers,
  !             for now, it does not since the layers have been interpolated
  !             from 32 layer to 40 depths in the netcdf files.


  WRITE(6,*) "phi0(1,1) = ", phi0(1,1)
  WRITE(6,*) "phi0(nlon,nlat) = ", phi0(nlon,nlat)

  ! Find the 'kmt' value for the depth of the ocean grid cells (phi0)
  ! sum the layer thicknesses where depth > 0
  kmt=0
!  do k=nlev,1,-1
    do j=1,nlat
      do i=1,nlon
!       if (height(i,j,k) > 0 .and. wet(i,j) > 0 .and. kmt(i,j) == 0 ) then
!        if (kmt(i,j) == 0 .and. sample_t(i,j,k) < ncundef-2 ) then
! JILI change undef to a smaller number
! JILI kmt=0 for land and nlev for ocean
!        if (kmt(i,j) == 0 .and. sample_t(i,j,k) < 100000.0 ) then
!!         phi0(i,j) = sum(height(i,j,1:k))
!          kmt(i,j) = k
         if (phi0(i,j) > 100000.0) then
            kmt(i,j)=0
         else
            kmt(i,j)=nlev
         endif
      enddo
    enddo
!  enddo
  !kmt=nlev
  kmt0 = REAL(kmt,r_size)

  !
  ! Corioris parameter
  !
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)
  do j=1,nlon
    do i=1,nlat
      fcori2d(i,j) = 2.0d0 * r_omega * sin(lat2d(i,j)*pi/180.0d0)
    enddo
  enddo


  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  !STEVE: for MOM4p1:
  wrapgap = 360.0d0 - abs(modulo(lon0,360.0)) - abs(modulo(lonf,360.0))

  if (dodebug) then
    WRITE(6,*) "common_hycom.f90::set_common_hycom:: lon0,lonf,lat0,latf = ", lon0,lonf,lat0,latf
  endif

END SUBROUTINE set_common_oceanmodel


!STEVE: ISSUE: need to update for HYCOM
SUBROUTINE read_etaclm
!===============================================================================
! Read in (e.g.) 1993-1999 model climatology to subtract as done by AVISO
!===============================================================================
  USE netcdf

  IMPLICIT NONE

  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: i,j
  INTEGER :: ncid, varid

  ! read the model SSH climatology netcdf file
  ! read into: SSHclm_m
  CALL check( NF90_OPEN(SSHclm_file,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file

  buf4=0.0
  CALL check( NF90_INQ_VARID(ncid,'ssh',varid) )
  CALL check( NF90_GET_VAR(ncid,varid,buf4) )
  DO j=1,nlat
    DO i=1,nlon
      !STEVE: Hopefully reading in meters here... (data might be in cm)
      SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
    END DO
  END DO

  CALL check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm


!STEVE: add this:
SUBROUTINE read_diag(infile,v3d,v2d,prec_in)
!===============================================================================
! Read a archive/restart file in binary format extracted from hycom ab-format
!===============================================================================
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by HYCOM
  !       The t,s,u,v,(ssh) fields are read in from the ab archive files so that the o-f data can be computed
  !       by the obsop.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  USE netcdf
  USE hycom_io, ONLY: read_hycom, hycom_undef,get_hycom
  USE params_model, ONLY: base,base_a,base_b
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: infile
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:),  INTENT(OUT) :: v2d !(nlon,nlat,nv2d)
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  INTEGER :: prec
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4 !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8 !(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .true.
  CHARACTER(slen) :: varname, binfile,infile_a,infile_b

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2d(nlon,nlat,nv2d))

  ALLOCATE(buf8(nlon,nlat,nlev))

  ! Check for input argument for precision
  if (PRESENT(prec_in)) then
    prec = prec_in
  else
    prec = 1 ! default is 1
  endif

  ! Use specified precision
  if (prec == 1) then
    WRITE(6,*) "read_diag::  single precision requested, but option not available at this time."
    WRITE(6,*) "read_diag::  using double precision read..."
  elseif (prec == 2) then
    WRITE(6,*) "read_diag::  using double precision read..."
  else
    WRITE(6,*) "read_diag:: unsupported option, prec = ", prec
    STOP
  endif

  ! STEVE: this is provided externally at the moment
  binfile = trim(infile)//trim(base)
  infile_a = trim(infile)//trim(base_a)
  infile_b = trim(infile)//trim(base_b)


  if (dodebug) WRITE(6,*) "read_diag:: Pre-read_hycom: Reading file: ", trim(binfile)
! JILI Read HYCOM a and b files
  CALL get_hycom(infile_a,infile_b,v3d,v2d)

  if (dodebug) WRITE(6,*) "read_diag:: Post-read_hycom: Just read file: ", trim(binfile)

  if (dodebug) then
    WRITE(6,*) "read_diag:: Post-read_hycom"
    buf8 = v3d(:,:,:,iv3d_t)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_t)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_t)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_s)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_s)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_s)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_u)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_u)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_u)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_v)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_v)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_v)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_h)
    where (buf8 == hycom_undef) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_h)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_h)) = ", MINVAL(buf8(:,:,1))
    WRITE(6,*)
  endif
  
  DEALLOCATE(buf8)

END SUBROUTINE read_diag


!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec_in)
  !STEVE: This subroutine reads the HYCOM restart file.
  !       The t,s,u,v,(ssh) fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsope.f90 program prior
  !       to running letkf.f90, so the archive (diagnostic) files should not have to be touched during letkf runtime.
  USE netcdf
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h

  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:),  INTENT(OUT) :: v2d !(nlon,nlat,nv2d)
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  INTEGER :: prec
  REAL(r_size), ALLOCATABLE :: v3d8(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d8(:,:,:)

  ! Check for input argument for precision
  if (PRESENT(prec_in)) then
    prec = prec_in
  else
    prec = 2 ! default is 2
  endif

! !STEVE: for now, only use the HYCOM archive format. The output ab file will be converted to a restart file externally.
  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
  ALLOCATE(v3d8(nlon,nlat,nlev,nv3d),v2d8(nlon,nlat,nv2d))
  CALL read_diag(infile,v3d8,v2d8,prec) !WARNING - default for read_diag may be prec=1
  v3d = REAL(v3d8,r_sngl)
  v2d = REAL(v2d8,r_sngl)
  DEALLOCATE(v3d8,v2d8)

END SUBROUTINE read_restart

!-----------------------------------------------------------------------
subroutine check(status)
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check


!-----------------------------------------------------------------------
! Write a set of MOM6 restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d,v2d)
  !STEVE: This writes out the analysis to a pre-existing template netcdf file.
  !       IN THE FUTURE: output directly as a HYCOM restart file
  USE netcdf
  USE hycom_io, ONLY: write_hycom, hycom_undef,put_hycom
  USE vars_model, ONLY:pmsk,umsk,vmsk
  !USE params_model, ONLY: base, base_a,base_b,do_physlimit,nlev,nlat,nlon
  USE params_model

  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  !REAL(r_sngl),DIMENSION(:,:,:,:),INTENT(IN) :: v3d !(nlon,nlat,nlev,nv3d)
  !REAL(r_sngl),DIMENSION(:,:,:),  INTENT(IN) :: v2d !(nlon,nlat,nv2d)
  REAL(r_sngl),DIMENSION(:,:,:,:) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),DIMENSION(:,:,:)   :: v2d !(nlon,nlat,nv2d)

  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  CHARACTER(slen) :: binfile,infile_a,infile_b
  CHARACTER(3) :: MEM3

  ! STEVE: this is provided externally at the moment
  binfile = trim(outfile)//trim(base)
  infile_a = "gs01"//outfile(5:7)//trim(base_a)
  infile_b = "gs01"//outfile(5:7)//trim(base_b)
  
  ! STEVE: for safety, clean up the variables for output:
  ! JILI for land grids, also set variables to undef
  ! JILI a potential probelm: bottom layers large spread-need to be very careful
  if (do_physlimit) then
  do k=1,nlev
    do j=1,nlat
      do i=1,nlon
        if (v3d(i,j,k,iv3d_t) < min_t) then

          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = min_t
        endif


        if (v3d(i,j,k,iv3d_t) > max_t) then
          WRITE(6,*) "WARNING: Bad temp value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
          v3d(i,j,k,iv3d_t) = max_t
        endif


        if (v3d(i,j,k,iv3d_s) < min_s ) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = min_s
        endif


        if (v3d(i,j,k,iv3d_s) > max_s) then
          WRITE(6,*) "WARNING: Bad salt value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
          v3d(i,j,k,iv3d_s) = max_s
        endif

        if (v3d(i,j,k,iv3d_u) < min_uv) then
          WRITE(6,*) "WARNING: Bad u-vel value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
          v3d(i,j,k,iv3d_u) = min_uv
        endif

        if (v3d(i,j,k,iv3d_u) > max_uv) then
          WRITE(6,*) "WARNING: Bad u-vel value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
          v3d(i,j,k,iv3d_u) = max_uv
        endif

        if (v3d(i,j,k,iv3d_v) < min_uv) then
          WRITE(6,*) "WARNING: Bad v-vel value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
          v3d(i,j,k,iv3d_v) = min_uv
        endif

        if (v3d(i,j,k,iv3d_v) > max_uv) then
          WRITE(6,*) "WARNING: Bad v-vel value in analysis output:"
          WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
          v3d(i,j,k,iv3d_v) = max_uv
        endif

        if (phi0(i,j) > 100000) then
          v3d(i,j,k,1:4)=ncundef
        endif

      enddo
    enddo
  enddo

    do j=1,nlat
      do i=1,nlon
        if (phi0(i,j) > 100000.0) then
          v2d(i,j,:)=ncundef
        endif
      end do
    end do

  endif


  WRITE(6,*) "common_hycom.f90::write_restart:: writing file: ", trim(binfile)

  WRITE(6,*) "common_hycom.f90::write_restart:: calling write_hycom..."
  !JILI HYCOM archive file output 
  CALL put_hycom(infile_a,infile_b,v3d,v2d,pmsk,umsk,vmsk)
  WRITE(6,*) "common_hycom.f90::write_restart:: Finished calling write_hycom."

END SUBROUTINE write_restart


!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------

!-----------------------------------------------------------------------
! Read a grid file with single precision
!-----------------------------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
  USE params_model, ONLY: nij0
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:), INTENT(OUT) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:),   INTENT(OUT) :: v2d !(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  iunit=11
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)


END SUBROUTINE read_grd


!-----------------------------------------------------------------------
! Write a grid file in single precision
!-----------------------------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  ! STEVE: binary native format for outputting ensemble mean and spread
  !        (may want to change to HYCOM ab archive format in the future)
  USE params_model, ONLY: nij0
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),DIMENSION(:,:,:,:), INTENT(IN) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),DIMENSION(:,:,:),   INTENT(IN) :: v2d !(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL, PARAMETER :: dodebug=.false.

  if (dodebug) print *, "write_grd4:: open filename = ",filename
  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_grd4:: nij0,iolength = ", nij0,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      if (dodebug) print *, "write_grd4:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    if (dodebug) print *, "write_grd4:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

END SUBROUTINE write_grd


SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
  use hycom_io, ONLY: hycom_undef
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: cnt3d(nij,nlev,nv3d)
  INTEGER :: cnt2d(nij,nv2d)
  INTEGER :: i,k,m,n

! ALLOCATE(v3dm(nij,nlev,nv3d),v2dm(nij,nv2d))


  !STEVE: For HYCOM, we have to check whether the member point
  !       is undefined, because some members can differ in the
  !       grid definition.

  cnt3d = 0
  v3dm  = 0.0d0
  do n=1,nv3d
    do k=1,nlev
      do i=1,nij
        do m=1,member
          if (v3d(i,k,m,n) < hycom_undef) then
            v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
            cnt3d(i,k,n) = cnt3d(i,k,n) + 1
          endif
        enddo
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(cnt3d(i,k,n),r_size)
      enddo
    enddo
  enddo

  cnt2d = 0
  v2dm = 0.0d0
  do n=1,nv2d
    do i=1,nij
      do m=1,member
        if (v2d(i,m,n) < hycom_undef) then
          v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
          cnt2d(i,n) = cnt2d(i,n) + 1
        endif
      enddo
      v2dm(i,n) = v2dm(i,n) / REAL(cnt2d(i,n),r_size)
    enddo
  enddo

END SUBROUTINE ensmean_grd


END MODULE common_oceanmodel

