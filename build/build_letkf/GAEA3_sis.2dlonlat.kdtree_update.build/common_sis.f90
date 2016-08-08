MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for SIS
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   04/06/2015 Steve Penny converted for use with SIS
!
!=======================================================================

  USE common
  IMPLICIT NONE
  PUBLIC

  !STEVE: for filtering undef values from netcdf file
  REAL(r_size), PARAMETER :: vmax = 1.0e18

CONTAINS


!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_oceanmodel
  USE netcdf
  USE params_model
  USE vars_model
! USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_MLD, DO_SLA

  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex

  WRITE(6,'(A)') 'Hello from set_common_sis'
  !
  ! Lon, Lat, f, orography
  !
!STEVE: GOAL: to utilize all netcdf grid data to completely define the grid and all grid-dependent operations
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile
    WRITE(6,*) "Exiting common_sis.f90..."
    STOP(2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile
  call check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )
  call check( NF90_INQ_VARID(ncid,'grid_x_T',varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon) )
  WRITE(6,*) "lon(1) = ", lon(1)
  WRITE(6,*) "lon(nlon) = ", lon(nlon)
  call check( NF90_INQ_VARID(ncid,'grid_y_T',varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat) )
  WRITE(6,*) "lat(1) = ", lat(1)
  WRITE(6,*) "lat(nlat) = ", lat(nlat)
  call check( NF90_INQ_VARID(ncid,'zt',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !
  ! dx and dy
  !
  call check( NF90_INQ_VARID(ncid,'ds_01_21_T',varid) )    ! width of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dx) )
  call check( NF90_INQ_VARID(ncid,'ds_10_12_T',varid) )    ! height of T_cell (meters)
  call check( NF90_GET_VAR(ncid,varid,dy) )
  call check( NF90_INQ_VARID(ncid,'area_T',varid) )        ! area of T_cell
  call check( NF90_GET_VAR(ncid,varid,area_t) )
  WRITE(6,*) "common_sis:: grid_spec.nc MIN(dx) = ", MINVAL(dx)
  WRITE(6,*) "common_sis:: grid_spec.nc MAX(dx) = ", MAXVAL(dx)
  WRITE(6,*) "common_sis:: grid_spec.nc MIN(dy) = ", MINVAL(dy)
  WRITE(6,*) "common_sis:: grid_spec.nc MAX(dy) = ", MAXVAL(dy)
  WRITE(6,*) "common_sis:: grid_spec.nc MIN(area_t) = ", MINVAL(area_t)
  WRITE(6,*) "common_sis:: grid_spec.nc MAX(area_t) = ", MAXVAL(area_t)

  !
  ! kmt data
  !
  call check( NF90_INQ_VARID(ncid,'num_levels',varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,kmt0) )
  WRITE(6,*) "kmt0(1,1) = ", kmt0(1,1)
  WRITE(6,*) "kmt0(nlon,nlat) = ", kmt0(nlon,nlat)
  kmt = NINT(kmt0)
  call check( NF90_INQ_VARID(ncid,'wet',varid) )        ! land/sea flag (0=land) for T-cell
  call check( NF90_GET_VAR(ncid,varid,wet) )
  WRITE(6,*) "wet(1,1) = ", wet(1,1)
  WRITE(6,*) "wet(nlon,nlat) = ", wet(nlon,nlat)

  WRITE(6,*) "Using dx and dy from netcdf file: ", gridfile
  WRITE(6,*) "dx(1,1) = ", dx(1,1)
  WRITE(6,*) "dx(nlon,nlat) = ", dx(nlon,nlat)
  WRITE(6,*) "dy(1,1) = ", dy(1,1)
  WRITE(6,*) "dy(nlon,nlat) = ", dy(nlon,nlat)

  !
  ! Corioris parameter
  !
  fcori(:) = 2.0d0 * r_omega * sin(lat(:)*pi/180.0d0)

  ! Close the grid_spec.nc file:
  call check( NF90_CLOSE(ncid) )

  ! STEVE: for (more) generalized (longitude) grid:
  lon0 = lon(1)
  lonf = lon(nlon)
  lat0 = lat(1)
  latf = lat(nlat)
  wrapgap = 360.0d0 - abs(lon0) - abs(lonf)

END SUBROUTINE set_common_oceanmodel

!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!STEVE: add this:
!-- Read a diagnostic file in sis netcdf format ---------------------------------------------------
SUBROUTINE read_diag(infile,v3d,v2d,prec)
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_t2, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by SIS
  !       The (hs, hi, t1, t2, cn, ui, uv, etc.) fields are read in from the files so that the o-f data can be computed
  !       by the obsop_XXX.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  !       Instead, the restart file at the analysis time is read in and used to generate the analysis itself.
  !       (see read_restart in this module)
  USE netcdf
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: infile
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1=single, 2=double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4 !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8 !(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .true.
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the T/S netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_diag :: just opened file ", infile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! snow thickness
  varname='h_snow'
  ivid=iv3d_hs

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )

  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select

  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-HS"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hs) = ",MAXVAL(v3d(:,:,k,iv3d_hs))
    enddo
    do k=1,nlev
      WRITE(6,*) "min val for level v3d(:,:,", k, ",iv3d_hs) = ",MINVAL(v3d(:,:,k,iv3d_hs))
    enddo
  endif
  ! !STEVE: end

  !!! ice thickness
  varname='h_ice'
  ivid=iv3d_hi

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-HI"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hi) = ",MAXVAL(v3d(:,:,k,iv3d_hi))
    enddo
    do k=1,nlev
      WRITE(6,*) "min val for level v3d(:,:,", k, ",iv3d_hi) = ",MINVAL(v3d(:,:,k,iv3d_hi))
    enddo
  endif
  ! !STEVE: end

  !!! layer 1 ice temperature
  varname='t_ice1'
  ivid=iv3d_t1

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)
  
  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-T1"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t1) = ",MAXVAL(v3d(:,:,k,iv3d_t1))
    enddo
    do k=1,nlev
      WRITE(6,*) "min val for level v3d(:,:,", k, ",iv3d_t1) = ",MINVAL(v3d(:,:,k,iv3d_t1))
    enddo
  endif
  ! !STEVE: end

  !!! v
  varname='t_ice2'
  ivid=iv3d_t2

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = REAL(buf4,r_size)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = buf8
  end select
  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-T2"
    WRITE(6,*) "read_diag:: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t2) = ",MAXVAL(v3d(:,:,k,iv3d_t2))
    enddo
    do k=1,nlev
      WRITE(6,*) "min val for level v3d(:,:,", k, ",iv3d_t2) = ",MINVAL(v3d(:,:,k,iv3d_t2))
    enddo
  endif
  ! !STEVE: end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D-Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case(prec)
    case(1)
      DEALLOCATE(buf4)
      ALLOCATE(buf4(nlon,nlat,nlev+1))
    case(2)
      DEALLOCATE(buf8)
      ALLOCATE(buf8(nlon,nlat,nlev+1))
  end select

  varname='part_size'
  ivid=iv3d_ps

  WRITE(6,*) "Reading part_size..."
  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,2:nlev+1)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,2:nlev+1),r_sngl)
  end select
  if (dodebug) WRITE(6,*) "read_diag :: just got data for variable part_size"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-part_size"
    WRITE(6,*) "read_diag :: infile = ", infile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t2) = ", MAXVAL(v3d(:,:,k,iv3d_ps))
    enddo
    do k=1,nlev
      WRITE(6,*) "min val for level v3d(:,:,", k, ",iv3d_t2) = ", MINVAL(v3d(:,:,k,iv3d_ps))
    enddo
  endif
! !STEVE: end

  if (dodebug) WRITE(6,*) "read_diag :: finished processing data for variable part_size"

  call check( NF90_CLOSE(ncid) )

  !STEVE: assign concentration as integral of part size across all categories:
  ivid=iv2d_cn
  WRITE(6,*) "Compute sum along categories (dimension 3)..."
! v2d(:,:,ivid) = 0.0d0
! do k=1,nlev
!   do j=1,nlat
!     do i=1,nlon
!       v2d(i,j,ivid) = v2d(i,j,ivid) + v3d(i,j,k,iv3d_ps)
!     enddo
!   enddo
! enddo
  v2d(:,:,ivid) = sum(v3d(:,:,:,iv3d_ps),dim=3)
  WRITE(6,*) "Assign values > 1 to be equal to 1..."
  where(v2d(:,:,ivid)>1) v2d(:,:,ivid)=1.0d0
  WRITE(6,*) "Assign values < 0 to be equal to 0..."
  where(v2d(:,:,ivid)<0) v2d(:,:,ivid)=0.0d0

  WRITE(6,*) "Finished read_diag..."

! STOP(13)
  if (ALLOCATED(buf4)) then
    DEALLOCATE(buf4)
  endif
  if (ALLOCATED(buf8)) then
    DEALLOCATE(buf8)
  endif
  
END SUBROUTINE read_diag

!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the ice_model.res.nc restart files produced by SIS
  !       The hs,hi,t1,t2 prognostic fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsope.f90 program prior
  !       to running letkf.f90, so the diagnostic files should not have to be touched during letkf runtime.
  USE netcdf
  USE params_model, ONLY: basefile
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_t2, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: nij0
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER, INTENT(IN) :: prec ! precision, 1 = single, 2 = double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8
  CHARACTER(slen) :: ncfile
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  !STEVE:
  REAL(r_size) :: meanSSH !STEVE: for temporary SSH estimate based on heat content
  REAL(r_size) :: videpth !STEVE: depth of vertically integrated heat content
  !STEVE: for debugging:
  CHARACTER(32) :: testfile
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .false.
  CHARACTER(3) :: MEM3
  CHARACTER(slen) :: varname
  INTEGER :: ivid

  ncfile = trim(infile)//'.'//trim(basefile)

  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

!! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the ice netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(ncfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart:: just opened file ", ncfile

  !!! thickness of the upper snow layer
  varname='h_snow'
  ivid=iv3d_hs

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable part_size"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable part_size"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-h_snow"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hs) = ",MAXVAL(v3d(:,:,k,iv3d_hs))
    enddo
  endif
! !STEVE: end

  !!! thickness of the lower ice layer
  varname='h_ice'
  ivid=iv3d_hi

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable salt"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable salt"

 !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-SALT"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hi) = ", MAXVAL(v3d(:,:,k,iv3d_hi))
    enddo 
  endif
! !STEVE: end

  !!! temperature of the top half ice layer
  varname='t_ice1'
  ivid=iv3d_t1

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable t_ice1"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable t_ice1"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-t_ice1"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t1) = ",MAXVAL(v3d(:,:,k,iv3d_t1))
    enddo
  endif
! !STEVE: end

  !!! temperature of the bottom half ice layer
  varname='t_ice2'
  ivid=iv3d_t2

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,:)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,:),r_sngl)
  end select

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable t_ice2"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable t_ice2"

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-t_ice2"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t2) = ", MAXVAL(v3d(:,:,k,iv3d_t2))
    enddo 
  endif
! !STEVE: end

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable t_ice2"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable t_ice2"

  !!! part_size, used to compute concentration
  varname='part_size'
  ivid=iv3d_ps

  select case(prec)
    case(1)
      DEALLOCATE(buf4)
      ALLOCATE(buf4(nlon,nlat,nlev+1))
    case(2)
      DEALLOCATE(buf8)
      ALLOCATE(buf8(nlon,nlat,nlev+1))
  end select

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )
  select case(prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      v3d(:,:,:,ivid) = buf4(:,:,2:nlev+1)
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      v3d(:,:,:,ivid) = REAL(buf8(:,:,2:nlev+1),r_sngl)
  end select

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-part_size"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t2) = ", MAXVAL(v3d(:,:,k,iv3d_ps))
    enddo 
  endif
! !STEVE: end

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable part_size"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable part_size"

  call check( NF90_CLOSE(ncid) )

  !STEVE: assign concentration as integral of part size across all categories:
  ivid=iv2d_cn
  v2d(:,:,ivid) = sum(v3d(:,:,:,iv3d_ps),dim=3)
  where(v2d(:,:,ivid)>1) v2d(:,:,ivid)=1.0d0

  ! For additional variables:
  ! E.g. surface fluxes, drifters

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !STEVE: debug test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.false.) then
    testfile = "test_read4.grd"
!   CALL write_grd(trim(testfile),v3d,v2d)

    iunit=55
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=testfile,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

    WRITE(6,*) "Writing to", testfile
    irec=1
    DO n=1,nv3d
      DO k=1,nlev
        WRITE(6,*) "n, k, irec = ", n, k, irec
        WRITE(6,*) "max v3d(n) = ", MAXVAL(v3d(:,:,k,n))
        WRITE(iunit,REC=irec) v3d(:,:,k,n)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      WRITE(iunit,REC=irec) v2d(:,:,n)
      irec = irec + 1
    END DO
    CLOSE(iunit)

    WRITE(6,*) "Initially read from file: ", infile
    WRITE(6,*) "STOP 10"
    STOP(10)
  endif
! !STEVE: debug end

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE
!
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
! Write a set of SIS restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d_in,v2d_in)
  USE netcdf
  USE params_model, ONLY: basefile
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_t2, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),INTENT(IN) :: v3d_in(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d_in(nlon,nlat,nv2d)
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:) !(nlon,nlat,nv2d)
  REAL(r_size), ALLOCATABLE :: t3d(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  CHARACTER(slen) :: ncfile
  INTEGER :: ncid,istat,varid
  INTEGER :: m,k,j,i !STEVE: for debugging
  REAL(r_size), DIMENSION(:), ALLOCATABLE :: mlon, mlat
  REAL(r_size), DIMENSION(:,:,:), ALLOCATABLE :: data4D
  CHARACTER(3) :: MEM3
  LOGICAL, PARAMETER :: do_physlimit=.false.  !STEVE (SIS) need to update for ice requirements

  !STEVE: this is the routine that writes out the individual analysis files for
  !       each esnsemble member in netcdf format.

  ncfile = trim(outfile)//'.'//trim(basefile)

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

  v3d = REAL(v3d_in,r_size)
  v2d = REAL(v2d_in,r_size)

  ! STEVE: for safety, clean up the variables for output:
  if (do_physlimit) then
  do k=1,nlev
    do j=1,nlat
      do i=1,nlon
!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_t) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
!         v3d(i,j,k,iv3d_t) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_s) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
!         v3d(i,j,k,iv3d_s) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_u) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
!         v3d(i,j,k,iv3d_u) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (kmt(i,j) .lt. k .and. v3d(i,j,k,iv3d_v) .ne. 0.0 ) then
!         WRITE(6,*) "WARNING: data on land point in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
!         v3d(i,j,k,iv3d_v) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (k .eq. 1 .and. kmt(i,j) .eq. 0) then 
!         if (v2d(i,j,iv2d_sst) .ne. 0.0 ) v2d(i,j,iv2d_sst) = 0.0 !NF90_FILL_FLOAT
!         if (v2d(i,j,iv2d_sss) .ne. 0.0 ) v2d(i,j,iv2d_sss) = 0.0 !NF90_FILL_FLOAT
!         if (v2d(i,j,iv2d_ssh) .ne. 0.0 ) v2d(i,j,iv2d_ssh) = 0.0 !NF90_FILL_FLOAT
!       endif

!       if (v3d(i,j,k,iv3d_t) < -4) then
!         WRITE(6,*) "WARNING: Bad temp value in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
!         v3d(i,j,k,iv3d_t) = -4.0
!       endif

!       if (v3d(i,j,k,iv3d_t) > 40.0) then
!         WRITE(6,*) "WARNING: Bad temp value in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
!         v3d(i,j,k,iv3d_t) = 40.0
!       endif

!       if (v3d(i,j,k,iv3d_s) < 0 ) then
!         WRITE(6,*) "WARNING: Bad salt value in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
!         v3d(i,j,k,iv3d_s) = 0.0
!       endif

!       if (v3d(i,j,k,iv3d_s) > 50.0) then
!         WRITE(6,*) "WARNING: Bad salt value in analysis output:"
!         WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
!         v3d(i,j,k,iv3d_s) = 50.0
!       endif

!       if (nv2d > 1 .and. k .eq. 1) then
!         if (v2d(i,j,iv2d_sst) < -4) v2d(i,j,iv2d_sst) = -4.0
!         if (v2d(i,j,iv2d_sst) > 40.0) v2d(i,j,iv2d_sst) = 40.0
!         if (v2d(i,j,iv2d_sss) < 0) v2d(i,j,iv2d_sss) = 0.0
!         if (v2d(i,j,iv2d_sss) > 50.0) v2d(i,j,iv2d_sss) = 50.0
!       endif

      enddo
    enddo
  enddo
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !STEVE: open netcdf diagnostic/restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(ncfile,NF90_WRITE,ncid) )

  !!! h_snow
  !STEVE: for debugging
  if (.false.) then
  do m=1,nv3d
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
!       if ( isnan( REAL(v3d(i,j,k,m),r_size) ) )then
!         WRITE(6,*) "common_sis.f90::write_grd:: ERROR: found NaN..."
!         WRITE(6,*) "v3d(i,j,k,m) contains NaN. i,j,k,m = ", i,j,k,m
!         STOP 1
!       endif
        enddo
      enddo
    enddo
  enddo
  endif

  call check( NF90_INQ_VARID(ncid,'h_snow',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_hs)) )

  !!! h_ice
  call check( NF90_INQ_VARID(ncid,'h_ice',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_hi)) )

  !!! t_ice1
  call check( NF90_INQ_VARID(ncid,'t_ice1',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_t1)) )

  !!! t_ice2
  call check( NF90_INQ_VARID(ncid,'t_ice2',varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_t2)) )

  call check( NF90_CLOSE(ncid) )

END SUBROUTINE write_restart


!-----------------------------------------------------------------------
! Read a grid file with single precision
!-----------------------------------------------------------------------
SUBROUTINE read_grd(filename,v3d,v2d)
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: nij0
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec

! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))

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

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

END SUBROUTINE read_grd

!-----------------------------------------------------------------------
! Write a grid file in single precision
!-----------------------------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: nij0
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: iunit,iolen
  INTEGER :: i,j,k,n,irec
  LOGICAL, PARAMETER :: dodebug=.false.

  if (dodebug) print *, "write_grd:: open filename = ",filename
  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  if (dodebug) print *, "write_grd:: nij0,iolength = ", nij0,iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)

  irec=1
  DO n=1,nv3d
    DO k=1,nlev
      if (dodebug) print *, "write_grd:: n,k,irec = ",n,k,irec
      WRITE(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO
  END DO

  DO n=1,nv2d
    if (dodebug) print *, "write_grd:: n,irec = ",n,irec
    WRITE(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
    irec = irec + 1
  END DO

  CLOSE(iunit)

END SUBROUTINE write_grd
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: element
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

END SUBROUTINE monit_grd

!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

! ALLOCATE(v3dm(nij,nlev,nv3d),v2dm(nij,nv2d))

  DO n=1,nv3d
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
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
