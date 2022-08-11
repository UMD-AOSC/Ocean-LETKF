MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for CICE
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   04/26/2011 Steve Penny, converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   04/06/2015 Steve Penny converted for use with SIS
!   03/14/2017 Steve Penny converted for use with CICE
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
  USE params_model, ONLY: nlon, nlat, nlev !STEVE: for debugging
  USE params_model, ONLY: gridfile
  USE params_model, ONLY: grid_lon_name, grid_lat_name, grid_lev_name
  USE params_model, ONLY: grid_lon2d_name, grid_lat2d_name
  USE params_model, ONLY: grid_wet_name, grid_kmt_name
  USE params_model, ONLY: initialize_params_model
  USE vars_model,   ONLY: initialize_vars_model
  USE vars_model,   ONLY: lon, lat, lev, lon2d, lat2d
  USE vars_model,   ONLY: kmt0, kmt, wet
  USE vars_model,   ONLY: set_vars_model

  IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname
  LOGICAL :: ex
  LOGICAL :: doverbose = .true.

  WRITE(6,'(A)') 'Hello from set_common_cice5'
  CALL initialize_params_model ! (checks to make sure it is initialized)
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  !STEVE: still using the grid_spec.nc file from GFDL MOM5 here for the grid definitions: (ISSUE)

  !-----------------------------------------------------------------------------
  ! Lon, Lat, f, orography
  !-----------------------------------------------------------------------------
  INQUIRE(FILE=trim(gridfile),EXIST=ex)
  IF(.not. ex) THEN
    WRITE(6,*) "The file does not exist: ", gridfile
    WRITE(6,*) "Exiting common_cice.f90..."
    STOP (2)
  ENDIF
  WRITE(6,'(A)') '  >> accessing file: ', gridfile

  !-----------------------------------------------------------------------------
  ! Get 1d and 2d longitude fields:
  !-----------------------------------------------------------------------------
  CALL check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )
  CALL check( NF90_INQ_VARID(ncid,grid_lon_name,varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lon) )

  if (doverbose) then
    WRITE(6,*) "lon(1) = ", lon(1)
    WRITE(6,*) "lon(nlon) = ", lon(nlon)
  endif

  CALL check( NF90_INQ_VARID(ncid,grid_lon2d_name,varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lon2d) )

  if (doverbose) then
    WRITE(6,*) "lon2d(1,1)       = ", lon2d(1,1)
    WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)
  endif

  !-----------------------------------------------------------------------------
  ! Get 1d and 2d latitude fields:
  !-----------------------------------------------------------------------------
  CALL check( NF90_INQ_VARID(ncid,grid_lat_name,varid) )   ! Latitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lat) )

  if (doverbose) then
    WRITE(6,*) "lat(1) = ", lat(1)
    WRITE(6,*) "lat(nlat) = ", lat(nlat)
  endif

  CALL check( NF90_INQ_VARID(ncid,grid_lat2d_name,varid) )   ! Longitude for T-cell
  CALL check( NF90_GET_VAR(ncid,varid,lat2d) )

  if (doverbose) then
    WRITE(6,*) "lat2d(1,1)       = ", lat2d(1,1)
    WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)
  endif

  ! 
  ! Lev data
  !
  call check( NF90_INQ_VARID(ncid,grid_lev_name,varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  !
  ! kmt data
  !
  call check( NF90_INQ_VARID(ncid,grid_kmt_name,varid) ) ! number of vertical T-cells
  call check( NF90_GET_VAR(ncid,varid,kmt0) )
  WRITE(6,*) "kmt0(1,1) = ", kmt0(1,1)
  WRITE(6,*) "kmt0(nlon,nlat) = ", kmt0(nlon,nlat)
  kmt = NINT(kmt0)
  call check( NF90_INQ_VARID(ncid,grid_wet_name,varid) )        ! land/sea flag (0=land) for T-cell
  call check( NF90_GET_VAR(ncid,varid,wet) )
  WRITE(6,*) "wet(1,1) = ", wet(1,1)
  WRITE(6,*) "wet(nlon,nlat) = ", wet(nlon,nlat)

  ! Set model variables that depend on initialization and further processing.
  ! (e.g. lon0, lat0, lonf, latf, wrapgap, ...)
  CALL set_vars_model

END SUBROUTINE set_common_oceanmodel

!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!STEVE: add this:
!-- Read a diagnostic file in cice5 netcdf format ---------------------------------------------------
SUBROUTINE read_diag(infile,v3d,v2d,prec_in)
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: diag_hs_name, diag_hi_name
  USE params_model, ONLY: diag_t1_name, diag_t2_name
  USE params_model, ONLY: diag_ps_name
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by CICE5
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
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:) :: buf4 !(nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:,:) :: buf8 !(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  INTEGER :: iunit,iolen,n,irec
  LOGICAL, PARAMETER :: dodebug = .true.
  CHARACTER(slen) :: varname
  INTEGER :: ivid
  INTEGER :: prec
  REAL(r_sngl), ALLOCATABLE :: v3d0(:,:,:,:)
  REAL(r_sngl), ALLOCATABLE :: v2d0(:,:,:)

  ! If prec is a provided argument, use indicated precision,
  if (PRESENT(prec_in)) then
    prec = prec_in 
  else
    ! otherwise default to single precision.
    prec = 1
    WRITE(6,*) "common_cice.f90::read_diag:: using default precision = ", prec
  endif

  !STEVE: for now, assume all input files are only restart files.
  !       Preferably, a set of minimal diagnostic files would be
  !       use instead, since this is only called for computing the
  !       model equivalent H(xb) in the observation operator in order
  !       to compute the obs-minus-forecast later in letkf.f90
  ALLOCATE(v3d0(nlon,nlat,nlev,nv3d), v2d0(nlon,nlat,nv2d))
  CALL read_restart(infile,v3d0,v2d0,prec)
  v3d = REAL(v3d0,r_size)
  v2d = REAL(v2d0,r_size)
  
END SUBROUTINE read_diag

!-----------------------------------------------------------------------
!-- Read a restart file at the analysis time  --------------------------
SUBROUTINE read_restart(infile,v3d,v2d,prec)
  !STEVE: This subroutine reads the cice5_model.res.nc restart files produced by CICE5
  !       The aicen,vicen,tsfcn,uvel,vvel prognostic fields are read in from a file so that the transform can be applied to them.
  !       It is assumed that the o-f data have already been computed by the obsope_XXXX.f90 program prior
  !       to running letkf.f90, so the diagnostic files should not have to be touched during letkf runtime.
  USE netcdf
  USE params_model, ONLY: basefile
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: nij0
  USE params_model, ONLY: rsrt_hs_name, rsrt_hi_name
  USE params_model, ONLY: rsrt_t1_name, rsrt_t2_name
  USE params_model, ONLY: rsrt_ps_name
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
  CHARACTER(3) :: MEM3
  CHARACTER(slen) :: varname
  INTEGER :: ivid
  LOGICAL, PARAMETER :: dodebug = .true.

  ncfile = trim(infile)//'.'//trim(basefile)

  select case(prec)
    case(1)
      if (dodebug) print *, "read_restart:: using single precision"
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      if (dodebug) print *, "read_restart:: using double precision"
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

!! ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the ice netcdf restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(ncfile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_restart:: just opened file ", ncfile

  !!! volume per unit area of snow (in category n) m
  varname=rsrt_hs_name
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

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable vsnon"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable vsnon"

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-h_snow"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hs) = ",MAXVAL(v3d(:,:,k,iv3d_hs))
    enddo
  endif
! !STEVE: end

  !!!  volume per unit area of ice (in category n) m
  varname=rsrt_hi_name
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

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable vicen"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable vicen"

 !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-h_ice"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_hi) = ", MAXVAL(v3d(:,:,k,iv3d_hi))
    enddo 
  endif
! !STEVE: end

  !!! temperature of ice/snow top surface (in category n) C
  varname=rsrt_t1_name
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

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable ", varname
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable ", varname

  ! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-t_ice1"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_t1) = ",MAXVAL(v3d(:,:,k,iv3d_t1))
    enddo
  endif
! !STEVE: end

  !!! total concentration of ice in grid cell (in category n)
  !   (in SIS part_size) used to compute concentration
  varname=rsrt_ps_name
  ivid=iv3d_ps

  select case(prec)
    case(1)
      DEALLOCATE(buf4)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      DEALLOCATE(buf8)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

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

! !STEVE: debug
  if (dodebug) then
    WRITE(6,*) "POST-part_size"
    WRITE(6,*) "read_restart :: ncfile = ", ncfile
    do k=1,nlev
      WRITE(6,*) "max val for level v3d(:,:,", k, ",iv3d_ps) = ", MAXVAL(v3d(:,:,k,iv3d_ps))
    enddo 
  endif
! !STEVE: end

  if (dodebug) WRITE(6,*) "read_restart :: just got data for variable aicen"
  if (dodebug) WRITE(6,*) "read_restart :: finished processing data for variable aicen"

  ! CLOSE NETCDF FILE
  call check( NF90_CLOSE(ncid) )

  !STEVE: assign concentration as integral of part size across all categories:
  ivid=iv2d_cn
  v2d(:,:,ivid) = sum(v3d(:,:,:,iv3d_ps),dim=3)
  where(v2d(:,:,ivid)>1) v2d(:,:,ivid)=1.0d0

! DEALLOCATE(v3d,v2d) !INTENT OUT, so no DEALLOCATE

  !STEVE: clean up undefined values:
  WHERE (ABS(v3d) .ge. vmax) v3d = 0.0
  WHERE (ABS(v2d) .ge. vmax) v2d = 0.0

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
! Write a set of CICE5 restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d_in,v2d_in)
  USE netcdf
  USE params_model, ONLY: basefile
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_ps, iv2d_cn
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_model, ONLY: rsrt_hs_name, rsrt_hi_name
  USE params_model, ONLY: rsrt_t1_name, rsrt_t2_name
  USE params_model, ONLY: rsrt_ps_name
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !STEVE: open netcdf diagnostic/restart file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call check( NF90_OPEN(ncfile,NF90_WRITE,ncid) )

  !!! vsnon
  call check( NF90_INQ_VARID(ncid,rsrt_hs_name,varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_hs)) )

  !!! vicen
  call check( NF90_INQ_VARID(ncid,rsrt_hi_name,varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_hi)) )

  !!! tsfcn
  call check( NF90_INQ_VARID(ncid,rsrt_t1_name,varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_t1)) )

  !!! aicen
  call check( NF90_INQ_VARID(ncid,rsrt_ps_name,varid) )
  call check( NF90_PUT_VAR(ncid,varid,v3d(:,:,:,iv3d_ps)) )

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
  LOGICAL, PARAMETER :: dodebug=.true.

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
