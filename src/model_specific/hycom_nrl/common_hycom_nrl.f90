MODULE common_oceanmodel
!=======================================================================
!
! [PURPOSE:] Common Information for HYCOM
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi created
!   01/23/2009 Takemasa Miyoshi modified
!   04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   06/08/2015 Steve Penny converted for use with HYCOM
!   12/12/2017 Steve Penny modified for use with NCODA at NRL
!
!=======================================================================
  USE common
  USE params_model, ONLY: glon, glat, glev
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_letkf, ONLY: DO_ALTIMETRY, DO_DRIFTERS, DO_SLA
  USE vars_model,   ONLY: lon, lat, lev, lon2d, lat2d
  USE vars_model,   ONLY: lon0, lonf, lat0, latf, wrapgap, kmt, kmt0
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
  USE params_model, ONLY: initialize_params_model
  USE vars_model,   ONLY: initialize_vars_model, set_vars_model
  USE params_model, ONLY: nlon,nlat,nlev,istart,iend,jstart,jend
  USE params_model, ONLY: glon, glat, glev
  USE vars_model,   ONLY: lev,z_lvl
  USE input_nml_oceanmodel, ONLY: read_ncoda_namelists
  USE hycom_io,     ONLY: read_grdfiles, read_kmt

  INTEGER :: i,j,k
  INTEGER :: ncid,ncid2,ncid3,istat,varid,dimid
  LOGICAL :: ex
  INTEGER :: fid=21
  INTEGER :: ierr

  ! For debugging:
  LOGICAL :: dodebug = .true.

  !-----------------------------------------------------------------------------
  ! Commence
  !-----------------------------------------------------------------------------
  WRITE(6,'(A)') 'Hello from set_common_oceanmodel'

  !-----------------------------------------------------------------------------
  ! Get grid dimensions
  !
  ! Using glon,glat,glev (specified in the input.nml namelist) for the global 
  ! grid and nlon,nlat,nlev for the LAM subgrid (specified in an input grid
  ! specification file)
  !-----------------------------------------------------------------------------

  ! Open NCODA namelists with grid dimensions and z-level definitions
  WRITE(6,*) "Read NCODA namelists for grid dimensions and z level coordinates."
  WRITE(6,*) "Calling read_ncoda_namelists..."
  CALL read_ncoda_namelists
  WRITE(6,*) "Finished calling read_ncoda_namelists."

  !-----------------------------------------------------------------------------
  ! Set subgrid grid dimensions using range indices (from input.nml)
  !-----------------------------------------------------------------------------
  ! Also read the LAM start and end indices: istart,iend,jstart,jend
  !STEVE: Set this in the input.nml or some other input file (ISSUE)
  ! WARNING: Assuming here that the grid does not overlap a boundary

  if (dodebug) then
    WRITE(6,*) "glon = ", glon
    WRITE(6,*) "glat = ", glat
    WRITE(6,*) "glev = ", glev
    WRITE(6,*) "istart,iend = ", istart, iend
    WRITE(6,*) "jstart,jend = ", jstart, jend
  endif

  ! Set subgrid to whole domain if they are not input via namelist:
  if (iend == 0) then
    istart = 1
    iend = glon
  endif
  if (jend == 0) then
    jstart = 1
    jend = glat
  endif

  ! Calculate the grid dimensions for the analysis domain
  nlon   = abs(iend - istart) + 1
  nlat   = abs(jend - jstart) + 1
  nlev   = glev

  if (dodebug) then
    WRITE(6,*) "nlon,nlat,nlev = ", nlon, nlat, nlev
  endif

  !-----------------------------------------------------------------------------
  ! Initialize the model parameters and variables
  !-----------------------------------------------------------------------------
  CALL initialize_params_model ! (checks to make sure it is initialized)
  CALL initialize_vars_model   ! (checks to make sure it is initialized)

  !-----------------------------------------------------------------------------
  ! Get the -local subgrid- 2d grid specifications
  !-----------------------------------------------------------------------------
  WRITE(6,*) "set_common_hycom:: Reading NCODA grd files for grid coordinates..."
  CALL read_grdfiles(lon2d,lat2d)
  WRITE(6,*) "set_common_hycom:: finished read_grdfiles."

  !-----------------------------------------------------------------------------
  ! Get the kmt information (2d land/sea mask specifying depth of ocean grid)
  !-----------------------------------------------------------------------------
  WRITE(6,*) "set_common_hycom:: Reading NCODA land/sea mask file for kmt data..."
  CALL read_kmt(kmt)
  WRITE(6,*) "set_common_hycom:: finished read_kmt."
  kmt0 = kmt  ! Convert to real (r_size)

  !-----------------------------------------------------------------------------
  ! Set model variables that depend on initialization and further processing.
  !-----------------------------------------------------------------------------
  ! (e.g. lon0, lat0, lonf, latf, wrapgap, ...)
  WRITE(6,*) "set_common_hycom:: Setting model variables..."
  CALL set_vars_model
  WRITE(6,*) "set_common_hycom:: finished set_vars_model."

  if (dodebug) then
    if (allocated(lon2d)) then
      WRITE(6,*) "pre-assignment:"
      WRITE(6,*) "SIZE(lon2d,1) = ", SIZE(lon2d,1)
      WRITE(6,*) "SIZE(lon2d,2) = ", SIZE(lon2d,2)
      WRITE(6,*) "lon2d(1,1) = ", lon2d(1,1)
      WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)
    else
      WRITE(6,*) "lon2d is not allocated"
    endif

    if (allocated(lat2d)) then
      WRITE(6,*) "pre-assignment:"
      WRITE(6,*) "SIZE(lat2d,1) = ", SIZE(lat2d,1)
      WRITE(6,*) "SIZE(lat2d,2) = ", SIZE(lat2d,2)
      WRITE(6,*) "lat2d(1,1) = ", lat2d(1,1)
      WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)
    else
      WRITE(6,*) "lat2d is not allocated"
    endif

    if (allocated(lev)) then
      WRITE(6,*) "pre-assignment:"
      WRITE(6,*) "SIZE(lev,1) = ", SIZE(lev,1)
      WRITE(6,*) "lev(1) = ", lev(1)
      WRITE(6,*) "lev(nlev) = ", lev(nlev)
    else
      WRITE(6,*) "lev is not allocated"
    endif
  endif

  WRITE(6,*) "set_common_hycom:: Finished."

END SUBROUTINE set_common_oceanmodel


!STEVE: ISSUE: need to update for HYCOM
SUBROUTINE read_etaclm
!===============================================================================
! Read in (e.g.) 1993-1999 model climatology to subtract as done by AVISO
!===============================================================================
!STEVE: commented out to compile without netcdf library support
! USE netcdf

  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: i,j
  INTEGER :: ncid, varid

  ! read the model SSH climatology netcdf file
  ! read into: SSHclm_m
! CALL check( NF90_OPEN(SSHclm_file,NF90_NOWRITE,ncid) )
! WRITE(6,*) "read_etaclm:: just opened file ", SSHclm_file

! buf4=0.0
! CALL check( NF90_INQ_VARID(ncid,'ssh',varid) )
! CALL check( NF90_GET_VAR(ncid,varid,buf4) )
! do j=1,nlat
!   do i=1,nlon
!     !STEVE: Hopefully reading in meters here... (data might be in cm)
!     SSHclm_m(i,j) = REAL(buf4(i,j),r_size)
!   enddo
! enddo

! CALL check( NF90_CLOSE(ncid) )

END SUBROUTINE read_etaclm


SUBROUTINE read_diag(infile,v3d,v2d,prec_in)
!===============================================================================
! Read a archive/restart file in binary format extracted from hycom ab-format
!===============================================================================
  !STEVE: This subroutine reads the hourly/daily diagnostic files produced by HYCOM
  !       The t,s,u,v,(ssh) fields are read in from the ab archive files so that the o-f data can be computed
  !       by the obsop.f90 program prior to running letkf.f90, so the diagnostic files do not have 
  !       to be touched during letkf runtime.
  USE hycom_io,     ONLY: hycom_undef,hycom_eps
  USE params_model, ONLY: base,base_a,base_b
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h

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
    STOP (126)  ! CDA: missing STOP number, temporarily set 126
  endif

  ! STEVE: this is provided externally at the moment
  binfile = trim(infile)//trim(base)
  infile_a = trim(infile)//trim(base_a)
  infile_b = trim(infile)//trim(base_b)


  if (dodebug) WRITE(6,*) "read_diag:: Pre-read_hycom: Reading file: ", trim(binfile)
! JILI Read HYCOM a and b files
! CALL get_hycom(infile_a,infile_b,v3d,v2d)

  if (dodebug) WRITE(6,*) "read_diag:: Post-read_hycom: Just read file: ", trim(binfile)

  if (dodebug) then
    WRITE(6,*) "read_diag:: Post-read_hycom"
    buf8 = v3d(:,:,:,iv3d_t)
    where (abs(buf8-hycom_undef) < hycom_eps) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_t)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_t)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_s)
    where (abs(buf8-hycom_undef) < hycom_eps) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_s)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_s)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_u)
    where (abs(buf8-hycom_undef) < hycom_eps) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_u)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_u)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_v)
    where (abs(buf8-hycom_undef) < hycom_eps) buf8 = 0.0d0
    WRITE(6,*) "MAXVAL(v3d(:,:,1,iv3d_v)) = ", MAXVAL(buf8(:,:,1))
    WRITE(6,*) "MINVAL(v3d(:,:,1,iv3d_v)) = ", MINVAL(buf8(:,:,1))

    buf8 = v3d(:,:,:,iv3d_h)
    where (abs(buf8-hycom_undef) < hycom_eps) buf8 = 0.0d0
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
  !       It is assumed that the o-f data have already been computed by the obsop_xxx.f90 program prior
  !       to running letkf.f90, so the archive (diagnostic) files should not have to be touched during letkf runtime.
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h
  USE params_model, ONLY: rsrt_tbase,rsrt_sbase,rsrt_ubase,rsrt_vbase,rsrt_sshbase
  USE hycom_io,     ONLY: read_hycom_ncoda

  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:),  INTENT(OUT) :: v2d !(nlon,nlat,nv2d)
  INTEGER, INTENT(IN), OPTIONAL :: prec_in ! precision, 1=single, 2=double
  INTEGER :: prec
  CHARACTER(slen) :: filename
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:,:)
  LOGICAL :: dodebug = .true.

  ! Check for input argument for precision
  !STEVE: (not curerntly supported)
  if (PRESENT(prec_in)) then
    prec = prec_in
  else
    prec = 1 ! default is 1
  endif

! !STEVE: using the HYCOM/NCODA bacgkround format
  ALLOCATE(v3d(nlon,nlat,nlev,nv3d),v2d(nlon,nlat,nv2d))
  CALL read_hycom_ncoda(infile,v3d,v2d)

END SUBROUTINE read_restart

!-----------------------------------------------------------------------
subroutine check(status)
!STEVE: commented out to compile without netcdf library support
! USE netcdf
  INTEGER, intent (in) :: status

! if(status /= nf90_noerr) then 
!   print *, trim(nf90_strerror(status))
!   stop "Stopped"
! end if
end subroutine check


!-----------------------------------------------------------------------
! Write a set of MOM6 restart files to initialize the next model run
!-----------------------------------------------------------------------
SUBROUTINE write_restart(outfile,v3d,v2d)
  !STEVE: This writes out the analysis to a pre-existing template
  !       as defined in the hycom_io module
  USE hycom_io,     ONLY: write_hycom_ncoda
  USE params_model, ONLY: do_physlimit, nlev, nlat, nlon
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_ssh
  USE params_model, ONLY: min_t, max_t, min_s, max_s, min_uv, max_uv, min_ssh, max_ssh

  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),DIMENSION(:,:,:,:),INTENT(INOUT) :: v3d !(nlon,nlat,nlev,nv3d) !STEVE: 'INOUT' allows for modifications to v3d and v2d
  REAL(r_sngl),DIMENSION(:,:,:),  INTENT(INOUT) :: v2d !(nlon,nlat,nv2d)

  INTEGER :: m,k,j,i !STEVE: for debugging

  LOGICAL :: dodebug = .false.

  if (dodebug) WRITE(6,*) "write_restart:: begin..."
  
  ! STEVE: for safety, clean up the variables for output:
  ! JILI for land grids, also set variables to undef
  ! JILI a potential probelm: bottom layers large spread-need to be very careful
  physlimit : if (do_physlimit) then
    do k=1,nlev
      do j=1,nlat
        do i=1,nlon
          ! Check allowable temperature bounds
          if (v3d(i,j,k,iv3d_t) < min_t) then
            WRITE(6,*) "WARNING: Bad temp value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
            v3d(i,j,k,iv3d_t) = min_t
          elseif (v3d(i,j,k,iv3d_t) > max_t) then
            WRITE(6,*) "WARNING: Bad temp value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_t)
            v3d(i,j,k,iv3d_t) = max_t
          endif
  
          ! Check allowable salinity bounds
          if (v3d(i,j,k,iv3d_s) < min_s) then
            WRITE(6,*) "WARNING: Bad salt value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
            v3d(i,j,k,iv3d_s) = min_s
          elseif (v3d(i,j,k,iv3d_s) > max_s) then
            WRITE(6,*) "WARNING: Bad salt value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_s)
            v3d(i,j,k,iv3d_s) = max_s
          endif
  
          ! Check allowable u current bounds
          if (v3d(i,j,k,iv3d_u) < min_uv) then
            WRITE(6,*) "WARNING: Bad u-vel value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
            v3d(i,j,k,iv3d_u) = min_uv
          elseif (v3d(i,j,k,iv3d_u) > max_uv) then
            WRITE(6,*) "WARNING: Bad u-vel value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_u)
            v3d(i,j,k,iv3d_u) = max_uv
          endif

          ! Check allowable v current bounds
          if (v3d(i,j,k,iv3d_v) < min_uv) then
            WRITE(6,*) "WARNING: Bad v-vel value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
            v3d(i,j,k,iv3d_v) = min_uv
          elseif (v3d(i,j,k,iv3d_v) > max_uv) then
            WRITE(6,*) "WARNING: Bad v-vel value in analysis output:"
            WRITE(6,*) "v3d(",i,",",j,",",k,") = ", v3d(i,j,k,iv3d_v)
            v3d(i,j,k,iv3d_v) = max_uv
          endif

          ! Check allowable ssh bounds
          if (k==1) then
            if (iv2d_ssh <= nv2d) then
              if (v2d(i,j,iv2d_ssh) < min_ssh) then
                WRITE(6,*) "WARNING: Bad ssh value in analysis output:"
                WRITE(6,*) "v2d(",i,",",j,") = ", v2d(i,j,iv2d_ssh)
                v2d(i,j,iv2d_ssh) = min_ssh
              elseif (v2d(i,j,iv2d_ssh) > max_ssh) then
                WRITE(6,*) "WARNING: Bad ssh value in analysis output:"
                WRITE(6,*) "v2d(",i,",",j,") = ", v2d(i,j,iv2d_ssh)
                v2d(i,j,iv2d_ssh) = max_ssh
              endif
            endif
          endif

        enddo
      enddo
    enddo

  endif physlimit

  ! Write local analysis
  if (dodebug) WRITE(6,*) "write_restart:: writing local analysis..."
  CALL write_hycom_ncoda(trim(outfile),v3d,v2d)
  if (dodebug) WRITE(6,*) "write_restart:: complete." 

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
  use hycom_io, ONLY: hycom_undef,hycom_eps

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
          if (abs(v3d(i,k,m,n)-hycom_undef) > hycom_eps )  then
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
        if (abs(v2d(i,m,n)-hycom_undef) > hycom_eps )  then
          v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
          cnt2d(i,n) = cnt2d(i,n) + 1
        endif
      enddo
      v2dm(i,n) = v2dm(i,n) / REAL(cnt2d(i,n),r_size)
    enddo
  enddo

END SUBROUTINE ensmean_grd


END MODULE common_oceanmodel

