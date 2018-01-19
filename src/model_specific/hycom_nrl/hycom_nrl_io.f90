MODULE hycom_io
  USE common, ONLY: r_size, r_sngl, slen
  USE params_model
! USE mod_za

  IMPLICIT NONE

  PUBLIC

  REAL(r_size), PARAMETER :: hycom_undef=2.0**100
  REAL(r_size), PARAMETER :: hycom_eps = EPSILON(1.0)  !(intentionally at real precision)              
  
CONTAINS 



SUBROUTINE read_hycom_ncoda(infile,v3d,v2d)
!===============================================================================
! Read in an ncoda-formatted hycom background file
!===============================================================================
  !STEVE: This subroutine reads the HYCOM/NCODA restart file.
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h
  USE params_model, ONLY: rsrt_tbase,rsrt_sbase,rsrt_ubase,rsrt_vbase,rsrt_sshbase

  CHARACTER(*),INTENT(IN) :: infile
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:),  INTENT(OUT) :: v2d !(nlon,nlat,nv2d)
  CHARACTER(slen) :: filename
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:,:)
  LOGICAL :: dodebug = .true.

  allocate( buf4(nlon,nlat,nlev) )
  allocate( v3d(nlon,nlat,nlev,nv3d), v2d(nlon,nlat,nv2d) )

  !
  ! Temperature
  !
  filename = trim(infile)//'.'//trim(rsrt_tbase)
  if (dodebug) write(6,*)'read_ens_mpi, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,nlev)
  v3d(:,:,:,iv3d_t) = buf4(:,:,:)
  if (dodebug) then
    write (6, *) "v3d(m/2,y/2,l,iv3d_t)= ", v3d(nlon/2,nlat/2,1,iv3d_t)
    write (6, *) " max/min v3d = ",  maxval(buf4), minval(buf4)
  endif

  !
  ! Salinity
  !
  filename = trim(infile)//'.'//trim(rsrt_sbase)
  if (dodebug) write(6,*)'read_ens_mpi,  reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,nlev)
  v3d(:,:,:,iv3d_s) = buf4(:,:,:)
  if (dodebug) then
    write (6, *) "v3d(m/2,y/2,l,iv3d_s)= ", v3d(nlon/2,nlat/2,1,iv3d_s)
    write (6, *) " max/min v3d = ",  maxval(buf4), minval(buf4)
  endif

  !
  ! u
  !
  filename = trim(infile)//'.'//trim(rsrt_ubase)
  if (dodebug) write(6,*)'read_ens_mpi, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,nlev)
  v3d(:,:,:,iv3d_u) = buf4(:,:,:)
  if (dodebug) then 
    write (6, *) "v3d(m/2,y/2,l,iv3d_u)= ", v3d(nlon/2,nlat/2,1,iv3d_u)
    write (6, *) " max/min v3d = ",  maxval(buf4), minval(buf4)
  endif

  !
  !  v
  ! 
  filename = trim(infile)//'.'//trim(rsrt_vbase)
  if (dodebug) write(6,*)'read_ens_mpi, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,nlev)
  v3d(:,:,:,iv3d_v) = buf4(:,:,:)
  if (dodebug) then
    write (6, *) "v3d(m/2,y/2,l,iv3d_v)= ", v3d(nlon/2,nlat/2,1,iv3d_v)
    write (6, *) " max/min v3d = ",  maxval(buf4), minval(buf4)
  endif

  !
  ! ssh
  !
  filename = trim(infile)//'.'//trim(rsrt_sshbase)
  if (dodebug) write(6,*)'read_ens_mpi, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4(:,:,1:1),1)
  v2d(:,:,iv2d_ssh) = buf4(:,:,1)
  if (dodebug) then
    write (6,*) "v2d(m/2,n/2,iv2d_ssh)= ", v2d(nlon/2,nlat/2,iv2d_ssh)
    write (6,*) " max/min v2d = ",  maxval(v2d(:,:,iv2d_ssh)), minval(v2d(:,:,iv2d_ssh))
  endif

END SUBROUTINE read_hycom_ncoda



SUBROUTINE read_ncoda_bckgrnd(filename, indata, kmax)
!===============================================================================
! Read the NCODA format background file for a specific slice defined in
! params_model
!===============================================================================
! Initial Installation - April 1994 -- Cummings, J.
! M.Wei, based on rd_ocn_anl.f
! STEVE: modify to read a segment of the total grid

  USE params_model, ONLY: glon, glat, glev
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: istart, iend, jstart, jend

  CHARACTER(*), INTENT(IN) :: filename
  REAL(r_sngl), DIMENSION(:,:,:), INTENT(OUT) :: indata
  INTEGER, INTENT(IN) :: kmax  ! Maximum number of levels to read
  REAL(r_sngl), DIMENSION(:,:,:), ALLOCATABLE :: buf4
  INTEGER :: is,ie,js,je
  INTEGER, PARAMETER :: fid=21
  LOGICAL   :: ex
  INTEGER :: mx, ny, lz, reclen, iolen
  LOGICAL :: dodebug = .false.

  !
  ! global  ..long integer array dimensions
  !
  is = istart
  ie = iend
  js = jstart
  je = jend

  mx = glon
  ny = glat
  lz = min(glev,kmax)

  INQUIRE(IOLENGTH=iolen) iolen
  reclen = mx * ny * lz * iolen
  allocate(buf4(mx,ny,lz)) !STEVE: TEMPORORAY - replace with direct access read of subgrid

  if (dodebug) then
    WRITE(6,*) "mx = ", mx
    WRITE(6,*) "ny = ", ny
    WRITE(6,*) "lz = ", lz
    WRITE(6,*) "reclen = ", reclen
  endif

  if (dodebug) then
    WRITE(6,*) "is = ", is
    WRITE(6,*) "ie = ", ie
    WRITE(6,*) "js = ", js
    WRITE(6,*) "je = ", je
    WRITE(6,*) "lz = ", lz
    WRITE(6,*) "SHAPE(buf4) = ", SHAPE(buf4)
    WRITE(6,*) "SHAPE(indata) = ", SHAPE(indata)
  endif

  !---------------------------------
  ! Read the NCODA background file
  !---------------------------------
  inquire (file=trim(filename), exist=ex)
  if (ex) then
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
!STEVE: (ISSUE) replace with stream I/O to read only desired section and cut down on I/O costs
!   open (unit=fid, file=trim(filename), status='old', access='stream', form='unformatted', recl=reclen)
    read (fid,rec=1) buf4
    close (fid)
    indata = buf4(is:ie,js:je,1:lz) 
    WRITE(6,*) "read_ncoda_bckgrnd:: Finished assigning buf4 to indata."
  else
    write (6,*) 'read_ncoda_bckgrnd:: Background file missing: ', trim(filename)
    STOP("read_ncoda_bckgrnd:: EXITING...")
  endif

  if (dodebug) then
    write (6, *) "buf4(m/2,y/2,1)= ", buf4(mx/2,ny/2,1)
    write (6, *) " max/min buf4 = ",  maxval(buf4), minval(buf4)
    write (6, *) "indata(m/2,y/2,1)= ", indata(mx/2,ny/2,1)
    write (6, *) " max/min indata = ",  maxval(indata), minval(indata)
  endif

  deallocate(buf4)

END SUBROUTINE read_ncoda_bckgrnd


SUBROUTINE read_ncoda_intfile(filename, indata, kmax)
!===============================================================================
! Read the NCODA format file using integers
!===============================================================================
! Initial Installation - April 1994 -- Cummings, J.
! STEVE: modify to read a segment of the total grid

  USE params_model, ONLY: glon, glat, glev
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: istart, iend, jstart, jend

  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: indata
  INTEGER, INTENT(IN) :: kmax
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: buf4
  INTEGER :: is,ie,js,je
  INTEGER, PARAMETER :: fid=20
  LOGICAL   :: ex
  INTEGER*8 :: mx, ny, lz, reclen

  !
  ! global  ..long integer array dimensions
  !
  is = istart
  ie = iend
  js = jstart
  je = jend

  mx = glon
  ny = glat
  lz = min(glev,kmax)

  reclen = mx * ny * lz * 4
  allocate(buf4(glon,glat,lz)) !STEVE: TEMPORORAY - replace with direct access read of subgrid

  !---------------------------------
  ! Read the NCODA integer (e.g. land/sea mask) file
  !---------------------------------
  inquire (file=trim(filename), exist=ex)
  if (ex) then
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    read (fid,rec=1) buf4
    close (fid)
    write (6, '( ''read_ncoda_intfile:: Finished reading background file: '', a)' ) trim(filename)
    !STEVE: (ISSUE) replace with stream I/O to read only desired section and cut
    !down on I/O costs
    indata = buf4(is:ie,js:je,1:lz) 
  else
    write (6, '( ''read_ncoda_intfile:: Background file missing: '',  a)' ) trim(filename)
    STOP("read_ncoda_intfile:: EXITING...")
  endif

  deallocate(buf4)

END SUBROUTINE read_ncoda_intfile


SUBROUTINE read_grdfiles(lon2d,lat2d) !,lev) !,lev3d)
!===============================================================================
! Read in an ncoda-formatted hycom grid specification file
!===============================================================================
  !STEVE: This subroutine reads the HYCOM/NCODA grid specification file
  USE params_model, ONLY: glon, glat, glev
  USE params_model, ONLY: gridfile_lon, gridfile_lat, gridfile_lev
  USE params_model, ONLY: istart,iend,jstart,jend

  REAL(r_size), DIMENSION(:,:), INTENT(OUT) :: lon2d, lat2d
  CHARACTER(slen) :: filename
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:,:)
  LOGICAL :: ex
  LOGICAL :: dodebug = .true.

  allocate( buf4(nlon,nlat,1) )

  ! Read the global grid specifications (glon x glat)
  ! and convert to a local grid (nlon x nlat)

  !
  ! lon
  !
  filename = gridfile_lon
  if (dodebug) write(6,'(2A)')'read_grdfiles, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,1)
  lon2d(1:nlon,1:nlat) = buf4(:,:,1)

  !
  ! lat
  !
  filename = gridfile_lat
  if (dodebug) write(6,'(2A)')'read_grdfiles, reading file: ',trim(filename)
  buf4 = 0.0
  CALL read_ncoda_bckgrnd(filename,buf4,1)
  lat2d(1:nlon,1:nlat) = buf4(:,:,1)

END SUBROUTINE read_grdfiles


SUBROUTINE read_kmt(kmt)
!===============================================================================
! Read in an ncoda-formatted land/sea mask specifying depth of ocean grid cells
!===============================================================================
  !STEVE: This subroutine reads the HYCOM/NCODA grid specification file
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: gridfile_kmt
  USE params_model, ONLY: istart,iend,jstart,jend

  INTEGER, DIMENSION(:,:), INTENT(OUT) :: kmt
  CHARACTER(slen) :: filename
  INTEGER, ALLOCATABLE :: buf4(:,:,:)
  LOGICAL :: dodebug = .true.

  allocate( buf4(nlon,nlat,1) )

  ! Read the global grid specifications (glon x glat)
  ! and convert to a local grid (nlon x nlat)

  !
  ! land sea mask (kmt)
  !
  filename = gridfile_kmt
  if (dodebug) write(6,*)'read_grdfiles, reading file',trim(filename)
  buf4 = 0
  CALL read_ncoda_intfile(filename,buf4,1)
  kmt(1:nlon,1:nlat) = buf4(:,:,1)

END SUBROUTINE read_kmt


SUBROUTINE write_hycom_ncoda(outfile,v3d,v2d)
!===============================================================================
! Write an ncoda-formatted hycom analysis file
!===============================================================================
  !STEVE: This subroutine write to the HYCOM/NCODA restart file.
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv3d_h
  USE params_model, ONLY: rsrt_tbase,rsrt_sbase,rsrt_ubase,rsrt_vbase,rsrt_sshbase

  CHARACTER(*),INTENT(IN) :: outfile
  REAL(r_sngl),DIMENSION(:,:,:,:),INTENT(IN) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),DIMENSION(:,:,:),  INTENT(IN) :: v2d !(nlon,nlat,nv2d)
  CHARACTER(slen) :: filename
  REAL(r_sngl), ALLOCATABLE :: buf4(:,:,:)
  LOGICAL :: dodebug = .true.

  allocate( buf4(nlon,nlat,nlev) )

  !
  ! Temperature
  !
  filename = trim(outfile)//'.'//trim(rsrt_tbase)
  if (dodebug) write(6,*)'write_hycom_ncoda, writing file: ',trim(filename)
  buf4 = v3d(:,:,:,iv3d_t)
  CALL write_ncoda_analysis(filename,buf4,nlev)

  !
  ! Salinity
  !
  filename = trim(outfile)//'.'//trim(rsrt_sbase)
  if (dodebug) write(6,*)'write_hycom_ncoda, writing file: ',trim(filename)
  buf4 = v3d(:,:,:,iv3d_s)
  CALL write_ncoda_analysis(filename,buf4,nlev)
  
  !
  ! u
  !
  filename = trim(outfile)//'.'//trim(rsrt_ubase)
  if (dodebug) write(6,*)'write_hycom_ncoda, writing file: ',trim(filename)
  buf4 = v3d(:,:,:,iv3d_u)
  CALL write_ncoda_analysis(filename,buf4,nlev)

  !
  !  v
  ! 
  filename = trim(outfile)//'.'//trim(rsrt_vbase)
  if (dodebug) write(6,*)'write_hycom_ncoda, writing file: ',trim(filename)
  buf4 = v3d(:,:,:,iv3d_v)
  CALL read_ncoda_bckgrnd(filename,buf4,nlev)

  !
  ! ssh
  !
  filename = trim(outfile)//'.'//trim(rsrt_sshbase)
  if (dodebug) write(6,*)'write_hycom_ncoda, writing file: ',trim(filename)
  buf4(:,:,1) = v2d(:,:,iv2d_ssh)
  CALL write_ncoda_analysis(filename,buf4(:,:,1:1),1)


END SUBROUTINE write_hycom_ncoda


SUBROUTINE write_ncoda_analysis(filename, outdata, kmax)
  USE params_model, ONLY: glon, glat, glev
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: istart, iend, jstart, jend

  CHARACTER(*), INTENT(IN) :: filename
  REAL(r_sngl), DIMENSION(nlon,nlat,nlev), INTENT(IN) :: outdata
  INTEGER, INTENT(IN) :: kmax
  REAL(r_sngl), DIMENSION(:,:,:), ALLOCATABLE :: buf4
  INTEGER :: is,ie,js,je
  INTEGER, PARAMETER :: fid=20
  LOGICAL   :: ex
  INTEGER*8 :: mx, ny, lz, reclen

!
! global  ..long integer array dimensions
!
  is = istart
  ie = iend
  js = jstart
  je = jend

  mx = glon
  ny = glat
  lz = min(glev,kmax)

  reclen = mx * ny * lz * 4
  allocate(buf4(glon,glat,lz)) !STEVE: TEMPORORAY - replace with stream access write of subgrid

!---------------------------------
! Write the NCODA analysis file
!---------------------------------
  inquire (file=trim(filename), exist=ex)
  if (ex) then
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
!STEVE: (ISSUE) replace with stream I/O to write only desired section and cut down on I/O costs
!   open (unit=fid, file=trim(filename), status='old', access='stream', form='unformatted', recl=reclen)
    read (fid,rec=1) buf4
    buf4(is:ie,js:je,1:lz) = outdata
    write (fid,rec=1) buf4
    close (fid)
    write (6, *) 'write_ncoda_analysis:: Finished writing analysis file: ', trim(filename)
  else
    write (6, *) 'write_ncoda_analysis:: Analysis file missing: ', trim(filename)
    STOP("write_ncoda_analysis:: EXITING...")
  endif

  deallocate(buf4)

END SUBROUTINE write_ncoda_analysis


SUBROUTINE read_blkdat(infile,dims,sigma)
!===============================================================================
! read in the blkdat.input file produced for hycom input to get
! the input parameters and grid definitions to use for letkf.
!===============================================================================

  CHARACTER(*), INTENT(IN) :: infile !='blkdat.input'
  INTEGER, DIMENSION(3), INTENT(OUT) :: dims
  REAL(r_sngl), ALLOCATABLE, INTENT(OUT) :: sigma(:)

! INTEGER, PARAMETER :: slen=512
  INTEGER      :: idm, jdm
  INTEGER      :: itest, jtest
  INTEGER      :: kdm
  INTEGER      :: nhybrd, nsigma
  REAL(r_sngl) :: dp00, dp00x, ds00, ds00x, dp00i
  REAL(r_sngl) :: isotop
  INTEGER      :: thflag
  REAL(r_sngl) :: thbase
  INTEGER      :: vsigma

  INTEGER, PARAMETER :: fid=12
  INTEGER            :: nh=4                    ! nh :: The number of header lines to skip
  CHARACTER(80)    :: dummy
  CHARACTER(70)    :: dummy2
  INTEGER :: i,j,k, ierr

  LOGICAL, PARAMETER :: dodebug = .true.

  !-----------------------------------------------------------------------------
  ! Open the blkdat.input file
  !-----------------------------------------------------------------------------
  OPEN( unit=fid, file=infile, status='old', &
        access='sequential', form='formatted', action='read' )

  !-----------------------------------------------------------------------------
  ! Read the 4 header lines as dummy arguments
  !-----------------------------------------------------------------------------
  do i=1,nh
    READ(fid,'(A)') dummy
    if (dodebug) print *, trim(dummy)
  enddo

  ! Skip 2 lines
  do i=1,2
    READ(fid,'(A)') dummy
    if (dodebug) print *, trim(dummy)
  enddo

! Read the following rows, restarting the count at 1:

  !-----------------------------------------------------------------------------
  ! Row 3: 1500      'idm   ' = longitudinal array size
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) idm
  if (dodebug) print *, "idm = ", idm
  if (ierr>0) then
    print *, "Error on read."
    STOP 3
  elseif (ierr<0) then
    print *, "End of file."
  endif

  dims(1) = idm

  !-----------------------------------------------------------------------------
  ! Row 4: 1100      'jdm   ' = latitudinal  array size
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) jdm
  if (dodebug) print *, "jdm = ", jdm
  if (ierr>0) then
    print *, "Error on read."
    STOP 4
  elseif (ierr<0) then
    print *, "End of file."
  endif

  dims(2) = jdm

  !-----------------------------------------------------------------------------
  ! Row 5   970      'itest ' = grid point where detailed diagnostics are desired
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) itest
  if (dodebug) print *, "itest = ", itest
  if (ierr>0) then
    print *, "Error on read."
    STOP 5
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 6    70      'jtest ' = grid point where detailed diagnostics are desired
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) jtest
  if (dodebug) print *, "jtest = ", jtest
  if (ierr>0) then
    print *, "Error on read."
    STOP 6
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 7    32      'kdm   ' = number of layers
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) kdm
  if (dodebug) print *, "kdm = ", kdm
  if (ierr>0) then
    print *, "Error on read."
    STOP 7
  elseif (ierr<0) then
    print *, "End of file."
  endif

  dims(3) = kdm

  !-----------------------------------------------------------------------------
  ! Row 8    32      'nhybrd' = number of hybrid levels (0=all isopycnal)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) nhybrd
  if (dodebug) print *, "nhybrd = ", nhybrd
  if (ierr>0) then
    print *, "Error on read."
    STOP 8
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 9     0      'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) nsigma
  if (dodebug) print *, "nsigma = ", nsigma
  if (ierr>0) then
    print *, "Error on read."
    STOP 9
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 10   3.0     'dp00  ' = deep    z-level spacing minimum thickness (m)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) dp00
  if (dodebug) print *, "dp00 = ", dp00
  if (ierr>0) then
    print *, "Error on read."
    STOP 10
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 11 450.0     'dp00x ' = deep    z-level spacing maximum thickness (m)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) dp00x
  if (dodebug) print *, "dp00x = ", dp00x
  if (ierr>0) then
    print *, "Error on read."
    STOP 11
  elseif (ierr<0) then
    print *, "End of file."
  endif

  ! Skip 1 line
  do i=1,1
    READ(fid,'(A)') dummy
    if (dodebug) print *, trim(dummy)
  enddo

  !-----------------------------------------------------------------------------
  ! Row 13   3.0    'ds00  ' = shallow z-level spacing minimum thickness (m)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) ds00
  if (dodebug) print *, "ds00 = ", ds00
  if (ierr>0) then
    print *, "Error on read."
    STOP 11
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 14 450.0    'ds00x ' = shallow z-level spacing maximum thickness (m)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) ds00x
  if (dodebug) print *, "ds00x = ", ds00x
  if (ierr>0) then
    print *, "Error on read."
    STOP 11
  elseif (ierr<0) then
    print *, "End of file."
  endif


  ! Skip 1 line
  do i=1,1
    READ(fid,'(A)') dummy
    if (dodebug) print *, trim(dummy)
  enddo

  !-----------------------------------------------------------------------------
  ! Row 16   1.0    'dp00i ' = deep iso-pycnal spacing minimum thickness (m)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) dp00i
  if (dodebug) print *, "dp00i = ", dp00i
  if (ierr>0) then
    print *, "Error on read."
    STOP 11
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 17   6.0    'isotop' = shallowest depth for isopycnal layers (m), <0 from file
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) isotop
  if (dodebug) print *, "isotop = ", isotop
  if (ierr>0) then
    print *, "Error on read."
    STOP 11
  elseif (ierr<0) then
    print *, "End of file."
  endif

  ! Skip 3 lines
  do i=1,3
    READ(fid,'(A)') dummy
    if (dodebug) print *, trim(dummy)
  enddo

  !-----------------------------------------------------------------------------
  ! Row 21   2      'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) thflag
  if (dodebug) print *, "thflag = ", thflag
  if (ierr>0) then
    print *, "Error on read."
    STOP 21
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 22  34.0    'thbase' = reference density (sigma units)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,F6.1)', iostat=ierr) thbase
  if (dodebug) print *, "thbase = ", thbase
  if (ierr>0) then
    print *, "Error on read."
    STOP 22
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 23   1      'vsigma' = spacially varying isopycnal target densities (0=F,1=T)
  !-----------------------------------------------------------------------------
  read(fid,'(T1,I4)', iostat=ierr) vsigma
  if (dodebug) print *, "vsigma = ", vsigma
  if (ierr>0) then
    print *, "Error on read."
    STOP 23
  elseif (ierr<0) then
    print *, "End of file."
  endif

  !-----------------------------------------------------------------------------
  ! Row 24  28.10   'sigma ' = layer  1 isopycnal target density (sigma units)
  ! Rows 25 to 23+kdm :: layers 2-32
  !-----------------------------------------------------------------------------
  ALLOCATE(sigma(kdm))

  do i=1,kdm
    read(fid,'(T1,F7.2)', iostat=ierr) sigma(i)
    if (dodebug) print *, "sigma(i) = ", i, sigma(i)
    if (ierr>0) then
      print *, "Error on read."
      STOP 24
    elseif (ierr<0) then
      print *, "End of file."
    endif
  enddo

END SUBROUTINE read_blkdat


END MODULE hycom_io
