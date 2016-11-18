MODULE hycom_io
! INTEGER,PARAMETER :: r_size=kind(0.0d0)
! INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
! INTEGER, PARAMETER :: slen=240
! INTEGER, PARAMETER :: nlon=1500,nlat=1100,nlev=32
! INTEGER, PARAMETER :: nv2d=3,nv3d=5
! REAL(r_size)   :: v2d(nlon,nlat,nv2d)
! REAL(r_size)   :: v3d(nlon,nlat,nlev,nv3d)
! INTEGER        :: k,nn,nv
! CHARACTER(slen) :: file_in != 'archv.2009_335_01_fsd.bin'

  USE common, ONLY: r_size, r_sngl, slen
  USE params_model

  IMPLICIT NONE

  PUBLIC

  REAL, PARAMETER :: hycom_undef=1e10
  
CONTAINS 

SUBROUTINE read_hycom(file_in,v3d,v2d)
!===============================================================================
! Read in a (big-endian) binary file extracted from hycom ab-format
!===============================================================================
  CHARACTER(slen), INTENT(IN) :: file_in
  REAL(r_size), INTENT(OUT)   :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size), INTENT(OUT)   :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl)                :: buf(nlon,nlat)
  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12
  INTEGER, PARAMETER, DIMENSION(nv2d) :: input_order_2d = (/ 1,2,3 /), output_order_2d = (/ 1,2,3 /)
  INTEGER, PARAMETER, DIMENSION(nv3d) :: input_order_3d = (/ 1,2,3,4,5 /), output_order_3d = (/ 1,2,5,3,4/)

  !***************************************************
  ! INPUT:
  !    v2d(i,j,1) = SSH
  !    v2d(i,j,2) = u_barotropic
  !    v2d(i,j,3) = v_barotropic
  !    v3d(i,j,k,1) = u_baroclinic (k=1 is top layer)
  !    v3d(i,j,k,2) = v_baroclinic  
  !    v2d(i,j,k,3) = layer thickness
  !    v2d(i,j,k,4) = temperature
  !    v2d(i,j,k,5) = salinity
  ! OUTPUT Changed to:
  !    v2d(i,j,1) = SSH
  !    v2d(i,j,2) = u_barotropic
  !    v2d(i,j,3) = v_barotropic
  !    v3d(i,j,k,1) = u_baroclinic (k=1 is top layer)
  !    v3d(i,j,k,2) = v_baroclinic  
  !    v2d(i,j,k,3) = temperature
  !    v2d(i,j,k,4) = salinity
  !    v2d(i,j,k,5) = layer thickness
  !***************************************************

  ! Direct Access:
  if (hycom_io_access==0) then
    open(fid,file=file_in,form='unformatted',access='direct',recl=nlon*nlat)

    do n=1,nv2d
      read(fid,rec=n) buf 
      v2d(:,:,output_order_2d(n)) = buf
    enddo

    do n=1,nv3d
      do k=1,nlev
        nrec=nv2d+(n-1)*nlev+k
        read(fid,rec=nrec) buf
        v3d(:,:,k,output_order_3d(n)) = buf
      enddo
    enddo

  ! Sequential Access:
  elseif (hycom_io_access==1) then
    open(fid,file=file_in,form='unformatted',access='sequential')

    do n=1,nv2d
      read(fid) v2d(:,:,output_order_2d(n))
    enddo

    do n=1,nv3d
      do k=1,nlev
        read(fid) v3d(:,:,k,output_order_3d(n))
      enddo
    enddo

  else
    WRITE(6,*) "read_hycom:: access not supported :: hycom_io_access = ", hycom_io_access
    STOP(894)
  endif

END SUBROUTINE read_hycom

SUBROUTINE write_hycom(file_out,v3d,v2d)
!===============================================================================
! Read in a (big-endian) binary file extracted from hycom ab-format
!===============================================================================
  CHARACTER(slen), INTENT(IN) :: file_out
  REAL(r_sngl), INTENT(IN)    :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), INTENT(IN)    :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl)                :: buf(nlon,nlat)
  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12
  INTEGER, PARAMETER, DIMENSION(nv2d) :: input_order_2d = (/ 1,2,3 /), output_order_2d = (/ 1,2,3 /)
  INTEGER, PARAMETER, DIMENSION(nv3d) :: input_order_3d = (/ 1,2,3,4,5 /), output_order_3d = (/ 1,2,5,3,4/)
  LOGICAL :: dodebug = .false.

  !***************************************************
  ! LETKF representation:
  !    v2d(i,j,1) = SSH
  !    v2d(i,j,2) = u_barotropic
  !    v2d(i,j,3) = v_barotropic
  !    v3d(i,j,k,1) = u_baroclinic (k=1 is top layer)
  !    v3d(i,j,k,2) = v_baroclinic  
  !    v2d(i,j,k,3) = temperature
  !    v2d(i,j,k,4) = salinity
  !    v2d(i,j,k,5) = layer thickness
  ! OUTPUT Changed to:
  !    v2d(i,j,1) = SSH
  !    v2d(i,j,2) = u_barotropic
  !    v2d(i,j,3) = v_barotropic
  !    v3d(i,j,k,1) = u_baroclinic (k=1 is top layer)
  !    v3d(i,j,k,2) = v_baroclinic  
  !    v2d(i,j,k,3) = layer thickness
  !    v2d(i,j,k,4) = temperature
  !    v2d(i,j,k,5) = salinity
  !***************************************************


  ! Direct Access:
  if (hycom_io_access==0) then
    open(fid,file=file_out,form='unformatted',access='direct',recl=nlon*nlat)

    do n=1,nv2d
      if (dodebug) WRITE(6,*) n
      if (dodebug) then
        WRITE(6,*) "output_order_2d(n) = ", output_order_2d(n)
        WRITE(6,*) "nlon = ", nlon
        WRITE(6,*) "nlat = ", nlat
        WRITE(6,*) "nv2d = ", nv2d
        WRITE(6,*) "v2d(1,1,output_order_2d(n)) = ", v2d(1,1,output_order_2d(n))
      endif
      buf = v2d(:,:,output_order_2d(n))
      write(fid,rec=n) buf
    enddo

    do n=1,nv3d
      do k=1,nlev
        nrec=nv2d+(n-1)*nlev+k
        if (dodebug) WRITE(6,*) nrec
        if (dodebug) then
          WRITE(6,*) "output_order_3d(n) = ", output_order_3d(n)
          WRITE(6,*) "v3d(1,1,k,output_order_3d(n)) = ", v3d(1,1,k,output_order_3d(n))
        endif
        buf = v3d(:,:,k,output_order_3d(n))
        write(fid,rec=nrec) buf
      enddo
    enddo

  ! Sequential Access:
  elseif (hycom_io_access==1) then
    open(fid,file=file_out,form='unformatted',access='sequential')

    do n=1,nv2d
      write(fid) v2d(:,:,output_order_2d(n))
    enddo

    do n=1,nv3d
      do k=1,nlev
        write(fid) v3d(:,:,k,output_order_3d(n))
      enddo
    enddo

  else
    WRITE(6,*) "write_hycom:: access not supported :: hycom_io_access = ", hycom_io_access
    STOP(892)
  endif

END SUBROUTINE write_hycom

SUBROUTINE read_blkdat(infile,sigma)
!===============================================================================
! read in the blkdat.input file produced for hycom input to get
! the input parameters and grid definitions to use for letkf.
!===============================================================================

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: infile !='blkdat.input'
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
