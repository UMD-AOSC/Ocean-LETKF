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
  USE mod_za

  IMPLICIT NONE

  PUBLIC

  REAL(r_size), PARAMETER :: hycom_undef=2.0**100
  REAL(r_size), PARAMETER :: hycom_eps = EPSILON(1.0)  !(intentionally at real precision)              
  
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

SUBROUTINE get_hycom(file_in_a,file_in_b,v3d,v2d)

!
! --- hycom/micom field extractor
!
! Added by Jili Dong Jan. 2017

  CHARACTER(slen), INTENT(IN) :: file_in_a,file_in_b
  REAL(r_size), INTENT(OUT)   :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size), INTENT(OUT)   :: v2d(nlon,nlat,nv2d)
  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12

  character cline*80
  character ctitle(4)*80
  real :: dummy_2d(nlon,nlat)
  integer :: jversn,iexpt,yrflag
  INTEGER :: MSK(nlon,nlat)
  real :: density(nlon,nlat,nlev)
  real :: HMINA,HMAXA,hminb,hmaxb,thbase
  integer :: i,nstep,nz,layer
  double precision :: time(3)

        call XCSPMD
        ! JILI skip ZAIOST
        !CALL ZAIOST

! Please make sure the variables in "a" file are in the right order
! Check "b" file to make sure
 
! open archive file
        CALL ZAIOPF(file_in_a,'OLD',21)
        open(110,file=file_in_b,form='formatted',&
         &status='old',action='read')

        read(110,116) ctitle,jversn,iexpt,yrflag,idm,jdm
116    format (a80/a80/a80/a80/&
        &i5,4x,'''iversn'' = hycom version number x10'/&
        &i5,4x,'''iexpt '' = experiment number x10'/&
        &i5,4x,'''yrflag'' = days in year flag'/&
        &i5,4x,'''idm   '' = longitudinal array size'/&
        &i5,4x,'''jdm   '' = latitudinal  array size'/&
        &'field       time step  model day',&
        &'  k  dens        min              max')

! read montg1
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read srfhgt
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("srfhgt  ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v2d(:,:,iv2d_ssh)=dummy_2d
        where (abs(dummy_2d-hycom_undef) > hycom_eps) v2d(:,:,iv2d_ssh)=dummy_2d/9.806  

! read steric
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read surflx
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read salflx
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read bl_dpth
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read mix_dpth
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read tmix
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read smix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read thmix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read umix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read vmix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read u_batro
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("u_btrop ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v2d(:,:,iv2d_ubt)=dummy_2d

! read v_batro
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("v_btrop ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v2d(:,:,iv2d_vbt)=dummy_2d


!***************************************************
!    v2d(i,j,1) = SSH
!    v2d(i,j,2) = x_barotropic velocity  (not eastward) 
!    v2d(i,j,3) = y_barotropic velocity (not northward)
!    u(i,j,k) = x_baroclinic velocity (k=1 is top layer)
!    v(i,j,k) = y_baroclinic velocity  
!    dp(i,j,k) = layer thickness
!    temp(i,j,k) = temperature
!    saln(i,j,k) = salinity
!***************************************************
        do k=1,nlev
! total u
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("u-vel.  ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
!        v3d(:,:,k,iv3d_u)=dummy_2d+v2d(:,:,iv2d_ubt)
! For now, use u perturbation
        v3d(:,:,k,iv3d_u)=dummy_2d

! total v   
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("v-vel.  ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
!        v3d(:,:,k,iv3d_v)=dummy_2d+v2d(:,:,iv2d_vbt)
! For now, use v perturbation
        v3d(:,:,k,iv3d_v)=dummy_2d



! dp
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("thknss  ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v3d(:,:,k,iv3d_h)=dummy_2d
        where (abs(dummy_2d-hycom_undef) > hycom_eps) v3d(:,:,k,iv3d_h)=dummy_2d/9806.                  

! temp    
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("temp    ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v3d(:,:,k,iv3d_t)=dummy_2d


! saln
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("salin   ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        v3d(:,:,k,iv3d_s)=dummy_2d


! dens 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        call check_ab("density ",cline(1:8),hmina,hminb,hmaxa,hmaxb)
        density(:,:,k)=dummy_2d

        enddo


      call zaiocl(21)
      close (110)


END SUBROUTINE get_hycom 

SUBROUTINE get_hycom_sample(sample_3d,pmsk,umsk,vmsk)

!
! --- hycom/micom field extractor
!
! Added by Jili Dong Jan. 2017

  REAL(r_size), INTENT(OUT)   :: sample_3d(nlon,nlat,nlev)
  INTEGER,INTENT(OUT)  :: pmsk(nlon,nlat),umsk(nlon,nlat),vmsk(nlon,nlat)
  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12

  character cline*80
  character ctitle(4)*80
  real :: dummy_2d(nlon,nlat)
  integer :: jversn,iexpt,yrflag
  INTEGER :: MSK(nlon,nlat)
  real :: HMINA,HMAXA,hminb,hmaxb,thbase
  integer :: i,nstep,nz,layer,ni,nj
  double precision :: time(3)

        call XCSPMD
! JILI skip zaiost
!        CALL ZAIOST
! Please make sure the variables in "a" file are in the right order
! Check "b" file to make sure

pmsk=1
umsk=1
vmsk=1

! open archive file
        CALL ZAIOPF("hycom_sample.a",'OLD',21)
        open(110,file="hycom_sample.b",form='formatted',&
         &status='old',action='read')

        read(110,116) ctitle,jversn,iexpt,yrflag,idm,jdm
116    format (a80/a80/a80/a80/&
        &i5,4x,'''iversn'' = hycom version number x10'/&
        &i5,4x,'''iexpt '' = experiment number x10'/&
        &i5,4x,'''yrflag'' = days in year flag'/&
        &i5,4x,'''idm   '' = longitudinal array size'/&
        &i5,4x,'''jdm   '' = latitudinal  array size'/&
        &'field       time step  model day',&
        &'  k  dens        min              max')

! read montg1
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read srfhgt
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read steric
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read surflx
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read salflx
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read bl_dpth
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read mix_dpth
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read tmix
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        do ni=1,nlon
          do nj=1,nlat
            if (dummy_2d(ni,nj) > 100000.0) pmsk(ni,nj)=0
          enddo
        enddo

! read smix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read thmix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! read umix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        do ni=1,nlon
          do nj=1,nlat
            if (dummy_2d(ni,nj) > 100000.0) umsk(ni,nj)=0
          enddo
        enddo


! read vmix 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        do ni=1,nlon
          do nj=1,nlat
            if (dummy_2d(ni,nj) > 100000.0) vmsk(ni,nj)=0
          enddo
        enddo



! read u_batro
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb


! read v_batro
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

!***************************************************
!    v2d(i,j,1) = SSH
!    v2d(i,j,2) = x_barotropic velocity  (not eastward) 
!    v2d(i,j,3) = y_barotropic velocity (not northward)
!    u(i,j,k) = x_baroclinic velocity (k=1 is top layer)
!    v(i,j,k) = y_baroclinic velocity  
!    dp(i,j,k) = layer thickness
!    temp(i,j,k) = temperature
!    saln(i,j,k) = salinity
!***************************************************
        do k=1,nlev
! total u
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! total v   
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! dp
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb


! temp    
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        sample_3d(:,:,k)=dummy_2d


! saln
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

! dens 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        enddo


      call zaiocl(21)
      close (110)


END SUBROUTINE get_hycom_sample



SUBROUTINE get_hycom_depth(phi0,lat2d,lon2d,dx,dy)

!
! --- hycom/micom grid depth field extractor
!
! Added by Jili Dong Jan. 2017

  REAL(r_size), INTENT(OUT)   :: phi0(nlon,nlat)
  REAL(r_size), INTENT(OUT)   :: lon2d(nlon,nlat)
  REAL(r_size), INTENT(OUT)   :: lat2d(nlon,nlat)
  REAL(r_size), INTENT(OUT)   :: dx(nlon,nlat)
  REAL(r_size), INTENT(OUT)   :: dy(nlon,nlat)
!  REAL(r_size), INTENT(OUT)   :: ulon2d(nlon,nlat)
!  REAL(r_size), INTENT(OUT)   :: ulat2d(nlon,nlat)
!  REAL(r_size), INTENT(OUT)   :: vlon2d(nlon,nlat)
!  REAL(r_size), INTENT(OUT)   :: vlat2d(nlon,nlat)


  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12

  character cline*80
  character ctitle(4)*80
  real :: dummy_2d(nlon,nlat)
  INTEGER :: MSK(nlon,nlat)
  real :: HMINA,HMAXA,hminb,hmaxb,thbase
  integer :: i
  character preambl(5)*79


        call XCSPMD
        CALL ZAIOST

! Please make sure the variables in "a" file are in the right order
! Check "b" file to make sure

      open (9,file="regional.depth.b",   &
           form='formatted',status='old',action='read')
      read (9, '(a79)') preambl
      read (9, '(a)')   cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(9)


      call zaiopf("regional.depth.a",'old', 9)
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      call zaiocl(9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),   &
                                     abs(hminb))*1.e-4 .or.   &
             abs(hmaxa-hmaxb).gt.max(abs(hmaxa),   &
                                     abs(hmaxb))*1.e-4     ) then
        write(*,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')   &
         'depth error - .a and .b files not consistent:',   &
         '.a,.b min = ',hmina,hminb,hmina-hminb,      &
         '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif

      phi0 = dummy_2d

!
!     grid location.
!
      open (unit=9,file='regional.grid.b',   &
           form='formatted',status='old',action='read')

      read(9,*) ! skip idm
      read(9,*) ! skip jdm
      read(9,*) ! skip mapflg

! read plon
      call zaiopf('regional.grid.a','old', 9)
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call check_grid_ab("plon",cline(1:4),hmina,hminb,hmaxa,hmaxb)
      lon2d=dummy_2d


! read plat
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call check_grid_ab("plat",cline(1:4),hmina,hminb,hmaxa,hmaxb)
      lat2d=dummy_2d

! skip qlon
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb


! skip qlat
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb


! skip ulon
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb


! skip ulat
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb


! skip vlon
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb


! skip vlat
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb

! skip pang 
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb

! read dx (pscx) 
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call check_grid_ab("pscx",cline(1:4),hmina,hminb,hmaxa,hmaxb)
      dx=dummy_2d

! read dy (pscy) 
      call zaiord(dummy_2d,MSK,.false., hmina,hmaxa, 9)
      read(  9,'(a)') cline
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call check_grid_ab("pscy",cline(1:4),hmina,hminb,hmaxa,hmaxb)
      dy=dummy_2d

      close(9)
      call zaiocl(9)

end subroutine get_hycom_depth 



subroutine check_ab(var_name_a,var_name_b,hmina,hminb,hmaxa,hmaxb)

! Check HYCOM a and b files consistency
! Added by Jili Dong Jan. 2017

real :: hmina,hminb,hmaxa,hmaxb
character(8) :: var_name_a,var_name_b 

!write(6,*),"read ",var_name_b," as ",var_name_a

if (var_name_a .ne. var_name_b) then
   write(6,*),var_name_a," not consistent with ",var_name_b 
   stop
endif

        if     (abs(hmina-hminb).gt.max(abs(hmina),     &
                                       abs(hminb))*1.e-4 .or.   &
               abs(hmaxa-hmaxb).gt.max(abs(hmaxa),     &
                                       abs(hmaxb))*1.e-4     ) then
          write(*,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')     &
           'error - .a and .b files not consistent:',      &
           '.a,.b min = ',hmina,hminb,hmina-hminb,         &
           '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif

end subroutine check_ab

subroutine check_grid_ab(var_name_a,var_name_b,hmina,hminb,hmaxa,hmaxb)

! Check HYCOM a and b files consistency
! Added by JILI DONG Jan. 2017


real :: hmina,hminb,hmaxa,hmaxb
character(4) :: var_name_a,var_name_b


write(6,*),"read ",var_name_b," as ",var_name_a

if (var_name_a .ne. var_name_b) then
   write(6,*),var_name_a," not consistent with ",var_name_b
   stop
endif

        if     (abs(hmina-hminb).gt.max(abs(hmina),     &
                                       abs(hminb))*1.e-4 .or.   &
               abs(hmaxa-hmaxb).gt.max(abs(hmaxa),     &
                                       abs(hmaxb))*1.e-4     ) then
          write(*,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')     &
           'grid error - .a and .b files not consistent:',      &
           '.a,.b min = ',hmina,hminb,hmina-hminb,         &
           '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif

end subroutine check_grid_ab


SUBROUTINE put_hycom(file_in_a,file_in_b,v3d,v2d)

!
! --- hycom output 
!
! Added by Jili Dong Apr. 2017
! Open an old archive as a template and write out updated T/S/U/V

  CHARACTER(slen), INTENT(IN) :: file_in_a,file_in_b
  REAL(r_sngl),DIMENSION(:,:,:,:) :: v3d !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),DIMENSION(:,:,:) :: v2d !(nlon,nlat,nv2d)
  INTEGER                     :: k,n,nrec
  INTEGER, PARAMETER          :: fid=12
  character cline*80
  character ctitle(4)*80
  real :: dummy_2d(nlon,nlat)
  integer :: jversn,iexpt,yrflag
  real(r_size) :: hycom_depth(nlon,nlat)
  INTEGER :: mskp(nlon,nlat),msku(nlon,nlat),mskv(nlon,nlat)
  INTEGER :: MSK(nlon,nlat)
  real :: HMINA,HMAXA,hminb,hmaxb,thbase
  real :: xmin, xmax
  integer :: i,i1,j1,nstep,nz,layer
  double precision :: time(3)

        call XCSPMD
        ! JILI skip ZAIOST
        !CALL ZAIOST

mskp=0
msku=0
mskv=0


! Please make sure the variables in "a" file are in the right order
! Check "b" file to make sure

! open old archive file
        CALL ZAIOPF(file_in_a,'OLD',21)
        open(110,file=file_in_b,form='formatted',&
         &status='old',action='read')

! open new archive file
        open (unit=24,file='anal'//file_in_b(5:7)//'.b',form='formatted',  &
     &          status='new',action='write')
        call zaiopf('anal'//file_in_b(5:7)//'.a','new', 24)

        read(110,116) ctitle,jversn,iexpt,yrflag,idm,jdm
116    format (a80/a80/a80/a80/&
        &i5,4x,'''iversn'' = hycom version number x10'/&
        &i5,4x,'''iexpt '' = experiment number x10'/&
        &i5,4x,'''yrflag'' = days in year flag'/&
        &i5,4x,'''idm   '' = longitudinal array size'/&
        &i5,4x,'''jdm   '' = latitudinal  array size'/&
        &'field       time step  model day',&
        &'  k  dens        min              max')

! write header
        write(24,116) ctitle,jversn,iexpt,yrflag,idm,jdm

! read montg1
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        where (abs(dummy_2d-hycom_undef) > hycom_eps) mskp=1


        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'montg1  ',nstep,time(1),layer,thbase,xmin,xmax

! read srfhgt 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        v2d(:,:,iv2d_ssh)=v2d(:,:,iv2d_ssh)*9.806

        CALL ZAIOWR(v2d(:,:,iv2d_ssh),mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'srfhgt  ',nstep,time(1),layer,thbase,xmin,xmax

! read steric 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'steric  ',nstep,time(1),layer,thbase,xmin,xmax

! read surflx 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'surflx  ',nstep,time(1),layer,thbase,xmin,xmax

! read salflx 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'salflx  ',nstep,time(1),layer,thbase,xmin,xmax

! read bl_dpth 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'bl_dpth ',nstep,time(1),layer,thbase,xmin,xmax

! read mix_dpth 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'mix_dpth',nstep,time(1),layer,thbase,xmin,xmax

! read tmix 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'tmix    ',nstep,time(1),layer,thbase,xmin,xmax

! read smix 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'smix    ',nstep,time(1),layer,thbase,xmin,xmax

! read thmix 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'thmix   ',nstep,time(1),layer,thbase,xmin,xmax

! read umix 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        where (abs(dummy_2d-hycom_undef) > hycom_eps) msku=1

        CALL ZAIOWR(dummy_2d,msku,.true.,xmin,xmax,24,.false.)
        write(24,117) 'umix    ',nstep,time(1),layer,thbase,xmin,xmax

! read vmix 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb
        where (abs(dummy_2d-hycom_undef) > hycom_eps) mskv=1

        CALL ZAIOWR(dummy_2d,mskv,.true.,xmin,xmax,24,.false.)
        write(24,117) 'vmix    ',nstep,time(1),layer,thbase,xmin,xmax

! read u_batro 
        CALL ZAIORD(dummy_2d,MSK,.false.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v2d(:,:,iv2d_ubt),msku,.true.,xmin,xmax,24,.false.)
        write(24,117) 'u_btrop ',nstep,time(1),layer,thbase,xmin,xmax

! read v_batro 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v2d(:,:,iv2d_vbt),mskv,.true.,xmin,xmax,24,.false.)
        write(24,117) 'v_btrop ',nstep,time(1),layer,thbase,xmin,xmax

! 3d variables

        do k=1,nlev

! u pert
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v3d(:,:,k,iv3d_u),msku,.true.,xmin,xmax,24,.false.)
        write(24,117) 'u-vel.  ',nstep,time(1),layer,thbase,xmin,xmax

! v pert
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v3d(:,:,k,iv3d_v),mskv,.true.,xmin,xmax,24,.false.)
        write(24,117) 'v-vel.  ',nstep,time(1),layer,thbase,xmin,xmax


! thickness (not updated for now) 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'thknss  ',nstep,time(1),layer,thbase,xmin,xmax

! temp 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v3d(:,:,k,iv3d_t),mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'temp    ',nstep,time(1),layer,thbase,xmin,xmax

! saln 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(v3d(:,:,k,iv3d_s),mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'salin    ',nstep,time(1),layer,thbase,xmin,xmax

! density 
        CALL ZAIORD(dummy_2d,MSK,.FALSE.,HMINA,HMAXA,21)
        read(110,'(a)') cline
        i= index(cline,'=')
        read (cline(i+1:),*) nstep,time(1),layer,thbase,hminb,hmaxb

        CALL ZAIOWR(dummy_2d,mskp,.true.,xmin,xmax,24,.false.)
        write(24,117) 'density ',nstep,time(1),layer,thbase,xmin,xmax

        enddo

 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)

      call zaiocl(21)
      close (110)
      call zaiocl(24)
      close (24)

end subroutine put_hycom 


END MODULE hycom_io
