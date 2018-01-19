program test_read_bin
!
! Sample compile line:
! pgf90 -byteswapio test_read_bin.f90 -o trb.x
!
! Matching NRL optionsi for LETKF:
! pgf90 -fast -byteswapio -Mextend -tp=nehalem,sandybridge -r4 -i4 -Mpreprocess test_read_bin.f90 -o trb.x
!

  IMPLICIT NONE

  INTEGER :: iolen, mx, ny, lz, reclen
  CHARACTER(256) :: filename
! REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: indata
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: indata
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: buf4
  INTEGER :: fid = 21
  INTEGER :: is,ie,js,je
  INTEGER :: nlon,nlat,nlev
  LOGICAL :: ex
  LOGICAL :: dodebug = .true.

  ! File details:
! filename = 'gues.t.dat'
  filename = 'anal.t.dat'
  mx = 187
  ny = 121
  lz = 32
  is = 85
  ie = 105
  js = 50
  je = 70
  nlon=ie-is+1
  nlat=je-js+1
  nlev=lz

  ! Allocate memory:
  allocate(buf4(mx,ny,lz))
  allocate(indata(nlon,nlat,nlev))

  ! Inquire io details
  INQUIRE(IOLENGTH=iolen) iolen
  reclen = mx * ny * lz * iolen ! * 4

  ! Debug 
  if (dodebug) then
    WRITE(6,*) "mx = ", mx
    WRITE(6,*) "ny = ", ny
    WRITE(6,*) "lz = ", lz
    WRITE(6,*) "reclen = ", reclen
  endif

  !---------------------------------
  ! Read the NCODA background file
  !---------------------------------
  inquire (file=trim(filename), exist=ex)
  if (ex) then
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    read (fid,rec=1) buf4
    close (fid)
    indata = buf4(is:ie,js:je,1:lz)
    WRITE(6,*) "Finished assigning buf4 to indata."
  else
    STOP("input file does not exists. EXITING...")
  endif

  if (dodebug) then
    write (6, *) "buf4(m/2,y/2,1)= ", buf4(mx/2,ny/2,1)
    write (6, *) " max/min buf4 = ",  maxval(buf4), minval(buf4)
    write (6, *) "indata(m/2,y/2,1)= ", indata(nlon/2,nlat/2,1)
    write (6, *) " max/min indata = ",  maxval(indata), minval(indata)
  endif


end program test_read_bin
