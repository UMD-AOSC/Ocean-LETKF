PROGRAM check_for_ocean
! This program reads in land-sea mask and checks for ocean points
! present in a subgrid

  IMPLICIT NONE

  INTEGER, PARAMETER :: slen = 64
  ! Default inputs:
  INTEGER :: ntiles = 308
  INTEGER :: glon = 4500
  INTEGER :: glat = 3298
  INTEGER :: glev = 41
  CHARACTER(slen) :: infile = 'kmt.dat'
  CHARACTER(slen) :: outfile = 'no_ocean.skip'
  !
  INTEGER :: nlon, nlat, nlev
  INTEGER :: i,j,k
  INTEGER :: is=0,ie=0,js=0,je=0
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: buf4
  INTEGER :: fid = 21
  LOGICAL :: ex
  LOGICAL :: dodebug = .true.
  INTEGER :: kmax = 1
  INTEGER :: reclen_mult = 4

  !-----------------------------------------------------------------------------
  ! Read command line options to override defaults
  !-----------------------------------------------------------------------------
  CALL process_command_line

  !-----------------------------------------------------------------------------
  ! Read tiled data and check for positive values
  !-----------------------------------------------------------------------------
  nlon = ie - is + 1
  nlat = je - js + 1
  nlev = 1
  allocate(buf4(nlon,nlat,nlev))

  ! Read subgrid
  CALL read_tile_int(infile,buf4,kmax)

  ! Print the max and min of the local tile:
  print *, "Max(buf4) = ", MAXVAL(buf4) 
  print *, "Min(buf4) = ", MINVAL(buf4) 
  print *, "Avg(buf4) = ", SUM(buf4)/(nlon*nlat*nlev)

  if (MAXVAL(buf4)<=0) then 
    !-----------------------------------------------------------------------------
    ! Create emtpy file flagging 'no ocean points'
    !-----------------------------------------------------------------------------
    OPEN(fid,file=outfile,status='new')
    CLOSE(fid)
  else
    print *, "max = ", MAXVAL(buf4)
    if (dodebug) print *, "buf4 = ", buf4
  endif

CONTAINS

!-------------------------------------------------------------------------------
! Read tiled data
!-------------------------------------------------------------------------------
SUBROUTINE read_tile_int(filename,buf,kmax)
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: buf
  INTEGER, INTENT(IN) :: kmax
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: buf_global
  INTEGER*8 :: reclen, iolen
  INTEGER :: mx,ny,lz, nlon,nlat,nlev
  LOGICAL :: ex
  INTEGER :: fid = 22
  LOGICAL :: dodebug = .true.

  nlon = SIZE(buf,1)
  nlat = SIZE(buf,2)
  nlev = SIZE(buf,3)

  mx = glon
  ny = glat
  lz = min(glev,nlev)
  lz = min(lz,kmax)

  allocate(buf_global(mx,ny,lz))
  reclen = mx * ny * lz * reclen_mult

  if (dodebug) then
    write(6,*) "nlon   = ", nlon
    write(6,*) "nlat   = ", nlat
    write(6,*) "nlev   = ", nlev
    write(6,*) "mx     = ", mx
    write(6,*) "ny     = ", ny
    write(6,*) "lz     = ", lz
    write(6,*) "_mult  = ", reclen_mult
    write(6,*) "mx * ny * lz * reclen_mult = ", mx * ny * lz * reclen_mult
    write(6,*) "reclen = ", reclen
    write(6,*) "iolen  = ", iolen
  endif

  inquire (file=trim(filename), exist=ex)
  if (ex) then
    if (dodebug) write(6,*) "read_tile_int:: Opening file: ", trim(filename)
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    read (fid,rec=1) buf_global
    close (fid)
    if (dodebug) then
      print *, "Max(buf_global) = ", MAXVAL(buf_global)
      print *, "Min(buf_global) = ", MINVAL(buf_global)
      print *, "Avg(buf_global) = ", SUM(buf_global)/(mx*ny*lz)
      print *, "Maxloc(buf_global) = ", MAXLOC(buf_global)
      print *, "Minloc(buf_global) = ", MINLOC(buf_global)
      print *, "is = ", is
      print *, "ie = ", ie
      print *, "js = ", js
      print *, "je = ", je
      print *, "lz = ", lz
      print *, "Max(buf_global(is:ie,js:je,1:lz)) = ", MAXVAL(buf_global(is:ie,js:je,1:lz))
      print *, "Min(buf_global(is:ie,js:je,1:lz)) = ", MINVAL(buf_global(is:ie,js:je,1:lz))
      print *, "Avg(buf_global(is:ie,js:je,1:lz)) = ", SUM(buf_global(is:ie,js:je,1:lz))/(nlon*nlat*nlev)
    endif
    buf(1:nlon,1:nlat,1:nlev) = buf_global(is:ie,js:je,1:lz)
    if (dodebug) write (6,*) 'read_tile_int:: Finished reading tile file: ', trim(filename)
  else
    write (6, *) 'read_tile_int:: tile file missing: ', trim(filename)
    STOP("merge_tiles.f90:: EXITING...")
  endif

  deallocate(buf_global)

END SUBROUTINE read_tile_int


SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In obsop_tprof.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)
  select case (arg1)
    case('-f')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      infile = arg2
    case('-o')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      outfile = arg2
    case('-is')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) is
    case('-ie')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) ie
    case('-js')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) js
    case('-je')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) je
    case('-kmax')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) kmax
    case('-ntiles')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) ntiles
    case('-glon')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) glon
    case('-glat')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) glat
    case('-glev')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) glev
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line


END PROGRAM check_for_ocean
