PROGRAM merge_tiles
! This program reads in tiled subgrid NCODA-HYCOM binary files 
! and writes to disk a single global binary file

  IMPLICIT NONE

  INTEGER, PARAMETER :: slen = 64
  ! Default inputs:
  INTEGER :: ntiles = 308
  INTEGER :: glon = 4500
  INTEGER :: glat = 3298
  INTEGER :: glev = 41
  CHARACTER(slen) :: infile = 'gs01001.t.dat'
  CHARACTER(slen) :: outfile = 'anal001.t.dat'
  CHARACTER(slen) :: filelist_file = 'filelist.txt'
  CHARACTER(slen) :: tiledef_file = 'tiledef.txt'
  CHARACTER(slen) :: lsmask_file = 'kmt.dat'
  LOGICAL :: DO_INC=.false.
  !
  INTEGER :: nlon, nlat, nlev
  INTEGER :: i,j,k
  INTEGER :: is,ie,js,je
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: buf4, buf4_global
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: bufint_global
  CHARACTER(slen), DIMENSION(:), ALLOCATABLE :: filelist
  INTEGER, DIMENSION(:), ALLOCATABLE :: is_array, ie_array, js_array, je_array
  INTEGER :: fid = 21
  INTEGER :: maxtiles = 99999
  LOGICAL :: ex
  LOGICAL :: dodebug = .true.
  INTEGER :: reclen_mult = 4
  REAL    :: undef = -999.0000  ! perhaps 0.0 is better for analysis increments
  INTEGER :: kmax = 32

  NAMELIST /grid_def/ glon,glat,glev,ntiles

  !-----------------------------------------------------------------------------
  ! Read namelist to override defaults
  !-----------------------------------------------------------------------------
  INQUIRE(FILE="input_mt.nml", EXIST=ex)
  if (ex) then
    OPEN(fid,file="input_mt.nml", status='OLD')
    READ(fid,nml=grid_def)
    CLOSE(fid)
  else
    WRITE(6,*) "Please supply namelist to override defaults."
  endif
  WRITE(6,grid_def)

  !-----------------------------------------------------------------------------
  ! Read command line options to override defaults
  !-----------------------------------------------------------------------------
  CALL process_command_line

  !-----------------------------------------------------------------------------
  ! Read filelist
  !-----------------------------------------------------------------------------
  allocate(filelist(ntiles))
  open(FILE=filelist_file,UNIT=fid,FORM='formatted')
  do i=1,ntiles
    read(fid,'(A)') filelist(i)
  enddo
  close(fid)

  !-----------------------------------------------------------------------------
  ! Read grid definitions
  !-----------------------------------------------------------------------------
  allocate(is_array(ntiles),ie_array(ntiles),js_array(ntiles),je_array(ntiles))
  open(FILE=tiledef_file,UNIT=fid,FORM='formatted')
  do i=1,ntiles
    read(fid,*) is_array(i), ie_array(i), js_array(i), je_array(i)
  enddo
  close(fid)

  !-----------------------------------------------------------------------------
  ! Read global land/sea mask
  !-----------------------------------------------------------------------------
  allocate(bufint_global(glon,glat,glev))
  CALL read_global_int(trim(lsmask_file),bufint_global,1)

  !-----------------------------------------------------------------------------
  ! Read global background field as the baseline
  !-----------------------------------------------------------------------------
  allocate(buf4_global(glon,glat,glev))
  CALL read_global(trim(infile),buf4_global,kmax)

  !-----------------------------------------------------------------------------
  ! Read tiled data and put into global array (add to -buf4_global to get an
  ! increment instead)
  !-----------------------------------------------------------------------------
  do i=1,ntiles
    is = is_array(i)
    ie = ie_array(i)
    js = js_array(i)
    je = je_array(i)
    nlon = ie - is + 1
    nlat = je - js + 1
    nlev = glev
    allocate(buf4(nlon,nlat,nlev))

    ! Read subgrid
    CALL read_tile(trim(filelist(i)),buf4,kmax)

    ! Assign subgrid to section of global output array
    ! If specified, compute increment instead of analysis field
    if (DO_INC) then
      where (bufint_global(is:ie,js:je,1:nlev)>0) 
        buf4_global(is:ie,js:je,1:nlev) = buf4(1:nlon,1:nlat,1:nlev) - buf4_global(is:ie,js:je,1:nlev)
      end where
    else
      where (bufint_global(is:ie,js:je,1:nlev)>0) 
        buf4_global(is:ie,js:je,1:nlev) = buf4(1:nlon,1:nlat,1:nlev)
      end where
    endif

    deallocate(buf4) 
  enddo

  !-----------------------------------------------------------------------------
  ! Double check land/sea mask to filter land points
  !-----------------------------------------------------------------------------
  if (DO_INC) undef = 0.0
  lsmask : where (bufint_global==0)
    buf4_global = undef
  end where lsmask

  !-----------------------------------------------------------------------------
  ! Write global data
  !-----------------------------------------------------------------------------
  CALL write_global(trim(outfile),buf4_global,kmax)

  print *, "Deallocating bufint_global..."
  if (allocated(bufint_global)) then
    deallocate(bufint_global)
  endif
  print *, "Deallocating buf4_global..."
  if (allocated(buf4_global)) then
    deallocate(buf4_global)
  endif

CONTAINS


!===============================================================================
! Read a global field of integers (e.g. land/sea mask)
!===============================================================================
SUBROUTINE read_global_int(filename,bufint,kmax)
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: bufint
  INTEGER, INTENT(IN) :: kmax
  INTEGER, PARAMETER :: fid=20
  LOGICAL   :: ex
  INTEGER*8 :: mx, ny, lz, reclen
  LOGICAL :: dodebug = .false.

  mx = SIZE(bufint,1)
  ny = SIZE(bufint,2)
  lz = MIN(SIZE(bufint,3),kmax)
  reclen = mx * ny * lz * reclen_mult

  inquire (file=trim(filename), exist=ex)
  if (ex) then
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    read (fid,rec=1) bufint(:,:,1:kmax)
    close (fid)
    write (6, '( ''read_ncoda_intfile:: Finished reading integer file: '',a)' ) trim(filename)
  else
    write (6, '( ''read_ncoda_intfile:: Integer file missing: '',  a)' ) trim(filename)
    STOP("read_ncoda_intfile:: EXITING...")
  endif

  if (dodebug) then
    print *, "Max(bufint) = ", MAXVAL(bufint)
    print *, "Min(bufint) = ", MINVAL(bufint)
    print *, "Avg(bufint) = ", SUM(bufint)/(mx*ny*lz)
    print *, "Maxloc(bufint) = ", MAXLOC(bufint)
    print *, "Minloc(bufint) = ", MINLOC(bufint)
  endif

END SUBROUTINE read_global_int


!===============================================================================
! Read global data
!===============================================================================
SUBROUTINE read_global(filename,buf4,kmax)
  CHARACTER(*), INTENT(IN) :: filename
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: buf4
  INTEGER, INTENT(IN) :: kmax
  INTEGER*8 :: mx,ny,lz
  INTEGER*8 :: reclen
  LOGICAL :: ex
  INTEGER :: fid = 23
  LOGICAL :: dodebug = .false.

  mx = SIZE(buf4,1)
  ny = SIZE(buf4,2)
  lz = MIN(SIZE(buf4,3),kmax)
  reclen = mx * ny * lz * reclen_mult

  if (dodebug) then
    write(6,*) "mx     = ", mx
    write(6,*) "ny     = ", ny
    write(6,*) "lz     = ", lz
    write(6,*) "_mult  = ", reclen_mult
    write(6,*) "mx * ny * lz * reclen_mult = ", mx * ny * lz * reclen_mult
    write(6,*) "reclen = ", reclen
  endif

  inquire (file=trim(filename), exist=ex)
  if (ex) then
    if (dodebug) write(6,*) "read_global: Opening file: ", trim(filename)
    if (dodebug) write(6,*) "with reclen: ", reclen
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    ! This may take a long time, so output a note of what is happening...
    write (6, *) 'read_global:: Reading background file: ', trim(filename)
    read (fid,rec=1) buf4(:,:,1:kmax)
    close (fid)
    write (6, *) 'read_global:: Finished reading background file: ', trim(filename)
  else
    write (6, *) 'read_global:: global background file missing: ', trim(filename)
    STOP("merge_tiles.f90:: EXITING...")
  endif

  if (dodebug) then
    print *, "Max(buf4) = ", MAXVAL(buf4)
    print *, "Min(buf4) = ", MINVAL(buf4)
    print *, "Avg(buf4) = ", SUM(buf4)/(mx*ny*lz)
    print *, "Maxloc(buf4) = ", MAXLOC(buf4)
    print *, "Minloc(buf4) = ", MINLOC(buf4)
  endif

END SUBROUTINE read_global


!===============================================================================
! Read tiled data
!===============================================================================
SUBROUTINE read_tile(filename,buf4,kmax)
  CHARACTER(*), INTENT(IN) :: filename
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: buf4
  INTEGER, INTENT(IN) :: kmax
  INTEGER :: reclen
  INTEGER :: mx,ny,lz
  LOGICAL :: ex
  INTEGER :: fid = 22
  LOGICAL :: dodebug = .false.

  mx = SIZE(buf4,1)
  ny = SIZE(buf4,2)
  lz = MIN(SIZE(buf4,3),kmax)

  reclen = mx * ny * lz * reclen_mult

  if (dodebug) then
    write(6,*) "mx     = ", mx
    write(6,*) "ny     = ", ny
    write(6,*) "lz     = ", lz
    write(6,*) "_mult  = ", reclen_mult
    write(6,*) "mx * ny * lz * reclen_mult = ", mx * ny * lz * reclen_mult
    write(6,*) "reclen = ", reclen
  endif

  inquire (file=trim(filename), exist=ex)
  if (ex) then
    if (dodebug) write(6,*) "read_tile:: Opening file: ", trim(filename)
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    read (fid,rec=1) buf4(:,:,1:kmax)
    close (fid)
    if (dodebug) write (6, *) 'read_tile:: Finished reading analysis tile file: ', trim(filename)
  else
    write (6, *) 'read_tile:: analysis tile file missing: ', trim(filename)
    STOP("merge_tiles.f90:: EXITING...")
  endif

  if (dodebug) then
    print *, "Max(buf4) = ", MAXVAL(buf4)
    print *, "Min(buf4) = ", MINVAL(buf4)
    print *, "Avg(buf4) = ", SUM(buf4)/(mx*ny*lz)
    print *, "Maxloc(buf4) = ", MAXLOC(buf4)
    print *, "Minloc(buf4) = ", MINLOC(buf4)
  endif

END SUBROUTINE read_tile


!===============================================================================
! Write global data
!===============================================================================
SUBROUTINE write_global(filename,buf4,kmax)
  CHARACTER(*), INTENT(IN) :: filename
  REAL, DIMENSION(:,:,:), INTENT(IN) :: buf4
  INTEGER, INTENT(IN) :: kmax
  INTEGER*8 :: mx,ny,lz
  INTEGER*8 :: reclen
  LOGICAL :: ex
  INTEGER :: fid = 23
  LOGICAL :: dodebug = .true.

  mx = SIZE(buf4,1)
  ny = SIZE(buf4,2)
  lz = MIN(SIZE(buf4,3),kmax)

  reclen = mx * ny * lz * reclen_mult

  if (dodebug) then
    write(6,*) "mx     = ", mx
    write(6,*) "ny     = ", ny
    write(6,*) "lz     = ", lz
    write(6,*) "_mult  = ", reclen_mult
    write(6,*) "mx * ny * lz * reclen_mult = ", mx * ny * lz * reclen_mult
    write(6,*) "reclen = ", reclen
  endif

  inquire (file=trim(filename), exist=ex)
  if (ex) then
    if (dodebug) write(6,*) "write_global: Opening file: ", trim(filename)
    if (dodebug) write(6,*) "with reclen: ", reclen
    open (unit=fid, file=trim(filename), status='old', access='direct', form='unformatted', recl=reclen)
    ! This may take a long time, so output a note of what is happening...
    write (6, *) 'write_global:: Writing analysis file: ', trim(filename)
    write (fid,rec=1) buf4(:,:,1:kmax)
    close (fid)
    write (6, *) 'write_global:: Finished writing analysis file: ', trim(filename)
  else
    write (6, *) 'write_global:: global analysis file missing: ', trim(filename)
    STOP("merge_tiles.f90:: EXITING...")
  endif

  if (dodebug) then
    print *, "Max(buf4) = ", MAXVAL(buf4)
    print *, "Min(buf4) = ", MINVAL(buf4)
    print *, "Avg(buf4) = ", SUM(buf4)/(mx*ny*lz)
!   print *, "Maxloc(buf4) = ", MAXLOC(buf4)
!   print *, "Minloc(buf4) = ", MINLOC(buf4)
  endif

END SUBROUTINE write_global

!===============================================================================
! Process command line arguments
!===============================================================================
SUBROUTINE process_command_line
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In merge_tiles.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)
  select case (arg1)
    case('-f')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      infile = arg2
    case('-flist')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      filelist_file = arg2
    case('-o')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      outfile = arg2
    case('-lsmask')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      lsmask_file = arg2
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
    case('-inc')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_INC
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug
    case('-kmax')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) kmax
    case('-undef')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) undef
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line


END PROGRAM merge_tiles
