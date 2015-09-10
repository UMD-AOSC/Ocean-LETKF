PROGRAM read_drifters
IMPLICIT NONE

INTERFACE
RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices)
    REAL,    DIMENSION(:), INTENT(INOUT) :: Item      ! array of values
    INTEGER,               INTENT(IN)    :: First,Last
    INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
END SUBROUTINE Quicksort
END INTERFACE

INTEGER :: np, i,j,k,ii,jj
INTEGER, PARAMETER :: slen = 256
CHARACTER(slen) :: infile, outfile
INTEGER :: it_id, nf, nd
REAL(kind=8), DIMENSION(:), ALLOCATABLE :: intime
REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: positions
REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: fields
INTEGER, DIMENSION(:), ALLOCATABLE :: index_time, ids
INTEGER :: use_flag
INTEGER, PARAMETER :: outdim = 10000
REAL(kind=8), DIMENSION(outdim) :: lon, lat, depth, temp, salt, time, work
INTEGER, DIMENSION(outdim) :: itime, dids
INTEGER, DIMENSION(outdim) :: didx
INTEGER :: num_drifters, num_times
INTEGER :: outidx = 0
INTEGER, PARAMETER :: fid = 33
INTEGER :: write_option = 0 ! 0 = write all, 1 = write last 1 observation per id

! Input: number of processors
! -np <number of processors> -o <outfile>
print *, "Input -np <number of processors> -o <outfile>"
CALL process_command_line

! Cycle over all processors and read in only files with it_id > 0
do ii=1,np
  WRITE(infile,'(A,I6.6)') 'drifters_out.nc.',ii
  print *, "Processing file: ", infile
   CALL ncread(trim(infile))

  ! If there's data, compile it into one array to be sorted next
  if (use_flag > 0) then
    lon(outidx+1:outidx+it_id) = positions(1,:) 
    lat(outidx+1:outidx+it_id) = positions(2,:) 
    depth(outidx+1:outidx+it_id) = positions(3,:) 
    temp(outidx+1:outidx+it_id) = fields(4,:)
    salt(outidx+1:outidx+it_id) = fields(5,:)
    time(outidx+1:outidx+it_id) = intime
    itime(outidx+1:outidx+it_id) = index_time
    dids(outidx+1:outidx+it_id) = ids
    outidx=outidx+it_id
  endif
enddo
print *, "outidx = ", outidx

CALL sort_drifters(REAL(dids(1:outidx)),1,outidx,outidx)

! Count unique ids
j=1
didx(j) = 1
do i=2,outidx
  if (dids(i) .ne. dids(i-1)) then
    j = j+1
    didx(j) = i
  endif
enddo
didx(j+1) = outidx+1 ! Use this to make the printing loop below simpler.
num_drifters=j
print *, "num_drifters = ", num_drifters

! Sort on the subset of unique ids, by index time
do k=1,j
  CALL sort_drifters(REAL(itime(1:outidx)),didx(k),didx(k+1)-1,outidx)
enddo
num_times=didx(2)-didx(1)
print *, "num_times = ", num_times

! Open outfile
open(fid,FILE=outfile)
! Write heading labels
write(fid,'(A16,I8,A16,I8)')  "num_drifters = ", num_drifters, "num_times = ", num_times
write(fid,'(8A12)') "ids","lon","lat","depth (m)","temp (ºC)","salt (psu)","time (s)","idx time"

! Write data
if (write_option .eq. 0) then
  ! Write all positions of each drifter:
  do k=1,j
    j = didx(k)
    jj = didx(k+1)-1
    do i=j,jj
      write(fid,'(I12,6F12.4,I12)') dids(i), lon(i), lat(i), depth(i), temp(i), salt(i), time(i), itime(i)
    enddo
  enddo
elseif (write_option .eq. 1) then
  ! Write the last position of each drifter:
  do k=1,j
    i = didx(k+1)-1
    write(fid,'(I12,6F12.4,I12)') dids(i), lon(i), lat(i), depth(i), temp(i), salt(i), time(i), itime(i)
  enddo
endif

close(fid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ncread(infile)
USE netcdf
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(*) :: infile
CHARACTER(32) :: dimname
INTEGER :: ncid,istat,varid,dimid,dimlen
INTEGER :: ndims_in, nvars_in, ngatts_in, unlimdimid_in
INTEGER :: i
it_id=0

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
call check( nf90_open(infile, NF90_NOWRITE, ncid) )

! NQ tells how many netCDF variables, dimensions, and global attributes are in the
! file; also the dimension id of the unlimited dimension, if there is one.
call check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )
print *, "ndims_in, nvars_in, ngatts_in, unlimdimid_in = ", ndims_in, nvars_in, ngatts_in, unlimdimid_in

! Inquire about dimension information
call check( nf90_inq_dimid(ncid, "it_id", dimid) )
call check( nf90_Inquire_Dimension(ncid, dimid, dimname, it_id) )
print *, "dimname = ", dimname
print *, "dimlen  = ", it_id

if (it_id < 1) then
  print *, "No data in ", infile
  use_flag = 0
  call check( nf90_close(ncid) )
  return
else
  use_flag = 1
endif

! Inquire about dimension information
call check( nf90_inq_dimid(ncid, "nf", dimid) )
call check( nf90_Inquire_Dimension(ncid, dimid, dimname, nf) )
print *, "dimname = ", dimname
print *, "dimlen  = ", nf

! Inquire about dimension information
call check( nf90_inq_dimid(ncid, "nd", dimid) )
call check( nf90_Inquire_Dimension(ncid, dimid, dimname, nd) )
print *, "dimname = ", dimname
print *, "dimlen  = ", nd

! Get the varid of the data variable, based on its name.
call check( nf90_inq_varid(ncid, "index_time", varid) )
! Read the data.
if (ALLOCATED(index_time)) then
  DEALLOCATE(index_time)
endif
ALLOCATE(index_time(it_id))
call check( nf90_get_var(ncid, varid, index_time) )
print *, "index_time = ", index_time

! Get the varid of the data variable, based on its name.
call check( nf90_inq_varid(ncid, "time", varid) )
! Read the data.
if (ALLOCATED(intime)) then
  DEALLOCATE(intime)
endif
ALLOCATE(intime(it_id))
call check( nf90_get_var(ncid, varid, intime) )
print *, "time = ", intime

! Get the varid of the data variable, based on its name.
call check( nf90_inq_varid(ncid, "ids", varid) )
! Read the data.
if (ALLOCATED(ids)) then
  DEALLOCATE(ids)
endif
ALLOCATE(ids(it_id))
call check( nf90_get_var(ncid, varid, ids) )
print *, "ids = ", ids

! Get the varid of the data variable, based on its name.
call check( nf90_inq_varid(ncid, "positions", varid) )
! Read the data.
if (ALLOCATED(positions)) then
  DEALLOCATE(positions)
endif
ALLOCATE (positions(nd,it_id))
call check( nf90_get_var(ncid, varid, positions) )
print *, "positions = "
print *, positions

! Get the varid of the data variable, based on its name.
call check( nf90_inq_varid(ncid, "fields", varid) )

! Read the data.
if (ALLOCATED(fields)) then
  DEALLOCATE(fields)
endif
ALLOCATE (fields(nf,it_id))
call check( nf90_get_var(ncid, varid, fields) )
print *, "fields(4:5,:) = "
do i=1,it_id
  print *, fields(4:5,i)
enddo

! Close the file, freeing all resources.
call check( nf90_close(ncid) )
print *,"*** SUCCESS reading example file ", infile, "! "

END SUBROUTINE ncread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check(status)
  USE netcdf
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, intent (in) :: status
    
  if(status /= NF90_NOERR) then 
    print *, trim(NF90_STRERROR(status))
    stop "Stopped"
  end if
end subroutine check 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE process_command_line

IMPLICIT NONE

CHARACTER(slen) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
DO i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In grd2cor.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case ('-np')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) np
    case('-f')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      infile = arg2
    case('-o')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      outfile = arg2
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
ENDDO

END SUBROUTINE process_command_line

SUBROUTINE sort_drifters(sort_array,start,finish,outidx)
IMPLICIT NONE
REAL, DIMENSION(*), INTENT(IN) :: sort_array
INTEGER, INTENT(IN) :: start,finish,outidx
REAL, DIMENSION(:), ALLOCATABLE :: sortarr
INTEGER, DIMENSION(:), ALLOCATABLE :: indarr ! for sorting
INTEGER :: i

! Sort by ids (dids) and by time
if (ALLOCATED(sortarr)) then
  DEALLOCATE(sortarr)
endif
ALLOCATE(sortarr(outidx))
!sortarr = dids(1:outidx)
sortarr = sort_array(1:outidx)

if (ALLOCATED(indarr)) then
  DEALLOCATE(indarr)
endif
ALLOCATE(indarr(outidx))
indarr = (/ (i,i=1,outidx) /)           ! indexical array

call Quicksort(sortarr,start,finish,indarr)

print *, "before: "
write(*,'(8A12)') "ids","lon","lat","depth (m)","temp (ºC)","salt (psu)","time (s)","idx time"
do i=1,outidx
  write(*,'(I12,6F12.4,I12)') dids(i), lon(i), lat(i), depth(i), temp(i), salt(i), time(i), itime(i)
enddo

! Sort ids
work = dids
do i=1,outidx
  dids(i) = work(indarr(i))
enddo

! Sort index_time
work = itime
do i=1,outidx
  itime(i) = work(indarr(i))
enddo

! Sort lon
work = lon
do i=1,outidx
  lon(i) = work(indarr(i))
enddo

! Sort lat
work = lat
do i=1,outidx
  lat(i) = work(indarr(i))
enddo

! Sort depth
work = depth
do i=1,outidx
  depth(i) = work(indarr(i))
enddo

! Sort temp
work = temp
do i=1,outidx
  temp(i) = work(indarr(i))
enddo

! Sort salt
work = salt
do i=1,outidx
  salt(i) = work(indarr(i))
enddo

! Sort time
work = time
do i=1,outidx
  time(i) = work(indarr(i))
enddo

print *, "after: "
write(*,'(8A12)') "ids","lon","lat","depth (m)","temp (ºC)","salt (psu)","time (s)","idx time"
do i=1,outidx
  write(*,'(I12,6F12.4,I12)') dids(i), lon(i), lat(i), depth(i), temp(i), salt(i), time(i), itime(i)
enddo

END SUBROUTINE sort_drifters


END PROGRAM read_drifters

!----------------------------------------------------------------------------
   !
   ! This file is based on the the routine in "Fortran 90 for Engineers
   ! & 
   ! Scientists" by Nyhoff and Leestma
   !
   ! Note: In the following subroutines, Item is an assumed-shape array
   !       so a program unit that calls these subroutines must:
   !       1. contain this subroutine as an internal subprogram,
   !       2. import this subroutine from a module, or
   !       3. contain an interface block for this subroutine.
   !
   !----------------------------------------------------------------------------

   !-Quicksort------------------------------------------------------------------
   !
   ! Subroutine to sort a list using the quicksort method. Call it with 
   ! First = the lower bound on the subscripts of the array and 
   ! Last  = the upper bound. 
   !
   ! Accepts : Array "Item", array "Indices"
   ! Returns : Array "Item"    (modified) with elements in ascending
   ! order
   !           array "Indices" (modified) with elements 
   !----------------------------------------------------------------------------
   RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices)
   !----------------------------------------------------------------------------
   ! This routine is based on a similar routine in "Fortran 90 for
   ! Engineers & 
   ! Scientists" by Nyhoff and Leestma.  I modified it to return an
   ! integer 
   ! array sorted based on the relationship of the real data in "Item".
   !
   ! Example:
   ! real,    dimension(100) :: randvals,randcopy
   ! integer, dimension(100) :: indarr = (/ (i, i=1,100) /)
   !  ...
   ! call random_number(randvals)               ! F90 intrinsic
   ! subroutine
   ! randcopy = randvals                        ! save for comparison
   ! call Quicksort(randvals,1,size(randvals),indarr)
   ! print *,'sorted - indexed original is ',SUM(randvals -
   ! randcopy(indarr))
   !
   ! TJH 21 Oct 1998
   !----------------------------------------------------------------------------

   REAL,    DIMENSION(:), INTENT(INOUT) :: Item      ! array of values
   INTEGER,               INTENT(IN)    :: First,Last
   INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices

   !--------------------------------------------------------------------
   ! Interface block(s) & Local Variables
   !--------------------------------------------------------------------

   INTERFACE
   SUBROUTINE Split(Item, Low, High, Mid, Indices)
       REAL,    DIMENSION(:), INTENT(INOUT) :: Item
       INTEGER,               INTENT(IN)    :: Low, High
       INTEGER,               INTENT(OUT)   :: Mid
       INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
   END SUBROUTINE Split
   END INTERFACE

   INTEGER   :: Mid

   !--------------------------------------------------------------------

   IF (First < Last) THEN                            ! IF list size >= 2
       CALL Split(Item, First, Last, Mid, Indices)    ! Split it
       CALL Quicksort(Item, First, Mid-1, Indices)    ! Sort left half
       CALL Quicksort(Item, Mid+1, Last,  Indices)    ! Sort right half
   END IF

   END SUBROUTINE Quicksort

   !-Split----------------------------------------------------------------------
   !
   ! Subroutine to split a list into two sublists, using the first
   ! element as a pivot, and return the position of the element about which the 
   ! list was divided. Local variables used are:
   ! Left       : position of the first element
   ! Right      : position of the last element
   ! Pivot      : pivot element
   ! Swap       : used to swap elements
   !
   ! Accepts:   Array Item and positions Low and High of the first and 
   !            last elements
   ! Returns:   Array Item (modified) with elements in ascending order
   !
   ! Note:      Item is an assumed-shape array so a program unit that
   ! calls
   !            this subroutine must:
   !            1. contain this subroutine as an internal subprogram,
   !            2. import this subroutine from a module
   !            3. contain an interface block for this subroutine.
   !----------------------------------------------------------------------------

   SUBROUTINE Split(Item, Low, High, Mid, Indices)

      REAL,    DIMENSION(:), INTENT(INOUT) :: Item
      INTEGER,               INTENT(IN)    :: Low, High
      INTEGER,               INTENT(OUT)   :: Mid
      INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices


      INTEGER ::   Left, Right
      REAL    ::  Pivot,  Swap
      INTEGER :: iPivot, iSwap

      Left   = Low
      Right  = High
      Pivot  = Item(Low)
      iPivot = Indices(Low)

      ! Repeat the following while Left and Right haven't met

      DO
         IF ( Left >= Right ) Exit

         ! Scan right to left to find element < Pivot

         DO
            IF ( Left >= Right .OR. Item(Right) < Pivot ) EXIT
            Right = Right - 1
         END DO

         ! Scan left to right to find element > Pivot

         DO
            IF (Item(Left) > Pivot) EXIT
            Left = Left + 1
         END DO

         ! If Left and Right haven't met, exchange the items

         IF (Left < Right) THEN
            Swap        = Item(Left)            ! EXCHANGE THE ARRAY ITEMS
            Item(Left)  = Item(Right)
            Item(Right) = Swap

            iSwap          = Indices(Left)      ! EXCHANGE THE INDICES ITEMS
            Indices(Left)  = Indices(Right)
            Indices(Right) = iSwap
         END IF

      END DO

      ! Switch element in split position with pivot

      Item(Low)   = Item(Right)                 ! SWITCH ARRAY ELEMS
      Item(Right) = Pivot
      Mid         = Right

      Indices(Low)   = Indices(Right)           ! SWITCH ARRAY ELEMS
      Indices(Right) = iPivot

   END SUBROUTINE Split
