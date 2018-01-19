PROGRAM convert_date

use calendar

implicit none

integer :: day, month, year, hour, minute, second
CHARACTER(4) :: YYYY
CHARACTER(2) :: MM, DD, HH, NN, SS
integer :: values(3)
integer :: julian = 0
integer :: v(9)
integer :: iy, ja, jm, jy, ierr
integer :: iday,imonth,iyear
INTEGER, PARAMETER :: slen = 256
INTEGER :: i,j,k

! Process the command line
CALL process_command_line

values = (/ year, month, day /)

julian = date_to_julian(day, month, year, values, ierr)
! julian = date_to_julian(day, month, year, ierr)

print *, "Date = ", year, month, day
print *, "julian = ", julian

! The Julian date for CE  1985 January  1 00:00:00.0 UT is
! JD 2446066.500000
print *, "Convert julian to <days since 1/1/1985> = ", julian - 2446066.500000 + 1

open(UNIT=20,FILE='aid.dat',STATUS='new')
write(20,*) FLOOR(julian - 2446066.500000 + 1)
write(20,*) year, month, day
close(20)

CONTAINS

SUBROUTINE process_command_line

IMPLICIT NONE

CHARACTER(slen) :: arg1,arg2
INTEGER :: ierr
INTEGER, DIMENSION(3) :: values
! STEVE: add input error handling!^M
! inputs are in the format "-x xxx"^M
DO i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case ('-y')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) year
      WRITE(YYYY,'(i4.4)') year
      print *, "year = ", arg2
      print *, "YYYY = ", YYYY
    case('-m')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) month
      WRITE(MM,'(i2.2)') month
      print *, "month = ", arg2
      print *, "MM = ", MM
    case('-d')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) day
      WRITE(DD,'(i2.2)') day
      print *, "day = ", arg2
      print *, "DD = ", DD
    case('-h')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      READ(arg2,*) hour
      WRITE(HH,'(i2.2)') hour
      print *, "hour = ", arg2
      print *, "HH = ", HH
  end select
END DO

END SUBROUTINE process_command_line




END PROGRAM convert_date

