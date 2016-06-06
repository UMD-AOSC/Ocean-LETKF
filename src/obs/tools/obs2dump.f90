PROGRAM obs2dump
!===============================================================================
! PROGRAM: obs2dump
! 
! USES: none
!
! USAGE:
!   Run to print out data from obs2-formatted observation file
! 
! !REVISION HISTORY:
!   04/03/2014 Steve Penny modified for use with OCEAN at NCEP.
!   04/03/2013 Takemasa Miyoshi created for SPEEDY atmospheric model.
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny $
!===============================================================================

  IMPLICIT NONE
  REAL(4) :: wk(8)
  !REAL(8) :: wk(8)
  INTEGER :: ios,n
  CHARACTER(1) :: S
  INTEGER, PARAMETER :: fid=3
  OPEN(fid,FORM='unformatted')
  n=0
  do
    n=n+1
    READ(3,IOSTAT=ios) wk
    if (ios /= 0) THEN
      print '(A)','END OF FILE'
      EXIT
    endif
!   if (NINT(wk(1)) .ne. 3073 .and. wk(5) > 30) then
!   print '(I6,2F7.2,F10.2,2ES12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
    print '(I6,2F7.2,F10.2,4F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8)
    if (wk(6) <= 0) then
      print *, "STEVE: oerr <= 0, must be > 0 ..."
      print *, "STEVE: oerr(n) = ", wk(6)
      print *, "STEVE: n = ", n
      print '(A)','PRESS "S" TO STOP'
      if(S == 'S')then
        EXIT
      else
        CONTINUE
      endif
    endif
!   print '(A)','PRESS "S" TO STOP'
!   READ(5,'(A1)') S
!   if (S == 'S') EXIT
!   else
!     continue
!   endif
  enddo
  CLOSE(fid)
END PROGRAM obs2dump
