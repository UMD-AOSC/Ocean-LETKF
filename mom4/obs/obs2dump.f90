PROGRAM obs2dump
  IMPLICIT NONE
  REAL(4) :: wk(8)
  !REAL(8) :: wk(8)
  INTEGER :: ios,n
  CHARACTER(1) :: S
  OPEN(3,FORM='unformatted')
  n=0
  DO
    n=n+1
    READ(3,IOSTAT=ios) wk
    IF(ios /= 0) THEN
      PRINT '(A)','END OF FILE'
      EXIT
    END IF
!   if (NINT(wk(1)) .ne. 3073 .and. wk(5) > 30) then
!   PRINT '(I6,2F7.2,F10.2,2ES12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
    PRINT '(I6,2F7.2,F10.2,4F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8)
    if (wk(6) <= 0) then
      print *, "STEVE: oerr <= 0, must be > 0 ..."
      print *, "STEVE: oerr(n) = ", wk(6)
      print *, "STEVE: n = ", n
      PRINT '(A)','PRESS "S" TO STOP'
      if(S == 'S')then
        EXIT
      else
        CONTINUE
      endif
    endif
!   PRINT '(A)','PRESS "S" TO STOP'
!   READ(5,'(A1)') S
!   IF(S == 'S') EXIT
!   else
!     continue
!   endif
  END DO
  CLOSE(3)
  STOP
END PROGRAM obs2dump
