PROGRAM obs2_verify
  IMPLICIT NONE
  REAL(4) :: wkb(8)
  REAL(4) :: wka(8)
  !REAL(8) :: wk(8)
  INTEGER :: ios,n
  CHARACTER(1) :: S
  OPEN(3,FORM='unformatted')
  OPEN(4,FORM='unformatted')
  n=0
  DO
    n=n+1
    READ(3,IOSTAT=ios) wkb
    READ(4,IOSTAT=ios) wka
    IF(ios /= 0) THEN
      PRINT '(A)','END OF FILE'
      EXIT
    END IF
!   if (NINT(wk(1)) .ne. 3073 .and. wk(5) > 30) then
!   PRINT '(I6,2F7.2,F10.2,2ES12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
    PRINT '(I6,2F7.2,F10.2,4F12.2)',NINT(wkb(1)),wkb(2),wkb(3),wkb(4),wkb(5),wkb(6),wkb(7),wkb(8)
    PRINT '(I6,2F7.2,F10.2,4F12.2)',NINT(wka(1)),wka(2),wka(3),wka(4),wka(5),wka(6),wka(7),wka(8)
    if (wkb(6) <= 0) then
      print *, "STEVE: oerr <= 0, must be > 0 ..."
      print *, "STEVE: oerr(n) = ", wkb(6)
      print *, "STEVE: n = ", n
    endif
    if (wkb(1) .ne. wka(1) .or. &
        wkb(2) .ne. wka(2) .or. &
        wkb(3) .ne. wka(3) .or. &
        wkb(4) .ne. wka(4) .or. &
        wkb(5) .ne. wka(5) .or. &
        wkb(6) .ne. wka(6) .or. &
!       wkb(7) .ne. wka(7) .or. &
        wkb(8) .ne. wka(8) ) then
      print *, "Records do not correspond :: rec n = ", n

      PRINT '(A)','PRESS "S" TO STOP'
      READ(5,'(A1)') S
      IF(S == 'S') EXIT
    endif
  
  END DO
  CLOSE(3)
  STOP
END PROGRAM obs2_verify
