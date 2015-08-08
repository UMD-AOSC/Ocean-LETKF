PROGRAM convert_infile

  IMPLICIT NONE

  INTEGER :: fid = 101
  INTEGER :: nn = 1000
  INTEGER :: n, i, ii
  
  INTEGER, ALLOCATABLE :: elem(:)
  REAL, ALLOCATABLE :: rlon(:)
  REAL, ALLOCATABLE :: rlat(:)
  REAL, ALLOCATABLE :: rlev(:)
  REAL, ALLOCATABLE :: oerr(:)
  REAL, ALLOCATABLE :: otime(:)
  REAL, ALLOCATABLE :: oqc(:)
  INTEGER, ALLOCATABLE :: obsid(:)
  CHARACTER(LEN=12) :: head_line(7)
  CHARACTER(LEN=12) :: head_line1(7)
  REAL :: wk(7)
  
  ALLOCATE(  elem(nn)  )
  ALLOCATE(  rlon(nn)  )
  ALLOCATE(  rlat(nn)  )
  ALLOCATE(  rlev(nn)  )
  ALLOCATE(  obsid(nn) )
  ALLOCATE(  oerr(nn)  )
  ALLOCATE(   oqc(nn)  ) 
  ALLOCATE(  otime(nn) )	

  open(fid,FILE='drifters_out.txt')
  read(fid,*) head_line
  read(fid,*) head_line1

  DO n=1,nn
    READ(fid,*) wk
    DO i = 1,3
      ii = 3 * (n-1) + i
      SELECT CASE(i)
        CASE(1)
          elem(ii) = 1111
        CASE(2)
          elem(ii) = 2222
        CASE(3)
          elem(ii) = 3333
      END SELECT
      obsid(ii) = INT(wk(1))
      rlon(ii) = REAL(wk(2))
      rlat(ii) = REAL(wk(3))
      rlev(ii) = REAL(wk(4))
      oerr(ii) = REAL(wk(5))
       oqc(ii) = REAL(wk(6))
     otime(ii) = REAL(wk(7))
     
    PRINT *, elem(ii), rlon(ii), rlat(ii), rlev(ii), obsid(ii), oerr(ii), oqc(ii), otime(ii)
    end DO  
  end DO
  close(fid)

END PROGRAM convert_infile
