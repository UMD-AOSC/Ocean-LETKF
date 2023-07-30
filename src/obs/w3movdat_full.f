!-----------------------------------------------------------------------
      subroutine w3movdat(rinc,idat,jdat)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3MOVDAT       RETURN A DATE FROM A TIME INTERVAL AND DATE
!   AUTHOR: MARK IREDELL     ORG: WP23       DATE: 98-01-05
!
! ABSTRACT: THIS SUBPROGRAM RETURNS THE DATE AND TIME THAT IS A GIVEN
!   NCEP RELATIVE TIME INTERVAL FROM AN NCEP ABSOLUTE DATE AND TIME.
!   THE OUTPUT IS IN THE NCEP ABSOLUTE DATE AND TIME DATA STRUCTURE.
!
! PROGRAM HISTORY LOG:
!   98-01-05  MARK IREDELL
!
! USAGE:  CALL W3MOVDAT(RINC,IDAT,JDAT)
!
!   INPUT VARIABLES:
!     RINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!     IDAT       INTEGER (8) NCEP ABSOLUTE DATE AND TIME
!                (YEAR, MONTH, DAY, TIME ZONE,
!                 HOUR, MINUTE, SECOND, MILLISECOND)
!
!   OUTPUT VARIABLES:
!     JDAT       INTEGER (8) NCEP ABSOLUTE DATE AND TIME
!                (YEAR, MONTH, DAY, TIME ZONE,
!                 HOUR, MINUTE, SECOND, MILLISECOND)
!                (JDAT IS LATER THAN IDAT IF TIME INTERVAL IS POSITIVE.)
!
! SUBPROGRAMS CALLED:
!     IW3JDN         COMPUTE JULIAN DAY NUMBER     
!     W3FS26         YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
!     W3REDDAT       REDUCE A TIME INTERVAL TO A CANONICAL FORM
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      real rinc(5)
      integer idat(8),jdat(8)
      real rinc1(5),rinc2(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  add the interval to the input time of day and put into reduced form
!  and then compute new date using julian day arithmetic.
      rinc1(1)=rinc(1)
      rinc1(2:5)=rinc(2:5)+idat(5:8)
      call w3reddat(-1,rinc1,rinc2)
      jldayn=iw3jdn(idat(1),idat(2),idat(3))+nint(rinc2(1))
      call w3fs26(jldayn,jdat(1),jdat(2),jdat(3),jdow,jdoy)
      jdat(4)=idat(4)
      jdat(5:8)=nint(rinc2(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endsubroutine


       FUNCTION IW3JDN(IYEAR,MONTH,IDAY)
C$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
C
C SUBPROGRAM: IW3JDN         COMPUTE JULIAN DAY NUMBER
C   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-29
C
C ABSTRACT: COMPUTES JULIAN DAY NUMBER FROM YEAR (4 DIGITS), MONTH,
C   AND DAY. IW3JDN IS VALID FOR YEARS 1583 A.D. TO 3300 A.D.
C   JULIAN DAY NUMBER CAN BE USED TO COMPUTE DAY OF WEEK, DAY OF
C   YEAR, RECORD NUMBERS IN AN ARCHIVE, REPLACE DAY OF CENTURY,
C   FIND THE NUMBER OF DAYS BETWEEN TWO DATES.
C
C PROGRAM HISTORY LOG:
C   87-03-29  R.E.JONES
C   89-10-25  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
C
C USAGE:   II = IW3JDN(IYEAR,MONTH,IDAY)
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IYEAR  ARG LIST  INTEGER   YEAR           ( 4 DIGITS)
C     MONTH  ARG LIST  INTEGER   MONTH OF YEAR   (1 - 12)
C     IDAY   ARG LIST  INTEGER   DAY OF MONTH    (1 - 31)
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IW3JDN FUNTION   INTEGER   JULIAN DAY NUMBER
C                      JAN. 1,1960 IS JULIAN DAY NUMBER 2436935
C                      JAN. 1,1987 IS JULIAN DAY NUMBER 2446797
C
C   REMARKS: JULIAN PERIOD WAS DEVISED BY JOSEPH SCALIGER IN 1582.
C     JULIAN DAY NUMBER #1 STARTED ON JAN. 1,4713 B.C. THREE MAJOR
C     CHRONOLOGICAL CYCLES BEGIN ON THE SAME DAY. A 28-YEAR SOLAR
C     CYCLE, A 19-YEAR LUNER CYCLE, A 15-YEAR INDICTION CYCLE, USED
C     IN ANCIENT ROME TO REGULATE TAXES. IT WILL TAKE 7980 YEARS
C     TO COMPLETE THE PERIOD, THE PRODUCT OF 28, 19, AND 15.
C     SCALIGER NAMED THE PERIOD, DATE, AND NUMBER AFTER HIS FATHER
C     JULIUS (NOT AFTER THE JULIAN CALENDAR). THIS SEEMS TO HAVE
C     CAUSED A LOT OF CONFUSION IN TEXT BOOKS. SCALIGER NAME IS
C     SPELLED THREE DIFFERENT WAYS. JULIAN DATE AND JULIAN DAY
C     NUMBER ARE INTERCHANGED. A JULIAN DATE IS USED BY ASTRONOMERS
C     TO COMPUTE ACCURATE TIME, IT HAS A FRACTION. WHEN TRUNCATED TO
C     AN INTEGER IT IS CALLED AN JULIAN DAY NUMBER. THIS FUNCTION
C     WAS IN A LETTER TO THE EDITOR OF THE COMMUNICATIONS OF THE ACM
C     VOLUME 11 / NUMBER 10 / OCTOBER 1968. THE JULIAN DAY NUMBER
C     CAN BE CONVERTED TO A YEAR, MONTH, DAY, DAY OF WEEK, DAY OF
C     YEAR BY CALLING SUBROUTINE W3FS26.
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY Y-MP8/864, CRAY Y-MP EL2/256
C
C$$$
C
       IW3JDN  =    IDAY - 32075
     &            + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4
     &            + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12
     &            - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4
       RETURN
       ENDFUNCTION


      subroutine w3reddat(it,rinc,dinc)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3REDDAT       REDUCE A TIME INTERVAL TO A CANONICAL FORM
!   AUTHOR: MARK IREDELL     ORG: WP23       DATE: 98-01-05
!
! ABSTRACT: THIS SUBPROGRAM REDUCES AN NCEP RELATIVE TIME INTERVAL
!   INTO ONE OF SEVEN CANONICAL FORMS, DEPENDING ON THE INPUT IT VALUE.
!
!   First reduced format type (IT=-1):
!        RINC(1) is an arbitrary integer.
!        RINC(2) is an integer between 00 and 23, inclusive.
!        RINC(3) is an integer between 00 and 59, inclusive.
!        RINC(4) is an integer between 00 and 59, inclusive.
!        RINC(5) is an integer between 000 and 999, inclusive.
!      If RINC(1) is negative, then the time interval is negative.
!    
!   Second reduced format type (IT=0):
!      If the time interval is not negative, then the format is:
!        RINC(1) is zero or a positive integer. 
!        RINC(2) is an integer between 00 and 23, inclusive.
!        RINC(3) is an integer between 00 and 59, inclusive.
!        RINC(4) is an integer between 00 and 59, inclusive.
!        RINC(5) is an integer between 000 and 999, inclusive.
!      Otherwise if the time interval is negative, then the format is:
!        RINC(1) is zero or a negative integer. 
!        RINC(2) is an integer between 00 and -23, inclusive.
!        RINC(3) is an integer between 00 and -59, inclusive.
!        RINC(4) is an integer between 00 and -59, inclusive.
!        RINC(5) is an integer between 000 and -999, inclusive.
!    
!   Days format type (IT=1):
!        RINC(1) is arbitrary.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is zero.
!    
!   Hours format type (IT=2):
!        RINC(1) is zero.
!        RINC(2) is arbitrary.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than 300 years.)
!    
!   Minutes format type (IT=3):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is arbitrary.
!        RINC(4) is zero.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than five years.)
!    
!   Seconds format type (IT=4):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is arbitrary.
!        RINC(5) is zero.
!      (This format should not express time intervals longer than one month.)
!    
!   Milliseconds format type (IT=5):
!        RINC(1) is zero.
!        RINC(2) is zero.
!        RINC(3) is zero.
!        RINC(4) is zero.
!        RINC(5) is arbitrary.
!     (This format should not express time intervals longer than one hour.)
!
! PROGRAM HISTORY LOG:
!   98-01-05  MARK IREDELL
!
! USAGE:  CALL W3REDDAT(IT,RINC,DINC)
!
!   INPUT VARIABLES:
!     IT         INTEGER RELATIVE TIME INTERVAL FORMAT TYPE
!                (-1 FOR FIRST REDUCED TYPE (HOURS ALWAYS POSITIVE),
!                 0 FOR SECOND REDUCED TYPE (HOURS CAN BE NEGATIVE),
!                 1 FOR DAYS ONLY, 2 FOR HOURS ONLY, 3 FOR MINUTES ONLY,
!                 4 FOR SECONDS ONLY, 5 FOR MILLISECONDS ONLY)
!     RINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!
!   OUTPUT VARIABLES:
!     DINC       REAL (5) NCEP RELATIVE TIME INTERVAL
!                (DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
!
! SUBPROGRAMS CALLED:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      real rinc(5),dinc(5)
!  parameters for number of units in a day
!  and number of milliseconds in a unit
!  and number of next smaller units in a unit, respectively
      integer,dimension(5),parameter:: itd=(/1,24,1440,86400,86400000/),
     &                                 itm=itd(5)/itd
      integer,dimension(4),parameter:: itn=itd(2:5)/itd(1:4)
      integer,parameter:: np=16
      integer iinc(4),jinc(5),kinc(5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  first reduce to the first reduced form
      iinc=floor(rinc(1:4))
!  convert all positive fractional parts to milliseconds
!  and determine canonical milliseconds
      jinc(5)=nint(dot_product(rinc(1:4)-iinc,real(itm(1:4)))+rinc(5))
      kinc(5)=modulo(jinc(5),itn(4))
!  convert remainder to seconds and determine canonical seconds
      jinc(4)=iinc(4)+(jinc(5)-kinc(5))/itn(4)
      kinc(4)=modulo(jinc(4),itn(3))
!  convert remainder to minutes and determine canonical minutes
      jinc(3)=iinc(3)+(jinc(4)-kinc(4))/itn(3)
      kinc(3)=modulo(jinc(3),itn(2))
!  convert remainder to hours and determine canonical hours
      jinc(2)=iinc(2)+(jinc(3)-kinc(3))/itn(2)
      kinc(2)=modulo(jinc(2),itn(1))
!  convert remainder to days and compute milliseconds of the day
      kinc(1)=iinc(1)+(jinc(2)-kinc(2))/itn(1)
      ms=dot_product(kinc(2:5),itm(2:5))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  next reduce to either single value canonical form
!  or to one of the two reduced forms
      if(it.ge.1.and.it.le.5) then
!  ensure that exact multiples of 1./np are expressed exactly
!  (other fractions may have precision errors)
        rp=(np*ms)/itm(it)+mod(np*ms,itm(it))/real(itm(it))
        dinc=0
        dinc(it)=real(kinc(1))*itd(it)+rp/np
      else
!  the reduced form is done except the second reduced form is modified
!  for negative time intervals with fractional days
        dinc=kinc
        if(it.eq.0.and.kinc(1).lt.0.and.ms.gt.0) then
          dinc(1)=dinc(1)+1
          dinc(2:5)=mod(ms-itm(1),itm(1:4))/itm(2:5)
        endif
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endsubroutine
       SUBROUTINE W3FS26(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
C$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
C
C SUBPROGRAM: W3FS26         YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
C   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-29
C
C ABSTRACT: COMPUTES YEAR (4 DIGITS), MONTH, DAY, DAY OF WEEK, DAY
C   OF YEAR FROM JULIAN DAY NUMBER. THIS SUBROUTINE WILL WORK
C   FROM 1583 A.D. TO 3300 A.D.
C
C PROGRAM HISTORY LOG:
C   87-03-29  R.E.JONES
C   89-10-25  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
C
C USAGE:  CALL W3FS26(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     JLDAYN ARG LIST  INTEGER   JULIAN DAY NUMBER
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IYEAR  ARG LIST  INTEGER   YEAR  (4 DIGITS)
C     MONTH  ARG LIST  INTEGER   MONTH
C     IDAY   ARG LIST  INTEGER   DAY
C     IDAYWK ARG LIST  INTEGER   DAY OF WEEK (1 IS SUNDAY, 7 IS SAT)
C     IDAYYR ARG LIST  INTEGER   DAY OF YEAR (1 TO 366)
C
C   REMARKS: A JULIAN DAY NUMBER CAN BE COMPUTED BY USING ONE OF THE
C     FOLLOWING STATEMENT FUNCTIONS. A DAY OF WEEK CAN BE COMPUTED
C     FROM THE JULIAN DAY NUMBER. A DAY OF YEAR CAN BE COMPUTED FROM
C     A JULIAN DAY NUMBER AND YEAR.
C
C      IYEAR (4 DIGITS)
C
C      JDN(IYEAR,MONTH,IDAY) = IDAY - 32075
C    &            + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4
C    &            + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12
C    &            - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4
C
C      IYR (4 DIGITS) , IDYR(1-366) DAY OF YEAR
C
C      JULIAN(IYR,IDYR) = -31739 + 1461 * (IYR + 4799) / 4
C    &                    -3 * ((IYR + 4899) / 100) / 4 + IDYR
C
C      DAY OF WEEK FROM JULIAN DAY NUMBER, 1 IS SUNDAY, 7 IS SATURDAY.
C
C      JDAYWK(JLDAYN) = MOD((JLDAYN + 1),7) + 1
C
C      DAY OF YEAR FROM JULIAN DAY NUMBER AND 4 DIGIT YEAR.
C
C      JDAYYR(JLDAYN,IYEAR) = JLDAYN -
C     &  (-31739+1461*(IYEAR+4799)/4-3*((IYEAR+4899)/100)/4)
C
C      THE FIRST FUNCTION WAS IN A LETTER TO THE EDITOR COMMUNICATIONS
C      OF THE ACM  VOLUME 11 / NUMBER 10 / OCTOBER, 1968. THE 2ND
C      FUNCTION WAS DERIVED FROM THE FIRST. THIS SUBROUTINE WAS ALSO
C      INCLUDED IN THE SAME LETTER. JULIAN DAY NUMBER 1 IS
C      JAN 1,4713 B.C. A JULIAN DAY NUMBER CAN BE USED TO REPLACE A
C      DAY OF CENTURY, THIS WILL TAKE CARE OF THE DATE PROBLEM IN
C      THE YEAR 2000, OR REDUCE PROGRAM CHANGES TO ONE LINE CHANGE
C      OF 1900 TO 2000. JULIAN DAY NUMBERS CAN BE USED FOR FINDING
C      RECORD NUMBERS IN AN ARCHIVE OR DAY OF WEEK, OR DAY OF YEAR.
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY Y-MP8/864
C
C$$$
C
       L      = JLDAYN + 68569
       N      = 4 * L / 146097
       L      = L - (146097 * N + 3) / 4
       I      = 4000 * (L + 1) / 1461001
       L      = L - 1461 * I / 4 + 31
       J      = 80 * L / 2447
       IDAY   = L - 2447 * J / 80
       L      = J / 11
       MONTH  = J + 2 - 12 * L
       IYEAR  = 100 * (N - 49) + I + L
       IDAYWK = MOD((JLDAYN + 1),7) + 1
       IDAYYR = JLDAYN -
     &  (-31739 +1461 * (IYEAR+4799) / 4 - 3 * ((IYEAR+4899)/100)/4)
       RETURN
       END

      SUBROUTINE W3FS21(IDATE, NMIN)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:   W3FS21       NUMBER OF MINUTES SINCE JAN 1, 1978
C   PRGMMR: REJONES          ORG: NMC421     DATE: 89-07-17
C
C ABSTRACT: CALCULATES THE NUMBER OF MINUTES SINCE 0000,
C   1 JANUARY 1978.
C
C PROGRAM HISTORY LOG:
C   84-06-21  A. DESMARAIS
C   89-07-14  R.E.JONES    CONVERT TO CYBER 205 FORTRAN 200,
C                          CHANGE LOGIC SO IT WILL WORK IN
C                          21 CENTURY.
C   89-11-02  R.E.JONES    CONVERT TO CRAY CFT77 FORTRAN
C
C USAGE:    CALL W3FS21 (IDATE, NMIN)
C   INPUT ARGUMENT LIST:
C     IDATE    - INTEGER  SIZE 5 ARRAY CONTAINING YEAR OF CENTURY,
C                MONTH, DAY, HOUR AND MINUTE.  IDATE(1) MAY BE
C                A TWO DIGIT YEAR OR 4. IF 2 DIGITS AND GE THAN 78
C                1900 IS ADDED TO IT. IF LT 78 THEN 2000 IS ADDED
C                TO IT. IF 4 DIGITS THE SUBROUTINE WILL WORK
C                CORRECTLY TO THE YEAR 3300 A.D.
C
C   OUTPUT ARGUMENT LIST:
C     NMIN     - INTEGER NUMBER OF MINUTES SINCE 1 JANUARY 1978
C
C   SUBPROGRAMS CALLED:
C     LIBRARY:
C       W3LIB    - IW3JDN
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY Y-MP8/832
C
C$$$
C
      INTEGER  IDATE(5)
      INTEGER  NMIN
      INTEGER  JDN78
C
      DATA  JDN78 / 2443510 /
C
C***   IDATE(1)       YEAR OF CENTURY
C***   IDATE(2)       MONTH OF YEAR
C***   IDATE(3)       DAY OF MONTH
C***   IDATE(4)       HOUR OF DAY
C***   IDATE(5)       MINUTE OF HOUR
C
      NMIN  = 0
C
      IYEAR = IDATE(1)
C
      IF (IYEAR.LE.99) THEN
        IF (IYEAR.LT.78) THEN
          IYEAR = IYEAR + 2000
        ELSE
          IYEAR = IYEAR + 1900
        ENDIF
      ENDIF
C
C     COMPUTE JULIAN DAY NUMBER FROM YEAR, MONTH, DAY
C
      IJDN  = IW3JDN(IYEAR,IDATE(2),IDATE(3))
C
C     SUBTRACT JULIAN DAY NUMBER OF JAN 1,1978 TO GET THE
C     NUMBER OF DAYS BETWEEN DATES
C
      NDAYS = IJDN - JDN78
C
C***  NUMBER OF MINUTES
C
      NMIN = NDAYS * 1440 + IDATE(4) * 60 + IDATE(5)
C
      RETURN
      END



