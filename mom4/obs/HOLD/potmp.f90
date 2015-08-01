SUBROUTINE potmp(PRESS,TEMP,S,RP,POTEMP)
!SOURCE: http://www.ncl.ucar.edu/Support/talk_archives/2011/att-1609/potmp.f
!
!     TITLE:
!     *****
!
!       POTMP  -- CALCULATE POTENTIAL TEMPERATURE FOR AN ARBITRARY
!                 REFERENCE PRESSURE
!
!     PURPOSE:
!     *******
!
!       TO CALCULATE POTENTIAL TEMPERATURE
!
!       REF: N.P. FOFONOFF
!            DEEP SEA RESEARCH
!            IN PRESS NOV 1976
!
!     PARAMETERS:
!     **********
!
!       PRESS  -> PRESSURE IN DECIBARS
!       TEMP   -> TEMPERATURE IN CELSIUS DEGREES
!       S      -> SALINITY PSS 78
!       RP     -> REFERENCE PRESSURE IN DECIBARS
!                 (0.0 FOR THE QUANTITY THETA)
!       POTEMP <- POTENTIAL TEMPERATURE (DEG C)
!
  REAL PRESS,TEMP,S,RP,POTEMP
!
!     VARIABLES:
!     *********
!
  INTEGER I,J,N
  REAL*4 DP,P,Q,R1,R2,R3,R4,R5,S1,T,X
!
!     CODE:
!     ****
!
  S1 = S-35.0
  P  = PRESS
  T  = TEMP

  DP = RP - P
  N  = IFIX(ABS(DP)/1000.) + 1
  DP = DP/FLOAT(N)

  outer : do 10 I=1,N
    inner : do 20 J=1,4

      R1 = ((-2.1687E-16*T+1.8676E-14)*T-4.6206E-13)*P
      R2 = (2.7759E-12*T-1.1351E-10)*S1
      R3 = ((-5.4481E-14*T+8.733E-12)*T-6.7795E-10)*T
      R4 = (R1+(R2+R3+1.8741E-8))*P+(-4.2393E-8*T+1.8932E-6)*S1
      R5 = R4+((6.6228E-10*T-6.836E-8)*T+8.5258E-6)*T+3.5803E-5

      X  = DP*R5

      GO TO (100,200,300,400),J

  100 CONTINUE
      T = T+.5*X
      Q = X
      P = P + .5*DP
      CYCLE inner

  200 CONTINUE
      T = T + .29298322*(X-Q)
      Q = .58578644*X + .121320344*Q
      CYCLE inner

  300 CONTINUE
      T = T + 1.707106781*(X-Q)
      Q = 3.414213562*X - 4.121320344*Q
      P = P + .5*DP
      CYCLE inner

  400 CONTINUE
      T = T + (X-2.0*Q)/6.0

    enddo inner
  enddo outer

      POTEMP = T
!
!       END POTMP
!

END SUBROUTINE potmp

