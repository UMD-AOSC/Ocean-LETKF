MODULE pzcon

IMPLICIT NONE

PRIVATE 

PUBLIC :: p_from_z

CONTAINS

PURE SUBROUTINE p_from_z(cnvrsn,grav,dynz,prsdb,depth)
!http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/pzcon.f
!     title:
!     *****
!
!       pzcon  -- convert pressure in decibars to depth in meters
!                 (or visa versa)
!
!     system:
!     ******
!
!       pacodf hydrographic data library
!
!     purpose:
!     *******
!
!       to calculate depth in meters (mantyla,saunders)
!       from pressure in decibars (or pressure from depth).
!
!       ref: journ. physical ocean.,vol 11 no. 4, april, 1981
!            (saunders)
!            private correspondence 1982-1983
!            (mantyla)
!
!     method:
!     ******
!
!       a standard ocean (salinity=35.0,t=0.0) is used plus a dynamic
!       height correction to account for deviations from the standard
!       ocean. pressure to depth conversion is effected as:
!
!       z = p/(b*c) + dynz/b
!
!       where:
!
!          p    = insitu pressure (decibars)
!          b    = insitu gravity as a function of latitude and pressure
!                 (decameters/sec/sec)
!          c    = insitu mean density rho(35.0,0.0,p)
!                 (grams/centimeter**3)
!          dynz = dynamic height in dynamic meters
!          z    = depth in meters
!
!     parameters:
!     **********
!
!       cnvrsn  -> conversion to be performed:
!                  0 = pressure to depth
!                  1 = depth to pressure
!       grav    -> acceleration of gravity meters/sec/sec
!                  at station latitude, pressure =0.0db
!       dynz    -> dynamic height in dynamic meters
!       prsdb   -> pressure in decibars (cnvrsn=0)
!               <- pressure in decibars (cnvrsn=1)
!       depth   <- depth in meters (cnvrsn=0)
!               -> depth in meters (cnvrsn=1)
!
  INTEGER, INTENT(IN) :: cnvrsn
  REAL, INTENT(IN)    :: grav,dynz
  REAL, INTENT(INOUT) :: prsdb,depth
!
!     variables:
!     *********
!
  REAL*4 :: a,b,c
  REAL*8 :: dd,dg,dh,da,db,dc
!
!     external functions:
!     ******** *********
!
  INTEGER :: isnan
!	  /* test for Not-A-Number */

!-------------------------------------------------------------------------------
!     (convert pressure to depth):
!-------------------------------------------------------------------------------
  if(cnvrsn.eq.0) then
    a = 2.2e-6*prsdb
! /* pressure correction to gravity */
    b = 0.1*(grav+0.5*a)
! /* insitu gravity decameters/sec/sec */
    c = 1.0285+a
! /* insitu mean density rho(35.0,0.0,p) */
    depth = prsdb/(b*c)
! /* pressure to depth conversion in standard ocean */
!   (dynamic height):
    if(isnan(dynz).eq.0) then
      depth = depth + dynz/b
    endif

  else
!-------------------------------------------------------------------------------
! (convert depth to pressure):
!-------------------------------------------------------------------------------
    dd = dble(depth)
    dg = dble(grav*0.1)
    dh = dble(dynz)
    da = 2.42d-13*dd
    db = dd*(1.13135d-7+2.2e-6*dg)-1.0
    dc = dg*dd
!   select
!   (dynamic height):
    if(isnan(dynz).eq.0) then
      db = db - 2.2d-6*dh
      dc = dc-dh
    else
      dc = 1.0285d0*dc
      prsdb  = 0.0
      if(da.ne.0.0) then
        prsdb  = sngl((-db-dsqrt(db*db-4.0d0*da*dc))/(da+da))
      endif
    endif

  endif

END SUBROUTINE p_from_z

PURE INTEGER FUNCTION isnan(r)
!http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/isnan.f
!
!	title:
!	*****
!
!	  isnan -- test for missing real*4 value
!
!	purpose:
!	*******
!
!	  to test for a real number flagged as missing
!
!	parameters:
!	**********
!
!	  r     -> real value to test
!	  isnan <- 0 iff number, 1 iff missing
!
  REAL, INTENT(IN) :: r

  if(r.eq.9.9e37) then
    isnan = 1
  else
    isnan = 0
  endif

END FUNCTION isnan

END MODULE pzcon
