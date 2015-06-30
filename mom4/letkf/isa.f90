MODULE isa

IMPLICIT NONE
!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
PRIVATE r_size, r_dble, r_sngl
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_dble=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

CONTAINS

LOGICAL FUNCTION isnan(a) 
real(r_size) ::  a 
if (a.ne.a) then 
  isnan = .true. 
else 
  isnan = .false. 
end if 
return 
END FUNCTION isnan 

LOGICAL FUNCTION isinf(a) 
real(r_size) ::  a 
if ((a*0).ne.0) then 
  isinf = .true. 
else 
  isinf = .false. 
end if 
return 
END FUNCTION isinf 

LOGICAL FUNCTION isnan4(a) 
real(r_sngl) ::  a 
if (a.ne.a) then 
  isnan4 = .true. 
else 
  isnan4 = .false. 
end if 
return 
END FUNCTION isnan4 

LOGICAL FUNCTION isinf4(a) 
real(r_sngl) ::  a 
if ((a*0).ne.0) then 
  isinf4 = .true. 
else 
  isinf4 = .false. 
end if 
return 
END FUNCTION isinf4 

END MODULE isa
