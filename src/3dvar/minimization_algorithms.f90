MODULE minimization_algorithms

USE common, ONLY: r_size, r_sngl

PRIVATE

PUBLIC :: pcg



CONTAINS

!===============================================================================
! Preconditioned conjugate gradient (pcg)
!===============================================================================
SUBROUTINE pcg(xb,xa)
USE Bmod, ONLY: B
USE Rmod, ONLY: R ! assumed diagonal
IMPLICIT NONE
REAL(r_size), INTENT(IN)  :: xb
REAL(r_size), INTENT(OUT) :: xa


! A = I + BH'inv(R)H
! b = BH'inv(R)(yo - Hxb)
! A(dx)=b

! Apply each matrix operation as a function call

! Compute cg method

! Call special function to get b = B*H'*inv(R)*(yo-H*xb)
CALL bvec(b)

! Call special function to get Ax0
CALL Amat_xvec(Ax0)

! r0 = b - Ax0
r0 = b - Ax0
! Initial search direction is the residual
p=r0
r2=DOT_PRODUCT(r0,r0)
do i=1,m
  !STEVE: call special function to get A*p
  Ap=A*p
  alpha=r2/DOT_PRODUCT(p,Ap)
  x=x+alpha*p
  r=r-alpha*Ap

  r2new=DOT_PRODUCT(r,r)
  if (r2new<tol) EXIT

  beta=r2new/r2
  p=r+beta*p

  r2=r2new
enddo


END SUBROUTINE pcg


END MODULE minimization_algorithms
