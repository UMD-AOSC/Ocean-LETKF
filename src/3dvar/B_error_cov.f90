MODULE B_error_cov




CONTAINS

SUBROUTINE bkg_err_cov_times_vec(x,Bx)
REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d), INTENT(IN) :: x
REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d), INTENT(OUT) :: Bx

idx=0
do k=1,nlev
  do j=1,nlat
    do i=1,nlon
      ! Compute the estimated B matrix times the input vector
      idx=idx+1
      CALL Berrcov(i,j,k,x,Bx(idx))
    enddo
  enddo
enddo


END SUBROUTINE bkg_err_cov_times_vec

!===============================================================================
! Call special function to get b = B*H'*inv(R)*(yo-H*xb)
!===============================================================================
SUBROUTINE bvec(b)
USE obs_vars, ONLY: obsdep
IMPLICIT NONE
REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d), INTENT(OUT) :: b
REAL(r_size), DIMENSION(:), ALLOCATABLE :: y
REAL(r_size), DIMENSION(nlon,nlat,nlev,nv3d) :: x

ALLOCATE(y(nobs))

! from module: array of innovations (yo-Hxb)

! function to compute R*v
y = obsdep(1:nobs) / rdiag(1:nobs)

! function to compute H*v
CALL H_YtoX(y,x)

! function to compute B*v by location:
CALL bkg_err_cov_times_vec(x,b)


END SUBROUTINE bvec

!===============================================================================
! Call special function to get Ax = (I + B*H'*inv(R)*H)*x
!===============================================================================
SUBROUTINE Amat_xvec(x,Ax)
REAL(r_size), DIMENSION(:), INTENT(IN) :: x
REAL(r_size), DIMENSION(:), INTENT(OUT) :: Ax
REAL(r_size), DIMENSION(:), ALLOCATABLE :: y

ALLOCATE(y(nobs))

! Compute Hx
CALL H_XtoY(x,y)

! Compute inv(r)*HX
y = y(1:nobs) / rdiag(1:nobs)

! Compute H'* <above> 
CALL H_YtoX(y,x)

! function to compute B*v by location:
CALL bkg_err_cov_times_vec(x,b)

! Add "1" to all diagonal entries: (I + Bv)

END SUBROUTINE

END MODULE B_error_cov
