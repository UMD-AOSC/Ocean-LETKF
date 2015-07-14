MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!  04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!
!=======================================================================
  USE common
  USE common_mtx
  USE params_letkf, ONLY: nbv

  IMPLICIT NONE

  PUBLIC

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!   OUTPUT
!     trans(nbv,nbv) : transformation matrix
!=======================================================================
SUBROUTINE letkf_core(nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  REAL(r_size),INTENT(IN) :: hdxb(1:nobs,1:nbv)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN) :: rloc(1:nobs)
  REAL(r_size),INTENT(IN) :: dep(1:nobs)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size),INTENT(OUT) :: trans(nbv,nbv)
  REAL(r_size) :: hdxb_rinv(nobsl,nbv)
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv,nobsl)
  REAL(r_size) :: work3(nbv)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.004d0 !0.04d0 !error stdev of parm_infl
  INTEGER :: i,j,k
  REAL(r_size) :: parm_in !STEVE (OCEAN)
  LOGICAL :: debug_hdxb_0 = .false.

  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,nbv
      trans(i,i) = SQRT(parm_infl)
    END DO
    RETURN
  ELSE

  !STEVE: store input inflation value  
  parm_in = parm_infl

!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  !STEVE: debug
  if ( MINVAL(rdiag(1:nobsl)) .le. 0.0 ) then
    WRITE(6,*) "common_letkf.f90:: ERROR: rdiag ≤ 0 (i.e. there is an obserr ≤ 0)"
    WRITE(6,*) "nbv = ", nbv
    WRITE(6,*) "rdiag = ",rdiag
    stop 1
  endif
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  CALL dgemm('t','n',nbv,nbv,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
    & nobsl,0.0d0,work1,nbv)
!DGEMM - Performs one of the matrix-matrix operations
!     C := alpha*op( A )*op( B ) + beta*C
!     where  op( X ) is one of
!        op( X ) = X   or   op( X ) = X',
!     alpha and beta are scalars, and A, B and C are matrices,
!     with op( A ) an m by k matrix,  op( B )  a  k by n matrix
!     and  C an m by n matrix.

!  DO j=1,nbv
!    DO i=1,nbv
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / parm_infl
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
    !STEVE: error check
!   if ( isnan(work1(i,i)) ) then
!     print *, "work1(i,i) = ", work1(i,i)
!     print *, "i = ", i
!     print *, "rho = ", rho
!     stop 2
!   endif
    !STEVE: end
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
  !STEVE: debug
  if ( MINVAL(eival) .le. 0.0 ) then
    WRITE(6,*) "common_letkf.f90:: ERROR: matrix eigenvalue ≤ 0"
    WRITE(6,*) "eival = ", eival
    WRITE(6,*) "eivec = ", eivec
    WRITE(6,*) "i = ", i
    WRITE(6,*) "nbv = ", nbv
    WRITE(6,*) "work1 = ", work1
    WRITE(6,*) "hdxb(1:nobsl,:) = ", hdxb(1:nobsl,:)
    WRITE(6,*) "hdxb_rinv(1:nobsl,:) = ", hdxb_rinv(1:nobsl,:)
    WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl) 
    WRITE(6,*) "rloc(1:nobsl) = ", rloc(1:nobsl)
    !STEVE: (1) This should not happen. (2) If it does, consider replacing non positive evals with mean of all positive evals...
    stop 1
  endif
  !STEVE: end
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,pa,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  CALL dgemm('n','t',nbv,nobsl,nbv,1.0d0,pa,nbv,hdxb_rinv,&
    & nobsl,0.0d0,work2,nbv)
!  DO j=1,nobsl
!    DO i=1,nbv
!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
!      DO k=2,nbv
!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,trans,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      trans(i,j) = trans(i,j) + work3(i)
    END DO
  END DO

  RETURN
  END IF
END SUBROUTINE letkf_core

END MODULE common_letkf
