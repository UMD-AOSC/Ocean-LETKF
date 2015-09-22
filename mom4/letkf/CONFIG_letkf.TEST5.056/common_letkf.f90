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

  PUBLIC :: letkf_core

  PRIVATE

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
  REAL(r_size),INTENT(INOUT) :: parm_infl   !STEVE: INOUT depricated, change to IN
  REAL(r_size),INTENT(OUT) :: trans(nbv,nbv)
  REAL(r_size) :: hdxb_rinv(nobsl,nbv)
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv,nobsl)
  REAL(r_size) :: work3(nbv)
  REAL(r_size) :: rho
  INTEGER :: i,j,k
  REAL(r_size) :: parm_in !STEVE (OCEAN)
  LOGICAL :: debug_hdxb_0 = .false.

  if (nobsl == 0) then

    trans = 0.0d0
    do i=1,nbv
      trans(i,i) = SQRT(parm_infl)
    enddo

    RETURN

  else

    !STEVE: store input inflation value  
    parm_in = parm_infl

    !---------------------------------------------------------------------------
    !  hdxb Rinv
    !---------------------------------------------------------------------------
    !STEVE: debug
    if ( MINVAL(rdiag(1:nobsl)) .le. 0.0 ) then
      WRITE(6,*) "common_letkf.f90:: ERROR: rdiag ≤ 0 (i.e. there is an obserr ≤ 0)"
      WRITE(6,*) "nbv = ", nbv
      WRITE(6,*) "rdiag = ",rdiag
      stop 1
    endif
    do j=1,nbv
      do i=1,nobsl
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
      enddo
    enddo

    !---------------------------------------------------------------------------
    !  hdxb^T Rinv hdxb
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    !  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
    !---------------------------------------------------------------------------
    rho = 1.0d0 / parm_infl
    do i=1,nbv
      work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
    enddo

    !---------------------------------------------------------------------------
    !  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
    !---------------------------------------------------------------------------
    CALL mtx_eigen(1,nbv,work1,eival,eivec,i)

    !!STEVE: LAPACK alternative: https://software.intel.com/en-us/node/469180
    !        "If the eigenvectors are requested, then this routine uses a divide 
    !         and conquer algorithm to compute eigenvalues and eigenvectors."
    !!call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info) !f77
    ! call dsyevd('V','U',nbv, a, lda, w, work, 2*nbv**2+6*n+1, iwork, liwork, info) !f77
    ! OR:
    !!call syevd(a, w [,jobz] [,uplo] [,info]) !f95
    ! call syevd(work1, eival ,'V') !f95
    ! eivec = work1
    !! include files: mkl.fi, lapack.f90 

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

    !---------------------------------------------------------------------------
    !  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
    !---------------------------------------------------------------------------
    do j=1,nbv
      do i=1,nbv
        work1(i,j) = eivec(i,j) / eival(j)
      enddo
    enddo
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

    !---------------------------------------------------------------------------
    !  Pa hdxb_rinv^T
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    !  Pa hdxb_rinv^T dep
    !---------------------------------------------------------------------------
    do i=1,nbv
      work3(i) = work2(i,1) * dep(1)
      do j=2,nobsl
        work3(i) = work3(i) + work2(i,j) * dep(j)
      enddo 
    enddo

    !---------------------------------------------------------------------------
    !  T = sqrt[(m-1)Pa]
    !---------------------------------------------------------------------------
    do j=1,nbv
      rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
      do i=1,nbv
        work1(i,j) = eivec(i,j) * rho
      enddo
    enddo 
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
    do j=1,nbv
      do i=1,nbv
        trans(i,j) = trans(i,j) + work3(i)
      enddo
    enddo

  endif
 
END SUBROUTINE letkf_core

END MODULE common_letkf
