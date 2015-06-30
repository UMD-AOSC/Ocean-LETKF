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
!  09/01/2013 Steve Penny converted adaptive inflation to SFC only
!  10/01/2014 Steve Penny Rewrite to optimize for accuracy and speed 
!                         using symmetric BLAS routines
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC
!=======================================================================
!  LEKF Model Independent Parameters
!=======================================================================
  INTEGER,PARAMETER :: nbv=56 !8 !28 !2 !28 !4 !ensemble size

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
! REAL(r_size),PARAMETER :: sigma_b = 0.01d0 !(0.001d0 is the 1st and 2nd runs) !error stdev of parm_infl
  !STEVE: I put this in common so I could output it to keep in the records
  INTEGER :: i,j,k
  !STEVE: turn off inflation, but still compute adaptive inflation
  LOGICAL,PARAMETER :: USE_INFL=.true. !STEVE: make this an input argument, and use only for 2d variables
  LOGICAL :: dodebug=.true.

  if (dodebug) then
    print *, "!-----------------------------------------------------------------------"
    print *, "!  letkf_core"
    print *, "!-----------------------------------------------------------------------"
  endif

  IF(nobsl == 0) THEN
    trans = 0.0d0
    ! Grow the unobserved area to maintain smoothness in analysis field
    if (USE_INFL) then
      DO i=1,nbv
        trans(i,i) = SQRT(parm_infl)
      END DO
    endif
    RETURN
  ELSE

!!-----------------------------------------------------------------------
!!  hdxb Rinv
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!  hdxb^T Rinv hdxb
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!!-----------------------------------------------------------------------
!
    if (dodebug) then
      print *, "!-----------------------------------------------------------------------"
      print *, "!  hdxb^T Rinv hdxb + (m-1) I / rho "
      print *, "!-----------------------------------------------------------------------"
    endif

    !STEVE: replacing previous 3 sections
    ! Split diagonal matrix into two: sqrt(Rinv)*sqrt(Rinv), and then merge with hdxb:
    DO j=1,nbv
      DO i=1,nobsl
        hdxb_rinv(i,j) = hdxb(i,j) * SQRT(rloc(i) / rdiag(i))
        if (dodebug) then
          print *, "hdxb_rinv(i,j) = ",
          print *, "hdxb(i,j) = ", hdxb(i,j) 
          print *, "rloc(i) = ", rloc(i)
          print *, "rdiag(i) = ", rdiag(i)
          print *, "SQRT(rloc(i) / rdiag(i)) = ", SQRT(rloc(i) / rdiag(i))
        endif
      END DO
    END DO

    ! Form identity matrix:
    work1 = 0.0d0
    DO i=1,nbv
      work1(i,i) = 1.0d0 
    END DO

    CALL dsyrk( 'L','T',nbv,nobsl,1.0d0,hdxb_rinv,nbv,1.0d0,work1,REAL(nbv-1,r_size)/parm_infl )

    if (dodebug) then
      print *, "work1 = ", work1
    endif

    ! Then fill in the strictly upper diagonal part with the entries from the lower diagonal:
    do ri=1,nbv
      do ci=ri+1,nbv
        if (dodebug) print *, "before, work1(ri,ci) = ", work1(ri,ci)
        work1(ri,ci) = work1(ci,ri)
        if (dodebug) print *, " after, work1(ri,ci) = ", work1(ri,ci)
      enddo
    enddo

!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
    CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
    if (dodebug) then
      print *, "eival = ", eival
      print *, "eivec = ", eivec
    endif

    if ( MINVAL(eival) .le. 0.0 ) then
      print *, "common_letkf.f90:: ERROR: matrix eigenvalue ≤ 0"
      print *, "eival = ", eival
      print *, "nbv = ", nbv
      print *, "i = ", i
      stop 1
    endif

!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
    if (dodebug) then
      print *, "!-----------------------------------------------------------------------"
      print *, "!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv"
      print *, "!-----------------------------------------------------------------------"
    endif
    eival = SQRT(eival)
    DO j=1,nbv
      DO i=1,nbv
        work1(i,j) = eivec(i,j) / eival(j)
        if (dodebug) then
          print *, "work1(",i,",",j,") = ", work1(i,j)
        endif
      END DO
    END DO
    ! CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,nbv,0.0d0,pa,nbv)
    !STEVE: ensure that Pa is symmetric:
    CALL dsyrk( 'L','T',nbv,nbv,1.0d0,work1,nbv,0.0d0,pa,nbv )
    ! Then fill in the strictly upper diagonal part with the entries from the lower diagonal:
    do ri=1,nbv
      do ci=ri+1,nbv
        pa(ri,ci) = pa(ci,ri)
        if (dodebug) then
          print *, "pa(",ri,",",ci") = " pa(ri,ci)
        endif
      enddo
    enddo

!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
    if (dodebug) then
      print *, "!-----------------------------------------------------------------------"
      print *, "!  T = sqrt[(m-1)Pa]"
      print *, "!-----------------------------------------------------------------------"
    endif
    eival = 1.0d0/SQRT(eival)
    DO j=1,nbv
      DO i=1,nbv
        work1(i,j) = eivec(i,j) * eival(j)
        if (dodebug) then
          print *, "work1(",i,",",j,") = ", work1(i,j)
        endif
      END DO
    END DO
    ! CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,nbv,0.0d0,trans,nbv)
    !STEVE: ensure that sqrt(Pa) is symmetric:
    CALL dsyrk( 'L','T',nbv,nbv,1.0d0,work1,nbv,0.0d0,trans,nbv )
    ! Then fill in the strictly upper diagonal part with the entries from the lower diagonal:
    do ri=1,nbv
      do ci=ri+1,nbv
        trans(ri,ci) = trans(ci,ri)
        if (dodebug) then
          print *, "trans(",ri,",",ci") = " trans(ri,ci)
        endif
      enddo
    enddo

!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
    !STEVE: better to multiply hdxb_rinv^T dep, then Pa (hdxb_rinv^T dep)
    !STEVE: added this in here to make up for earlier change to hdxb_rinv
    if (dodebug) then
      print *, "!-----------------------------------------------------------------------"
      print *, "!  Pa hdxb_rinv^T dep"
      print *, "!-----------------------------------------------------------------------"
    endif
    DO j=1,nbv
      DO i=1,nobsl
        hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
      END DO
    END DO
    CALL dgemv ('t', nobsl, nbv, 1.0d0, hdxb_rinv, LDA, dep, 1, 0.0d0, work3, 1)
    CALL dsymv ('L', nbv, 1.0d0, pa, LDA, nbv, 1, 0.0d0, work3, 1)

!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
    if (dodebug) then
      print *, "!-----------------------------------------------------------------------"
      print *, "!  T + Pa hdxb_rinv^T dep"
      print *, "!-----------------------------------------------------------------------"
    endif
    DO j=1,nbv
      trans(:,j) = trans(:,j) + work3(:)
    END DO

    !STEVE: testing...
    STOP(1984)

    return

  ENDIF
END SUBROUTINE letkf_core

END MODULE common_letkf
