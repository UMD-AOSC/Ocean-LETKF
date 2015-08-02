MODULE so_vertical_covariance



CONTAINS


!===============================================================================
! <SUBROUTINE NAME="vertical_covariance"
!
! <DESCRIPTION>
! Compute the vertical covariance matrix.
! A smoother is used.  The scale is proportional to layer thickness.
! </DESCRIPTION>
!
!===============================================================================
  SUBROUTINE vertical_covariance

  INTEGER, PARAMETER                :: kex = 10, nitz = 100
  INTEGER                           :: k, kmp, kmpm1, kssm1, ksg
  INTEGER                           :: ks, kk, n, ni2
  REAL                              :: scl, ascl
  REAL, ALLOCATABLE, DIMENSION(:)   :: wm, wp, w, scz, scz2
  LOGICAL :: dodebug=.false.
!
  kmp = kass + 2*kex
  kmpm1 = kmp - 1
  kssm1 = kass - 1
!
  allocate ( wm(kmp) )
  allocate ( wp(kmp) )
  allocate ( w(kmp) )
  allocate ( scz(kmp) )
  allocate ( scz2(kmp) )

  !STEVE:
  !if (dodebug) print *, "allocate (21), kmp = ", kmp

!
  if (dodebug) print *, "kmp = ", kmp
  if (dodebug) print *, "kass = ", kass
  if (dodebug) print *, "kex = ", kex
  if (dodebug) print *, "kssm1 = ", kssm1
  do k=1,kssm1
    if (dodebug) print *, "k,kssm1 = ", k,kssm1
    scl = vsclf*Grd%dzt(k)
    ascl = scl*scl / float(nitz);
    if (k .eq. 1) then
      wm(k+kex) = ascl / (2.0*Grd%dzw(k-1)*Grd%dzt(k))
    else
      wm(k+kex) = ascl / (Grd%dzw(k-1)*Grd%dzt(k))
    endif
    wp(k+kex) = ascl / (Grd%dzw(k)*Grd%dzt(k))
    w(k+kex) = 1.0 - wm(k+kex) - wp(k+kex)
  enddo
  do k=1,kex
    if (dodebug) print *, "k,kex = ", k,kex
    wm(k) = wm(kex+1)
    wp(k) = wp(kex+1)
    w(k) = w(kex+1)
  enddo
  do k=kass+kex,kmp
    if (dodebug) print *, "k,kmp = ", k,kmp
    wm(k) = wm(kex+kssm1)
    wp(k) = wp(kex+kssm1)
    w(k) = w(kex+kssm1)
  enddo
!
  ni2 = nitz / 2
!
  do ks=1,kass
    if (dodebug) print *, "ks,kass = ", ks,kass
    do k=1,kmp
      scz(k) = 0.0
    enddo
    ksg = ks + kex
    scz(ksg) = 1.0
    do n=1,ni2
      do k=2,kmpm1
        scz2(k) = w(k)*scz(k) + wm(k)*scz(k-1) + wp(k)*scz(k+1)
      enddo
      scz2(1) = scz2(2)
      scz2(kmp) = scz2(kmpm1)

      do k=2,kmpm1
       scz(k) = w(k)*scz2(k) + wp(k-1)*scz2(k-1) + wm(k+1)*scz2(k+1)
      enddo
      scz(1) = scz(2)
      scz(kmp) = scz(kmpm1)
    enddo
!
    ascl = scz(ksg)
    do k=1,kmp
      scz(k) = scz(k) / ascl
    enddo
!
    do k=1,kass
      cvn(ks,k) = scz(k+kex)
    enddo
  enddo

  do k=1,kass
    if (dodebug) print *, "k,kass = ", ks,kass
    cvn(k,k) = vcvn
    do kk=k+1,kass
      ascl = 0.5*(cvn(k,kk) + cvn(kk,k))*vcvn
      cvn(k,kk) = ascl
      cvn(kk,k) = ascl
    enddo
  enddo
!
  cvnsalt = cvn
!
  deallocate ( wm )
  deallocate ( wp )
  deallocate ( w )
  deallocate ( scz )
  deallocate ( scz2 )
!
  if (dodebug) print *, "DONE SUBROUTINE vertical_covariance"

  END SUBROUTINE vertical_covariance



END MODULE so_vertical_covariance
