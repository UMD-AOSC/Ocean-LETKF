MODULE H_operator
!=======================================================================
!
! [PURPOSE:]  Linear H Operator
!
! [HISTORY:]
!   09/04/2016 Steve Penny created
!
!=======================================================================
  USE common
  USE common_oceanmodel
  USE params_obs

  IMPLICIT NONE
  PUBLIC H_XtoY, H_YtoX

CONTAINS


!-----------------------------------------------------------------------
! Transformation from model space to observation space (i.e. H-operator)
!-----------------------------------------------------------------------
SUBROUTINE H_XtoY(elm,ri,rj,rk,v3d,v2d,yobs)        !(OCEAN)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_eta, iv2d_sst, iv2d_sss
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER :: i,j,k,n
  INTEGER :: intelm

  intelm = NINT(elm)
  SELECT CASE (intelm)
  CASE(id_u_obs)  ! U
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ri,rj,rk,yobs)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,rj,rk,yobs)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,yobs)
  CASE(id_s_obs)  ! S                             !(OCEAN)
    CALL itpl_3d(v3d(:,:,:,iv3d_s),ri,rj,rk,yobs) !(OCEAN)
  CASE(id_eta_obs) ! SSH                          !(OCEAN)
    !STEVE: use mom6 surface height to form Yb (i.e. hdxf)
    CALL itpl_2d(v2d(:,:,iv2d_eta),ri,rj,yobs)    !(OCEAN)
  CASE(id_sst_obs) ! SST                          !(OCEAN)
    CALL itpl_2d(v2d(:,:,iv2d_sst),ri,rj,yobs)    !(OCEAN)
  CASE(id_sss_obs) ! SSS                          !(OCEAN)
    CALL itpl_2d(v2d(:,:,iv2d_sss),ri,rj,yobs)    !(OCEAN)
  CASE DEFAULT
    print *, "ERROR::H_XtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_u_obs, id_v_obs, id_t_obs, id_s_obs, id_ssh_obs, id_eta_obs, id_sst_obs, id_sss_obs, id_x_obs, id_y_obs, id_z_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

END SUBROUTINE H_XtoY

!-----------------------------------------------------------------------
! Transformation from model space to observation space (i.e. H-operator)
!-----------------------------------------------------------------------
!SUBROUTINE H_YtoX(elm,ri,rj,rk,yval,ci,cj,ck,xval)        !(OCEAN)
SUBROUTINE H_YtoX(elm,ri,rj,rk,yval,xval)        !(OCEAN)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_eta, iv2d_sst, iv2d_sss
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(in) :: yval
! REAL(r_size),DIMENSION(2),INTENT(OUT) :: ci,cj,ck
  REAL(r_size),INTENT(OUT) :: xval
  INTEGER :: i,j,k,n
  INTEGER :: intelm

  intelm = NINT(elm)
  SELECT CASE (intelm)
  CASE(id_u_obs)  ! U
    CALL itpl_3d(yval,ri,rj,rk,xval)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(yval,ri,rj,rk,xval)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(yval,ri,rj,rk,xval)
  CASE(id_s_obs)  ! S                             !(OCEAN)
    CALL itpl_3d(yval,ri,rj,rk,xval)
  CASE(id_eta_obs) ! SSH                          !(OCEAN)
    !STEVE: use mom6 surface height to form Yb (i.e. hdxf)
    CALL itpl_2d(yval,ri,rj,xval)
  CASE(id_sst_obs) ! SST                          !(OCEAN)
    CALL itpl_2d(yval,ri,rj,xval)
  CASE(id_sss_obs) ! SSS                          !(OCEAN)
    CALL itpl_2d(yval,ri,rj,xval)
  CASE DEFAULT
    print *, "ERROR::H_XtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_u_obs, id_v_obs, id_t_obs, id_s_obs, id_ssh_obs, id_eta_obs, id_sst_obs, id_sss_obs, id_x_obs, id_y_obs, id_z_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

END SUBROUTINE H_YtoX

!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  if (i <= nlon) then
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj

!   if (.true.) then
!     print *, "common_obs_mom6.f90::itpl_2d::"
!     print *, "var(i-1,j-1) = ", var(i-1,j-1)
!     print *, "(1-ai) = ", (1-ai)
!     print *, "(1-aj) = ", (1-aj)
!     print *, "var(i  ,j-1) = ", var(i  ,j-1)
!     print *, "ai = ", ai
!     print *, "var(i-1,j  ) = ", var(i-1,j  )
!     print *, "aj = ", aj
!     print *, "var(i  ,j  ) = ", var(i  ,j  )
!     print *, "var5 = ", var5
!   endif

  else
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  endif

END SUBROUTINE itpl_2d


SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k
  LOGICAL, PARAMETER :: dodebug = .false.

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  if (dodebug) WRITE(6,*) "i,j,k,ai,aj,ak = ", i,j,k,ai,aj,ak

  if (i <= nlon) then
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  else
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  endif

  if (dodebug) WRITE(6,*) "var5 = ", var5

END SUBROUTINE itpl_3d



END MODULE H_operator
