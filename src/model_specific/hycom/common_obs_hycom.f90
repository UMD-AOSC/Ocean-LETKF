MODULE common_obs_oceanmodel
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   04/26/2011 Steve PENNY converted to OCEAN for use with MOM4. grep '(OCEAN)' for changes.
!   01/18/2015 Steve Penny converted for use with MOM6
!   06/05/2017 Jili Dong modified for use with HYCOM 
!
!=======================================================================
  USE common
  USE common_oceanmodel
  USE params_obs
  USE kdtree

  IMPLICIT NONE
  PUBLIC Trans_XtoY, phys2ijk, phys2ij,read_obs, get_nobs, read_obs2, write_obs2, itpl_2d, itpl_3d, monit_dep
  PUBLIC center_obs_coords

  PRIVATE
  TYPE(KD_ROOT), SAVE :: kdtree_root
  INTEGER, SAVE       :: initialized = 0
  REAL(r_size), SAVE  :: prev_lon, prev_lat      !STEVE: for checking if this longitude was last searched
  INTEGER             :: k_sought=1
  INTEGER             :: k_found

  INTEGER, SAVE, ALLOCATABLE      :: idx(:)      !STEVE: index of the observations that are found by kd_search_radius
  REAL(r_size), SAVE, ALLOCATABLE :: dist(:)     !STEVE: distance from the center grid point
  INTEGER, SAVE :: nn                            !STEVE: total number of local observations found by kd_search_radius

  ! For debugging kdtree:
  REAL(r_size), ALLOCATABLE, DIMENSION(:) :: lon2ij, lat2ij

CONTAINS


!-----------------------------------------------------------------------
! Transformation from model space to observation space (i.e. H-operator)
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rlev,v3d,v2d,yobs)        !(OCEAN)
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: iv3d_u, iv3d_v, iv3d_t, iv3d_s, iv2d_ssh, iv2d_eta, iv2d_sst, iv2d_sss
  USE vars_model,   ONLY: lat2d,phi0

  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rlev
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  REAL(r_size)  :: temp_obs(2,2)
  INTEGER :: i,j,k,n
  INTEGER :: intelm,obs_id_hycom

  intelm = NINT(elm)
  SELECT CASE (intelm)
  CASE(id_u_obs)  ! U
    obs_id_hycom=1
    CALL hycom_intrp(v3d(:,:,:,:),ri,rj,rlev,temp_obs,obs_id_hycom)
    CALL itpl_2d_4pts(temp_obs,ri,rj,yobs)
  CASE(id_v_obs)  ! V
    obs_id_hycom=2
    CALL hycom_intrp(v3d(:,:,:,:),ri,rj,rlev,temp_obs,obs_id_hycom)
    CALL itpl_2d_4pts(temp_obs,ri,rj,yobs)
  CASE(id_t_obs)  ! T
    obs_id_hycom=3
    CALL hycom_intrp(v3d(:,:,:,:),ri,rj,rlev,temp_obs,obs_id_hycom)
    CALL itpl_2d_4pts(temp_obs,ri,rj,yobs)
  CASE(id_s_obs)  ! S                             !(OCEAN)
    obs_id_hycom=4
    CALL hycom_intrp(v3d(:,:,:,:),ri,rj,rlev,temp_obs,obs_id_hycom)
    CALL itpl_2d_4pts(temp_obs,ri,rj,yobs)
  CASE(id_eta_obs) ! SSH                          !(OCEAN)
    !STEVE: use mom6 surface height to form Yb (i.e. hdxf)
    CALL itpl_2d(v2d(:,:,iv2d_ssh),ri,rj,yobs)    !(OCEAN)
  CASE(id_sst_obs) ! SST                          !(OCEAN)
    CALL itpl_2d(v3d(:,:,1,iv3d_t),ri,rj,yobs)    !(OCEAN)
  CASE(id_sss_obs) ! SSS                          !(OCEAN)
    CALL itpl_2d(v2d(:,:,iv2d_sss),ri,rj,yobs)    !(OCEAN)
  CASE DEFAULT
    print *, "ERROR::Trans_XtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_u_obs, id_v_obs, id_t_obs, id_s_obs, id_ssh_obs, id_eta_obs, id_sst_obs, id_sss_obs, id_x_obs, id_y_obs, id_z_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

END SUBROUTINE Trans_XtoY


!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(elem,rlon,rlat,rlev,ri,rj,rk)     !(OCEAN)
  USE params_model, ONLY: nlon, nlat, nlev
  USE vars_model,   ONLY: lon2d, lat2d, lev2d, lon, lat, lev
  USE vars_model,   ONLY: lon0, lonf, lat0, latf, wrapgap
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: aj,ak,ai, rrlon, glon,glat,xlon,xlat
! REAL(r_size) :: lnps(nlon,nlat)
  REAL(r_size) :: plev(nlev)
  INTEGER :: i,j,k,n, ni,nj, ii,jj
  LOGICAL, PARAMETER :: dodebug = .false.

  !-----------------------------------------------------------------------------
  ! Initialize the KD-tree for the model-grid
  !-----------------------------------------------------------------------------
  !STEVE: NOTE: this is replicated for each process
  if (initialized == 0) then
    WRITE(6,*) "Initializing the obs_local()"
    initialized = 1
    call kd_init( kdtree_root, RESHAPE(lon2d(:,:),(/nlon*nlat/)), RESHAPE(lat2d(:,:),(/nlon*nlat/)) )
    WRITE(6,*) "Done constructing KD search tree."
    WRITE(6,*) "nlon*nlat = ", nlon*nlat
    ALLOCATE(dist(k_sought))
    ALLOCATE( idx(k_sought))
    prev_lon = -1e10
    prev_lat = -1e10

    ! To debug kdtree:
    if (dodebug) then
      ALLOCATE(lon2ij(nlon*nlat),lat2ij(nlon*nlat))
      lon2ij = RESHAPE(lon2d(:,:),(/nlon*nlat/))
      lat2ij = RESHAPE(lat2d(:,:),(/nlon*nlat/))
    endif
  endif

  !-----------------------------------------------------------------------------
  ! query the KD-tree
  !-----------------------------------------------------------------------------
  if (rlon .ne. prev_lon .or. rlat .ne. prev_lat) then
    CALL kd_search_nnearest(kdtree_root,  (/rlon, rlat/), k_sought, idx, dist, k_found, .false.)
    prev_lon = rlon
    prev_lat = rlat
    if (dodebug) then
      WRITE(6,*) "common_obs_mom6.f90::phys2ijk ---------------------------- "
      WRITE(6,*) "base observation point (lon/lat) :: ", rlon, rlat
      WRITE(6,*) "k_sought = ", k_sought
      WRITE(6,*) "k_found  = ", k_found
      WRITE(6,*) "Grid points found by kd_search_nnearest = "
      WRITE(6,*) "idx,ni,nj,lon,lat,dist"
      do n=1,k_found
        ni = MODULO(idx(n)-1,nlon)+1
        nj = idx(n)/nlon+1 ! Integer division
        WRITE(6,*) idx(n), ni,nj, lon2d(ni,nj), lat2d(ni,nj), dist(n)
      enddo
    endif
  endif

  !-----------------------------------------------------------------------------
  ! rlon -> ri, rlat -> rj
  !-----------------------------------------------------------------------------

  ! Find the appropriate grid box containing this observation

  ! First, find the index of the nearest grid point
  ni = MODULO(idx(1)-1,nlon)+1
  nj = FLOOR(REAL(idx(1))/REAL(nlon))+1 ! Integer division
  glon = lon2d(ni,nj)
  glat = lat2d(ni,nj)

  ! For debugging:
  if (dodebug) then
    xlon = lon2ij(idx(1))
    xlat = lat2ij(idx(1))
  endif

  !STEVE: initialize to something that will throw errors if it's not changed within
  ri = -1
  rj = -1
  rk = -1

  ! Second, check to see whether the observation is to the left/right, up/down from this point

  !-----------------------------------------------------------------------------
  ! Get longitude scaling
  !-----------------------------------------------------------------------------
  if (glon == rlon) then
    ai = 0.0
  elseif (glon < rlon) then
    ai = (rlon - lon2d(ni,nj)) / (lon2d(ni+1,nj) - lon2d(ni,nj))
  elseif (rlon < glon) then
    ai = (rlon - lon2d(ni,nj)) / (lon2d(ni,nj) - lon2d(ni-1,nj))
  endif
  ri = REAL(ni,r_size) + ai

  !-----------------------------------------------------------------------------
  ! Get latitude scaling
  !-----------------------------------------------------------------------------
  if (glat == rlat) then
    aj = 0.0
  elseif (glat < rlat) then
    aj = (rlat - lat2d(ni,nj)) / (lat2d(ni,nj+1) - lat2d(ni,nj))
  elseif (rlat < glat) then
    aj = (rlat - lat2d(ni,nj)) / (lat2d(ni,nj) - lat2d(ni,nj-1))
    !STEVE: should be negative
  endif
  rj = REAL(nj,r_size) + aj

  ! Set flags to mark rlon outside of map range:
  if (ri < 1) then
    ri = 0
  elseif(ri >= nlon+1) then
    ri = nlon+1
  endif

  if (rj < 1) then
    rj = 0
  elseif (rj >= nlat+1) then
    rj = nlat+1
  endif

  ! if (dodebug .and. (ri > nlon .or. ri < 1)) then
  if (dodebug) then
    WRITE(6,*) "================="
    WRITE(6,*) "In common_obs_mom6.f90::phys2ijk,"
    WRITE(6,*) "-----------------"
    WRITE(6,*) "rlon = ", rlon
    if (dodebug) WRITE(6,*) "xlon = ", xlon
    WRITE(6,*) "glon = ", glon
    WRITE(6,*) "ni = ", ni
    WRITE(6,*) "ai = ", ai
    WRITE(6,*) "ri = ", ri
    WRITE(6,*) "lon0 = ", lon0
    WRITE(6,*) "lonf = ", lonf
    WRITE(6,*) "-----------------"
    WRITE(6,*) "rlat = ", rlat
    if (dodebug) WRITE(6,*) "xlat = ", xlat
    WRITE(6,*) "glat = ", glat
    WRITE(6,*) "nj = ", nj
    WRITE(6,*) "aj = ", aj
    WRITE(6,*) "rj = ", rj
    WRITE(6,*) "lat0 = ", lat0
    WRITE(6,*) "latf = ", latf
    WRITE(6,*) "================="
  endif

  if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) RETURN
  if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) RETURN

  !-----------------------------------------------------------------------------
  ! rlev -> rk
  !-----------------------------------------------------------------------------
  if (NINT(elem) == id_ssh_obs) then     ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif(NINT(elem) == id_eta_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif(NINT(elem) == id_sst_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif(NINT(elem) == id_sss_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  else
    !
    ! vertical interpolation
    !
    !
    ! find rk
    !
!   if (dodebug) WRITE(6,*) "rlev = ", rlev
    do k=1,nlev
!     if (dodebug) WRITE(6,*) "lev(",k,") = ", lev(k)
      if (rlev < lev(k)) EXIT
      if (k .eq. nlev .and. rlev .eq. lev(nlev)) EXIT !STEVE: added this case for simulated obs that reach the lowest model levels.
                                                      !       Otherwise, k iterates to nlev+1 before exiting loop.
      if (k>1 .and. lev(k)==0.0) then
        WRITE(6,*) "common_obs_mom4.f90::phys2ijk:: lev(k)==0.0 but k>1. ERROR, EXITING..."
        WRITE(6,*) "(model levels array lev() was not appropriately set)"
        STOP (934)
      endif
    enddo

    if (rk > nlev) then
      WRITE(6,*) "common_obs_mom4.f90::phys2ijk:: rk>nlev. ERROR, EXITING..."
      WRITE(6,*) "(appropriate model leve not found for rlev = ",rlev,")"
      STOP (935)
    endif

    if (k .eq. 1) then
      rk = 1 !ak
    else
      !STEVE: now apply the interpolation at the identified model level:
      ak = (rlev - lev(k-1)) / (lev(k) - lev(k-1))
      rk = REAL(k-1,r_size) + ak
    endif
    if (dodebug) WRITE(6,*) "phys2ijk:: rlon,rlat,rlev,ri,rj,rk = ", rlon,rlat,rlev,ri,rj,rk

  endif

  if (dodebug) then
    WRITE(6,*) "-----------------"
    WRITE(6,*) "rlev = ", rlev
    WRITE(6,*) "k  = ", k
    WRITE(6,*) "lev(k-1) = ", lev(k-1)
    WRITE(6,*) "lev(k)   = ", lev(k)
    WRITE(6,*) "ak = ", ak
    WRITE(6,*) "rk = ", rk
    WRITE(6,*) "lev0 = ", lev(1)
    WRITE(6,*) "levf = ", lev(nlev)
    WRITE(6,*) "================="
  endif 
  
END SUBROUTINE phys2ijk

!-----------------------------------------------------------------------
! Coordinate conversion just for 2D
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(elem,rlon,rlat,rlev,ri,rj)     !(OCEAN)
  USE params_model, ONLY: nlon, nlat
  USE vars_model,   ONLY: lon2d, lat2d, lev2d, lon, lat
  USE vars_model,   ONLY: lon0, lonf, lat0, latf, wrapgap
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size) :: aj,ai, rrlon, glon,glat,xlon,xlat
! REAL(r_size) :: lnps(nlon,nlat)
  INTEGER :: i,j,n, ni,nj, ii,jj
  LOGICAL, PARAMETER :: dodebug = .false.

  !-----------------------------------------------------------------------------
  ! Initialize the KD-tree for the model-grid
  !-----------------------------------------------------------------------------
  !STEVE: NOTE: this is replicated for each process
  if (initialized == 0) then
    WRITE(6,*) "Initializing the obs_local()"
    initialized = 1
    call kd_init( kdtree_root, RESHAPE(lon2d(:,:),(/nlon*nlat/)), RESHAPE(lat2d(:,:),(/nlon*nlat/)) )
    WRITE(6,*) "Done constructing KD search tree."
    WRITE(6,*) "nlon*nlat = ", nlon*nlat
    ALLOCATE(dist(k_sought))
    ALLOCATE( idx(k_sought))
    prev_lon = -1e10
    prev_lat = -1e10

    ! To debug kdtree:
    if (dodebug) then
      ALLOCATE(lon2ij(nlon*nlat),lat2ij(nlon*nlat))
      lon2ij = RESHAPE(lon2d(:,:),(/nlon*nlat/))
      lat2ij = RESHAPE(lat2d(:,:),(/nlon*nlat/))
    endif
  endif

  !-----------------------------------------------------------------------------
  ! query the KD-tree
  !-----------------------------------------------------------------------------
  if (rlon .ne. prev_lon .or. rlat .ne. prev_lat) then
    CALL kd_search_nnearest(kdtree_root,  (/rlon, rlat/), k_sought, idx, dist, k_found, .false.)
    prev_lon = rlon
    prev_lat = rlat
    if (dodebug) then
      WRITE(6,*) "common_obs_mom6.f90::phys2ijk ---------------------------- "
      WRITE(6,*) "base observation point (lon/lat) :: ", rlon, rlat
      WRITE(6,*) "k_sought = ", k_sought
      WRITE(6,*) "k_found  = ", k_found
      WRITE(6,*) "Grid points found by kd_search_nnearest = "
      WRITE(6,*) "idx,ni,nj,lon,lat,dist"
      do n=1,k_found
        ni = MODULO(idx(n)-1,nlon)+1
        nj = idx(n)/nlon+1 ! Integer division
        WRITE(6,*) idx(n), ni,nj, lon2d(ni,nj), lat2d(ni,nj), dist(n)
      enddo
    endif
  endif

  !-----------------------------------------------------------------------------
  ! rlon -> ri, rlat -> rj
  !-----------------------------------------------------------------------------

  ! Find the appropriate grid box containing this observation

  ! First, find the index of the nearest grid point
  ni = MODULO(idx(1)-1,nlon)+1
  nj = FLOOR(REAL(idx(1))/REAL(nlon))+1 ! Integer division
  glon = lon2d(ni,nj)
  glat = lat2d(ni,nj)

  ! For debugging:
  if (dodebug) then
    xlon = lon2ij(idx(1))
    xlat = lat2ij(idx(1))
  endif

  !STEVE: initialize to something that will throw errors if it's not changed within
  ri = -1
  rj = -1

  ! Second, check to see whether the observation is to the left/right, up/down from this point

  !-----------------------------------------------------------------------------
  ! Get longitude scaling
  !-----------------------------------------------------------------------------
  if (glon == rlon) then
    ai = 0.0
  elseif (glon < rlon) then
    ai = (rlon - lon2d(ni,nj)) / (lon2d(ni+1,nj) - lon2d(ni,nj))
  elseif (rlon < glon) then
    ai = (rlon - lon2d(ni,nj)) / (lon2d(ni,nj) - lon2d(ni-1,nj))
  endif
  ri = REAL(ni,r_size) + ai

  !-----------------------------------------------------------------------------
  ! Get latitude scaling
  !-----------------------------------------------------------------------------
  if (glat == rlat) then
    aj = 0.0
  elseif (glat < rlat) then
    aj = (rlat - lat2d(ni,nj)) / (lat2d(ni,nj+1) - lat2d(ni,nj))
  elseif (rlat < glat) then
    aj = (rlat - lat2d(ni,nj)) / (lat2d(ni,nj) - lat2d(ni,nj-1))
    !STEVE: should be negative
  endif
  rj = REAL(nj,r_size) + aj

  ! Set flags to mark rlon outside of map range:
  if (ri < 1) then
    ri = 0
  elseif(ri >= nlon+1) then
    ri = nlon+1
  endif

  if (rj < 1) then
    rj = 0
  elseif (rj >= nlat+1) then
    rj = nlat+1
  endif

  ! if (dodebug .and. (ri > nlon .or. ri < 1)) then
  if (dodebug) then
    WRITE(6,*) "================="
    WRITE(6,*) "In common_obs_hycom.f90::phys2ij,"
    WRITE(6,*) "-----------------"
    WRITE(6,*) "rlon = ", rlon
    if (dodebug) WRITE(6,*) "xlon = ", xlon
    WRITE(6,*) "glon = ", glon
    WRITE(6,*) "ni = ", ni
    WRITE(6,*) "ai = ", ai
    WRITE(6,*) "ri = ", ri
    WRITE(6,*) "lon0 = ", lon0
    WRITE(6,*) "lonf = ", lonf
    WRITE(6,*) "-----------------"
    WRITE(6,*) "rlat = ", rlat
    if (dodebug) WRITE(6,*) "xlat = ", xlat
    WRITE(6,*) "glat = ", glat
    WRITE(6,*) "nj = ", nj
    WRITE(6,*) "aj = ", aj
    WRITE(6,*) "rj = ", rj
    WRITE(6,*) "lat0 = ", lat0
    WRITE(6,*) "latf = ", latf
    WRITE(6,*) "================="
  endif

  if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) RETURN
  if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) RETURN


  
END SUBROUTINE phys2ij

!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  USE params_model, ONLY: nlon, nlat
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

!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d_4pts(var,ri,rj,var5)
  USE params_model, ONLY: nlon, nlat
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(2,2)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

    var5 = var(1,1) * (1-ai) * (1-aj) &
       & + var(2  ,1) *    ai  * (1-aj) &
       & + var(1,2  ) * (1-ai) *    aj  &
       & + var(2  ,2  ) *    ai  *    aj

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


END SUBROUTINE itpl_2d_4pts




SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  USE params_model, ONLY: nlon, nlat, nlev
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


!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_s,rmse_ssh,rmse_eta,rmse_sst,rmse_sss !(OCEAN)
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_ssh,bias_eta,bias_sst,bias_sss !(OCEAN)
  INTEGER :: n,iu,iv,it,is,issh,ieta,isst,isss                                !(OCEAN)

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_s = 0.0d0    !(OCEAN)
  rmse_ssh = 0.0d0  !(OCEAN)
  rmse_eta = 0.0d0  !(OCEAN)
  rmse_sst = 0.0d0  !(OCEAN)
  rmse_sss = 0.0d0  !(OCEAN)
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_s = 0.0d0    !(OCEAN)
  bias_ssh = 0.0d0  !(OCEAN)
  bias_eta = 0.0d0  !(OCEAN)
  bias_sst = 0.0d0  !(OCEAN)
  bias_sss = 0.0d0  !(OCEAN)
  iu = 0
  iv = 0
  it = 0
  is = 0            !(OCEAN)
  issh = 0          !(OCEAN)
  ieta = 0          !(OCEAN)
  isst = 0          !(OCEAN)
  isss = 0          !(OCEAN)
  DO n=1,nn
    if (qc(n) /= 1) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep(n)**2
      bias_u = bias_u + dep(n)
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep(n)**2
      bias_v = bias_v + dep(n)
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep(n)**2
      bias_t = bias_t + dep(n)
      it = it + 1
    CASE(id_s_obs)                    !(OCEAN)
      rmse_s = rmse_s + dep(n)**2     !(OCEAN)
      bias_s = bias_s + dep(n)        !(OCEAN)
      is = is + 1                     !(OCEAN)
    CASE(id_ssh_obs)                  !(OCEAN)
      rmse_ssh = rmse_ssh + dep(n)**2 !(OCEAN)
      bias_ssh = bias_ssh + dep(n)    !(OCEAN)
      issh = issh + 1                 !(OCEAN)
    CASE(id_eta_obs)                  !(OCEAN)
      rmse_eta = rmse_eta + dep(n)**2 !(OCEAN)
      bias_eta = bias_eta + dep(n)    !(OCEAN)
      ieta = ieta + 1                 !(OCEAN)
    CASE(id_sst_obs)                  !(OCEAN)
      rmse_sst = rmse_sst + dep(n)**2 !(OCEAN)
      bias_sst = bias_sst + dep(n)    !(OCEAN)
      isst = isst + 1                 !(OCEAN)
    CASE(id_sss_obs)                  !(OCEAN)
      rmse_sss = rmse_sss + dep(n)**2 !(OCEAN)
      bias_sss = bias_sss + dep(n)    !(OCEAN)
      isss = isss + 1                 !(OCEAN)
    END SELECT
  enddo
  if (iu == 0) then
    rmse_u = undef
    bias_u = undef
  else
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  endif
  if (iv == 0) then
    rmse_v = undef
    bias_v = undef
  else
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  endif
  if (it == 0) then
    rmse_t = undef
    bias_t = undef
  else
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  endif
  if (is == 0) then                                   !(OCEAN)
    rmse_s = undef                                   !(OCEAN)
    bias_s = undef                                   !(OCEAN)
  else
    rmse_s = SQRT(rmse_s / REAL(is,r_size))          !(OCEAN)
    bias_s = bias_s / REAL(is,r_size)                !(OCEAN)
  endif
  if (issh == 0) then                                 !(OCEAN)
    rmse_ssh = undef                                 !(OCEAN)
    bias_ssh = undef                                 !(OCEAN)
  else
    rmse_ssh = SQRT(rmse_ssh / REAL(issh,r_size))    !(OCEAN)
    bias_ssh = bias_ssh / REAL(issh,r_size)          !(OCEAN)
  endif
  if (ieta == 0) then                                 !(OCEAN)
    rmse_eta = undef                                 !(OCEAN)
    bias_eta = undef                                 !(OCEAN)
  else
    rmse_eta = SQRT(rmse_eta / REAL(ieta,r_size))    !(OCEAN)
    bias_eta = bias_eta / REAL(ieta,r_size)          !(OCEAN)
  endif
  if (isst == 0) then                                 !(OCEAN)
    rmse_sst = undef                                 !(OCEAN)
    bias_sst = undef                                 !(OCEAN)
  else
    rmse_sst = SQRT(rmse_sst / REAL(isst,r_size))    !(OCEAN)
    bias_sst = bias_sst / REAL(isst,r_size)          !(OCEAN)
  endif
  if (isss == 0) then                                 !(OCEAN)
    rmse_sss = undef                                 !(OCEAN)
    bias_sss = undef                                 !(OCEAN)
  else
    rmse_sss = SQRT(rmse_sss / REAL(isss,r_size))    !(OCEAN)
    bias_sss = bias_sss / REAL(isss,r_size)          !(OCEAN)
  endif

  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ============================================='
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','eta','SST','SSS'                                       !(OCEAN)
  WRITE(6,'(7ES12.3)') bias_u,bias_v,bias_t,bias_s,bias_ssh,bias_eta,bias_sst,bias_sss               !(OCEAN)
  WRITE(6,'(7ES12.3)') rmse_u,rmse_v,rmse_t,rmse_s,rmse_ssh,bias_eta,rmse_sst,bias_sss               !(OCEAN)
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ============================'
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','eta','SST','SSS'                                       !(OCEAN)
  WRITE(6,'(7I12)') iu,iv,it,is,issh,ieta,isst,isss                                              !(OCEAN)
  WRITE(6,'(A)') '========================================================================'

END SUBROUTINE monit_dep


!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:) 
  INTEGER :: ios
  INTEGER :: iu,iv,it,is,issh,ieta,isst,isss,ix,iy,iz !(OCEAN)
  INTEGER :: nprof !(OCEAN)
  REAL(r_sngl) :: lon_m1, lat_m1
  INTEGER :: iunit
  LOGICAL :: ex
  LOGICAL, PARAMETER :: dodebug=.false.

  ALLOCATE(wk(nrec))
  nn = 0
  iu = 0
  iv = 0
  it = 0
  is = 0    !(OCEAN)
  issh = 0  !(OCEAN)
  ieta = 0  !(OCEAN)
  isst = 0  !(OCEAN)
  isss = 0  !(OCEAN)
  ix = 0    !(OCEAN)
  iy = 0    !(OCEAN)
  iz = 0    !(OCEAN)
  nprof = 0    !(OCEAN)
  lon_m1 = 0    !(OCEAN)
  lat_m1 = 0    !(OCEAN)
  iunit=91
  if (dodebug) print *, "get_nobs::"
  INQUIRE(FILE=cfile,EXIST=ex)
  if (ex) then
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      if (dodebug .and. nrec.eq.6) then
        PRINT '(I6,2F7.2,F10.2,2F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
      elseif (dodebug .and. nrec .eq. 8) then
        PRINT '(I6,2F7.2,F10.2,4F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8)
      elseif (dodebug .and. nrec .eq. 9) then
        PRINT '(I6,2F7.2,F10.2,5F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
      endif
      if (ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        iu = iu + 1
      CASE(id_v_obs)
        iv = iv + 1
      CASE(id_t_obs)
        it = it + 1
      CASE(id_s_obs)     !(OCEAN)
        is = is + 1      !(OCEAN)
      CASE(id_ssh_obs)   !(OCEAN)
        issh = issh + 1  !(OCEAN)
      CASE(id_eta_obs)   !(OCEAN)
        ieta = ieta + 1  !(OCEAN)
      CASE(id_sst_obs)   !(OCEAN)
        isst = isst + 1  !(OCEAN)
      CASE(id_sss_obs)   !(OCEAN)
        isss = isss + 1  !(OCEAN)
      CASE(id_x_obs)     !(OCEAN)
        ix = ix + 1      !(OCEAN)
      CASE(id_y_obs)     !(OCEAN)
        iy = iy + 1      !(OCEAN)
      CASE(id_z_obs)     !(OCEAN)
        iz = iz + 1      !(OCEAN)
      END SELECT
      if (wk(2) .ne. lon_m1 .and. wk(3) .ne. lat_m1) then
        nprof=nprof+1
        lon_m1 = wk(2)
        lat_m1 = wk(3)
      endif
      nn = nn + 1
    enddo
    WRITE(6,'(I10,A)') nprof ,' PROFILES INPUT'
    WRITE(6,'(I10,A)') nn,' TOTAL OBSERVATIONS INPUT (in get_nobs)'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '          S:',is   !(OCEAN)
    WRITE(6,'(A12,I10)') '        SSH:',issh !(OCEAN)
    WRITE(6,'(A12,I10)') '        eta:',ieta !(OCEAN)
    WRITE(6,'(A12,I10)') '        SST:',isst !(OCEAN)
    WRITE(6,'(A12,I10)') '        SSS:',isss !(OCEAN)
    WRITE(6,'(A12,I10)') '          X:',ix   !(OCEAN)
    WRITE(6,'(A12,I10)') '          Y:',iy   !(OCEAN)
    WRITE(6,'(A12,I10)') '          Z:',iz   !(OCEAN)
    CLOSE(iunit)
  else
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  endif
  DEALLOCATE(wk)

  if (nn .eq. 0) then
    WRITE(6,*) "get_nobs:: No observations have been found. Exiting..."
    !STOP(60)
  endif

END SUBROUTINE get_nobs


SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,obhr)
  USE vars_model, ONLY: lon0, lonf, lat0, latf, wrapgap
! USE vars_model, ONLY: lon, lat, lev
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: obhr(nn)
! REAL(r_sngl) :: wk(6)
  REAL(r_sngl) :: wk(9)
  !REAL(r_size) :: wk(6) !(OCEAN) STEVE: I changed this because the netcdf observation files are stored as double precision
  INTEGER :: n,iunit
  ! STEVE: for general grid
  LOGICAL, PARAMETER :: dodebug = .false.
  LOGICAL :: process_obs = .true.

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    READ(iunit) wk
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    obhr(n) = REAL(wk(7),r_size)

    !STEVE: error check
    if (oerr(n) .le. 0) then
      print *, "common_obs_mom6.f90::read_obs:: WARNING!"
      print *, "STEVE: oerr <= 0, must be > 0 ..." 
      print *, "STEVE: oerr(",n,") = ", oerr(n)
      PRINT '(I6,2F7.2,F10.2,2F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
!     stop 9
    endif
    
    ! Special processing for obs:
    !STEVE: this handles the fact that the observations are typically on a 0 to 360º grid, while
    !       the NCEP mom4p1 grid configuration is on a ~ -285 to 75º grid
    if (process_obs) then
      CALL center_obs_coords(rlon,oerr,nn)
    endif 

  enddo

  CLOSE(iunit)

  if (MAXVAL(rlon) > lonf) then
    WRITE(6,*) "read_obs:: Error: MAX(observation lon, i.e. rlon) > lonf"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lonf = ", lonf
    STOP (2)
  endif
  if (MINVAL(rlon) < lon0) then
    WRITE(6,*) "read_obs:: Error: MIN(observation lon, i.e. rlon) < lon0"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lon0 = ", lon0
    STOP (2)
  endif

END SUBROUTINE read_obs


SUBROUTINE read_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: ohx(nn)
  REAL(r_size),INTENT(OUT) :: obhr(nn)
  INTEGER,INTENT(OUT) :: oqc(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  do n=1,nn
    READ(iunit) wk
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    ohx(n)  = REAL(wk(7),r_size)
    oqc(n)  = NINT(wk(8))
    obhr(n) = REAL(wk(9),r_size)
  enddo

  CLOSE(iunit)

END SUBROUTINE read_obs2


!STEVE: adding for support
SUBROUTINE write_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,obhr)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: odat(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: obhr(nn)
  REAL(r_sngl) :: wk(9)
! REAL(r_size) :: wk(6) !(OCEAN) STEVE: I changed this because the netcdf observation files are stored as double precision
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  do n=1,nn
    wk(1) = REAL(elem(n),r_sngl)
    wk(2) = REAL(rlon(n),r_sngl)
    wk(3) = REAL(rlat(n),r_sngl)
    wk(4) = REAL(rlev(n),r_sngl)
    wk(5) = REAL(odat(n),r_sngl)
    wk(6) = REAL(oerr(n),r_sngl)
    wk(7) = REAL(obhr(n),r_sngl)
    wk(8) = -1
    wk(9) = -1
    WRITE(iunit) wk
  enddo

  CLOSE(iunit)

END SUBROUTINE write_obs


SUBROUTINE write_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,qcflag_in)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: odat(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: ohx(nn)
  REAL(r_size),INTENT(IN) :: obhr(nn)
  LOGICAL, INTENT(IN), OPTIONAL :: qcflag_in
  LOGICAL :: qcflag
  INTEGER,INTENT(IN) :: oqc(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit
  LOGICAL, PARAMETER :: dodebug=.false.

  if (PRESENT(qcflag_in)) then
    qcflag = qcflag_in
  else
    qcflag = .false.
  endif

  iunit=92
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  do n=1,nn
    if (qcflag .and. oqc(n)==0) CYCLE
    wk(1) = REAL(elem(n),r_sngl)  ! ID for observation type
    wk(2) = REAL(rlon(n),r_sngl)  ! Ob lon
    wk(3) = REAL(rlat(n),r_sngl)  ! Ob lat
    wk(4) = REAL(rlev(n),r_sngl)  ! Ob level
    wk(5) = REAL(odat(n),r_sngl)  ! Observed data quantity
    wk(6) = REAL(oerr(n),r_sngl)  ! Estimated observation error
    wk(7) = REAL(ohx(n),r_sngl)   ! Model forecast transformed to observation space: H(xb)
    wk(8) = REAL(oqc(n),r_sngl)   ! Quality control ID (1==keep, 0==discard) for use in assimilation
    wk(9) = REAL(obhr(n),r_sngl)  ! Quality control ID (1==keep, 0==discard) for use in assimilation
    if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
    WRITE(iunit) wk
  enddo

  CLOSE(iunit)

END SUBROUTINE write_obs2


SUBROUTINE center_obs_coords(rlon,oerr,nn)
!===============================================================================
! Center all observations within the longitudes defined by the model grid
!===============================================================================
  USE vars_model,   ONLY: lon0, lonf, wrapgap
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: rlon(nn)
  REAL(r_size),INTENT(INOUT) :: oerr(nn)
  INTEGER, INTENT(IN) :: nn
  INTEGER :: n
  LOGICAL :: dodebug = .false.
  LOGICAL :: dodebug1 = .false.

  if (dodebug) then
    WRITE(6,*) "center_obs_coords::"
    WRITE(6,*) "wrapgap = ", wrapgap
  endif

  do n=1,nn

    if (rlon(n) >= lonf) then
      if (dodebug1) WRITE(6,*) "(a)"

      if (dodebug) then
        WRITE(6,*) "letkf_obs.f90:: Wrapping large lon obs to model grid: n = ", n
        WRITE(6,*) "pre : rlon(n) = ", rlon(n)
        WRITE(6,*) "360 - rlon(n) = ", 360.0 - rlon(n)
        WRITE(6,*) "abs(rlon(n) - lonf) = ", abs(rlon(n) - lonf)
      endif

      ! Update the coordinate
      if ( abs(rlon(n) - lonf) < wrapgap ) then
        if (dodebug1) WRITE(6,*) "(b)"
          ! First, handle observations that are just outside of the model grid
          !STEVE: shift it if it's just outside grid

        if (abs(rlon(n) - lonf) < wrapgap/2) then
            rlon(n) = lonf
        else
            rlon(n) = lon0
        endif

        ! Increase error to compensate
        oerr(n) = oerr(n)*2
      else
        if (dodebug1) WRITE(6,*) "(c)"
        ! Otherwise, wrap the observation coordinate to be inside of the defined model grid coordinates
        !Wrap the coordinate
        if (dodebug) then
          WRITE(6,*) "lon0    = ", lon0
          WRITE(6,*) "lonf    = ", lonf
          WRITE(6,*) "rlon(n) = ", rlon(n)
          WRITE(6,*) "wrapgap = ", wrapgap
        endif
        rlon(n) = REAL(lon0 + abs(rlon(n) - lonf) - wrapgap,r_size)
        if (dodebug) then
          WRITE(6,*) "After update:"
          WRITE(6,*) "rlon(n) = REAL(lon0 + abs(rlon(n) - lonf) - wrapgap,r_size)"
          WRITE(6,*) "rlon(n) = ", rlon(n)
        endif
      endif

      if (dodebug) WRITE(6,*) "post : rlon(n) = ", rlon(n)

    elseif (rlon(n) < lon0) then
      if (dodebug1) WRITE(6,*) "(d)"

      if (dodebug) then
          WRITE(6,*) "letkf_obs.f90:: Wrapping small lon obs to model grid: n = ", n
          WRITE(6,*) "pre  : rlon(n) = ", rlon(n)
          WRITE(6,*) "360 - rlon(n) = ", 360.0 - rlon(n)
      endif

      ! Update the coordinate
      if (abs(lon0 - rlon(n)) < wrapgap ) then
        if (dodebug1) WRITE(6,*) "(e)"
        !STEVE: shift it if it's just outside grid
        if (abs(lon0 - rlon(n)) < wrapgap/2) then
            rlon(n) = lon0
        else
            rlon(n) = lonf
        endif
        ! Increase error to compensate                                 
        oerr(n) = oerr(n)*2
      else
        if (dodebug1) WRITE(6,*) "(f)"
        !Wrap the coordinate
        rlon(n) = REAL(lonf - abs(lon0 - rlon(n)) + wrapgap,r_size)
      endif

      if (dodebug) WRITE(6,*) "post : rlon(n) = ", rlon(n)

    endif
  enddo

  WRITE(6,*) "center_obs_coords:: finished recentering observations."

  if (MAXVAL(rlon) > lonf) then
    WRITE(6,*) "read_obs:: Error: MAX(observation lon, i.e. rlon) > lonf"
    WRITE(6,*) "MAXVAL(rlon) = ", MAXVAL(rlon)
    WRITE(6,*) "lonf = ", lonf
!   WRITE(6,*) "lon(nlon) = ", lon(nlon)
    STOP (22)
  endif
  if (MINVAL(rlon) < lon0) then
    WRITE(6,*) "read_obs:: Error: MIN(observation lon, i.e. rlon) < lon0"
    WRITE(6,*) "MINVAL(rlon) = ", MINVAL(rlon)
    WRITE(6,*) "lon0 = ", lon0
!   WRITE(6,*) "lon(1) = ", lon(1)
    STOP (23)
  endif

END SUBROUTINE center_obs_coords


END MODULE common_obs_oceanmodel
