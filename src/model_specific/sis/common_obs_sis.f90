MODULE common_obs_oceanmodel
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   04/26/2011 Steve PENNY converted to OCEAN for use with MOM4. grep '(OCEAN)' for changes.
!   01/18/2015 Steve Penny converted for use with MOM6
!   04/06/2015 Steve Penny converted for use with SIS
!
!=======================================================================
  USE common
  USE common_oceanmodel
! USE params_obs
  USE kdtree

  IMPLICIT NONE

  PUBLIC Trans_XtoY, phys2ijk, read_obs, get_nobs, read_obs2, write_obs2, itpl_2d, itpl_3d, monit_dep
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
SUBROUTINE Trans_XtoY(elm,ri,rj,rk,v3d,v2d,yobs)        !(OCEAN)
  USE params_model, ONLY: iv3d_hs, iv3d_hi, iv3d_t1, iv3d_t2, iv3d_ps, iv2d_cn, iv2d_ui, iv2d_vi
  USE params_model, ONLY: nlon,nlat,nlev,nv3d,nv2d
  USE params_obs,   ONLY: id_hs_obs, id_hi_obs, id_t1_obs, id_t2_obs, id_cn_obs, id_ui_obs, id_vi_obs
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER :: i,j,k
  INTEGER :: intelm

  intelm = NINT(elm)
  SELECT CASE (intelm)
  CASE(id_cn_obs)  ! concentration of sea-ice         !(SIS)
    !STEVE: need a conversion from modeled to observed field, unless the diagnostic field has it...
!   CALL itpl_3d(v3d(:,:,:,iv3d_hs),ri,rj,rk,yobs)    !(SIS)
    CALL itpl_2d(v2d(:,:,iv2d_cn),ri,rj,yobs)         !(SIS)
  CASE(id_ui_obs)  ! zonal drift of sea-ice           !(SIS)
    CALL itpl_2d(v2d(:,:,iv2d_ui),ri,rj,yobs)         !(SIS)
  CASE(id_vi_obs)  ! meridional drift of sea-ice      !(SIS)
    CALL itpl_2d(v2d(:,:,iv2d_vi),ri,rj,yobs)         !(SIS)
  CASE DEFAULT
    print *, "ERROR::Trans_XtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_hs_obs, id_hi_obs, id_t1_obs, id_t1_obs, id_cn_obs, id_ui_obs, id_vi_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

END SUBROUTINE Trans_XtoY


!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ijk(elem,rlon,rlat,rlev,ri,rj,rk)     !(OCEAN)
  USE vars_model,   ONLY: lon, lat, lev
  USE vars_model,   ONLY: lon0, lonf, lat0, latf
  USE vars_model,   ONLY: lon2d, lat2d !, lev2d
  USE params_model, ONLY: nlon, nlat, nlev
! USE params_model, ONLY: iv2d_eta, iv2d_sst, iv2d_sss, iv2d_ssh
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: aj,ak,ai, rrlon, glon,glat,xlon,xlat
  INTEGER :: i,j,k, n, ni,nj, ii,jj
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

  ! No vertical dimension expected
  rk = 0.0  

END SUBROUTINE phys2ijk


!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
PURE SUBROUTINE itpl_2d(var,ri,rj,var5)
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

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj

!   if (.true.) then
!     print *, "common_obs_sis.f90::itpl_2d::"
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

  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

END SUBROUTINE itpl_2d

PURE SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  USE params_model, ONLY: nlon,nlat,nlev
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

! if (dodebug) WRITE(6,*) "i,j,k,ai,aj,ak = ", i,j,k,ai,aj,ak

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  ELSE
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  END IF

! if (dodebug) WRITE(6,*) "var5 = ", var5

END SUBROUTINE itpl_3d

!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc)
  USE params_obs,   ONLY: id_hs_obs, id_hi_obs, id_t1_obs, id_t2_obs, id_cn_obs, id_ui_obs, id_vi_obs
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  REAL(r_size) :: rmse_hs,rmse_hi,rmse_t1,rmse_t2,rmse_cn,rmse_ui,rmse_vi
  REAL(r_size) :: bias_hs,bias_hi,bias_t1,bias_t2,bias_cn,bias_ui,bias_vi
  INTEGER :: n,ihs,ihi,it1,it2,icn,iui,ivi

  rmse_hs = 0.0d0
  rmse_hi = 0.0d0
  rmse_t1 = 0.0d0
  rmse_t2 = 0.0d0    !(OCEAN)
  rmse_cn = 0.0d0  !(OCEAN)
  rmse_ui = 0.0d0  !(OCEAN)
  rmse_vi = 0.0d0  !(OCEAN)
  bias_hs = 0.0d0
  bias_hi = 0.0d0
  bias_t1 = 0.0d0
  bias_t2 = 0.0d0    !(OCEAN)
  bias_cn = 0.0d0  !(OCEAN)
  bias_ui = 0.0d0  !(OCEAN)
  bias_vi = 0.0d0  !(OCEAN)
  ihs = 0
  ihi = 0
  it1 = 0
  it2 = 0            !(OCEAN)
  icn = 0          !(OCEAN)
  iui = 0          !(OCEAN)
  ivi = 0          !(OCEAN)
  DO n=1,nn
    IF(qc(n) /= 1) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_hs_obs)
      rmse_hs = rmse_hs + dep(n)**2
      bias_hs = bias_hs + dep(n)
      ihs = ihs + 1
    CASE(id_hi_obs)
      rmse_hi = rmse_hi + dep(n)**2
      bias_hi = bias_hi + dep(n)
      ihi = ihi + 1
    CASE(id_t1_obs)
      rmse_t1 = rmse_t1 + dep(n)**2
      bias_t1 = bias_t1 + dep(n)
      it1 = it1 + 1
    CASE(id_t2_obs)                    !(OCEAN)
      rmse_t2 = rmse_t2 + dep(n)**2     !(OCEAN)
      bias_t2 = bias_t2 + dep(n)        !(OCEAN)
      it2 = it2 + 1                     !(OCEAN)
    CASE(id_cn_obs)                  !(OCEAN)
      rmse_cn = rmse_cn + dep(n)**2 !(OCEAN)
      bias_cn = bias_cn + dep(n)    !(OCEAN)
      icn = icn + 1                 !(OCEAN)
    CASE(id_ui_obs)                  !(OCEAN)
      rmse_ui = rmse_ui + dep(n)**2 !(OCEAN)
      bias_ui = bias_ui + dep(n)    !(OCEAN)
      iui = iui + 1                 !(OCEAN)
    CASE(id_vi_obs)                  !(OCEAN)
      rmse_vi = rmse_vi + dep(n)**2 !(OCEAN)
      bias_vi = bias_vi + dep(n)    !(OCEAN)
      ivi = ivi + 1                 !(OCEAN)
    END SELECT
  END DO
  IF(ihs == 0) THEN
    rmse_hs = undef
    bias_hs = undef
  ELSE
    rmse_hs = SQRT(rmse_hs / REAL(ihs,r_size))
    bias_hs = bias_hs / REAL(ihs,r_size)
  END IF
  IF(ihi == 0) THEN
    rmse_hi = undef
    bias_hi = undef
  ELSE
    rmse_hi = SQRT(rmse_hi / REAL(ihi,r_size))
    bias_hi = bias_hi / REAL(ihi,r_size)
  END IF
  IF(it1 == 0) THEN
    rmse_t1 = undef
    bias_t1 = undef
  ELSE
    rmse_t1 = SQRT(rmse_t1 / REAL(it1,r_size))
    bias_t1 = bias_t1 / REAL(it1,r_size)
  END IF
  IF(it2 == 0) THEN                                   !(OCEAN)
    rmse_t2 = undef                                   !(OCEAN)
    bias_t2 = undef                                   !(OCEAN)
  ELSE
    rmse_t2 = SQRT(rmse_t2 / REAL(it2,r_size))          !(OCEAN)
    bias_t2 = bias_t2 / REAL(it2,r_size)                !(OCEAN)
  END IF
  IF(icn == 0) THEN                                 !(OCEAN)
    rmse_cn = undef                                 !(OCEAN)
    bias_cn = undef                                 !(OCEAN)
  ELSE
    rmse_cn = SQRT(rmse_cn / REAL(icn,r_size))    !(OCEAN)
    bias_cn = bias_cn / REAL(icn,r_size)          !(OCEAN)
  END IF
  IF(iui == 0) THEN                                 !(OCEAN)
    rmse_ui = undef                                 !(OCEAN)
    bias_ui = undef                                 !(OCEAN)
  ELSE
    rmse_ui = SQRT(rmse_ui / REAL(iui,r_size))    !(OCEAN)
    bias_ui = bias_ui / REAL(iui,r_size)          !(OCEAN)
  END IF
  IF(ivi == 0) THEN                                 !(OCEAN)
    rmse_vi = undef                                 !(OCEAN)
    bias_vi = undef                                 !(OCEAN)
  ELSE
    rmse_vi = SQRT(rmse_vi / REAL(ivi,r_size))    !(OCEAN)
    bias_vi = bias_vi / REAL(ivi,r_size)          !(OCEAN)
  END IF

  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ============================================='
  WRITE(6,'(7A12)') 'HS','HI','T1','T2','CN','UI','VI'
  WRITE(6,'(7ES12.3)') bias_hs,bias_hi,bias_t1,bias_t2,bias_cn,bias_ui,bias_vi
  WRITE(6,'(7ES12.3)') rmse_hs,rmse_hi,rmse_t1,rmse_t2,rmse_cn,bias_ui,rmse_vi
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ============================'
  WRITE(6,'(7A12)') 'HS','HI','T1','T2','CN','UI','VI'
  WRITE(6,'(7I12)') ihs,ihi,it1,it2,icn,iui,ivi
  WRITE(6,'(A)') '========================================================================'

END SUBROUTINE monit_dep
!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nrec,nn)
  USE params_obs,   ONLY: id_hs_obs, id_hi_obs, id_t1_obs, id_t2_obs, id_cn_obs, id_ui_obs, id_vi_obs
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:) 
  INTEGER :: ios
  INTEGER :: ihs,ihi,it1,it2,icn,iui,ivi
  INTEGER :: iunit
  LOGICAL :: ex
  LOGICAL, PARAMETER :: dodebug=.false.
  INTEGER :: cnt_ice, cnt_water, cnt_mix

  ALLOCATE(wk(nrec))
  nn = 0
  ihs = 0
  ihi = 0
  it1 = 0
  it2 = 0       !(OCEAN)
  icn = 0       !(OCEAN)
  iui = 0       !(OCEAN)
  ivi = 0       !(OCEAN)
  cnt_ice = 0
  cnt_water = 0
  cnt_mix = 0
  iunit=91
  if (dodebug) print *, "get_nobs::"
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
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
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_hs_obs)
        ihs = ihs + 1
      CASE(id_hi_obs)
        ihi = ihi + 1
      CASE(id_t1_obs)
        it1 = it1 + 1
      CASE(id_t2_obs)     !(OCEAN)
        it2 = it2 + 1      !(OCEAN)
      CASE(id_cn_obs)   !(OCEAN)
        icn = icn + 1  !(OCEAN)
        if (wk(5) > 0) then
          if (wk(5) < 1) then
            cnt_mix = cnt_mix + 1
          else
            cnt_ice = cnt_ice + 1
          endif
        else
          cnt_water = cnt_water + 1
        endif
      CASE(id_ui_obs)   !(OCEAN)
        iui = iui + 1  !(OCEAN)
      CASE(id_vi_obs)   !(OCEAN)
        ivi = ivi + 1  !(OCEAN)
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,'(I10,A)') nn,' TOTAL OBSERVATIONS INPUT (in get_nobs)'
    WRITE(6,'(A12,I10)') '         HS:',ihs
    WRITE(6,'(A12,I10)') '         HI:',ihi
    WRITE(6,'(A12,I10)') '         T1:',it1
    WRITE(6,'(A12,I10)') '         T2:',it2   !(OCEAN)
    WRITE(6,'(A12,I10)') '         CN:',icn !(OCEAN)
    WRITE(6,'(A12,I10)') '         UI:',iui !(OCEAN)
    WRITE(6,'(A12,I10)') '         VI:',ivi !(OCEAN)
    WRITE(6,'(A27)') 'Concentration Distribution:'
    WRITE(6,'(A13,I10,A11,F5.2,A1)') ' water/land: ',cnt_water, ' cells, or ', 100*REAL(cnt_water)/REAL(icn),'%'
    WRITE(6,'(A13,I10,A11,F5.2,A1)') '  mixed ice: ',cnt_mix,   ' cells, or ', 100*REAL(cnt_mix)/REAL(icn),'%'
    WRITE(6,'(A13,I10,A11,F5.2,A1)') '  solid ice: ',cnt_ice,   ' cells, or ', 100*REAL(cnt_ice)/REAL(icn),'%'
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  if (nn .eq. 0) then
    WRITE(6,*) "get_nobs:: WARNING: No observations have been found! Ensure that at least one timeslot has observations."
!   STOP(60)
  endif

END SUBROUTINE get_nobs


SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,obhr)
  USE params_model, ONLY: nlon
  USE vars_model, ONLY: lon0, lonf
  USE vars_model, ONLY: lon,lat
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
      print *, "common_obs_sis.f90::read_obs:: WARNING!"
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

  END DO
  CLOSE(iunit)

  if (MAXVAL(rlon) > lonf) then
    WRITE(6,*) "read_obs:: Error: MAX(observation lon, i.e. rlon) > lonf"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lonf = ", lonf
    WRITE(6,*) "lon(nlon) = ", lon(nlon)
    STOP (2)
  endif
  if (MINVAL(rlon) < lon0) then
    WRITE(6,*) "read_obs:: Error: MIN(observation lon, i.e. rlon) < lon0"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lon0 = ", lon0
    WRITE(6,*) "lon(1) = ", lon(1)
    STOP (2)
  endif

END SUBROUTINE read_obs


SUBROUTINE read_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
  USE params_model, ONLY: obs_noise_coeff
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
  REAL(r_size) :: rand(nn)   !(SIS)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit

  !STEVE: trying to handle numerical issues with ice concentration data type: (SIS)
  CALL com_randn(nn,rand) !(SIS)

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    READ(iunit) wk
!   SELECT CASE(NINT(wk(1)))
!   CASE(id_u_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_v_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_t_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!   CASE(id_q_obs)
!     wk(4) = wk(4) * 100.0 ! hPa -> Pa
!  END SELECT
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    ohx(n)  = REAL(wk(7),r_size)
    oqc(n)  = NINT(wk(8))
    obhr(n) = REAL(wk(9),r_size)

    !STEVE: add some 'inflation' (noise) to the ohx so the ensemble spread is not 0, (SIS)
    !       which is particularly important in high ice or open water areas...
    if ( odat(n) >= 0.99 .or. odat(n) <= 0.01 ) then
      ohx(n) = ohx(n) + rand(n) * (obs_noise_coeff * oerr(n)) !STEVE: scale the added noise to a percentage of the obs error
    endif
  END DO
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
  DO n=1,nn
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
  END DO
  CLOSE(iunit)

END SUBROUTINE write_obs

SUBROUTINE write_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
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
  INTEGER,INTENT(IN) :: oqc(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit
  LOGICAL, PARAMETER :: dodebug=.false.

  iunit=92
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    wk(1) = REAL(elem(n),r_sngl)  ! ID for observation type
    wk(2) = REAL(rlon(n),r_sngl)  ! Ob lon
    wk(3) = REAL(rlat(n),r_sngl)  ! Ob lat
    wk(4) = REAL(rlev(n),r_sngl)  ! Ob level
    wk(5) = REAL(odat(n),r_sngl)  ! Observed data quantity
    wk(6) = REAL(oerr(n),r_sngl)  ! Estimated observation error
    wk(7) = REAL(ohx(n),r_sngl)   ! Model forecast transformed to observation space: H(xb)
    wk(8) = REAL(oqc(n),r_sngl)   ! Quality control ID (1==keep, 0==discard) for use in assimilation
    wk(9) = REAL(obhr(n),r_sngl)   ! Quality control ID (1==keep, 0==discard) for use in assimilation
    if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
    WRITE(iunit) wk
  END DO
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
