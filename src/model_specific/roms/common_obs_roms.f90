MODULE common_obs_oceanmodel
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   02/03/2009 Takemasa MIYOSHI  modified for ROMS
!
!=======================================================================
  USE common
  USE params_model
  USE vars_model,     ONLY: phi0
  USE params_obs

  IMPLICIT NONE
  PUBLIC

CONTAINS
!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rlev,v3d,v2d,yobs)
  USE common_oceanmodel, ONLY: calc_depth
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rlev
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  REAL(r_size)             :: wk1(1),wk2(1),depth(nlev)
  INTEGER :: k

  SELECT CASE (NINT(elm))
  CASE(id_u_obs) ! U
    yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_u) ! only surface
  CASE(id_v_obs) ! V
    yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_v) ! only surface
  CASE(id_t_obs) ! T
    wk1(1) = rlev
    CALL calc_depth(v2d(NINT(ri),NINT(rj),iv2d_z),phi0(NINT(ri),NINT(rj)),depth) 
    CALL com_interp_spline(nlev,depth,v3d(NINT(ri),NINT(rj),:,iv3d_t),1,wk1,wk2)
    yobs = wk2(1) 
!    write(*,*),'wk1(1)_t>>',wk1(1),'wk2(1)_t>>',wk2
  CASE(id_s_obs) ! S
    wk1(1) = rlev
    CALL calc_depth(v2d(NINT(ri),NINT(rj),iv2d_z),phi0(NINT(ri),NINT(rj)),depth) 
    CALL com_interp_spline(nlev,depth,v3d(NINT(ri),NINT(rj),:,iv3d_s),1,wk1,wk2)
    yobs = wk2(1)
!    write(*,*),'wk1(1)_s>>',wk1(1),'wk2(1)_s>>',wk2
!    write(*,*),'next>>'
  CASE(id_z_obs) ! Z
    yobs = v2d(NINT(ri),NINT(rj),iv2d_z)
  END SELECT

END SUBROUTINE Trans_XtoY


!(UPDATE) with lon2d/lat2d search for observation (LON2D/LAT2D)
SUBROUTINE phys2ijk(elem,rlon,rlat,rlev,ri,rj,rk)     !(OCEAN)
!===============================================================================
! Coordinate conversion
!===============================================================================
  USE vars_model,   ONLY: lon, lat, lev, lon0, lonf, lon2d, lat2d, lev2d
  USE params_model, ONLY: nlon, nlat, nlev
  USE params_model, ONLY: iv2d_eta, iv2d_sst, iv2d_sss, iv2d_ssh
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: aj,ak,ai, rrlon
! REAL(r_size) :: lnps(nlon,nlat)
  REAL(r_size) :: plev(nlev)
  INTEGER :: i,j,k
  LOGICAL :: dodebug = .false.
  ! To speed up 2d lon/lat search:
  INTEGER, PARAMETER :: ird = 10 ! gridpoint radius
  INTEGER, PARAMETER :: jrd = 10 ! gridpoint radius
  REAL(r_size) :: tripolar_lat = 65.0 ! cutoff for tripolar grid
  LOGICAL :: ijset=.false.
  REAL(r_size) :: lonA,lonB,latA,latB
  INTEGER :: i0,j0,i1,j1
!STEVE: probably want to replace this with a kd-tree lookup for general search (KDTREE)

  !STEVE: initialize to something that will throw errors if it's not changed within
  ri = -1
  rj = -1
  rk = -1
  rrlon=rlon

  if (rlon > tripolar_lat) then

    !---------------------------------------------------------------------------
    ! Get a guess of the starting longitude
    !---------------------------------------------------------------------------
    do i=1,nlon
      if (rrlon < lon(i)) EXIT
    enddo
    i0 = max(1,   i-ird)
    i1 = min(nlon,i+ird)

    !---------------------------------------------------------------------------
    ! Get a guess of the starting latitude
    !---------------------------------------------------------------------------
    do j=1,nlat
      if (rlat < lat(j)) EXIT
    enddo
    j0 = max(1,  j-jrd)
    j1 = min(nlat,j+jrd)

    lonA = lon(i-1)
    lonB = lon(i)
    latA = lat(j-1)
    latB = lat(j)

  endif

  if (i < 2 .OR. nlon < i) RETURN
  if (j < 2 .OR. nlat < j) RETURN

  ! Interpolate lons
  ai = (rrlon - lonA) / (lonB - lonA)
  ri = REAL(i-1,r_size) + ai

  ! Interpolate lats
  aj = (rlat - latA) / (latB - latA)
  rj = REAL(j-1,r_size) + aj

  if (dodebug .and. (ri > nlon .or. ri < 1)) then
    WRITE(6,*) "In common_obs_mom4.f90::phys2ijk,"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "rrlon = ", rrlon
    WRITE(6,*) "ai = ", ai
    WRITE(6,*) "ri = ", ri
    WRITE(6,*) "lon0 = ", lon0
    WRITE(6,*) "lonf = ", lonf
    WRITE(6,*) "i = ", i
  endif

  !-----------------------------------------------------------------------------
  ! rlev -> rk
  !-----------------------------------------------------------------------------
  if (NINT(elem) == id_ssh_obs) then     ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif (NINT(elem) == id_eta_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif (NINT(elem) == id_sst_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  elseif (NINT(elem) == id_sss_obs) then ! surface observation !(OCEAN)
    rk = 0.0d0
  else

    !
    ! do vertical interpolation
    !

    !
    ! find rk
    !
    do k=1,nlev
      if (rlev < lev(k)) EXIT
      if (k .eq. nlev .and. rlev .eq. lev(nlev)) EXIT !STEVE: added this case for simulated obs that reach the lowest model levels.
                                                     !       Otherwise, k iterates to nlev+1 before exiting loop.
    enddo

    if (k .eq. 1) then
!     print *, "k = 1, rlev = ", rlev
!     print *, "STEVE: STOPPING ON PURPOSE in common_obs_mom4.f90..."
!     print *, "We can't have k=1"
!     stop 1
      !STEVE: interpolate from SFC (0) to model level 1
      !print *, "NOTICE: observation is above model SFC level => i,j,k = ",i,j,k 
      !ak = (rlev - 0) / (lev(k) - 0)
      rk = 1 !ak
    else
      !STEVE: now apply the interpolation at the identified model level:
      ak = (rlev - lev(k-1)) / (lev(k) - lev(k-1))
      rk = REAL(k-1,r_size) + ak
    endif

  endif

END SUBROUTINE phys2ijk


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

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

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

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

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

END SUBROUTINE itpl_3d

!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nn,len)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  INTEGER,INTENT(IN),OPTIONAL :: len
  REAL(r_sngl) :: wk(6)
  INTEGER :: ios
  INTEGER :: iu,iv,it,is,iz
  INTEGER :: iunit
  LOGICAL :: ex

  nn = 0
  iu = 0
  iv = 0
  it = 0
  is = 0
  iz = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)

  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',ACTION='read')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs)
        iu = iu + 1
      CASE(id_v_obs)
        iv = iv + 1
      CASE(id_t_obs)
        it = it + 1
      CASE(id_s_obs)
        is = is + 1
      CASE(id_z_obs)
        iz = iz + 1
      END SELECT
      nn = nn + 1
    END DO
!    WRITE(*,*) nn,' OBSERVATIONS INPUT - PASSEI'
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '       SALT:',is
    WRITE(6,'(A12,I10)') '       ZETA:',iz
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

END SUBROUTINE get_nobs


SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn) ! longitude
  REAL(r_size),INTENT(OUT) :: rlat(nn) ! latitude
  REAL(r_size),INTENT(OUT) :: rlev(nn) ! depth [meters]
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_sngl) :: wk(6)
  INTEGER :: n,iunit

  iunit=91
  WRITE(6,*) '==  OBSERVATION  =='
  WRITE(6,*) '== Reading cfile ==',cfile

  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  DO n=1,nn
    READ(iunit) wk
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
  END DO
  write(*,*),'elem(1) --- Leo',elem(1)
  write(*,*),'rlon(1) --- Leo',rlon(1)  
  write(*,*),'rlat(1) --- Leo',rlat(1)
  write(*,*),'rlev(1) --- Leo',rlev(1)
  write(*,*),'odat(1) --- Leo',odat(1)
  write(*,*),'oerr(1) --- Leo',oerr(1)

  CLOSE(iunit)

END SUBROUTINE read_obs


SUBROUTINE read_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
!===============================================================================
! Read in observations with appended H(xb) for each ob
!===============================================================================
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
  REAL(r_size),OPTIONAL,INTENT(OUT) :: obhr(nn)
  INTEGER,INTENT(OUT) :: oqc(nn)
  REAL(r_sngl),ALLOCATABLE :: wk(:)
  INTEGER :: n,iunit

  ALLOCATE(wk(obs2nrec))
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
    if (obs2nrec>8) then
      obhr(n) = REAL(wk(9),r_size)
    endif
  enddo

  CLOSE(iunit)

END SUBROUTINE read_obs2


!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,qc)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  INTEGER,INTENT(IN) :: qc(nn)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_s,rmse_z
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_z
  INTEGER :: n,iu,iv,it,is,iz

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_s = 0.0d0
  rmse_z = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_s = 0.0d0
  bias_z = 0.0d0
  iu = 0
  iv = 0
  it = 0
  is = 0
  iz = 0
  DO n=1,nn
    IF(qc(n) /= 1) CYCLE
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
    CASE(id_s_obs)
      rmse_s = rmse_s + dep(n)**2
      bias_s = bias_s + dep(n)
      is = is + 1
    CASE(id_z_obs)
      rmse_z = rmse_z + dep(n)**2
      bias_z = bias_z + dep(n)
      iz = iz + 1
    END SELECT
  END DO
  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(is == 0) THEN
    rmse_s = undef
    bias_s = undef
  ELSE
    rmse_s = SQRT(rmse_s / REAL(is,r_size))
    bias_s = bias_s / REAL(is,r_size)
  END IF
  IF(iz == 0) THEN
    rmse_z = undef
    bias_z = undef
  ELSE
    rmse_z = SQRT(rmse_z / REAL(iz,r_size))
    bias_z = bias_z / REAL(iz,r_size)
  END IF
  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ================================='
  WRITE(6,'(5A12)') 'U','V','T','SALT','ZETA'
  WRITE(6,'(5ES12.3)') bias_u,bias_v,bias_t,bias_s,bias_z
  WRITE(6,'(5ES12.3)') rmse_u,rmse_v,rmse_t,rmse_s,rmse_z
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ================'
  WRITE(6,'(5A12)') 'U','V','T','SALT','ZETA'
  WRITE(6,'(5I12)') iu,iv,it,is,iz
  WRITE(6,'(A)') '============================================================'

  RETURN
END SUBROUTINE monit_dep


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


  WRITE(6,*) "ROMS does not require recentering, returning..."
  RETURN

  !STEVE: This routine can be used to recenter the longitude coorindates of the observations,
  !       in case the model is set up on a grid that does not match the observations.

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


SUBROUTINE write_obs2(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
!===============================================================================
! Write out observations with appended H(xb) for each ob
!===============================================================================
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
  REAL(r_size),OPTIONAL,INTENT(IN) :: obhr(nn)
  INTEGER,INTENT(IN) :: oqc(nn)
  REAL(r_sngl),ALLOCATABLE :: wk(:)
  INTEGER :: n,iunit
  LOGICAL, PARAMETER :: dodebug=.false.

  ALLOCATE(wk(obs2nrec))
  iunit=92
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  do n=1,nn
    wk(1) = REAL(elem(n),r_sngl)   ! ID for observation type
    wk(2) = REAL(rlon(n),r_sngl)   ! Ob lon
    wk(3) = REAL(rlat(n),r_sngl)   ! Ob lat
    wk(4) = REAL(rlev(n),r_sngl)   ! Ob level
    wk(5) = REAL(odat(n),r_sngl)   ! Observed data quantity
    wk(6) = REAL(oerr(n),r_sngl)   ! Estimated observation error
    wk(7) = REAL(ohx(n),r_sngl)    ! Model forecast transformed to observation space: H(xb)
    wk(8) = REAL(oqc(n),r_sngl)    ! Quality control ID (1==keep, 0==discard) for use in assimilation
    if (obs2nrec>8) then
      wk(9) = REAL(obhr(n),r_sngl) ! Time in (noninteger) hours, assuming the data is organized by day
      if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
    else
      if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8)
    endif
    WRITE(iunit) wk
  enddo

  CLOSE(iunit)

END SUBROUTINE write_obs2



END MODULE common_obs_oceanmodel
