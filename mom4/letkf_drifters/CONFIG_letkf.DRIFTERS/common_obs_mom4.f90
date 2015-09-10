MODULE common_obs_mom4
!===============================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   04/26/2011 Steve PENNY converted to OCEAN for use with MOM4. grep '(OCEAN)' for changes.
!
!===============================================================================
  USE common
  USE common_mom4
  USE params_obs

  IMPLICIT NONE
  PUBLIC

CONTAINS

SUBROUTINE Trans_XtoY(elm,ri,rj,rk,v3d,v2d,yobs)        !(OCEAN)
!===============================================================================
! Transformation from model space to observation space (i.e. H-operator)
!===============================================================================
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
  CASE(id_u_obs)  ! U
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ri,rj,rk,yobs)
  CASE(id_v_obs)  ! V
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,rj,rk,yobs)
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,yobs)
  CASE(id_s_obs)  ! S                             !(OCEAN)
    CALL itpl_3d(v3d(:,:,:,iv3d_s),ri,rj,rk,yobs) !(OCEAN)
  CASE(id_eta_obs) ! SSH                          !(OCEAN)
    !STEVE: use mom4 surface height to form Yb (i.e. hdxf)
    CALL itpl_2d(v2d(:,:,iv2d_eta),ri,rj,yobs)    !(OCEAN)
  CASE(id_sst_obs) ! SST                          !(OCEAN)
    CALL itpl_2d(v2d(:,:,iv2d_sst),ri,rj,yobs)    !(OCEAN)
  CASE(id_sss_obs) ! SSS                          !(OCEAN)
    CALL itpl_2d(v2d(:,:,iv2d_sss),ri,rj,yobs)    !(OCEAN)
  CASE DEFAULT
    print *, "ERROR::Trans_XtoY:: observation type not recognized."
    print *, "element id = ", intelm
    print *, "available id's = ", id_u_obs, id_v_obs, id_t_obs, id_s_obs, id_ssh_obs, id_eta_obs, id_sst_obs, id_sss_obs, id_x_obs, id_y_obs, id_z_obs
    print *, "STEVE: STOPPING ON PURPOSE..."
    STOP 1
  END SELECT

  RETURN
END SUBROUTINE Trans_XtoY


SUBROUTINE phys2ijk(elem,rlon,rlat,rlev,ri,rj,rk)     !(OCEAN)
!===============================================================================
! Coordinate conversion
!===============================================================================
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

!STEVE: probably want to replace this with a kd-tree lookup for general search (KDTREE)

!
! rlon -> ri
!
  rrlon=rlon
  ! STEVE: Since this is mapping physical coordinates to the model grid, this is
  ! my preferred/ideal place to do the conversion. However, to prevent potential unknown
  ! issues within the code, I converted the observation coordinates immediately
  ! upon reading the data.
  !
  ! If the observations and map are on a different coordinate grid, e.g. NCEP
  ! mom4p1 on -285 to 75, and obs on 0 to 360, then adjust obs to map
  ! coordinates (attempting to do it in a general way)
! IF(rlon >= lonf) THEN
!   rrlon = REAL(lon0 + abs(rlon - lonf) - wrapgap,r_size)
! ELSE IF(rlon < lon0) THEN
!   rrlon = REAL(lonf - abs(lon0 - rlon) - wrapgap,r_size)
! ENDIF


  !STEVE: initialize to something that will throw errors if it's not changed within
  ri = -1
  rj = -1
  rk = -1

  ! Find the correct longitude:
  DO i=1,nlon
    IF(rrlon < lon(i)) EXIT
  END DO

  ! Set flags to mark rlon outside of map range:
  IF(i == 1) THEN
    ri = 0
  ELSE IF(i == nlon+1) THEN
    ri = nlon+1
  ELSE
    ! Otherwise, Interpolate
    ai = (rrlon - lon(i-1)) / (lon(i) - lon(i-1))
    ri = REAL(i-1,r_size) + ai
  END IF

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

  IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) RETURN

!
! rlat -> rj
!
  DO j=1,nlat
    IF(rlat < lat(j)) EXIT
  END DO

  IF(j == 1) THEN
    rj = (rlat + 90.0d0) / (lat(1) + 90.0d0)
  ELSE IF(j == nlat+1) THEN
    aj = (rlat - lat(nlat)) / (90.0d0 - lat(nlat))
    rj = REAL(nlat,r_size) + aj
  ELSE
    !STEVE: Interpolate
    aj = (rlat - lat(j-1)) / (lat(j) - lat(j-1))
    rj = REAL(j-1,r_size) + aj
  END IF
  IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) RETURN

!
! rlev -> rk
!
  IF(NINT(elem) == id_ssh_obs) THEN     ! surface observation !(OCEAN)
    rk = 0.0d0
  ELSEIF(NINT(elem) == id_eta_obs) THEN ! surface observation !(OCEAN)
    rk = 0.0d0
  ELSEIF(NINT(elem) == id_sst_obs) THEN ! surface observation !(OCEAN)
    rk = 0.0d0
  ELSEIF(NINT(elem) == id_sss_obs) THEN ! surface observation !(OCEAN)
    rk = 0.0d0
  ELSE
    !
    ! vertical interpolation
    !
    !
    ! find rk
    !
    DO k=1,nlev
      IF(rlev < lev(k)) EXIT
      IF(k .eq. nlev .and. rlev .eq. lev(nlev)) EXIT !STEVE: added this case for simulated obs that reach the lowest model levels.
                                                     !       Otherwise, k iterates to nlev+1 before exiting loop.
    END DO

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

  END IF

END SUBROUTINE phys2ijk

PURE SUBROUTINE itpl_2d(var,ri,rj,var5)
!===============================================================================
! Bilinear nterpolation in 2D
!===============================================================================
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
!     print *, "common_obs_mom4.f90::itpl_2d::"
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

  RETURN
END SUBROUTINE itpl_2d

PURE SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
!===============================================================================
! Interpolation in 3D
!===============================================================================
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

  RETURN
END SUBROUTINE itpl_3d


SUBROUTINE monit_dep(nn,elm,dep,qc)
!===============================================================================
! Monitor departure
!===============================================================================
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
  IF(is == 0) THEN                                   !(OCEAN)
    rmse_s = undef                                   !(OCEAN)
    bias_s = undef                                   !(OCEAN)
  ELSE
    rmse_s = SQRT(rmse_s / REAL(is,r_size))          !(OCEAN)
    bias_s = bias_s / REAL(is,r_size)                !(OCEAN)
  END IF
  IF(issh == 0) THEN                                 !(OCEAN)
    rmse_ssh = undef                                 !(OCEAN)
    bias_ssh = undef                                 !(OCEAN)
  ELSE
    rmse_ssh = SQRT(rmse_ssh / REAL(issh,r_size))    !(OCEAN)
    bias_ssh = bias_ssh / REAL(issh,r_size)          !(OCEAN)
  END IF
  IF(ieta == 0) THEN                                 !(OCEAN)
    rmse_eta = undef                                 !(OCEAN)
    bias_eta = undef                                 !(OCEAN)
  ELSE
    rmse_eta = SQRT(rmse_eta / REAL(ieta,r_size))    !(OCEAN)
    bias_eta = bias_eta / REAL(ieta,r_size)          !(OCEAN)
  END IF
  IF(isst == 0) THEN                                 !(OCEAN)
    rmse_sst = undef                                 !(OCEAN)
    bias_sst = undef                                 !(OCEAN)
  ELSE
    rmse_sst = SQRT(rmse_sst / REAL(isst,r_size))    !(OCEAN)
    bias_sst = bias_sst / REAL(isst,r_size)          !(OCEAN)
  END IF
  IF(isss == 0) THEN                                 !(OCEAN)
    rmse_sss = undef                                 !(OCEAN)
    bias_sss = undef                                 !(OCEAN)
  ELSE
    rmse_sss = SQRT(rmse_sss / REAL(isss,r_size))    !(OCEAN)
    bias_sss = bias_sss / REAL(isss,r_size)          !(OCEAN)
  END IF

  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ============================================='
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','eta','SST','SSS'                                       !(OCEAN)
  WRITE(6,'(7ES12.3)') bias_u,bias_v,bias_t,bias_s,bias_ssh,bias_eta,bias_sst,bias_sss               !(OCEAN)
  WRITE(6,'(7ES12.3)') rmse_u,rmse_v,rmse_t,rmse_s,rmse_ssh,bias_eta,rmse_sst,bias_sss               !(OCEAN)
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ============================'
  WRITE(6,'(7A12)') 'U','V','T','S','SSH','eta','SST','SSS'                                       !(OCEAN)
  WRITE(6,'(7I12)') iu,iv,it,is,issh,ieta,isst,isss                                              !(OCEAN)
  WRITE(6,'(A)') '========================================================================'

  RETURN
END SUBROUTINE monit_dep

SUBROUTINE get_nobs(cfile,nrec,nn)
!===============================================================================
! Read in observation file to count the number of observations
!===============================================================================
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl),ALLOCATABLE :: wk(:) 
  INTEGER :: ios
  INTEGER :: iu,iv,it,is,issh,ieta,isst,isss,ix,iy,iz !(OCEAN)
  INTEGER :: iunit
  LOGICAL :: ex
  LOGICAL :: dodebug=.false.

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
      nn = nn + 1
    END DO
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT (in get_nobs)'
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
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF
  DEALLOCATE(wk)

  RETURN
END SUBROUTINE get_nobs


SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,obhr)
!===============================================================================
! Read in observations
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
  REAL(r_size), OPTIONAL, INTENT(OUT) :: obhr(nn)
  REAL(r_sngl), ALLOCATABLE :: wk(:) !obs1nrec
  !REAL(r_size) :: wk(6) !(OCEAN) STEVE: I changed this because the netcdf observation files are stored as double precision
  INTEGER :: n,iunit
  ! STEVE: for general grid
  LOGICAL :: dodebug = .false.
  LOGICAL :: process_obs = .true.

  ALLOCATE(wk(obs1nrec))
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
    if (obs1nrec > 6) then
      obhr(n) = REAL(wk(7),r_size)
    endif
    !STEVE: error check
    if (oerr(n) .le. 0) then
      print *, "common_obs_mom4.f90::read_obs:: WARNING!"
      print *, "STEVE: oerr <= 0, must be > 0 ..." 
      print *, "STEVE: oerr(",n,") = ", oerr(n)
      PRINT '(I6,2F7.2,F10.2,2F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6)
      stop 9
    endif
    
    ! Special processing for obs:
    !STEVE: this handles the fact that the observations are typically on a 0 to 360º grid, while
    !       the NCEP mom4p1 grid configuration is on a ~ -285 to 75º grid
    if (process_obs) then
      IF(rlon(n) >= lonf) THEN
        if (dodebug) then
          WRITE(6,*) "letkf_obs.f90:: Wrapping large lon obs to model grid: n = ", n
          WRITE(6,*) "pre  - rlon(n) = ", rlon(n)
          WRITE(6,*) "360 - rlon(n) = ", 360.0 - rlon(n)
        endif
        ! Update the coordinate
        if (abs(rlon(n) - lonf) < wrapgap ) then
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
          ! Otherwise, wrap the observation coordinate to be inside of the defined model grid coordinates
          !Wrap the coordinate
          rlon(n) = REAL(lon0 + abs(rlon(n) - lonf) - wrapgap,r_size)
        endif
        if (dodebug) WRITE(6,*) "post - rlon(n) = ", rlon(n)
        ELSE IF(rlon(n) < lon0) THEN
        if (dodebug) then
          WRITE(6,*) "letkf_obs.f90:: Wrapping small lon obs to model grid: n = ", n
          WRITE(6,*) "pre  - rlon(n) = ", rlon(n)
          WRITE(6,*) "360 - rlon(n) = ", 360.0 - rlon(n)
        endif
        ! Update the coordinate
        if (abs(lon0 - rlon(n)) < wrapgap ) then
          !STEVE: shift it if it's just outside grid
          if (abs(lon0 - rlon(n)) < wrapgap/2) then
            rlon(n) = lon0
          else
            rlon(n) = lonf 
          endif
          ! Increase error to compensate                                 
          oerr(n) = oerr(n)*2
        else
          !Wrap the coordinate
          rlon(n) = REAL(lonf - abs(lon0 - rlon(n)) + wrapgap,r_size)
        endif
        if (dodebug) WRITE(6,*) "post - rlon(n) = ", rlon(n)
      ENDIF
    endif 
  END DO
  CLOSE(iunit)

  if (MAXVAL(rlon) > lonf) then
    WRITE(6,*) "read_obs:: Error: MAX(observation lon, i.e. rlon) > lonf"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lonf = ", lonf
    WRITE(6,*) "lon(nlon) = ", lon(nlon)
    STOP(2)
  endif
  if (MINVAL(rlon) < lon0) then
    WRITE(6,*) "read_obs:: Error: MIN(observation lon, i.e. rlon) < lon0"
    WRITE(6,*) "rlon = ", rlon
    WRITE(6,*) "lon0 = ", lon0
    WRITE(6,*) "lon(1) = ", lon(1)
    STOP(2)
  endif

  RETURN
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


!STEVE: adding for support
SUBROUTINE write_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,obhr)
!===============================================================================
! Write out observations
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
  REAL(r_size),OPTIONAL,INTENT(IN) :: obhr(nn)
  REAL(r_sngl),ALLOCATABLE :: wk(:)
! REAL(r_size) :: wk(6) !(OCEAN) STEVE: I changed this because the netcdf observation files are stored as double precision
  INTEGER :: n,iunit

  ALLOCATE(wk(obs1nrec))
  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  do n=1,nn
    wk(1) = REAL(elem(n),r_sngl)
    wk(2) = REAL(rlon(n),r_sngl)
    wk(3) = REAL(rlat(n),r_sngl)
    wk(4) = REAL(rlev(n),r_sngl)
    wk(5) = REAL(odat(n),r_sngl)
    wk(6) = REAL(oerr(n),r_sngl)
    if (obs1nrec>6) then
      wk(7) = REAL(obhr(n),r_sngl)
    endif
    WRITE(iunit) wk
  enddo

  CLOSE(iunit)

END SUBROUTINE write_obs


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
    wk(1) = REAL(elem(n),r_sngl)  ! ID for observation type
    wk(2) = REAL(rlon(n),r_sngl)  ! Ob lon
    wk(3) = REAL(rlat(n),r_sngl)  ! Ob lat
    wk(4) = REAL(rlev(n),r_sngl)  ! Ob level
    wk(5) = REAL(odat(n),r_sngl)  ! Observed data quantity
    wk(6) = REAL(oerr(n),r_sngl)  ! Estimated observation error
    wk(7) = REAL(ohx(n),r_sngl)   ! Model forecast transformed to observation space: H(xb)
    wk(8) = REAL(oqc(n),r_sngl)   ! Quality control ID (1==keep, 0==discard) for use in assimilation
    if (obs2nrec>8) then
      wk(9) = REAL(obhr(n),r_sngl)   ! Quality control ID (1==keep, 0==discard) for use in assimilation
      if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8),wk(9)
    else
      if (dodebug) PRINT '(I6,2F7.2,F10.2,6F12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7),wk(8)
    endif
    WRITE(iunit) wk
  enddo

  CLOSE(iunit)

END SUBROUTINE write_obs2

END MODULE common_obs_mom4
