MODULE common_obs_roms
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   02/03/2009 Takemasa MIYOSHI  modified for ROMS
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_roms

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332

CONTAINS
!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,ri,rj,rlev,v3d,v2d,yobs)
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

  RETURN
END SUBROUTINE Trans_XtoY
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

  RETURN
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

  RETURN
END SUBROUTINE itpl_3d

!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
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

  RETURN
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

  RETURN
END SUBROUTINE read_obs

END MODULE common_obs_roms
