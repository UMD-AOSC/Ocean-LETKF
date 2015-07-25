PROGRAM gsw_pot_to_insitu

USE common
USE params_model

IMPLICIT NONE

REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)

! Loop through all temperature observations
! already interpolated to the obs location
! and convert from potential temperature
! to the in situ temperature
do n=1,nobs
  if ( obselm(i) .ne. t_obs_id) CYCLE
  
  ! ri,rj,rk :=  phys2ijk(olon,olat,olev)
  CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk)

  ! pt := TransXtoY v3d(:,:,:,iv3d_t) 
  ! sp := TransXtoY v3d(:,:,:,iv3d_s) 
  CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,pt)
  CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,sp)

  !STEVE: eventually use pressure from model output,
  !       for now, use simple computation:
  !  p := g*h*rho
  !       Use average g and average rho
  !       For h, use observed depth from argo data
! rho = 1025 ! 1027[kg/m^3]
  !STEVE: slightly more complex version:
  p = p_from_z(obslev(n),obslat(n))

  !STEVE: for some reason, realtime GTS data is
  !       distributed using depth as the vertical
  !       coordinate rather than pressure:
  !       http://www.argo.ucsd.edu/Argo_date_guide.html#gtsusers

  ohx(n) = t_from_pt(pt,sp,p,obslon(n),obslat(n))

enddo

CONTAINS

FUNCTION t_from_pt(pt,sp,p,lon,lat)

IMPLICIT NONE

INTEGER, PARAMETER :: r14 = selected_real_kind(14,30)
REAL(r14), INTENT(IN) :: pt, sp, p, lon, lat
REAL(r14) :: t_from_pt
REAL(r14) :: sa, ct, T

! Compute in situ temperature from potential temperature by computing
! the conservative temperature from potential temperature then computing
! the in situ temperature from the conservative temperature.
! sa     : Absolute Salinity                   [g/kg]
! sp     : Practical Salinity                  [unitless]
! ct     : Conservative Temperature            [deg C]
! pt     : potential temperature with          [deg C]
!           reference pressure of 0 dbar
! c      : conductivity                        [mS/cm]
! t      : in-situ temperature [ITS-90]        [deg C]
! p      : sea pressure                        [dbar]
! z      : depth                               [m]

! Compute the absolute salinity from the practical salinity
sa = gsw_sa_from_sp(sp,p,long,lat)

! Compute the conservative temperature
ct = gsw_ct_from_pt(sa,pt)

!! Compute the pressure
!gsw_rho(sa,ct,p)

! Compute the in situ temperature
t_from_pt  = gsw_t_from_ct(sa,ct,p)

END FUNCTION t_from_pt

FUNCTION p_from_z(dpth,xlat)
! pressure from depth from saunder's formula with eos80.
! reference: saunders,peter m., practical conversion of pressure
!            to depth., j.p.o. , april 1981.
! r millard
! march 9, 1983
! check value: p_from_z=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
! http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/ppsw.f
  IMPLICIT NONE
  REAL(r14), INTENT(IN) :: dpth
  REAL(r14), INTENT(IN) :: xlat
  PARAMETER :: pi=3.141592654
  REAL(r14) :: p_from_z
  REAL(r14) :: plat, d, c1

  plat=abs(xlat*pi/180.)
  d=sin(plat)
  c1=5.92e-3+d**2*5.25e-3
  p_from_z=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6

END FUNCTION p_from_z

END PROGRAM
