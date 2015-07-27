MODULE gsw_pot_to_insitu
!===============================================================================
!
! MODULE:
!  gsw_pot_to_insitu
!
! PURPOSE:
!  This module provides additional subroutines needed to compute the in situ
!  temperature from a numerical model's potential temperature (ºC) and practical 
!  salinity (psu). It utilizes the TEOS GSW Fortran package, version 3.03,
!  available at http://www.teos-10.org
!
!
!-------------------------------------------------------------------------------
! Author:  Steve Penny
! Contact: Steve.Penny@noaa.gov
!===============================================================================

IMPLICIT NONE

PUBLIC :: t_from_pt, p_from_z

PRIVATE

INTEGER, PARAMETER :: r14 = selected_real_kind(14,30)

CONTAINS

!===============================================================================
FUNCTION t_from_pt(pt_in,sp_in,p_in,lon_in,lat_in)
!===============================================================================
  USE common, ONLY: r_size
! USE gsw_oceanographic_toolbox, ONLY: gsw_sa_from_sp, gsw_ct_from_pt, gsw_t_from_ct !, r14

  IMPLICIT NONE
  REAL(r14) :: gsw_sa_from_sp, gsw_ct_from_pt, gsw_t_from_ct
  REAL(r_size) :: t_from_pt
  REAL(r_size), INTENT(IN) :: pt_in, sp_in, p_in, lon_in, lat_in
  REAL(r14) :: pt, sp, p, lon, lat
  REAL(r14) :: sa, ct, T
 
  pt = pt_in
  sp = sp_in
  p  = p_in
  lon = lon_in
  lat = lat_in

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
  sa = gsw_sa_from_sp(sp,p,lon,lat)

  ! Compute the conservative temperature
  ct = gsw_ct_from_pt(sa,pt)

  ! Compute the in situ temperature
  t_from_pt  = gsw_t_from_ct(sa,ct,p)

END FUNCTION t_from_pt

!===============================================================================
PURE FUNCTION p_from_z(dpth,xlat)
!===============================================================================
! pressure from depth from saunder's formula with eos80.
! reference: saunders,peter m., practical conversion of pressure
!            to depth., j.p.o. , april 1981.
! r millard
! march 9, 1983
! check value: p_from_z=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
! http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/ppsw.f
  IMPLICIT NONE
  REAL(r14) :: p_from_z
  REAL(r14), INTENT(IN) :: dpth
  REAL(r14), INTENT(IN) :: xlat
  REAL(r14), PARAMETER :: pi=3.141592654
  REAL(r14) :: plat, d, c1

  plat=abs(xlat*pi/180.)
  d=sin(plat)
  c1=5.92e-3+d**2*5.25e-3
  p_from_z=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6

END FUNCTION p_from_z

END MODULE gsw_pot_to_insitu
