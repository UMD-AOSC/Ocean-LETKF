MODULE letkf_local
!===============================================================================
! MODULE: letkf_local
! 
! USES:
!   use common
!   use common_mpi
!   use common_mom4
!   use common_mpi_mom4
!   use common_letkf
!   use letkf_obs
!   use params_letkf
!   use vars_letkf
!
! PUBLIC TYPES:
!                 implicit none
!                 [save]
!
!                 <TYPE declaration>
!     
! PUBLIC MEMBER FUNCTIONS:
!           <function>                     ! Description      
!
! PUBLIC DATA MEMBERS:
!           <TYPE> :: <variable>           ! Variable description

! DESCRIPTION: 
!   This module contains the subroutines necessary to perform localization.
!   It is an offshoot of letkf_tools, which previously contained
!   all localization routines. The localization algorithm has been significantly
!   updated versus Miyoshi's original approach, and special modifications
!   have been made specific to the ocean domain.
!
!   The longer-term purpose of isolating these routines in
!   an independent module is to allow for further development of localization
!   approaches, in particular those utilizing computational geometry tools
!   such as:
!   kd-tree, kNN search and range search
!   Graph representation and A* search
!   
!   In this implementation, the kd-tree has been applied by Travis Sluka
!   for range-based searches.
!
! !REVISION HISTORY:
!   03/30/2016 Steve Penny updated with kdtree code developed by Travis Sluka
!   04/03/2014 Steve Penny created for use with OCEAN at NCEP.
!
! Originally designed by Prof. Stephen G. Penny
! University of Maryland, College Park
! 
!-------------------------------------------------------------------------------
! $Authors: Steve Penny, Travis Sluka $
!===============================================================================

  USE common
  USE common_obs_mom4
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE letkf_obs !contains debug_hdxf_0, and nobsgrd
  USE vars_letkf,   ONLY: var_local_n2n
  USE params_obs,   ONLY: nid_obs, id_sst_obs, id_ssh_obs, id_sss_obs
  USE vars_obs,     ONLY: obslon, obslat, obslev, obselm, obsdep, obshdxf, obserr
  USE params_letkf, ONLY: sigma_obs, sigma_obs0, DO_NO_VERT_LOC
  USE kdtree

  IMPLICIT NONE

  PUBLIC :: obs_local
  
  PRIVATE

  TYPE(KD_ROOT), SAVE :: kdtree_root
  INTEGER, SAVE       :: initialized = 0
  REAL(r_size), SAVE  :: prev_lon, prev_lat      !STEVE: for checking if this longitude was last searched

  INTEGER, SAVE, ALLOCATABLE      :: idx(:)      !STEVE: index of the observations that are found by kd_search
  REAL(r_size), SAVE, ALLOCATABLE :: dist(:)     !STEVE: distance from the center grid point
  INTEGER, SAVE :: nn                            !STEVE: total number of local observations found by kd_search

  !! Localization
  !! ------------------------------------------------------------
  REAL(r_size), DIMENSION(2) :: sigma_ocn_h
  REAL(r_size) :: sigma_ocn_v
  
  REAL(r_size), DIMENSION(2) :: sigma_atm_h
  REAL(r_size) :: sigma_atm_v

CONTAINS

SUBROUTINE obs_local(ij,ilev,mlev,var_local,hdxf,rdiag,rloc,dep,nobsl,nobstotal)
  !===============================================================================
  ! Project global observations to local
  !     (hdxf_global,dep_global,rdiag_global) -> (hdxf,dep,rdiag)
  !===============================================================================
  INTEGER,      INTENT(IN)  :: ij,ilev,mlev,nobstotal
  REAL(r_size), INTENT(IN)  :: var_local(nid_obs)
  REAL(r_size), INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size), INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size), INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size), INTENT(OUT) :: dep(nobstotal)
  INTEGER,      INTENT(OUT) :: nobsl

  !! things used internally for hz loclzlization
  REAL(r_size) :: sigma_atm_h_ij, sigma_ocn_h_ij
  REAL(r_size) :: sigma_atm_v_ij, sigma_ocn_v_ij
  REAL(r_size) :: dist_zero_h, sigma_h
  REAL(r_size) :: dist_zero_v, sigma_v
  REAL(r_size) :: dist_zero_ocn_h, dist_zero_ocn_v
  REAL(r_size) :: dist_zero_atm_h, dist_zero_atm_v
  REAL(r_size) :: sigma_max_h_d0

  REAL(r_size) :: loc_h, loc_v
  INTEGER :: n,i,j,k 
  REAL(r_size):: dlev
  REAL(r_size) :: cutoff_scaling = SQRT(10.0d0/3.0d0) * 2.0d0

  !STEVE: debug
  LOGICAL :: dodebug = .false.

  !-----------------------------------------------------------------------------
  ! Initialize the KD-tree
  !-----------------------------------------------------------------------------
  !STEVE: NOTE: this is replicated for each process
  if (initialized == 0) then
    WRITE(6,*) "Initializing the obs_local()"
    initialized = 1
    call kd_init(kdtree_root, obslon, obslat)
    WRITE(6,*) "Done constructing KD search tree."
    WRITE(6,*) "nobstotal = ", nobstotal
    ALLOCATE(dist(nobstotal))
    ALLOCATE(idx(nobstotal))
    prev_lon = -1e10
    prev_lat = -1e10

    ! Enable initialization via namelist:
    sigma_ocn_h = (/ sigma_obs, sigma_obs0 /)
    sigma_ocn_v = sigma_obsv
    sigma_atm_h   = (/ sigma_obs, sigma_obs0 /)
    sigma_atm_v   = sigma_obsv
    if (DO_NO_VERT_LOC) then
      sigma_ocn_v = -1
      sigma_atm_v = -1
    endif
  endif

  !-----------------------------------------------------------------------------
  ! determine the maximum search radius for the initial obs search 
  ! (NOTE: this is gridpoint-dependent, so must be updated every time ij changes)
  !-----------------------------------------------------------------------------
  sigma_ocn_h_ij = (1.0d0-(abs(lat1(ij))/90.0d0))*(sigma_ocn_h(1)-sigma_ocn_h(2))+sigma_ocn_h(2)
  sigma_ocn_v_ij = sigma_ocn_v
  sigma_atm_h_ij = (1.0d0-(abs(lat1(ij))/90.0d0))*(sigma_atm_h(1)-sigma_atm_h(2))+sigma_atm_h(2)    
  sigma_atm_v_ij = sigma_atm_v
  dist_zero_ocn_h = sigma_ocn_h_ij * cutoff_scaling
  dist_zero_ocn_v = sigma_ocn_v_ij * cutoff_scaling
  dist_zero_atm_h = sigma_atm_h_ij * cutoff_scaling
  dist_zero_atm_v = sigma_atm_v_ij * cutoff_scaling
  sigma_max_h_d0 = max(sigma_ocn_h_ij, sigma_atm_h_ij)  * cutoff_scaling

  !-----------------------------------------------------------------------------
  ! query the KD-tree
  !-----------------------------------------------------------------------------
  if (lon1(ij) /= prev_lon .or. lat1(ij) /= prev_lat) then
    call kd_search(kdtree_root, obslon, obslat, (/lon1(ij),lat1(ij)/),sigma_max_h_d0, idx, dist, nn)
    prev_lon = lon1(ij)
    prev_lat = lat1(ij)
    if (dodebug) then
      WRITE(6,*) "letkf_local.f90 ------------------------------------ ij = ", ij
      WRITE(6,*) "sigma_max_h_d0 = ", sigma_max_h_d0
      WRITE(6,*) "dist_zero_ocn_h = ", dist_zero_ocn_h
      WRITE(6,*) "dist_zero_ocn_v = ", dist_zero_ocn_v
      WRITE(6,*) "Observations found by kd_search = "
      WRITE(6,*) "idx,lon,lat,lev,dist"
      do n=1,nn
        WRITE(6,*) idx(n), obslon(idx(n)), obslat(idx(n)), obslev(idx(n)), dist(n)
      enddo
    endif
  endif

  !-----------------------------------------------------------------------------
  !  For each observation found in the local radius, 
  !  scale the weight of the observation error based on distance from gridpoint.
  !  Outside of a maximum radius, set the obs impact to zero (i.e. remove the ob)
  !-----------------------------------------------------------------------------
  nobsl=0        
  do n = 1, nn
    loc_h = 0.0
    loc_v = 0.0
!   if (dodebug) WRITE(6,*) "n = ", n
       
    if(nobsl >= nobstotal) EXIT

    !!------------------------------          
    !! atmospheric observation
    !!------------------------------
    !STEVE: put this back in when necessary

    !!------------------------------          
    !! ocean observation
    !!------------------------------
    sigma_h = sigma_ocn_h_ij
    dist_zero_h = dist_zero_ocn_h
    dist_zero_v = dist_zero_ocn_v

    if (sigma_ocn_v < 0) then
      !! vertical localization is off
      dlev = 0
      sigma_v = HUGE(1.0d0)
      dist_zero_v= HUGE(1.0d0)
    else
      dlev = abs(lev(ilev)-obslev(idx(n)))
      sigma_v = sigma_ocn_v
    endif             

    !! horizontal localization cutoff
    if (dist(n) > dist_zero_h) CYCLE

    !! vertical localization cutoff
    if (dlev > dist_zero_v) CYCLE

    !! Skip surface obs below the mixed layer:
    if (ilev > mlev .and. obselm(idx(n))==id_sst_obs) CYCLE
!   if (ilev > mlev .and. obselm(idx(n))==id_ssh_obs) CYCLE
!   if (ilev > mlev .and. obselm(idx(n))==id_sss_obs) CYCLE
       
    !! Else, use this observation!
    nobsl = nobsl+1              
    if (dodebug) WRITE(6,*) "nobsl = ", nobsl

    if (DO_NO_VERT_LOC) then
      rloc(nobsl) = EXP( -0.5d0 * (dist(n)/sigma_h)**2 )
    else
      rloc(nobsl) = EXP(-0.5d0 * ((dist(n)/sigma_h)**2 + (dlev/sigma_v)**2))
    endif
    hdxf(nobsl,:) = obshdxf(idx(n),:)
    dep(nobsl)  = obsdep(idx(n))
    rdiag(nobsl) = obserr(idx(n))**2

    if (dodebug) WRITE(6,*) "rloc(nobsl)  = ", rloc(nobsl)
    if (dodebug) WRITE(6,*) "rdiag(nobsl) = ", rdiag(nobsl)
    if (dodebug) WRITE(6,*) "dep(nobsl)   = ", dep(nobsl)
  enddo

END SUBROUTINE obs_local


!(OCEAN) STEVE: add checks for atlantic/pacific basin boundary
PURE SUBROUTINE atlpac (xlat, xlon_in, lxap)
REAL(r_size), INTENT(IN) :: xlat, xlon_in
REAL(r_size), INTENT(OUT) :: lxap
REAL(r_size) :: xlon

! Ensure the longitude is specified on a 0-360 grid:
xlon = modulo(xlon_in,360.0)

! STEVE: Stolen from SODA: 
! ISSUE: use until we have a general method for managing land-blocked ocean basins...
!=================================================================
! X. Cao 12/9/99
!
!   to make a mark to the location of a point in Caribbean area
! (xlat.gt.-2..and.xlat.lt.32..and.xlon.gt.245..and.xlon.lt.295.)
! to indicate if it is in Atlantic ocean (1) or in Pacific ocean (2)
! or out of Caribbean area (0)
!=================================================================
!
  lxap=0
!
! -- Atlantic? Pacific?
!
  if(xlat.gt.-2.0.and.xlat.le.8.5) then
!   if(xlon.lt.285.) then
    if(xlon.lt.285. .AND. xlon .gt. 30.0) then !STEVE: this bug was identified by Guillaume V., causing a discontinuity at 0ºlongitude.fixed 3/16/16
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.8.5.and.xlat.le.15.5) then
    if(xlon.lt.276.) then
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.15.5.and.xlat.le.19.5) then
    if(xlon.lt.270.) then
      lxap=2
    else
      lxap=1
    endif
  endif

  if(xlat.gt.19.5.and.xlat.le.32.0) then
    if(xlon.lt.258.) then
      lxap=2
    else
      lxap=1
    endif
  endif

END SUBROUTINE atlpac

END MODULE letkf_local
