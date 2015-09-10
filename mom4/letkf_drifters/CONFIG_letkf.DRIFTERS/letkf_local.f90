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
!
! PUBLIC TYPES:
!                 implicit none
!                 [save]
!
!                 <type declaration>
!     
! PUBLIC MEMBER FUNCTIONS:
!           <function>                     ! Description      
!
! PUBLIC DATA MEMBERS:
!           <type> :: <variable>           ! Variable description

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
! !REVISION HISTORY:
!   04/03/2014 Steve Penny created for use with OCEAN at NCEP.
!
! Designed by Prof. Stephen G. Penny
! University of Maryland, College Park
! 
!-------------------------------------------------------------------------------
! $Author: Steve Penny $
!===============================================================================

  USE common
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE letkf_obs !contains debug_hdxf_0, and nobsgrd
  USE params_letkf, ONLY: nbv, DO_NO_VERT_LOC, localization_method

  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,nid_obs) = 1.0d0

!  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,nid_obs) = RESHAPE( &      !(OCEAN) !(DO_SFCFLUXES)
!!           U      V      T      S    SSH    SST    SSS     uflx   vflx  tflux  qflux    u10    v10    t2m    q2m PRESmsl prate    dlw    dsw  longwv shrtwv
!   & (/ 1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! U !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! V !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! T !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! S !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SSH !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SST !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /) & ! SSS !(OCEAN)
!  & ,(/nv3d+nv2d,nid_obs/))
! STEVE: note: this was a sort of 'trick' to make sure that the SFC variables
! were analyzed separately. All that was needed for this to occur was to have
! the var_local vector associated with each SFC variable to be different than
! the vector associated with the ocean interior variables. !(DO_SFCFLUXES)
! As a result, I can compute adaptive inflation separately for the SFC (and
! apply none to the ocean interior).

!  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,nid_obs) = RESHAPE( &      !(OCEAN)
!!           U      V      T      S    SSH    SST    SSS     uflx   vflx  tflux  qflux    u10    v10    t2m    q2m PRESmsl prate    dlw    dsw  longwv shrtwv
!   & (/ 1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! U !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! V !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! T !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! S !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SSH !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SST !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) & ! SSS !(OCEAN)
!  & ,(/nv3d+nv2d,nid_obs/))
!!  &
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! uflx !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! vflx !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! tflux !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! qflux !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! u10 !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! v10 !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! t2m !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! q2m !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! PRESmsl !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! prate !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! dlw !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! dsw !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! longwv !(OCEAN)
!  &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)& ! shrtwv !(OCEAN)
!  & ,(/nv3d+nv2d,nid_obs/))
   !NOTE: the obs are the rows and the model variables the columns

!  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,nid_obs) = RESHAPE( &      !(OCEAN)
!!           U      V      T      S    SSH    SST    SSS                      !(OCEAN)
!   & (/ 1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! U             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! V             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! T             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! S             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! SSH           !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! SST           !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0 /)& ! SSS           !(OCEAN)
!   & ,(/nv3d+nv2d,nid_obs/))
   !NOTE: the obs are the rows and the model variables the columns
  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)

CONTAINS

SUBROUTINE obs_local(ij,ilev,var_local,hdxf,rdiag,rloc,dep,nobsl,nobstotal)
!===============================================================================
! Project global observations to local
!     (hdxf_global,dep_global,rdiag_global) -> (hdxf,dep,rdiag)
!===============================================================================
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nobstotal
  REAL(r_size),INTENT(IN) :: var_local(nid_obs)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,tvnn,iobs
  !STEVE: for (OCEAN):
  LOGICAL :: blocked_by_land     !Use these three lines to identify gulf vs. pacific points in localization
  REAL(r_size) :: xlat,xlon,lxpa
  REAL(r_size) :: olat,olon,olev,oelm,lopa
  REAL(r_size) :: f1lon,f1lat,f2lon,f2lat
  REAL(r_size) :: dist1,dist2,a,b,u,v,lu,lv
  REAL(r_size) :: ecc,ecc2,fcc,fcc2,theta,theta_rotate,sigma_min
  REAL(r_size) :: dist_zero_a,dist_zero_b,dist_min
  REAL(r_size) :: sigma_a,sigma_b
  REAL(r_size) :: minr,maxr,xdis,fcc0,fcc1,maxdN,maxdS,xrad
  REAL(r_size), PARAMETER :: cmpersec2kmperday=0.864d0 !, days = 5.0d0
  !STEVE: for initializing splits of model grid:
  LOGICAL :: splitlon = .false. ! initialize
  LOGICAL :: splitlonL = .false., splitlonR = .false. ! initialize
  LOGICAL :: splitlat = .false. ! initialize
  LOGICAL :: fulllon = .false. ! initialize
  !STEVE:
  LOGICAL :: dodebug = .false.
  !debug TEST:
  INTEGER :: nn_old,imin_old,imax_old,jmin_old,jmax_old
  REAL(r_size) :: minlon_old,maxlon_old,minlat_old,maxlat_old

  !-----------------------------------------------------------------------------
  ! INITIALIZE
  !-----------------------------------------------------------------------------
  if( nobs > 0 ) then
    ALLOCATE(nobs_use(nobs))
  endif
  nn = 0

  !-----------------------------------------------------------------------------
  ! TEST: For comparison (and debugging) call old Miyoshi version:
  !-----------------------------------------------------------------------------
  if (dodebug) CALL obs_local_setup_old(ij,nn_old,minlon_old,maxlon_old,minlat_old,maxlat_old,imin_old,imax_old,jmin_old,jmax_old)

  !-----------------------------------------------------------------------------
  ! Set model coordinates
  !-----------------------------------------------------------------------------
  xlat = lat1(ij)
  xlon = lon1(ij)
  if (dodebug) WRITE(6,*) "---------------------------------- ij = ", ij
  
  !
  ! data search (NOTE: this only wraps in the longitude direction)
  ! Use bounding box first to reduce obs considered.
  !
  minlon = xlon - dlon_zero(ij)
  maxlon = xlon + dlon_zero(ij)
  minlat = xlat - dlat_zero !STEVE: this should be changed to be (ij) dependent
  maxlat = xlat + dlat_zero

  if (dodebug) then 
    WRITE(6,*) "---minlon = ", minlon
    WRITE(6,*) "---maxlon = ", maxlon
    WRITE(6,*) "---minlat = ", minlat
    WRITE(6,*) "---maxlat = ", maxlat
    
    WRITE(6,*) "+++minlon = ", minlon_old
    WRITE(6,*) "+++maxlon = ", maxlon_old
    WRITE(6,*) "+++minlat = ", minlat_old
    WRITE(6,*) "+++maxlat = ", maxlat_old

    WRITE(6,*) "nn_old   = ", nn_old
    WRITE(6,*) "imin_old = ", imin_old
    WRITE(6,*) "imax_old = ", imax_old
    WRITE(6,*) "jmin_old = ", jmin_old
    WRITE(6,*) "jmax_old = ", jmax_old
  endif

  !-----------------------------------------------------------------------------
  ! Handle the case when the bounding box is split of the edge of the grid
  !-----------------------------------------------------------------------------
  !STEVE: ISSUE
  !STEVE: we need handling of the fold in the arctic due to the Murray tripolar grid.
  fulllon = .false.
  splitlonL = .false.
  splitlonR = .false.
  splitlat = .false.
  if (maxlon - minlon > lonf - lon0) then 
    fulllon = .true.
    if (dodebug) then 
      WRITE(6,*) "FullLon."
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "xlat = ", xlat
    endif
    minlon = lon0
    maxlon = lonf
  else

    if (minlon < lon0) then ! Local region overlaps on left
      if (dodebug) WRITE(6,*) "split minlon = ", minlon
      splitlonL = .true.
      minlon = REAL(lonf - abs(lon0 - minlon) + wrapgap,r_size)
    endif

    if (maxlon > lonf) then ! Local region overlaps on right
      if (dodebug) WRITE(6,*) "split maxlon = ", maxlon
      if (splitlonL) WRITE(6,*) "WARNING: lon already split. ij = ", ij
      splitlonR = .true.
      maxlon = REAL(lon0 + abs(maxlon - lonf) - wrapgap,r_size)
    endif 

    if (minlat < lat0) then ! Local region overlaps on bottom
      splitlat= .true.
      minlat = lat0
    endif

    if (maxlat > latf) then ! Local region overlaps on top
      if (splitlat) WRITE(6,*) "WARNING: lat already split. ij = ", ij
      splitlat = .true.
      maxlat = latf
    endif 

  endif

  if (splitlat) then
    ! Too close to pole, need to work out the details of this
    ! For now, cut off the local region at the pole 
    WRITE(6,*) "letkf_local.f90::obs_local: Grid point is close to pole."
  endif

  ! Find jmin and jmax
  do jmin=1,nlat-2
    if(minlat < lat(jmin+1)) exit
  enddo
  do jmax=1,nlat-2
    if(maxlat < lat(jmax+1)) exit
  enddo
  nn = 1

  if ( (splitlonL .or. splitlonR) .and. .not. fulllon) then
    if (dodebug) WRITE(6,*) "splitlonL = ", splitlonL
    if (dodebug) WRITE(6,*) "splitlonR = ", splitlonR

    ! Do left
    imin = 1
    do imax=1,nlon-1
      if(maxlon < lon(imax+1)) exit
    enddo
    if (dodebug) then
      WRITE(6,*) "splitlon left: nn = ", nn
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif
    if( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)

    ! Do right
    do imin=1,nlon-1
      if(minlon < lon(imin+1)) exit
    enddo
    imax = nlon
    if (dodebug) then
      WRITE(6,*) "splitlon right: nn = ", nn
      WRITE(6,*) "xlon = ", xlon
      WRITE(6,*) "minlon = ", minlon
      WRITE(6,*) "maxlon = ", maxlon
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
    endif
    if( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)

  elseif ( .not. fulllon ) then
    ! The local range is entirely contained within the global domain
    ! Find the min and max latitudes
    do jmin=1,nlat-2
      if(minlat < lat(jmin+1)) exit
    enddo
    do jmax=jmin,nlat-2
      if(maxlat < lat(jmax+1)) exit
    enddo

    ! Find the min and max longitudes
    do imin=1,nlon-1
      if(minlon < lon(imin+1)) exit
    enddo
    do imax=imin,nlon-1
      if(maxlon < lon(imax+1)) exit
    enddo

    if( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)

  else
    ! The domain is the entire circle around the globe (zonally)
    ! Find the min and max latitudes
    do jmin=1,nlat-2
      if(minlat < lat(jmin+1)) exit
    enddo
    do jmax=jmin,nlat-2
      if(maxlat < lat(jmax+1)) exit
    enddo

    if (dodebug) then
      WRITE(6,*) "Doing FullLon..."
    endif
    imin=1
    imax=nlon
    if( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  endif
  nn = nn-1
  !STEVE: end rewrite

  ! If no observations remain, then we're done.
  if(nn < 1) then
    nobsl = 0
    RETURN
  endif

  !----------------------------------------------------------------------------
  ! Cycle through all observations to identify which should be kept
  !----------------------------------------------------------------------------
  !STEVE: most of the localization section has been completely edited for (OCEAN)
  !STEVE: This should eventually be replaced with a tree-based search algorithm (e.g. R-Tree, kd-tree)
  !STEVE: (future) use custom localization with CGAL/BOOST algorithms
  nobsl = 0
  do n=1,nn
    !---------------------------------------------------------------------------
    ! Observational localization Distance Cutoff
    !---------------------------------------------------------------------------
    olat = obslat(nobs_use(n))
    olon = obslon(nobs_use(n))
    olev = obslev(nobs_use(n))
    oelm = obselm(nobs_use(n))

    !---------------------------------------------------------------------------
    ! variable localization
    !---------------------------------------------------------------------------
    !STEVE: ISSUE: this should be handled in a more generic way so new
    !              observations don't requrie editing here.
    SELECT CASE(NINT(oelm))
    CASE(id_u_obs)
        iobs=1
    CASE(id_v_obs)
        iobs=2
    CASE(id_t_obs)
        iobs=3
    CASE(id_s_obs)   !(OCEAN)
        iobs=4
    CASE(id_ssh_obs) !(OCEAN)
        iobs=5
    CASE(id_sst_obs) !(OCEAN)
        iobs=6
    CASE(id_sss_obs) !(OCEAN)
        iobs=7
    CASE(id_eta_obs) !(OCEAN)
        iobs=8
    CASE DEFAULT
        WRITE(6,*) "letkf_local.f90 :: there is no variable localization for obs-type :: ", oelm
        WRITE(6,*) "letkf_local.f90 :: exitING..."
        STOP(95)
    END SELECT
    if(var_local(iobs) < TINY(var_local)) CYCLE  ! Skip obs that are set to "0" impact on this model variable.

    !---------------------------------------------------------------------------
    ! vertical localization
    !---------------------------------------------------------------------------
    !STEVE: make vertical localization depth dependent
    !STEVE: Could alternatively treat surface obs as if they occur at depth 0, regardless of model's top level
    !STEVE: ISSUE: need to make this more general, define vertical localization types for specific obs elsewhere
    if(NINT(oelm) == id_ssh_obs) then
      dlev = 0.0d0 ! allow all levels to be influenced by ssh
    else if(NINT(oelm) == id_sst_obs) then
      dlev = ABS( olev - lev(ilev) )
      !WRITE(6,*) "obs_local:: SST ob, dlev = ",dlev
    else if(NINT(oelm) == id_sss_obs) then
      dlev = ABS( olev - lev(ilev) )
    else if(NINT(oelm) == id_u_obs .or. &
            NINT(oelm) == id_v_obs .or. &
            NINT(oelm) == id_t_obs .or. &
            NINT(oelm) == id_s_obs   ) then
      dlev = ABS( olev - lev(ilev) )
    else
      dlev = 0.0d0
    endif
    if(dlev > dist_zerov) CYCLE

    !---------------------------------------------------------------------------
    ! horizontal localization
    !---------------------------------------------------------------------------
    horizontal_localization : if (localization_method .eq. 1) then
      !STEVE: make horizontal localization latitude dependent
      ! STEVE: make sigma_obs a linear/lookup function of latitude
      dist_min = sigma_obs0 * (SQRT(10.0d0/3.0d0) * 2.0d0)
      minr = dist_min
      maxr = dist_zero
      maxdN = abs( (90.0d0)*(pi/180.0d0)*re)
      maxdS = abs((-90.0d0)*(pi/180.0d0)*re)

      ! Shrink radius far from equator
      ! Shrink foci far from equator
      xdis = abs(xlat)*(pi/180.0d0)*re
      if (xlat >= 0) then
        xrad = (1 - xdis/maxdN)*maxr + (xdis/maxdN)*minr
      else
        xrad = (1 - xdis/maxdS)*maxr + (xdis/maxdS)*minr
      endif

      sigma_a = xrad / (SQRT(10.0d0/3.0d0) * 2.0d0)
      sigma_b = sigma_a
      dist_zero_a = xrad

      !-------------------------------------------------------------------------
      ! Compute the great circle distance between the grid point and the observation 
      !-------------------------------------------------------------------------
      CALL com_distll_1(olon,olat,xlon,xlat,dist)

      if (.false.) then !(modulo(n,10000) .eq. 0) then
        !STEVE: DEBUG: print out the important values for a selection of obs
        print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
        print *, "DEBUG letkf_local.f90..."
        print *, "n = ", n
        print *, "dist_min = ", dist_min
        print *, "sigma_obs0 = ", sigma_obs0
        print *, "minr = ", minr
        print *, "maxr = ", maxr
        print *, "maxdN = ", maxdN
        print *, "maxdS = ", maxdS
        print *, "xlat = ", xlat
        print *, "xdis = ", xdis
        print *, "xrad = ", xrad
        print *, "sigma_a = ", sigma_a
        print *, "sigma_b = ", sigma_b
        print *, "sigma_obs = ", sigma_obs
        print *, "dist_zero_a = ", dist_zero_a
        print *, "dist = ", dist
!       print *, "dist > dist_zero_a = ", dist > dist_zero_a
        print *, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
      endif

      if(dist > dist_zero_a) CYCLE  !Points are outside of radius

    else !STEVE: use the original (default) localization

      CALL com_distll_1(olon,olat,xlon,xlat,dist)
      if(dist > dist_zero ) CYCLE

    endif horizontal_localization

    ! STEVE: ADD CHECK FOR OBSERVATIONS THAT ARE OCCLUDED BY LAND!
    ! STOLE atlpac from SODA, in the future, implement a general method
    ! for occluding points that can't be reached by a SPM in the local range
    blocked_by_land = .false.
    CALL atlpac(xlat,xlon,lxpa)
    CALL atlpac(olat,olon,lopa)
    IF (  lxpa .eq. 1 .and. lopa .eq. 2 &        !in caribbean, pacific and atlantic
       .or. lxpa .eq. 2 .and. lopa .eq. 1 &      !in caribbean, atlantic and pacific
         ) blocked_by_land = .true.
    IF (blocked_by_land) CYCLE                   !STEVE: the two points are in atl and pac

    !STEVE: ISSUE: update this to read in a table with basin ids to 
    !              only assimilate data in current or adjacent basin

!--------------------------- End of Observation Culling -----------------------!

    !---------------------------------------------------------------------------
    ! Collect all identifed local observations into data structures for output
    !---------------------------------------------------------------------------
    nobsl = nobsl + 1
    hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
    dep(nobsl)    = obsdep(nobs_use(n))
    tmperr=obserr(nobs_use(n))
    rdiag(nobsl) = tmperr * tmperr
    if (ALLOCATED(obs_useidx)) then
      obs_useidx(nobsl) = nobs_use(n)
    endif

    !---------------------------------------------------------------------------
    ! Observational localization (weighting)
    !---------------------------------------------------------------------------
    !
    ! Note: var_local scales the localization weighting based on the parameter
    ! set in letkf_tools.f90. A row corresponding to the model parameter is
    ! input to this subroutine, and the column indicates the proportion of that
    ! type of observation to use.
    observation_localization : if (localization_method .eq. 1 ) then
      if (DO_NO_VERT_LOC) then
        rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
      else
        rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                   & * var_local(iobs)
      endif

    elseif (localization_method .eq. 0) then
      !STEVE: this is R^2 and the R localization gaussian function
      rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                                                   & * var_local(iobs)
    else

      print *, "ERROR:: Localization method not supported. localization_method = ", localization_method 
      STOP(3)

    endif observation_localization

    !STEVE: debugging error... keeping this around to prevent it 
    !                          from occuring again and going unnoticed
    if (rloc(nobsl) > 1) then
      WRITE(6,*) "rloc(nobsl) > 1 !"
      WRITE(6,*) "localization_method = ", localization_method
      WRITE(6,*) "nobsl = ", nobsl
      WRITE(6,*) "dist = ", dist
      WRITE(6,*) "sigma_a = ", sigma_a
      WRITE(6,*) "sigma_obs = ", sigma_obs
      WRITE(6,*) "var_local(iobs) = ", var_local(iobs)
      WRITE(6,*) "rloc(nobsl) = ", rloc(nobsl)
      WRITE(6,*) "letkf_local.f90:: exitING.."
      STOP 93 
    endif

  enddo

  !----------------------------------------------------------------------------
  !STEVE: this should never happen, if it does something went wrong
  !----------------------------------------------------------------------------
  if( nobsl > nobstotal ) then
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN,TVNN=', ij, nn, tvnn
    STOP 99
  endif
 
  if( nobs > 0 ) then
    DEALLOCATE(nobs_use)
  endif

END SUBROUTINE obs_local

PURE SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!===============================================================================
! Identify the observations within the local region
!===============================================================================
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  ! Cycle through each latitude
  do j=jmin,jmax
    ! Find the number of accumulated obs at the bottom of the range
    if(imin > 1) then
      ib = nobsgrd(imin-1,j)+1
    else
      ! Wrap around to the previous latitude at the last longitude
      if(j > 1) then
        ib = nobsgrd(nlon,j-1)+1
      else
        ib = 1
      endif
    endif
    ! Find the number of accumulated obs at the top of the range
    ie = nobsgrd(imax,j)

    ! Subtract to get the number of obs in this region
    n = ie - ib + 1

    if(n == 0) CYCLE !there are no obs here

    do ip=ib,ie
      if(nn > nobs) then
!       WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!       stop 1  !STEVE: (added)
      endif
      ! Index for observation used
      nobs_use(nn) = ip
      ! Count up the total obs used so far
      nn = nn + 1
    enddo
  enddo

  RETURN
END SUBROUTINE obs_local_sub

SUBROUTINE obs_local_setup_old(ij,nn,minlon,maxlon,minlat,maxlat,imin,imax,jmin,jmax)
!===============================================================================
! Original Miyoshi code
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!===============================================================================
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij
  INTEGER,INTENT(OUT) :: nn
  REAL(r_size), INTENT(OUT) :: minlon,maxlon,minlat,maxlat
  INTEGER, INTENT(OUT) :: imin,imax,jmin,jmax
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: n,nobsl,im
  !STEVE: for (OCEAN):
  REAL(r_size) :: xlat,xlon
  LOGICAL :: dodebug = .false.

!
! INITIALIZE
!
  if( nobs > 0 ) then
    ALLOCATE(nobs_use(nobs))
  endif
!
! data search
!
  minlon = lon1(ij) - dlon_zero(ij)
  maxlon = lon1(ij) + dlon_zero(ij)
  minlat = lat1(ij) - dlat_zero
  maxlat = lat1(ij) + dlat_zero

  do jmin=1,nlat-2
    if(minlat < lat(jmin+1)) exit
  enddo
  do jmax=1,nlat-2
    if(maxlat < lat(jmax+1)) exit
  enddo
  nn = 1
  if(minlon >= lon0 .AND. maxlon <= lonf) then
    do imin=1,nlon-1
      if(minlon < lon(imin+1)) exit
    enddo
    do imax=1,nlon-1
      if(maxlon < lon(imax+1)) exit
    enddo
    if( nobs > 0 ) &
    & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  else if(minlon >= lon0 .AND. maxlon > lonf) then
    do imin=1,nlon-1
      if(minlon < lon(imin+1)) exit
    enddo
    maxlon = maxlon - 360.0d0
    if(maxlon > lonf) then
      imin = 1
      imax = nlon
      if( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    else
      do imax=1,nlon-1
        if(maxlon < lon(imax+1)) exit
      enddo
      if(imax > imin) then
        imin = 1
        imax = nlon
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      else
        imin = 1
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        do imin=1,nlon-1
          if(minlon < lon(imin+1)) exit
        enddo
        imax = nlon
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      endif
    endif
  else if(minlon < lon0 .AND. maxlon <= lonf) then
    do imax=1,nlon-1
      if(maxlon < lon(imax+1)) exit
  enddo
    minlon = minlon + 360.0d0
    if(minlon < lon0) then
      imin = 1
      imax = nlon
      if( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    else
      do imin=1,nlon-1
        if(minlon < lon(imin+1)) exit
      enddo
      if(imin < imax) then
        imin = 1
        imax = nlon
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      else
        imin = 1
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
        do imin=1,nlon-1
          if(minlon < lon(imin+1)) exit
        enddo
        imax = nlon
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      endif
    endif
  else
    maxlon = maxlon - 360.0d0
    minlon = minlon + 360.0d0
    if(maxlon > lonf .OR. minlon < lon0) then
      imin = 1
      imax = nlon
      if( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    else
      do imin=1,nlon-1
        if(minlon < lon(imin+1)) exit
      enddo
      do imax=1,nlon-1
        if(maxlon < lon(imax+1)) exit
      enddo
      if(imin > imax) then
        imin = 1
        imax = nlon
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      else
        if( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
      endif
    endif
  endif
  nn = nn-1
  if(nn < 1) then
    nobsl = 0
    RETURN
  endif
END SUBROUTINE obs_local_setup_old

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
    if(xlon.lt.285.) then
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
