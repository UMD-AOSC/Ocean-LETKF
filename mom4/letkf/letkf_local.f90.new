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
  USE params_letkf, ONLY: nbv, DO_NO_VERT_LOC, DO_MLD, localization_method, DO_IRREG_GRID
  USE vars_letkf,   ONLY: var_local, var_local_n2n, get_iobs

  PUBLIC :: obs_local

  PRIVATE

  TYPE local_panel
    INTEGER:: panel_index
    REAL(r_size) :: w !minlon
    REAL(r_size) :: e !maxlon
    REAL(r_size) :: s !minlat
    REAL(r_size) :: n !maxlat
    CHARACTER(8) :: panel_label !e.g. 'low left'
    CHARACTER(1) :: split_label !=('NSEW')
  END TYPE local_panel

  REAL(r_size), PARAMETER :: npole=90.0d0
  REAL(r_size), PARAMETER :: spole=-90.0d0

  INTEGER, PARAMETER :: npmax=4  ! Maximum number of panels
  REAL(r_size), SAVE :: wbdy, ebdy

CONTAINS


SUBROUTINE obs_local(ij,ilev,mlev,var_local,hdxf,rdiag,rloc,dep,nobsl,nobstotal)
!===============================================================================
! Project global observations to local
!     (hdxf_global,dep_global,rdiag_global) -> (hdxf,dep,rdiag)
!===============================================================================
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,mlev,nobstotal
  REAL(r_size),INTENT(IN) :: var_local(nid_obs)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size), DIMENSION(4) :: minlons,maxlons,minlats,maxlats
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
  !STEVE:
  LOGICAL :: dodebug = .false.

  TYPE (local_panel), DIMENSION(npmax) :: panel
  INTEGER :: npanels = 0         ! Split the Minimum Bounding Box into as many panels as necessary
  INTEGER :: ip

  !-----------------------------------------------------------------------------
  ! INITIALIZE
  !-----------------------------------------------------------------------------
  if( nobs > 0 ) then
    ALLOCATE(nobs_use(nobs))
  else
    WRITE(6,*) "letkf_local.f90::obs_local NO OBSERVATIONS!"
    STOP(695)
  endif
  wbdy=lon0
  ebdy=lonf

  !-----------------------------------------------------------------------------
  ! Set model coordinates - this is the gridpoint we're analyzing at the moment
  !-----------------------------------------------------------------------------
  xlat = lat1(ij)
  xlon = lon1(ij)
  if (dodebug) WRITE(6,*) "obs_local::---------------------------------- ij = ", ij
  if (dodebug) WRITE(6,*) "obs_local::-------------------------------- xlat = ", xlat
  if (dodebug) WRITE(6,*) "obs_local::-------------------------------- xlon = ", xlon
  
  !
  ! data search
  ! Use minimum bounding box (MBB) first to reduce obs considered.
  !
  ! Split bounding box into 2-4 panels if it overlaps boundaries
  npanels    = 1
  panel(1)%e = xlon - dlon_zero(ij)
  panel(1)%w = xlon + dlon_zero(ij)
  panel(1)%s = xlat - dlat_zero !STEVE: this should be changed to be (ij) dependent
  panel(1)%n = xlat + dlat_zero
  panel(1)%panel_label = 'global'

  ! Identify whether panels should be split
  CALL split_panels(panel,npanels,'tripolar')

  ! Loop through panels and find local observations in each panel
  nn=0
  do ip=1,npanels
    if (dodebug) WRITE(6,*) "ip/npanels = ", ip, "/", npanels

    ! Get indices of each panel for identifying local observations
    CALL find_mbb_idx(panel(ip)%e,panel(ip)%w,panel(ip)%s,panel(ip)%n,imin,imax,jmin,jmax)

    ! Find local observations
    if( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)

    if (dodebug) then
      WRITE(6,*) "imin = ", imin
      WRITE(6,*) "imax = ", imax
      WRITE(6,*) "jmin = ", jmin
      WRITE(6,*) "jmax = ", jmax
      WRITE(6,*) "nn   = ", nn
    endif
  enddo
  ! (nn now holds the number of local obs)

  ! If no observations remain, then we're done.
  if(nn < 1) then
    if (dodebug) then
      WRITE(6,*) "ip/npanels = ", ip, "/", npanels
      WRITE(6,*) "WARNING! no observations remaining after call to obs_local_sub (using panels). RETURNING..."
    endif
    nobsl = 0
    RETURN
  endif

  !----------------------------------------------------------------------------
  ! Cycle through all observations that are remaining after Minimum Bounding Box
  ! culling to identify which should be kept based on a geodetic metric
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
    CALL get_iobs(oelm,iobs) !STEVE: get the the index in the var_local for this ob type (in vars_letkf.f90)
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
      if (DO_MLD) then
        if (ilev <= mlev) then
          if (DO_NO_VERT_LOC) dlev = 0.0d0 !STEVE: this will apply SST obs with equal weight throughout mixed layer
        else
          dlev = HUGE(1.0d0) !STEVE: this will skip all SST obs below the mlixed layer depth
        endif
      endif
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
    observation_localization : if (localization_method .eq. 1 ) then ! Latitude-dependent localization
      if (DO_NO_VERT_LOC) then
        rloc(nobsl) =EXP(-0.5d0 * (dist/sigma_a)**2) * var_local(iobs)
      else
        rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_a)**2 + (dlev/sigma_obsv)**2)) &
                                                   & * var_local(iobs)
      endif

    elseif (localization_method .eq. 0) then  ! Constant localization
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
      WRITE(6,*) "letkf_local.f90:: EXITING..."
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


SUBROUTINE split_panels(panel,np,gtype)
  TYPE(local_panel), DIMENSION(npmax), INTENT(INOUT) :: panel
  INTEGER, INTENT(INOUT) :: np                  ! number of panels
  CHARACTER(*), INTENT(IN) :: gtype          ! Grid type (e.g. 'tripolar')
  INTEGER :: npnew,npin
  LOGICAL :: dodebug = .false.

  ! if gtype == 'tripolar':
  
! STEVE: almost all of them are this:
! if (dodebug) then
!   WRITE(6,*) "split_panels:: init"
!   CALL print_panel(panel,np)
! endif
 
  ! If the box intersects with BOTH the north and south boundaries, then 
  ! this isn't really local anymore, but still, we'll cover it as a case
  npin=np
  npnew=np
  if (panel(1)%n > npole .and. panel(1)%s < spole) then
    ! Keep same number of panels as input
    npnew=np

    do ip=1,np
      panel(ip)%n = npole
      panel(ip)%s = spole
    enddo

    if (dodebug) then
      WRITE(6,*) "split_panels:: BOTH north and south split"
      CALL print_panel(panel,np)
    endif
 
  ! If the box intersects with the northern boundary, split the panel into two
  elseif (panel(1)%n > npole) then ! Local region overlaps on top
    ! split the panel into two
    npnew=np*2
    do ip=1,np  !(assuming max 4)
      if (np+ip>npmax) then
        WRITE(6,*) "Max number of panels is: ", npmax
        WRITE(6,*) "You have: ", np+ip
        STOP(98) 
      endif
      panel(np+ip) = panel(ip) 
      
      ! Update the minimum latitude:
      panel(np+ip)%s = npole - (panel(ip)%n - npole)

      ! Update the maximum latitude:
      panel(np+ip)%n = npole

      ! Update the minimum longitude:
      !STEVE: this is not quite exact...
      rw = abs(lon(NINT(nlon/2.0)) - panel(ip)%e)
      panel(np+ip)%w = lon(NINT(nlon/2.0)) + rw

      ! Update the maximum longitude:
      lw = abs(panel(ip)%w - wbdy)
      panel(np+ip)%e = lon(nlon) - lw

      ! Update the panel at the fold:
      panel(ip)%n    = npole

    enddo
  
    if (dodebug) then
      WRITE(6,*) "split_panels:: north split"
      CALL print_panel(panel,npnew)
    endif
 
  ! If the box intersects with the southern boundary, assuming it is a true pole,
  elseif (panel(1)%s < spole) then
    ! split the panel into two
    npnew=np*2
    do ip=1,np  !(assuming max 4)
      if (np+ip>npmax) then
        WRITE(6,*) "Max number of panels is: ", npmax
        WRITE(6,*) "You have: ", np+ip
        STOP(99) 
      endif
      panel(np+ip) = panel(ip) 
      
      ! Update the minimum latitude:
      panel(np+ip)%s = spole

      ! Update the maximum latitude:
      panel(np+ip)%n = spole - (panel(ip)%s - spole)

      ! Update the minimum longitude:
      rw = abs(ebdy - panel(ip)%e)
      panel(np+ip)%w = wbdy + rw

      ! Update the maximum longitude:
      lw = abs(panel(ip)%w - wbdy)
      panel(np+ip)%e = lon(nlon) - lw

      ! Update the original panel:
      panel(ip)%s    = spole

    enddo
    
    if (dodebug) then
      WRITE(6,*) "split_panels:: south split"
      CALL print_panel(panel,npnew)
    endif
 
  endif 
  np=npnew !np is updated with the new number of panels

  ! If the box interesects with BOTH the east and west boundaries, then just update the existing panels
  if (panel(1)%w < wbdy .and. panel(1)%e > ebdy) then 
    ! Keep the same number of panels
    npnew = np

    do ip=1,np
      panel(ip)%w = wbdy
      panel(ip)%e = ebdy
    enddo

    if (dodebug) then
      WRITE(6,*) "split_panels:: BOTH west and east split"
      CALL print_panel(panel,np)
    endif
 
  ! If the box intersects with the east/west boundary, split the panel(s) into two (again)
  elseif (panel(1)%w < wbdy) then ! Local region overlaps on western side
    ! split the panel into two
    npnew=np*2
    do ip=1,np  !(assuming max 4)
      if (np+ip>npmax) then
        WRITE(6,*) "Max number of panels is: ", npmax
        WRITE(6,*) "You have: ", np+ip
        STOP(100) 
      endif
      panel(np+ip) = panel(ip) 
      
      ! (The minimum latitude will be the same)

      ! (The maximum latitude will be the same)

      ! Update the minimum longitude:
      lw = abs(wbdy - panel(ip)%w)
      panel(np+ip)%w = ebdy - lw

      ! Update the maximum longitude:
      panel(np+ip)%e = ebdy

      ! Update the original panel:
      panel(ip)%w = wbdy

    enddo

    if (dodebug) then
      WRITE(6,*) "split_panels:: west split"
      CALL print_panel(panel,npnew)
    endif
 
  elseif (panel(1)%e > ebdy) then ! Local region overlaps on eastern side

    ! split the panel into two
    npnew=np*2
    do ip=1,np  !(assuming max 4)
      if (np+ip>npmax) then
        WRITE(6,*) "Max number of panels is: ", npmax
        WRITE(6,*) "You have: ", np+ip
        STOP(101) 
      endif
      panel(np+ip) = panel(ip) 
      
      ! (The minimum latitude will be the same)

      ! (The maximum latitude will be the same)

      ! Update the minimum longitude:
      panel(np+ip)%w = wbdy

      ! Update the maximum longitude:
      rw = abs(panel(ip)%e - ebdy)
      panel(np+ip)%e = wbdy + rw

      ! Update the original panel:
      panel(ip)%e = ebdy

    enddo

    if (dodebug) then
      WRITE(6,*) "split_panels:: east split"
      CALL print_panel(panel,npnew)
    endif
 
  endif 
  np=npnew ! set np equal to the new number of panels

  ! Do a check to make sure the panel wasn't eliminated
  if (npin > np) then
    WRITE(6,*) "ERROR:: problem in letkf_local::split_panels. npin > np, but we must have np >= npin. EXITING..."
    STOP(60)
  endif

END SUBROUTINE split_panels


SUBROUTINE print_panel(panel,np)
  TYPE (local_panel), DIMENSION(npmax), INTENT(IN) :: panel
  INTEGER, INTENT(IN) :: np
  INTEGER :: ip
 
  WRITE(6,*) "=letkf_local.f90==========================="
  do ip=1,np 
    WRITE(6,*) "-panel-type--------------------------------"
    WRITE(6,*) "ip = ", ip
    WRITE(6,*) "n = ", panel(ip)%n
    WRITE(6,*) "s = ", panel(ip)%s
    WRITE(6,*) "e = ", panel(ip)%e
    WRITE(6,*) "w = ", panel(ip)%w
    WRITE(6,*) "-------------------------------------------"
  enddo
  WRITE(6,*) "==========================================="

END SUBROUTINE print_panel


PURE SUBROUTINE find_mbb_idx(minlon,maxlon,minlat,maxlat,imin,imax,jmin,jmax)
  REAL(r_size), INTENT(IN)  :: minlon,maxlon,minlat,maxlat
  INTEGER, INTENT(OUT) :: imin,imax,jmin,jmax

  ! The local range is entirely contained within the global domain

  !-----------------------------------------------------------------------------
  ! Find jmin and jmax
  !-----------------------------------------------------------------------------
  ! Find the min and max latitudes
  do jmin=1,nlat-2
    if (DO_IRREG_GRID) then
      if (minlat < lat2d(ii,jmin+1)) exit
    else
      if (minlat < lat(jmin+1)) exit        !STEVE: (ORIGINAL)
    endif
  enddo
  do jmax=jmin,nlat-2
    if (DO_IRREG_GRID) then
      if (maxlat < lat2d(ii,jmax+1)) exit
    else
      if (maxlat < lat(jmax+1)) exit        !STEVE: (ORIGINAL)
    endif
  enddo

  !-----------------------------------------------------------------------------
  ! Find imin and imax
  !-----------------------------------------------------------------------------
  ! Find the min and max longitudes
  do imin=1,nlon-1
    if (DO_IRREG_GRID) then
      if (minlon < lon2d(imin+1,jj)) exit
    else
      if (minlon < lon(imin+1)) exit        !STEVE: (ORIGINAL)
    endif
  enddo
  do imax=imin,nlon-1
    if (DO_IRREG_GRID) then
      if (maxlon < lon2d(imax+1,jj)) exit
    else
      if (maxlon < lon(imax+1)) exit        !STEVE: (ORIGINAL)
    endif
  enddo

END SUBROUTINE find_mbb_idx


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
      if(nn >= nobs) then
!       WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
!       stop 1  !STEVE: (added)
      endif
      ! Index for observation used
      nn = nn + 1  ! Count of total obs used so far
      nobs_use(nn) = ip
    enddo
  enddo

END SUBROUTINE obs_local_sub


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
