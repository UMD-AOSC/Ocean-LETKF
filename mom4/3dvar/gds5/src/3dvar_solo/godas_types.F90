
! The include file fms_platform.h will handle the conversion of POINTER to ALLOCATABLE arrays
! for derived type members. The conversion affects performance only and should not change
! any numeric result. It is limited to member arrays that are used within MOM4 only
! and to arrays that are never associated (=>) with another array.

! Fortran 90 requires arrays members of derived type to have the POINTER attribute.
! However, most Fortran 95 compilers now also support ALLOCATABLE array components
! (a Fortran 2003 feature). This avoids the aliasing problem afflicting pointers.
! Some compilers may require an additional switch (e.g. -fno-alias) to fully exploit
! the performance benefit of the conversion.
!
! Macros used from fms_platform.h:
! _ALLOCATABLE maps to either POINTER  or ALLOCATABLE
! _NULL        maps to either =>NULL() or "nothing"

#include <fms_platform.h>

MODULE godas_types_mod
!
!<CONTACT EMAIL="Steve.Penny@noaa.gov"> Steve Penny
!</CONTACT>
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module contains type declarations and default values for godas.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains type declarations and default values for godas.
!</DESCRIPTION>
!
  USE fms_mod,          only: write_version_number
  USE mpp_domains_mod,  only: domain2d
  USE mpp_mod,          only: FATAL, mpp_error
  USE time_manager_mod, only: time_type

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: module_is_initialized=.false.
  CHARACTER(len=128) :: version = &
     '$Id: godas_types.F90,v 1.0 2006/11/21 08:47:00 gtn Exp $'
  CHARACTER (len=128) :: tagname = &
     '$Name: gds2p0d $'

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
!################################################################################################

  TYPE, PUBLIC :: ocean_prog_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                                       ! to determine if ready to do mpp updates

     REAL, dimension(isd:ied,jsd:jed,nk,3) :: field            ! tracer field

     REAL                                  :: min_range        ! min value used for calls to diagnostic manager
     REAL                                  :: max_range        ! max value used for calls to diagnostic manager
     LOGICAL                               :: init             ! true if the input restart file is an initial condition file
     CHARACTER(len=128)                    :: file_in          ! name for input restart file
     CHARACTER(len=128)                    :: file_out         ! name for output restart file
     CHARACTER(len=32)                     :: name_in          ! name of variable to use from the input restart file
  END TYPE ocean_prog_tracer_type

  TYPE, PUBLIC :: ocean_external_mode_type
     REAL, dimension(isd:ied,jsd:jed)   :: eta_t        ! surface height on tracer cell center (m)
  END TYPE ocean_external_mode_type

  TYPE, PUBLIC :: ocean_cor_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                                       ! to determine if ready to do mpp updates

     REAL, dimension(isd:ied,jsd:jed,nk) :: fcor               ! correction for tracer field

     REAL                                  :: min_range        ! min value used for calls to diagnostic manager
     REAL                                  :: max_range        ! max value used for calls to diagnostic manager
     LOGICAL                               :: init             ! true if the input restart file is an initial condition file
     CHARACTER(len=128)                    :: file_in          ! name for input restart file
     CHARACTER(len=128)                    :: file_out         ! name for output restart file
     CHARACTER(len=32)                     :: name_in          ! name of variable to use from the input restart file
  END TYPE ocean_cor_tracer_type

  TYPE, PUBLIC :: ocean_obsz_type
     CHARACTER(len=32)  :: name                                ! e.g. temp, salt
     CHARACTER(len=32)  :: units
     CHARACTER(len=32)  :: platform                            ! e.g. XBT, TAO, Argo
     CHARACTER(len=16)  :: platid                              ! e.g. Q1901234, 52123, V60F12

     INTEGER                               :: code             ! T=1, S=11
     TYPE(time_type)                       :: obs_time
     REAL                                  :: olat             ! obs latitude
     REAL                                  :: olng             ! obs longitude
     INTEGER                               :: kd               ! number of levels
     INTEGER                               :: io               ! model grid point, sw corner
     INTEGER                               :: jo               ! model grid point, sw corner
     REAL(kind=4)                          :: a00              ! interpolation factor, sw
     REAL(kind=4)                          :: a01              ! interpolation factor, nw
     REAL(kind=4)                          :: a11              ! interpolation factor, ne
     REAL(kind=4)                          :: a10              ! interpolation factor, se
     REAL,   dimension(nk)                 :: val              ! obs value
     REAL,   dimension(nk)                 :: err              ! obs error estimate
     REAL,   dimension(nk)                 :: aerr             ! obs error adjusted in init_grad
     REAL,   dimension(nk)                 :: inv              ! innovation at minimum time offset
     REAL,   dimension(nk)                 :: inc              ! increment at minimum time offset
     REAL,   dimension(nk)                 :: bke              ! background error estimate
     INTEGER                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     LOGICAL                               :: win              ! true if within obs window
  END TYPE ocean_obsz_type

  TYPE, PUBLIC :: ocean_obs0_type
     CHARACTER(len=32)  :: name                                ! e.g. SST, Altimetry
     CHARACTER(len=32)  :: units
     CHARACTER(len=32)  :: platform                            ! e.g. SST-OI, Jason-1
     CHARACTER(len=16)  :: platid                              ! e.g. Reynolds, Aquarius, Jason2

     INTEGER                               :: code             ! SST=21, ALTM=24; do not conflict w. obsz
     TYPE(time_type)                       :: obs_time
     REAL                                  :: olat             ! obs latitude
     REAL                                  :: olng             ! obs longitude
     INTEGER                               :: kd               ! number of levels=1, ALWAYS!
     INTEGER                               :: io               ! model grid point, sw corner
     INTEGER                               :: jo               ! model grid point, sw corner
     REAL(kind=4)                          :: a00              ! interpolation factor, sw
     REAL(kind=4)                          :: a01              ! interpolation factor, nw
     REAL(kind=4)                          :: a11              ! interpolation factor, ne
     REAL(kind=4)                          :: a10              ! interpolation factor, se
     REAL(kind=4)                          :: val              ! obs value
     REAL(kind=4)                          :: err              ! obs error estimate
     REAL(kind=4)                          :: aerr             ! obs error adjusted in init_grad
     REAL(kind=4)                          :: inv              ! innovation at minimum time offset
     REAL(kind=4)                          :: inc              ! increment at minimum time offset
     REAL(kind=4)                          :: bke              ! background error estimate
     INTEGER                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     LOGICAL                               :: win              ! true if within obs window
  END TYPE ocean_obs0_type

  TYPE, PUBLIC :: ocean_rstr_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                             ! to determine if ready to do mpp updates

     REAL, dimension(isd:ied,jsd:jed) :: frstr       ! restoration for tracer field

     REAL                        :: min_range        ! min value for diagnostic manager
     REAL                        :: max_range        ! max value for to diagnostic manager
     LOGICAL                     :: init             ! true -> start is initial condition file
     CHARACTER(len=128)          :: file_in          ! name for input restart file
     CHARACTER(len=128)          :: file_out         ! name for output restart file
     CHARACTER(len=32)           :: name_in          ! name of variable in the restart file
  END TYPE ocean_rstr_tracer_type

#else
!############################################################################################
! not STATIC_MEMORY

  TYPE, PUBLIC :: ocean_prog_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                                       ! to determine if ready to do mpp updates

     REAL, _ALLOCATABLE, dimension(:,:,:,:) :: field _NULL        ! tracer field

     REAL                                  :: min_range        ! min value used for calls to diagnostic manager
     REAL                                  :: max_range        ! max value used for calls to diagnostic manager
     LOGICAL                               :: init             ! true if the input restart file is an initial condition file
     CHARACTER(len=128)                    :: file_in          ! name for input restart file
     CHARACTER(len=128)                    :: file_out         ! name for output restart file
     CHARACTER(len=32)                     :: name_in          ! name of variable to use from the input restart file
  END TYPE ocean_prog_tracer_type

  TYPE, PUBLIC :: ocean_external_mode_type
     REAL, _ALLOCATABLE, dimension(:,:) :: eta_t _NULL        ! surface height on tracer cell center (m)
  END TYPE ocean_external_mode_type

  TYPE, PUBLIC :: ocean_cor_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                                       ! to determine if ready to do mpp updates

     REAL, _ALLOCATABLE, dimension(:,:,:) :: fcor _NULL        ! correction for tracer field

     REAL(kind=4)                          :: min_range        ! min value used for calls to diagnostic manager
     REAL(kind=4)                          :: max_range        ! max value used for calls to diagnostic manager
     LOGICAL                               :: init             ! true if the input restart file is an initial condition file
     CHARACTER(len=128)                    :: file_in          ! name for input restart file
     CHARACTER(len=128)                    :: file_out         ! name for output restart file
     CHARACTER(len=32)                     :: name_in          ! name of variable to use from the input restart file
  END TYPE ocean_cor_tracer_type

  TYPE, PUBLIC :: ocean_obsz_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=32)  :: platform                            ! e.g. XBT, TAO, Argo
     CHARACTER(len=16)  :: platid                              ! e.g. Q1901234, 52123, V60F12

     INTEGER                               :: code             ! T=1, S=2
     TYPE(time_type)                       :: obs_time
     REAL                                  :: olat             ! obs latitude
     REAL                                  :: olng             ! obs longitude
     INTEGER                               :: kd               ! number of levels
     INTEGER                               :: io               ! model grid point, sw corner
     INTEGER                               :: jo               ! model grid point, sw corner
     REAL(kind=4)                          :: a00              ! interpolation factor, sw
     REAL(kind=4)                          :: a01              ! interpolation factor, nw
     REAL(kind=4)                          :: a11              ! interpolation factor, ne
     REAL(kind=4)                          :: a10              ! interpolation factor, se
     REAL, _ALLOCATABLE, dimension(:)      :: val _NULL        ! obs value
     REAL, _ALLOCATABLE, dimension(:)      :: err _NULL        ! obs error estimate
     REAL, _ALLOCATABLE, dimension(:)      :: aerr _NULL       ! obs error adjusted in init_grad
     REAL, _ALLOCATABLE, dimension(:)      :: inv _NULL        ! innovation at minimum time offset
     REAL, _ALLOCATABLE, dimension(:)      :: inc _NULL        ! increment at minimum time offset
     REAL, _ALLOCATABLE, dimension(:)      :: bke _NULL        ! background error estimate
     INTEGER                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     LOGICAL                               :: win              ! true if within obs window
  END TYPE ocean_obsz_type

  TYPE, PUBLIC :: ocean_obs0_type
     CHARACTER(len=32)  :: name                                ! e.g. SST, Altimetry
     CHARACTER(len=32)  :: units
     CHARACTER(len=32)  :: platform                            ! e.g. SST-OI, Jason-1
     CHARACTER(len=16)  :: platid                              ! e.g. Reynolds, Aquarius, Jason2

     INTEGER                               :: code             ! SST=21, ALTM=24; do not conflict w. obsz
     TYPE(time_type)                       :: obs_time
     REAL                                  :: olat             ! obs latitude
     REAL                                  :: olng             ! obs longitude
     INTEGER                               :: kd               ! number of levels=1, ALWAYS!
     INTEGER                               :: io               ! model grid point, sw corner
     INTEGER                               :: jo               ! model grid point, sw corner
     REAL(kind=4)                          :: a00              ! interpolation factor, sw
     REAL(kind=4)                          :: a01              ! interpolation factor, nw
     REAL(kind=4)                          :: a11              ! interpolation factor, ne
     REAL(kind=4)                          :: a10              ! interpolation factor, se
     REAL(kind=4)                          :: val              ! obs value
     REAL(kind=4)                          :: err              ! obs error estimate
     REAL(kind=4)                          :: aerr             ! obs error adjusted in init_grad
     REAL(kind=4)                          :: inv              ! innovation at minimum time offset
     REAL(kind=4)                          :: inc              ! increment at minimum time offset
     REAL(kind=4)                          :: bke              ! background error estimate
     INTEGER                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     LOGICAL                               :: win              ! true if within obs window
  END TYPE ocean_obs0_type

  TYPE, PUBLIC :: ocean_rstr_tracer_type
     CHARACTER(len=32)  :: name
     CHARACTER(len=32)  :: units
     CHARACTER(len=128) :: longname

     LOGICAL :: complete                    ! to determine if ready to do mpp updates

     REAL, _ALLOCATABLE, dimension(:,:) :: frstr _NULL        ! restoration for tracer field

     REAL(kind=4)       :: min_range        ! min value used for calls to diagnostic manager
     REAL(kind=4)       :: max_range        ! max value used for calls to diagnostic manager
     LOGICAL            :: init             ! true -> start is an initial condition file
     CHARACTER(len=128)    :: file_in       ! name for input restart file
     CHARACTER(len=128)    :: file_out      ! name for output restart file
     CHARACTER(len=32)     :: name_in       ! name of variable in the restart file
  END TYPE ocean_rstr_tracer_type

#endif
!############################################################################################
! end of STATIC_MEMORY

PUBLIC godas_types_init

CONTAINS

!#######################################################################
! <SUBROUTINE NAME="godas_types_init">
!
! <DESCRIPTION>
! Initialize the godas types.
! </DESCRIPTION>
!
  SUBROUTINE godas_types_init()

    if (module_is_initialized) then
       call mpp_error( FATAL, '==>Error: godas_types_init: module already initialized')
    endif
    module_is_initialized = .true.

    call write_version_number(version, tagname)

    return

  END SUBROUTINE godas_types_init
! </SUBROUTINE> NAME="godas_types_init"


END MODULE godas_types_mod

