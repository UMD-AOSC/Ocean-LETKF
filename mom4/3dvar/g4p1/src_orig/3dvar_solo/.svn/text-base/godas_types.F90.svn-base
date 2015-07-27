
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

module godas_types_mod
!
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
  use fms_mod,          only: write_version_number
  use mpp_domains_mod,  only: domain2d
  use mpp_mod,          only: FATAL, mpp_error
  use time_manager_mod, only: time_type

  implicit none

  private

  logical :: module_is_initialized=.false.
  character(len=128) :: version = &
     '$Id: godas_types.F90,v 1.0 2006/11/21 08:47:00 gtn Exp $'
  character (len=128) :: tagname = &
     '$Name: gds2p0d $'

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
!################################################################################################

  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                                       ! to determine if ready to do mpp updates

     real, dimension(isd:ied,jsd:jed,nk,3) :: field              ! tracer field

     real                                  :: min_range        ! min value used for calls to diagnostic manager
     real                                  :: max_range        ! max value used for calls to diagnostic manager
     logical                               :: init             ! true if the input restart file is an initial condition file
     character(len=128)                    :: file_in          ! name for input restart file
     character(len=128)                    :: file_out         ! name for output restart file
     character(len=32)                     :: name_in          ! name of variable to use from the input restart file
  end type ocean_prog_tracer_type

  type, public :: ocean_external_mode_type
     real, dimension(isd:ied,jsd:jed)   :: eta_t        ! surface height on tracer cell center (m)
  end type ocean_external_mode_type

  type, public :: ocean_cor_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                                       ! to determine if ready to do mpp updates

     real, dimension(isd:ied,jsd:jed,nk) :: fcor               | correction for tracer field

     real                                  :: min_range        ! min value used for calls to diagnostic manager
     real                                  :: max_range        ! max value used for calls to diagnostic manager
     logical                               :: init             ! true if the input restart file is an initial condition file
     character(len=128)                    :: file_in          ! name for input restart file
     character(len=128)                    :: file_out         ! name for output restart file
     character(len=32)                     :: name_in          ! name of variable to use from the input restart file
  end type ocean_cor_tracer_type

  type, public :: ocean_obsz_type
     character(len=32)  :: name                                ! e.g. temp, salt
     character(len=32)  :: units
     character(len=32)  :: platform                            ! e.g. XBT, TAO, Argo
     character(len=16)  :: platid                              ! e.g. Q1901234, 52123, V60F12

     integer                               :: code             ! T=1, S=11
     type(time_type)                       :: obs_time
     real                                  :: olat             ! obs latitude
     real                                  :: olng             ! obs longitude
     integer                               :: kd               ! number of levels
     integer                               :: io               ! model grid point, sw corner
     integer                               :: jo               ! model grid point, sw corner
     real(kind=4)                          :: a00              ! interpolation factor, sw
     real(kind=4)                          :: a01              ! interpolation factor, nw
     real(kind=4)                          :: a11              ! interpolation factor, ne
     real(kind=4)                          :: a10              ! interpolation factor, se
     real,   dimension(nk)                 :: val              ! obs value
     real,   dimension(nk)                 :: err              ! obs error estimate
     real,   dimension(nk)                 :: aerr             ! obs error adjusted in init_grad
     real,   dimension(nk)                 :: inv              ! innovation at minimum time offset
     real,   dimension(nk)                 :: inc              ! increment at minimum time offset
     real,   dimension(nk)                 :: bke              ! background error estimate
     integer                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     logical                               :: win              ! true if within obs window
  end type ocean_obsz_type

  type, public :: ocean_obs0_type
     character(len=32)  :: name                                ! e.g. SST, Altimetry
     character(len=32)  :: units
     character(len=32)  :: platform                            ! e.g. SST-OI, Jason-1
     character(len=16)  :: platid                              ! e.g. Reynolds, Aquarius, Jason2

     integer                               :: code             ! SST=21, ALTM=24; do not conflict w. obsz
     type(time_type)                       :: obs_time
     real                                  :: olat             ! obs latitude
     real                                  :: olng             ! obs longitude
     integer                               :: kd               ! number of levels=1, ALWAYS!
     integer                               :: io               ! model grid point, sw corner
     integer                               :: jo               ! model grid point, sw corner
     real(kind=4)                          :: a00              ! interpolation factor, sw
     real(kind=4)                          :: a01              ! interpolation factor, nw
     real(kind=4)                          :: a11              ! interpolation factor, ne
     real(kind=4)                          :: a10              ! interpolation factor, se
     real(kind=4)                          :: val              ! obs value
     real(kind=4)                          :: err              ! obs error estimate
     real(kind=4)                          :: aerr             ! obs error adjusted in init_grad
     real(kind=4)                          :: inv              ! innovation at minimum time offset
     real(kind=4)                          :: inc              ! increment at minimum time offset
     real(kind=4)                          :: bke              ! background error estimate
     integer                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     logical                               :: win              ! true if within obs window
  end type ocean_obs0_type

  type, public :: ocean_rstr_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                             ! to determine if ready to do mpp updates

     real, dimension(isd:ied,jsd:jed) :: frstr       ! restoration for tracer field

     real                        :: min_range        ! min value for diagnostic manager
     real                        :: max_range        ! max value for to diagnostic manager
     logical                     :: init             ! true -> start is initial condition file
     character(len=128)          :: file_in          ! name for input restart file
     character(len=128)          :: file_out         ! name for output restart file
     character(len=32)           :: name_in          ! name of variable in the restart file
  end type ocean_rstr_tracer_type

#else
!############################################################################################
! not STATIC_MEMORY

  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                                       ! to determine if ready to do mpp updates

     real, _ALLOCATABLE, dimension(:,:,:,:) :: field _NULL        ! tracer field

     real                                  :: min_range        ! min value used for calls to diagnostic manager
     real                                  :: max_range        ! max value used for calls to diagnostic manager
     logical                               :: init             ! true if the input restart file is an initial condition file
     character(len=128)                    :: file_in          ! name for input restart file
     character(len=128)                    :: file_out         ! name for output restart file
     character(len=32)                     :: name_in          ! name of variable to use from the input restart file
  end type ocean_prog_tracer_type

  type, public :: ocean_external_mode_type
     real, _ALLOCATABLE, dimension(:,:) :: eta_t _NULL        ! surface height on tracer cell center (m)
  end type ocean_external_mode_type

  type, public :: ocean_cor_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                                       ! to determine if ready to do mpp updates

     real, _ALLOCATABLE, dimension(:,:,:) :: fcor _NULL        ! correction for tracer field

     real(kind=4)                          :: min_range        ! min value used for calls to diagnostic manager
     real(kind=4)                          :: max_range        ! max value used for calls to diagnostic manager
     logical                               :: init             ! true if the input restart file is an initial condition file
     character(len=128)                    :: file_in          ! name for input restart file
     character(len=128)                    :: file_out         ! name for output restart file
     character(len=32)                     :: name_in          ! name of variable to use from the input restart file
  end type ocean_cor_tracer_type

  type, public :: ocean_obsz_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=32)  :: platform                            ! e.g. XBT, TAO, Argo
     character(len=16)  :: platid                              ! e.g. Q1901234, 52123, V60F12

     integer                               :: code             ! T=1, S=2
     type(time_type)                       :: obs_time
     real                                  :: olat             ! obs latitude
     real                                  :: olng             ! obs longitude
     integer                               :: kd               ! number of levels
     integer                               :: io               ! model grid point, sw corner
     integer                               :: jo               ! model grid point, sw corner
     real(kind=4)                          :: a00              ! interpolation factor, sw
     real(kind=4)                          :: a01              ! interpolation factor, nw
     real(kind=4)                          :: a11              ! interpolation factor, ne
     real(kind=4)                          :: a10              ! interpolation factor, se
     real, _ALLOCATABLE, dimension(:)      :: val _NULL        ! obs value
     real, _ALLOCATABLE, dimension(:)      :: err _NULL        ! obs error estimate
     real, _ALLOCATABLE, dimension(:)      :: aerr _NULL       ! obs error adjusted in init_grad
     real, _ALLOCATABLE, dimension(:)      :: inv _NULL        ! innovation at minimum time offset
     real, _ALLOCATABLE, dimension(:)      :: inc _NULL        ! increment at minimum time offset
     real, _ALLOCATABLE, dimension(:)      :: bke _NULL        ! background error estimate
     integer                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     logical                               :: win              ! true if within obs window
  end type ocean_obsz_type

  type, public :: ocean_obs0_type
     character(len=32)  :: name                                ! e.g. SST, Altimetry
     character(len=32)  :: units
     character(len=32)  :: platform                            ! e.g. SST-OI, Jason-1
     character(len=16)  :: platid                              ! e.g. Reynolds, Aquarius, Jason2

     integer                               :: code             ! SST=21, ALTM=24; do not conflict w. obsz
     type(time_type)                       :: obs_time
     real                                  :: olat             ! obs latitude
     real                                  :: olng             ! obs longitude
     integer                               :: kd               ! number of levels=1, ALWAYS!
     integer                               :: io               ! model grid point, sw corner
     integer                               :: jo               ! model grid point, sw corner
     real(kind=4)                          :: a00              ! interpolation factor, sw
     real(kind=4)                          :: a01              ! interpolation factor, nw
     real(kind=4)                          :: a11              ! interpolation factor, ne
     real(kind=4)                          :: a10              ! interpolation factor, se
     real(kind=4)                          :: val              ! obs value
     real(kind=4)                          :: err              ! obs error estimate
     real(kind=4)                          :: aerr             ! obs error adjusted in init_grad
     real(kind=4)                          :: inv              ! innovation at minimum time offset
     real(kind=4)                          :: inc              ! increment at minimum time offset
     real(kind=4)                          :: bke              ! background error estimate
     integer                               :: stat             ! if (1,2,3) fill (inv,inc & bke,outfile)

     logical                               :: win              ! true if within obs window
  end type ocean_obs0_type

  type, public :: ocean_rstr_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: complete                    ! to determine if ready to do mpp updates

     real, _ALLOCATABLE, dimension(:,:) :: frstr _NULL        ! restoration for tracer field

     real(kind=4)       :: min_range        ! min value used for calls to diagnostic manager
     real(kind=4)       :: max_range        ! max value used for calls to diagnostic manager
     logical            :: init             ! true -> start is an initial condition file
     character(len=128)    :: file_in       ! name for input restart file
     character(len=128)    :: file_out      ! name for output restart file
     character(len=32)     :: name_in       ! name of variable in the restart file
  end type ocean_rstr_tracer_type

#endif
!############################################################################################
! end of STATIC_MEMORY

public godas_types_init

contains

!#######################################################################
! <SUBROUTINE NAME="godas_types_init">
!
! <DESCRIPTION>
! Initialize the godas types.
! </DESCRIPTION>
!
  subroutine godas_types_init()

    if (module_is_initialized) then
       call mpp_error( FATAL, '==>Error: godas_types_init: module already initialized')
    endif
    module_is_initialized = .true.

    call write_version_number(version, tagname)

    return

  end subroutine godas_types_init
! </SUBROUTINE> NAME="godas_types_init"


end module godas_types_mod

