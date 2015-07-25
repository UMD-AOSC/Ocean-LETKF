MODULE godas_rstr_mod
!
!<CONTACT EMAIL="Steve.Penny@noaa.gov"> Steve Penny
!</CONTACT>
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module manages restoration fields for GODAS.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module manages restoration fields for GODAS.
! Initialization should be called after GODAS initialization.
!</DESCRIPTION>
!
use fms_mod,                    ONLY: FATAL, WARNING, NOTE, stdout, stdlog
use fms_mod,                    ONLY: file_exist, read_data
use mpp_mod,                    ONLY: mpp_error, mpp_pe
use mpp_io_mod,                 ONLY: mpp_open, mpp_close
use mpp_io_mod,                 ONLY: MPP_WRONLY, MPP_RDONLY
use ocean_types_mod,            ONLY: ocean_grid_type, ocean_domain_type
use ocean_types_mod,            ONLY: ocean_time_type
use ocean_parameters_mod,       ONLY: missing_value
use ocean_domains_mod,          ONLY: get_local_indices
use time_interp_external_mod,   ONLY: time_interp_external, init_external_field
! use diag_manager_mod,           ONLY: register_diag_field
use time_manager_mod,           ONLY: get_date

use godas_types_mod,            ONLY: ocean_prog_tracer_type
use godas_types_mod,            ONLY: ocean_rstr_tracer_type
use godas_data_mod,             ONLY: num_rstr_tracers, id_rstr, rstr_time
use godas_data_mod,             ONLY: sst_damp, sss_damp, scl_incr
use godas_data_mod,             ONLY: single_incr, godas_at_end
!
implicit none

private

LOGICAL :: godas_rstr_module_initialized = .false.

CHARACTER(len=256) :: version = '$Id: godas_rstr.F90,v 1.0 2008/01/05 12:28:00 gtn Exp $'
CHARACTER(len=256) :: tagname = 'Tag $Name: gds1p0d $'
CHARACTER(len=48), parameter          :: mod_name = 'godas_rstr_mod'

#include <ocean_memory.h>

TYPE(ocean_grid_type), pointer   :: Grd =>NULL()
TYPE(ocean_domain_type), pointer :: Dom =>NULL()

INTEGER, allocatable, dimension(:) :: id_tie_rstr
INTEGER         :: index_temp=-1
INTEGER         :: index_salt=-1

TYPE data_type
   CHARACTER(len=3) :: gridname
   CHARACTER(len=128) :: fieldname_code ! used in user's code (e.g. temp, salt)
   CHARACTER(len=128) :: fieldname_file ! fieldname used in the data file (not used)
   CHARACTER(len=128) :: file_name      ! name of data file
   LOGICAL :: ongrid                    ! false, not relevant, here for compatibility
   real :: factor                       ! For unit conversion, default=1
END TYPE data_type

INTEGER, parameter :: max_table=10
real, allocatable, dimension(:,:) :: rdata, opn_ocn_msk

TYPE(data_type), dimension(max_table) :: data_table

! for ascii output
INTEGER :: unit=6

public  godas_rstr_init
public  godas_rstr_comp

CONTAINS

!ISSUE: change to 'godas_restore_init'
!#######################################################################
! <FUNCTION NAME="godas_rstr_init">
!
! <DESCRIPTION>
! Initialization code for restoration fields
! </DESCRIPTION>
!
FUNCTION godas_rstr_init (Grid, Domain, Time, T_prog) RESULT (T_rstr)

  TYPE(ocean_grid_type), intent(in), target   :: Grid
  TYPE(ocean_domain_type), intent(in), target :: Domain
  TYPE(ocean_time_type), intent(in)           :: Time
  TYPE(ocean_prog_tracer_type), intent(in)    :: T_prog(:)

  ! return value
  TYPE(ocean_rstr_tracer_type), dimension(:), pointer :: T_rstr

  INTEGER               :: n, nf, ntable, num_rstr, num_rstr_files
  INTEGER               :: nu, pe, num_prog_tracers
  CHARACTER(len=256)    :: record
  TYPE(data_type)       :: default_table, data_entry

  CHARACTER(len=48),  parameter :: sub_name = 'godas_rstr_init'
  CHARACTER(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) // &
                                                  '(' // trim(sub_name) // '): '
  CHARACTER(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) // &
                                                 '(' // trim(sub_name) // '): '
  CHARACTER(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) // &
                                                 '(' // trim(sub_name) // '): '

  real, dimension(2)                    :: range_array

  if (godas_rstr_module_initialized) then
    CALL mpp_error(FATAL, trim(error_header) // ' GODAS RSTR already initialized')
  endif

  pe = mpp_pe()

  nullify(T_rstr)

  if (num_rstr_tracers < 1) return

  ! allocate T_rstr
  ALLOCATE( T_rstr  (num_rstr_tracers) )

  do n=1,num_rstr_tracers-1
    T_rstr(n)%complete=.false.
  enddo
  T_rstr(num_rstr_tracers)%complete=.true.

  Grd => Grid
  Dom => Domain

  CALL get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)

  do n=1,num_rstr_tracers
#ifndef STATIC_MEMORY
    ALLOCATE( T_rstr(n)%frstr(isd:ied,jsd:jed) )
#endif
    T_rstr(n)%frstr(:,:)        = 0.0
  enddo

  num_prog_tracers = size(T_prog)
  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  ALLOCATE( rdata (isd:ied,jsd:jed) )
  ALLOCATE( opn_ocn_msk (isd:ied,jsd:jed) )
  ALLOCATE( id_rstr (num_prog_tracers) )
  ALLOCATE( id_tie_rstr (num_prog_tracers) )
  id_tie_rstr = -1

! initialize data table
  default_table%gridname = 'none'
  default_table%fieldname_code = 'none'
  default_table%fieldname_file = 'none'
  default_table%file_name = 'none'
  default_table%ongrid = .FALSE.
  default_table%factor = 1.
  do n = 1,max_table
    data_table(n) = default_table
  enddo

! read restoration table
  CALL mpp_open(nu, 'data_table', action=MPP_RDONLY)
  ntable = 0
  num_rstr_files = 0
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'ORS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 CALL mpp_error(FATAL,'error reading data_table')
9 continue
  CALL mpp_close(nu)
  if (ntable .eq. 0) then
    CALL mpp_error(FATAL,'no ORS entry in data_table')
  endif
  num_rstr_files = 0
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') num_rstr_files = num_rstr_files + 1
    if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') num_rstr_files = num_rstr_files + 1
  enddo

  num_rstr = 0
  if (num_rstr_files .gt. 0) then
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') then
        if (file_exist(trim(data_table(nf)%file_name))) then
            id_tie_rstr(index_temp) = init_external_field(trim(data_table(nf)%file_name), &
                 trim(T_prog(index_temp)%name), domain=Dom%domain2d)
            write(stdout(),*) '==>Note from GODAS: applying surface temp restoring'
            if (id_tie_rstr(index_temp) == -1) CALL mpp_error(FATAL,'==>Error in GODAS: cannot find restore field')
        endif
        num_rstr = num_rstr + 1
      else if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') then
        if (file_exist(trim(data_table(nf)%file_name))) then
            id_tie_rstr(index_salt) = init_external_field(trim(data_table(nf)%file_name), &
                 trim(T_prog(index_salt)%name), domain=Dom%domain2d)
            write(stdout(),*) '==>Note from GODAS: applying surface salt restoring'
            if (id_tie_rstr(index_salt) == -1) CALL mpp_error(FATAL,'==>Error in GODAS: cannot find restore field')
        endif
        num_rstr = num_rstr + 1
      endif
    enddo
  endif

  if (num_rstr /= num_rstr_tracers) CALL mpp_error(FATAL,'==>Error in GODAS: number of restore fields in namelist and data_table disagree')

! ----------------------------------------------
! register diagnostics
! ----------------------------------------------
!
!ISSUE: these should be able to be deleted
  do n=1,num_rstr_tracers
    if (n == index_temp) then
      T_rstr(n)%name='trstr'
      T_rstr(n)%units='Deg_C'
      T_rstr(n)%longname='potential temperature restoration'
      T_rstr(n)%min_range=-10.0
      T_rstr(n)%max_range=100.0
      T_rstr(n)%init=.false.
      T_rstr(n)%file_in='INPUT/ocean_rstr.res.nc'
      T_rstr(n)%file_out='RESTART/ocean_rstr.res.nc'
      T_rstr(n)%name_in='trstr'
    else if (n == index_salt) then
      T_rstr(n)%name='srstr'
      T_rstr(n)%units='psu'
      T_rstr(n)%longname='salinity restoration'
      T_rstr(n)%min_range=-10.0
      T_rstr(n)%max_range=100.0
      T_rstr(n)%init=.false.
      T_rstr(n)%file_in='INPUT/ocean_rstr.res.nc'
      T_rstr(n)%file_out='RESTART/ocean_rstr.res.nc'
      T_rstr(n)%name_in='srstr'
    endif
  enddo

  !STEVE: dave said this can be commented out:
! if (godas_at_end) then
!   do n=1,num_rstr_tracers
!     if (.not. T_rstr(n)%init) then
!       write (stdout(),'(/a,a)') 'Expecting to read a GODAS restart file, ', T_rstr(n)%file_in
!     endif

!     T_rstr(n)%frstr(:,:) = 0.0

!     if (file_exist(T_rstr(n)%file_in)) then
!       CALL read_data(T_rstr(n)%file_in, T_rstr(n)%name_in, T_rstr(n)%frstr(:,:), Dom%domain2d, timelevel=1)

!       if (.not. single_incr) then
!         T_rstr(n)%frstr(:,:) = scl_incr * T_rstr(n)%frstr(:,:)
!         write (stdout(),'(/a,1pe12.3)') 'GODAS restart increment rescaled: ',scl_incr
!       endif

!     else
!       write (stdout(),'(/a)')'GODAS restart not found, increments set to zero.'
!     endif
!   enddo
! endif

  godas_rstr_module_initialized = .true.

END FUNCTION godas_rstr_init
! </FUNCTION> NAME="godas_rstr_init">


! <SUBROUTINE NAME="godas_rstr_comp">
!
! <DESCRIPTION>
! Computes the surface restoration fields.
! </DESCRIPTION>
!
SUBROUTINE godas_rstr_comp (Time, T_prog, T_rstr)

  TYPE(ocean_time_type), intent(in)            :: Time
  TYPE(ocean_prog_tracer_type), intent(in)     :: T_prog(:)
  TYPE(ocean_rstr_tracer_type), intent(inout)  :: T_rstr(:)

  INTEGER         :: i, j, n, taup1
  INTEGER         :: year, month, day, hour, minute, second

  taup1   = Time%taup1

! use near-frazil condition as a proxy for where sea-ice is present

  CALL get_date(Time%model_time, year, month, day, hour, minute, second)

  if ((rstr_time(1) == hour .and. rstr_time(2) == minute .and. rstr_time(3) == second) .or. sum(rstr_time) < 0) then

    opn_ocn_msk(isd:ied,jsd:jed) = Grd%tmask(isd:ied,jsd:jed,1)
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,1) == 1.0) then
          if (T_prog(index_temp)%field(i,j,1,taup1) <= -0.0539*T_prog(index_salt)%field(i,j,1,taup1)) then
            opn_ocn_msk(i,j) = 0.0
          endif
        endif
      enddo
    enddo
  
! SST restoration
  
    if (id_tie_rstr(index_temp) > 0) then
      CALL time_interp_external(id_tie_rstr(index_temp), Time%model_time, rdata)
      do j = jsc,jec
        do i = isc,iec
          T_rstr(index_temp)%frstr(i,j) = &
              sst_damp*opn_ocn_msk(i,j)*(rdata(i,j)-T_prog(index_temp)%field(i,j,1,taup1))
        enddo
      enddo
    endif
  
! SSS restoration
  
    if (id_tie_rstr(index_salt) > 0) then
      CALL time_interp_external(id_tie_rstr(index_salt), Time%model_time, rdata)
      do j = jsc,jec
        do i = isc,iec
          T_rstr(index_salt)%frstr(i,j) = &
              sss_damp*opn_ocn_msk(i,j)*(rdata(i,j)-T_prog(index_salt)%field(i,j,1,taup1))
        enddo
      enddo
    endif
  else
    T_rstr(index_temp)%frstr(:,:) = 0.0
    T_rstr(index_salt)%frstr(:,:) = 0.0
  endif

  return

END SUBROUTINE godas_rstr_comp
! </SUBROUTINE> NAME="godas_rstr_comp"

END MODULE godas_rstr_mod
