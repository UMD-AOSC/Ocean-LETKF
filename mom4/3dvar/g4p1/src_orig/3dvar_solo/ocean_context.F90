module ocean_context_mod

use fms_mod,                  only: open_namelist_file, close_file, check_nml_error
use fms_mod,                  only: FATAL, NOTE, WARNING, file_exist
use fms_io_mod,               only: set_domain, nullify_domain
use mpp_mod,              only: mpp_error, mpp_pe, mpp_npes, stdlog, stdout, mpp_sync
use mpp_io_mod,           only: mpp_open, mpp_read, mpp_close
use mpp_io_mod,           only: MPP_RDONLY, MPP_SINGLE, MPP_NETCDF, MPP_MULTI
use mpp_io_mod,           only: mpp_get_info, mpp_get_fields, mpp_get_atts, mpp_get_axes
use mpp_io_mod,           only: axistype, atttype, fieldtype

use time_manager_mod,         only: JULIAN, get_date, get_time
use time_manager_mod,         only: time_type, operator( /= ), operator( < ), operator ( / )
use time_manager_mod,         only: set_time, operator(-), operator( + ), operator( == )
use time_manager_mod,         only: operator(*)
use time_interp_external_mod, only: time_interp_external_init

use ocean_types_mod,              only: ocean_domain_type
use ocean_types_mod,              only: ocean_grid_type
use ocean_types_mod,              only: ocean_time_type

use ocean_grids_mod,              only: ocean_grids_init, set_ocean_grid_size
use ocean_grids_mod,              only: set_ocean_hgrid_arrays, set_ocean_vgrid_arrays
use ocean_domains_mod,            only: set_ocean_domain
use ocean_domains_mod,            only: get_local_indices, get_global_indices
use ocean_workspace_mod,          only: ocean_workspace_init, wrk1
use ocean_topog_mod,              only: ocean_topog_init
use ocean_util_mod,               only: ocean_util_init, write_timestamp

use ocean_parameters_mod,         only: ZSTAR
use ocean_parameters_mod,         only: DEPTH_BASED
use ocean_parameters_mod,         only: QUASI_HORIZONTAL

use godas_types_mod,            only: ocean_prog_tracer_type, ocean_external_mode_type
use godas_types_mod,            only: ocean_cor_tracer_type, ocean_rstr_tracer_type
use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
use godas_obs_mod,              only: godas_obsz_init, godas_obs0_init, godas_obsa_init
use godas_mod,                  only: godas_init, godas_increment, godas_end
use godas_rstr_mod,             only: godas_rstr_init

implicit none

private

#include <ocean_memory.h>

  ! domain layout for parallel processors. for npe=1, layout(2)=(/1,1/)
  integer :: layout(2)=(/1,1/)

  ! IO domain layout for parallel processors.
  integer :: io_layout(2)=(/0,0/)

  integer :: num_prog_tracers=2 ! (e.g., temp, salt, age)
  integer :: num_cor_tracers=2
  integer :: num_obs_z=-1        ! (e.g. total of T(z), S(z))
  integer :: num_obs_0=-1        ! (e.g. total of SST, SSS))
  integer :: num_obs_a=-1        ! (e.g. total of Altimetry))

  character(len=32) :: vertical_coordinate='zstar'
  integer :: vert_coordinate
  integer :: vert_coordinate_class
  integer :: vert_coordinate_type

  type(ocean_domain_type),        target, save   :: Domain
  type(ocean_grid_type),          target, save   :: Grid
  type(ocean_time_type),          target, save   :: Time

  type(ocean_external_mode_type), save           :: Ext_mode
  type(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
  type(ocean_cor_tracer_type), dimension(:), pointer, save :: T_cor =>NULL()
  type(ocean_obsz_type), dimension(:), pointer, save :: obs_Z =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_0 =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_A =>NULL()
  type(ocean_rstr_tracer_type), dimension(:), pointer, save :: T_rstr =>NULL()

  public ocean_context

  ! to print various cksums for debugging purposes
  logical :: debug = .false.

  logical :: module_is_initialized =.false.
  logical :: have_obc              =.false.

  character(len=30) :: rfile = 'INPUT/ocean_temp_salt.res.nc'
  type(fieldtype), allocatable, dimension(:), target :: fields
  type(fieldtype), pointer :: src_temp => NULL(), src_salt => NULL()
  character(len=50) :: fldname
  logical :: found_src_temp, found_src_salt
  integer :: unit, ndim, nvar, natt, ntime, nt
  real, allocatable, dimension(:,:,:) :: temp, salt


  namelist /ocean_context_nml/ vertical_coordinate

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_context">
!
! </DESCRIPTION>
!
subroutine ocean_context(Time_init, Time_in)
    type(time_type),       intent(in)           :: Time_init
    type(time_type),       intent(in)           :: Time_in

    integer :: i, j, k, n, taup1, lay_out(2)
    integer :: ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit
  stdoutunit=stdout();stdlogunit=stdlog()

    if (module_is_initialized) then
      call mpp_error(FATAL, &
      '==>Error in ocean_context_mod: module already initialized')
    endif

    module_is_initialized = .true.

    call time_interp_external_init()

    ! provide for namelist over-ride of defaults
    ioun = open_namelist_file()
    read  (ioun, ocean_context_nml,iostat=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_context_nml)
    write (stdlogunit, ocean_context_nml)
    ierr = check_nml_error(io_status,'ocean_context_nml')
    call close_file (ioun)

    write (stdoutunit,'(/a,i6,a/)') ' ==>Note: Running godas_solo on',mpp_npes(),'processors.'

    ! initialize ocean time type information

    Time%calendar  = JULIAN
    Time%init       = (Time_in == Time_init)
    Time%Time_init  = Time_init
    Time%model_time = Time_in
    Time%itt        = 0
    Time%taum1      = 1
    Time%tau        = 2
    Time%taup1      = 3

    write(stdoutunit,'(/a)') &
    ' ==>Note: Time%model_time = time stamp at start of this analysis:'
    call write_timestamp(Time%model_time)

    if(vertical_coordinate=='zstar') then
      vert_coordinate       = ZSTAR
      vert_coordinate_class = DEPTH_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
    else
       call mpp_error(FATAL, &
       '==>Error from ocean_context_mod: setup only for zstar vertical coordinate.')
    endif

    ! initialize grid and domain information
    call ocean_grids_init(debug, vert_coordinate, vert_coordinate_class)
    call set_ocean_grid_size(Grid, 'INPUT/grid_spec.nc')
    call set_ocean_domain(Domain, Grid, layout=layout, io_layout=io_layout)
    call set_domain(Domain%domain2d)
    call ocean_workspace_init(Domain, Grid)
    call set_ocean_hgrid_arrays(Domain, Grid)
    call ocean_topog_init(Domain, Grid, 'INPUT/grid_spec.nc', vert_coordinate_type)
    call set_ocean_vgrid_arrays(Time, Domain, Grid, have_obc)
    call ocean_util_init(Domain)

    call nullify_domain()

    write(stdoutunit,'(/12x,a/)') '======== COMPLETED OCEAN CONTEXT SETUP ========'

  ! read temperature and salinity restart

    allocate( T_prog  (num_prog_tracers) )
    T_prog(1)%name = 'temp'
    T_prog(2)%name = 'salt'

  ! set local array indices
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    ni=Grid%ni
    nj=Grid%nj
    nk=Grid%nk

    do n=1,num_prog_tracers
      allocate( T_prog(n)%field(isd:ied,jsd:jed,nk,3) )
      T_prog(n)%field(:,:,:,:)        = 0.0
    enddo

    call mpp_open(unit, trim(rfile), action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)

    found_src_temp = .false.
    found_src_salt = .false.
    do i=1,nvar
      call mpp_get_atts(fields(i),name=fldname)
      if ('temp' == trim(fldname)) then
         found_src_temp = .true.
         src_temp => fields(i)
      else if ('salt' == trim(fldname)) then
         found_src_salt = .true.
         src_salt => fields(i)
      endif
    enddo
    if(.not. found_src_temp .or. .not. found_src_salt) call mpp_error(FATAL, 'temp and/or salt not in the file '//trim(rfile) )
    allocate(temp(ni,nj,nk), salt(ni,nj,nk))
    nt = 1
    call mpp_read(unit,src_temp, temp, tindex=nt)
    call mpp_read(unit,src_salt, salt, tindex=nt)
    call mpp_close(unit)

    taup1 = Time%taup1
    do k=1,nk
      do j=jsd,jed
        do i=isd,ied
          T_prog(1)%field(i,j,k,taup1) = temp(i,j,k)
          T_prog(2)%field(i,j,k,taup1) = salt(i,j,k)
        enddo
      enddo
    enddo

    deallocate(temp, salt)

    write(stdoutunit,'(/12x,a/)') '======== TEMP & SALT RESTART READ ========'

    T_cor => godas_init(Grid, Domain, Time, T_prog(:), num_cor_tracers, debug)
    obs_Z => godas_obsz_init(Grid, Domain, Time, num_obs_z, debug)
    obs_0 => godas_obs0_init(Grid, Domain, Time, num_obs_0, debug)
    obs_A => godas_obsa_init(Grid, Domain, Time, num_obs_a, debug)
    T_rstr => godas_rstr_init(Grid, Domain, Time, T_prog(:))
    call godas_increment(Time, T_prog(:), Ext_mode, T_cor(:), obs_Z(:), obs_0(:), obs_A(:), T_rstr(:))

    call godas_end(Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)

end subroutine ocean_context

! </SUBROUTINE NAME="ocean_context">

end module ocean_context_mod
