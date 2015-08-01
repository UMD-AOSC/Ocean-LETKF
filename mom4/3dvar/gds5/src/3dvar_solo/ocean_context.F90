MODULE ocean_context_mod

USE fms_mod,                  ONLY: open_namelist_file, close_file, check_nml_error
USE fms_mod,                  ONLY: FATAL, NOTE, WARNING, file_exist
USE fms_io_mod,               ONLY: set_domain, nullify_domain
USE mpp_mod,              ONLY: mpp_error, mpp_pe, mpp_npes, stdlog, stdout, mpp_sync
USE mpp_io_mod,           ONLY: mpp_open, mpp_read, mpp_close
USE mpp_io_mod,           ONLY: MPP_RDONLY, MPP_SINGLE, MPP_NETCDF, MPP_MULTI
USE mpp_io_mod,           ONLY: mpp_get_info, mpp_get_fields, mpp_get_atts, mpp_get_axes
USE mpp_io_mod,           ONLY: axistype, atttype, fieldtype

USE time_manager_mod,         ONLY: JULIAN, get_date, get_time
USE time_manager_mod,         ONLY: time_type, operator( /= ), operator( < ), operator ( / )
USE time_manager_mod,         ONLY: set_time, operator(-), operator( + ), operator( == )
USE time_manager_mod,         ONLY: operator(*)
USE time_interp_external_mod, ONLY: time_interp_external_init

USE ocean_types_mod,              ONLY: ocean_domain_type
USE ocean_types_mod,              ONLY: ocean_grid_type
USE ocean_types_mod,              ONLY: ocean_time_type

USE ocean_grids_mod,              ONLY: ocean_grids_init, set_ocean_grid_size
USE ocean_grids_mod,              ONLY: set_ocean_hgrid_arrays, set_ocean_vgrid_arrays
USE ocean_domains_mod,            ONLY: set_ocean_domain
USE ocean_domains_mod,            ONLY: get_local_indices, get_global_indices
USE ocean_workspace_mod,          ONLY: ocean_workspace_init, wrk1
USE ocean_topog_mod,              ONLY: ocean_topog_init
USE ocean_util_mod,               ONLY: ocean_util_init, write_timestamp

USE ocean_parameters_mod,         ONLY: ZSTAR
USE ocean_parameters_mod,         ONLY: DEPTH_BASED
USE ocean_parameters_mod,         ONLY: QUASI_HORIZONTAL

USE godas_types_mod,            ONLY: ocean_prog_tracer_type, ocean_external_mode_type
USE godas_types_mod,            ONLY: ocean_cor_tracer_type, ocean_rstr_tracer_type
USE godas_types_mod,            ONLY: ocean_obsz_type, ocean_obs0_type
USE godas_obs_mod,              ONLY: godas_obsz_init, godas_obs0_init, godas_obsa_init
USE godas_mod,                  ONLY: godas_init, godas_increment, godas_end
USE godas_rstr_mod,             ONLY: godas_rstr_init

IMPLICIT NONE

PRIVATE

#include <ocean_memory.h>

  ! domain layout for parallel processors. for npe=1, layout(2)=(/1,1/)
  INTEGER :: layout(2)=(/1,1/)

  ! IO domain layout for parallel processors.
  INTEGER :: io_layout(2)=(/0,0/)

  INTEGER :: num_prog_tracers=2 ! (e.g., temp, salt, age)
  INTEGER :: num_cor_tracers=2
  INTEGER :: num_obs_z=-1        ! (e.g. total of T(z), S(z))
  INTEGER :: num_obs_0=-1        ! (e.g. total of SST, SSS))
  INTEGER :: num_obs_a=-1        ! (e.g. total of Altimetry))

  CHARACTER(len=32) :: vertical_coordinate='zstar'
  INTEGER :: vert_coordinate
  INTEGER :: vert_coordinate_class
  INTEGER :: vert_coordinate_type

  TYPE(ocean_domain_type),        target, save   :: Domain
  TYPE(ocean_grid_type),          target, save   :: Grid
  TYPE(ocean_time_type),          target, save   :: Time

  TYPE(ocean_external_mode_type), save           :: Ext_mode
  TYPE(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
  TYPE(ocean_cor_tracer_type), dimension(:), pointer, save :: T_cor =>NULL()
  TYPE(ocean_obsz_type), dimension(:), pointer, save :: obs_Z =>NULL()
  TYPE(ocean_obs0_type), dimension(:), pointer, save :: obs_0 =>NULL()
  TYPE(ocean_obs0_type), dimension(:), pointer, save :: obs_A =>NULL()
  TYPE(ocean_rstr_tracer_type), dimension(:), pointer, save :: T_rstr =>NULL()

  public ocean_context

  ! to print various cksums for debugging purposes
  LOGICAL :: debug = .false.

  LOGICAL :: module_is_initialized =.false.
  LOGICAL :: have_obc              =.false.

  CHARACTER(len=30) :: rfile = 'INPUT/ocean_temp_salt.res.nc'
  TYPE(fieldtype), allocatable, dimension(:), target :: fields
  TYPE(fieldtype), pointer :: src_temp => NULL(), src_salt => NULL()
  CHARACTER(len=50) :: fldname
  LOGICAL :: found_src_temp, found_src_salt
  INTEGER :: unit, ndim, nvar, natt, ntime, nt
  real, allocatable, dimension(:,:,:) :: temp, salt


  namelist /ocean_context_nml/ vertical_coordinate

CONTAINS


!#######################################################################
! <SUBROUTINE NAME="ocean_context">
!
! </DESCRIPTION>
!
SUBROUTINE ocean_context(Time_init, Time_in)
    TYPE(time_type),       intent(in)           :: Time_init
    TYPE(time_type),       intent(in)           :: Time_in

    INTEGER :: i, j, k, n, taup1, lay_out(2)
    INTEGER :: ioun, io_status, ierr

  INTEGER :: stdoutunit,stdlogunit
  stdoutunit=stdout();stdlogunit=stdlog()

    if (module_is_initialized) then
      CALL mpp_error(FATAL, &
      '==>Error in ocean_context_mod: module already initialized')
    endif

    module_is_initialized = .true.

    CALL time_interp_external_init()

    ! provide for namelist over-ride of defaults
    ioun = open_namelist_file()
    read  (ioun, ocean_context_nml,iostat=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_context_nml)
    write (stdlogunit, ocean_context_nml)
    ierr = check_nml_error(io_status,'ocean_context_nml')
    CALL close_file (ioun)

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
    CALL write_timestamp(Time%model_time)

    if(vertical_coordinate=='zstar') then
      vert_coordinate       = ZSTAR
      vert_coordinate_class = DEPTH_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
    else
       CALL mpp_error(FATAL, &
       '==>Error from ocean_context_mod: setup only for zstar vertical coordinate.')
    endif

    ! initialize grid and domain information
    CALL ocean_grids_init(debug, vert_coordinate, vert_coordinate_class)
    CALL set_ocean_grid_size(Grid, 'INPUT/grid_spec.nc')
    CALL set_ocean_domain(Domain, Grid, layout=layout, io_layout=io_layout)
    CALL set_domain(Domain%domain2d)
    CALL ocean_workspace_init(Domain, Grid)
    CALL set_ocean_hgrid_arrays(Domain, Grid)
    CALL ocean_topog_init(Domain, Grid, 'INPUT/grid_spec.nc', vert_coordinate_type)
    CALL set_ocean_vgrid_arrays(Time, Domain, Grid, have_obc)
    CALL ocean_util_init(Domain)

    CALL nullify_domain()

    write(stdoutunit,'(/12x,a/)') '======== COMPLETED OCEAN CONTEXT SETUP ========'

  ! read temperature and salinity restart

    allocate( T_prog  (num_prog_tracers) )
    T_prog(1)%name = 'temp'
    T_prog(2)%name = 'salt'

  ! set local array indices
    CALL get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    ni=Grid%ni
    nj=Grid%nj
    nk=Grid%nk

    do n=1,num_prog_tracers
      allocate( T_prog(n)%field(isd:ied,jsd:jed,nk,3) )
      T_prog(n)%field(:,:,:,:)        = 0.0
    enddo

    !STEVE: read the restart file
    CALL mpp_open(unit, trim(rfile), action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    !STEVE: get file info from netcdf restart file
    CALL mpp_get_info(unit, ndim, nvar, natt, ntime)

    !STEVE: get fields from netcdf restart file
    allocate(fields(nvar))
    CALL mpp_get_fields(unit, fields)

    !STEVE: get attributes from netcdf restart file
    found_src_temp = .false.
    found_src_salt = .false.
    do i=1,nvar
      CALL mpp_get_atts(fields(i),name=fldname)
      if ('temp' == trim(fldname)) then
         found_src_temp = .true.
         src_temp => fields(i)
      else if ('salt' == trim(fldname)) then
         found_src_salt = .true.
         src_salt => fields(i)
      endif
    enddo
    if(.not. found_src_temp .or. .not. found_src_salt) CALL mpp_error(FATAL, 'temp and/or salt not in the file '//trim(rfile) )

    !STEVE: read the temperature and salinity from the netcdf restart file
    allocate(temp(ni,nj,nk), salt(ni,nj,nk))
    nt = 1
    CALL mpp_read(unit,src_temp, temp, tindex=nt)
    CALL mpp_read(unit,src_salt, salt, tindex=nt)
    CALL mpp_close(unit)

    !STEVE: assign retart values to the T_prog data structure
    taup1 = Time%taup1
    do k=1,nk
      do j=jsd,jed
        do i=isd,ied
          !STEVE: editing this to fix the halo issue:
!         T_prog(1)%field(i,j,k,taup1) = temp(i,j,k)
!         T_prog(2)%field(i,j,k,taup1) = salt(i,j,k)
          T_prog(1)%field(i,j,k,taup1) = temp(MODULO(i-1,ni)+1,MODULO(j-1,nj)+1,k)
          T_prog(2)%field(i,j,k,taup1) = salt(MODULO(i-1,ni)+1,MODULO(j-1,nj)+1,k)
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

    print *, "Calling godas_increment ..."
    CALL godas_increment(Time, T_prog(:), Ext_mode, T_cor(:), obs_Z(:), obs_0(:), obs_A(:), T_rstr(:))

    print *, "Calling godas_end ..."
    CALL godas_end(Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)

END SUBROUTINE ocean_context

! </SUBROUTINE NAME="ocean_context">

END MODULE ocean_context_mod
