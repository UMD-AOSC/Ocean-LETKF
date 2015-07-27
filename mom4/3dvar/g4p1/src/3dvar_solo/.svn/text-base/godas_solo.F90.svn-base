program main
!
  use constants_mod,            only: constants_init
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: mpp_error, FATAL, NOTE, mpp_pe, mpp_npes, mpp_set_current_pelist
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id, mpp_sync
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
  use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
  use time_manager_mod,         only: date_to_string

  use ocean_domains_mod,          only: get_local_indices, get_global_indices
  use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type
  use ocean_types_mod,            only: ocean_time_type

  use godas_types_mod,            only: ocean_prog_tracer_type, ocean_external_mode_type
  use godas_types_mod,            only: ocean_cor_tracer_type, ocean_rstr_tracer_type
  use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
  use godas_obs_mod,              only: godas_obsz_init, godas_obs0_init, godas_obsa_init
  use godas_mod,                  only: godas_init, godas_increment, godas_end
  use godas_rstr_mod,             only: godas_rstr_init

  use ocean_context_mod,      only: ocean_context

  implicit none

  type(ocean_grid_type), pointer   :: Grd =>NULL()
  type(ocean_domain_type), pointer :: Dom =>NULL()

  type(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
  type(ocean_cor_tracer_type), dimension(:), pointer, save :: T_cor =>NULL()
  type(ocean_obsz_type), dimension(:), pointer, save :: obs_Z =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_0 =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_A =>NULL()
  type(ocean_rstr_tracer_type), dimension(:), pointer, save :: T_rstr =>NULL()

  ! define some time types 
  type(time_type) :: Time_init    ! initial time for experiment
  type(time_type) :: Time_start   ! start time for experiment
  type(time_type) :: Time_end     ! end time for experiment (as determined by dtts)
  type(time_type) :: Run_len      ! length of experiment 
  type(time_type) :: Time_step_coupled
  type(time_type) :: Time_restart_init
  type(time_type) :: Time_restart
  type(time_type) :: Time_restart_current
  type(time_type) :: Time        

  character(len=17) :: calendar = 'julian'
  character(len=50) :: rfile = 'ocean_temp_salt.res.nc'
  character(len=20) :: name, src_temp, src_salt, fields(20)
  logical :: found_src_temp = .false., found_src_salt = .false.
  integer :: ndim, nvar, natt, ntime, i, nt

  integer :: num_prog_tracers=2
  integer :: num_cor_tracers=-1
  integer :: num_obs_z=2        ! (e.g. total of T(z), S(z))
  integer :: num_obs_0=0        ! (e.g. total of SST, SSS))
  integer :: num_obs_a=0        ! (e.g. total of Altimetry))

  integer :: nc
  integer :: calendar_type=-1

  integer :: date_init(6)=0, date(6)
  integer :: date_restart(6)
  integer :: years=0, months=0, days=0, hours=0, minutes=0, seconds=0

  integer :: isc,iec,jsc,jec
  integer :: unit, io_status, ierr

  integer :: flags=0
  integer :: nfields 
  
  character(len=256) :: version = ''
  character(len=256) :: tag = ''

  character(len=9) :: month
  character(len=64):: timestamp

  integer :: n, m 
  integer :: layout_mask(2) = (/0 , 0/)
  integer :: n_mask = 0
  integer :: mask_list(2,MAXPES)
  integer, parameter :: mp = 2*MAXPES
  data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/
  integer :: restart_interval(6) = (/0,0,0,0,0,0/)
  integer :: mpi_comm_mom
  integer ::  stdoutunit,stdlogunit

  namelist /godas_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds

  call fms_init()

  call constants_init

  flags = MPP_CLOCK_SYNC

  stdoutunit=stdout();stdlogunit=stdlog()

  ! provide for namelist over-ride
  unit = open_namelist_file('input.nml')
  read  (unit, godas_solo_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,'(/12x,a/)') '======== MODEL BEING DRIVEN BY GODAS_SOLO_MOD ========'
  write (stdoutunit, godas_solo_nml)  
  write (stdlogunit, godas_solo_nml)
  ierr = check_nml_error(io_status,'godas_solo_nml')
  call close_file (unit)

  if ( uppercase(trim(calendar)) == 'JULIAN' ) then
    calendar_type = JULIAN
  else
    call mpp_error(FATAL,&
      '==>Error from godas_solo_mod: only the Julian calendar is supported')
  endif

  ! get ocean_solo.restart : this can override settings from namelist
  if (file_exist('INPUT/ocean_solo.res')) then
      call mpp_open(unit,'INPUT/ocean_solo.res',form=MPP_ASCII,action=MPP_RDONLY)
      read(unit,*) calendar_type 
      read(unit,*) date_init
      read(unit,*) date
      call mpp_close(unit)
  endif

  call set_calendar_type (calendar_type)

  call time_interp_external_init()

  if (sum(date_init) <= 0) then
      call mpp_error(FATAL,&
      '==>Error from godas_solo_mod: date_init must be set either in ocean_solo.res or in godas_solo_nml')
  else
      Time_init  = set_date(date_init(1),date_init(2), date_init(3), &
           date_init(4),date_init(5),date_init(6))
  endif

  if (file_exist('INPUT/ocean_solo.res')) then
      Time_start =  set_date(date(1),date(2),date(3),date(4),date(5),date(6))
  else
      Time_start = Time_init
      date = date_init
  endif

  call ocean_context(Time_init, Time_start)

  call fms_io_exit

  call fms_end

  write (stdoutunit,'(/12x,a/)') 'GODAS4: --- completed ---'
  write (stdlogunit,'(/12x,a/)') 'GODAS4: --- completed ---'

end program main

