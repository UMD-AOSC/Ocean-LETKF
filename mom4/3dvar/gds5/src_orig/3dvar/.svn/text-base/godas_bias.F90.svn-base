module godas_bias_mod
!
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module conducts a bias reduction cycle.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module conducts a bias reduction cycle.
! Initialization for godas_bias is done as well.
!</DESCRIPTION>
!
!<NAMELIST NAME="godas_bias_nml">
!
!  <DATA NAME="debug_godas_bias" TYPE="logical">
!  For debugging the godas_bias module.
!  </DATA>
!</NAMELIST>
!
use fms_mod,           only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_mod,           only: read_data, write_data
use fms_mod,           only: FATAL, WARNING, NOTE, stdout, stdlog
use mpp_mod,           only: mpp_error, mpp_pe, mpp_npes, mpp_root_pe
use mpp_mod,           only: mpp_sync, ALL_PES
use mpp_mod,           only: mpp_broadcast, mpp_transmit
use mpp_domains_mod,   only: mpp_update_domains
use mpp_domains_mod,   only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_domains_mod,   only: mpp_get_compute_domains
use mpp_io_mod,        only: mpp_open, mpp_close
use mpp_io_mod,        only: MPP_WRONLY, MPP_RDONLY, MPP_IEEE32, MPP_DIRECT, MPP_SEQUENTIAL
use mpp_io_mod,        only: MPP_SINGLE, MPP_MULTI
use time_manager_mod,  only: time_type, set_date, get_date, set_time, get_time, print_date
use time_manager_mod,  only: increment_time, decrement_time, repeat_alarm
use time_manager_mod,  only: operator(-), operator(>), operator(<), operator(<=)
use diag_manager_mod,  only: register_diag_field, send_data
use constants_mod,     only: pi

use ocean_domains_mod,          only: get_local_indices, get_global_indices
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,            only: ocean_prog_tracer_type
use ocean_types_mod,            only: ocean_external_mode_type
use ocean_types_mod,            only: ocean_time_type
use ocean_parameters_mod,       only: missing_value

use godas_types_mod,            only: ocean_cor_tracer_type
use godas_types_mod,            only: ocean_obsz_type
use godas_data_mod,             only: apply_bias_correction
use godas_data_mod,             only: id_bias_cor
use godas_data_mod,             only: num_cor_tracers
use godas_data_mod,             only: kass, kass2, ksalt, nsgobs, nsgsobs, maxits, npits, jemx
use godas_data_mod,             only: no_asm_rep, asm_bias_cnt, obs_bias_trk_cnt
use godas_data_mod,             only: gds_step, scl_incr, tovrf_b, sovrf_b
use godas_data_mod,             only: tbvrf_b, sbvrf_b, hrscl_b, vcvn_b, vsclf_b, no_lat_mx, yscl_b, ys2
use godas_data_mod,             only: xcb, xce, xcsz, ycb, yce, ycsz
use godas_data_mod,             only: s2, s1, wgns, elipt
use godas_data_mod,             only: wcn, wea, wwe, wso, wno, wgta
use godas_data_mod,             only: cvn, cvnsalt, vtmp, vsal, ev, wrkk
use godas_data_mod,             only: num_obsbiasz
use godas_data_mod,             only: temp_code, salt_code, ts_code
use godas_data_mod,             only: dtemp_max, dtemp_elm, dsalt_max, dsalt_elm
use godas_data_mod,             only: tz_wndw_fwd, tz_wndw_bwd, rtzw
use godas_data_mod,             only: sz_wndw_fwd, sz_wndw_bwd, rszw
use godas_data_mod,             only: wndw_secs
use godas_data_mod,             only: g_cg, d_cg, f_cg, e_cg, t_cg, h_cg
use godas_data_mod,             only: gds_freq, alrm_dur
use godas_data_mod,             only: assrestrt, save_all_inv, debug_godas, ovr_alrm
use godas_data_mod,             only: single_incr, apply_bias_incr, asm_ts_seq, asm_sfc_split, godas_at_end
use godas_data_mod,             only: asm_bTz, asm_bSz
use godas_data_mod,             only: spd
use godas_obs_bias_mod,         only: godas_obs_bias_track
!
implicit none

private

logical :: godas_bias_module_initialized = .false.

character(len=256) :: version = '$Id: godas_bias.F90,v 1.0 2006/11/21 08:47:00 gtn Exp $'
character(len=256) :: tagname = 'Tag $Name: gds2p0d $'
character(len=48), parameter          :: mod_name = 'godas_bias_mod'

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

integer         :: index_temp
integer         :: index_salt
integer         :: asm_code

type data_type
   character(len=3) :: gridname
   character(len=128) :: fieldname_code ! used in user's code (e.g. temp, salt, etc.)
   character(len=128) :: fieldname_file ! fieldname used in the data file (not used)
   character(len=128) :: file_name      ! name of data file
   logical :: ongrid                    ! false, not relevant, here for compatibility
   real :: factor                       ! For unit conversion, default=1
end type data_type

integer, parameter :: max_table=10

type(data_type), dimension(max_table) :: data_table

real          :: aeval, dbsq

! for diagnostics
logical :: used

! for ascii output
! integer :: unit=6

public  godas_bias_init
public  godas_bias_increment
public  godas_bias_end

! logical :: dbg = .true.

namelist /godas_bias_nml/ apply_bias_correction, tovrf_b, sovrf_b, tbvrf_b, sbvrf_b, &
                          hrscl_b, vcvn_b, vsclf_b, yscl_b

contains


! #######################################################################
! <FUNCTION NAME="godas_bias_init">
!
! <DESCRIPTION>
! Initialization code for godas_bias, returning a pointer to
! the T_bias_cor array.
! </DESCRIPTION>
!
function godas_bias_init (Grid, Domain, Time, T_prog, num_cor, debug) &
                    result (T_bias_cor)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
  integer, intent(out)                        :: num_cor

  logical, intent(in), optional               :: debug

  integer :: n

  ! return value
  type(ocean_cor_tracer_type), dimension(:), pointer :: T_bias_cor

  integer               :: pe, ierr, mtss, mtsd
  integer               :: year, month, day, hour, minute, second
  integer               :: i, j, k, ig, jg, kg
  integer               :: num_prog_tracers
  integer               :: ioun, io_status
  character(len=256)    :: record
  type(data_type)       :: default_table, data_entry
  integer               :: nu, nf, ntable, num_obs_bias_types
  real(kind=4), dimension(:,:,:), allocatable   :: buf
  real(kind=4), dimension(:,:), allocatable   :: buf2
  real(kind=4), dimension(:), allocatable   :: buf1

  character(len=48),  parameter :: sub_name = 'godas_bias_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  real, dimension(2)                    :: range_array

  if (godas_bias_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' GODAS already initialized')
  endif

  nullify(T_bias_cor)

  write( stdlog(),'(/a/)') trim(version)

  num_prog_tracers = size(T_prog)
  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  pe = mpp_pe()

! set namelist defaults (see godas_data for descriptions)

  apply_bias_correction      = .false.
  tovrf_b             = 1.0
  sovrf_b            = 1.0
  tbvrf_b             = 0.01
  sbvrf_b             = 0.01
  hrscl_b             = 3.99
  vsclf_b             = 0.5
  vcvn_b              = 1.0
  yscl_b              = 4.0

! provide for namelist over-ride

  ioun = open_namelist_file()
  read  (ioun, godas_bias_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'godas_bias_nml')
  call close_file (ioun)

! return if apply_bias_correction = .false.

  if ( apply_bias_correction ) then
    write (stdout(),'(a)') 'Bias correction is being applied.'
  else
    write (stdout(),'(a)') 'NO bias correction is being applied.'
    return
  endif

  num_cor = num_cor_tracers

! do some namelist based adjustments

  write (stdout(),'(/)')
  write (stdout(), godas_bias_nml)
  write (stdlog(), godas_bias_nml)

  asm_bias_cnt = no_asm_rep + 1
  obs_bias_trk_cnt = 0

  ! allocate T_bias_cor
  allocate( T_bias_cor  (num_cor_tracers) )
  allocate( id_bias_cor (num_cor_tracers) )

  id_bias_cor(:) = -1

  do n=1,num_cor_tracers-1
    T_bias_cor(n)%complete=.false.
  enddo
  T_bias_cor(num_cor_tracers)%complete=.true.

  ! set local array indices
  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Dom, isg, ieg, jsg, jeg)
  nk=Grd%nk
  nj=Grd%nj

  do j=1,nj
    if (Grd%grid_y_t(j) < no_lat_mx) jemx = j
  enddo
  ys2 = yscl_b * yscl_b

  do n=1,num_cor_tracers
#ifndef STATIC_MEMORY
    allocate( T_bias_cor(n)%fcor(isd:ied,jsd:jed,nk) )
#endif
    T_bias_cor(n)%fcor(:,:,:)        = 0.0
  enddo
!
! ----------------------------------------------
! read in time-constant but geo-varying error variance
! ----------------------------------------------
!
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

! read observations table
  call mpp_open(nu, 'data_table', action=MPP_RDONLY)
  ntable = 0
  asm_bTz = .false.
  asm_bSz = .false.
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'BAS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no MODEL BIAS CORRECTION entry in data_table')
  endif
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') then
      asm_bTz = .true.
    else if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') then
      asm_bSz = .true.
    endif
  enddo
  num_obs_bias_types = 0
  if (asm_bTz) num_obs_bias_types = num_obs_bias_types + 1
  if (asm_bSz) num_obs_bias_types = num_obs_bias_types + 1
  if (num_obs_bias_types .ne. num_cor_tracers) then
    call mpp_error(FATAL,'ERROR: # BAS types (T/S) in data_file does not match num_cor_tracers')
  endif
!
! ----------------------------------------------
! register diagnostics
! ----------------------------------------------
!
  do n=1,num_cor_tracers
    if (n == index_temp) then
      T_bias_cor(n)%name='tbcor'
      T_bias_cor(n)%units='Deg_C'
      T_bias_cor(n)%longname='potential temperature bias correction'
      T_bias_cor(n)%min_range=-10.0
      T_bias_cor(n)%max_range=100.0
      T_bias_cor(n)%init=.false.
      T_bias_cor(n)%file_in='INPUT/ocean_bcor.res.nc'
      T_bias_cor(n)%file_out='RESTART/ocean_bcor.res.nc'
      T_bias_cor(n)%name_in='tbcor'
    else if (n == index_salt) then
      T_bias_cor(n)%name='sbcor'
      T_bias_cor(n)%units='psu'
      T_bias_cor(n)%longname='salinity bias correction'
      T_bias_cor(n)%min_range=-10.0
      T_bias_cor(n)%max_range=100.0
      T_bias_cor(n)%init=.false.
      T_bias_cor(n)%file_in='INPUT/ocean_bcor.res.nc'
      T_bias_cor(n)%file_out='RESTART/ocean_bcor.res.nc'
      T_bias_cor(n)%name_in='sbcor'
    endif
  enddo

  ! register diagnostics  (only if godas_at_end=.false.)

  if (.not. godas_at_end) then
    do n=1,num_cor_tracers
      range_array(1) = T_bias_cor(n)%min_range
      range_array(2) = T_bias_cor(n)%max_range
      id_bias_cor(n) = register_diag_field ('ocean_model', trim(T_bias_cor(n)%name), &
           Grd%tracer_axes(1:3),                                             &
           Time%model_time, trim(T_bias_cor(n)%longname), trim(T_bias_cor(n)%units), &
           missing_value=missing_value, range=range_array)
    enddo
  endif

  if (godas_at_end) then
    do n=1,num_cor_tracers
      if (.not. T_bias_cor(n)%init) then
        write (stdout(),'(/a,a)') 'Expecting to read a GODAS restart file, ', T_bias_cor(n)%file_in
      endif

      T_bias_cor(n)%fcor(:,:,:) = 0.0
  
      if (file_exist(T_bias_cor(n)%file_in)) then
        call read_data(T_bias_cor(n)%file_in, T_bias_cor(n)%name_in, T_bias_cor(n)%fcor(:,:,:), Dom%domain2d, timelevel=1)

        if (.not. single_incr) then
          T_bias_cor(n)%fcor(:,:,:) = scl_incr * T_bias_cor(n)%fcor(:,:,:)
          write (stdout(),'(/a,1pe12.3)') 'GODAS restart increment rescaled: ',scl_incr
        endif

      else
        write (stdout(),'(/a)')'GODAS restart not found, increments set to zero.'
      endif
    enddo
  endif

  godas_bias_module_initialized = .true.

end function godas_bias_init
! </FUNCTION> NAME="godas_bias_init">


!#######################################################################
! <SUBROUTINE NAME="godas_bias_increment">
!
! <DESCRIPTION>
! Apply bias corrections from analysis.  Update analysis at specified interval.
! </DESCRIPTION>
!
subroutine godas_bias_increment (Time, T_prog, T_bias_cor, obs_bias_Z) 

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_cor_tracer_type), intent(inout)           :: T_bias_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_bias_Z(:)

  integer         :: n, taup1, pe
  integer         :: year, month, day, hour, minute, second
integer :: i, j, k

  if (.not.apply_bias_correction) then
    return
  endif

! if godas_at_end=.true. no action is taken here. instead a single analysis is called
!  in subroutine godas_bias_end

  if (.not. godas_at_end) then

    pe = mpp_pe()

    call get_date(Time%model_time, year, month, day, hour, minute, second)

    if (repeat_alarm(Time%model_time, gds_freq, alrm_dur) .or. ovr_alrm) then
      if (asm_ts_seq) then
        asm_code = temp_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
        asm_code = salt_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
      else
        asm_code = ts_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
      endif
      apply_bias_incr = .true.
      call godas_obs_bias_track(Time, obs_bias_Z)

!  scl_incr = (time step) / (time between analyses)
      if (.not. single_incr) then
        do n=1,num_cor_tracers
          T_bias_cor(n)%fcor(:,:,:) = scl_incr * T_bias_cor(n)%fcor(:,:,:)
        enddo
      endif

      if (ovr_alrm) then
        asm_bias_cnt = 2
        ovr_alrm = .false.
      else
        asm_bias_cnt = 1
      endif
    else if (asm_bias_cnt .lt. no_asm_rep) then
      if (asm_ts_seq) then
        asm_code = temp_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
        asm_code = salt_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
      else
        asm_code = ts_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z) 
      endif

!  scl_incr = (time step) / (time between analyses)
      if (.not. single_incr) then
        do n=1,num_cor_tracers
          T_bias_cor(n)%fcor(:,:,:) = scl_incr * T_bias_cor(n)%fcor(:,:,:)
        enddo
      endif

      asm_bias_cnt = asm_bias_cnt + 1
    else
      asm_bias_cnt = asm_bias_cnt + 1
    endif

! send increments/corrections to diag_manager

    do n=1,num_cor_tracers
      if (id_bias_cor(n) > 0) used = send_data (id_bias_cor(n), T_bias_cor(n)%fcor(:,:,:), &
                                  Time%model_time,rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    enddo

  endif

! apply increments 
!  increments may be applied once after each analysis cycle or applied in equal parts
!  at each time step between analyses. 
!  setting the namelist logical single_incr=.true. selects the single application and
!  and sets no_asm_rep=1.
!
  taup1   = Time%taup1
  do n=1,num_cor_tracers
    where (T_prog(index_temp)%field(:,:,:,taup1) > 0.0) &
      T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,taup1) + T_bias_cor(n)%fcor(:,:,:)
  enddo

! if using the single application zero out the increments after one use
!
  if (single_incr) then
    do n=1,num_cor_tracers
       T_bias_cor(n)%fcor(:,:,:) = 0.0
    enddo
  endif

end subroutine godas_bias_increment
! </SUBROUTINE> NAME="godas_bias_increment">


!#######################################################################
! <SUBROUTINE NAME="godas_bias_analysis">
!
! <DESCRIPTION>
! Perform an analysis.
! </DESCRIPTION>
!
subroutine godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_cor_tracer_type), intent(inout)           :: T_bias_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_bias_Z(:)

  integer         :: pe, n
  integer         :: i, ip, j, jp, k, iter
  real            :: alpha, beta, gh_old, gh_new, df_val
integer :: ii, jj

  ni    = Grd%ni
  nj    = Grd%nj
  nk    = Grd%nk
  pe = mpp_pe()

!-----------------------------------------------------------------------
!  Find the first interation of the gradient of the functional (g^1)
!  by setting the intial guess for the correction field to zero (T^1 = 0),
!  comparing the model with the observations, weighting their difference
!  with the inverse of the observation error covariance (F) and projecting
!  this onto the model grid.
!       T^1 = 0
!       g^1 = -trnsD invF To
!-----------------------------------------------------------------------
!
  t_cg = 0.0

  call init_grad_bias (Time, T_prog, obs_bias_Z)

!-----------------------------------------------------------------------
!  Do the first multiplication of the gradient by the background
!  error covariance matrix (E).
!       h^1 = E g^1
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

  call eg_lpsmthr_bias ()

!-----------------------------------------------------------------------
!  Set the initial search directions to zero.
!       d^0 = 0
!       e^0 = 0
!-----------------------------------------------------------------------

  d_cg = 0.0
  e_cg = 0.0

!-----------------------------------------------------------------------
!  Set the initial value of beta to zero
!-----------------------------------------------------------------------

  beta = 0.0

!-----------------------------------------------------------------------
!  Begin the iteration loop
!-----------------------------------------------------------------------

  do iter=1,maxits

!-----------------------------------------------------------------------
!  Update the search directions
!-----------------------------------------------------------------------

    d_cg = beta * d_cg - h_cg
    e_cg = beta * e_cg - g_cg
    do k=1,kass2
      call mpp_update_domains (d_cg(:,:,k), Dom%domain2d)
      call mpp_update_domains (e_cg(:,:,k), Dom%domain2d)
    enddo

!-----------------------------------------------------------------------
!  Compute f
!      f^n = e^n + trnsD invF D d^n
!-----------------------------------------------------------------------

    call comp_f_bias (obs_bias_Z)

!-----------------------------------------------------------------------
!  Compute the inner products <g,h  and <d,f and update alpha
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------

    gh_new = mpp_global_sum(Dom%domain2d,g_cg(:,:,:)*h_cg(:,:,:)*Grd%tmask(:,:,:),BITWISE_EXACT_SUM)
    df_val = mpp_global_sum(Dom%domain2d,d_cg(:,:,:)*f_cg(:,:,:)*Grd%tmask(:,:,:),BITWISE_EXACT_SUM)

    alpha = gh_new / df_val

!-----------------------------------------------------------------------
!  Update the field correction (T) and the gradient (g)
!      T^(n+1) = T^n + alpha d^n
!      g^(n+1) = g^n + alpha f^n
!-----------------------------------------------------------------------

    t_cg = t_cg + alpha * d_cg
    do k=1,kass2
      call mpp_update_domains (t_cg(:,:,k), Dom%domain2d)
    enddo

    if (iter .lt. maxits) then
      g_cg = g_cg + alpha * f_cg
      do k=1,kass2
        call mpp_update_domains (g_cg(:,:,k), Dom%domain2d)
      enddo

!-----------------------------------------------------------------------
!  Update h by multiplying the new gradient ( g^(n+1) ) by the
!  background error covariance E.
!       h^(n+1) = E g^(n+1)
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

      call eg_lpsmthr_bias ()

!-----------------------------------------------------------------------
!  Compute a new inner product <g,h and update beta
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------

      gh_old = gh_new
      gh_new = mpp_global_sum(Dom%domain2d,g_cg(:,:,:)*h_cg(:,:,:)*Grd%tmask(:,:,:),BITWISE_EXACT_SUM)

      beta = gh_new / gh_old

    endif
  enddo

  do n=1,num_cor_tracers
    if (n .eq. index_temp .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      T_bias_cor(n)%fcor(:,:,:) = 0.0
      do k=1,kass
        T_bias_cor(n)%fcor(:,:,k) = t_cg(:,:,k)
      enddo
    else if (n .eq. index_salt .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      T_bias_cor(n)%fcor(:,:,:) = 0.0
      do k=1,kass
        T_bias_cor(n)%fcor(:,:,k) = t_cg(:,:,k+ksalt)
      enddo
    endif
  enddo

! put increments and background error variance into obs_bias_Z

  do n=1,num_obsbiasz
    if (obs_bias_Z(n)%stat .eq. 2) then
      i = obs_bias_Z(n)%io
      ip = i + 1
      j = obs_bias_Z(n)%jo
      jp = j + 1
      if (obs_bias_Z(n)%code .eq. temp_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
        do k=1,obs_bias_Z(n)%kd
          obs_bias_Z(n)%inc(k) = obs_bias_Z(n)%a00 * T_bias_cor(index_temp)%fcor(i,j,k) + &
                            obs_bias_Z(n)%a01 * T_bias_cor(index_temp)%fcor(i,jp,k) + &
                            obs_bias_Z(n)%a11 * T_bias_cor(index_temp)%fcor(ip,jp,k) + &
                            obs_bias_Z(n)%a10 * T_bias_cor(index_temp)%fcor(ip,j,k)
          obs_bias_Z(n)%bke(k) = obs_bias_Z(n)%a00 * vtmp(i,j,k) + &
                            obs_bias_Z(n)%a01 * vtmp(i,jp,k) + &
                            obs_bias_Z(n)%a11 * vtmp(ip,jp,k) + &
                            obs_bias_Z(n)%a10 * vtmp(ip,j,k)
        enddo
        obs_bias_Z(n)%stat = 3
      else if (obs_bias_Z(n)%code .eq. salt_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
        do k=1,obs_bias_Z(n)%kd
          obs_bias_Z(n)%inc(k) = obs_bias_Z(n)%a00 * T_bias_cor(index_salt)%fcor(i,j,k) + &
                            obs_bias_Z(n)%a01 * T_bias_cor(index_salt)%fcor(i,jp,k) + &
                            obs_bias_Z(n)%a11 * T_bias_cor(index_salt)%fcor(ip,jp,k) + &
                            obs_bias_Z(n)%a10 * T_bias_cor(index_salt)%fcor(ip,j,k)
          obs_bias_Z(n)%bke(k) = obs_bias_Z(n)%a00 * vsal(i,j,k) + &
                            obs_bias_Z(n)%a01 * vsal(i,jp,k) + &
                            obs_bias_Z(n)%a11 * vsal(ip,jp,k) + &
                            obs_bias_Z(n)%a10 * vsal(ip,j,k)
        enddo
        obs_bias_Z(n)%stat = 3
      endif
    endif
  enddo

end subroutine godas_bias_analysis
! </SUBROUTINE> NAME="godas_bias_analysis"


!#######################################################################
! <SUBROUTINE NAME="godas_bias_end">
!
! <DESCRIPTION>
! Write GODAS restarts
! </DESCRIPTION>
!
subroutine godas_bias_end(Time, T_prog, T_bias_cor, obs_bias_Z)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_cor_tracer_type), intent(inout)           :: T_bias_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_bias_Z(:)
! logical :: ens_ocean

  integer :: num_prog_tracers
  integer :: i, len, m, n, taup1, pe
  character(len=128) :: filename

  if (.not.apply_bias_correction) then
    return
  endif

  pe = mpp_pe()
  write(stdout(), '(a)') 'Later a GODAS restart will be written'

  if (godas_at_end) then

! next do a single GODAS analysis.  increments will be applied after next restart.

    if (asm_ts_seq .or. num_cor_tracers .eq. 1) then
      if (asm_bTz) then
        asm_code = temp_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z)
      endif
      if (asm_bSz) then
        asm_code = salt_code
        call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z)
      endif
    else
      asm_code = ts_code
      call godas_bias_analysis (Time, T_prog, T_bias_cor, obs_bias_Z)
    endif
    call godas_obs_bias_track(Time, obs_bias_Z)

! write a restart file of GODAS increments

    do n=1,num_cor_tracers
!     call write_data(trim(T_bias_cor(n)%file_out), trim(T_bias_cor(n)%name), T_bias_cor(n)%fcor(:,:,:), domain = Dom%domain2d, append_pelist_name = ens_ocean)
      call write_data(trim(T_bias_cor(n)%file_out), trim(T_bias_cor(n)%name),&
                T_bias_cor(n)%fcor(:,:,:), domain = Dom%domain2d)
    enddo

  endif

end subroutine godas_bias_end
! </SUBROUTINE> NAME="godas_bias_end"

!!
!!#######################################################################
!! <SUBROUTINE NAME="init_grad_bias">
!!
!! <DESCRIPTION>
!! Compute the initial estimate of the gradient of the functional (g)
!! </DESCRIPTION>
!!
subroutine init_grad_bias (Time, T_prog, obs_bias_Z)
!
  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
  type(ocean_obsz_type), intent(inout)              :: obs_bias_Z(:)
!
  integer         :: i, ip, j, jp, k, kks, n, taup1, pe
  integer         :: year, month, day, hour, minute, second
  real            :: ov, aerr, aov
  type(time_type) :: diff_time, wndw_fwd, wndw_bwd
  integer         :: dsec, dday
  real            :: time_sep, time_adj
 integer :: tobc, sobc
 real :: tsum, ssum

!-----------------------------------------------------------------------
!  data types are encoded in obs%code
!      T(z)        1 <= code <= 10
!      S(z)       11 <= code <= 20
!  these are set in godas_data_mod, and can be modified if needed
!-----------------------------------------------------------------------
!
  pe = mpp_pe()

  taup1   = Time%taup1

!-----------------------------------------------------------------------
!   Set g to zero.
!-----------------------------------------------------------------------
!
  g_cg = 0.0
!
!-----------------------------------------------------------------------
!  For each observation
!   1) interpolate model forecast to observation position
!   2) compute innovation
!   3) adjust error inverse assigned to observation
!       i) adjust for time separation
!      ii) adjust for obs too far from model
!   4) multiply error inverse times inovation
!   5) project result back onto model grid
!   6) sum the gridded result in g_cg
!-----------------------------------------------------------------------
!  The profile observations
!-----------------------------------------------------------------------
!
 tobc = 0
 tsum = 0.0
 sobc = 0
 ssum = 0.0
  do n=1,num_obsbiasz
    i = obs_bias_Z(n)%io
    ip = i + 1
    j = obs_bias_Z(n)%jo
    jp = j + 1
    if (obs_bias_Z(n)%code .eq. temp_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, tz_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, tz_wndw_bwd)
      if (obs_bias_Z(n)%obs_time < wndw_fwd .and. obs_bias_Z(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_bias_Z(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rtzw)
        obs_bias_Z(n)%win = .true.
        if (save_all_inv) then
          obs_bias_Z(n)%stat = 1
        else
          if (obs_bias_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_bias_Z(n)%stat = 1
        endif
        do k=1,obs_bias_Z(n)%kd
          ov = obs_bias_Z(n)%a00 * T_prog(index_temp)%field(i,j,k,taup1) + &
               obs_bias_Z(n)%a01 * T_prog(index_temp)%field(i,jp,k,taup1) + &
               obs_bias_Z(n)%a11 * T_prog(index_temp)%field(ip,jp,k,taup1) + &
               obs_bias_Z(n)%a10 * T_prog(index_temp)%field(ip,j,k,taup1)
          ov = ov - obs_bias_Z(n)%val(k)
          if (obs_bias_Z(n)%stat .eq. 1) then
            obs_bias_Z(n)%inv(k) = ov
            if (k .eq. obs_bias_Z(n)%kd) obs_bias_Z(n)%stat = 2
          endif
          aerr = obs_bias_Z(n)%err(k)*time_adj
          aov = abs(ov)
          if (aov .lt. dtemp_max) then
            if (aov .gt. dtemp_elm) then
              aerr = aerr/(1.0+aov-dtemp_elm)**2
            endif
            g_cg(i,j,k) = g_cg(i,j,k) + ov*aerr*obs_bias_Z(n)%a00
            g_cg(i,jp,k) = g_cg(i,jp,k) + ov*aerr*obs_bias_Z(n)%a01
            g_cg(ip,jp,k) = g_cg(ip,jp,k) + ov*aerr*obs_bias_Z(n)%a11
            g_cg(ip,j,k) = g_cg(ip,j,k) + ov*aerr*obs_bias_Z(n)%a10
            obs_bias_Z(n)%aerr(k) = aerr
          else
            obs_bias_Z(n)%aerr(k) = 0.0
          endif
        enddo
      else
        time_adj = 0.0
        obs_bias_Z(n)%win = .false.
      endif
    else if (obs_bias_Z(n)%code .eq. salt_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, sz_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, sz_wndw_bwd)
      if (obs_bias_Z(n)%obs_time < wndw_fwd .and. obs_bias_Z(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_bias_Z(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rszw)
        obs_bias_Z(n)%win = .true.
        if (save_all_inv) then
          obs_bias_Z(n)%stat = 1
        else
          if (obs_bias_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_bias_Z(n)%stat = 1
        endif
        do k=1,obs_bias_Z(n)%kd
          ov = obs_bias_Z(n)%a00 * T_prog(index_salt)%field(i,j,k,taup1) + &
               obs_bias_Z(n)%a01 * T_prog(index_salt)%field(i,jp,k,taup1) + &
               obs_bias_Z(n)%a11 * T_prog(index_salt)%field(ip,jp,k,taup1) + &
               obs_bias_Z(n)%a10 * T_prog(index_salt)%field(ip,j,k,taup1)
          ov = ov - obs_bias_Z(n)%val(k)
          if (obs_bias_Z(n)%stat .eq. 1) then
            obs_bias_Z(n)%inv(k) = ov
            if (k .eq. obs_bias_Z(n)%kd) obs_bias_Z(n)%stat = 2
          endif
          aerr = obs_bias_Z(n)%err(k)*time_adj
          aov = abs(ov)
          if (aov .lt. dsalt_max) then
            if (aov .gt. dsalt_elm) then
              aerr = aerr/(1.0+aov-dsalt_elm)**2
            endif
            g_cg(i,j,k+ksalt) = g_cg(i,j,k+ksalt) + ov*aerr*obs_bias_Z(n)%a00
            g_cg(i,jp,k+ksalt) = g_cg(i,jp,k+ksalt) + ov*aerr*obs_bias_Z(n)%a01
            g_cg(ip,jp,k+ksalt) = g_cg(ip,jp,k+ksalt) + ov*aerr*obs_bias_Z(n)%a11
            g_cg(ip,j,k+ksalt) = g_cg(ip,j,k+ksalt) + ov*aerr*obs_bias_Z(n)%a10
            obs_bias_Z(n)%aerr(k) = aerr
          else
            obs_bias_Z(n)%aerr(k) = 0.0
          endif
        enddo
      else
        time_adj = 0.0
        obs_bias_Z(n)%win = .false.
      endif
    endif
  enddo
!
  do k=1,kass2
    call mpp_update_domains (g_cg(:,:,k), Dom%domain2d)
  enddo
!
  end subroutine init_grad_bias
! </SUBROUTINE> NAME="init_grad_bias"

!!
!!#######################################################################
!! <SUBROUTINE NAME="comp_f_bias">
!!
!! <DESCRIPTION>
!! Compute the 
!! </DESCRIPTION>
!!
subroutine comp_f_bias (obs_bias_Z)
!
  type(ocean_obsz_type), intent(in)        :: obs_bias_Z(:)
!
  integer :: i, ip, j, jp, k, ks1, n, pe
  real    :: ov, aerr
integer :: noot, noos
!-----------------------------------------------------------------------
!  data types are encoded in obs%code
!      T(z)        1 <= code <= 10
!      S(z)       11 <= code <= 20
!  these are set in godas_data_mod, and can be modified if needed
!-----------------------------------------------------------------------
!
  pe = mpp_pe()

!-----------------------------------------------------------------------
!   Set f to zero.
!-----------------------------------------------------------------------
!
  f_cg = 0.0
!
!-----------------------------------------------------------------------
!  For each observation
!   1) interpolate d_cg to observation position
!   2) multiply error inverse times interpolated d_cg
!   3) project result back onto model grid
!   4) sum the gridded result in f_cg
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  Profile observations
!-----------------------------------------------------------------------
!
  do n=1,num_obsbiasz
    i = obs_bias_Z(n)%io
    ip = i + 1
    j = obs_bias_Z(n)%jo
    jp = j + 1
    if (obs_bias_Z(n)%code .eq. temp_code .and. obs_bias_Z(n)%win .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      do k=1,obs_bias_Z(n)%kd
        ov = obs_bias_Z(n)%a00 * d_cg(i,j,k) + &
             obs_bias_Z(n)%a01 * d_cg(i,jp,k) + &
             obs_bias_Z(n)%a11 * d_cg(ip,jp,k) + &
             obs_bias_Z(n)%a10 * d_cg(ip,j,k)
        aerr = obs_bias_Z(n)%aerr(k)
        f_cg(i,j,k) = f_cg(i,j,k) + ov*aerr*obs_bias_Z(n)%a00
        f_cg(i,jp,k) = f_cg(i,jp,k) + ov*aerr*obs_bias_Z(n)%a01
        f_cg(ip,jp,k) = f_cg(ip,jp,k) + ov*aerr*obs_bias_Z(n)%a11
        f_cg(ip,j,k) = f_cg(ip,j,k) + ov*aerr*obs_bias_Z(n)%a10
      enddo
    else if (obs_bias_Z(n)%code .eq. salt_code .and. obs_bias_Z(n)%win .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      do k=1,obs_bias_Z(n)%kd
        ov = obs_bias_Z(n)%a00 * d_cg(i,j,k+ksalt) + &
             obs_bias_Z(n)%a01 * d_cg(i,jp,k+ksalt) + &
             obs_bias_Z(n)%a11 * d_cg(ip,jp,k+ksalt) + &
             obs_bias_Z(n)%a10 * d_cg(ip,j,k+ksalt)
        aerr = obs_bias_Z(n)%aerr(k)
        f_cg(i,j,k+ksalt) = f_cg(i,j,k+ksalt) + ov*aerr*obs_bias_Z(n)%a00
        f_cg(i,jp,k+ksalt) = f_cg(i,jp,k+ksalt) + ov*aerr*obs_bias_Z(n)%a01
        f_cg(ip,jp,k+ksalt) = f_cg(ip,jp,k+ksalt) + ov*aerr*obs_bias_Z(n)%a11
        f_cg(ip,j,k+ksalt) = f_cg(ip,j,k+ksalt) + ov*aerr*obs_bias_Z(n)%a10
      enddo
    endif
  enddo
!
!-----------------------------------------------------------------------
!   Add e to f
!-----------------------------------------------------------------------
!
  f_cg = f_cg + e_cg
!
  do k=1,kass2
    call mpp_update_domains (f_cg(:,:,k), Dom%domain2d)
  enddo
!
  end subroutine comp_f_bias
! </SUBROUTINE> NAME="comp_f_bias"


!!
!#######################################################################
! <SUBROUTINE NAME="eg_lpsmthr_bias">
!
! <DESCRIPTION>
! This subroutine multiplies g by an approximation to the first guess
! error covariance matrix [e] to get the vector h. The approximation to
! [e] is made by a series of multiplications by 1+laplacian.
! </DESCRIPTION>
!
subroutine eg_lpsmthr_bias

  integer         :: nit, n, i, j, k, kk, ka, kp, kkp
  integer         :: np, npid2
  integer         :: jbg, jfn, jgbg, jgfn
  real            :: con, col
real :: cs
integer :: pe

  ni    = Grd%ni
  nj    = Grd%nj
  nk    = Grd%nk
pe = mpp_pe()

  npid2=npits/2

  jbg = jsc
  if (jbg .eq. 1) jbg = 3
  jfn = jec
  if (jfn .ge. jemx) jfn = jemx - 1

!
!-----------------------------------------------------------------------
!   multiply g by the square root of the local vertical background
!   error covariance
!-----------------------------------------------------------------------
!
! do j=jsc,jec
!   do i=isc,iec
  do j=jsd,jed
    do i=isd,ied

      if (asm_code .eq. temp_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vtmp(i,j,k))
          enddo
        endif

! The computation of wrkk(k) alterred 6/10/08 - dwb

        do k=1,kass
          wrkk(k) = 0.0
! DBG     cs = 0.0
          do kk=1,kass
            wrkk(k) = wrkk(k) + g_cg(i,j,kk) * cvn(k,kk)
! DBG       cs = cs + cvn(k,kk)
          enddo
! DBG     wrkk(k) = ev(k) * wrkk(k) / cs
          wrkk(k) = ev(k) * wrkk(k)
        enddo

      endif

      if (asm_code .eq. salt_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vsal(i,j,k))
          enddo
        endif

! The computation of wrkk(k) alterred 6/10/08 - dwb

        do k=1,kass
          kp = k+ksalt
          wrkk(kp) = 0.0
! DBG     cs = 0.0
          do kk=1,kass
            kkp = kk+ksalt
            wrkk(kp) = wrkk(kp) + g_cg(i,j,kkp) * cvn(k,kk)
! DBG       cs = cs + cvn(k,kk)
          enddo
! DBG     wrkk(kp) = ev(k) * wrkk(kp) / cs
          wrkk(kp) = ev(k) * wrkk(kp)
        enddo

      endif

      do k=1,kass2
        g_cg(i,j,k) = wrkk(k)
      enddo

    enddo
  enddo

  do k=1,kass2
    do j=jsd,jed
      do i=isd,ied
        wcn(i,j) = 1.0 - wso(i,j) - wno(i,j) - wea(i,j) - wwe(i,j)
      enddo
    enddo

    if (k .le. kass) then
      kk = k                           ! temperature
    else
      kk = k - ksalt                   ! salinity
    endif

    do j=jsc,jec
      do i=isc,iec
        con = wso(i,j)
        if (Grd%kmt(i,j+1) .lt. kk) con = wno(i,j)
        col = con*con/((con+aeval)*con+dbsq)
        if (Grd%kmt(i,j-1) .lt. kk) wcn(i,j) = wcn(i,j)+con*col
        if (Grd%kmt(i,j+1) .lt. kk) wcn(i,j) = wcn(i,j)+con*col

        con = wwe(i,j)
        col = con*con*con/((con+aeval)*con+dbsq)
        if (Grd%kmt(i-1,j) .lt. kk) wcn(i,j) = wcn(i,j)+col
        if (Grd%kmt(i+1,j) .lt. kk) wcn(i,j) = wcn(i,j)+col
      enddo
    enddo
    call mpp_update_domains (wcn, Dom%domain2d)

    s1(:,:) = g_cg(:,:,k) * wgta(:,:)
    s2(:,:) = 0.0

    do np=1,npid2
      do j=jbg,jfn
        do i=isc,iec
          s2(i,j) = ( wcn(i,j) * s1(i,j) + wso(i,j) * s1(i,j-1) + wno(i,j) * s1(i,j+1) &
                        + wwe(i,j) * s1(i-1,j) + wea(i,j) * s1(i+1,j) ) * Grd%tmask(i,j,kk)
        enddo
      enddo
      call mpp_update_domains (s2, Dom%domain2d)
      do j=jbg,jfn
        do i=isc,iec
          s1(i,j) = ( wcn(i,j) * s2(i,j) + wno(i,j-1) * s2(i,j-1) + wso(i,j+1) * s2(i,j+1) &
                        + wwe(i,j) * s2(i-1,j) + wea(i,j) * s2(i+1,j) ) * Grd%tmask(i,j,kk)
        enddo
      enddo
      call mpp_update_domains (s1, Dom%domain2d)
    enddo

    h_cg(:,:,k) = s1(:,:) * wgta(:,:)

  enddo
!
!-----------------------------------------------------------------------
!   multiply h by the square root of the local vertical background
!   error covariance
!-----------------------------------------------------------------------
!
! do j=jsc,jec
!   do i=isc,iec
  do j=jsd,jed
    do i=isd,ied

      if (asm_code .eq. temp_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vtmp(i,j,k))
          enddo
        endif

! The computation of wrkk(k) alterred 6/10/08 - dwb

        do k=1,kass
          wrkk(k) = 0.0
! DBG     cs = 0.0
          do kk=1,kass
            wrkk(k) = wrkk(k) + h_cg(i,j,kk) * cvn(k,kk)
! DBG       cs = cs + cvn(k,kk)
          enddo
! DBG     wrkk(k) = ev(k) * wrkk(k) / cs
          wrkk(k) = ev(k) * wrkk(k)
        enddo

      endif

      if (asm_code .eq. salt_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vsal(i,j,k))
          enddo
        endif

! The computation of wrkk(k) alterred 6/10/08 - dwb

        do k=1,kass
          kp = k+ksalt
          wrkk(kp) = 0.0
! DBG     cs = 0.0
          do kk=1,kass
            kkp = kk+ksalt
            wrkk(kp) = wrkk(kp) + h_cg(i,j,kkp) * cvn(k,kk)
! DBG       cs = cs + cvn(k,kk)
          enddo
! DBG     wrkk(kp) = ev(k) * wrkk(kp) / cs
          wrkk(kp) = ev(k) * wrkk(kp)
        enddo

      endif

      do k=1,kass2
        h_cg(i,j,k) = wrkk(k)
      enddo

    enddo
  enddo

  do k=1,kass2
    call mpp_update_domains (h_cg(:,:,k), Dom%domain2d)
  enddo

end subroutine eg_lpsmthr_bias
! </SUBROUTINE> NAME="eg_lpsmthr_bias"

end module godas_bias_mod
