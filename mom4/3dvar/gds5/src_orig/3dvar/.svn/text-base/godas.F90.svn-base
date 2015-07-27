module godas_mod
!
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module conducts an assimilation cycle.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module conducts an assimilation cycle.
! Initialization for godas is done as well.
!</DESCRIPTION>
!
!<NAMELIST NAME="godas_nml">
!
!  <DATA NAME="debug_godas" TYPE="logical">
!  For debugging the godas module.
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
use godas_types_mod,            only: ocean_cor_tracer_type, ocean_rstr_tracer_type
use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
use godas_data_mod,             only: id_cor, id_rstr
use godas_data_mod,             only: num_cor_tracers, num_rstr_tracers, rstr_time
use godas_data_mod,             only: sst_damp, sss_damp
use godas_data_mod,             only: kass, kass2, ksalt, nsgobs, nsgsobs, maxits, npits, jemx
use godas_data_mod,             only: no_asm_rep, asm_cnt, obs_trk_cnt
use godas_data_mod,             only: gds_step, scl_incr, tov0f, sov0f, tbv0f, sbv0f, tovrf, sovrf
use godas_data_mod,             only: tbvrf, sbvrf, hrscl, hrscl0, vcvn, vsclf, no_lat_mx, yscl, ys2
use godas_data_mod,             only: xcb, xce, xcsz, ycb, yce, ycsz
use godas_data_mod,             only: s2, s1, wgns, elipt
use godas_data_mod,             only: wcn, wea, wwe, wso, wno, wgta
use godas_data_mod,             only: wcn_s, wea_s, wwe_s, wso_s, wno_s, wgta_s
use godas_data_mod,             only: cvn, cvnsalt, vtmp, vsal, ev, wrkk, vtmp_s, vsal_s
use godas_data_mod,             only: eta_clm, cdnz, cdnzs
use godas_data_mod,             only: num_obsz, num_obs0, num_obsa
use godas_data_mod,             only: temp_code, salt_code, sst_code, sss_code, altm_code, ts_code
use godas_data_mod,             only: dtemp_max, dtemp_elm, dsalt_max, dsalt_elm
use godas_data_mod,             only: dsst_max, dsst_elm, dsss_max, dsss_elm, daltm_max, daltm_elm
use godas_data_mod,             only: tz_wndw_fwd, tz_wndw_bwd, rtzw
use godas_data_mod,             only: sz_wndw_fwd, sz_wndw_bwd, rszw
use godas_data_mod,             only: t0_wndw_fwd, t0_wndw_bwd, rt0w
use godas_data_mod,             only: s0_wndw_fwd, s0_wndw_bwd, rs0w
use godas_data_mod,             only: al_wndw_fwd, al_wndw_bwd, ralw, wndw_secs
use godas_data_mod,             only: g_cg, d_cg, f_cg, e_cg, t_cg, h_cg
use godas_data_mod,             only: g_cg_s, d_cg_s, f_cg_s, e_cg_s, t_cg_s, h_cg_s
use godas_data_mod,             only: gds_freq, alrm_dur
use godas_data_mod,             only: assrestrt, rstrestrt, restore_sfc, save_all_inv, debug_godas, ovr_alrm
use godas_data_mod,             only: single_incr, apply_incr, asm_ts_seq, asm_sfc_split, godas_at_end
use godas_data_mod,             only: asm_Tz, asm_Sz, asm_T0, asm_S0, asm_Al
use godas_data_mod,             only: spd

use godas_obs_mod,              only: godas_obs_track
use godas_rstr_mod,             only: godas_rstr_comp
!
implicit none

private

logical :: godas_module_initialized = .false.

character(len=256) :: version = '$Id: godas.F90,v 1.0 2006/11/21 08:47:00 gtn Exp $'
character(len=256) :: tagname = 'Tag $Name: gds2p0d $'
character(len=48), parameter          :: mod_name = 'godas_mod'

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

integer         :: index_temp
integer         :: index_salt
integer         :: asm_code

type data_type
   character(len=3) :: gridname
   character(len=128) :: fieldname_code ! used in user's code (e.g. mdl_tvv, mdl_svv, etc.)
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

public  godas_init
public  godas_increment
public  godas_end

! logical :: dbg = .true.

namelist /godas_nml/ num_cor_tracers, kass, nsgobs, nsgsobs, maxits, npits, &
                     gds_step, no_asm_rep, single_incr, tov0f, sov0f, tbv0f, &
                     sbv0f, tovrf, sovrf, tbvrf, sbvrf, hrscl, hrscl0, vcvn, vsclf, &
                     no_lat_mx, yscl, tz_wndw_fwd, tz_wndw_bwd, sz_wndw_fwd, &
                     sz_wndw_bwd, t0_wndw_fwd, t0_wndw_bwd, s0_wndw_fwd, &
                     s0_wndw_bwd, al_wndw_fwd, al_wndw_bwd, assrestrt, &
                     save_all_inv, asm_ts_seq, asm_sfc_split, godas_at_end, &
                     restore_sfc, num_rstr_tracers, rstr_time, sst_damp, sss_damp, &
                     rstrestrt, debug_godas

contains


!#######################################################################
! <FUNCTION NAME="godas_init">
!
! <DESCRIPTION>
! Initialization code for godas, returning a pointer to
! the T_cor array.
! </DESCRIPTION>
!
function godas_init (Grid, Domain, Time, T_prog, num_cor, debug) &
                    result (T_cor)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
  integer, intent(out)                        :: num_cor

  logical, intent(in), optional               :: debug

  integer :: n

  ! return value
  type(ocean_cor_tracer_type), dimension(:), pointer :: T_cor

  integer               :: pe, ierr, mtss, mtsd
  integer               :: year, month, day, hour, minute, second
  integer               :: i, j, k, ig, jg, kg
  integer               :: num_prog_tracers
  integer               :: ioun, io_status
  character(len=256)    :: record
  type(data_type)       :: default_table, data_entry
  integer               :: nu, nf, ntable, num_obs_types, num_var_files, num_vr0_files
  integer               :: num_alt_files
  real(kind=4), dimension(:,:,:), allocatable   :: buf
  real(kind=4), dimension(:,:), allocatable   :: buf2
  real(kind=4), dimension(:), allocatable   :: buf1
  logical               :: Tzv, Szv, T0v, S0v, Aave, Acdnz

  character(len=48),  parameter :: sub_name = 'godas_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  real, dimension(2)                    :: range_array

  if (godas_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' GODAS already initialized')
  endif

  nullify(T_cor)

  write( stdlog(),'(/a/)') trim(version)

  num_prog_tracers = size(T_prog)
  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  pe = mpp_pe()

! set namelist defaults (see godas_data for descriptions)

  num_cor_tracers     = 2
  kass                = 30         ! standard for 40-level MOM; use 35 for deep
  nsgobs              = 1          ! standard
  nsgsobs             = 1          ! standard
  maxits              = 3
  npits               = 200
  gds_step            = 43200
  no_asm_rep          = 1
  single_incr         = .false.
  tov0f               = 1.0
  sov0f               = 0.1
  tbv0f               = 1.0
  sbv0f               = 0.1
  tovrf               = 1.0
  sovrf               = 1.0
  tbvrf               = 0.01
  sbvrf               = 0.01
  hrscl               = 3.99
  hrscl0              = 3.99
  vsclf               = 0.5
  vcvn                = 1.0
  no_lat_mx           = 63.0
  yscl                = 4.0
  tz_wndw_fwd         = 14
  tz_wndw_bwd         = 14
  sz_wndw_fwd         = 14
  sz_wndw_bwd         = 14
  t0_wndw_fwd         = 7
  t0_wndw_bwd         = 7
  s0_wndw_fwd         = 7
  s0_wndw_bwd         = 7
  al_wndw_fwd         = 14
  al_wndw_bwd         = 14
  assrestrt           = .true.
  restore_sfc         = .false.
  num_rstr_tracers    = 2
  rstr_time(:)        = -1
  rstrestrt           = .true.
  sst_damp            = 0.1
  sss_damp            = 0.1
  save_all_inv        = .false.
  asm_ts_seq          = .false.
  asm_sfc_split       = .false.
  godas_at_end        = .false.
  debug_godas         = .false.

! provide for namelist over-ride

  ioun = open_namelist_file()
  read  (ioun, godas_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'godas_nml')
  call close_file (ioun)

  num_cor = num_cor_tracers

! do some namelist based adjustments

  if (asm_ts_seq) then
    ksalt = 0
  else
    ksalt = kass
  endif
  kass2 = kass + ksalt

  if (.not. restore_sfc) num_rstr_tracers = 0
  if (restore_sfc) asm_sfc_split = .true.

  if (.not. asm_sfc_split) then
    tov0f = tovrf
    tbv0f = tbvrf
    sov0f = sovrf
    sbv0f = sbvrf
  endif

  if (single_incr) no_asm_rep = 1

  if (godas_at_end) then
    wndw_secs = 43200
    rtzw = 0.0
    rszw = 0.0
    rt0w = 0.0
    rs0w = 0.0
    ralw = 0.0
  else
    wndw_secs = 0
    if (tz_wndw_fwd > tz_wndw_bwd) then
      rtzw = 1.0 / real(tz_wndw_fwd)
    else
      rtzw = 1.0 / real(tz_wndw_bwd)
    endif
    if (sz_wndw_fwd > sz_wndw_bwd) then
      rszw = 1.0 / real(sz_wndw_fwd)
    else
      rszw = 1.0 / real(sz_wndw_bwd)
    endif
    if (t0_wndw_fwd > t0_wndw_bwd) then
      rt0w = 1.0 / real(t0_wndw_fwd)
    else
      rt0w = 1.0 / real(t0_wndw_bwd)
    endif
    if (s0_wndw_fwd > s0_wndw_bwd) then
      rs0w = 1.0 / real(s0_wndw_fwd)
    else
      rs0w = 1.0 / real(s0_wndw_bwd)
    endif
    if (al_wndw_fwd > al_wndw_bwd) then
      ralw = 1.0 / real(al_wndw_fwd)
    else
      ralw = 1.0 / real(al_wndw_bwd)
    endif
  endif

  if (sst_damp .lt. 0.0) sst_damp = 0.0
  if (sst_damp .gt. 1.0) sst_damp = 1.0
  if (sss_damp .lt. 0.0) sss_damp = 0.0
  if (sss_damp .gt. 1.0) sss_damp = 1.0

  write (stdout(),'(/)')
  write (stdout(), godas_nml)
  write (stdlog(), godas_nml)

  ovr_alrm = assrestrt
  asm_cnt = no_asm_rep + 1
  obs_trk_cnt = 0

  if (PRESENT(debug) .and. .not. debug_godas) then
    debug_godas = debug
  endif

  gds_freq = set_time(gds_step, 0)
  call get_time(Time%Time_step, mtss, mtsd)
  alrm_dur = set_time(mtss, 0)

! for now "godas_at_end" depends on gds_step = length_of_run in seconds

  if (godas_at_end) then
    scl_incr = float(mtss) / float(gds_step)
  else
    scl_incr = float(mtss) / float(gds_step)
  endif

  ! allocate T_cor
  allocate( T_cor  (num_cor_tracers) )
  allocate( id_cor (num_cor_tracers) )

  id_cor(:) = -1

  do n=1,num_cor_tracers-1
    T_cor(n)%complete=.false.
  enddo
  T_cor(num_cor_tracers)%complete=.true.

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
  ys2 = yscl * yscl

  do n=1,num_cor_tracers
#ifndef STATIC_MEMORY
    allocate( T_cor(n)%fcor(isd:ied,jsd:jed,nk) )
#endif
    T_cor(n)%fcor(:,:,:)        = 0.0
  enddo
!
! ----------------------------------------------
! allocate various computational arrays
! ----------------------------------------------
!
  allocate (d_cg(isd:ied,jsd:jed,kass2))
  allocate (f_cg(isd:ied,jsd:jed,kass2))
  allocate (g_cg(isd:ied,jsd:jed,kass2))
  allocate (e_cg(isd:ied,jsd:jed,kass2))
  allocate (t_cg(isd:ied,jsd:jed,kass2))
  allocate (h_cg(isd:ied,jsd:jed,kass2))
!
  allocate (d_cg_s(isd:ied,jsd:jed))
  allocate (f_cg_s(isd:ied,jsd:jed))
  allocate (g_cg_s(isd:ied,jsd:jed))
  allocate (e_cg_s(isd:ied,jsd:jed))
  allocate (t_cg_s(isd:ied,jsd:jed))
  allocate (h_cg_s(isd:ied,jsd:jed))
!
! ----------------------------------------------
! begin calculating weights for horizontal smoother
! ----------------------------------------------
!
  call init_wghts()
!
! ----------------------------------------------
! compute vertical covariance matrix
!-----------------------------------------------
!
  allocate (cvn(kass,kass))
  allocate (cvnsalt(kass,kass))
!
  call vertical_covariance()
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
  asm_Tz = .false.
  asm_Sz = .false.
  asm_T0 = .false.
  asm_S0 = .false.
  asm_Al = .false.
  Tzv = .false.
  Szv = .false.
  T0v = .false.
  S0v = .false.
  Aave = .false.
  Acdnz = .false.
  num_var_files = 0
  num_vr0_files = 0
  num_alt_files = 0
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'OBS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no MODEL VARIANCE entry in data_table')
  endif
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') then
      asm_Tz = .true.
    else if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') then
      asm_Sz = .true.
    else if (data_table(nf)%fieldname_code(1:3) .eq. 'sst') then
      asm_T0 = .true.
    else if (data_table(nf)%fieldname_code(1:3) .eq. 'sss') then
      asm_S0 = .true.
    else if (data_table(nf)%fieldname_code(1:4) .eq. 'altm') then
      asm_Al = .true.
    else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_tvv') then
      Tzv = .true.
      num_var_files = num_var_files + 1
    else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_svv') then
      Szv = .true.
      num_var_files = num_var_files + 1
    else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_tv0') then
      T0v = .true.
      num_vr0_files = num_vr0_files + 1
    else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_sv0') then
      S0v = .true.
      num_vr0_files = num_vr0_files + 1
    else if (data_table(nf)%fieldname_code(1:4) .eq. 'etac') then
      Aave = .true.
      num_alt_files = num_alt_files + 1
    else if (data_table(nf)%fieldname_code(1:4) .eq. 'cdnz') then
      Acdnz = .true.
      num_alt_files = num_alt_files + 1
    endif
  enddo
  num_obs_types = 0
  if (asm_Tz .or. asm_T0) num_obs_types = num_obs_types + 1
  if (asm_Sz .or. asm_S0) num_obs_types = num_obs_types + 1
  if (asm_Al) num_obs_types = 2
  if (asm_Tz .and. .not.Tzv) then
    call mpp_error(FATAL,'ERROR: in data_table found T(z) OBS file, but no variance file')
  else if (Tzv .and. .not.asm_Tz) then
    call mpp_error(FATAL,'ERROR: in data_table found T(z) variance file, but no OBS file')
  else if (asm_Sz .and. .not.Szv) then
    call mpp_error(FATAL,'ERROR: in data_table found S(z) OBS file, but no variance file')
  else if (Szv .and. .not.asm_Sz) then
    call mpp_error(FATAL,'ERROR: in data_table found S(z) variance file, but no OBS file')
  else if (asm_T0 .and. .not.T0v) then
    call mpp_error(FATAL,'ERROR: in data_table found SST OBS file, but no variance file')
  else if (T0v .and. .not.asm_T0) then
    call mpp_error(FATAL,'ERROR: in data_table found SST variance file, but no OBS file')
  else if (asm_S0 .and. .not.S0v) then
    call mpp_error(FATAL,'ERROR: in data_table found SSS OBS file, but no variance file')
  else if (S0v .and. .not.asm_S0) then
    call mpp_error(FATAL,'ERROR: in data_table found SSS variance file, but no OBS file')
  else if (asm_Al .and. .not.Aave) then
    call mpp_error(FATAL,'ERROR: in data_table found Alt OBS file, but no Alt average file')
  else if (Aave .and. .not.asm_Al) then
    call mpp_error(FATAL,'ERROR: in data_table found Alt average file, but no Alt OBS file')
  else if (asm_Al .and. .not.Acdnz) then
    call mpp_error(FATAL,'ERROR: in data_table found Alt OBS file, but no Alt lin-den coef file')
  else if (Acdnz .and. .not.asm_Al) then
    call mpp_error(FATAL,'ERROR: in data_table found Alt lin-den coef file, but no Alt OBS file')
  else if (num_obs_types .ne. num_cor_tracers) then
    call mpp_error(FATAL,'ERROR: # OBS types (T/S) in data_file does not match num_cor_tracers in namelist')
  endif
!
  if (num_var_files .ne. 0) then
    allocate (vtmp(isd:ied,jsd:jed,kass))
    allocate (vsal(isd:ied,jsd:jed,kass))
    allocate (buf(isg:ieg,jsg:jeg,kass))
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_tvv') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        read (nu) year, month, day
        read (nu) ig, jg, kg
        read (nu) buf
        close(nu)
!
        do k=1,kass
          do j=jsd,jed
            do i=isd,ied
              vtmp(i,j,k) = tbvrf * buf(i,j,k)
!             vtmp(i,j,k) = buf(i,j,k)
            enddo
          enddo
        enddo
!
      else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_svv') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        read (nu) year, month, day
        read (nu) ig, jg, kg
        read (nu) buf
        close(nu)
!
        do k=1,kass
          do j=jsd,jed
            do i=isd,ied
              vsal(i,j,k) = sbvrf * buf(i,j,k)
!             vsal(i,j,k) = buf(i,j,k)
            enddo
          enddo
        enddo
!
      endif
    enddo
    deallocate (buf)
  endif
!
  if (num_vr0_files .gt. 0) then
    if (num_vr0_files .ne. num_cor_tracers) then
      call mpp_error(FATAL,'no. of var0 files != no. of tracers to correct')
    endif
    allocate (vtmp_s(isd:ied,jsd:jed))
    allocate (vsal_s(isd:ied,jsd:jed))
    allocate (buf2(isg:ieg,jsg:jeg))
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_tv0') then
        if (len_trim(data_table(nf)%file_name) .eq. 0) then
          buf2 = data_table(nf)%factor
        else
          call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                          access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
          read (nu) year, month, day
          read (nu) ig, jg
          read (nu) buf2
          close(nu)
        endif
!
        do j=jsd,jed
          do i=isd,ied
            vtmp_s(i,j) = tbv0f * buf2(i,j)
!           vtmp_s(i,j) = buf2(i,j)
          enddo
        enddo
!
      else if (data_table(nf)%fieldname_code(1:7) .eq. 'mdl_sv0') then
        if (len_trim(data_table(nf)%file_name) .eq. 0) then
          buf2 = data_table(nf)%factor
        else
          call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                          access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
          read (nu) year, month, day
          read (nu) ig, jg
          read (nu) buf2
          close(nu)
        endif
!
        do j=jsd,jed
          do i=isd,ied
            vsal_s(i,j) = sbv0f * buf2(i,j)
!           vsal_s(i,j) = buf2(i,j)
          enddo
        enddo
!
      endif
    enddo
    deallocate (buf2)
  else
    if (asm_sfc_split .and. .not.restore_sfc) then
      call mpp_error(FATAL,'surface assimilation to be split and var0 files are not provided')
    endif
  endif
!
  if (num_alt_files .eq. 2) then
    allocate (cdnz(kass))
    allocate (cdnzs(kass))
    allocate (eta_clm(isd:ied,jsd:jed))
    allocate (buf1(kass))
    allocate (buf2(isg:ieg,jsg:jeg))
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:4) .eq. 'cdnz') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        read (nu) buf1
        do k=1,kass
          cdnz(k) = buf1(k)
        enddo
        read (nu) buf1
        do k=1,kass
          cdnzs(k) = buf1(k)
        enddo
        close(nu)
        deallocate (buf1)
      else if (data_table(nf)%fieldname_code(1:4) .eq. 'etac') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        read (nu) buf2
        close(nu)
!
        do j=jsd,jed
          do i=isd,ied
            eta_clm(i,j) = buf2(i,j)
          enddo
        enddo
        deallocate (buf2)
      endif
    enddo
  endif
!
  allocate (ev(kass))
  allocate (wrkk(kass2))

! ----------------------------------------------
! register diagnostics
! ----------------------------------------------
!
  do n=1,num_cor_tracers
    if (n == index_temp) then
      T_cor(n)%name='tcor'
      T_cor(n)%units='Deg_C'
      T_cor(n)%longname='potential temperature correction'
      T_cor(n)%min_range=-10.0
      T_cor(n)%max_range=100.0
      T_cor(n)%init=.false.
      T_cor(n)%file_in='INPUT/ocean_cor.res.nc'
      T_cor(n)%file_out='RESTART/ocean_cor.res.nc'
      T_cor(n)%name_in='tcor'
    else if (n == index_salt) then
      T_cor(n)%name='scor'
      T_cor(n)%units='psu'
      T_cor(n)%longname='salinity correction'
      T_cor(n)%min_range=-10.0
      T_cor(n)%max_range=100.0
      T_cor(n)%init=.false.
      T_cor(n)%file_in='INPUT/ocean_cor.res.nc'
      T_cor(n)%file_out='RESTART/ocean_cor.res.nc'
      T_cor(n)%name_in='scor'
    endif
  enddo

  ! register diagnostics  (only if godas_at_end=.false.)

  if (.not. godas_at_end) then
    do n=1,num_cor_tracers
      range_array(1) = T_cor(n)%min_range
      range_array(2) = T_cor(n)%max_range
      id_cor(n) = register_diag_field ('ocean_model', trim(T_cor(n)%name), &
           Grd%tracer_axes(1:3),                                             &
           Time%model_time, trim(T_cor(n)%longname), trim(T_cor(n)%units), &
           missing_value=missing_value, range=range_array)
    enddo
  endif

  if (godas_at_end) then
    do n=1,num_cor_tracers
      if (.not. T_cor(n)%init) then
        write (stdout(),'(/a,a)') 'Expecting to read a GODAS restart file, ', T_cor(n)%file_in
      endif

      T_cor(n)%fcor(:,:,:) = 0.0
  
      if (file_exist(T_cor(n)%file_in)) then
        call read_data(T_cor(n)%file_in, T_cor(n)%name_in, T_cor(n)%fcor(:,:,:), Dom%domain2d, timelevel=1)

        if (.not. single_incr) then
          T_cor(n)%fcor(:,:,:) = scl_incr * T_cor(n)%fcor(:,:,:)
          write (stdout(),'(/a,1pe12.3)') 'GODAS restart increment rescaled: ',scl_incr
        endif

      else
        write (stdout(),'(/a)')'GODAS restart not found, increments set to zero.'
      endif
    enddo
  endif

  godas_module_initialized = .true.

end function godas_init
! </FUNCTION> NAME="godas_init">


!#######################################################################
! <SUBROUTINE NAME="godas_increment">
!
! <DESCRIPTION>
! Apply corrections from analysis.  Update analysis at specified interval.
! </DESCRIPTION>
!
subroutine godas_increment (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode
  type(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_Z(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_0(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_A(:)
  type(ocean_rstr_tracer_type), intent(inout)          :: T_rstr(:)

  integer         :: n, taup1, pe
  integer         :: year, month, day, hour, minute, second
integer :: i, j, k

! if godas_at_end=.true. no action is taken here. instead a single analysis is called
!  in subroutine godas_end

  if (.not. godas_at_end) then

    pe = mpp_pe()

    call get_date(Time%model_time, year, month, day, hour, minute, second)

    if (repeat_alarm(Time%model_time, gds_freq, alrm_dur) .or. ovr_alrm) then
      if (asm_ts_seq) then
        asm_code = temp_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
        asm_code = salt_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      else
        asm_code = ts_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      endif
      if (asm_sfc_split) then
        if (restore_sfc) then
          do n=1,num_cor_tracers
            T_cor(n)%fcor(:,:,1) = 0.0
          enddo
          call godas_rstr_comp (Time, T_prog, T_rstr)
        else
          asm_code = temp_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
          asm_code = salt_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
        endif
      endif
      apply_incr = .true.
      call godas_obs_track(Time, obs_Z, obs_0, obs_A)

!  scl_incr = (time step) / (time between analyses)
      if (.not. single_incr) then
        do n=1,num_cor_tracers
          T_cor(n)%fcor(:,:,:) = scl_incr * T_cor(n)%fcor(:,:,:)
        enddo
        if (restore_sfc) then
          do n=1,1,num_rstr_tracers
            T_rstr(n)%frstr(:,:) = scl_incr * T_rstr(n)%frstr(:,:)
          enddo
        endif
      endif

      if (ovr_alrm) then
        asm_cnt = 2
        ovr_alrm = .false.
      else
        asm_cnt = 1
      endif
    else if (asm_cnt .lt. no_asm_rep) then
      if (asm_ts_seq) then
        asm_code = temp_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
        asm_code = salt_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      else
        asm_code = ts_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      endif
      if (asm_sfc_split) then
        if (restore_sfc) then
          call godas_rstr_comp (Time, T_prog, T_rstr)
        else
          asm_code = temp_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
          asm_code = salt_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
        endif
      endif

!  scl_incr = (time step) / (time between analyses)
      if (.not. single_incr) then
        do n=1,num_cor_tracers
          T_cor(n)%fcor(:,:,:) = scl_incr * T_cor(n)%fcor(:,:,:)
        enddo
        if (restore_sfc) then
          do n=1,1,num_rstr_tracers
            T_rstr(n)%frstr(:,:) = scl_incr * T_rstr(n)%frstr(:,:)
          enddo
        endif
      endif

      asm_cnt = asm_cnt + 1
    else
      asm_cnt = asm_cnt + 1
    endif

! send increments/corrections to diag_manager

    do n=1,num_cor_tracers
      if (id_cor(n) > 0) used = send_data (id_cor(n), T_cor(n)%fcor(:,:,:), &
                                  Time%model_time,rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    enddo

! send restorations to diag_manager

    do n=1,num_rstr_tracers
      if (id_rstr(n) > 0) used = send_data (id_rstr(n), T_rstr(n)%frstr(:,:), &
                                  Time%model_time,rmask=Grd%tmask(:,:,1), &
                                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
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
      T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,taup1) + T_cor(n)%fcor(:,:,:)
  enddo

  if (restore_sfc) then
    do n=1,num_rstr_tracers
      T_prog(n)%field(:,:,1,taup1) = T_prog(n)%field(:,:,1,taup1) + T_rstr(n)%frstr(:,:)
    enddo
  endif

! if using the single application zero out the increments after one use
!
  if (single_incr) then
    do n=1,num_cor_tracers
       T_cor(n)%fcor(:,:,:) = 0.0
    enddo
    if (restore_sfc) then
      do n=1,num_rstr_tracers
        T_rstr(n)%frstr(:,:) = 0.0
      enddo
    endif
  endif

end subroutine godas_increment
! </SUBROUTINE> NAME="godas_increment">


!#######################################################################
! <SUBROUTINE NAME="godas_analysis">
!
! <DESCRIPTION>
! Perform an analysis.
! </DESCRIPTION>
!
subroutine godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode
  type(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_Z(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_0(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_A(:)

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

  call init_grad (Time, T_prog, Ext_mode, obs_Z, obs_0, obs_A)

!-----------------------------------------------------------------------
!  Do the first multiplication of the gradient by the background
!  error covariance matrix (E).
!       h^1 = E g^1
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

  call eg_lpsmthr ()

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

    call comp_f (obs_Z, obs_0, obs_A)

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

      call eg_lpsmthr ()

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
      T_cor(n)%fcor(:,:,:) = 0.0
      do k=1,kass
        T_cor(n)%fcor(:,:,k) = t_cg(:,:,k)
      enddo
    else if (n .eq. index_salt .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      T_cor(n)%fcor(:,:,:) = 0.0
      do k=1,kass
        T_cor(n)%fcor(:,:,k) = t_cg(:,:,k+ksalt)
      enddo
    endif
  enddo

! put increments and background error variance into obs_Z, obs_0, and obs_A

  do n=1,num_obsz
    if (obs_Z(n)%stat .eq. 2) then
      i = obs_Z(n)%io
      ip = i + 1
      j = obs_Z(n)%jo
      jp = j + 1
      if (obs_Z(n)%code .eq. temp_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
        do k=1,obs_Z(n)%kd
          obs_Z(n)%inc(k) = obs_Z(n)%a00 * T_cor(index_temp)%fcor(i,j,k) + &
                            obs_Z(n)%a01 * T_cor(index_temp)%fcor(i,jp,k) + &
                            obs_Z(n)%a11 * T_cor(index_temp)%fcor(ip,jp,k) + &
                            obs_Z(n)%a10 * T_cor(index_temp)%fcor(ip,j,k)
          obs_Z(n)%bke(k) = obs_Z(n)%a00 * vtmp(i,j,k) + &
                            obs_Z(n)%a01 * vtmp(i,jp,k) + &
                            obs_Z(n)%a11 * vtmp(ip,jp,k) + &
                            obs_Z(n)%a10 * vtmp(ip,j,k)
        enddo
        obs_Z(n)%stat = 3
      else if (obs_Z(n)%code .eq. salt_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
        do k=1,obs_Z(n)%kd
          obs_Z(n)%inc(k) = obs_Z(n)%a00 * T_cor(index_salt)%fcor(i,j,k) + &
                            obs_Z(n)%a01 * T_cor(index_salt)%fcor(i,jp,k) + &
                            obs_Z(n)%a11 * T_cor(index_salt)%fcor(ip,jp,k) + &
                            obs_Z(n)%a10 * T_cor(index_salt)%fcor(ip,j,k)
          obs_Z(n)%bke(k) = obs_Z(n)%a00 * vsal(i,j,k) + &
                            obs_Z(n)%a01 * vsal(i,jp,k) + &
                            obs_Z(n)%a11 * vsal(ip,jp,k) + &
                            obs_Z(n)%a10 * vsal(ip,j,k)
        enddo
        obs_Z(n)%stat = 3
      endif
    endif
  enddo

  if (.not. asm_sfc_split) then
    do n=1,num_obs0
      if (obs_0(n)%stat .eq. 2) then
        i = obs_0(n)%io
        ip = i + 1
        j = obs_0(n)%jo
        jp = j + 1
        if (obs_0(n)%code .eq. sst_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
          obs_0(n)%inc = obs_0(n)%a00 * T_cor(index_temp)%fcor(i,j,1) + &
                         obs_0(n)%a01 * T_cor(index_temp)%fcor(i,jp,1) + &
                         obs_0(n)%a11 * T_cor(index_temp)%fcor(ip,jp,1) + &
                         obs_0(n)%a10 * T_cor(index_temp)%fcor(ip,j,1)
          obs_0(n)%bke = obs_0(n)%a00 * vtmp(i,j,1) + &
                         obs_0(n)%a01 * vtmp(i,jp,1) + &
                         obs_0(n)%a11 * vtmp(ip,jp,1) + &
                         obs_0(n)%a10 * vtmp(ip,j,1)
          obs_0(n)%stat = 3
        else if (obs_0(n)%code .eq. sss_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
          obs_0(n)%inc = obs_0(n)%a00 * T_cor(index_salt)%fcor(i,j,1) + &
                         obs_0(n)%a01 * T_cor(index_salt)%fcor(i,jp,1) + &
                         obs_0(n)%a11 * T_cor(index_salt)%fcor(ip,jp,1) + &
                         obs_0(n)%a10 * T_cor(index_salt)%fcor(ip,j,1)
          obs_0(n)%bke = obs_0(n)%a00 * vsal(i,j,1) + &
                         obs_0(n)%a01 * vsal(i,jp,1) + &
                         obs_0(n)%a11 * vsal(ip,jp,1) + &
                         obs_0(n)%a10 * vsal(ip,j,1)
          obs_0(n)%stat = 3
        endif
      endif
    enddo
  endif

end subroutine godas_analysis
! </SUBROUTINE> NAME="godas_analysis"


!#######################################################################
! <SUBROUTINE NAME="godas_analysis_sfc">
!
! <DESCRIPTION>
! Perform a surface analysis. Temperature and salinity must be done individually.
! </DESCRIPTION>
!
subroutine godas_analysis_sfc (Time, T_prog, T_cor, obs_0)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
  type(ocean_obs0_type), intent(inout)                 :: obs_0(:)

  integer         :: pe, n
  integer         :: i, ip, j, jp, k, iter
  real            :: alpha, beta, gh_old, gh_new, df_val
integer :: ii, jj

  ni    = Grd%ni
  nj    = Grd%nj
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
  t_cg_s = 0.0

  call init_grad_sfc (Time, T_prog, obs_0)

!-----------------------------------------------------------------------
!  Do the first multiplication of the gradient by the background
!  error covariance matrix (E).
!       h^1 = E g^1
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

  call eg_lpsmthr_sfc ()

!-----------------------------------------------------------------------
!  Set the initial search directions to zero.
!       d^0 = 0
!       e^0 = 0
!-----------------------------------------------------------------------

  d_cg_s = 0.0
  e_cg_s = 0.0

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

    d_cg_s = beta * d_cg_s - h_cg_s
    e_cg_s = beta * e_cg_s - g_cg_s
    call mpp_update_domains (d_cg_s, Dom%domain2d)
    call mpp_update_domains (e_cg_s, Dom%domain2d)

!-----------------------------------------------------------------------
!  Compute f
!      f^n = e^n + trnsD invF D d^n
!-----------------------------------------------------------------------

    call comp_f_sfc (obs_0)

!-----------------------------------------------------------------------
!  Compute the inner products <g,h  and <d,f and update alpha
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------

    gh_new = mpp_global_sum(Dom%domain2d,g_cg_s(:,:)*h_cg_s(:,:)*Grd%tmask(:,:,1),BITWISE_EXACT_SUM)
    df_val = mpp_global_sum(Dom%domain2d,d_cg_s(:,:)*f_cg_s(:,:)*Grd%tmask(:,:,1),BITWISE_EXACT_SUM)

    alpha = gh_new / df_val

!-----------------------------------------------------------------------
!  Update the field correction (T) and the gradient (g)
!      T^(n+1) = T^n + alpha d^n
!      g^(n+1) = g^n + alpha f^n
!-----------------------------------------------------------------------

    t_cg_s = t_cg_s + alpha * d_cg_s
    call mpp_update_domains (t_cg_s, Dom%domain2d)

    if (iter .lt. maxits) then
      g_cg_s = g_cg_s + alpha * f_cg_s
      call mpp_update_domains (g_cg_s, Dom%domain2d)

!-----------------------------------------------------------------------
!  Update h by multiplying the new gradient ( g^(n+1) ) by the
!  background error covariance E.
!       h^(n+1) = E g^(n+1)
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

      call eg_lpsmthr_sfc ()

!-----------------------------------------------------------------------
!  Compute a new inner product <g,h and update beta
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------

      gh_old = gh_new
      gh_new = mpp_global_sum(Dom%domain2d,g_cg_s(:,:)*h_cg_s(:,:)*Grd%tmask(:,:,1),BITWISE_EXACT_SUM)

      beta = gh_new / gh_old

    endif
  enddo

!-----------------------------------------------------------------------
! Put increments into T_cor
!-----------------------------------------------------------------------

  if (asm_code .eq. temp_code) then
    T_cor(index_temp)%fcor(:,:,1) = t_cg_s(:,:)
  else if (asm_code .eq. salt_code) then
    T_cor(index_salt)%fcor(:,:,1) = t_cg_s(:,:)
  endif

!-----------------------------------------------------------------------
! put increments and background error variance into obs_0
!-----------------------------------------------------------------------

  do n=1,num_obs0
    if (obs_0(n)%stat .eq. 2) then
      i = obs_0(n)%io
      ip = i + 1
      j = obs_0(n)%jo
      jp = j + 1
      if (obs_0(n)%code .eq. sst_code .and. asm_code .eq. temp_code) then
        obs_0(n)%inc = obs_0(n)%a00 * T_cor(index_temp)%fcor(i,j,1) + &
                       obs_0(n)%a01 * T_cor(index_temp)%fcor(i,jp,1) + &
                       obs_0(n)%a11 * T_cor(index_temp)%fcor(ip,jp,1) + &
                       obs_0(n)%a10 * T_cor(index_temp)%fcor(ip,j,1)
        obs_0(n)%bke = obs_0(n)%a00 * vtmp_s(i,j) + &
                       obs_0(n)%a01 * vtmp_s(i,jp) + &
                       obs_0(n)%a11 * vtmp_s(ip,jp) + &
                       obs_0(n)%a10 * vtmp_s(ip,j)
        obs_0(n)%stat = 3
      else if (obs_0(n)%code .eq. sss_code .and. asm_code .eq. salt_code) then
        obs_0(n)%inc = obs_0(n)%a00 * T_cor(index_salt)%fcor(i,j,1) + &
                       obs_0(n)%a01 * T_cor(index_salt)%fcor(i,jp,1) + &
                       obs_0(n)%a11 * T_cor(index_salt)%fcor(ip,jp,1) + &
                       obs_0(n)%a10 * T_cor(index_salt)%fcor(ip,j,1)
        obs_0(n)%bke = obs_0(n)%a00 * vsal_s(i,j) + &
                       obs_0(n)%a01 * vsal_s(i,jp) + &
                       obs_0(n)%a11 * vsal_s(ip,jp) + &
                       obs_0(n)%a10 * vsal_s(ip,j)
        obs_0(n)%stat = 3
      endif
    endif
  enddo

end subroutine godas_analysis_sfc
! </SUBROUTINE> NAME="godas_analysis_sfc"


!#######################################################################
! <SUBROUTINE NAME="godas_end">
!
! <DESCRIPTION>
! Write GODAS restarts
! </DESCRIPTION>
!
subroutine godas_end(Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode
  type(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
  type(ocean_obsz_type), intent(inout)                 :: obs_Z(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_0(:)
  type(ocean_obs0_type), intent(inout)                 :: obs_A(:)
  type(ocean_rstr_tracer_type), intent(inout)          :: T_rstr(:)
! logical :: ens_ocean

  integer :: num_prog_tracers
  integer :: i, len, m, n, taup1, pe
  character(len=128) :: filename

  pe = mpp_pe()
  write(stdout(), '(a)') 'Later a GODAS restart will be written'

  if (godas_at_end) then

! next do a single GODAS analysis.  increments will be applied after next restart.

    if (asm_ts_seq .or. num_cor_tracers .eq. 1) then
      if (asm_Tz) then
        asm_code = temp_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      endif
      if (asm_Sz) then
        asm_code = salt_code
        call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
      endif
    else
      asm_code = ts_code
      call godas_analysis (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
    endif
    if (asm_sfc_split) then
      if (restore_sfc) then
        do n=1,num_cor_tracers
          T_cor(n)%fcor(:,:,1) = 0.0
        enddo
        call godas_rstr_comp (Time, T_prog, T_rstr)
      else
        if (asm_T0) then
          asm_code = temp_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
        endif
        if (asm_S0) then
          asm_code = salt_code
          call godas_analysis_sfc (Time, T_prog, T_cor, obs_0)
        endif
      endif
    endif
    call godas_obs_track(Time, obs_Z, obs_0, obs_A)

! write a restart file of GODAS increments

    do n=1,num_cor_tracers
!     call write_data(trim(T_cor(n)%file_out), trim(T_cor(n)%name), T_cor(n)%fcor(:,:,:), domain = Dom%domain2d, append_pelist_name = ens_ocean)
      call write_data(trim(T_cor(n)%file_out), trim(T_cor(n)%name),&
                T_cor(n)%fcor(:,:,:), domain = Dom%domain2d)
    enddo

    if (restore_sfc) then
      do n=1,num_rstr_tracers
!       call write_data(trim(T_rstr(n)%file_out), trim(T_rstr(n)%name), T_rstr(n)%frstr(:,:), domain = Dom%domain2d, append_pelist_name = ens_ocean)
        call write_data(trim(T_rstr(n)%file_out), trim(T_rstr(n)%name),&
                T_rstr(n)%frstr(:,:), domain = Dom%domain2d)
      enddo
    endif

  endif

end subroutine godas_end
! </SUBROUTINE> NAME="godas_end"


!!#######################################################################
!! <SUBROUTINE NAME="init_wghts">
!!
!! <DESCRIPTION>
!! Compute the the weights for the lp-smoother
!! </DESCRIPTION>
!!
subroutine init_wghts

  integer       :: pe, pes, npes, len
  integer       :: i, j, n, ii, jj, np, npid2
  integer       :: jbg, jfn, jgbg, jgfn
  real          :: cuj, cujm, ctjr
  real          :: re, re2, con, col, acon, wsnd(2)
  real          :: b2                    !  = (hrscl * pi / 360.0)**2


  pe = mpp_pe()
  npes = mpp_npes()
  ni    = Grd%ni
  nj    = Grd%nj
  nk    = Grd%nk

  allocate (xcb(npes),xce(npes),xcsz(npes))
  allocate (ycb(npes),yce(npes),ycsz(npes))
  call mpp_get_compute_domains(Dom%domain2d,xcb,xce,xcsz,ycb,yce,ycsz)

  aeval = 0.0
  do n=1,npits-1
    aeval=1./float(n)+aeval
  enddo
  aeval=aeval/float(npits)
  dbsq = 0.5*aeval*aeval

  re = 6370.0e3
  re2 = re**2

! hrscl = hrscl * 111.324e3
! b2 = 0.25*hrscl**2
  b2 = (hrscl * pi / 360.0)**2
  acon = b2/float(npits)

  allocate( s1(isd:ied,jsd:jed) )
  allocate( s2(isd:ied,jsd:jed) )
  allocate( wgns(isd:ied,jsd:jed) )
  s1 = 0.0
  s2 = 0.0
  wgns = 1.0

  allocate( elipt(isd:ied,jsd:jed) )
  elipt = 1.0
  do j = jsd, jed
    do i = isd,ied
      if (Grd%yt(i,j) .lt. -10.0) then
        elipt(i,j) = 1.0 - 0.75*exp(-((Grd%yt(i,j)+10.0)**2/900.0))
      else if (Grd%yt(i,j) .gt. 10.0) then
        elipt(i,j) = 1.0 - 0.75*exp(-((Grd%yt(i,j)-10.0)**2/900.0))
      else
        elipt(i,j) = 0.25
      endif
    enddo
  enddo

  allocate( wgta(isd:ied,jsd:jed) )
  allocate( wcn(isd:ied,jsd:jed) )
  allocate( wea(isd:ied,jsd:jed) )
  allocate( wwe(isd:ied,jsd:jed) )
  allocate( wno(isd:ied,jsd:jed) )
  allocate( wso(isd:ied,jsd:jed) )
  wgta(:,:) = 0.0
  wcn(:,:) = 0.0
  wea(:,:) = 0.0
  wwe(:,:) = 0.0
  wno(:,:) = 0.0
  wso(:,:) = 0.0
  do j = jsc,jec
    do i = isc,iec
      cuj = cos(Grd%phiu(i,j))
      cujm = cos(Grd%phiu(i,j-1))
      ctjr = 1.0/cos(Grd%phit(i,j))
      wea(i,j) = re2*cuj*cuj*ctjr*Grd%dxter(i,j)*Grd%dxtr(i,j)*acon/elipt(i,j)
      wno(i,j) = re2*cuj*cuj*ctjr*Grd%dytnr(i,j)*Grd%dytr(i,j)*acon
      wwe(i,j) = re2*cuj*cuj*ctjr*Grd%dxter(i-1,j)*Grd%dxtr(i,j)*acon/elipt(i,j)
      wso(i,j) = re2*cujm*cuj*ctjr*Grd%dytnr(i,j-1)*Grd%dytr(i,j)*acon
    enddo
  enddo
  call mpp_update_domains (wea(:,:), Dom%domain2d)
  call mpp_update_domains (wno(:,:), Dom%domain2d)
  call mpp_update_domains (wwe(:,:), Dom%domain2d)
  call mpp_update_domains (wso(:,:), Dom%domain2d)
  do j = jsd, jed
    do i = isd,ied
      wcn(i,j) = 1.0 - wso(i,j) - wno(i,j) - wea(i,j) - wwe(i,j)
    enddo
  enddo
  if (jsc .eq. jsg) then
    do i=isc,iec
      con = wso(i,jsc)
      col = con*con/((con+aeval)*con+dbsq)
      wcn(i,jsc) = wcn(i,jsc)+con*col
    enddo
  endif
  if (jec .eq. jeg) then
    do i=isc,iec
      con = wno(i,jec)
      col = con*con/((con+aeval)*con+dbsq)
      wcn(i,jec) = wcn(i,jec)+con*col
    enddo
  endif

  npid2=npits/2
  
  jgbg = 2
  jgfn = jeg - 1
  jbg = jsc
  if (jbg .eq. 1) jbg = 2
  jfn = jec
  if (jfn .ge. jemx) jfn = jemx - 1
 ii = (isg+ieg)/2
  do jj=jgbg,jgfn
!   do ii=isg,ieg
      s1(:,:) = 0.0
      s2(:,:) = 0.0
      pes = -1
      do n=1,npes
        if (ii .ge. xcb(n) .and. ii .le. xce(n) .and. jj .ge. ycb(n) .and. jj .le. yce(n)) then
          pes = n - 1
        endif
      enddo
      if (pes .lt. 0) call mpp_error(FATAL,'ERROR in INIT_WGHTS')
      if (pe .eq. pes) then
        s1(ii,jj) = 1.0
      endif
      call mpp_update_domains (s1, Dom%domain2d)
      do np=1,npid2
        do j=jbg,jfn
          do i=isc,iec
            s2(i,j) = wcn(i,j) * s1(i,j) + wso(i,j) * s1(i,j-1) + wno(i,j) * s1(i,j+1) &
                                           + wwe(i,j) * s1(i-1,j) + wea(i,j) * s1(i+1,j)
          enddo
        enddo
        call mpp_update_domains (s2, Dom%domain2d)
        do j=jbg,jfn
          do i=isc,iec
            s1(i,j) = wcn(i,j) * s2(i,j) + wno(i,j-1) * s2(i,j-1) + wso(i,j+1) * s2(i,j+1) &
                                           + wwe(i,j) * s2(i-1,j) + wea(i,j) * s2(i+1,j)
          enddo
        enddo
        call mpp_update_domains (s1, Dom%domain2d)
      enddo
      wsnd = 0.0
      if (pe .eq. pes) then
        wsnd(1) = s1(ii,jj)
      endif
      len = 2
!     call mpp_transmit (wsnd,len,ALL_PES,wsnd,len,pes)
      call mpp_broadcast(wsnd,len,pes)
      if (jj .ge. jbg .and. jj .le. jfn) then
        do i=isc,iec
          wgns(i,jj) = wsnd(1)
        enddo
      endif
!   enddo
  enddo
  call mpp_sync()

  do j=jbg,jfn
    do i=isc,iec
      if (wgns(i,j) .lt. 0.0) then
        wgns(i,j) = 0.00001
      endif
      wgta(i,j)=sqrt(1.0/wgns(i,j))
!     wgta(i,j)=sqrt(tbvrf/wgns(i,j))
    enddo
  enddo
  call mpp_update_domains (wgta, Dom%domain2d)

  if (asm_sfc_split) then

    allocate( wgta_s(isd:ied,jsd:jed) )
    allocate( wcn_s(isd:ied,jsd:jed) )
    allocate( wea_s(isd:ied,jsd:jed) )
    allocate( wwe_s(isd:ied,jsd:jed) )
    allocate( wno_s(isd:ied,jsd:jed) )
    allocate( wso_s(isd:ied,jsd:jed) )

    if (abs(hrscl - hrscl0) .lt. 0.05) then
      wgta_s = wgta
      wcn_s = wcn
      wea_s = wea
      wwe_s = wwe
      wno_s = wno
      wso_s = wso

    else

      b2 = (hrscl0 * pi / 360.0)**2
      acon = b2/float(npits)

      wgta_s(:,:) = 0.0
      wcn_s(:,:) = 0.0
      wea_s(:,:) = 0.0
      wwe_s(:,:) = 0.0
      wno_s(:,:) = 0.0
      wso_s(:,:) = 0.0
      do j = jsc,jec
        do i = isc,iec
          cuj = cos(Grd%phiu(i,j))
          cujm = cos(Grd%phiu(i,j-1))
          ctjr = 1.0/cos(Grd%phit(i,j))
          wea_s(i,j) = re2*cuj*cuj*ctjr*Grd%dxter(i,j)*Grd%dxtr(i,j)*acon/elipt(i,j)
          wno_s(i,j) = re2*cuj*cuj*ctjr*Grd%dytnr(i,j)*Grd%dytr(i,j)*acon
          wwe_s(i,j) = re2*cuj*cuj*ctjr*Grd%dxter(i-1,j)*Grd%dxtr(i,j)*acon/elipt(i,j)
          wso_s(i,j) = re2*cujm*cuj*ctjr*Grd%dytnr(i,j-1)*Grd%dytr(i,j)*acon
        enddo
      enddo
      call mpp_update_domains (wea_s(:,:), Dom%domain2d)
      call mpp_update_domains (wno_s(:,:), Dom%domain2d)
      call mpp_update_domains (wwe_s(:,:), Dom%domain2d)
      call mpp_update_domains (wso_s(:,:), Dom%domain2d)
      do j = jsd, jed
        do i = isd,ied
          wcn_s(i,j) = 1.0 - wso_s(i,j) - wno_s(i,j) - wea_s(i,j) - wwe_s(i,j)
        enddo
      enddo
      if (jsc .eq. jsg) then
        do i=isc,iec
          con = wso_s(i,jsc)
          col = con*con/((con+aeval)*con+dbsq)
          wcn_s(i,jsc) = wcn_s(i,jsc)+con*col
        enddo
      endif
      if (jec .eq. jeg) then
        do i=isc,iec
          con = wno_s(i,jec)
          col = con*con/((con+aeval)*con+dbsq)
          wcn_s(i,jec) = wcn_s(i,jec)+con*col
        enddo
      endif

      npid2=npits/2

      jgbg = 2
      jgfn = jeg - 1
      jbg = jsc
      if (jbg .eq. 1) jbg = 2
      jfn = jec
      if (jfn .ge. jemx) jfn = jemx - 1
     ii = (isg+ieg)/2
      do jj=jgbg,jgfn
!       do ii=isg,ieg
          s1(:,:) = 0.0
          s2(:,:) = 0.0
          pes = -1
          do n=1,npes
            if (ii .ge. xcb(n) .and. ii .le. xce(n) .and. jj .ge. ycb(n) .and. jj .le. yce(n)) then
              pes = n - 1
            endif
          enddo
          if (pes .lt. 0) call mpp_error(FATAL,'ERROR in INIT_WGHTS')
          if (pe .eq. pes) then
            s1(ii,jj) = 1.0
          endif
          call mpp_update_domains (s1, Dom%domain2d)
          do np=1,npid2
            do j=jbg,jfn
              do i=isc,iec
                s2(i,j) = wcn_s(i,j) * s1(i,j) + wso_s(i,j) * s1(i,j-1) + wno_s(i,j) * s1(i,j+1) &
                                               + wwe_s(i,j) * s1(i-1,j) + wea_s(i,j) * s1(i+1,j)
              enddo
            enddo
            call mpp_update_domains (s2, Dom%domain2d)
            do j=jbg,jfn
              do i=isc,iec
                s1(i,j) = wcn_s(i,j) * s2(i,j) + wno_s(i,j-1) * s2(i,j-1) + wso_s(i,j+1) * s2(i,j+1) &
                                               + wwe_s(i,j) * s2(i-1,j) + wea_s(i,j) * s2(i+1,j)
              enddo
            enddo
            call mpp_update_domains (s1, Dom%domain2d)
          enddo
          wsnd = 0.0
          if (pe .eq. pes) then
            wsnd(1) = s1(ii,jj)
          endif
          len = 2
!         call mpp_transmit (wsnd,len,ALL_PES,wsnd,len,pes)
          call mpp_broadcast(wsnd,len,pes)
          if (jj .ge. jbg .and. jj .le. jfn) then
            do i=isc,iec
              wgns(i,jj) = wsnd(1)
            enddo
          endif
!       enddo
      enddo
      call mpp_sync()

      do j=jbg,jfn
        do i=isc,iec
          if (wgns(i,j) .lt. 0.0) then
            wgns(i,j) = 0.00001
          endif
          wgta_s(i,j)=sqrt(1.0/wgns(i,j))
!         wgta_s(i,j)=sqrt(tbvrf/wgns(i,j))
        enddo
      enddo
      call mpp_update_domains (wgta_s, Dom%domain2d)
    endif
  endif

 deallocate (elipt)
  deallocate(wgns)

end subroutine init_wghts
! </SUBROUTINE> NAME="init_wghts"


!!
!!#######################################################################
!! <SUBROUTINE NAME="vertical_covariance"
!!
!! <DESCRIPTION>
!! Compute the vertical covariance matrix.
!! A smoother is used.  The scale is proportional to layer thickness.
!! </DESCRIPTION>
!!
  subroutine vertical_covariance

  integer, parameter                :: kex = 10, nitz = 100
  integer                           :: k, kmp, kmpm1, kssm1, ksg
  integer                           :: ks, kk, n, ni2
  real                              :: scl, ascl 
  real, allocatable, dimension(:)   :: wm, wp, w, scz, scz2
!
  kmp = kass + 2*kex
  kmpm1 = kmp - 1
  kssm1 = kass - 1
!
  allocate ( wm(kmp) )
  allocate ( wp(kmp) )
  allocate ( w(kmp) )
  allocate ( scz(kmp) )
  allocate ( scz2(kmp) )
!
  do k=1,kssm1
    scl = vsclf*Grd%dzt(k)
    ascl = scl*scl / float(nitz);
    if (k .eq. 1) then
      wm(k+kex) = ascl / (2.0*Grd%dzw(k-1)*Grd%dzt(k))
    else
      wm(k+kex) = ascl / (Grd%dzw(k-1)*Grd%dzt(k))
    endif
    wp(k+kex) = ascl / (Grd%dzw(k)*Grd%dzt(k))
    w(k+kex) = 1.0 - wm(k+kex) - wp(k+kex)
  enddo
  do k=1,kex
    wm(k) = wm(kex+1)
    wp(k) = wp(kex+1)
    w(k) = w(kex+1)
  enddo
  do k=kass+kex,kmp
    wm(k) = wm(kex+kssm1)
    wp(k) = wp(kex+kssm1)
    w(k) = w(kex+kssm1)
  enddo
!
  ni2 = nitz / 2
!
  do ks=1,kass
    do k=1,kmp
      scz(k) = 0.0
    enddo
    ksg = ks + kex
    scz(ksg) = 1.0
    do n=1,ni2
      do k=2,kmpm1
        scz2(k) = w(k)*scz(k) + wm(k)*scz(k-1) + wp(k)*scz(k+1)
      enddo
      scz2(1) = scz2(2)
      scz2(kmp) = scz2(kmpm1)

      do k=2,kmpm1
       scz(k) = w(k)*scz2(k) + wp(k-1)*scz2(k-1) + wm(k+1)*scz2(k+1)
      enddo
      scz(1) = scz(2)
      scz(kmp) = scz(kmpm1)
    enddo
!
    ascl = scz(ksg)
    do k=1,kmp
      scz(k) = scz(k) / ascl
    enddo
!
    do k=1,kass
      cvn(ks,k) = scz(k+kex)
    enddo
  enddo
!
  do k=1,kass
    cvn(k,k) = vcvn
    do kk=k+1,kass
      ascl = 0.5*(cvn(k,kk) + cvn(kk,k))*vcvn
      cvn(k,kk) = ascl
      cvn(kk,k) = ascl
    enddo
  enddo
!
  cvnsalt = cvn
!
  deallocate ( wm )
  deallocate ( wp )
  deallocate ( w )
  deallocate ( scz )
  deallocate ( scz2 )
!
  end subroutine vertical_covariance
! </SUBROUTINE> NAME="vertical_covariance"

!!
!!#######################################################################
!! <SUBROUTINE NAME="init_grad">
!!
!! <DESCRIPTION>
!! Compute the initial estimate of the gradient of the functional (g)
!! </DESCRIPTION>
!!
subroutine init_grad (Time, T_prog, Ext_mode, obs_Z, obs_0, obs_A)
!
  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
  type(ocean_external_mode_type), intent(in)        :: Ext_mode
  type(ocean_obsz_type), intent(inout)              :: obs_Z(:)
  type(ocean_obs0_type), intent(inout)              :: obs_0(:)
  type(ocean_obs0_type), intent(inout)              :: obs_A(:)
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
!      SST        21 <= code <= 22
!      SSS        23 <= code <= 23
!      ALTM       24 <= code <= 26
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
  do n=1,num_obsz
    i = obs_Z(n)%io
    ip = i + 1
    j = obs_Z(n)%jo
    jp = j + 1
    if (obs_Z(n)%code .eq. temp_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, tz_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, tz_wndw_bwd)
      if (obs_Z(n)%obs_time < wndw_fwd .and. obs_Z(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_Z(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rtzw)
        obs_Z(n)%win = .true.
        if (save_all_inv) then
          obs_Z(n)%stat = 1
        else
          if (obs_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_Z(n)%stat = 1
        endif
        do k=1,obs_Z(n)%kd
          ov = obs_Z(n)%a00 * T_prog(index_temp)%field(i,j,k,taup1) + &
               obs_Z(n)%a01 * T_prog(index_temp)%field(i,jp,k,taup1) + &
               obs_Z(n)%a11 * T_prog(index_temp)%field(ip,jp,k,taup1) + &
               obs_Z(n)%a10 * T_prog(index_temp)%field(ip,j,k,taup1)
          ov = ov - obs_Z(n)%val(k)
          if (obs_Z(n)%stat .eq. 1) then
            obs_Z(n)%inv(k) = ov
            if (k .eq. obs_Z(n)%kd) obs_Z(n)%stat = 2
          endif
          aerr = obs_Z(n)%err(k)*time_adj
          aov = abs(ov)
          if (aov .lt. dtemp_max) then
            if (aov .gt. dtemp_elm) then
              aerr = aerr/(1.0+aov-dtemp_elm)**2
            endif
            g_cg(i,j,k) = g_cg(i,j,k) + ov*aerr*obs_Z(n)%a00
            g_cg(i,jp,k) = g_cg(i,jp,k) + ov*aerr*obs_Z(n)%a01
            g_cg(ip,jp,k) = g_cg(ip,jp,k) + ov*aerr*obs_Z(n)%a11
            g_cg(ip,j,k) = g_cg(ip,j,k) + ov*aerr*obs_Z(n)%a10
            obs_Z(n)%aerr(k) = aerr
          else
            obs_Z(n)%aerr(k) = 0.0
          endif
        enddo
      else
        time_adj = 0.0
        obs_Z(n)%win = .false.
      endif
    else if (obs_Z(n)%code .eq. salt_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, sz_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, sz_wndw_bwd)
      if (obs_Z(n)%obs_time < wndw_fwd .and. obs_Z(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_Z(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rszw)
        obs_Z(n)%win = .true.
        if (save_all_inv) then
          obs_Z(n)%stat = 1
        else
          if (obs_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_Z(n)%stat = 1
        endif
        do k=1,obs_Z(n)%kd
          ov = obs_Z(n)%a00 * T_prog(index_salt)%field(i,j,k,taup1) + &
               obs_Z(n)%a01 * T_prog(index_salt)%field(i,jp,k,taup1) + &
               obs_Z(n)%a11 * T_prog(index_salt)%field(ip,jp,k,taup1) + &
               obs_Z(n)%a10 * T_prog(index_salt)%field(ip,j,k,taup1)
          ov = ov - obs_Z(n)%val(k)
          if (obs_Z(n)%stat .eq. 1) then
            obs_Z(n)%inv(k) = ov
            if (k .eq. obs_Z(n)%kd) obs_Z(n)%stat = 2
          endif
          aerr = obs_Z(n)%err(k)*time_adj
          aov = abs(ov)
          if (aov .lt. dsalt_max) then
            if (aov .gt. dsalt_elm) then
              aerr = aerr/(1.0+aov-dsalt_elm)**2
            endif
            g_cg(i,j,k+ksalt) = g_cg(i,j,k+ksalt) + ov*aerr*obs_Z(n)%a00
            g_cg(i,jp,k+ksalt) = g_cg(i,jp,k+ksalt) + ov*aerr*obs_Z(n)%a01
            g_cg(ip,jp,k+ksalt) = g_cg(ip,jp,k+ksalt) + ov*aerr*obs_Z(n)%a11
            g_cg(ip,j,k+ksalt) = g_cg(ip,j,k+ksalt) + ov*aerr*obs_Z(n)%a10
            obs_Z(n)%aerr(k) = aerr
          else
            obs_Z(n)%aerr(k) = 0.0
          endif
        enddo
      else
        time_adj = 0.0
        obs_Z(n)%win = .false.
      endif
    endif
  enddo
!
!-----------------------------------------------------------------------
!  The surface observations (exlcuding altimetry)
!-----------------------------------------------------------------------
!
  if (.not. asm_sfc_split) then
    do n=1,num_obs0
      i = obs_0(n)%io
      ip = i + 1
      j = obs_0(n)%jo
      jp = j + 1
      if (obs_0(n)%code .eq. sst_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
        wndw_fwd = increment_time(Time%model_time, wndw_secs, t0_wndw_fwd)
        wndw_bwd = decrement_time(Time%model_time, wndw_secs, t0_wndw_bwd)
        if (obs_0(n)%obs_time < wndw_fwd .and. obs_0(n)%obs_time > wndw_bwd) then
          diff_time = Time%model_time - obs_0(n)%obs_time
          call get_time (diff_time, dsec, dday)
          time_sep = real(dday) + real(dsec)/real(spd)
          time_adj = (1.0-time_sep*rt0w)
          obs_0(n)%win = .true.
          if (save_all_inv) then
            obs_0(n)%stat = 1
          else
            if (obs_0(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_0(n)%stat = 1
          endif
          ov = obs_0(n)%a00 * T_prog(index_temp)%field(i,j,1,taup1) + &
               obs_0(n)%a01 * T_prog(index_temp)%field(i,jp,1,taup1) + &
               obs_0(n)%a11 * T_prog(index_temp)%field(ip,jp,1,taup1) + &
               obs_0(n)%a10 * T_prog(index_temp)%field(ip,j,1,taup1)
          ov = ov - obs_0(n)%val
          if (obs_0(n)%stat .eq. 1) then
            obs_0(n)%inv = ov
            obs_0(n)%stat = 2
          endif
          aerr = obs_0(n)%err*time_adj
          aov = abs(ov)
          if (aov .lt. dtemp_max) then
            if (aov .gt. dtemp_elm) then
              aerr = aerr/(1.0+aov-dtemp_elm)**2
            endif
            g_cg(i,j,1) = g_cg(i,j,1) + ov*aerr*obs_0(n)%a00
            g_cg(i,jp,1) = g_cg(i,jp,1) + ov*aerr*obs_0(n)%a01
            g_cg(ip,jp,1) = g_cg(ip,jp,1) + ov*aerr*obs_0(n)%a11
            g_cg(ip,j,1) = g_cg(ip,j,1) + ov*aerr*obs_0(n)%a10
            obs_0(n)%aerr = aerr
          else
            obs_0(n)%aerr = 0.0
          endif
        else
          time_adj = 0.0
          obs_0(n)%win = .false.
        endif
      else if (obs_0(n)%code .eq. sss_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
        wndw_fwd = increment_time(Time%model_time, wndw_secs, s0_wndw_fwd)
        wndw_bwd = decrement_time(Time%model_time, wndw_secs, s0_wndw_bwd)
        if (obs_0(n)%obs_time < wndw_fwd .and. obs_0(n)%obs_time > wndw_bwd) then
          diff_time = Time%model_time - obs_0(n)%obs_time
          call get_time (diff_time, dsec, dday)
          time_sep = real(dday) + real(dsec)/real(spd)
          time_adj = (1.0-time_sep*rs0w)
          obs_0(n)%win = .true.
          if (save_all_inv) then
            obs_0(n)%stat = 1
          else
            if (obs_0(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_0(n)%stat = 1
          endif
          ov = obs_0(n)%a00 * T_prog(index_salt)%field(i,j,1,taup1) + &
               obs_0(n)%a01 * T_prog(index_salt)%field(i,jp,1,taup1) + &
               obs_0(n)%a11 * T_prog(index_salt)%field(ip,jp,1,taup1) + &
               obs_0(n)%a10 * T_prog(index_salt)%field(ip,j,1,taup1)
          ov = ov - obs_0(n)%val
          if (obs_0(n)%stat .eq. 1) then
            obs_0(n)%inv = ov
            obs_0(n)%stat = 2
          endif
          aerr = obs_0(n)%err*time_adj
          aov = abs(ov)
          if (aov .lt. dsalt_max) then
            if (aov .gt. dsalt_elm) then
              aerr = aerr/(1.0+aov-dsalt_elm)**2
            endif
            g_cg(i,j,1+ksalt) = g_cg(i,j,1+ksalt) + ov*aerr*obs_0(n)%a00
            g_cg(i,jp,1+ksalt) = g_cg(i,jp,1+ksalt) + ov*aerr*obs_0(n)%a01
            g_cg(ip,jp,1+ksalt) = g_cg(ip,jp,1+ksalt) + ov*aerr*obs_0(n)%a11
            g_cg(ip,j,1+ksalt) = g_cg(ip,j,1+ksalt) + ov*aerr*obs_0(n)%a10
            obs_0(n)%aerr = aerr
          else
            obs_0(n)%aerr = 0.0
          endif
        else
          time_adj = 0.0
          obs_0(n)%win = .false.
        endif
      endif
    enddo
  endif
!
!-----------------------------------------------------------------------
!  The altimeter observations
!-----------------------------------------------------------------------
!
  do n=1,num_obsa
    i = obs_A(n)%io
    ip = i + 1
    j = obs_A(n)%jo
    jp = j + 1
    if (obs_A(n)%code .eq. altm_code) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, al_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, al_wndw_bwd)
      if (obs_A(n)%obs_time < wndw_fwd .and. obs_A(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_A(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*ralw)
        obs_A(n)%win = .true.
        if (save_all_inv) then
          obs_A(n)%stat = 1
        else
          if (obs_A(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_A(n)%stat = 1
        endif
        ov = obs_A(n)%a00 * (Ext_mode%eta_t(i,j,taup1) - eta_clm(i,j)) + &
             obs_A(n)%a01 * (Ext_mode%eta_t(i,jp,taup1) - eta_clm(i,jp)) + &
             obs_A(n)%a11 * (Ext_mode%eta_t(ip,jp,taup1) - eta_clm(ip,jp)) + &
             obs_A(n)%a10 * (Ext_mode%eta_t(ip,j,taup1) - eta_clm(ip,j))
        ov = ov - obs_A(n)%val
        if (obs_A(n)%stat .eq. 1) then
          obs_A(n)%inv = ov
          obs_A(n)%stat = 2
        endif
        aerr = obs_A(n)%err*time_adj
        aov = abs(ov)
        if (aov .lt. daltm_max) then
          if (aov .gt. daltm_elm) then
            aerr = aerr/(1.0+aov-daltm_elm)**2
          endif
          do k=1,kass
            g_cg(i,j,k) = g_cg(i,j,k) + ov*aerr*obs_A(n)%a00*cdnz(k)
            g_cg(i,jp,k) = g_cg(i,jp,k) + ov*aerr*obs_A(n)%a01*cdnz(k)
            g_cg(ip,jp,k) = g_cg(ip,jp,k) + ov*aerr*obs_A(n)%a11*cdnz(k)
            g_cg(ip,j,k) = g_cg(ip,j,k) + ov*aerr*obs_A(n)%a10*cdnz(k)
            kks = k+ksalt
            g_cg(i,j,kks) = g_cg(i,j,kks) + ov*aerr*obs_A(n)%a00*cdnzs(k)
            g_cg(i,jp,kks) = g_cg(i,jp,kks) + ov*aerr*obs_A(n)%a01*cdnzs(k)
            g_cg(ip,jp,kks) = g_cg(ip,jp,kks) + ov*aerr*obs_A(n)%a11*cdnzs(k)
            g_cg(ip,j,kks) = g_cg(ip,j,kks) + ov*aerr*obs_A(n)%a10*cdnzs(k)
          enddo
          obs_A(n)%aerr = aerr
        else
          obs_A(n)%aerr = 0.0
        endif
      else
        time_adj = 0.0
        obs_A(n)%win = .false.
      endif
    endif
  enddo
!
  do k=1,kass2
    call mpp_update_domains (g_cg(:,:,k), Dom%domain2d)
  enddo
!
  end subroutine init_grad
! </SUBROUTINE> NAME="init_grad"

!!
!!#######################################################################
!! <SUBROUTINE NAME="init_grad_sfc">
!!
!! <DESCRIPTION>
!! Compute the initial estimate of the gradient of the functional (g)
!! </DESCRIPTION>
!!
subroutine init_grad_sfc (Time, T_prog, obs_0)
!
  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
  type(ocean_obs0_type), intent(inout)              :: obs_0(:)
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
!      SST        21 <= code <= 22
!      SSS        23 <= code <= 23
!  these are set in godas_data_mod, and can be modified if needed
!-----------------------------------------------------------------------
!
  pe = mpp_pe()

  taup1   = Time%taup1

!-----------------------------------------------------------------------
!   Set g to zero.
!-----------------------------------------------------------------------
!
  g_cg_s = 0.0
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
!   6) sum the gridded result in g_cg_s
!-----------------------------------------------------------------------
!  The profile observations
!-----------------------------------------------------------------------
!
 tobc = 0
 tsum = 0.0
 sobc = 0
 ssum = 0.0
!
!-----------------------------------------------------------------------
!  The surface observations
!-----------------------------------------------------------------------
!
  do n=1,num_obs0
    i = obs_0(n)%io
    ip = i + 1
    j = obs_0(n)%jo
    jp = j + 1
    if (obs_0(n)%code .eq. sst_code .and. asm_code .eq. temp_code) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, t0_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, t0_wndw_bwd)
      if (obs_0(n)%obs_time < wndw_fwd .and. obs_0(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_0(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rt0w)
        obs_0(n)%win = .true.
        if (save_all_inv) then
          obs_0(n)%stat = 1
        else
          if (obs_0(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_0(n)%stat = 1
        endif
        ov = obs_0(n)%a00 * T_prog(index_temp)%field(i,j,1,taup1) + &
             obs_0(n)%a01 * T_prog(index_temp)%field(i,jp,1,taup1) + &
             obs_0(n)%a11 * T_prog(index_temp)%field(ip,jp,1,taup1) + &
             obs_0(n)%a10 * T_prog(index_temp)%field(ip,j,1,taup1)
        ov = ov - obs_0(n)%val
        if (obs_0(n)%stat .eq. 1) then
          obs_0(n)%inv = ov
          obs_0(n)%stat = 2
        endif
        aerr = obs_0(n)%err*time_adj
        aov = abs(ov)
        if (aov .lt. dsst_max) then
          if (aov .gt. dsst_elm) then
            aerr = aerr/(1.0+aov-dsst_elm)**2
          endif
          g_cg_s(i,j) = g_cg_s(i,j) + ov*aerr*obs_0(n)%a00
          g_cg_s(i,jp) = g_cg_s(i,jp) + ov*aerr*obs_0(n)%a01
          g_cg_s(ip,jp) = g_cg_s(ip,jp) + ov*aerr*obs_0(n)%a11
          g_cg_s(ip,j) = g_cg_s(ip,j) + ov*aerr*obs_0(n)%a10
          obs_0(n)%aerr = aerr
        else
          obs_0(n)%aerr = 0.0
        endif
      else
        time_adj = 0.0
        obs_0(n)%win = .false.
      endif
    else if (obs_0(n)%code .eq. sss_code .and. asm_code .eq. salt_code) then
      wndw_fwd = increment_time(Time%model_time, wndw_secs, s0_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, s0_wndw_bwd)
      if (obs_0(n)%obs_time < wndw_fwd .and. obs_0(n)%obs_time > wndw_bwd) then
        diff_time = Time%model_time - obs_0(n)%obs_time
        call get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rs0w)
        obs_0(n)%win = .true.
        if (save_all_inv) then
          obs_0(n)%stat = 1
        else
          if (obs_0(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_0(n)%stat = 1
        endif
        ov = obs_0(n)%a00 * T_prog(index_salt)%field(i,j,1,taup1) + &
             obs_0(n)%a01 * T_prog(index_salt)%field(i,jp,1,taup1) + &
             obs_0(n)%a11 * T_prog(index_salt)%field(ip,jp,1,taup1) + &
             obs_0(n)%a10 * T_prog(index_salt)%field(ip,j,1,taup1)
        ov = ov - obs_0(n)%val
        if (obs_0(n)%stat .eq. 1) then
          obs_0(n)%inv = ov
          obs_0(n)%stat = 2
        endif
        aerr = obs_0(n)%err*time_adj
        aov = abs(ov)
        if (aov .lt. dsss_max) then
          if (aov .gt. dsss_elm) then
            aerr = aerr/(1.0+aov-dsss_elm)**2
          endif
          g_cg_s(i,j) = g_cg_s(i,j) + ov*aerr*obs_0(n)%a00
          g_cg_s(i,jp) = g_cg_s(i,jp) + ov*aerr*obs_0(n)%a01
          g_cg_s(ip,jp) = g_cg_s(ip,jp) + ov*aerr*obs_0(n)%a11
          g_cg_s(ip,j) = g_cg_s(ip,j) + ov*aerr*obs_0(n)%a10
          obs_0(n)%aerr = aerr
        else
          obs_0(n)%aerr = 0.0
        endif
      else
        time_adj = 0.0
        obs_0(n)%win = .false.
      endif
    endif
  enddo
!
  call mpp_update_domains (g_cg_s, Dom%domain2d)
!
  end subroutine init_grad_sfc
! </SUBROUTINE> NAME="init_grad_sfc"


!!
!!#######################################################################
!! <SUBROUTINE NAME="comp_f">
!!
!! <DESCRIPTION>
!! Compute the 
!! </DESCRIPTION>
!!
subroutine comp_f (obs_Z, obs_0, obs_A)
!
  type(ocean_obsz_type), intent(in)        :: obs_Z(:)
  type(ocean_obs0_type), intent(in)        :: obs_0(:)
  type(ocean_obs0_type), intent(in)        :: obs_A(:)
!
  integer :: i, ip, j, jp, k, ks1, n, pe
  real    :: ov, aerr
integer :: noot, noos
!-----------------------------------------------------------------------
!  data types are encoded in obs%code
!      T(z)        1 <= code <= 10
!      S(z)       11 <= code <= 20
!      SST        21 <= code <= 22
!      SSS        23 <= code <= 23
!      ALTM       24 <= code <= 26
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
  do n=1,num_obsz
    i = obs_Z(n)%io
    ip = i + 1
    j = obs_Z(n)%jo
    jp = j + 1
    if (obs_Z(n)%code .eq. temp_code .and. obs_Z(n)%win .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      do k=1,obs_Z(n)%kd
        ov = obs_Z(n)%a00 * d_cg(i,j,k) + &
             obs_Z(n)%a01 * d_cg(i,jp,k) + &
             obs_Z(n)%a11 * d_cg(ip,jp,k) + &
             obs_Z(n)%a10 * d_cg(ip,j,k)
        aerr = obs_Z(n)%aerr(k)
        f_cg(i,j,k) = f_cg(i,j,k) + ov*aerr*obs_Z(n)%a00
        f_cg(i,jp,k) = f_cg(i,jp,k) + ov*aerr*obs_Z(n)%a01
        f_cg(ip,jp,k) = f_cg(ip,jp,k) + ov*aerr*obs_Z(n)%a11
        f_cg(ip,j,k) = f_cg(ip,j,k) + ov*aerr*obs_Z(n)%a10
      enddo
    else if (obs_Z(n)%code .eq. salt_code .and. obs_Z(n)%win .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      do k=1,obs_Z(n)%kd
        ov = obs_Z(n)%a00 * d_cg(i,j,k+ksalt) + &
             obs_Z(n)%a01 * d_cg(i,jp,k+ksalt) + &
             obs_Z(n)%a11 * d_cg(ip,jp,k+ksalt) + &
             obs_Z(n)%a10 * d_cg(ip,j,k+ksalt)
        aerr = obs_Z(n)%aerr(k)
        f_cg(i,j,k+ksalt) = f_cg(i,j,k+ksalt) + ov*aerr*obs_Z(n)%a00
        f_cg(i,jp,k+ksalt) = f_cg(i,jp,k+ksalt) + ov*aerr*obs_Z(n)%a01
        f_cg(ip,jp,k+ksalt) = f_cg(ip,jp,k+ksalt) + ov*aerr*obs_Z(n)%a11
        f_cg(ip,j,k+ksalt) = f_cg(ip,j,k+ksalt) + ov*aerr*obs_Z(n)%a10
      enddo
    endif
  enddo
!
!-----------------------------------------------------------------------
!  Surface observations (excluding altimetry)
!-----------------------------------------------------------------------
!
noot = 0
noos = 0
  if (.not. asm_sfc_split) then
    do n=1,num_obs0
      i = obs_0(n)%io
      ip = i + 1
      j = obs_0(n)%jo
      jp = j + 1
      if (obs_0(n)%code .eq. sst_code .and. obs_0(n)%win .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
  noot = noot+1
        ov = obs_0(n)%a00 * d_cg(i,j,1) + &
             obs_0(n)%a01 * d_cg(i,jp,1) + &
             obs_0(n)%a11 * d_cg(ip,jp,1) + &
             obs_0(n)%a10 * d_cg(ip,j,1)
        aerr = obs_0(n)%aerr
        f_cg(i,j,1) = f_cg(i,j,1) + ov*aerr*obs_0(n)%a00
        f_cg(i,jp,1) = f_cg(i,jp,1) + ov*aerr*obs_0(n)%a01
        f_cg(ip,jp,1) = f_cg(ip,jp,1) + ov*aerr*obs_0(n)%a11
        f_cg(ip,j,1) = f_cg(ip,j,1) + ov*aerr*obs_0(n)%a10
      else if (obs_0(n)%code .eq. sss_code .and. obs_0(n)%win .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
  noos = noos+1
        ks1 = ksalt + 1
        ov = obs_0(n)%a00 * d_cg(i,j,ks1) + &
             obs_0(n)%a01 * d_cg(i,jp,ks1) + &
             obs_0(n)%a11 * d_cg(ip,jp,ks1) + &
             obs_0(n)%a10 * d_cg(ip,j,ks1)
        aerr = obs_0(n)%aerr
        f_cg(i,j,ks1) = f_cg(i,j,ks1) + ov*aerr*obs_0(n)%a00
        f_cg(i,jp,ks1) = f_cg(i,jp,ks1) + ov*aerr*obs_0(n)%a01
        f_cg(ip,jp,ks1) = f_cg(ip,jp,ks1) + ov*aerr*obs_0(n)%a11
        f_cg(ip,j,ks1) = f_cg(ip,j,ks1) + ov*aerr*obs_0(n)%a10
      endif
    enddo
  endif
!
!-----------------------------------------------------------------------
!  Altimeter observations
!-----------------------------------------------------------------------
!
  do n=1,num_obsa
    i = obs_A(n)%io
    ip = i + 1
    j = obs_A(n)%jo
    jp = j + 1
    if (obs_A(n)%code .eq. altm_code .and. obs_A(n)%win) then
      ov = 0.0
      do k=1,kass
        ov = ov + cdnz(k) * &
                ( obs_A(n)%a00 * d_cg(i,j,k) + &
                  obs_A(n)%a01 * d_cg(i,jp,k) + &
                  obs_A(n)%a11 * d_cg(ip,jp,k) + &
                  obs_A(n)%a10 * d_cg(ip,j,k) ) + &
                  cdnzs(k) * &
                ( obs_A(n)%a00 * d_cg(i,j,k+ksalt) + &
                  obs_A(n)%a01 * d_cg(i,jp,k+ksalt) + &
                  obs_A(n)%a11 * d_cg(ip,jp,k+ksalt) + &
                  obs_A(n)%a10 * d_cg(ip,j,k+ksalt) )
      enddo
      aerr = obs_A(n)%aerr
      do k=1,kass
        f_cg(i,j,k) = f_cg(i,j,k) + ov*aerr*obs_A(n)%a00*cdnz(k)
        f_cg(i,jp,k) = f_cg(i,jp,k) + ov*aerr*obs_A(n)%a01*cdnz(k)
        f_cg(ip,jp,k) = f_cg(ip,jp,k) + ov*aerr*obs_A(n)%a11*cdnz(k)
        f_cg(ip,j,k) = f_cg(ip,j,k) + ov*aerr*obs_A(n)%a10*cdnz(k)
        f_cg(i,j,k+ksalt) = f_cg(i,j,k+ksalt) + ov*aerr*obs_A(n)%a00*cdnzs(k)
        f_cg(i,jp,k+ksalt) = f_cg(i,jp,k+ksalt) + ov*aerr*obs_A(n)%a01*cdnzs(k)
        f_cg(ip,jp,k+ksalt) = f_cg(ip,jp,k+ksalt) + ov*aerr*obs_A(n)%a11*cdnzs(k)
        f_cg(ip,j,k+ksalt) = f_cg(ip,j,k+ksalt) + ov*aerr*obs_A(n)%a10*cdnzs(k)
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
  end subroutine comp_f
! </SUBROUTINE> NAME="comp_f"


!!
!!#######################################################################
!! <SUBROUTINE NAME="comp_f_sfc">
!!
!! <DESCRIPTION>
!! Compute the
!! </DESCRIPTION>
!!
subroutine comp_f_sfc (obs_0)
!
  type(ocean_obs0_type), intent(in)        :: obs_0(:)
!
  integer :: i, ip, j, jp, k, ks1, n, pe
  real    :: ov, aerr
integer :: noot, noos
!-----------------------------------------------------------------------
!  data types are encoded in obs%code
!      SST        21 <= code <= 22
!      SSS        23 <= code <= 23
!  these are set in godas_data_mod, and can be modified if needed
!-----------------------------------------------------------------------
!
  pe = mpp_pe()

!-----------------------------------------------------------------------
!   Set f to zero.
!-----------------------------------------------------------------------
!
  f_cg_s = 0.0
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
!  Surface observations
!-----------------------------------------------------------------------
!
noot = 0
noos = 0
  do n=1,num_obs0
    i = obs_0(n)%io
    ip = i + 1
    j = obs_0(n)%jo
    jp = j + 1
    if (obs_0(n)%code .eq. sst_code .and. obs_0(n)%win .and. asm_code .eq. temp_code) then
noot = noot+1
      ov = obs_0(n)%a00 * d_cg_s(i,j) + &
           obs_0(n)%a01 * d_cg_s(i,jp) + &
           obs_0(n)%a11 * d_cg_s(ip,jp) + &
           obs_0(n)%a10 * d_cg_s(ip,j)
      aerr = obs_0(n)%aerr
      f_cg_s(i,j) = f_cg_s(i,j) + ov*aerr*obs_0(n)%a00
      f_cg_s(i,jp) = f_cg_s(i,jp) + ov*aerr*obs_0(n)%a01
      f_cg_s(ip,jp) = f_cg_s(ip,jp) + ov*aerr*obs_0(n)%a11
      f_cg_s(ip,j) = f_cg_s(ip,j) + ov*aerr*obs_0(n)%a10
    else if (obs_0(n)%code .eq. sss_code .and. obs_0(n)%win .and. asm_code .eq. salt_code) then
noos = noos+1
      ov = obs_0(n)%a00 * d_cg_s(i,j) + &
           obs_0(n)%a01 * d_cg_s(i,jp) + &
           obs_0(n)%a11 * d_cg_s(ip,jp) + &
           obs_0(n)%a10 * d_cg_s(ip,j)
      aerr = obs_0(n)%aerr
      f_cg_s(i,j) = f_cg_s(i,j) + ov*aerr*obs_0(n)%a00
      f_cg_s(i,jp) = f_cg_s(i,jp) + ov*aerr*obs_0(n)%a01
      f_cg_s(ip,jp) = f_cg_s(ip,jp) + ov*aerr*obs_0(n)%a11
      f_cg_s(ip,j) = f_cg_s(ip,j) + ov*aerr*obs_0(n)%a10
    endif
  enddo
!
!-----------------------------------------------------------------------
!   Add e to f
!-----------------------------------------------------------------------
!
  f_cg_s = f_cg_s + e_cg_s
!
  call mpp_update_domains (f_cg_s, Dom%domain2d)
!
  end subroutine comp_f_sfc
! </SUBROUTINE> NAME="comp_f_sfc"



!#######################################################################
! <SUBROUTINE NAME="eg_lpsmthr">
!
! <DESCRIPTION>
! This subroutine multiplies g by an approximation to the first guess
! error covariance matrix [e] to get the vector h. The approximation to
! [e] is made by a series of multiplications by 1+laplacian.
! </DESCRIPTION>
!
subroutine eg_lpsmthr

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

end subroutine eg_lpsmthr
! </SUBROUTINE> NAME="eg_lpsmthr"


!#######################################################################
! <SUBROUTINE NAME="eg_lpsmthr_sfc">
!
! <DESCRIPTION>
! This subroutine multiplies g by an approximation to the first guess
! error covariance matrix [e] to get the vector h. The approximation to
! [e] is made by a series of multiplications by 1+laplacian.
! </DESCRIPTION>
!
subroutine eg_lpsmthr_sfc

  integer         :: nit, n, i, j, k, kk, ka, kp, kkp
  integer         :: np, npid2
  integer         :: jbg, jfn, jgbg, jgfn
  real            :: con, col
real :: cs, ev1
integer :: pe

  ni    = Grd%ni
  nj    = Grd%nj
pe = mpp_pe()

  npid2=npits/2

  jbg = jsc
  if (jbg .eq. 1) jbg = 3
  jfn = jec
  if (jfn .ge. jemx) jfn = jemx - 1

!
!-----------------------------------------------------------------------
!   multiply g by the square root of the local background error variance
!-----------------------------------------------------------------------
!

  where (Grd%kmt .eq. 0) g_cg_s = 0.0
  if (asm_code .eq. temp_code) then
    g_cg_s = g_cg_s * sqrt(vtmp_s)
  else if (asm_code .eq. salt_code) then
    g_cg_s = g_cg_s * sqrt(vsal_s)
  endif

  do j=jsd,jed
    do i=isd,ied
      wcn_s(i,j) = 1.0 - wso_s(i,j) - wno_s(i,j) - wea_s(i,j) - wwe_s(i,j)
    enddo
  enddo

  s1(:,:) = g_cg_s(:,:) * wgta_s(:,:)
  s2(:,:) = 0.0

  do np=1,npid2
    do j=jbg,jfn
      do i=isc,iec
        s2(i,j) = ( wcn_s(i,j) * s1(i,j) + wso_s(i,j) * s1(i,j-1) + wno_s(i,j) * s1(i,j+1) &
                        + wwe_s(i,j) * s1(i-1,j) + wea_s(i,j) * s1(i+1,j) ) * Grd%tmask(i,j,1)
      enddo
    enddo
    call mpp_update_domains (s2, Dom%domain2d)
    do j=jbg,jfn
      do i=isc,iec
        s1(i,j) = ( wcn_s(i,j) * s2(i,j) + wno_s(i,j-1) * s2(i,j-1) + wso_s(i,j+1) * s2(i,j+1) &
                      + wwe_s(i,j) * s2(i-1,j) + wea_s(i,j) * s2(i+1,j) ) * Grd%tmask(i,j,1)
      enddo
    enddo
    call mpp_update_domains (s1, Dom%domain2d)
  enddo

  h_cg_s(:,:) = s1(:,:) * wgta_s(:,:)

  where (Grd%kmt .eq. 0) h_cg_s = 0.0

  if (asm_code .eq. temp_code) then
    h_cg_s = h_cg_s * sqrt(vtmp_s)
  else if (asm_code .eq. salt_code) then
    h_cg_s = h_cg_s * sqrt(vsal_s)
  endif

  call mpp_update_domains (h_cg_s, Dom%domain2d)

end subroutine eg_lpsmthr_sfc
! </SUBROUTINE> NAME="eg_lpsmthr_sfc"

end module godas_mod
