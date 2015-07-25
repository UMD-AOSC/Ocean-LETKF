!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_model_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Time step the ocean model using either a twolevel staggered scheme
! (the default) or threelevel leap-frog scheme (the older approach).
! Threelevel scheme remains only for legacy purposes and is not 
! recommended for normal use.   
!</OVERVIEW>
!
!<DESCRIPTION>
! Top level module for ocean model.  Contains routines for 
! initialization, termination, and update of ocean model state.
!
! Design consideration: declarations of top level ocean variables
! are private to this module and hence are only available to other routines
! through argument lists.  For instance, timestep information is passed to
! the various modules on the initialization call and stored internally
! in the respective modules.  This is a crucial design consideration sinces
! it maintains modularity and hence maintainability of the code.  (mjh)
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_model_nml">
!  <DATA NAME="layout" TYPE="integer">
!  Processor domain layout for ocean model. 
!  </DATA> 
!
!  <DATA NAME="io_layout" TYPE="integer, dimension(2)">
!  Processor IO domain layout for ocean model. The default value is (0,0).
!  If either io_layout(1) or (2) is 0, it will default to the number of
!  processors in the computational layout, except restart file will default 
!  to single file if fms_io_nml fileset_write is set to 'single'.  When
!  both entry of io_layout is positive, io_domain will be defined(a pointer in domain2d)
!  and number of distributed files will be layout(1)*layout(2). For example, assume 
!  the restart file is ocean_velocity.res.nc and the diagnostics file is ocean_daily.nc, 
!  if the layout = (1,2), the restart files will be ocean_velocity.res.nc.0000 and 
!  ocean_veloicity.res.nc.0001, the diagnostics files will be ocean_daily.res.nc.0000 
!  and ocean_daily.res.nc.0001. When the io_domain is defined, restart file and 
!  diagnostics file name will be controlled by the io_domain (ignoring fms_io_nml fileset_write). 
!  </DATA> 
!
!   <DATA NAME="dt_ocean"  TYPE="integer"  DEFAULT="-1">
!     Ocean model time step in seconds. 
!   </DATA>
!
!  <DATA NAME="time_tendency" TYPE="character">
!
!  Possible time stepping schemes are the following. 
!
!  1. "threelevel" has the following characteristics
!
!     leap-frog for the time tendency which means the 
!     inviscid/nondissipative processes are at time tau.  
!
!     forward for lateral mixing processes (dissipation at taum1)
!
!     implicit for vertical dissipative (with aidif = 1.0)
!
!     semi-implicit for Coriolis (with acor>0) 
!
!     Because of the need to apply time filters to suppress 
!     leap-frog splitting, the threelevel time stepping scheme
!     does not conserve total tracer content in the model.  
!
!  2. "twolevel" has the following characteristics: 
!
!     staggered 2nd order forward time tendency, which means 
!     that tracer advection, lateral tracer and velocity mixing, 
!     are at time tau. Pressure gradients are at taup1.  
!
!     Adams-Bashforth (either 2nd or 3rd order) for velocity advection  
!     Third order is default as it is more stable.  
!
!     implicit vertical mixing (with aidif = 1.0)
!
!     semi-implicit for Coriolis (with acor > 0) 
!
!     This scheme conserves total volume and tracer in the ocean model.  
!
!  </DATA> 
!
!  <DATA NAME="impose_init_from_restart" TYPE="logical">
!  Consider the following situation:  We have run the model for many years
!  and generated restarts. Time%init is then .false.  Then, we wish to start
!  a series of perturbation experiments from this restart file.  The generic
!  situation is for Time%init to then be .true. However, we need it to be 
!  .false. in mom4 in order to have a proper reading of the full restart
!  information. Setting impose_init_from_restart=.true. will facilitate
!  this setup.  The default is impose_init_from_restart=.false., in which case
!  the model will run through its normal start/stop segments using restarts. 
!  </DATA> 
!
!  <DATA NAME="baroclinic_split" TYPE="integer">
!  baroclinic_split = dtts/dtuv 
!                   = (tracer time step)/(baroclinic time step)
!                   = (ocean model time step)/(baroclinic time step)
!  Transients corrupted if baroclinic_split > 1, so it is recommended
!  to use baroclinic_split=1. 
!  </DATA> 
!  <DATA NAME="barotropic_split" TYPE="integer">
!  Ratio barotropic_split = dtuv/dtbt
!                         = (baroclinic time step)/(barotropic time step). 
!  Must be large enough to resolve the barotropic gravity waves 
!  captured by the barotropic part of the model. 
!  Barotropic waves are dissipated when this splitting 
!  is greater than unity. Model algorithm is not fully 
!  implemented when barotropic_split=1, so user beware
!  if wishing to run an unsplit model simulation. 
!  </DATA> 
!  <DATA NAME="surface_height_split" TYPE="integer">
!  Ratio surface_height_split = dtts/dteta
!                             = (tracer time step)/(surface height time step)
!                             = (tracer time step)/(bottom pressure time step) 
!  Typically this split is set to unity for models where baroclinic_split=1,
!  but something larger when baroclinic_split is order 10.  dteta is the time 
!  step used for update of eta_t or pbot_t. If surface_height_split is 
!  not equal to unity, then tracer conservation properties are compromised. 
!  </DATA> 
!
!  <DATA NAME="reinitialize_thickness" TYPE="logical">
!  When initialized with a nontrivial eta field, it is 
!  necessary to reinitialize the thickness arrays.
!  </DATA> 
!
!  <DATA NAME="cmip_units" TYPE="logical">
!  For CMIP output, we need to have temperature in deg K and 
!  mass transport in kg/s.  The flag cmip_units=.true. will 
!  diagnose CMIP5-related fields with the CMIP units for sending
!  to the diagnostic manager. 
!  Default cmip_units=.false. 
!  </DATA> 
!
!  <DATA NAME="debug" TYPE="logical">
!  For overall model debugging. Set true to print cksums at 
!  each timestep for debugging purposes.
!  </DATA> 
!</NAMELIST>

use fms_mod,                  only: FATAL, NOTE, WARNING, file_exist
use fms_mod,                  only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,                  only: clock_flag_default
use fms_io_mod,               only: set_domain, nullify_domain
use mpp_domains_mod,          only: mpp_global_sum, domain2d, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_domains_mod,          only: mpp_update_domains, BGRID_NE
use mpp_mod,                  only: mpp_error, mpp_pe, mpp_npes, stdlog, stdout
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE, CLOCK_ROUTINE
use stock_constants_mod,      only: ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT
use time_interp_external_mod, only: time_interp_external_init
use time_manager_mod,         only: JULIAN, get_date, get_time
use time_manager_mod,         only: time_type, operator( /= ), operator( < ), operator ( / )
use time_manager_mod,         only: set_time, operator(-), operator( + ), operator( == )
use time_manager_mod,         only: operator(*)

use ocean_advection_velocity_mod, only: ocean_advection_velocity_init, ocean_advection_velocity
use ocean_barotropic_mod,         only: ocean_barotropic_init, ocean_barotropic_end
use ocean_barotropic_mod,         only: update_ocean_barotropic
use ocean_barotropic_mod,         only: ocean_barotropic_forcing, ocean_mass_forcing
use ocean_barotropic_mod,         only: ocean_eta_smooth, ocean_pbot_smooth
use ocean_barotropic_mod,         only: eta_and_pbot_update, eta_and_pbot_diagnose, eta_and_pbot_tendency
use ocean_barotropic_mod,         only: ocean_barotropic_restart
use ocean_bbc_mod,                only: ocean_bbc_init, get_ocean_bbc
use ocean_bih_friction_mod,       only: ocean_bih_friction_init, ocean_bih_friction_end
use ocean_bih_friction_mod,       only: ocean_bih_friction_restart
use ocean_convect_mod,            only: ocean_convect_init
use ocean_coriolis_mod,           only: ocean_coriolis_init
use ocean_density_mod,            only: ocean_density_init, ocean_density_end, update_ocean_density
use ocean_density_mod,            only: ocean_density_restart, ocean_density_diag
use ocean_diagnostics_mod,        only: ocean_diag_init, ocean_diagnostics
use ocean_domains_mod,            only: set_ocean_domain, get_local_indices
use ocean_form_drag_mod,          only: ocean_form_drag_init, compute_visc_form_drag
use ocean_grids_mod,              only: ocean_grids_init, set_ocean_grid_size 
use ocean_grids_mod,              only: set_ocean_hgrid_arrays, set_ocean_vgrid_arrays
use ocean_increment_eta_mod,      only: ocean_increment_eta_init, ocean_increment_eta_source
use ocean_increment_tracer_mod,   only: ocean_increment_tracer_init, ocean_increment_tracer_source
use ocean_increment_velocity_mod, only: ocean_increment_velocity_init, ocean_increment_velocity_source
use ocean_lap_tracer_mod,         only: ocean_lap_tracer_init
use ocean_bih_tracer_mod,         only: ocean_bih_tracer_init
use ocean_lap_friction_mod,       only: ocean_lap_friction_init, ocean_lap_friction_end
use ocean_lap_friction_mod,       only: ocean_lap_friction_restart
use ocean_mixdownslope_mod,       only: ocean_mixdownslope_init, mixdownslope
use ocean_momentum_source_mod,    only: ocean_momentum_source_init
use ocean_nphysics_mod,           only: ocean_nphysics_init, ocean_nphysics_end, neutral_physics
use ocean_nphysics_mod,           only: ocean_nphysics_restart
use ocean_obc_mod,                only: ocean_obc_init, ocean_obc_end
use ocean_obc_mod,                only: ocean_obc_update_boundary, ocean_obc_prepare
use ocean_obc_mod,                only: ocean_obc_restart
use ocean_operators_mod,          only: ocean_operators_init
use ocean_overexchange_mod,       only: ocean_overexchange_init, overexchange
use ocean_overflow_mod,           only: ocean_overflow_init, overflow
use ocean_passive_mod,            only: ocean_passive_tracer_init
use ocean_polar_filter_mod,       only: ocean_polar_filter_init
use ocean_pressure_mod,           only: ocean_pressure_init
use ocean_rivermix_mod,           only: ocean_rivermix_init, rivermix
use ocean_riverspread_mod,        only: ocean_riverspread_init
use ocean_parameters_mod,         only: TWO_LEVEL, THREE_LEVEL
use ocean_parameters_mod,         only: GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, PSIGMA
use ocean_parameters_mod,         only: DEPTH_BASED, PRESSURE_BASED
use ocean_parameters_mod,         only: QUASI_HORIZONTAL, TERRAIN_FOLLOWING
use ocean_parameters_mod,         only: VERTMIX_GOTM
use ocean_parameters_mod,         only: cp_ocean 
use ocean_sbc_mod,                only: ocean_sbc_init, initialize_ocean_sfc, ocean_sfc_end
use ocean_sbc_mod,                only: sum_ocean_sfc, avg_ocean_sfc, get_ocean_sbc, flux_adjust
use ocean_sbc_mod,                only: ocean_sfc_restart
use ocean_shortwave_mod,          only: ocean_shortwave_init, sw_source
use ocean_sigma_transport_mod,    only: ocean_sigma_transport_init, sigma_transport, ocean_sigma_transport_end
use ocean_sigma_transport_mod,    only: ocean_sigma_transport_restart
use ocean_sponges_eta_mod,        only: ocean_sponges_eta_init, sponge_eta_source
use ocean_sponges_tracer_mod,     only: ocean_sponges_tracer_init, sponge_tracer_source
use ocean_sponges_velocity_mod,   only: ocean_sponges_velocity_init, sponge_velocity_source
use ocean_submesoscale_mod,       only: ocean_submesoscale_init, submeso_restrat
use ocean_tempsalt_mod,           only: ocean_tempsalt_ideal_reinit 
use ocean_thickness_mod,          only: ocean_thickness_init, ocean_thickness_init_adjust
use ocean_thickness_mod,          only: ocean_thickness_end, rho_dzt_tendency
use ocean_thickness_mod,          only: update_tcell_thickness, update_ucell_thickness
use ocean_thickness_mod,          only: ocean_thickness_restart
use ocean_time_filter_mod,        only: ocean_time_filter_init, time_filter 
use ocean_topog_mod,              only: ocean_topog_init
use ocean_tracer_advect_mod,      only: ocean_tracer_advect_init, ocean_tracer_advect_end
use ocean_tracer_advect_mod,      only: ocean_tracer_advect_restart
use ocean_tracer_mod,             only: ocean_prog_tracer_init, ocean_diag_tracer_init
use ocean_tracer_mod,             only: update_ocean_tracer, ocean_tracer_end, compute_tmask_limit
use ocean_tracer_mod,             only: ocean_tracer_restart
use ocean_tracer_util_mod,        only: ocean_tracer_util_init
use ocean_tpm_mod,                only: ocean_tpm_source, ocean_tpm_bbc, ocean_tpm_tracer
use ocean_tpm_mod,                only: ocean_tpm_end, ocean_tpm_start
use ocean_tpm_mod,                only: ocean_tpm_flux_init, ocean_tpm_init_sfc
use ocean_types_mod,              only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,              only: ocean_grid_type, ocean_thickness_type
use ocean_types_mod,              only: ocean_domain_type, ice_ocean_boundary_type, ocean_public_type
use ocean_types_mod,              only: ocean_external_mode_type, ocean_adv_vel_type
use ocean_types_mod,              only: ocean_velocity_type, ocean_density_type, ocean_options_type
use ocean_types_mod,              only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,              only: ocean_types_init
use ocean_util_mod,               only: ocean_util_init, write_timestamp
use ocean_velocity_advect_mod,    only: ocean_velocity_advect_init
use ocean_velocity_diag_mod,      only: pressure_conversion, energy_analysis
use ocean_velocity_mod,           only: ocean_velocity_init, update_ocean_velocity, ocean_velocity_end
use ocean_velocity_mod,           only: ocean_explicit_accel_a, ocean_explicit_accel_b, ocean_implicit_accel
use ocean_velocity_mod,           only: ocean_velocity_restart
use ocean_vert_mix_mod,           only: ocean_vert_mix_init, ocean_vert_mix_end, vert_mix_coeff
use ocean_vert_mix_mod,           only: ocean_vert_mix_restart
use ocean_vert_gotm_mod,          only: advect_gotm_compute
use ocean_workspace_mod,          only: ocean_workspace_init, wrk1
use ocean_xlandinsert_mod,        only: ocean_xlandinsert_init, xlandinsert
use ocean_xlandmix_mod,           only: ocean_xlandmix_init, xlandmix
use ocean_drifters_mod,           only : ocean_drifters_init, update_ocean_drifters, ocean_drifters_end

#ifdef ENABLE_GDS    
use godas_types_mod,            only: ocean_cor_tracer_type, ocean_rstr_tracer_type
use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
use godas_obs_mod,              only: godas_obsz_init, godas_obs0_init, godas_obsa_init
use godas_mod,                  only: godas_init, godas_increment, godas_end
use godas_rstr_mod,             only: godas_rstr_init

#endif

implicit none

private


#include <ocean_memory.h>

#ifdef MOM4_STATIC_ARRAYS

  real, dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt        ! diffusion coefficient at base of tracer cells (m^2/sec)  
                                                           ! Note: separate coefficient for temp and salinity 
  real, dimension(isd:ied,jsd:jed,nk)   :: visc_cbu        ! viscosity at base of momentum cells (m^2/sec)
  real, dimension(isd:ied,jsd:jed,nk)   :: gm_diffusivity  ! diffusivity from GM on T-cells (m^2/sec)

  real, dimension(isd:ied,jsd:jed,nk,2) :: visc_cbu_form_drag   ! viscosity (m^2/sec) from Greatbatch form drag
  real, dimension(isd:ied,jsd:jed)      :: pme           ! mass flux per horz area of precip-evap (kg/(s*m^2))
  real, dimension(isd:ied,jsd:jed)      :: melt          ! mass flux per horz area of ice melt water (kg/(s*m^2))
  real, dimension(isd:ied,jsd:jed,2)    :: upme          ! horizontal velocity of precip minus evap (m/s)
  real, dimension(isd:ied,jsd:jed)      :: river         ! mass flux of river (runoff+calving) per horz area (kg/(s*m^2))
  real, dimension(isd:ied,jsd:jed)      :: runoff        ! mass flux of river runoff (liquid) per horz area from (kg/(s*m^2))
  real, dimension(isd:ied,jsd:jed)      :: calving       ! mass flux of calving land ice per horz area from (kg/(s*m^2))
  real, dimension(isd:ied,jsd:jed,2)    :: uriver        ! horizontal velocity from river runoff+calving
  real, dimension(isd:ied,jsd:jed)      :: patm          ! pressure at ocean top from atmosphere and/or ice (Pa) 
  real, dimension(isd:ied,jsd:jed)      :: swflx         ! short wave radiation flux (W/m^2)
  real, dimension(isd:ied,jsd:jed)      :: swflx_vis     ! visible short wave radiation flux (W/m^2)
  real, dimension(isd:ied,jsd:jed,nk)   :: sw_frac_zt    ! short wave radiation flux fraction
  real, dimension(isd:ied,jsd:jed,nk)   :: opacity       ! attenuation of visible light (1/metre)
  real, dimension(isd:ied,jsd:jed)      :: surf_blthick  ! surf boundary layer depth from vertical mixing scheme (m)
  real, dimension(isd:ied,jsd:jed)      :: bott_blthick  ! bottom boundary layer depth from sigma transport (m)
  real, dimension(isd:ied,jsd:jed)      :: rossby_radius ! rossby radius (m)
  real, dimension(isd:ied,jsd:jed,nk)   :: swheat        ! external shortwave heating source W/m^2

#else

  real, pointer, dimension(:,:,:,:) :: diff_cbt        =>NULL() ! diffusion coefficient at base of tracer cells (m^2/sec) 
                                                                ! Note: separate coefficient for temp and salinity 
  real, pointer, dimension(:,:,:)   :: visc_cbu        =>NULL() ! viscosity at base of momentum cells (m^2/sec)
  real, pointer, dimension(:,:,:)   :: gm_diffusivity  =>NULL() ! diffusivity for GM on T-cells (m^2/sec)

  real, pointer, dimension(:,:,:,:) :: visc_cbu_form_drag  =>NULL() ! viscosity (m^2/sec) from Greatbatch form drag
  real, pointer, dimension(:,:)     :: pme        =>NULL()   ! mass flux per horz area from precip-evap (kg/(s*m^2))
  real, pointer, dimension(:,:)     :: melt       =>NULL()   ! mass flux per horz area of ice melt water (kg/(s*m^2))
  real, pointer, dimension(:,:,:)   :: upme       =>NULL()   ! horizontal velocity of precip minus evap (m/s)
  real, pointer, dimension(:,:)     :: river      =>NULL()   ! mass flux of river (runoff+calving) per horz area (kg/(s*m^2)) 
  real, pointer, dimension(:,:)     :: runoff     =>NULL()   ! mass flux of river runoff (liquid) per horz area from (kg/(s*m^2)) 
  real, pointer, dimension(:,:)     :: calving    =>NULL()   ! mass flux of calving land ice per horz area (kg/(s*m^2)) 
  real, pointer, dimension(:,:,:)   :: uriver     =>NULL()   ! horizontal velocity from river (m/s)
  real, pointer, dimension(:,:)     :: patm       =>NULL()   ! pressure at ocean top from sea ice and/or atmosphere (Pa) 
  real, pointer, dimension(:,:)     :: swflx      =>NULL()   ! short wave radiation flux (W/m^2) 
  real, pointer, dimension(:,:)     :: swflx_vis  =>NULL()   ! short wave radiation flux (W/m^2) 
  real, pointer, dimension(:,:,:)   :: sw_frac_zt =>NULL()   ! short wave radiation flux fraction
  real, pointer, dimension(:,:,:)   :: opacity    =>NULL()   ! attenuation of visible light (1/metre)
  real, pointer, dimension(:,:)     :: surf_blthick =>NULL() ! surf boundary layer depth from vertical mixing scheme (m)
  real, pointer, dimension(:,:)     :: bott_blthick =>NULL() ! bottom boundary layer depth from sigma transport (m)
  real, pointer, dimension(:,:)     :: rossby_radius=>NULL() ! rossby radius (m) 
  real, pointer, dimension(:,:,:)   :: swheat    =>NULL()    ! external shortwave heating source W/m^2

#endif

  ! for running Time%init=.true. yet using a restart file 
  logical :: impose_init_from_restart=.false. 

  ! to initialize with a nontrivial eta or pbot (from ocean_barotropic)
  ! we then need to reinitialize the thickness arrays, since they originally
  ! assumed eta_t=0 and pbot_t=pbot0. 
  logical :: reinitialize_thickness=.false.  
  ! for running with an externally provided shortwave heating source
  logical :: ext_swheat_is_set=.false.

  ! time step related variables 
  real :: dtts=0   ! tracer timestep (seconds)
  real :: dtuv=0   ! internal mode timestep (seconds)
  real :: dtbt=0   ! external mode timestep (seconds)
  real :: dteta=0  ! ocean volume time step (eta_t) (seconds)
  real :: dtime_t  ! 2*dtts  for threelevel and dtts  for twolevel
  real :: dtime_u  ! 2*dtuv  for threelevel and dtvu  for twolevel
  real :: dtime_e  ! 2*dteta for threelevel and dteta for twolevel

  ! setting the number of prognostic and diagnostic tracers  
  ! both determined by reading field_table 
  integer :: num_prog_tracers=-1 ! (e.g., temp, salt, age)
  integer :: num_diag_tracers=-1 ! (e.g., frazil, pH) 
#ifdef ENABLE_GDS
  integer :: num_cor_tracers=-1
  integer :: num_obs_z=-1        ! (e.g. total of T(z), S(z))
  integer :: num_obs_0=-1        ! (e.g. total of SST, SSS))
  integer :: num_obs_a=-1        ! (e.g. total of Altimetry))
#endif

  ! number of ocean calls 
  integer :: num_ocean_calls 
  logical :: first_ocn_call=.true. 

  ! for setting model time steps 
  integer :: baroclinic_split    =1  ! ratio of tracer timestep dtts to baroclinic timestep dtuv
  integer :: surface_height_split=1  ! ratio of tracer timestep dtts to the eta_t (or pbot_t) timestep dteta 
  integer :: barotropic_split    =30 ! ratio of baroclinic timestep dtuv to barotropic timestep dtbt

  ! for setting how terms in the equations are time stepped
  ! select time tendency ('threelevel' or 'twolevel")
  ! tendency is an integer corresponding to choice of time tendency 
  character(len=32) :: time_tendency='twolevel' 
  integer :: tendency=0                         

  ! valid vertical coordinate options are the following:
  ! geopotential, zstar, zsigma, pressure, pstar, psigma
  character(len=32) :: vertical_coordinate='geopotential' 

  ! integers corresponding to choice of vertical coordinate
  integer :: vert_coordinate         
  integer :: vert_coordinate_class
  integer :: vert_coordinate_type

  character(len=128) :: version = '$Id: ocean_model.F90,v 16.0.2.3.2.1.2.5.10.4.2.9.6.6.6.1.4.2 2009/12/10 22:20:47 smg Exp $'
  character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

  type(ocean_external_mode_type), save           :: Ext_mode
  type(ocean_adv_vel_type),       save           :: Adv_vel
  type(ocean_density_type),       target, save   :: Dens
  type(ocean_domain_type),        target, save   :: Domain

  type(ocean_grid_type),          target, save   :: Grid
  type(ocean_thickness_type),     target, save   :: Thickness

  type(ocean_time_type),          target, save   :: Time
  type(ocean_time_steps_type),    target, save   :: Time_steps
  type(ocean_options_type),       target, save   :: Ocean_options
  type(ocean_velocity_type),      target, save   :: Velocity

  type(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
  type(ocean_diag_tracer_type), dimension(:), pointer, save :: T_diag =>NULL() 
#ifdef ENABLE_GDS
  type(ocean_cor_tracer_type), dimension(:), pointer, save :: T_cor =>NULL()
  type(ocean_obsz_type), dimension(:), pointer, save :: obs_Z =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_0 =>NULL()
  type(ocean_obs0_type), dimension(:), pointer, save :: obs_A =>NULL()
  type(ocean_rstr_tracer_type), dimension(:), pointer, save :: T_rstr =>NULL()

#endif


  ! identification numbers for mpp clocks
  integer :: id_init
  integer :: id_advect
  integer :: id_vmix
  integer :: id_neutral
  integer :: id_compute_visc_form_drag
  integer :: id_submesoscale
  integer :: id_sw
  integer :: id_sponges_tracer
  integer :: id_sponges_eta
  integer :: id_sponges_velocity
  integer :: id_sigma
  integer :: id_tracer
  integer :: id_bbc
  integer :: id_sbc
  integer :: id_flux_adjust
  integer :: id_explicit_accel_a
  integer :: id_explicit_accel_b
  integer :: id_implicit_accel
  integer :: id_bottom_smooth
  integer :: id_surface_smooth
  integer :: id_eta_and_pbot_tendency
  integer :: id_eta_and_pbot_update
  integer :: id_eta_and_pbot_diagnose
  integer :: id_rho_dzt_tendency
  integer :: id_mass_forcing 
  integer :: id_barotropic_forcing
  integer :: id_barotropic_update
  integer :: id_velocity
  integer :: id_xlandinsert
  integer :: id_xlandmix
  integer :: id_overflow
  integer :: id_overexchange
  integer :: id_mixdownslope
  integer :: id_rivermix
  integer :: id_density
  integer :: id_density_diag
  integer :: id_otpm_source
  integer :: id_otpm_bbc
  integer :: id_otpm_tracer
  integer :: id_diagnostics
  integer :: id_time_filter
  integer :: id_tcell_thickness
  integer :: id_ucell_thickness
  integer :: id_update_halo_tracer
  integer :: id_update_halo_velocity
  integer :: id_godas
  integer :: id_godas_init
  integer :: id_ocean_sfc
  integer :: id_ocean_seg_end
  integer :: id_tmask_limit
  integer :: id_ocean
  integer :: id_advect_gotm
  integer :: id_increment_tracer
  integer :: id_increment_eta
  integer :: id_increment_velocity

  public ocean_model_init
  public ocean_model_end
  public update_ocean_model
  public get_ocean_domain
  public get_ocean_grid_size
  public ocean_public_type
  public ice_ocean_boundary_type
  public ocean_model_init_sfc
  public ocean_model_flux_init
  public ocean_stock_pe
  public ocean_model_restart

  ! routines for interfacing to other models 
  public mom4_get_Tsurf
  public mom4_get_Ssurf
  public mom4_get_UVsurf
  public mom4_get_thickness
  public mom4_get_density
  public mom4_get_prog_tracer
  public mom4_get_temperature_index
  public mom4_get_salinity_index
  public mom4_get_UV
  public mom4_get_dimensions
  public mom4_get_diag_axes
  public mom4_get_num_diag_tracers
  public mom4_get_num_prog_tracers
  public mom4_get_ocean_data
  public mom4_get_surface_tmask 
  public mom4_get_latlon_UV
  public mom4_set_swheat

  public    ocean_model_data_get
  interface ocean_model_data_get
    module procedure ocean_model_data1D_get 
    module procedure ocean_model_data2D_get 
  end interface

  ! for the temperature variable 
  integer :: temp_variable=0 

  ! index for temperature (degC) and salinity (psu) 
  integer :: index_temp=-1    
  integer :: index_salt=-1    

  ! domain layout for parallel processors. for npe=1, layout(2)=(/1,1/)
  integer :: layout(2)=(/1,1/) 

  ! IO domain layout for parallel processors. 
  integer :: io_layout(2)=(/0,0/)   

  ! to print various cksums for debugging purposes
  logical :: debug = .false. 

  logical :: module_is_initialized =.false.
  logical :: have_obc              =.false.   
  logical :: cmip_units            =.false.
 
  type, public ::  ocean_state_type; private
     ! This type is private, and can therefore vary between different ocean models.
     ! All information entire ocean state may be contained here, although it is not
     ! necessary that this is implemented with all models.
     logical       :: is_ocean_pe = .false.       ! .true. on processors that run the ocean model.
  end type ocean_state_type

  integer :: dt_ocean = -1  ! ocean tracer timestep

  namelist /ocean_model_nml/ time_tendency, impose_init_from_restart, reinitialize_thickness,    &
                             baroclinic_split, barotropic_split, surface_height_split,           &
                             layout, io_layout, debug, vertical_coordinate, dt_ocean, cmip_units

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_model_init">
!
! <DESCRIPTION>
! Initialize the ocean model. 
! Arguments: 
!  Ocean (inout)  - A structure containing various publicly visible ocean
!                    surface properties after initialization.
!  OS    (pointer)- A structure whose internal contents are private
!                    to ocean_model_mod that may be used to contain all
!                    information about the ocean's interior state.
!  Time_init (in) - The start time for the coupled model's calendar.
!  Time_in   (in) - The time at which to initialize the ocean model.
! </DESCRIPTION>
!
subroutine ocean_model_init(Ocean, Ocean_state, Time_init, Time_in)   
    type(ocean_public_type), intent(inout)      :: Ocean
    type(ocean_state_type),          pointer    :: Ocean_state
    type(time_type),       intent(in)           :: Time_init
    type(time_type),       intent(in)           :: Time_in 
    
    type(time_type) :: Time_step_ocean
    integer :: secs, days, secs0, days0 
    integer :: n, lay_out(2)
    integer :: ioun, io_status, ierr
 
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if (module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_model_mod (ocean_model_init): module already initialized')
    endif 

    module_is_initialized = .true.

    if (associated(Ocean_state)) then
       call mpp_error(WARNING, "ocean_model_init called with an associated "// &
            "ocean_state_type structure. Model is already initialized.")
       return
    endif
    allocate(Ocean_state)
    Ocean_state%is_ocean_pe = Ocean%is_ocean_pe !This is Not utilized in MOM currently
    
    ! set clock ids
    id_ocean               = mpp_clock_id( 'Ocean', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    id_init                = mpp_clock_id('(Ocean initialization) '         ,grain=CLOCK_SUBCOMPONENT)
    id_godas               = mpp_clock_id('GODAS', flags=clock_flag_default ,grain=CLOCK_COMPONENT)
    id_godas_init          = mpp_clock_id('(GODAS initialization) '         ,grain=CLOCK_SUBCOMPONENT)
    id_advect              = mpp_clock_id('(Ocean advection velocity) '     ,grain=CLOCK_MODULE)
    id_density_diag        = mpp_clock_id('(Ocean density diag) '           ,grain=CLOCK_MODULE)    
    id_density             = mpp_clock_id('(Ocean update density) '         ,grain=CLOCK_MODULE)    
    id_vmix                = mpp_clock_id('(Ocean vertical mixing coeff) '  ,grain=CLOCK_MODULE)
    id_neutral             = mpp_clock_id('(Ocean neutral physics) '        ,grain=CLOCK_MODULE)
    id_submesoscale        = mpp_clock_id('(Ocean submesoscale restrat)'    ,grain=CLOCK_MODULE)
    id_sw                  = mpp_clock_id('(Ocean shortwave) '              ,grain=CLOCK_MODULE)
    id_sponges_eta         = mpp_clock_id('(Ocean sponges_eta) '            ,grain=CLOCK_MODULE)
    id_sponges_tracer      = mpp_clock_id('(Ocean sponges_tracer) '         ,grain=CLOCK_MODULE)
    id_sponges_velocity    = mpp_clock_id('(Ocean sponges_velocity) '       ,grain=CLOCK_MODULE)
    id_xlandinsert         = mpp_clock_id('(Ocean xlandinsert) '            ,grain=CLOCK_MODULE)
    id_xlandmix            = mpp_clock_id('(Ocean xlandmix) '               ,grain=CLOCK_MODULE)
    id_rivermix            = mpp_clock_id('(Ocean rivermix) '               ,grain=CLOCK_MODULE)
    id_overexchange        = mpp_clock_id('(Ocean overexchange) '           ,grain=CLOCK_MODULE)
    id_mixdownslope        = mpp_clock_id('(Ocean mixdownslope) '           ,grain=CLOCK_MODULE)
    id_overflow            = mpp_clock_id('(Ocean overflow) '               ,grain=CLOCK_MODULE)
    id_sigma               = mpp_clock_id('(Ocean sigma transport) '        ,grain=CLOCK_MODULE)
    id_tracer              = mpp_clock_id('(Ocean tracer update) '          ,grain=CLOCK_MODULE)
    id_sbc                 = mpp_clock_id('(Ocean surface flux) '           ,grain=CLOCK_MODULE)
    id_bbc                 = mpp_clock_id('(Ocean bottom flux) '            ,grain=CLOCK_MODULE)
    id_flux_adjust         = mpp_clock_id('(Ocean restoring flux) '         ,grain=CLOCK_MODULE)
    id_otpm_source         = mpp_clock_id('(Ocean TPM source) '             ,grain=CLOCK_MODULE)
    id_otpm_bbc            = mpp_clock_id('(Ocean TPM bbc) '                ,grain=CLOCK_MODULE)
    id_otpm_tracer         = mpp_clock_id('(Ocean TPM tracer) '             ,grain=CLOCK_MODULE)
    id_explicit_accel_a    = mpp_clock_id('(Ocean explicit accel_a) '       ,grain=CLOCK_MODULE)
    id_explicit_accel_b    = mpp_clock_id('(Ocean explicit accel_b) '       ,grain=CLOCK_MODULE)
    id_implicit_accel      = mpp_clock_id('(Ocean implicit accel) '         ,grain=CLOCK_MODULE)
    id_eta_and_pbot_tendency= mpp_clock_id('(Ocean eta and pbot tendency)'  ,grain=CLOCK_MODULE)
    id_eta_and_pbot_update  = mpp_clock_id('(Ocean eta and pbot update)'    ,grain=CLOCK_MODULE)
    id_eta_and_pbot_diagnose= mpp_clock_id('(Ocean eta and pbot diagnose)'  ,grain=CLOCK_MODULE)
    id_rho_dzt_tendency    = mpp_clock_id('(Ocean rho_dzt tendency)'        ,grain=CLOCK_MODULE)
    id_surface_smooth      = mpp_clock_id('(Ocean surface height smooth)'   ,grain=CLOCK_MODULE)
    id_bottom_smooth       = mpp_clock_id('(Ocean bottom pressure smooth)'  ,grain=CLOCK_MODULE)
    id_mass_forcing        = mpp_clock_id('(Ocean mass forcing) '           ,grain=CLOCK_MODULE)
    id_barotropic_forcing  = mpp_clock_id('(Ocean barotropic forcing) '     ,grain=CLOCK_MODULE)
    id_barotropic_update   = mpp_clock_id('(Ocean barotropic dynamics)'     ,grain=CLOCK_MODULE)
    id_velocity            = mpp_clock_id('(Ocean velocity update) '        ,grain=CLOCK_MODULE)
    id_diagnostics         = mpp_clock_id('(Ocean diagnostics)'             ,grain=CLOCK_MODULE)
    id_time_filter         = mpp_clock_id('(Ocean time filter for 3-level)' ,grain=CLOCK_MODULE)
    id_tcell_thickness     = mpp_clock_id('(Ocean update T-cell thickness)' ,grain=CLOCK_MODULE)
    id_ucell_thickness     = mpp_clock_id('(Ocean update U-cell thickness)' ,grain=CLOCK_MODULE)
    id_update_halo_tracer  = mpp_clock_id('(Ocean tracer halo updates)'     ,grain=CLOCK_MODULE)
    id_update_halo_velocity= mpp_clock_id('(Ocean velocity halo update)'    ,grain=CLOCK_MODULE)
    id_ocean_sfc           = mpp_clock_id('(Ocean sum ocean surface)'       ,grain=CLOCK_MODULE)
    id_ocean_seg_end       = mpp_clock_id('(Ocean average state)'           ,grain=CLOCK_MODULE)
    id_tmask_limit         = mpp_clock_id('(Ocean tracer tmask limit)'      ,grain=CLOCK_MODULE)
    id_advect_gotm         = mpp_clock_id('(Ocean gotm: advection)'         ,grain=CLOCK_ROUTINE)
    id_increment_eta       = mpp_clock_id('(Ocean increment eta)'           ,grain=CLOCK_MODULE)
    id_increment_tracer    = mpp_clock_id('(Ocean increment tracer)'        ,grain=CLOCK_MODULE)
    id_increment_velocity  = mpp_clock_id('(Ocean increment velocity)'      ,grain=CLOCK_MODULE)


    call mpp_clock_begin(id_init)
    call write_version_number( version, tagname )

    write(stdoutunit,'(/54x,a/)') '======== STARTING MOM4 INITIALIZATION ========'

#ifdef STATIC_MEMORY
    write(stdoutunit,*) ' '
    write(stdoutunit,*)'==>Error: mom4p0 compiler option "STATIC_MEMORY" is called "MOM4_STATIC_ARRAYS" in mom4p1.'
    write(stdoutunit,*)'          Recompile code with this new name. Apologies for the inconvenience.' 
    call mpp_error(FATAL, &
     '==>Error: Change compiler option "STATIC_MEMORY" to "MOM4_STATIC_ARRAYS" and then recompile.')
#endif 

#ifdef MOM4_STATIC_ARRAYS
    write(stdoutunit,*) ' '
    write(stdoutunit,*)'==>NOTE: Using MOM4_STATIC_ARRAYS cpp option in mom4.'
#else
    write(stdoutunit,*) ' '
    write(stdoutunit,*)'==>NOTE: Using dynamically allocated array option in mom4'
#endif 

    call time_interp_external_init()

    ! provide for namelist over-ride of defaults 
    ioun = open_namelist_file()
    read  (ioun, ocean_model_nml,iostat=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_model_nml)
    write (stdlogunit, ocean_model_nml)
    ierr = check_nml_error(io_status,'ocean_model_nml')
    call close_file (ioun)

    write (stdoutunit,'(/a,i6,a/)') ' ==>Note: Running mom4 using',mpp_npes(),' computer processors.'  

    ! initialize ocean time type information

    Time%calendar  = JULIAN

    if(time_tendency=='threelevel') then 
      tendency        = THREE_LEVEL 
      Time%taum1      = 1
      Time%tau        = 2
      Time%taup1      = 3
      Time%tau_m2     = 1
      Time%tau_m1     = 2
      Time%tau_m0     = 3
      write (stdoutunit,*) ' ' 
      write (stdoutunit,*) '==>Note: Running mom4 with a leap frog to discretize the time tendency.'
      write (stdoutunit,*) '         Unfortunately, this method does not conserve volume and tracer because'
      write (stdoutunit,*) '         it is necessary to use time filtering.  Use the "twolevel" scheme to conserve.'
      write (stdoutunit,*) ' '
    elseif(time_tendency=='twolevel') then  
      tendency        = TWO_LEVEL 
      Time%taum1      = 2
      Time%tau        = 2
      Time%taup1      = 3
      Time%tau_m2     = 1
      Time%tau_m1     = 2
      Time%tau_m0     = 3
      write (stdoutunit,*) ' ' 
      write (stdoutunit,*) &
      '==>Note: Running mom4 with staggered twotime level scheme to compute time tendencies.'
      write (stdoutunit,*) &
      '         This is the default. Mass/volume and tracer are conserved with this scheme.'
      write (stdoutunit,*) ' '
    else
      call mpp_error(FATAL,&
      '==>Error from ocean_model_mod: time_tendency must be "twolevel" or "threelevel".')
    endif 

    if(dt_ocean .lt. 0.0) call mpp_error(FATAL,&
      '==>Error from ocean_model_mod: dt_ocean must be set to a positive integer in ocean_model_nml.')

    Time_step_ocean = set_time (dt_ocean,0)

    Time%init       = (Time_in == Time_init)
    Time%Time_init  = Time_init
    Time%Time_step  = Time_step_ocean
    Time%model_time = Time_in
    Time%itt        = 0

    call get_time(Time_step_ocean, secs, days)
    call get_time(Time_in-Time_init, secs0, days0)
    Time%itt0 = nint((days0+secs0/86400.0)/(days+secs/86400.0))

    write(stdoutunit,'(/a)')&
    ' ==>Note: Time%Time_init = time stamp at very start of the mom4 experiment is given by'
    call write_timestamp(Time%Time_init)

    write(stdoutunit,'(/a)') &
    ' ==>Note: Time%model_time = time stamp at start of this leg of the mom4 experiment is'
    call write_timestamp(Time%model_time)

    if(impose_init_from_restart) then 
        Time%init=.false.
        write(stdoutunit,'(/a)')&
        ' ==>Note: impose_init_from_restart=.true., so will initialize model from a restart file.'
    endif

    if(Time%init) then 
      write(stdoutunit,'(/a)')&
      ' ==>Note: Time%init=.true. =>mom4 will start from user specified initial conditions.' 
    endif 
    if(.not. Time%init) then 
      write(stdoutunit,'(/a)') &
      ' ==>Note: Time%init=.false. =>mom4 will start from restart conditions from previous leg of experiment.' 
    endif 

    dtts    = secs + days*86400
    dtuv    = dtts/baroclinic_split
    dteta   = dtts/surface_height_split
    dtbt    = dtuv/barotropic_split

    if(tendency==THREE_LEVEL) then 
      dtime_t = 2.0*dtts 
      dtime_u = 2.0*dtuv 
      dtime_e = 2.0*dteta
    elseif(tendency==TWO_LEVEL) then  
      dtime_t = dtts 
      dtime_u = dtuv
      dtime_e = dteta
    endif 

    Time_steps%time_tendency = time_tendency 
    Time_steps%tendency      = tendency 
    Time_steps%dtts          = dtts
    Time_steps%dtuv          = dtuv
    Time_steps%dtbt          = dtbt
    Time_steps%dteta         = dteta
    Time_steps%dtime_t       = dtime_t
    Time_steps%dtime_u       = dtime_u
    Time_steps%dtime_e       = dtime_e

    write(stdoutunit,'(/a)')' ==> Note: time steps (seconds) used for mom4' 
    write(stdoutunit,'(a,f10.2)')'  dtts  (tracer)                            = ',dtts
    write(stdoutunit,'(a,f10.2)')'  dtuv  (baroclinic)                        = ',dtuv
    write(stdoutunit,'(a,f10.2)')'  dteta (surface height or bottom pressure) = ',dteta
    write(stdoutunit,'(a,f10.2)')'  dtbt  (barotropic)                        = ',dtbt

    if(baroclinic_split < 1) then 
       call mpp_error(FATAL,&
       '==>Error from ocean_model_mod(ocean_model_init): baroclinic_split must be an integer >= 1')
    endif 
    if(baroclinic_split > 1) then 
       call mpp_error(NOTE,&
       '==>ocean_model_mod: baroclinic_split > 1 corrupts transients & can be unstable. Use with caution.')
       write(stdoutunit,'(/a/)') &
       '==>ocean_model_mod: baroclinic_split > 1 corrupts transients & can be unstable. Use with caution.'
    endif 
    if(surface_height_split < 1) then 
       call mpp_error(FATAL,&
       '==>Error from ocean_model_mod(ocean_model_init): surface_height_split must be an integer >= 1')
    endif 
    if(surface_height_split > 1) then 
       call mpp_error(NOTE, &
       '==>ocean_model_mod: surface_height_split > 1 corrupts transients. Use with caution.')
       write(stdoutunit,'(/a/)') &
       '==>ocean_model_mod: surface_height_split > 1 corrupts transients. Use with caution.'
    endif  

    if(barotropic_split < 1) then 
       call mpp_error(FATAL, &
       '==>Error from ocean_model_mod: barotropic_split must be an integer >= 1')
    endif 
    if(barotropic_split == 1) then 
       call mpp_error(WARNING, &
       '==>Warning from ocean_model_mod: barotropic_split=1 is NOT well tested in mom4p1.')
    endif 

    if (nint(dtuv) /= nint(dtbt)) then
        write (stdoutunit,'(/1x,a/)') &
         '==> Note: The velocity equations will be split into baroclinic and barotropic pieces.'
    endif

    ! vertical coordinate choice: 
    ! if-tests on character strings are slower than  
    ! if-tests on integers. hence it is useful to 
    ! introduce the following integers 
    ! GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, and PSIGMA 
    if(vertical_coordinate=='geopotential') then 
      vert_coordinate       = GEOPOTENTIAL 
      vert_coordinate_class = DEPTH_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note: Using mom4 with geopotential vertical coordinate.'
      write (stdoutunit,'(1x,a/)') &
      '     Beware of vanishing top model grid cells. '
      write (stdoutunit,'(a)'    ) &
      '     The equations are Boussinesq, and so conserve volume rather than mass.'     
      write (stdoutunit,'(a/)'    ) &
      '     Use one of the pressure-like coordinates to get non-Boussinesq effects.'
    elseif(vertical_coordinate=='zstar') then    
      vert_coordinate       = ZSTAR
      vert_coordinate_class = DEPTH_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
      write (stdoutunit,'(/1x,a/)') &
      ' ==> Note: Using mom4 with zstar vertical coordinate.'
      write (stdoutunit,'(a)'    )  &
      '     The equations are Boussinesq, and so conserve volume rather than mass.'     
      write (stdoutunit,'(a/)'    ) &
      '     Use one of the pressure-like coordinates to get non-Boussinesq effects.'
    elseif(vertical_coordinate=='zsigma') then    
      vert_coordinate       = ZSIGMA
      vert_coordinate_class = DEPTH_BASED
      vert_coordinate_type  = TERRAIN_FOLLOWING
      write (stdoutunit,'(/1x,a/)') &
      ' ==> Note: Using mom4 with zsigma terrain following vertical coordinate.'
      write (stdoutunit,'(a)'    )  &
      '     The equations are Boussinesq, and so conserve volume rather than mass.'     
      write (stdoutunit,'(a)'    )  &
      '     Use one of the pressure-like coordinates to get non-Boussinesq effects.'
    elseif(vertical_coordinate=='pressure') then    
      vert_coordinate       = PRESSURE
      vert_coordinate_class = PRESSURE_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note: Using mom4 with pressure vertical coordinate.'
      write (stdoutunit,'(1x,a/)') &
      '     Beware of vanishing bottom and top grid cells.'
      write (stdoutunit,'(a)'    ) &
      '     The equations are non-Boussinesq.'
    elseif(vertical_coordinate=='pstar') then    
      vert_coordinate       = PSTAR
      vert_coordinate_class = PRESSURE_BASED
      vert_coordinate_type  = QUASI_HORIZONTAL
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note: Using mom4 with pstar vertical coordinate.'
      write (stdoutunit,'(a)'    ) &
      '     The equations are non-Boussinesq.'
    elseif(vertical_coordinate=='psigma') then    
      vert_coordinate       = PSIGMA
      vert_coordinate_class = PRESSURE_BASED
      vert_coordinate_type  = TERRAIN_FOLLOWING
      write (stdoutunit,'(/1x,a)') &
      ' ==> Note: Using mom4 with psigma terrain following vertical coordinate.'
      write (stdoutunit,'(a)'    ) &
      '     The equations are non-Boussinesq.'
    else
       call mpp_error(FATAL, &
       '==>Error from ocean_model_mod: no valid vertical coordinate chosen.')
    endif 

    if(vert_coordinate_type==TERRAIN_FOLLOWING) then 
      write (stdoutunit,'(a)'    ) '     '
      write (stdoutunit,'(a)'    ) &
      '     WARNING: TERRAIN_FOLLOWING vertical coordinates are implemented in mom4p1'
      write (stdoutunit,'(a)'    ) &
      '     for process studies and/or regional models. They HAVE NOT been tested for'
      write (stdoutunit,'(a)'    ) &
      '     global climate simulations. The following issues have not been addressed'
      write (stdoutunit,'(a)'    ) &
      '     in the mom4p1 implementation:'
      write (stdoutunit,'(a)'    ) &
      '     (1) Grd%tripolar=.true. has large (wrong) wrho_bt on bottom of T-cell column.'
      write (stdoutunit,'(a)'    ) &
      '     (2) neutral physics is not generalized to terrain following coordinates.'
      write (stdoutunit,'(a)'    ) &
      '     (3) The pressure gradient computation is not sophisticated, and '
      write (stdoutunit,'(a)'    ) &
      '     thus the model likely will suffer from large pressure gradient errors when '
      write (stdoutunit,'(a)'    ) &
      '     using realistic topography and coarse resolutions.'
      write (stdoutunit,'(a)'    ) &
      '     (4) The KPP scheme is unstable with non_local_kpp=.true. and it exhibits a'
      write (stdoutunit,'(a)'    ) &
      '     spurious tracer balance.'
      write (stdoutunit,'(a)'    ) &
      '     The user should thus beware that terrain following coordinates '
      write (stdoutunit,'(a)'    ) &
      '     are still largely in the development stages with mom4p1.'
    endif 

    if(vert_coordinate /= GEOPOTENTIAL .and. tendency==THREE_LEVEL) then 
         write(stdoutunit,'(a)') &
         '==>Error:threelevel tendency only implemented for geopotential vert_coordinate.'
         call mpp_error(FATAL, &
         'ocean_model_init: threelevel tendency only implemented for geopotential')  
    endif 

   
    ! initialize grid and domain information
    call ocean_grids_init(debug, vert_coordinate, vert_coordinate_class)
    call set_ocean_grid_size(Grid, 'INPUT/grid_spec.nc') 
    if(ASSOCIATED(Ocean%maskmap)) then
       lay_out(1) = size(Ocean%maskmap,1)
       lay_out(2) = size(Ocean%maskmap,2)
       call set_ocean_domain(Domain, Grid, layout=lay_out, maskmap=Ocean%maskmap)
    else
       call set_ocean_domain(Domain, Grid, layout=layout, io_layout=io_layout)
    end if
    call set_domain(Domain%domain2d)

    call ocean_workspace_init(Domain, Grid)
    call set_ocean_hgrid_arrays(Domain, Grid)
    call ocean_topog_init(Domain, Grid, 'INPUT/grid_spec.nc', vert_coordinate_type)
    call ocean_obc_init(have_obc, Time, Time_steps, Domain, Grid, Ocean_options, &
          vert_coordinate, debug)
    call set_ocean_vgrid_arrays(Time, Domain, Grid, have_obc)
    call ocean_util_init(Domain)

    ! saved for use in the FMS coupler
    Ocean%axes(:) = Grid%tracer_axes(:) 

#ifndef MOM4_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk

    allocate(gm_diffusivity(isd:ied,jsd:jed,nk))
    allocate(visc_cbu_form_drag(isd:ied,jsd:jed,nk,2))
    allocate(visc_cbu(isd:ied,jsd:jed,nk))
    allocate(diff_cbt(isd:ied,jsd:jed,nk,2))
    allocate(pme(isd:ied,jsd:jed))
    allocate(melt(isd:ied,jsd:jed))
    allocate(upme(isd:ied,jsd:jed,2))
    allocate(river(isd:ied,jsd:jed))
    allocate(runoff(isd:ied,jsd:jed))
    allocate(calving(isd:ied,jsd:jed))
    allocate(uriver(isd:ied,jsd:jed,2))
    allocate(patm(isd:ied,jsd:jed))
    allocate(swflx(isd:ied,jsd:jed))    
    allocate(swflx_vis(isd:ied,jsd:jed))    
    allocate(sw_frac_zt(isd:ied,jsd:jed,nk))    
    allocate(opacity(isd:ied,jsd:jed,nk))    
    allocate(surf_blthick(isd:ied,jsd:jed))    
    allocate(bott_blthick(isd:ied,jsd:jed))    
    allocate(rossby_radius(isd:ied,jsd:jed))    
    allocate(swheat(isd:ied,jsd:jed,nk))
#endif
    visc_cbu_form_drag          = 0.0
    gm_diffusivity              = 0.0
    visc_cbu                    = 0.0
    diff_cbt                    = 0.0
    pme                         = 0.0
    melt                        = 0.0
    upme                        = 0.0
    river                       = 0.0
    runoff                      = 0.0
    calving                     = 0.0
    uriver                      = 0.0
    patm                        = 0.0
    swflx                       = 0.0    
    swflx_vis                   = 0.0    
    sw_frac_zt                  = 0.0    
    opacity                     = 0.0    
    surf_blthick                = 0.0
    bott_blthick                = 0.0
    rossby_radius               = 0.0
    swheat                      = 0.0

    call ocean_types_init()
    call ocean_tracer_util_init(Grid, Domain)
    call ocean_coriolis_init(Grid, Domain, Time, Time_steps)
    call ocean_velocity_init(Grid, Domain, Time, Time_steps, Ocean_options, Velocity, &
                             have_obc, debug )
    call ocean_barotropic_init(Grid, Domain, Time, Time_steps, Ocean_options, Ext_mode, have_obc, &
                               vert_coordinate, vert_coordinate_class, cmip_units, debug)
    call ocean_thickness_init(Time, Time_steps, Domain, Grid, Ext_mode, Thickness,          &
                              vert_coordinate, vert_coordinate_class, vert_coordinate_type, &
                              debug )
    call ocean_operators_init(Grid, Domain, Thickness)

    ! initialize prognostic tracers 
    T_prog => ocean_prog_tracer_init(Grid, Thickness, Ocean, Ocean_options, Domain, Time, Time_steps, &
                                     num_prog_tracers, vert_coordinate_type, have_obc, cmip_units, debug)

    ! initialize diagnostic tracers 
    T_diag => ocean_diag_tracer_init(Time, Thickness, vert_coordinate_type, num_diag_tracers, cmip_units, debug)

    call ocean_advection_velocity_init(Grid, Domain, Time, Time_steps, Thickness, Adv_vel, &
                                       vert_coordinate_class, have_obc, debug)

    call ocean_density_init(Grid, Domain, Time, Time_steps, Thickness, T_prog(:), T_diag(:), Ocean_options, &
                            Dens, vert_coordinate, debug)


    call ocean_thickness_init_adjust(Grid, Time, Dens, Ext_mode, Thickness)

    ! For some experiments, we may wish to initialize eta or pbot
    ! with nontrivial values.  In this case we need to reinitialize the 
    ! thickness arrays, since the original initialization assumed eta=0
    ! and pbot=pbot0.
    if(have_obc .or. reinitialize_thickness) then
      write(stdoutunit,'(a)') &
       '==>Note: Reinitializing thickness to allow for nontrivial initial eta or pbot.'
      call update_tcell_thickness(Time, Grid, Ext_mode, Dens, Thickness) 
      call update_ucell_thickness(Time, Grid, Ext_mode, Thickness) 
    endif 

    call ocean_pressure_init(Grid, Domain, Time, vert_coordinate, vert_coordinate_class,&
                             have_obc, debug)

    call ocean_vert_mix_init(Grid, Domain, Time, Time_steps, Ocean_options, T_prog(:), T_diag(:), &
                             vert_coordinate, vert_coordinate_class, have_obc, debug)
    call ocean_bih_tracer_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtime_t, have_obc)
    call ocean_lap_tracer_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtime_t, have_obc)
    call ocean_sigma_transport_init(Grid, Domain, Time, T_prog(:), Ocean_options, &
                                    dtime_t, vert_coordinate, vert_coordinate_type, debug)
    call ocean_nphysics_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog(:), Ocean_options, &
                             vert_coordinate_type, vert_coordinate_class, debug)
    call ocean_submesoscale_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtime_t, &
                                 vert_coordinate_class, debug)
    call ocean_lap_friction_init(Grid, Domain, Time, Ocean_options, dtime_u, have_obc, debug)
    call ocean_bih_friction_init(Grid, Domain, Time, Ocean_options, dtime_u, have_obc, debug)    
    call ocean_momentum_source_init(Grid, Domain, Time, Ocean_options, debug)    
    call ocean_form_drag_init(Grid, Domain, Time, Time_steps, Ocean_options, debug)    
    call ocean_tracer_advect_init(Grid, Domain, Time, T_prog(:), dtime_t, have_obc, debug)
    call ocean_velocity_advect_init(Grid, Domain, Time, have_obc, debug)
    call ocean_convect_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtime_t)
    call ocean_sbc_init(Grid, Domain, Time, T_prog(:), T_diag(:), Velocity, Ocean, &
                        time_tendency, dtime_t)
    call ocean_bbc_init(Grid, Domain, Time, T_prog(:), Velocity, Ocean_options, vert_coordinate_type)
    call ocean_shortwave_init(Grid, Domain, Time, vert_coordinate, Ocean_options)
    call ocean_sponges_tracer_init(Grid, Domain, Time, T_prog(:), dtime_t, Ocean_options)
    call ocean_sponges_velocity_init(Grid, Domain, Time, Velocity, dtime_u, Ocean_options)
    call ocean_sponges_eta_init(Grid, Domain, Time, Ext_mode, dtime_t, Ocean_options)
    call ocean_xlandmix_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtime_t)    
    call ocean_xlandinsert_init(Grid, Domain, Time, T_prog(:), Ocean_options, dtts)    
    call ocean_riverspread_init(Grid, Domain, Ocean_options)
    call ocean_rivermix_init(Grid, Domain, Time, Time_steps, T_prog(:), Ocean_options, debug)    
    call ocean_overexchange_init(Grid, Domain, Time, T_prog(:), Ocean_options, &
                                 vert_coordinate_type, dtime_t, debug)
    call ocean_mixdownslope_init(Grid, Domain, Time, T_prog(:), Ocean_options, &
                                 vert_coordinate_type, dtime_t, debug)
    call ocean_overflow_init(Grid, Domain, Time, T_prog(:), Ocean_options, &
                             vert_coordinate_type, dtime_t, debug)    
    call ocean_time_filter_init(Grid, Domain, Time, Time_steps, T_prog(:), &
                                vert_coordinate, time_tendency)
    call ocean_polar_filter_init(Grid, Domain, Time, T_prog(:), dtime_t)
    call initialize_ocean_sfc(Time, Thickness, T_prog(:), T_diag(:), Velocity, Ocean)
    call ocean_tpm_start(Domain, Grid, T_prog(:), T_diag(:), Time, Thickness)
    call ocean_diag_init(Grid, Domain, Time, Time_steps, T_prog(:), T_diag(:), Dens, &
                         vert_coordinate_class, have_obc, cmip_units)
    call ocean_increment_eta_init(Grid, Domain, Time)
    call ocean_increment_tracer_init(Grid, Domain, Time, T_prog(:))
    call ocean_increment_velocity_init(Grid, Domain, Time)


#ifdef ENABLE_GDS
! Must call godas_init before godas_obsz_init, godas_obs0_init and godas_obsa_init
! Must call godas_init before godas_rstr_init, it returns without initializing if
!   num_rstr_tracers is zero (num_rstr_tracers is set by namelist in godas_init).
!
    call mpp_clock_begin(id_godas_init)
    T_cor => godas_init(Grid, Domain, Time, T_prog(:), num_cor_tracers, debug)
    obs_Z => godas_obsz_init(Grid, Domain, Time, num_obs_z, debug)
    obs_0 => godas_obs0_init(Grid, Domain, Time, num_obs_0, debug)
    obs_A => godas_obsa_init(Grid, Domain, Time, num_obs_a, debug)
    T_rstr => godas_rstr_init(Grid, Domain, Time, T_prog(:))
    call mpp_clock_end(id_godas_init)

#endif

    call ocean_drifters_init(Domain, Grid, Time, T_prog(:), Velocity, Adv_vel)
    
    ! set index_temp and index_salt for this module 
    do n= 1, num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo

    ! re-initialize tempsalt according to idealized profile 
    call ocean_tempsalt_ideal_reinit(Thickness, index_temp, index_salt, T_prog(:))

    ! initialize passive tracers set according to neutral density or temperature.
    ! also initialize passive tracers to be restored to their initial values. 
    ! these tracers must be initialized after temperature and density are 
    ! initialized.  
    call ocean_passive_tracer_init(T_prog(:), T_diag(:), Dens%neutralrho(:,:,:), &
      Time%init, Time%tau, num_prog_tracers, index_temp, index_salt)

    call nullify_domain()
    call mpp_clock_end(id_init) 

    write(stdoutunit,'(/52x,a/)') '======== COMPLETED MOM4 INITIALIZATION ========'

  end subroutine ocean_model_init
! </SUBROUTINE> NAME="ocean_model_init"


 !#######################################################################
! <SUBROUTINE NAME="update_ocean_model">
!
! <DESCRIPTION>
! Update in time the ocean model fields. 
!   This subroutine uses the forcing in Ice_ocean_boundary to advance the
! ocean model's state from the input value of Ocean_state (which must be for
! time time_start_update) for a time interval of Ocean_coupling_time_step,
! returning the publicly visible ocean surface properties in Ocean_sfc and
! storing the new ocean properties in Ocean_state.
!
! Arguments: 
!  Ice_ocean_boundary - A structure containing the various forcing
!                                 fields coming from the ice. It is intent in.
!  Ocean_state - A structure containing the internal ocean state.
!  Ocean_sfc - A structure containing all the publicly visible ocean
!                        surface fields after a coupling time step.
!  time_start_update - The time at the beginning of the update step.
!  Ocean_coupling_time_step - The amount of time over which to advance
!                                       the ocean.

! Note: although several types are declared intent(inout), this is to allow for
!   the possibility of halo updates and to keep previously allocated memory.
!   In practice, Ice_ocean_boundary is intent in, Ocean_state is private to
!   this module and intent inout, and Ocean_sfc is intent out.
! </DESCRIPTION>
!
  subroutine update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, &
       time_start_update, Ocean_coupling_time_step)
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
    type(ocean_state_type),        pointer       :: Ocean_state
    type(ocean_public_type),       intent(inout) :: Ocean_sfc
    type(time_type), intent(in)                  :: time_start_update
    type(time_type), intent(in)                  :: Ocean_coupling_time_step
    
    integer :: num_ocn

    integer :: taum1, tau, taup1
    integer :: i, j, k, n

    call mpp_clock_begin(id_ocean)

    if(first_ocn_call) then 
      num_ocean_calls =  Ocean_coupling_time_step / Time%Time_step 
      if ( num_ocean_calls * Time%Time_step /= Ocean_coupling_time_step ) then 
       call mpp_error(FATAL, &
       '==>Error from ocean_model_mod(update_ocean_model): cpld time step is not a multiple of the ocean time step', FATAL)
      endif
      first_ocn_call=.false.
    endif 


    ! Loop over num_ocean_calls, moved here from the coupler due to interface changes
    do num_ocn = 1,num_ocean_calls

       ! increment ocean time and time labels 
       Time%model_time = Time%model_time + Time%Time_step
       Time%itt        = Time%itt+1
       Time%itt0       = Time%itt0+1

       Time%taum1      = mod(Time%itt+0,3)+1
       Time%tau        = mod(Time%itt+1,3)+1
       Time%taup1      = mod(Time%itt+2,3)+1

       Time%tau_m2     = mod(Time%itt+0,3)+1
       Time%tau_m1     = mod(Time%itt+1,3)+1
       Time%tau_m0     = mod(Time%itt+2,3)+1

       if(tendency==TWO_LEVEL) then  
          Time%taum1 = Time%tau
       endif

       taum1 = Time%taum1
       tau   = Time%tau
       taup1 = Time%taup1


       ! initialize some fields to zero at start of the time step

       ! tracer tendency th_tendency     [units rho_dzt*tracer concentration per time]
       ! tracer source has no weighting  [units tracer concentration per time] 
       do n=1,num_prog_tracers
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   T_prog(n)%th_tendency(i,j,k) = 0.0
                   T_prog(n)%source(i,j,k)      = 0.0
                enddo
             enddo
          enddo
       enddo

       ! sources of mass and form drag vertical viscosity and gm_diffusivity 
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Thickness%mass_source(i,j,k) = 0.0
                visc_cbu_form_drag(i,j,k,1)  = 0.0
                visc_cbu_form_drag(i,j,k,2)  = 0.0
                gm_diffusivity(i,j,k)        = 0.0
             enddo
          enddo
       enddo

       ! boundary layer thicknesses
       do j=jsd,jed
          do i=isd,ied
             surf_blthick(i,j) = 0.0
             bott_blthick(i,j) = 0.0
          enddo
       enddo


       ! calculate tracer tmask_limit based on tracer values at time tau
       call mpp_clock_begin(id_tmask_limit)
       call compute_tmask_limit(Time, T_prog(1:num_prog_tracers))
       call mpp_clock_end(id_tmask_limit)

       ! obtain surface boundary fluxes using boundary field from coupler
       call mpp_clock_begin(id_sbc)
       call get_ocean_sbc(Time, Ice_ocean_boundary, Thickness, Ext_mode, T_prog(1:num_prog_tracers), &
            Velocity, pme, melt, river, runoff, calving, upme, uriver, swflx, swflx_vis, patm)
       call mpp_clock_end(id_sbc)

       ! compute "flux adjustments" (e.g., surface tracer restoring, flux correction)
       call mpp_clock_begin(id_flux_adjust)
       call flux_adjust(Time, T_diag(1:num_diag_tracers), T_prog(1:num_prog_tracers), Velocity, river, melt, pme)
       call mpp_clock_end(id_flux_adjust)

       ! calculate bottom momentum fluxes and bottom tracer fluxes
       call mpp_clock_begin(id_bbc)
       call get_ocean_bbc(Time, Thickness, Velocity, T_prog(1:num_prog_tracers))
       call mpp_clock_end(id_bbc)

       ! add shortwave heating to T_prog%th_tendency 
       call mpp_clock_begin(id_sw)
       if(ext_swheat_is_set) then
          ! apply external sw heating source if it is set
          do k=1,nk-1
             do j=jsc,jec
                do i=isc,iec
                   T_prog(index_temp)%th_tendency(i,j,k) = &
                        T_prog(index_temp)%th_tendency(i,j,k) + swheat(i,j,k) / cp_ocean
                enddo
             enddo
          enddo
       else
        ! compute the sw penetration
        call sw_source(Time, Thickness, Dens, T_diag(:), swflx, swflx_vis, T_prog(index_temp), sw_frac_zt, opacity)
       endif
       call mpp_clock_end(id_sw)

       ! compute vertical mixing coefficients. 
       ! if use kpp, also add nonlocal tendency to T_prog%th_tendency 
       call mpp_clock_begin(id_vmix)    
       call vert_mix_coeff(Time, Thickness, Velocity, T_prog(1:num_prog_tracers),    &
            T_diag(1:num_diag_tracers), Dens, swflx, sw_frac_zt, pme, &
            river, visc_cbu, diff_cbt, surf_blthick)
       call mpp_clock_end(id_vmix)

       ! calculate some diagnostic arrays for ocean density
       call mpp_clock_begin(id_density_diag) 
       call ocean_density_diag(Time, T_prog(index_temp), &
            T_prog(index_salt), Dens, T_diag(:))   
       call mpp_clock_end(id_density_diag) 

       ! set the ocean source terms for the tracer packages
       call mpp_clock_begin(id_otpm_source)
       call ocean_tpm_source(isd, ied, jsd, jed, Domain, Grid, T_prog(:), T_diag(:),   &
            Time, Thickness, Dens, opacity, surf_blthick, dtts)
       call mpp_clock_end(id_otpm_source)

       ! set ocean surface boundary conditions for the tracer packages
       call mpp_clock_begin(id_otpm_bbc)
       call ocean_tpm_bbc(Domain, Grid, T_prog(:))
       call mpp_clock_end(id_otpm_bbc)

       ! add sponges to T_prog%th_tendency 
       call mpp_clock_begin(id_sponges_tracer)
       call sponge_tracer_source(Time, Thickness, T_prog(1:num_prog_tracers)) 
       call mpp_clock_end(id_sponges_tracer)

       ! add increments to T_prog%th_tendency as used in Australian OFAM/Bluelink
       call mpp_clock_begin(id_increment_tracer)
       call ocean_increment_tracer_source(Time, Thickness, T_prog(1:num_prog_tracers)) 
       call mpp_clock_end(id_increment_tracer)

       ! add cross land mixing to T_prog%th_tendency and Thickness%mass_source  
       call mpp_clock_begin(id_xlandmix)
       call xlandmix (Time, Ext_mode, Dens, Thickness, T_prog(1:num_prog_tracers))
       call mpp_clock_end(id_xlandmix)

       ! add cross land insertion to T_prog%th_tendency and Thickness%mass_source
       call mpp_clock_begin(id_xlandinsert)
       call xlandinsert (Time, Ext_mode, Dens, Thickness, T_prog(1:num_prog_tracers))
       call mpp_clock_end(id_xlandinsert)

       ! add river discharge to T_prog%th_tendency and/or enhance diff_cbt next to river mouths 
       call mpp_clock_begin(id_rivermix)
       call rivermix (Time, Thickness, Dens, T_prog(1:num_prog_tracers), river, runoff, calving, &
                      diff_cbt, index_temp, index_salt)
       call mpp_clock_end(id_rivermix)

       ! add discharge of dense shelf water into abyss to T_prog%th_tendency
       call mpp_clock_begin(id_overflow)
       call overflow (Time, Thickness, T_prog(1:num_prog_tracers), Dens, index_temp, index_salt)
       call mpp_clock_end(id_overflow)

       ! add exchange of dense shelf water properties into abyss to T_prog%th_tendency
       call mpp_clock_begin(id_overexchange)
       call overexchange (Time, Thickness, T_prog(1:num_prog_tracers), Dens, index_temp, index_salt)
       call mpp_clock_end(id_overexchange)

       ! add mixing of dense shelf water properties into abyss to T_prog%th_tendency
       call mpp_clock_begin(id_mixdownslope)
       call mixdownslope (Time, Thickness, T_prog(1:num_prog_tracers), Dens, index_temp, index_salt)
       call mpp_clock_end(id_mixdownslope)

       ! add surface height smoother to T_prog%surface_smooth and Thickness%mass_source
       if(vert_coordinate_class==DEPTH_BASED) then 
          call mpp_clock_begin(id_surface_smooth)
          call ocean_eta_smooth(Time, Thickness, Ext_mode, T_prog(1:num_prog_tracers))
          call mpp_clock_end(id_surface_smooth)
       endif

       ! add bottom pressure smoother to T_prog%bottom_smooth and Thickness%mass_source
       if(vert_coordinate_class==PRESSURE_BASED) then 
          call mpp_clock_begin(id_bottom_smooth)
          call ocean_pbot_smooth(Time, Thickness, Ext_mode, T_prog(1:num_prog_tracers))
          call mpp_clock_end(id_bottom_smooth)
       endif
       ! get prescribed OBC data from files
       if (have_obc) call ocean_obc_prepare(Time, Thickness, Ext_mode, T_prog(1:num_prog_tracers))

       ! vertical integral of mass forcing used for eta and pbot update
       call mpp_clock_begin(id_mass_forcing)
       call ocean_mass_forcing(Time, Thickness, Ext_mode) 
       call mpp_clock_end(id_mass_forcing)

       ! As used for the Australian OFAM/Bluelink eta forcing
       call mpp_clock_begin(id_increment_eta)
       call ocean_increment_eta_source(Time, Ext_mode) 
       call mpp_clock_end(id_increment_eta)

       ! add eta sponges
       call mpp_clock_begin(id_sponges_eta)
       call sponge_eta_source(Time, Ext_mode)
       call mpp_clock_end(id_sponges_eta)

       ! accumulate terms for surface height tendency or bottom pressure tendency 
       call mpp_clock_begin(id_eta_and_pbot_tendency)
       call eta_and_pbot_tendency(Time, pme, river, Ext_mode)
       call mpp_clock_end(id_eta_and_pbot_tendency)

       ! compute rho_dzt tendency (a function of eta_and_pbot_tendency)
       call mpp_clock_begin(id_rho_dzt_tendency)
       call rho_dzt_tendency(Time, Grid, Ext_mode, Thickness)
       call mpp_clock_end(id_rho_dzt_tendency)

       ! compute advective velocity components on faces of T-cells and U-cells.
       ! included are Thickness%mass_source and Thickness%rho_dzt_tendency
       call mpp_clock_begin(id_advect)
       call ocean_advection_velocity(Velocity, Time, Thickness, Dens, pme, river, Adv_vel)
       call mpp_clock_end(id_advect)

       ! update ocean free surface height or bottom pressure using "big time step"
       call mpp_clock_begin(id_eta_and_pbot_update)
       call eta_and_pbot_update(Time, Ext_mode)
       call mpp_clock_end(id_eta_and_pbot_update)

       ! update taup1 value of the tracer grid cell thickness.  note that after 
       ! this update, dzt and dzwt (and other Thickness%T-arrays with no 
       ! time indices) are now at time taup1.
       call mpp_clock_begin(id_tcell_thickness)
       call update_tcell_thickness(Time, Grid, Ext_mode, Dens, Thickness) 
       call mpp_clock_end(id_tcell_thickness)

       ! advect tke and diss for GOTM scheme in preparation for the next time step 
       if(Ocean_options%vertmix==VERTMIX_GOTM) then 
          call mpp_clock_begin(id_advect_gotm) 
          call advect_gotm_compute(Time, Adv_vel, Thickness, pme, river)
          call mpp_clock_end(id_advect_gotm) 
       endif

       ! add tendency from sigma transport to Tracer%th_tendency
       call mpp_clock_begin(id_sigma)
       call sigma_transport(Time, Thickness, Dens, T_prog(1:num_prog_tracers), Adv_vel, bott_blthick) 
       call mpp_clock_end(id_sigma)

       ! add tendency from neutral physics to T_prog%th_tendency.
       ! also compute T_prog%K33_implicit for use in time-implicit update.
       ! note that surf_blthick is needed from vert_mix_coeff (kpp) to 
       ! define "neutral physics surface boundary layer". 
       ! also, bott_blthick is needed from sigma_transport to 
       ! define "neutral physics bottom boundary layer". 
       call mpp_clock_begin(id_neutral)
       call neutral_physics(Time, Thickness, Dens, Dens%rho(:,:,:,taum1), &
         T_prog(1:num_prog_tracers), gm_diffusivity, surf_blthick, bott_blthick, rossby_radius)
       call mpp_clock_end(id_neutral)

       ! compute form drag vertical viscosity for vertical momentum transport
       call mpp_clock_begin(id_compute_visc_form_drag)
       call compute_visc_form_drag(Time, Thickness, Velocity, Dens, &
            gm_diffusivity, surf_blthick, visc_cbu_form_drag)
       call mpp_clock_end(id_compute_visc_form_drag)

       ! add tendency from submesoscale param to T_prog%th_tendency.
       call mpp_clock_begin(id_submesoscale)
       call submeso_restrat(Time, Thickness, Dens, T_prog, surf_blthick)
       call mpp_clock_end(id_submesoscale)

       ! update to time=taup1 the value of tracer concentrations
       call mpp_clock_begin(id_tracer)
       call update_ocean_tracer(Time, Dens, Adv_vel, Thickness, pme, diff_cbt, Dens%pressure_at_depth, &
            T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers))
       call mpp_clock_end(id_tracer)

       ! diagnose time=taup1 ocean free surface height or bottom pressure.
       ! also diagnose geodepth_zt.
       call mpp_clock_begin(id_eta_and_pbot_diagnose)
       call eta_and_pbot_diagnose(Time, Dens, Thickness, patm, pme, river, Ext_mode)
       call mpp_clock_end(id_eta_and_pbot_diagnose)

       ! perform extra calculations for the ocean tracer packages
       call mpp_clock_begin(id_otpm_tracer)
       call ocean_tpm_tracer(Domain, T_prog(:), T_diag(:), Grid, Time, Thickness, Dens, dtts, &
            surf_blthick, swflx_vis, opacity, diff_cbt, river)
       call mpp_clock_end(id_otpm_tracer)

       ! fill processor halos for tracers(taup1)
       ! halo values are needed prior to rho(taup1) computation 
       call mpp_clock_begin(id_update_halo_tracer)
       do n=1,num_prog_tracers
          call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), &
               Domain%domain2d, complete=T_prog(n)%complete)
       enddo
       do n=1,num_prog_tracers 
          if(have_obc)call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
       enddo
       call mpp_clock_end(id_update_halo_tracer)

       ! update to time=taup1 Dens%pressure_at_depth and Dens%rho
       call mpp_clock_begin(id_density) 
       call update_ocean_density(Time, Thickness, T_prog(index_temp), &
            T_prog(index_salt), Ext_mode, Dens)   
       call mpp_clock_end(id_density) 

       ! add time explicit contributions to Velocity%accel.  
       ! compute just those pieces needed to force barotropic dynamics 
       Velocity%rossby_radius(:,:) = Grid%umask(:,:,1)*rossby_radius(:,:)
       call mpp_clock_begin(id_explicit_accel_a)
       if(tendency == TWO_LEVEL) then 
          call ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, Dens,  &
               Dens%rho(:,:,:,taup1), pme, river, upme, uriver)
       elseif(tendency == THREE_LEVEL) then 
          call ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, Dens,  &
               Dens%rho(:,:,:,tau), pme, river, upme, uriver)
       endif
       call mpp_clock_end(id_explicit_accel_a)

    ! vertical integral of forcing used for barotropic dynamics 
    call mpp_clock_begin(id_barotropic_forcing)
    call ocean_barotropic_forcing(Time, Thickness, Velocity, Ext_mode) 
    call mpp_clock_end(id_barotropic_forcing)

    ! update (udrho,vdrho) and eta_t_bar or pbot_t_bar using barotropic timesteps 
    call mpp_clock_begin(id_barotropic_update)
    call update_ocean_barotropic (Time, Dens, Thickness, Velocity, Adv_vel, &
                                  T_prog(index_temp), T_prog(index_salt),   &
                                  Ext_mode, patm, pme, river)
    call mpp_clock_end(id_barotropic_update)

    ! remaining time explicit contributions to rho*dz*acceleration
    call mpp_clock_begin(id_explicit_accel_b)
    call ocean_explicit_accel_b(visc_cbu, Time, Thickness, Velocity)
    call mpp_clock_end(id_explicit_accel_b)

    ! compute pressure conversion diagnostic prior to updating ucell thickness 
    if(tendency == TWO_LEVEL) then 
        call pressure_conversion (Time, Thickness, Dens%rho(:,:,:,taup1), &
                                  Dens, Ext_mode, Adv_vel, Velocity, vert_coordinate_class)
    elseif(tendency == THREE_LEVEL) then 
        call pressure_conversion (Time, Thickness, Dens%rho(:,:,:,tau),   &
                                  Dens, Ext_mode, Adv_vel, Velocity, vert_coordinate_class)
    endif

    ! update time=taup1 value of the velocity grid cell thickness
    ! Thickness%U-arrays with no time index are now at taup1
    call mpp_clock_begin(id_ucell_thickness)
    call update_ucell_thickness(Time, Grid, Ext_mode, Thickness) 
    call mpp_clock_end(id_ucell_thickness)

    ! As used for the Australian OFAM/Bluelink increment to velocity
    call mpp_clock_begin(id_increment_velocity)
    call ocean_increment_velocity_source(Time, Thickness, Velocity)
    call mpp_clock_end(id_increment_velocity)

    ! add velocity sponges
    call mpp_clock_begin(id_sponges_velocity)
    call sponge_velocity_source(Time, Thickness, Velocity)
    call mpp_clock_end(id_sponges_velocity)

    ! density and thickness weighted acceleration from
    ! implicit vertical friction and implicit coriolis 
    call mpp_clock_begin(id_implicit_accel)
    call ocean_implicit_accel(visc_cbu, visc_cbu_form_drag, Time, Thickness, Velocity) 
    call mpp_clock_end(id_implicit_accel)

    ! update to time=taup1 the ocean velocity
    call mpp_clock_begin(id_velocity)
    call update_ocean_velocity(Time, Thickness, barotropic_split, vert_coordinate_class, &
                               Ext_mode, Velocity) 
    call mpp_clock_end(id_velocity)

    ! compute energy analysis diagnostic 
    call energy_analysis (Time, Thickness, Ext_mode, Adv_vel, Dens,               &
                          pme, river, upme, uriver, visc_cbu, visc_cbu_form_drag, &
                          Velocity)

    ! perform some numerical diagnostics (e.g., tracer and mass conservation, CFL checks, etc.)
    call mpp_clock_begin(id_diagnostics)
    call ocean_diagnostics(Time, Thickness, T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers), &
                           Adv_vel, Ext_mode, Dens, Velocity, pme, melt, runoff, calving, visc_cbu)
    call mpp_clock_end(id_diagnostics)

    ! fill halo values for the velocity field 
    call mpp_clock_begin(id_update_halo_velocity)
    call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1), &
                            Domain%domain2d,gridtype=BGRID_NE)
    if(have_obc) then
       call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','n')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','t')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','t')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','n')
    endif
    call mpp_clock_end(id_update_halo_velocity)
    
    ! apply time filter when using leap-frog time tendency 
    call mpp_clock_begin(id_time_filter)
    if(tendency==THREE_LEVEL .and. vert_coordinate == GEOPOTENTIAL) then 
      call time_filter(Time, Grid, Thickness, T_prog(1:num_prog_tracers), Velocity, Ext_mode)
    endif 
    call mpp_clock_end(id_time_filter)

    ! modifications to prognostic variables using ocean data assimilation 
#ifdef ENABLE_GDS


    call mpp_clock_begin(id_godas)
    call godas_increment(Time, T_prog(:), Ext_mode, T_cor(:), obs_Z(:), obs_0(:), obs_A(:), T_rstr(:))
    call mpp_clock_end(id_godas)
#endif

    call update_ocean_drifters(Velocity, Adv_vel, T_prog(:), Grid, Time)
    
    ! sum ocean sfc state over coupling interval
    call mpp_clock_begin(id_ocean_sfc)
    call sum_ocean_sfc(Time, Thickness, T_prog(1:num_prog_tracers), &
                       T_diag(1:num_diag_tracers), Dens, Velocity, Ocean_sfc)
    call mpp_clock_end(id_ocean_sfc)

    enddo ! {end of do no=1,num_ocean_calls}

    ! at end of coupling interval, pass averaged ocean state to other component models 
    !
       call mpp_clock_begin(id_ocean_seg_end)
       call avg_ocean_sfc(Time, Thickness, T_prog(1:num_prog_tracers), &
                          T_diag(1:num_diag_tracers), Velocity, Ocean_sfc)
       call mpp_clock_end(id_ocean_seg_end)

    call mpp_clock_end(id_ocean)

    return

  end subroutine update_ocean_model
! </SUBROUTINE> NAME="update_ocean_model"



!#######################################################################
! <SUBROUTINE NAME="get_ocean_grid_size">
!
! <DESCRIPTION>
! Obtain the ocean grid size. 
! </DESCRIPTION>
!
  subroutine get_ocean_grid_size(num_lon, num_lat, num_z)

    integer,           intent(out) :: num_lon, num_lat
    integer, optional, intent(out) :: num_z

    num_lon = Grid%ni
    num_lat = Grid%nj
    if (PRESENT(num_z)) num_z   = Grid%nk

    return
    
  end subroutine get_ocean_grid_size
! </SUBROUTINE> NAME="get_ocean_grid_size"

subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc

  !Output array2D is allocated on the coupled model compute domain
  !Grid% and T_prog% arrays below are allocated on
  !data   domain for dynamic case
  !or
  !global domain for static  case
  !Hence we need to start copying them from Domain%isc and Domain%jsc

  select case(name)
  case('area')
     array2D(isc:,jsc:) = Ocean%area(isc:,jsc:)
  case('t_surf')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)
  case('mask')
     array2D(isc:,jsc:) = Grid%tmask(Domain%isc:,Domain%jsc:,1)
  case('t_pme')
     array2D(isc:,jsc:) = T_prog(index_temp)%tpme(Domain%isc:,Domain%jsc:)
  case('t_runoff')
     array2D(isc:,jsc:) = T_prog(index_temp)%trunoff(Domain%isc:,Domain%jsc:)
  case('t_calving')
     array2D(isc:,jsc:) = T_prog(index_temp)%tcalving(Domain%isc:,Domain%jsc:)
  case('btfHeat')
     array2D(isc:,jsc:) = T_prog(index_temp)%conversion * T_prog(index_temp)%btf(Domain%isc:,Domain%jsc:)
  case default
     call mpp_error(FATAL,'get_ocean_grid_data2D: unknown argument name='//name)
  end select
  

end subroutine ocean_model_data2D_get

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value
  
  select case(name)
  case('c_p')
     value = cp_ocean
  case default
     call mpp_error(FATAL,'get_ocean_grid_data1D: unknown argument name='//name)
  end select
  

end subroutine ocean_model_data1D_get

!#######################################################################
! <SUBROUTINE NAME="get_ocean_domain">
!
! <DESCRIPTION>
! Obtain the ocean domain size. 
! </DESCRIPTION>
!
  subroutine get_ocean_domain(Ocean_domain)

    type(domain2d), intent(out) :: Ocean_domain

    if (.NOT. module_is_initialized) then
       call mpp_error(FATAL, &
       '==>Error in ocean_model_mod (get_ocean_domain): module is not initialized')
    endif

    Ocean_domain = Domain%domain2d

    return

  end subroutine get_ocean_domain
! </SUBROUTINE> NAME="get_ocean_domain"


!#######################################################################
! <SUBROUTINE NAME="ocean_model_init_sfc">
!
! <DESCRIPTION>
! Call ocean_tpm_init_sfc and pass it the needed arguments, most of which
! are local to the ocean model.
! </DESCRIPTION>
!
  subroutine ocean_model_init_sfc(Ocean_state, Ocean)
    type(ocean_state_type),  pointer       :: Ocean_state    
    type(ocean_public_type), intent(inout) :: Ocean

    call ocean_tpm_init_sfc(Domain, T_prog(:), Dens, Ocean, Time, Grid)
    
    return

  end subroutine ocean_model_init_sfc
! </SUBROUTINE> NAME="ocean_model_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_model_flux_init">
!
! <DESCRIPTION>
! Call ocean_tpm_flux_init and pass it the needed arguments, most of which
! are local to the ocean model. 
!
! Currently, no arguments are passed.
! </DESCRIPTION>
!
  subroutine ocean_model_flux_init(Ocean_state)
    type(ocean_state_type), pointer       :: Ocean_state
    
    call ocean_tpm_flux_init
    
    return

  end subroutine ocean_model_flux_init
! </SUBROUTINE> NAME="ocean_model_flux_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_model_end">
!
! <DESCRIPTION>
! Close down the ocean model 
!   This subroutine terminates the model run, saving the ocean state in a
! restart file and deallocating any data associated with the ocean.
!
! NOTE from nnz: This module keeps its own Time and does not need the Time_in argument.
! Arguments: 
!   Ocean_state (type(ocean_state_type), pointer) - A structure containing the internal ocean state.
!   Time_in     (type(time_type), intent(in))     - The model time, used for writing restarts.
!   Ocean_sfc   (type(ocean_public_type), optional, intent(inout))- An ocean_public_type structure that is to be
!                   deallocated upon termination.
! </DESCRIPTION>
!
  subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time_in)
    type(ocean_state_type),            pointer       :: Ocean_state
    type(time_type),                   intent(in)    :: Time_in
    type(ocean_public_type), optional, intent(inout) :: Ocean_sfc

  integer :: stdoutunit 
  stdoutunit=stdout() 

#ifdef ENABLE_GDS
    call godas_end(Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)
#endif
    call ocean_tracer_end(Time, T_prog(:), T_diag(:))
    call ocean_tracer_advect_end(Time, T_prog(:))
    call ocean_nphysics_end(Time)
    call ocean_bih_friction_end(Time)
    call ocean_lap_friction_end(Time)
    call ocean_sigma_transport_end(Time)
    call ocean_tpm_end(Domain, Grid, T_prog(:), T_diag(:), Time, Thickness)
    call ocean_velocity_end(Time, Velocity)
    call ocean_barotropic_end(Time, Ext_mode)
    call ocean_thickness_end(Time, Grid, Thickness)
    call ocean_density_end(Time, Dens)
    if(have_obc) call ocean_obc_end(Time, have_obc)
    call ocean_sfc_end(Ocean_sfc)
    call ocean_vert_mix_end(Time)
    call ocean_drifters_end(Grid)
    
    write (stdoutunit,'(//,1x,a)') &
    '==================Summary of completed MOM4 integration======================='

    if(vert_coordinate==GEOPOTENTIAL) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using GEOPOTENTIAL as the vertical coordinate.' 
    elseif(vert_coordinate==ZSTAR) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using ZSTAR as the vertical coordinate.' 
    elseif(vert_coordinate==ZSIGMA) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using ZSIGMA as the vertical coordinate.' 
    elseif(vert_coordinate==PRESSURE) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using PRESSURE as the vertical coordinate.' 
    elseif(vert_coordinate==PSTAR) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using PSTAR as the vertical coordinate.' 
    elseif(vert_coordinate==PSIGMA) then 
       write(stdoutunit,'(1x,a)') &
       ' Finished MOM4 integration using PSIGMA as the vertical coordinate.' 
    endif 
 
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of time steps                             = ',Time%itt
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of prog-tracers                           = ',num_prog_tracers
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of diag-tracers                           = ',num_diag_tracers
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of i-points(ni)                           = ',Grid%ni
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of j-points(nj)                           = ',Grid%nj
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of k-points(nk)                           = ',Grid%nk
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of computed ocean tracer points(ni*nj*nk) = ',Grid%total_t_points
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of wet ocean tracer points                = ',Grid%wet_t_points
    write (stdoutunit,'(1x,a,1x,i12)')  &
    ' number of wet ocean velocity points              = ',Grid%wet_u_points

    write (stdoutunit,'(1x,a,1x,a)')     &
    ' slow motion time tendency computed using         =       ',time_tendency
    write (stdoutunit,'(1x,a,1x,a)')     &
    ' barotropic motion computed using                 = ',Time_steps%barotropic_scheme
    write (stdoutunit,'(1x,a,1x,f16.6)') &
    ' tracer time step dtts (secs)                     = ',Time_steps%dtts
    write (stdoutunit,'(1x,a,1x,f16.6)') &
    ' baroclinic velocity time step dtuv (secs)        = ',Time_steps%dtuv
    write (stdoutunit,'(1x,a,1x,f16.6)') &
    ' barotropic time step dtbt (secs)                 = ',Time_steps%dtbt
    write (stdoutunit,'(1x,a,1x,f16.6)') &
    ' implicit Coriolis parameter acor                 = ',Time_steps%acor
    write (stdoutunit,'(1x,a,1x,f16.6)') &
    ' implicit vertical mixing parameter aidif         = ',Time_steps%aidif

    write (stdoutunit,'(a)')  ' '
    write (stdoutunit,'(2x,a)')             Ocean_options%OBC
    write (stdoutunit,'(2x,a)')             Ocean_options%baroclinic_tendency
    write (stdoutunit,'(2x,a)')             Ocean_options%barotropic_tendency
    write (stdoutunit,'(2x,a)')             Ocean_options%tracer_tendency
    write (stdoutunit,'(2x,a)')             Ocean_options%equation_of_state
    write (stdoutunit,'(2x,a)')             Ocean_options%temperature_variable
    write (stdoutunit,'(2x,a)')             Ocean_options%frazil_ice
    write (stdoutunit,'(2x,a)')             Ocean_options%vert_mix
    write (stdoutunit,'(2x,a)')             Ocean_options%tidal_wave_mix
    write (stdoutunit,'(2x,a)')             Ocean_options%tidal_drag_mix
    write (stdoutunit,'(2x,a)')             Ocean_options%bryan_lewis_mix
    write (stdoutunit,'(2x,a)')             Ocean_options%hwf_mix
    write (stdoutunit,'(2x,a)')             Ocean_options%tanh_diff_cbt
    write (stdoutunit,'(2x,a)')             Ocean_options%horz_bih_tracer
    write (stdoutunit,'(2x,a)')             Ocean_options%horz_lap_tracer
    write (stdoutunit,'(2x,a)')             Ocean_options%horz_lap_friction
    write (stdoutunit,'(2x,a)')             Ocean_options%horz_bih_friction
    write (stdoutunit,'(2x,a)')             Ocean_options%momentum_source
    write (stdoutunit,'(2x,a)')             Ocean_options%form_drag
    write (stdoutunit,'(2x,a)')             Ocean_options%bottom_roughness
    write (stdoutunit,'(2x,a)')             Ocean_options%geothermal_heating
    write (stdoutunit,'(2x,a)')             Ocean_options%neutral_physics
    write (stdoutunit,'(2x,a)')             Ocean_options%submesoscale
    write (stdoutunit,'(2x,a)')             Ocean_options%sigma_transport
    write (stdoutunit,'(2x,a)')             Ocean_options%overflow
    write (stdoutunit,'(2x,a)')             Ocean_options%overexchange
    write (stdoutunit,'(2x,a)')             Ocean_options%mixdownslope
    write (stdoutunit,'(2x,a)')             Ocean_options%shortwave
    write (stdoutunit,'(2x,a)')             Ocean_options%xlandmix
    write (stdoutunit,'(2x,a)')             Ocean_options%xlandinsert
    write (stdoutunit,'(2x,a)')             Ocean_options%rivermix
    write (stdoutunit,'(2x,a)')             Ocean_options%riverspread
    write (stdoutunit,'(2x,a)')             Ocean_options%passive_tracers
    write (stdoutunit,'(2x,a)')             Ocean_options%ocean_sponges_eta
    write (stdoutunit,'(2x,a)')             Ocean_options%ocean_sponges_tracer
    write (stdoutunit,'(2x,a)')             Ocean_options%ocean_sponges_velocity
    write (stdoutunit,'(2x,a/)')   '===================================================='
    
    return

  end subroutine ocean_model_end
! </SUBROUTINE> NAME="ocean_model_end"

!#######################################################################
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
  subroutine ocean_model_restart(Ocean_state, timestamp)
     type(ocean_state_type),    pointer     :: Ocean_state
     character(len=*), intent(in), optional :: timestamp

     call ocean_tracer_restart(Time, T_prog, timestamp)
     call ocean_tracer_advect_restart(T_prog, timestamp)
     call ocean_nphysics_restart(timestamp)
     call ocean_bih_friction_restart(timestamp)
     call ocean_lap_friction_restart(timestamp)
     call ocean_sigma_transport_restart(timestamp)
     call ocean_velocity_restart(Time, Velocity, timestamp)
     call ocean_barotropic_restart(Time, Ext_mode, timestamp)
     call ocean_thickness_restart(Time, Thickness, timestamp)
     call ocean_density_restart(Time, Dens, timestamp)
     if(have_obc) call ocean_obc_restart(timestamp)
     call ocean_sfc_restart(timestamp)
     call ocean_vert_mix_restart(timestamp)

  end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_stock_pe">
!
! <DESCRIPTION>
! Returns stocks of total ocean heat and water water for conservation 
! checks.  Report here values just on a single PE. Global sums 
! are done in the coupler. 
!
! This routine is part of a group of similar routines in other 
! FMS component models that aims to quantify the conservation of 
! scalar properties between the component models when running 
! coupled models. 
! 
! </DESCRIPTION>
!
subroutine ocean_stock_pe(Ocean_state, index, value, time_index)
  type(ocean_state_type),pointer     :: Ocean_state
  integer,               intent(in)  :: index
  real,                  intent(out) :: value
  integer, optional,     intent(in)  :: time_index ! -1=previous, 0=now, or +1=next

  integer  :: i,j,k
  integer  :: itime, t_index 
  integer, parameter :: INDEX_PREV=-1, INDEX_NOW=0, INDEX_NEXT=+1

  value = 0.0

  if (.not.associated(Ocean_state)) return
  if(.not. Ocean_state%is_ocean_pe) return

  !The default t_index has to be INDEX_NEXT because stocks are evaluated after 
  !ocean update (values at taup1 are set) but before we are in the next time step.
  t_index = INDEX_NEXT 
  if(present(time_index)) t_index = time_index


  select case(t_index)

  ! note that for time_tendency=='twolevel', taum1=tau
  case (INDEX_PREV)
     itime = Time % taum1

  case (INDEX_NOW)
     itime = Time % tau

  case default
     itime = Time % taup1

  end select


  select case (index)

  ! units of kg  
  case (ISTOCK_WATER)
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               value = value + Grid%tmask(i,j,k)*Grid%dat(i,j)*Thickness%rho_dzt(i,j,k,itime)
            enddo
         enddo
      enddo

  ! units of Joule 
  case (ISTOCK_HEAT)
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               value = value + Grid%tmask(i,j,k)*Grid%dat(i,j)*Thickness%rho_dzt(i,j,k,itime) &
                               *T_prog(index_temp)%field(i,j,k,itime)
            enddo
         enddo
      enddo
      value = value*T_prog(index_temp)%conversion

   case (ISTOCK_SALT)
      ! Return the mass of the salt in the ocean on this PE in kg.
      ! The T_prog(index_salt)%conversion = 1000 converts salinity in PSU to saltinity in kg kg-1.
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               value = value + Grid%tmask(i,j,k)*Grid%dat(i,j)*Thickness%rho_dzt(i,j,k,itime) &
                               *T_prog(index_salt)%field(i,j,k,itime)
            enddo
         enddo
      enddo
      value = value*T_prog(index_salt)%conversion


  case default

      value = 0.0

  end select 


end subroutine ocean_stock_pe
! </SUBROUTINE> NAME="ocean_stock_pe"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_Tsurf">
!
! <DESCRIPTION>
! Return the surface temperature in degrees K
! </DESCRIPTION>
!
subroutine mom4_get_Tsurf(Ocean, res)

  type(Ocean_public_type) :: Ocean
  real, intent(out)     :: res(Domain%isc:, Domain%jsc:)

  integer :: isc, iec, jsc, jec, i, j

  isc = Domain%isc ; iec = Domain%iec
  jsc = Domain%jsc ; jec = Domain%jec

  do j=jsc,jec
     do i=isc,iec
        res(i,j) = Ocean%t_surf(i,j)
     enddo
  enddo

end subroutine mom4_get_Tsurf
! </SUBROUTINE> NAME="mom4_get_Tsurf"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_Ssurf">
!
! <DESCRIPTION>
! Return the surface salinity in psu
! </DESCRIPTION>
!
subroutine mom4_get_Ssurf(Ocean, res)

  type(Ocean_public_type) :: Ocean
  real, intent(out)     :: res(Domain%isc:, Domain%jsc:)

  integer :: isc, iec, jsc, jec, i, j

  isc = Domain%isc ; iec = Domain%iec
  jsc = Domain%jsc ; jec = Domain%jec

  do j = jsc, jec
     do i = isc, iec
        res(i,j) = Ocean%s_surf(i,j)
     enddo
  enddo

end subroutine mom4_get_Ssurf
! </SUBROUTINE> NAME="mom4_get_Ssurf"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_thickness">
!
! <DESCRIPTION>
! Return thickness (in meters) of each layer.
! </DESCRIPTION>
!
subroutine mom4_get_thickness(fld)

  real, intent(out) :: fld(isc:,jsc:,:)

  integer :: i, j, k

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           fld(i,j,k) = Thickness%dzt(i,j,k)
        enddo
     enddo
  enddo
  
end subroutine mom4_get_thickness
! </SUBROUTINE> NAME="mom4_get_thickness"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_density">
!
! <DESCRIPTION>
! Return density (in kg/m^3).
! </DESCRIPTION>
!
subroutine mom4_get_density(fld)

  real, intent(out) :: fld(isc:,jsc:,:)

  integer :: i, j, k, tau

  tau = Time%tau

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           fld(i,j,k) = Dens%rho(i,j,k,tau)
        enddo
     enddo
  enddo
  
end subroutine mom4_get_density
! </SUBROUTINE> NAME="mom4_get_density"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_prog_tracer">
!
! <DESCRIPTION>
! Return prognostic tracer data.
! </DESCRIPTION>
!
subroutine mom4_get_prog_tracer(index, fld, units, longname)

  integer, intent(in) :: index 
  real, intent(out)   :: fld(Domain%isc:,Domain%jsc:,:)
  character(len=*), optional, intent(out)  :: units
  character(len=*), optional, intent(out)  :: longname

  integer :: i, j, k, tau

  tau = Time%tau

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           fld(i,j,k) = T_prog(index)%field(i,j,k,tau)
        enddo
     enddo
  enddo
  
  if(present(units   )) units    = T_prog(index)%units
  if(present(longname)) longname = T_prog(index)%longname
  
end subroutine mom4_get_prog_tracer
! </SUBROUTINE> NAME="mom4_get_prog_tracer"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_temperature_index">
!
! <DESCRIPTION>
! Return temperature index from prognostic tracer table, which can 
! then be used to extract data.
! </DESCRIPTION>
!
subroutine mom4_get_temperature_index(index)

  integer, intent(out) :: index

  index = INDEX_TEMP
  
end subroutine mom4_get_temperature_index
! </SUBROUTINE> NAME="mom4_get_temperature_index"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_salinity_index">
!
! <DESCRIPTION>
! Return salt index from prognostic tracer table, which can 
! then be used to extract data.
! </DESCRIPTION>
!
subroutine mom4_get_salinity_index(index)

  integer, intent(out) :: index

  index = INDEX_SALT
  
end subroutine mom4_get_salinity_index
! </SUBROUTINE> NAME="mom4_get_salinity_index"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_dimensions">
!
! <DESCRIPTION>
! Return dimensions of data in compute domain
! </DESCRIPTION>
!
subroutine mom4_get_dimensions(isc, iec, jsc, jec, isd, ied, jsd, jed, nk_out)

  integer, optional, intent(out) :: isc ! starting x-index of compute domain
  integer, optional, intent(out) :: iec ! ending x-index of compute domain
  integer, optional, intent(out) :: jsc ! starting y-index of compute domain
  integer, optional, intent(out) :: jec ! ending y-index of compute domain

  integer, optional, intent(out) :: isd ! starting x-index of data domain
  integer, optional, intent(out) :: ied ! ending x-index of data domain
  integer, optional, intent(out) :: jsd ! starting y-index of data domain
  integer, optional, intent(out) :: jed ! ending y-index of data domain

  integer, optional, intent(out) :: nk_out  ! number of vertical levels

  if(present(isc)) isc = Domain%isc
  if(present(iec)) iec = Domain%iec
  if(present(jsc)) jsc = Domain%jsc
  if(present(jec)) jec = Domain%jec

  if(present(isd)) isd = Domain%isd
  if(present(ied)) ied = Domain%ied
  if(present(jsd)) jsd = Domain%jsd
  if(present(jed)) jed = Domain%jed

  if(present(nk_out)) nk_out = NK
  
end subroutine mom4_get_dimensions
! </SUBROUTINE> NAME="mom4_get_dimensions"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_UVsurf">
!
! <DESCRIPTION>
! Return horizontal velocity vector components (u,v) on the 
! A grid (tracer-points).
!
! Note that these velocity components are oriented according to the 
! grid lines (i-lines and j-lines).  They are generally NOT mapped 
! to latitude-longitude lines, unless using a spherical coordinate 
! grid specification.  
! </DESCRIPTION>
!
subroutine mom4_get_UVsurf(Ocean, ua, va, ier)

  type(ocean_public_type)    :: Ocean
  real, dimension(Domain%isc:,Domain%jsc:), intent(out)        :: ua, va
  integer, intent(out)     :: ier  ! error code (0=ok)
  
  
  call mom4_U_to_T_2d(ub=Ocean%u_surf, &
       &              vb=Ocean%v_surf, &
       &              ua=ua, va=va, ier=ier)

end subroutine mom4_get_UVsurf
! </SUBROUTINE> NAME="mom4_get_UVsurf"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_UV">
!
! <DESCRIPTION>
! Return horizontal velocity vector (u,v) (in m/s) on T points (A mesh).
!
! Note that these velocity components are oriented according to the 
! grid lines (i-lines and j-lines).  They are generally NOT mapped 
! to latitude-longitude lines, unless using a spherical coordinate 
! grid specification. 
! </DESCRIPTION>
!
subroutine mom4_get_UV(ua, va, ier)

  real, dimension(Domain%isc:, Domain%jsc:, :), intent(out) :: ua, va
  integer, intent(out) :: ier ! error code (0=ok)

  real, dimension(Domain%isc:Domain%iec, Domain%jsc:Domain%jec, nk) :: ub_tmp, vb_tmp

  integer :: i, j, k, tau

  tau   = Time%tau

  do k=1,nk

     do j=jsc,jec
        do i=isc,iec
           ub_tmp(i,j,k) = Velocity%u(i,j,k,1,tau)
           vb_tmp(i,j,k) = Velocity%u(i,j,k,2,tau)
        enddo
     enddo
     
     call mom4_U_to_T_2d(ub=ub_tmp(Domain%isc:, Domain%jsc:, k), &
          &              vb=vb_tmp(Domain%isc:, Domain%jsc:, k), &
          &              ua=ua(Domain%isc:, Domain%jsc:, k), &
          &              va=va(Domain%isc:, Domain%jsc:, k), ier=ier)

  enddo
  
end subroutine mom4_get_UV
! </SUBROUTINE> NAME="mom4_get_UV"


!#######################################################################
! <SUBROUTINE NAME="mom4_U_to_T_2d">
!
! <DESCRIPTION>
! Interpolate (u,v) velocity components from U (B-grid) to 
! T points (A-grid).
! </DESCRIPTION>
!
subroutine mom4_U_to_T_2d(ub, vb, ua, va, ier)

  real, dimension(Domain%isc:, Domain%jsc:), intent(in)  :: ub, vb
  real, dimension(Domain%isc:, Domain%jsc:), intent(out) :: ua, va
  integer, intent(out) :: ier

  integer :: isc, iec, jsc, jec, i1, j1, i, j
  real, dimension(Domain%isd:Domain%ied, Domain%jsd:Domain%jed) :: ub_tmp, vb_tmp

  ier = 0

  if( Domain%yhalo < 1 .or. Domain%xhalo < 1 ) then
     ! error; need at least one halo point
     ier = 1
     return
  endif

  isc = Domain%isc ; iec = Domain%iec
  jsc = Domain%jsc ; jec = Domain%jec

  ! For the time being we assume the input arrays to be on the compute 
  ! domain. Ideally would want to handle input arrays both on the compute
  ! and data domain. In the latter case, a copy + mpp_update can be avoided.
  if(size(ub, dim=1) /= iec-isc+1 .or. size(vb, dim=2) /= jec-jsc+1) then
     ! input arrays must be on the compute domain
     ier = 2
     return
  endif

  ub_tmp(isc:iec, jsc:jec) = ub(isc:iec, jsc:jec)
  vb_tmp(isc:iec, jsc:jec) = vb(isc:iec, jsc:jec)

  call mpp_update_domains(ub_tmp, vb_tmp, Domain%domain2d)

  do j = jsc-1, jec-1
     j1 = j + 1
     do i = isc-1, iec-1
        i1 = i + 1
        ua(i1,j1) = ( &
             & ub_tmp(i ,j ) * Grid%dtn(i1,j1) * Grid%dte(i1,j1) + &
             & ub_tmp(i1,j ) * Grid%dtn(i1,j1) * Grid%dtw(i1,j1) + &
             & ub_tmp(i1,j1) * Grid%dts(i1,j1) * Grid%dtw(i1,j1) + &
             & ub_tmp(i ,j1) * Grid%dts(i1,j1) * Grid%dte(i1,j1) &
             & ) / Grid%dat(i1,j1)
        va(i1,j1) = ( &
             & vb_tmp(i ,j ) * Grid%dtn(i1,j1) * Grid%dte(i1,j1) + &
             & vb_tmp(i1,j ) * Grid%dtn(i1,j1) * Grid%dtw(i1,j1) + &
             & vb_tmp(i1,j1) * Grid%dts(i1,j1) * Grid%dtw(i1,j1) + &
             & vb_tmp(i ,j1) * Grid%dts(i1,j1) * Grid%dte(i1,j1) &
             & ) / Grid%dat(i1,j1)
     enddo
  enddo

end subroutine mom4_U_to_T_2d
! </SUBROUTINE> NAME="mom4_U_to_T_2d"

!#######################################################################
! <SUBROUTINE NAME="mom4_get_latlon_UV">
!
!   <OVERVIEW>
!     Gets horizontal velocity components (u,v) (in m/s) on T points (A mesh) 
!     along geographical (latlon) directions in compute domain.
!   </OVERVIEW>
!   <DESCRIPTION>
!    Note that these velocity components are oriented along the 
!    geographical latitude-longitude lines.  
!
!           im,j    i,j
!  B-------B-------B-------B       y
!  |       |       |       |       ^
!  |       |    i,j|       |       |   /lon
!  |---A---|---A---|---A---|       |  /   
!  |       |       |       |     \ | /
!  |       |im,jm  |i,jm   |      \|/ rot angle
!  B-------B-------B-------B    ---X-------------> x
!  |       |       |       |      /|\ 
!  |       |       |       |     / | \
!  |---A---|---A---|---A---|       |  \lat  
!  |       |       |       |       |   \   
!  |       |       |       |
!  B-------B-------B-------B
!
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_model_mod 
!     real, dimension(isc:, jsc:, :) :: u,v
!     integer :: ierr
!     call mom4_get_latlon_UV(ua, va, ierr)
!   </TEMPLATE>
!   <OUT NAME="ua"  TYPE="real, dimension(isc:, jsc:, :)" >
!     array will contain velocity component along x direction upon return
!   </OUT>
!   <OUT NAME="va"  TYPE="real, dimension(isc:, jsc:, :)" >
!     array will contain velocity component along y direction upon return
!   </OUT>
!   <OUT NAME="ierr"  TYPE="ineger" >
!     error status will be zero for success and nonzero for failure 
!   </OUT>
!
subroutine mom4_get_latlon_UV(ua, va, ier)

  real, dimension(isc:, jsc:, :), intent(out) :: ua, va
  integer, intent(out) :: ier ! error code (0=ok)

  real, dimension(isd:ied, jsd:jed) :: ub_tmp, vb_tmp

  integer :: i, j, k, im, jm, tau

  ier = 0
  tau   = Time%tau

  do k=1,nk

     ! Rotate the coordinates system to be along geographical (lat-lon) coordinates.
     ! (Passive transformation with rotation matrix).
     !
     do j=jsc,jec
        do i=isc,iec
           ub_tmp(i,j) = Grid%cos_rot(i,j)*Velocity%u(i,j,k,1,tau) &
           &           - Grid%sin_rot(i,j)*Velocity%u(i,j,k,2,tau) 
             
           vb_tmp(i,j) = Grid%sin_rot(i,j)*Velocity%u(i,j,k,1,tau) &
           &           + Grid%cos_rot(i,j)*Velocity%u(i,j,k,2,tau) 
        enddo
     enddo

     ! Fill in the halos
     call mpp_update_domains(ub_tmp, vb_tmp, Domain%domain2d)

     ! Take the average of fields at corners (B Grid) to estimate the field at center (A Grid).
     ! Try an area weighted average.

     do j = jsc, jec
        jm = j - 1
        do i = isc, iec
           im = i - 1
           ua(i,j,k) = ( &
                & ub_tmp(i ,j)  * Grid%dtn(i,j) * Grid%dte(i,j) + &
                & ub_tmp(im,j)  * Grid%dtn(i,j) * Grid%dtw(i,j) + &
                & ub_tmp(im,jm) * Grid%dts(i,j) * Grid%dtw(i,j) + &
                & ub_tmp(i ,jm) * Grid%dts(i,j) * Grid%dte(i,j) ) / Grid%dat(i,j)

           va(i,j,k) = ( &
                & vb_tmp(i ,j)  * Grid%dtn(i,j) * Grid%dte(i,j) + &
                & vb_tmp(im,j)  * Grid%dtn(i,j) * Grid%dtw(i,j) + &
                & vb_tmp(im,jm) * Grid%dts(i,j) * Grid%dtw(i,j) + &
                & vb_tmp(i ,jm) * Grid%dts(i,j) * Grid%dte(i,j) ) / Grid%dat(i,j)
        enddo
     enddo
     
  enddo !end k loop
  
end subroutine mom4_get_latlon_UV
! </SUBROUTINE> NAME="mom4_get_latlon_UV"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_diag_axes">
!
! <DESCRIPTION>
! Return axes indices for diag manager.
! </DESCRIPTION>
!
subroutine mom4_get_diag_axes(axes)

  integer, intent(out) :: axes(:)

  axes(:) = Grid%tracer_axes(:)

end subroutine mom4_get_diag_axes
! </SUBROUTINE> NAME="mom4_get_diag_axes"


!#######################################################################
! <FUNCTION NAME="mom4_get_num_diag_tracers">
!   <OVERVIEW>
!     Returns the module variable num_diag_tracers   
!   </OVERVIEW>
!   <DESCRIPTION>
!     This function returns the number of ocean diagnostic tracers if not -1.
!     It send a FATAL message if num_diag_tracers is not set (i.e. is -1) 
!     before this function call.
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_model_mod
!     mom4_get_num_diag_tracers()
!   </TEMPLATE>
!   <IN NAME=""  TYPE="" >
!     No inputs needed.
!   </IN>
!   <OUT NAME=""  TYPE="integer" >
!     This function returns an integer. 
!   </OUT>

function mom4_get_num_diag_tracers()
  integer :: mom4_get_num_diag_tracers

  if(num_diag_tracers == -1) then
     call mpp_error(FATAL, &
     'mom4_get_num_diag_tracers: num_diag_tracers is -1, should call ocean_diag_tracer_init before this function! ')
  endif

  mom4_get_num_diag_tracers = num_diag_tracers

end function mom4_get_num_diag_tracers
! </FUNCTION> NAME="mom4_get_num_diag_tracers"


!#######################################################################
! <FUNCTION NAME="mom4_get_num_prog_tracers">
!   <OVERVIEW>
!     Returns the module variable num_prog_tracers   
!   </OVERVIEW>
!   <DESCRIPTION>
!     This function returns the number of ocean prognostic tracers if not -1.
!     It send a FATAL message if num_prog_tracers is not set (i.e. is -1) 
!     before this function call.
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_model_mod
!     mom4_get_num_prog_tracers()
!   </TEMPLATE>
!   <IN NAME=""  TYPE="" >
!     No inputs needed.
!   </IN>
!   <OUT NAME=""  TYPE="integer" >
!     This function returns an integer. 
!   </OUT>

function mom4_get_num_prog_tracers()
  integer :: mom4_get_num_prog_tracers

  if(num_prog_tracers == -1) then
     call mpp_error(FATAL, &
     'mom4_get_num_prog_tracers: num_prog_tracers is -1, should call ocean_prog_tracer_init before this function! ')
  endif

  mom4_get_num_prog_tracers = num_prog_tracers

end function mom4_get_num_prog_tracers
! </FUNCTION> NAME="mom4_get_num_prog_tracers"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_surface_tmask">
!   <OVERVIEW>
!     Gets the pointer to 2D array Grid%tmask(:,:,1)   
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine gets the pointer to a 2D array with values of the tmask 
!     (land/sea mask for T cells based on s-coordinate) at the ocean surface .  
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_model_mod
!     real, dimension(:,:), pointer :: temp
!     call mom4_get_surface_tmask(temp)
!   </TEMPLATE>
!   <OUT NAME="surfaceTmask"  TYPE="real, dimension(:,:),pointer" >
!      pointer to 2 dimensional array of tmask at the ocean surface. 
!   </OUT>

subroutine mom4_get_surface_tmask(surfaceTmask)
  real, dimension(:,:),pointer  :: surfaceTmask

  surfaceTmask => Grid%tmask(isc:,jsc:,1)

end subroutine mom4_get_surface_tmask
! </SUBROUTINE> NAME="mom4_get_surface_tmask"


!#######################################################################
! <SUBROUTINE NAME="mom4_get_ocean_data">
!   <OVERVIEW>
!     Gets one of the 2D array data of ocean type   
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine gets one of the following data arrays of ocean_public_type
!     depending on the passed "name" argument, it sends a FATAL signal otherwise
!     Ocean%t_surf when name='t_surf' 
!     Ocean%s_surf when name='s_surf' 
!     Ocean%u_surf when name='u_surf' 
!     Ocean%v_surf when name='v_surf' 
!     Ocean%sea_lev when name='sea_lev' 
!     Ocean%frazil when name='frazil' 
!   </DESCRIPTION>
!   <TEMPLATE>
!     use ocean_model_mod
!     real, dimension(:,:), pointer :: temp
!     call mom4_get_ocean_data(Ocean,'s_surf',temp)
!   </TEMPLATE>
!   <IN NAME="Ocean"  TYPE="type(ocean_public_type)" >
!      ocean type
!   </IN>
!   <IN NAME="name"  TYPE="character(len=*)" >
!      one of 't_surf','s_surf','u_surf','v_surf','sea_lev','frazil'
!   </IN>
!   <OUT NAME="dataArrayPointer"  TYPE="real, dimension(:,:), pointer" >
!      pointer to 2 dimensional array corresponding to "name" argument, at the ocean surface. 
!   </OUT>

subroutine mom4_get_ocean_data(Ocean,name,dataArrayPointer)
  type(ocean_public_type), intent(in)    :: Ocean
  character(len=*), intent(in)  :: name
  real, dimension(:,:),pointer  :: dataArrayPointer
  

  select case(name)
  case('t_surf')
     if(.not.associated(Ocean%t_surf)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%t_surf is not associated!')  
     dataArrayPointer => Ocean%t_surf
  case('s_surf')
     if(.not.associated(Ocean%s_surf)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%s_surf is not associated!')  
     dataArrayPointer => Ocean%s_surf
  case('u_surf')
     if(.not.associated(Ocean%u_surf)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%u_surf is not associated!')  
     dataArrayPointer => Ocean%u_surf
  case('v_surf')
     if(.not.associated(Ocean%v_surf)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%v_surf is not associated!')  
     dataArrayPointer => Ocean%v_surf
  case('sea_lev')
     if(.not.associated(Ocean%sea_lev)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%sea_lev is not associated!') 
     dataArrayPointer => Ocean%sea_lev
  case('frazil')
     if(.not.associated(Ocean%frazil)) &
     call mpp_error(FATAL,'mom4_get_ocean_data: Ocean%frazil is not associated!')  
     dataArrayPointer => Ocean%frazil
  case default
     call mpp_error(FATAL,'mom4_get_ocean_data: unknown argument name='//name)
  end select
  

end subroutine mom4_get_ocean_data
! </SUBROUTINE> NAME="mom4_get_ocean_data"

subroutine mom4_set_swheat(swheat_input)
  real, dimension(isc:, jsc:,:), intent(in)  :: swheat_input
  
  ext_swheat_is_set=.true.
  swheat = swheat_input

end subroutine mom4_set_swheat


end module ocean_model_mod
  
