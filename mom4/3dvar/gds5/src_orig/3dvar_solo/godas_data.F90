module godas_data_mod
!
use time_manager_mod, only: time_type
!
public
!
!   Assimilation data for MOM
!
!   Basic paramters set by namelist (see godas_init)
!
!   num_cor_tracers  = number of tracers to be corrected (def:2, T&S)
!   kass             = number of vertical levels in assimilation
!   nsgobs           = number of segments of subsurface observations used per run
!   nsgsobs          = number of segments of surface observations used per run
!   maxits           = maximum iterations in oi analysis scheme
!   npits            = number of iterations used in lsmth to define [e]
!   asm_Tz           = will be set .true. if T(z) OBS file found in data_table
!   asm_Sz           = will be set .true. if S(z) OBS file found in data_table
!   asm_T0           = will be set .true. if SST OBS file found in data_table
!   asm_S0           = will be set .true. if SSS OBS file found in data_table
!   asm_Al           = will be set .true. if Altimetry OBS file found in data_table
!   asm_bTz          = will be set .true. if T(z) BAS file found in data_table
!   asm_bSz          = will be set .true. if S(z) BAS file found in data_table
!   asm_ts_seq       = if .true. temperature and salinity are assimilated sequentially
!   asm_sfc_split    = if .true. SST/SSS and T(z)/S(z) assimilations are split (SST & SSS done sequentially)
!   godas_at_end     = if .true. a single analysis is done just before the run concludes, some other
!                      parameters (e.g. gds_step, scl_incr, single_incr, no_asm_rep) become irrelevant
!   gds_step         = interval between GODAS analyses (sec)
!   scl_incr         = model_step / gds_step, fraction of increment added per time step
!   single_incr      = if .true. increment is added once per analysis cycle, not once per time step
!   no_asm_rep       = number of repeat assimilations at each triggering 
!   tovrf            = global scaling for observed vertical temperature error variance
!   sovrf            = global scaling for observed vertical salinity error variance
!   tbvrf            = global scaling for background vertical temperature error variance
!   sbvrf            = global scaling for background vertical salinity error variance
!   tov0f            = global scaling for observed surface temperature error variance
!   sov0f            = global scaling for observed surface salinity error variance
!   tbv0f            = global scaling for background surface temperature error variance
!   sbv0f            = global scaling for background surface salinity error variance
!   hrscl            = horizontal covariance scale (deg)
!   hrscl0           = horizontal covariance scale (deg) for the surface if asm_sfc_split=.true.
!   vcvn             = vertical covariance normalization factor
!   vsclf            = vertical covariance scale factor (scale is vsclf*(layer thickness))

! for bias analysis
!   tovrf_b          = global scaling for observed vertical temperature error variance
!   sovrf_b          = global scaling for observed vertical salinity error variance
!   tbvrf_b          = global scaling for background vertical temperature error variance
!   sbvrf_b          = global scaling for background vertical salinity error variance
!   hrscl_b          = horizontal covariance scale (deg)
!   vcvn_b           = vertical covariance normalization factor
!   vsclf_b          = vertical covariance scale factor (scale is vsclf*(layer thickness))
!   yscl_b           = scale for northern tapering factor : 1.0 - exp[ -(no_lat_mx - y)^2 / yscl^2 ]

!   no_lat_mx        = northern latitude maximum for assimilation
!   jemx             = northern imdex maximum for assimilation
!   yscl             = scale for northern tapering factor : 1.0 - exp[ -(no_lat_mx - y)^2 / yscl^2 ]
!   tz_wndw_fwd      = forward time window for temperature profiles (days)
!   tz_wndw_bwd      = backward time window for temperature profiles (days)
!   sz_wndw_fwd      = forward time window for salinity profiles (days)
!   sz_wndw_bwd      = backward time window for salinity profiles (days)
!   t0_wndw_fwd      = forward time window for SST (days)
!   t0_wndw_bwd      = backward time window for SST (days)
!   s0_wndw_fwd      = forward time window for SSS (days)
!   s0_wndw_bwd      = backward time window for SSS (days)
!   al_wndw_fwd      = forward time window for altimetry (days)
!   al_wndw_bwd      = backward time window for altimetry (days)
!   wndw_secs        = seconds added to all forward/backward windows (currently used only when
!                      godas_at_end=.true., it is then 43200, otherwise 0)
!   assrestrt        = if .true. an assimilation is done at restart
!   restore_sfc      = have GODAS do surface restoration rather than MOM4
!   num_rstr_tracers = number of tracers to be restored (def:0, SST, SSS)
!   rstr_time        = limits surface restore to one time of day (hour,minute,second)
!                      if defaults are set (-1,-1,-1) restore is done each time assimilation is done
!   sst_damp         = damping factor for SST (0: no damping, 1: set to restore field)
!   sss_damp         = damping factor for SSS (0: no damping, 1: set to restore field)
!   rstrestrt        = if .true. a restoration is done at restart
!   save_all_inv     = if .true. save all innovations etc. for post processing,
!                      if .false. only innovations etc. within 1 timestep of analysis are saved 
!   debug_godas      = if .true. print debug information

integer :: num_cor_tracers, kass, nsgobs, nsgsobs, maxits, npits, gds_step, no_asm_rep, jemx
real    :: scl_incr, tovrf, sovrf, tbvrf, sbvrf, hrscl, hrscl0, vcvn, vsclf, no_lat_mx, yscl, ys2
real    :: tov0f, sov0f, tbv0f, sbv0f
real    :: tovrf_b, sovrf_b, tbvrf_b, sbvrf_b, hrscl_b, vcvn_b, vsclf_b, yscl_b
integer :: tz_wndw_fwd, tz_wndw_bwd, sz_wndw_fwd, sz_wndw_bwd
integer :: t0_wndw_fwd, t0_wndw_bwd, s0_wndw_fwd, s0_wndw_bwd
integer :: al_wndw_fwd, al_wndw_bwd, wndw_secs
logical :: single_incr, assrestrt, asm_ts_seq, asm_sfc_split, save_all_inv, godas_at_end, debug_godas
logical :: asm_Tz, asm_Sz, asm_T0, asm_S0, asm_Al, asm_bTz, asm_bSz
integer :: num_rstr_tracers, rstr_time(3)
real    :: sst_damp, sss_damp
logical :: restore_sfc, rstrestrt

logical :: apply_bias_correction

!   ksalt  = k offset for assimilating salinity
!   kass2  = kass + ksalt
!   kass3  = 3*kass; used for multivariate analysis of u and v
!   kass4  = 4*kass; used for multivariate analysis of u and v

integer           :: ksalt, kass2, kass3, kass4, asm_cnt, obs_trk_cnt
integer           :: asm_bias_cnt, obs_bias_trk_cnt
type(time_type)   :: gds_freq, alrm_dur
logical           :: ovr_alrm, apply_incr, apply_bias_incr
integer, allocatable, dimension(:) :: id_cor
integer, allocatable, dimension(:) :: id_rstr
integer, allocatable, dimension(:) :: id_bias_cor

!   rtzw   = inverse of the larger of tz_wndw_fwd and tz_wndw_bwd
!   rszw   = inverse of the larger of sz_wndw_fwd and sz_wndw_bwd
!   rt0w   = inverse of the larger of t0_wndw_fwd and t0_wndw_bwd
!   rs0w   = inverse of the larger of s0_wndw_fwd and s0_wndw_bwd
!   ralw   = inverse of the larger of al_wndw_fwd and al_wndw_bwd

real :: rtzw, rszw, rt0w, rs0w, ralw

! Some constants

integer, parameter   :: spd = 86400

! The observations
!
!   All data sets are organized in arbitrary segments (formerly weeks).
!   Profile data are organized into obs_Z structures and surface data
!   (altimetry data) are organized into obs_0 structures.
!   Climate profile data for bias control are organized into obs_Z structures.
!   Data can be separated into various types, depending on the type
!   of obsevation and the type of platform.  These are managed 
!   externally by data_table entries and internally by integer codes.
!
!     temp_code  =  1, generic T(z) code
!     xbt_code   =  2, xbt T(z) code
!     tao_code   =  3, tao, triton, pirata T(z) code
!     argo_code  =  4, argo T(z) code
!     salt_code  = 11, generic S(z) code
!     stao_code  = 13, tao, triton, pirata S(z) code
!     sargo_code = 14, argo S(z) code
!     sst_code   = 21, generic SST code
!     sss_code   = 23, generic SSS code
!     altm_code  = 24, generic altimetry code
!     tp_code    = 25, topex/poseidon code
!     j1_code    = 26, jason-1 code
!
!     ts_code    = 51, T-S code
!
!     num_obsz     = total number of profiles of all kinds
!     num_obs0     = total number of surface obs of all kinds (excl. alt.)
!     num_obsa     = total number of altimetry obs of all knids
!     num_obsbiasz = total number of profiles of all kinds for bias correction
!
integer, parameter     ::     temp_code  =  1
integer, parameter     ::     xbt_code   =  2
integer, parameter     ::     tao_code   =  3
integer, parameter     ::     argo_code  =  4
integer, parameter     ::     salt_code  = 11
integer, parameter     ::     stao_code  = 13
integer, parameter     ::     sargo_code = 14
integer, parameter     ::     sst_code   = 21
integer, parameter     ::     sss_code   = 23
integer, parameter     ::     altm_code  = 24
integer, parameter     ::     tp_code    = 25
integer, parameter     ::     j1_code    = 26
integer, parameter     ::     ts_code    = 51
integer                ::     num_obsz, num_obs0, num_obsa, num_obsbiasz
!
! Limiting the size of the innovations
!
!   dtemp_max = maximum size of temperature innovations
!   dtemp_elm = increase error of temperature innovations larger than this limit
!   dsalt_max = maximum size of salinity innovations
!   dsalt_elm = increase error of salinity innovations larger than this limit
!   daltm_max = maximum size of altimetry innovations
!   daltm_elm = increase error of altimetry innovations larger than this limit
!
real(kind=4), save     ::  dtemp_max = 10.0
real(kind=4), save     ::  dtemp_elm = 5.0
real(kind=4), save     ::  dsalt_max = 5.0
real(kind=4), save     ::  dsalt_elm = 3.0
real(kind=4), save     ::  daltm_max = 5.0
real(kind=4), save     ::  daltm_elm = 3.0
!
! Limiting the size of the innovations at the surface when the surface anlaysis
!  is done separately from the 3D analysis. These settings effectively place
!  no limit on the size of innovations based on the assumption that QC'd surface
!  OIs are being assimilated.  The limiting code remains in "init_grad_sfc" for
!  possible future adjustment.
!
real(kind=4), save     ::  dsst_max = 100.0
real(kind=4), save     ::  dsst_elm = 100.0
real(kind=4), save     ::  dsss_max = 100.0
real(kind=4), save     ::  dsss_elm = 100.0
!
!  Only the variable part of the altimetry is assimilated. A surface height 
!  climatology must be subtracted from the model surface height. This is typically
!  bootstrapped from a model analysis that assimilates only T(z) and S(z).
!
!  eta_clm = model surface height climatology
!
real, save, dimension(:,:), allocatable :: eta_clm
!
!  The initial surface height innovation is assumed to represent an integral
!  error in the baroclinic field.  Coefficients (based on a linearized
!  dynamic height calculation) are used to distribute the innovation between
!  T and S down through the water column.
!
real, save, dimension(:), allocatable :: cdnz, cdnzs
!
!  Arrays used in specifying the background error covariance
!
real, save, dimension(:), allocatable    :: ev, wrkk       !  ev(kass), wrkk(kass2)
!
!   cvn  = local vertical covariance matrix for temperature
!   vtmp   = vertical variance for temperature
!   vtmp_s = local variance for surface temperature
!
real, save, dimension(:,:), allocatable :: cvn      !  cvn(kass,kass)
real, save, dimension(:,:,:), allocatable :: vtmp
real, save, dimension(:,:), allocatable :: vtmp_s
!
!   cvnsalt = local vertical covariance matrix for salinity
!   vsal    = vertical variance for salinity
!   vsal_s  = local variance for surface salinity
!
real, save, dimension(:,:), allocatable :: cvnsalt      !  cvnsalt(kass,kass)
real, save, dimension(:,:,:), allocatable :: vsal
real, save, dimension(:,:), allocatable :: vsal_s
!
!   Arrays used in minimizing the cost function. 
!
!   In the earlier MOM3 version, the Laplace smoother routine was decomposed
!   in the vertical. This was in contrast to the rest of the GODAS and MOM3
!   code which decomposed in the horizontal.  Here in the MOM4 version, all
!   of GODAS retains the 2D horizontal decomposition of MOM4.
!
!    elipt controls anisotropy of background error covariance (i,j)
!    wgta is normalization for Laplace smoother (lpwghts and lpsmthr) (i,j)
!    xcb etc. are used to hold compute domain decomposition data
!
real, save, dimension(:,:), allocatable :: wgta, wgta_s, elipt
integer, dimension(:), allocatable :: xcb, xce, xcsz, ycb, yce, ycsz
!
!     wcn, wea, wwe, wno, wso are C, E, W, N and S coefficients, allocated as (i,j)
!     wgns, s1, s2, dpth are arrays allocated as (i,j)
!
real, save, dimension(:,:), allocatable :: wcn, wea, wwe, wno, wso
real, save, dimension(:,:), allocatable :: wcn_s, wea_s, wwe_s, wno_s, wso_s
real, save, dimension(:,:), allocatable :: wgns, s1, s2, dpth
!
!    The following arrays represent T, d, e, f, g, and h in Derber 
!    and Rosati (1989) and will be dimensioned (imt,jmt,kass)
!    They are tagged with "_cg" to identify them with the congugate
!    gradient algorithm of the assimilation. The _s extension versions
!    are used with a surface only analysis.
!
real, save, dimension(:,:,:), allocatable :: t_cg, d_cg, e_cg, f_cg, g_cg, h_cg
real, save, dimension(:,:), allocatable :: t_cg_s, d_cg_s, e_cg_s, f_cg_s, g_cg_s, h_cg_s
!
!  For saving innovations
!
!   invdys = save innovations within rsdys of the model time
!   ioexti = file unit for saving extracted temperature innovations
!   ioexsi = file unit for saving extracted salinity innovations
!   ioexei = file unit for saving extracted eta (ALTM) innovations
!   tinvd  = holds filename for temperature innovations
!   sinvd  = holds filename for salinity innovations
!   einvd  = holds filename for eta (ALTM) innovations
!
real, save :: invdys = 0.25
integer, save :: ioexti, ioexsi, ioexei
character(len=10), save :: tinvd, sinvd, einvd
!
! Section for multivariate assimilation
!
!   geofrac   =  accounts for inexactness of geostrophic balance, geofrac < 1
!   yl0, yl1  =  latitudes for averaging geostrophic velocities across equator
!   jl0, jl1  =  indices for averaging geostrophic velocities across equator
!   rl0, rl1  =  interpolation factors for averaging across equator
!   kmva(imt,jstask:jetask) = a "kmt" for u, v multi-variable
!
! application of the geo correction near the equatorial undercurrent is
! controlled by the following variables
!   ylso, ylno = equatorial latitudes to exclude geo correction of undercurrent
!   jlso, jlno = j-indices to exclude geo correction of undercurrent
!   xlwe, xlea = longitudes to exclude geo correction of undercurrent
!   ilwe, ilea = i-indices to exclude geo correction of undercurrent
!  within the bounds set by ylso and ylno, geo correction is controlled by
!  twEq and kwEq (precidence is given to whichever of these 2 conditions 
!  is deeper) and by tcEq.
!   twEq       = sfc currents above this temperature will get geo correction
!   kwEq       = sfc currents k <= kwEq will get geo correction
!   tcEq       = sfc currents below this temperature will get geo correction
!
integer, save :: pneq, jseq, jeeq
integer, save, dimension(:,:), allocatable :: kmva
real, save, dimension(:), allocatable :: rl0, rl1
real, parameter :: geofrac = 0.9
real, parameter :: yl0 = -1.2, yl1 = 1.2
integer, save :: jl0, jl1
real, parameter :: ylso = -2.2, ylno = 2.2, twEq = 25.0, tcEq = 13.0
real, parameter :: xlwe = 39.0, xlea = 103.0
! integer, parameter :: kwEq = 2
integer, parameter :: kwEq = 3
! the following settings should allow geostrophy applied everywhere
! real, parameter :: ylso = -2.2, ylno = 2.2, twEq = 1.0, tcEq = 1.0
! real, parameter :: xlwe = 39.0, xlea = 103.0
! integer, parameter :: kwEq = 35
integer, save :: jlso, jlno, ilwe, ilea
real, save, dimension(:,:), allocatable :: awrk
!
! End of section for multivariate assimilation
!
end module godas_data_mod
