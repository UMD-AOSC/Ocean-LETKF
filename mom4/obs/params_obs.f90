! MODULE: params_obs
!
! This module contains all of the necessary parameters related to the
! observations, and observation operators.
!
! Author: Prof. Stephen G. Penny
!         University of Maryland, College Park
!         Department of Atmospheric and Oceanic Science
!
! 2016.4.7


MODULE params_obs

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

INTEGER,SAVE :: nobs
INTEGER,PARAMETER :: nid_obs=8
INTEGER,PARAMETER :: id_u_obs=2819
INTEGER,PARAMETER :: id_v_obs=2820
INTEGER,PARAMETER :: id_t_obs=3073
INTEGER,PARAMETER :: id_s_obs=5521      !(OCEAN)
INTEGER,PARAMETER :: id_ssh_obs=5526    !(OCEAN)
INTEGER,PARAMETER :: id_sst_obs=5525    !(OCEAN)
INTEGER,PARAMETER :: id_sss_obs=5522    !(OCEAN)
INTEGER,PARAMETER :: id_eta_obs=5351    !(OCEAN)
INTEGER,PARAMETER :: id_sic_obs=6282    !(SEAICE)
INTEGER,PARAMETER :: id_x_obs=1111   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters
INTEGER,PARAMETER :: id_y_obs=2222   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters
INTEGER,PARAMETER :: id_z_obs=3333   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters

!!------------------------------------------------------------
!! unique ID's for observations
!! STEVE: the following will replace what is above:
!!------------------------------------------------------------
!! atmosphere obs
INTEGER, PARAMETER :: obsid_atm_min  = 1000
INTEGER, PARAMETER :: obsid_atm_max  = 1999
INTEGER, PARAMETER :: obsid_atm_num  = 8
INTEGER, PARAMETER :: obsid_atm_offset = 0

INTEGER, PARAMETER :: obsid_atm_ps   = 1100
INTEGER, PARAMETER :: obsid_atm_rain = 1110
INTEGER, PARAMETER :: obsid_atm_t    = 1210
INTEGER, PARAMETER :: obsid_atm_tv   = 1211
INTEGER, PARAMETER :: obsid_atm_q    = 1220
INTEGER, PARAMETER :: obsid_atm_rh   = 1221
INTEGER, PARAMETER :: obsid_atm_u    = 1250
INTEGER, PARAMETER :: obsid_atm_v    = 1251

!! ocean obs
INTEGER, PARAMETER :: obsid_ocn_min  = 2000
INTEGER, PARAMETER :: obsid_ocn_max  = 2999
INTEGER, PARAMETER :: obsid_ocn_num  = 8
INTEGER, PARAMETER :: obsid_ocn_offset = obsid_atm_offset + obsid_atm_num

INTEGER, PARAMETER :: obsid_ocn_ssh  = 2100
INTEGER, PARAMETER :: obsid_ocn_eta  = 2101
INTEGER, PARAMETER :: obsid_ocn_sst  = 2110
INTEGER, PARAMETER :: obsid_ocn_sss  = 2120
INTEGER, PARAMETER :: obsid_ocn_t    = 2210
INTEGER, PARAMETER :: obsid_ocn_s    = 2220
INTEGER, PARAMETER :: obsid_ocn_u    = 2250
INTEGER, PARAMETER :: obsid_ocn_v    = 2251
!->
INTEGER, PARAMETER :: obsid_ocn_x    = 2301
INTEGER, PARAMETER :: obsid_ocn_y    = 2302
INTEGER, PARAMETER :: obsid_ocn_z    = 2303

!! sea-ice obs
INTEGER, PARAMETER :: obsid_sic_min  = 3000
INTEGER, PARAMETER :: obsid_sic_max  = 3999
INTEGER, PARAMETER :: obsid_sic_num  = 1
INTEGER, PARAMETER :: obsid_sic_offset = obsid_ocn_offset + obsid_ocn_num

INTEGER, PARAMETER :: obsid_sic_con  = 3100

!! land obs
INTEGER, PARAMETER :: obsid_lnd_min  = 4000
INTEGER, PARAMETER :: obsid_lnd_max  = 4999
INTEGER, PARAMETER :: obsid_lnd_num  = 1
INTEGER, PARAMETER :: obsid_lnd_offset = obsid_sic_offset + obsid_sic_num

INTEGER, PARAMETER :: obsid_lnd_wat  = 4100

!! wave obs
INTEGER, PARAMETER :: obsid_wav_min  = 5000
INTEGER, PARAMETER :: obsid_wav_max  = 5999
INTEGER, PARAMETER :: obsid_wav_num  = 1
INTEGER, PARAMETER :: obsid_wav_offset = obsid_lnd_offset + obsid_lnd_num

INTEGER, PARAMETER :: obsid_wav_hgt  = 5100

!! aerosol obs
INTEGER, PARAMETER :: obsid_aer_min  = 6000
INTEGER, PARAMETER :: obsid_aer_max  = 6999
INTEGER, PARAMETER :: obsid_aer_num  = 1
INTEGER, PARAMETER :: obsid_aer_offset = obsid_wav_offset + obsid_wav_num

INTEGER, PARAMETER :: obsid_aer_aod  = 6100

!-------------------------------------------------------------------------------
! arrays holding all observation id's and names, for easy iteration
! in loops that want to print stats for obs
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: obsid_num = 16
INTEGER, PARAMETER, DIMENSION(obsid_num) :: obsid_array = (/&
       obsid_atm_ps, obsid_atm_rain, obsid_atm_t, obsid_atm_tv, &
       obsid_atm_q, obsid_atm_rh, obsid_atm_u, obsid_atm_v, &
       obsid_ocn_ssh, obsid_ocn_eta, obsid_ocn_sst, obsid_ocn_sss, &
       obsid_ocn_t, obsid_ocn_s, obsid_ocn_u, obsid_ocn_v/)
CHARACTER (len=10) :: obsid_names(obsid_num) = (/&
       "ATM_PS  ", "ATM_RAIN", "ATM_T   ", "ATM_TV  ", &
       "ATM_Q   ", "ATM_RH  ", "ATM_U   ", "ATM_V   ", &
       "OCN_SSH ", "OCN_ETA ", "OCN_SST ", "OCN_SSS ",&
       "OCN_T   ", "OCN_S   ", "OCN_U   ", "OCN_V   "/)

!-------------------------------------------------------------------------------
! Number of records in obs1 or obs2 formatted observation input binary files. 
! ISSUE: make these namelist controllable:
!-------------------------------------------------------------------------------
INTEGER :: obs1nrec = 6                     ! The number of records in the obs1-formatted file (previous 6, 7 adds a time record).
INTEGER :: obs2nrec = 9                     ! The number of records in the obs2-formatted file (previous 8, 9 adds a time record).

!-------------------------------------------------------------------------------
! Remove all observations above 65ºN due to tripolar grid
!-------------------------------------------------------------------------------
LOGICAL :: DO_REMOVE_65N = .false.

!-------------------------------------------------------------------------------
! Temperature conversion method for compting OMFs
!-------------------------------------------------------------------------------
LOGICAL            :: DO_POTTEMP_to_INSITU = .true.  ! Conversion to observation space. This is needed if the
                                                     ! observations aren't converted to potential temperature
                                                     ! (as is done by most - NCEP, SODA, NASA/GMAO, etc.). But
                                                     ! unlike that approach, this does not require synthetic salinity
                                                     ! observations to be constructed from climatologies.
                                                     ! This approach is theoretically better, but investigation must
                                                     ! be done to ensure model biases do not cause significant errors.
                                                     ! (a warning from J. Carton of potential difficulty)
                                                     !
                                                     ! Only one can be true, this one takes prioirty
                                                     !
LOGICAL            :: DO_INSITU_to_POTTEMP = .false. ! Technically, this would require matching an observed salinity
                                                     ! measurement with each observed in situ temperature measurement
                                                     ! and using it to compute the potential temperature. The opposite
                                                     ! process is quite a bit easier.

END MODULE params_obs
