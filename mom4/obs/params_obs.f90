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
INTEGER,PARAMETER :: nid_dobs=3
INTEGER,PARAMETER :: id_x_obs=1111   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters
INTEGER,PARAMETER :: id_y_obs=2222   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters
INTEGER,PARAMETER :: id_z_obs=3333   !(OCEAN) (DRIFTERS) !STEVE: may want to change this depending on type of drifters
INTEGER,PARAMETER :: id_q_obs=3330
INTEGER,PARAMETER :: id_rh_obs=3331
INTEGER,PARAMETER :: id_ps_obs=14593
INTEGER,PARAMETER :: id_rain_obs=9999
INTEGER,PARAMETER :: nid_sfcflxobs=14       !(DO_SFCFLUXES)
INTEGER,PARAMETER :: id_atm_u_obs=5287      !(OCEAN) (ATMOS)
INTEGER,PARAMETER :: id_atm_v_obs=5284      !(OCEAN) (ATMOS)
INTEGER,PARAMETER :: id_atm_t_obs=5285      !(OCEAN) (ATMOS)
INTEGER,PARAMETER :: id_atm_q_obs=5281      !(OCEAN) (ATMOS)
INTEGER,PARAMETER :: id_atm_ps_obs=5280     !(OCEAN) (ATMOS)

!ISSUE: make these namelist controllable:
INTEGER :: obs1nrec = 6                     ! The number of records in the obs1-formatted file (previous 6, 7 adds a time record).
INTEGER :: obs2nrec = 9                     ! The number of records in the obs2-formatted file (previous 8, 9 adds a time record).

! Remove all observations above 65ºN due to tripolar grid
!LOGICAL :: DO_REMOVE_65N = .true.
!STEVE: I made this a command line adjustable argument in obsop.f90 :: -rm65N .true.

LOGICAL, PARAMETER :: DO_POTTEMP_to_INSITU = .false. ! Conversion to observation space. This is needed if the
                                                     ! observations aren't converted to potential temperature
                                                     ! (as is done by most - NCEP, SODA, NASA/GMAO, etc.). But
                                                     ! unlike that approach, this does not require synthetic salinity
                                                     ! observations to be constructed from climatologies.
                                                     ! This approach is theoretically better, but investigation must
                                                     ! be done to ensure model biases to not cause significant errors.
                                                     ! (a warning from J. Carton of potential difficulty)
                                                     !
                                                     ! Only one can be true, this one takes prioirty
                                                     !
LOGICAL, PARAMETER :: DO_INSITU_to_POTTEMP = .false. ! Technically, this would require matching an observed salinity
                                                     ! measurement with each observed in situ temperature measurement
                                                     ! and using it to compute the potential temperature. The opposite
                                                     ! process is quite a bit easier.

END MODULE params_obs
