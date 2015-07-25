MODULE params_obs

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

INTEGER,SAVE :: nobs
!STEVE: making these namelist accessible:
INTEGER,SAVE :: nslots=5                  ! number of time slots for 4D-LETKF
INTEGER,SAVE :: nbslot=5 !1               !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
REAL(r_size),SAVE :: sigma_obs=720.0d3    !3x Rossby Radius of Deformation at Equ. according to Chelton
REAL(r_size),SAVE :: sigma_obs0=200.0d3   !20x Rossby Radius of Deformation at Pole, according to Chelton
REAL(r_size),SAVE :: sigma_obsv=1000.0d0  !STEVE: doesn't matter if using option "DO_NO_VERT_LOC"
REAL(r_size),SAVE :: sigma_obst=5.0d0     ! Not using this at the moment
REAL(r_size),SAVE :: gross_error=3.0d0    ! number of standard deviations
! REAL(r_size),PARAMETER :: gross_error=10.0d0 !3.0d0 ! number of standard deviations   (Use for OSSEs)
                                                      ! used to filter out observations
INTEGER,PARAMETER :: nid_obs=7
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

INTEGER :: obs1nrec = 6                     ! The number of records in the obs1-formatted file (previous 6, 7 adds a time record).
INTEGER :: obs2nrec = 9                     ! The number of records in the obs2-formatted file (previous 8, 9 adds a time record).

END MODULE params_obs
