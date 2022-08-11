MODULE vars_letkf
  USE common,       ONLY: r_size
  USE params_obs,   ONLY: nid_obs
  USE params_model, ONLY: nv3d, nv2d, nv4d

  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d+nv4d)
  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d+nv4d,nid_obs) = 1.0d0

!  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d+nv4d,nid_obs) = RESHAPE( &      !(OCEAN) !(DO_SFCFLUXES)
!!           U      V      T      S    SSH    SST    SSS     uflx   vflx  tflux  qflux    u10    v10    t2m    q2m PRESmsl prate    dlw    dsw  longwv shrtwv
!   & (/ 1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! U !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, & ! V !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! T !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! S !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SSH !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, & ! SST !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /) & ! SSS !(OCEAN)
!  & ,(/nv3d+nv2d+nv4d,nid_obs/))
! STEVE: note: this was a sort of 'trick' to make sure that the SFC variables
! were analyzed separately. All that was needed for this to occur was to have
! the var_local vector associated with each SFC variable to be different than
! the vector associated with the ocean interior variables. !(DO_SFCFLUXES)
! As a result, I can compute adaptive inflation separately for the SFC (and
! apply none to the ocean interior).

!  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d+nv4d,nid_obs) = RESHAPE( &      !(OCEAN)
!!           U      V      T      S    SSH    SST    SSS                      !(OCEAN)
!   & (/ 1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! U             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! V             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! T             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! S             !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! SSH           !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0,  & ! SST           !(OCEAN)
!   &    1.0d0, 1.0d0, 1.0d0, 1.00d0, 1.0d0, 1.0d0, 1.0d0 /)& ! SSS           !(OCEAN)
!   & ,(/nv3d+nv2d+nv4d,nid_obs/))
   !NOTE: the obs are the rows and the model variables the columns


  CONTAINS

  !STEVE: I think this would be better in "vars_obs.f90", maybe. (ISSUE)
  SUBROUTINE get_iobs(oelm,iobs)
    USE params_obs, ONLY: id_u_obs, id_v_obs, id_t_obs, id_s_obs, &
                          id_ssh_obs, id_ssh_obs, id_sst_obs, id_sss_obs, id_eta_obs
    REAL(r_size), INTENT(IN) :: oelm
    INTEGER, INTENT(OUT) :: iobs
    !---------------------------------------------------------------------------
    ! Used for variable localization
    !---------------------------------------------------------------------------
    SELECT CASE(NINT(oelm))
    CASE(id_u_obs)
        iobs=1
    CASE(id_v_obs)
        iobs=2
    CASE(id_t_obs)
        iobs=3
    CASE(id_s_obs)   !(OCEAN)
        iobs=4
    CASE(id_ssh_obs) !(OCEAN)
        iobs=5
    CASE(id_sst_obs) !(OCEAN)
        iobs=6
    CASE(id_sss_obs) !(OCEAN)
        iobs=7
    CASE(id_eta_obs) !(OCEAN)
        iobs=8
    CASE DEFAULT
        WRITE(6,*) "vars_letkf.f90 :: there is no variable localization for obs-type :: ", oelm
        WRITE(6,*) "vars_letkf.f90 :: FATAL ERROR, exiting..."
        STOP (95)
    END SELECT
  END SUBROUTINE get_iobs

END MODULE vars_letkf
