MODULE so_pcg

IMPLICIT NONE
PRIVATE


CONTAINS

SUBROUTINE setup_pcg
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

END SUBROUTINE setup_pcg


SUBROUTINE pcg (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A)
!===============================================================================
! SUBROUTINE:
!  pcg
!
! PURPOSE:
!  Preconditioned Conjugate Gradient (PCG) algorithm designed by
!  Derber and Rosati (1989), further developed by Behringer et al.
!  to expand to the global domain and include salinity and altimetry
!  observations.
!
!
! Authors: Steve Penny, Dave Behringer
!===============================================================================

  TYPE(ocean_time_type), intent(in)                    :: Time
  TYPE(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  TYPE(ocean_external_mode_type), intent(inout)        :: Ext_mode
  TYPE(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
  TYPE(ocean_obsz_type), intent(inout)                 :: obs_Z(:)
  TYPE(ocean_obs0_type), intent(inout)                 :: obs_0(:)
  TYPE(ocean_obs0_type), intent(inout)                 :: obs_A(:)

  INTEGER         :: pe, n
  INTEGER         :: i, ip, j, jp, k, iter
  REAL            :: alpha, beta, gh_old, gh_new, df_val
  INTEGER :: ii, jj

  ni    = Grd%ni
  nj    = Grd%nj
  nk    = Grd%nk
  pe = mpp_pe()

!-----------------------------------------------------------------------
!  Find the first iteration of the gradient of the functional (g^1)
!  by setting the intial guess for the correction field to zero (T^1 = 0),
!  comparing the model with the observations, weighting their difference
!  with the inverse of the observation error covariance (R) and projecting
!  this onto the model grid.
!       T^1 = 0
!       g^1 = -trnsD invR To
!-----------------------------------------------------------------------
!
  t_cg = 0.0
  CALL init_grad (Time, T_prog, Ext_mode, obs_Z, obs_0, obs_A)

!-----------------------------------------------------------------------
!  Do the first multiplication of the gradient by the background
!  error covariance matrix (B).
!       h^1 = B g^1
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------

  CALL eg_lpsmthr ()  !ISSUE: change to Bg_laplacian_smoother

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
      CALL mpp_update_domains (d_cg(:,:,k), Dom%domain2d)
      CALL mpp_update_domains (e_cg(:,:,k), Dom%domain2d)
    enddo

    !-----------------------------------------------------------------------
    !  Compute f
    !      f^n = e^n + trnsD invR D d^n
    !-----------------------------------------------------------------------
    CALL comp_f (obs_Z, obs_0, obs_A)

    !-----------------------------------------------------------------------
    !  Compute the inner products <g,h  and <d,f and update alpha
    !  (only over the computational part of the processor domain)
    !-----------------------------------------------------------------------
    !ISSUE: use the BLAS
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
      CALL mpp_update_domains (t_cg(:,:,k), Dom%domain2d)
    enddo

    if (iter .lt. maxits) then
      g_cg = g_cg + alpha * f_cg
      do k=1,kass2
        CALL mpp_update_domains (g_cg(:,:,k), Dom%domain2d)
      enddo

      !-----------------------------------------------------------------------
      !  Update h by multiplying the new gradient ( g^(n+1) ) by the
      !  background error covariance B.
      !       h^(n+1) = B g^(n+1)
      !  In this version a laplace smoother is used.
      !-----------------------------------------------------------------------
      CALL eg_lpsmthr ()  !ISSUE: again, change to Bg_laplacian_smoother

      !-----------------------------------------------------------------------
      !  Compute a new inner product <g,h and update beta
      !  (only over the computational part of the processor domain)
      !-----------------------------------------------------------------------
      !ISSUE: again, use the BLAS
      gh_old = gh_new
      gh_new = mpp_global_sum(Dom%domain2d,g_cg(:,:,:)*h_cg(:,:,:)*Grd%tmask(:,:,:),BITWISE_EXACT_SUM)
      beta = gh_new / gh_old

    endif
  enddo


  !STEVE: the following code blocks store the pcg output to FMS format data structures.
  !ISSUE: rewrite to simpler data structures
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

END SUBROUTINE pcg


SUBROUTINE finish_pcg



END SUBROUTINE finish_pcg


END MODULE so_pcg
