MODULE model_tools


CONTAINS

!STEVE: borrowed from mom4p1, heavily edited
!#######################################################################
! <SUBROUTINE NAME="calc_mixed_layer_depth">
!
! <DESCRIPTION>
!
! Calculate the mixed layer depth (m), which is defined as the depth ( > 0 )
! where the buoyancy difference with respect to the surface level is
! equal to buoyancy_crit (m/s2). 
!
! Note that the mixed layer depth is taken with respect to the ocean surface
! at z=eta_t, so the mixed layer depth is always positive. That is, the mld 
! is here defined as a thickness of water.
!            
! </DESCRIPTION>
!
subroutine calc_mixed_layer_depth(Thickness, salinity, theta, rho, pressure, hmxl, smooth_mld_input)

  type(ocean_thickness_type),   intent(in)  :: Thickness
  real, dimension(isd:,jsd:,:), intent(in)  :: salinity
  real, dimension(isd:,jsd:,:), intent(in)  :: theta
  real, dimension(isd:,jsd:,:), intent(in)  :: rho
  real, dimension(isd:,jsd:,:), intent(in)  :: pressure
  real, dimension(isd:,jsd:),   intent(out) :: hmxl
  logical, optional,            intent(in)  :: smooth_mld_input

  real, parameter :: epsln=1.0e-20  ! for divisions 
  integer         :: i, j, k, km1, kb
  logical         :: smooth_mld_routine

  if (.not.module_is_initialized) then
    call mpp_error(FATAL, &
    '==>Error from ocean_tracer_diag_mod (calc_mixed_layer_depth): module needs initialization ')
  endif

   if (present(smooth_mld_input)) then
    smooth_mld_routine = smooth_mld_input
  else
    smooth_mld_routine = smooth_mld
  endif

  hmxl(:,:)   = 0.0
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0

  wrk1(:,:,:) = density_delta_sfc( rho(:,:,:), salinity(:,:,:), theta(:,:,:), pressure(:,:,:))
  do k=2,nk
     do j=jsc,jec
        do i=isc,iec
           wrk2(i,j,k) = -grav*Grd%tmask(i,j,k)*wrk1(i,j,k-1)/(epsln+rho(i,j,k))
        enddo
     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        kb=Grd%kmt(i,j)
        if(kb==0) then
            hmxl(i,j) = 0.0
        else
            hmxl(i,j) = Thickness%depth_zwt(i,j,kb)
        endif
     enddo
  enddo

  do k=2,nk
     km1 = k-1
     do j=jsc,jec
        do i=isc,iec
        kb=Grd%kmt(i,j)
           if (kb == 0) then
               hmxl(i,j) = 0.0
           else
               if ( wrk2(i,j,k) >= buoyancy_crit .and. hmxl(i,j)==Thickness%depth_zwt(i,j,kb)) then
                   hmxl(i,j) = Thickness%depth_zt(i,j,km1)                                &
                             - (Thickness%depth_zt(i,j,km1) - Thickness%depth_zt(i,j,k))  &
                             * (buoyancy_crit-wrk2(i,j,km1))                              &
                             / (wrk2(i,j,k) - wrk2(i,j,km1) + epsln)
               endif
               hmxl(i,j) = hmxl(i,j) * Grd%tmask(i,j,1)
           endif
        enddo
     enddo
  enddo

 ! smooth mld
  if(smooth_mld_routine) then
      call mpp_update_domains(hmxl(:,:), Dom%domain2d)
      hmxl(:,:) = S2D(hmxl(:,:))
  endif


end subroutine calc_mixed_layer_depth
! </SUBROUTINE>  NAME="calc_mixed_layer_depth"

END MODULE model_tools
