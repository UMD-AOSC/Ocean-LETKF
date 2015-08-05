MODULE laplace_smoother
!===============================================================================
! MODULE:
!  laplace_smoother
!
! PURPOSE:
!  This module is used to compute the background error covariance spatial distribution
!
!
!
! 
!-------------------------------------------------------------------------------
! AUTHOR: Steve Penny, Dave Behringer
!===============================================================================




CONTAINS


SUBROUTINE init_weights
!===============================================================================
! SUBROUTINE:
!   init_weights
!
! PURPOSE:
!  Compute the the weights for the Laplacian smoother
!
!===============================================================================

  USE common, ONLY: re

  IMPLICIT NONE

  INTEGER       :: pe, pes, npes, len
  INTEGER       :: i, j, n, ii, jj, np, npid2
  INTEGER       :: jbg, jfn, jgbg, jgfn
  REAL          :: cuj, cujm, ctjr
  REAL          :: re, re2, con, col, acon !, wsnd(2)
  REAL(kind=8)  :: wsnd(2)
  REAL          :: b2                    !  = (hrscl * pi / 360.0)**2


  pe = mpp_pe()
  npes = mpp_npes()

  allocate (xcb(npes),xce(npes))
  allocate (ycb(npes),yce(npes))

! CALL mpp_get_compute_domains(Dom%domain2d,xcb,xce,xcsz,ycb,yce,ycsz)
  xcb = 1
  xce = nlon
  ycb = 1
  yce = nlat
  isd = 1
  ied = nlon
  jsd = 1
  jed = nlat

  aeval = 0.0
  do n=1,npits-1
    aeval=1./float(n)+aeval
  enddo
  aeval=aeval/float(npits)
  dbsq = 0.5*aeval*aeval
  
  ! Earth's radius squared
  re2 = re**2
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

  !STEVE:
  !if (dodebug) print *, "allocate (19), isd:ied,jsd:jed = ", isd,ied,jsd,jed

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
  CALL mpp_update_domains (wea(:,:), Dom%domain2d)
  CALL mpp_update_domains (wno(:,:), Dom%domain2d)
  CALL mpp_update_domains (wwe(:,:), Dom%domain2d)
  CALL mpp_update_domains (wso(:,:), Dom%domain2d)
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
      s1(:,:) = 0.0
      s2(:,:) = 0.0
      pes = -1
      do n=1,npes
        if (ii .ge. xcb(n) .and. ii .le. xce(n) .and. jj .ge. ycb(n) .and. jj .le. yce(n)) then
          pes = n - 1
        endif
      enddo
      if (pes .lt. 0) CALL mpp_error(FATAL,'ERROR in INIT_WGHTS')
      if (pe .eq. pes) then
        s1(ii,jj) = 1.0
      endif
      CALL mpp_update_domains (s1, Dom%domain2d)
      do np=1,npid2
        do j=jbg,jfn
          do i=isc,iec
            s2(i,j) = wcn(i,j) * s1(i,j) + wso(i,j) * s1(i,j-1) + wno(i,j) * s1(i,j+1) &
                                           + wwe(i,j) * s1(i-1,j) + wea(i,j) * s1(i+1,j)
          enddo
        enddo
        CALL mpp_update_domains (s2, Dom%domain2d)
        do j=jbg,jfn
          do i=isc,iec
            s1(i,j) = wcn(i,j) * s2(i,j) + wno(i,j-1) * s2(i,j-1) + wso(i,j+1) * s2(i,j+1) &
                                           + wwe(i,j) * s2(i-1,j) + wea(i,j) * s2(i+1,j)
          enddo
        enddo
        CALL mpp_update_domains (s1, Dom%domain2d)
      enddo
      wsnd = 0.0
      if (pe .eq. pes) then
        wsnd(1) = s1(ii,jj)
      endif
      len = 2
      !STEVE: wsnd == real, len == integer, pes == integer
      CALL mpp_broadcast(wsnd,len,pes)
      if (jj .ge. jbg .and. jj .le. jfn) then
        do i=isc,iec
          wgns(i,jj) = wsnd(1)
        enddo
      endif
  enddo
  CALL mpp_sync()

  do j=jbg,jfn
    do i=isc,iec
      if (wgns(i,j) .lt. 0.0) then
        wgns(i,j) = 0.00001
      endif
      wgta(i,j)=sqrt(1.0/wgns(i,j))
    enddo
  enddo
  CALL mpp_update_domains (wgta, Dom%domain2d)

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

      CALL mpp_update_domains (wea_s(:,:), Dom%domain2d)
      CALL mpp_update_domains (wno_s(:,:), Dom%domain2d)
      CALL mpp_update_domains (wwe_s(:,:), Dom%domain2d)
      CALL mpp_update_domains (wso_s(:,:), Dom%domain2d)
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
          s1(:,:) = 0.0
          s2(:,:) = 0.0
          pes = -1
          do n=1,npes
            if (ii .ge. xcb(n) .and. ii .le. xce(n) .and. jj .ge. ycb(n) .and. jj .le. yce(n)) then
              pes = n - 1
            endif
          enddo
          if (pes .lt. 0) CALL mpp_error(FATAL,'ERROR in INIT_WGHTS')
          if (pe .eq. pes) then
            s1(ii,jj) = 1.0
          endif
          CALL mpp_update_domains (s1, Dom%domain2d)
          do np=1,npid2
            do j=jbg,jfn
              do i=isc,iec
                s2(i,j) = wcn_s(i,j) * s1(i,j) + wso_s(i,j) * s1(i,j-1) + wno_s(i,j) * s1(i,j+1) &
                                               + wwe_s(i,j) * s1(i-1,j) + wea_s(i,j) * s1(i+1,j)
              enddo
            enddo
            CALL mpp_update_domains (s2, Dom%domain2d)
            do j=jbg,jfn
              do i=isc,iec
                s1(i,j) = wcn_s(i,j) * s2(i,j) + wno_s(i,j-1) * s2(i,j-1) + wso_s(i,j+1) * s2(i,j+1) &
                                               + wwe_s(i,j) * s2(i-1,j) + wea_s(i,j) * s2(i+1,j)
              enddo
            enddo
            CALL mpp_update_domains (s1, Dom%domain2d)
          enddo
          wsnd = 0.0
          if (pe .eq. pes) then
            wsnd(1) = s1(ii,jj)
          endif
          len = 2
          CALL mpp_broadcast(wsnd,len,pes)
          if (jj .ge. jbg .and. jj .le. jfn) then
            do i=isc,iec
              wgns(i,jj) = wsnd(1)
            enddo
          endif
      enddo
      CALL mpp_sync()

      do j=jbg,jfn
        do i=isc,iec
          if (wgns(i,j) .lt. 0.0) then
            wgns(i,j) = 0.00001
          endif
          wgta_s(i,j)=sqrt(1.0/wgns(i,j))
        enddo
      enddo
      CALL mpp_update_domains (wgta_s, Dom%domain2d)
    endif
  endif

 deallocate (elipt)
 deallocate(wgns)

END SUBROUTINE init_weights

!===============================================================================
! <SUBROUTINE NAME="init_grad">
!
! <DESCRIPTION>
! Compute the initial estimate of the gradient of the functional (g)
! </DESCRIPTION>
!
!===============================================================================
SUBROUTINE init_grad (Time, T_prog, Ext_mode, obs_Z, obs_0, obs_A)
!
  TYPE(ocean_time_type), intent(in)                 :: Time
  TYPE(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
  TYPE(ocean_external_mode_type), intent(in)        :: Ext_mode
  TYPE(ocean_obsz_type), intent(inout)              :: obs_Z(:)
  TYPE(ocean_obs0_type), intent(inout)              :: obs_0(:)
  TYPE(ocean_obs0_type), intent(inout)              :: obs_A(:)
!
  INTEGER         :: i, ip, j, jp, k, kks, n, taup1, pe
  INTEGER         :: year, month, day, hour, minute, second
  REAL            :: ov, aerr, aov
  TYPE(time_type) :: diff_time, wndw_fwd, wndw_bwd
  INTEGER         :: dsec, dday, dticks
  REAL            :: time_sep, time_adj
  INTEGER :: tobc, sobc
  REAL :: tsum, ssum

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
!      STEVE: this will be done externally in obsop.f90
!   2) compute innovation
!      STEVE: this can be done immediately upon reading obs2 format file
!   3) adjust error inverse assigned to observation
!       i) adjust for time separation
!      ii) adjust for obs too far from model
!   4) multiply error inverse times innovation
!   5) project result back onto model grid
!      STEVE: this may not always be possible with 3DVar (for nonlinear obs)
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
    ! The model grid coordinates containing the observation:
    i = obs_Z(n)%io
    ip = i + 1
    j = obs_Z(n)%jo
    jp = j + 1

    ! Apply the following to temperature profiles
    if (obs_Z(n)%code .eq. temp_code .and. (asm_code .eq. temp_code .or. asm_code .eq. ts_code)) then
      ! Calculate the forward and backward time windows
      ! (STEVE: this can be done externally more easily)
      wndw_fwd = increment_time(Time%model_time, wndw_secs, tz_wndw_fwd)
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, tz_wndw_bwd)

      ! If the observation is within the window,
      if (obs_Z(n)%obs_time < wndw_fwd .and. obs_Z(n)%obs_time > wndw_bwd) then
        
        ! Compute the time difference
        diff_time = Time%model_time - obs_Z(n)%obs_time
        CALL get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rtzw)
        obs_Z(n)%win = .true.
        if (save_all_inv) then
          obs_Z(n)%stat = 1
        else
          if (obs_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_Z(n)%stat = 1
        endif

        ! Loop through all levels of the observation
        do k=1,obs_Z(n)%kd

          ! Interpolate the model to the observation location
          ov = obs_Z(n)%a00 * T_prog(index_temp)%field(i,j,k,taup1) + &
               obs_Z(n)%a01 * T_prog(index_temp)%field(i,jp,k,taup1) + &
               obs_Z(n)%a11 * T_prog(index_temp)%field(ip,jp,k,taup1) + &
               obs_Z(n)%a10 * T_prog(index_temp)%field(ip,j,k,taup1)

          ! Compute the observation innovation
          ov = ov - obs_Z(n)%val(k)

          if (obs_Z(n)%stat .eq. 1) then
            obs_Z(n)%inv(k) = ov
            if (k .eq. obs_Z(n)%kd) obs_Z(n)%stat = 2
          endif

          ! Adjust the obs error based on the time
          aerr = obs_Z(n)%err(k)*time_adj

          ! aov == absolute value of the innovation
          aov = abs(ov)

          ! If the innovation is in an allowable range, then continue
          if (aov .lt. dtemp_max) then
           
            ! If the innovation is still too big, then decrease(?) the error 
            if (aov .gt. dtemp_elm) then
              aerr = aerr/(1.0+aov-dtemp_elm)**2
            endif

            ! Update the conjugate vector g by adding: the innovation * the inverse obs error weighted to each grid point
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

    ! Apply the following to salinity profiles
    else if (obs_Z(n)%code .eq. salt_code .and. (asm_code .eq. salt_code .or. asm_code .eq. ts_code)) then
      
      ! Compute forward window
      wndw_fwd = increment_time(Time%model_time, wndw_secs, sz_wndw_fwd)
      ! Compute backward window
      wndw_bwd = decrement_time(Time%model_time, wndw_secs, sz_wndw_bwd)

      ! If the ob is inside the window, then
      if (obs_Z(n)%obs_time < wndw_fwd .and. obs_Z(n)%obs_time > wndw_bwd) then

        ! Adjust the scaling for the ob
        diff_time = Time%model_time - obs_Z(n)%obs_time
        CALL get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*rszw)

        ! The ob is in the window
        obs_Z(n)%win = .true.
        if (save_all_inv) then
          obs_Z(n)%stat = 1
        else
          if (obs_Z(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_Z(n)%stat = 1
        endif

        ! for each level:
        do k=1,obs_Z(n)%kd

          ! Interpolate the model state to the obs location
          ov = obs_Z(n)%a00 * T_prog(index_salt)%field(i,j,k,taup1) + &
               obs_Z(n)%a01 * T_prog(index_salt)%field(i,jp,k,taup1) + &
               obs_Z(n)%a11 * T_prog(index_salt)%field(ip,jp,k,taup1) + &
               obs_Z(n)%a10 * T_prog(index_salt)%field(ip,j,k,taup1)

          ! Compute the innovation
          ov = ov - obs_Z(n)%val(k)
          if (obs_Z(n)%stat .eq. 1) then
            obs_Z(n)%inv(k) = ov
            if (k .eq. obs_Z(n)%kd) obs_Z(n)%stat = 2
          endif

          ! Apply the time adjustment computed above
          aerr = obs_Z(n)%err(k)*time_adj

          ! aov == the absolute value of the obs innovation
          aov = abs(ov)
          
          ! If the ob is close enough to the background
          if (aov .lt. dsalt_max) then
            ! but still not close enough, then scale the error: 
            if (aov .gt. dsalt_elm) then
              aerr = aerr/(1.0+aov-dsalt_elm)**2
            endif
            ! Next, apply the innovation scaled by the error, weighted to each model grid point
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
          CALL get_time (diff_time, dsec, dday)
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
          CALL get_time (diff_time, dsec, dday)
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
        CALL get_time (diff_time, dsec, dday)
        time_sep = real(dday) + real(dsec)/real(spd)
        time_adj = (1.0-time_sep*ralw)
        obs_A(n)%win = .true.
        if (save_all_inv) then
          obs_A(n)%stat = 1
        else
          if (obs_A(n)%obs_time <= Time%model_time .and. diff_time < gds_freq) obs_A(n)%stat = 1
        endif
        !STEVE: ISSUE: subtract the model climatology immediately when read in, not here
        ov = obs_A(n)%a00 * (Ext_mode%eta_t(i,j) - eta_clm(i,j)) + &
             obs_A(n)%a01 * (Ext_mode%eta_t(i,jp) - eta_clm(i,jp)) + &
             obs_A(n)%a11 * (Ext_mode%eta_t(ip,jp) - eta_clm(ip,jp)) + &
             obs_A(n)%a10 * (Ext_mode%eta_t(ip,j) - eta_clm(ip,j))
        ov = ov - obs_A(n)%val
        if (obs_A(n)%stat .eq. 1) then
          obs_A(n)%inv = ov
          obs_A(n)%stat = 2
        endif
        !STEVE: if aerr is made depth-dependent, then this whole section could potentially
        !       be combined with the above T/S computations, conditionally multiplying aerr
        !       with cdnz and cdnzs if it is an altimetry ob.
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
    CALL mpp_update_domains (g_cg(:,:,k), Dom%domain2d)
  enddo
!
  END SUBROUTINE init_grad


!===============================================================================
! <SUBROUTINE NAME="Bg_laplace_smoother">
!
! <DESCRIPTION>
! This SUBROUTINE multiplies g by an approximation to the first guess
! error covariance matrix [B] to get the vector h. The approximation to
! [B] is made by a series of multiplications by 1+laplacian.
! </DESCRIPTION>
!
!===============================================================================
SUBROUTINE Bg_laplace_smoother

  INTEGER         :: nit, n, i, j, k, kk, ka, kp, kkp
  INTEGER         :: np, npid2
  INTEGER         :: jbg, jfn, jgbg, jgfn
  REAL            :: con, col
  REAL :: cs
  INTEGER :: pe

  pe    = mpp_pe()

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
  do j=jsd,jed
    do i=isd,ied

      !-------------------------------------------------------------------------
      ! Do for temperature
      !-------------------------------------------------------------------------
      if (asm_code .eq. temp_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ! ka is the minimum of the assimilation depth and the water column depth
          ka = min(kass,Grd%kmt(i,j))
          ! specify the background error at each depth as the sqrt of the variance
          do k=1,ka
            ev(k) = sqrt(vtmp(i,j,k))
          enddo
        endif

        do k=1,kass
          wrkk(k) = 0.0
          do kk=1,kass
            wrkk(k) = wrkk(k) + g_cg(i,j,kk) * cvn(k,kk)
          enddo
          wrkk(k) = ev(k) * wrkk(k)
        enddo

      endif

      !-------------------------------------------------------------------------
      ! Do for salinity
      !-------------------------------------------------------------------------
      if (asm_code .eq. salt_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vsal(i,j,k))
          enddo
        endif

        do k=1,kass
          kp = k+ksalt
          wrkk(kp) = 0.0
          do kk=1,kass
            kkp = kk+ksalt
            wrkk(kp) = wrkk(kp) + g_cg(i,j,kkp) * cvn(k,kk)
          enddo
          wrkk(kp) = ev(k) * wrkk(kp)
        enddo

      endif
      
      ! Assign the working directory to the conjugate gradient vector g
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
    CALL mpp_update_domains (wcn, Dom%domain2d)

    s1(:,:) = g_cg(:,:,k) * wgta(:,:)
    s2(:,:) = 0.0

    do np=1,npid2
      do j=jbg,jfn
        do i=isc,iec
          s2(i,j) = ( wcn(i,j) * s1(i,j) + wso(i,j) * s1(i,j-1) + wno(i,j) * s1(i,j+1) &
                        + wwe(i,j) * s1(i-1,j) + wea(i,j) * s1(i+1,j) ) * Grd%tmask(i,j,kk)
        enddo
      enddo
      CALL mpp_update_domains (s2, Dom%domain2d)
      do j=jbg,jfn
        do i=isc,iec
          s1(i,j) = ( wcn(i,j) * s2(i,j) + wno(i,j-1) * s2(i,j-1) + wso(i,j+1) * s2(i,j+1) &
                        + wwe(i,j) * s2(i-1,j) + wea(i,j) * s2(i+1,j) ) * Grd%tmask(i,j,kk)
        enddo
      enddo
      CALL mpp_update_domains (s1, Dom%domain2d)
    enddo

    h_cg(:,:,k) = s1(:,:) * wgta(:,:)

  enddo
!
!-----------------------------------------------------------------------
!   multiply h by the square root of the local vertical background
!   error covariance
!-----------------------------------------------------------------------
!
  do j=jsd,jed
    do i=isd,ied

      !-------------------------------------------------------------------------
      ! Do for temperature
      !-------------------------------------------------------------------------
      if (asm_code .eq. temp_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vtmp(i,j,k))
          enddo
        endif

        do k=1,kass
          wrkk(k) = 0.0
          do kk=1,kass
            wrkk(k) = wrkk(k) + h_cg(i,j,kk) * cvn(k,kk)
          enddo
          wrkk(k) = ev(k) * wrkk(k)
        enddo

      endif

      !-------------------------------------------------------------------------
      ! Do for salinity
      !-------------------------------------------------------------------------
      if (asm_code .eq. salt_code .or. asm_code .eq. ts_code) then

        ev = 0.0
        if (Grd%kmt(i,j) .gt. 0) then
          ka = min(kass,Grd%kmt(i,j))
          do k=1,ka
            ev(k) = sqrt(vsal(i,j,k))
          enddo
        endif

        do k=1,kass
          kp = k+ksalt
          wrkk(kp) = 0.0
          do kk=1,kass
            kkp = kk+ksalt
            wrkk(kp) = wrkk(kp) + h_cg(i,j,kkp) * cvn(k,kk)
          enddo
          wrkk(kp) = ev(k) * wrkk(kp)
        enddo

      endif

      do k=1,kass2
        h_cg(i,j,k) = wrkk(k)
      enddo

    enddo
  enddo

  do k=1,kass2
    CALL mpp_update_domains (h_cg(:,:,k), Dom%domain2d)
  enddo

END SUBROUTINE Bg_laplace_smoother



END MODULE laplace smoother
