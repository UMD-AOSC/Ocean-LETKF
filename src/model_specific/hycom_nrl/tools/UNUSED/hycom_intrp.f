      SUBROUTINE hycom_intrp(var_3d,ri,     
     &    rj,zz,utilz,obs_id_hycom)

      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80

      USE params_model, ONLY: nlon, nlat, nlev 
      USE vars_model, ONLY: lat2d, phi0 




c
c --- hycom/micom to 3-d z-level diagnostic field extractor
c
      real(kind(0.0d0)) :: var_3d(nlon,nlat,nlev,5)
      real,    allocatable, dimension (:)     ::
     &   ttk,ssk
      real,    allocatable, dimension (:,:,:)  :: p,dp,utilk
      real,    allocatable, dimension (:,:,:,:)  :: var_hycom 
      real(kind(0.0d0)) ::  utilz(2,2)
      real ::  utilz_temp(2,2)

c
      real      onecm
      data      onecm/.01/
      real(kind(0.0d0))      zz,ri,rj
      real zz_temp
      real, parameter :: flag = 2.0**100
      integer i,j,k,ii,jj,itype
      integer obs_id_hycom
      real dpth,pmid,phi,plo,sum_tem,bot


c
      ii = 2
      jj = 2
      itype = 1
      dpth=0.5*onecm
      bot = 0.0
      utilz=0.0
      utilz_temp=0.0
      zz_temp=zz
      allocate (p(ii,jj,nlev+1))
      allocate (dp(ii,jj,nlev))
      allocate (var_hycom(ii,jj,nlev,4))
      allocate(  utilk(ii,jj,nlev+1) )




c     HYCOM input
      dp(:,:,:) = var_3d (int(ri):int(ri)+1,int(rj):int(rj)+1,:,5)
      var_hycom  = var_3d (int(ri):int(ri)+1,int(rj):int(rj)+1,:,1:4)

      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo

      do 3 k=1,nlev
      do 3 j=1,jj
      do 3 i=1,ii

c
c --- convert layer thickness to meters
      if (phi0(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
      else
        var_hycom(i,j,k,:)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
      endif
 3    continue
 


c
c --- put vertically averaged values into massless layers
c
      allocate( ttk(nlev) )
c
      do 70 n=1,4
      do 70 j=1,jj
      do 70 i=1,ii
c
      if (phi0(i,j).gt.0.) then
        do k= 1,nlev
          ttk(k)=0.
          pmid=.5*(p(i,j,k)+p(i,j,k+1))
          phi=pmid+dpth
          plo=pmid-dpth
c
          sum_tem=0.
          do k1=1,nlev
            delp=max(0.,min(p(i,j,k1+1),phi)-max(p(i,j,k1),plo))
            sum_tem=sum_tem+delp

            ttk(k)=ttk(k)+var_hycom(i,j,k1,n)*delp
          enddo !k1
c
          ttk(k)=ttk(k)/sum_tem
        enddo !k
        do k= 1,nlev
          var_hycom(i,j,k,n)=ttk(k)
        enddo !k
      end if !ip
 70   continue


c --- -------------------
c --- u-velocity
c --- -------------------
c
      if (obs_id_hycom.eq.1) then
        do k= 1,nlev
          do j=1,jj
            do i=1,ii
              if (phi0(i,j).gt.0.) then  !may be approximate at i=1 and i=ii    
                if     (dp(i,j,k).gt.0.1 .and.
     &                  p(i,j,k).lt.p(i,j,nlev+1)-bot) then
                    utilk(i,j,k)=var_hycom(i,  j,k,1)
                else
                  utilk(i,j,k)=utilk(i,j,k-1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
        enddo
        call layer2z(utilk,p,utilz_temp,zz_temp,flag,ii,jj,nlev,1,itype)
      endif !u-velocity

c
c --- -------------------
c --- v-velocity
c --- -------------------
c
      if (obs_id_hycom.eq.2) then
        do k= 1,nlev
          do j=1,jj
            do i=1,ii
              if (phi0(i,j).gt.0.) then  !may be approximate at j=1 and j=jj      
                if     (dp(i,j,k).gt.0.1 .and.
     &                  p(i,j,k).lt.p(i,j,nlev+1)-bot) then
                    utilk(i,j,k)=var_hycom(i,  j,k,2)
                else
                  utilk(i,j,k)=utilk(i,j,k-1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
        enddo
        call layer2z(utilk,p,utilz_temp,zz_temp,flag,ii,jj,nlev,1,itype)
      endif !v-velocity

c
c --- ----------------
c --- temperature
c --- ----------------
c
      if (obs_id_hycom.eq.3) then
        do k= 1,nlev
          do j=1,jj
            do i=1,ii
              if (phi0(i,j).gt.0.) then
                if     (p(i,j,k).lt.p(i,j,nlev+1)-bot) then
                  utilk(i,j,k)=var_hycom(i,j,k,3)
c Put pt<->t insitu conversion to later in obsop
c                  dbar =
c     &              PPSW_p80(0.5*(p(i,j,k)+p(i,j,k+1)), lat2d(i,j))
c                  utilk(i,j,k)=
c     &              PPSW_theta(var_hycom(i,j,k,4),
c     &              var_hycom(i,j,k,3), 0.0,dbar)
                else
                  utilk(i,j,k)=utilk(i,j,k-1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
        enddo
        call layer2z(utilk,p,utilz_temp,zz_temp,flag,ii,jj,nlev,1,itype)

      endif !temperature

c
c --- -------------
c --- salinity
c --- -------------
c
      if (obs_id_hycom.eq.4) then
        do k= 1,nlev
          do j=1,jj
            do i=1,ii
              if (phi0(i,j).gt.0.) then
                if     (p(i,j,k).lt.p(i,j,nlev+1)-bot) then
                  utilk(i,j,k)=var_hycom(i,j,k,4)
                else
                  utilk(i,j,k)=utilk(i,j,k-1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo
          enddo
        enddo
        call layer2z(utilk,p,utilz_temp,zz_temp,flag,ii,jj,nlev,1,itype)
      endif


      utilz=utilz_temp

      deallocate (p)
      deallocate (dp)
      deallocate (var_hycom)
      deallocate(  utilk)
      deallocate(  ttk )



      end

