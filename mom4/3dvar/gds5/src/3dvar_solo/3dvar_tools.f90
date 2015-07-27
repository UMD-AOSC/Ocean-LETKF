MODULE 3dvar_tools
!===============================================================================
! MODULE: 3dvar_tools
! 
! USES:
!  use common
!  use common_mpi
!  use common_mom4
!  use common_mpi_mom4
!  use common_letkf
!  use letkf_obs !contains debug_hdxf_0
!  use letkf_local
!  use params_letkf, ONLY: nbv, cov_infl_mul, sp_infl_add, DO_INFL_RESET
!
! PUBLIC TYPES:
!                 das_3dvar
!
! SAVE:
!                 nobstotal (private)
!     
! PUBLIC MEMBER FUNCTIONS:
!           <function>                     ! Description      
!
! PUBLIC DATA MEMBERS:
!           <type> :: <variable>           ! Variable description
!
! DESCRIPTION: 
!   This module performs the main loop for the 3dvar data assimilation.
!   The 3dvar_core is called as the main routine. Note that the 3dvar
!   computation is global.
!
!   The goal here is to make 3dvar a callable routine at the end of letkf
!   execution to avoid performing the entire 3DVar-GODAS execution as an
!   additional step.
!
! REVISION HISTORY:
!   07/27/2015 Steve Penny created for use with OCEAN at NCEP.
!              based on GODAS code by D. Behringer and LETKF
!              code by T. Miyoshi.
! 
!-------------------------------------------------------------------------------
! $Authors: Steve Penny, Takemasa Miyoshi $
!===============================================================================

  USE common
  USE common_mpi
  USE common_mom4
  USE common_mpi_mom4
  USE common_letkf
  USE letkf_obs !contains debug_hdxf_0
  USE letkf_local !STEVE: separating localization functions
  USE params_letkf, ONLY: nbv, cov_infl_mul, sp_infl_add, DO_INFL_RESET 

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_3dvar

  INTEGER,SAVE :: nobstotal

CONTAINS

!-------------------------------------------------------------------------------
! Data Assimilation
!-------------------------------------------------------------------------------
SUBROUTINE das_3dvar(gues3d,gues2d,anal3d,anal2d)
!===============================================================================
! The main control structure for the data assimilation algorithm.
! Calls the 3dvar_core algorithm.
!===============================================================================
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: gues3d(nlon,nlat,nlev,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nlon,nlat,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nlon,nlat,nlev,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d)
  LOGICAL :: ex
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr
  !STEVE: for debugging
  LOGICAL :: debug_sfckmt = .false.
  LOGICAL :: dodebug = .false.
  INTEGER :: nn
  REAL(r_size) :: maxdep_val
  INTEGER :: maxdep_nn
  REAL(r_size) :: mindep_val
  INTEGER :: mindep_nn
  !STEVE: for DO_NO_VERT_LOC
  INTEGER :: klev

  WRITE(6,'(A)') 'Hello from das_3dvar'
  nobstotal = nobs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs

  !-----------------------------------------------------------------------------
  ! In case of no obs
  !-----------------------------------------------------------------------------
  if(nobstotal == 0) then
    WRITE(6,'(A)') 'No observation assimilated'
    anal3d = gues3d
    anal2d = gues2d
    RETURN ! <- return / exit the subroutine
  else                   !(OCEAN)
    anal3d = 0.0d0       !(OCEAN)
    anal2d = 0.0d0       !(OCEAN)
  endif
  
  !-----------------------------------------------------------------------------
  ! MAIN ASSIMILATION LOOP
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
! ALLOCATE(obs_useidx(1:nobs)) !STEVE: for debugging...
! obs_useidx=0
  WRITE(6,*) "... done."
  do ilev=1,nlev
    WRITE(6,'(A,I3)') 'ilev = ',ilev

    if (DO_NO_VERT_LOC .and. ilev > 1) CYCLE
          
    do ij=1,nij1 !STEVE: go through every possible coordinate of the grid in list form...
      if (dodebug) WRITE(6,*) "ij = ", ij

      !(OCEAN) STEVE: it's on land, so just assign undef values and CYCLE
      !STEVE: NEED to define kmt1 as ocean depth
      if (kmt1(ij) < ilev) then
        anal3d(ij,ilev,:,:) = 0.0
        work3d(ij,ilev,:) = 0.0
        if (ilev == 1) then
          anal2d(ij,:,:) = 0.0
          work2d(ij,:) = 0.0 !STEVE: added
        endif
        !STEVE: debug
        !WRITE(6,*) "CYCLE: kmt1(ij) = ", kmt1(ij)
        CYCLE
      endif

      !-------------------------------------------------------------------------
      ! Loop through all prognostic variables (e.g. temp, salt, u, v, etc.)
      !-------------------------------------------------------------------------
      do n=1,nv3d
        if (var_local_n2n(n) < n) then
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        else
          CALL obs_local(ij,ilev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)

          parm = work3d(ij,ilev,n)

          debug_local_obs : if ( .false. ) then !NINT(i1(ij)) .eq. 456 .and. NINT(j1(ij)) .eq. 319 .and. ilev .eq. 5) then
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "nobsl = ", nobsl
            WRITE(6,*) "i1(ij), j1(ij), ilev = ", i1(ij), j1(ij), ilev
            WRITE(6,*) "ij,n,var_local_n2n(n) = ", ij,n,var_local_n2n(n)
            maxdep_val = 0
            maxdep_nn = 0
            mindep_val = 0
            mindep_nn = 0
            do nn = 1 ,nobsl
             !if (obselm(obs_useidx(nn)) .eq. id_t_obs .OR. obselm(obs_useidx(nn)) .eq. id_s_obs .and. ALLOCATED(obs_useidx)) then
                WRITE(6,*) "---------- (letkf_tools.f90)"
                WRITE(6,*) "nn = ", nn
                if (NINT(obselm(obs_useidx(nn))) .eq. id_t_obs) WRITE(6,*) "TEMPOB"
                if (NINT(obselm(obs_useidx(nn))) .eq. id_s_obs) WRITE(6,*) "SALTOB"
                if (NINT(obselm(obs_useidx(nn))) .eq. id_sst_obs) WRITE(6,*) "SSTOB"
                WRITE(6,*) "obslon/obslat/obslev(obs_useidx(nn)) = ", obslon(obs_useidx(nn)), obslat(obs_useidx(nn)), obslev(obs_useidx(nn))
                WRITE(6,*) "obs_useidx(nn) = ", obs_useidx(nn)
                WRITE(6,*) "obselm(obs_useidx(nn)) = ", obselm(obs_useidx(nn))
                if (obsdat(obs_useidx(nn)) < 0) WRITE(6,*) "NEGOB!"
                WRITE(6,*) "obsdat(obs_useidx(nn)) = ", obsdat(obs_useidx(nn))
                WRITE(6,*) "dep(nn)    = ", dep(nn)
                WRITE(6,*) "obserr(obs_useidx(nn)) = ", obserr(obs_useidx(nn))
                WRITE(6,*) "rdiag(nn)    = ", rdiag(nn)
                WRITE(6,*) "rloc(nn)    = ", rloc(nn)
!               WRITE(6,*) "obsdep(obs_useidx(nn))    = ", obsdep(obs_useidx(nn))
                if (dep(nn) > maxdep_val) then
                  maxdep_val = dep(nn)
                  maxdep_nn = nn
                endif
                if (dep(nn) < mindep_val) then
                  mindep_val = dep(nn)
                  mindep_nn = nn
                endif
                WRITE(6,*) "hdxf(nn,:) = "
                WRITE(6,*) hdxf(nn,1:nbv)
!               WRITE(6,*) "obshdxf(obs_useidx(nn),:) = "
!               WRITE(6,*) obshdxf(obs_useidx(nn),1:nbv)
             !endif
            enddo
            WRITE(6,*) "maxdep_val = ", maxdep_val
            WRITE(6,*) "maxdep_nn  = ", maxdep_nn
            WRITE(6,*) "mindep_val = ", mindep_val
            WRITE(6,*) "mindep_nn  = ", mindep_nn
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
            WRITE(6,*) "------------------------------------------------------------"
          endif debug_local_obs

          !-------------------------------------------------------------------------
          ! Call 3DVar MAIN subroutine
          !-------------------------------------------------------------------------
          CALL 3dvar_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))   !STEVE: need to change for RIP

        endif

        !STEVE: Use the trans matrix computed in letkf_core to form the analysis
        do m=1,nbv
          anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)

          !STEVE: reset analysis to mean for all levels
          if(DO_NO_VERT_LOC .and. ilev .eq. 1) then
            do klev=2,nlev
              anal3d(ij,klev,m,n) = mean3d(ij,klev,n)
            enddo
          endif

          do k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) + gues3d(ij,ilev,k,n) * trans(k,m,n)

            if (DO_NO_VERT_LOC .and. ilev .eq. 1) then
              !STEVE: match up ij to ij at other vertical levels
              do klev=2,nlev
                if (kmt1(ij) < klev) then
                  anal3d(ij,klev,:,:) = 0.0
                  work3d(ij,klev,:) = 0.0
                else
                  anal3d(ij,klev,m,n) = anal3d(ij,klev,m,n) + gues3d(ij,klev,k,n) * trans(k,m,n)
                  work3d(ij,klev,:) = work3d(ij,ilev,:)
                endif
              enddo
            endif
          enddo

        enddo

      enddo ! n=1,nv3d

      !-------------------------------------------------------------------------
      ! Go through the 2d variables
      !-------------------------------------------------------------------------
      if (ilev == 1) then !update 2d variable at ilev=1
        do n=1,nv2d
          if (var_local_n2n(nv3d+n) <= nv3d) then
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n))
          elseif (var_local_n2n(nv3d+n) < nv3d+n) then
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)
          else
            CALL obs_local(ij,ilev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)
            parm = work2d(ij,n)
            !STEVE: check rdiag for > 0
            if (MINVAL(rdiag(1:nobsl)) .le. 0) then
              WRITE(6,*) "letkf.f90:: after obs_local, before letkf_core, for 2D, MINVAL(rdiag(1:nobsl)) ≤ 0"
              WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
            endif

            CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,nv3d+n)) !STEVE: change for RIP

            !Debugging:
!           print *, "Debugging SFC 2D adaptive inflation:"
!           print *, "pre letkf_core 2D, ilev=1: ij, n, work2d(ij,n) (parm out) = ", ij, n, work2d(ij,n)
            work2d(ij,n) = parm
!           print *, "post letkf_core 2D, ilev=1: ij, n, work2d(ij,n) (parm out) = ", ij, n, work2d(ij,n)
          endif

          !STEVE: process 2D SFC variables here:
          do m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            do k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m,nv3d+n)
            enddo
          enddo
        enddo

      endif !(ilev == 1)

    enddo !ij
  enddo !ilev

  DEALLOCATE(hdxf,rdiag,rloc,dep)
! DEALLOCATE(obs_useidx) !STEVE: for debugging...

  ! STEVE: I'd like to output the mean field in the gues spot
! gues3d(:,:,1,:) = mean3d(:,:,:)
! gues2d(:,1,:)   = mean2d(:,:)
  DEALLOCATE(mean3d,mean2d)

END SUBROUTINE das_3dvar


END MODULE 3dvar_tools
