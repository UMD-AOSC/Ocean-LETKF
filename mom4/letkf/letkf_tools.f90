MODULE letkf_tools
!===============================================================================
! MODULE: letkf_tools
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
!                 das_letkf, adapt_obserr, create_oer_init
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
!   This module performs the main loop for the letkf data assimilation.
!   The letkf_core is called for each grid point. In this routine,
!   the grid points have already been distributed across processes.
!   Note that the letkf_core computes the transformation matrix
!   to be applied to the background ensemble.
!
! REVISION HISTORY:
!   01/26/2009 Takemasa Miyoshi  created
!   04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!   04/03/2014 Steve Penny created for use with OCEAN at NCEP.
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
  PUBLIC ::  das_letkf, adapt_obserr, create_oer_init

  INTEGER,SAVE :: nobstotal

CONTAINS

!-------------------------------------------------------------------------------
! Data Assimilation
!-------------------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
!===============================================================================
! The main control structure for the data assimilation algorithm.
! Calls the letkf_core algorithm independently for each model grid point.
!===============================================================================
  IMPLICIT NONE
  CHARACTER(12) :: inflinfile='infl_mul.grd'
  CHARACTER(12) :: infloutfile='infl_out.grd'
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
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

  WRITE(6,'(A)') 'Hello from das_letkf'
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
  ! Variable localization
  !-----------------------------------------------------------------------------
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  WRITE(6,*) "var_local_n2n = ", var_local_n2n

  !-----------------------------------------------------------------------------
  ! Forecast perturbations
  !-----------------------------------------------------------------------------
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)

  DO n=1,nv3d
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! multiplicative inflation
  !-----------------------------------------------------------------------------
  IF(cov_infl_mul > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = cov_infl_mul
    work2d = cov_infl_mul
  END IF
  IF(cov_infl_mul <= 0.0d0) THEN ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    INQUIRE(FILE=inflinfile,EXIST=ex)
    IF(ex) THEN
      IF(myrank == 0) THEN
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflinfile
        CALL read_bingrd4(inflinfile,work3dg,work2dg)
      END IF
      CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
    ELSE
      WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflinfile
      work3d = -1.0d0 * cov_infl_mul
      work2d = -1.0d0 * cov_infl_mul
    END IF
  END IF

  !-----------------------------------------------------------------------------
  ! Reset inflation, if desired
  !-----------------------------------------------------------------------------
  if( DO_INFL_RESET ) then
    work3d = 1.0d0
    !work2d = 1.0d0 !STEVE: don't reset the 2D field
    !STEVE: will this mess up SST, SSH and SSS? maybe should reset them also.
    !       not sure if it matters since they are not prognostic in the model.

    ! Eventually, I'd like to back up the adaptive inflation from the depths to
    ! the corresponding surface level, following the water mass generation, but
    ! this would be a lot more complicated.

    ! If using the Hybrid-LETKF, multiplicative inflation is not necessary
  endif

  !-----------------------------------------------------------------------------
  ! MAIN ASSIMILATION LOOP
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
! ALLOCATE(obs_useidx(1:nobs)) !STEVE: for debugging...
! obs_useidx=0
  WRITE(6,*) "... done."
  DO ilev=1,nlev
    WRITE(6,'(A,I3)') 'ilev = ',ilev

    IF(DO_NO_VERT_LOC .and. ilev > 1) CYCLE
          
    DO ij=1,nij1 !STEVE: go through every possible coordinate of the grid in list form...
      if (dodebug) WRITE(6,*) "ij = ", ij

      !STEVE: debug
!     if (.false. .AND. ilev == 1 .AND. NINT(i1(ij)) == 131 .AND. NINT(j1(ij)) == 76) then
!       do m=1,nbv
!         WRITE(6,*) "letkf_tools.f90"
!         WRITE(6,*) "kmt1(ij) = ", kmt1(ij)
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "i1(ij) = ", i1(ij)
!         WRITE(6,*) "j1(ij) = ", j1(ij)
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!         WRITE(6,*) "should be > 20 in this region"
!       enddo
!     endif

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

      !WRITE(6,*) "ASSIM: kmt1(ij) = ", kmt1(ij)
      !(OCEAN) STEVE:end
       
!     if (ilev == 1) then 
!     do m=1,nbv
!       if (NINT(i1(ij)) == 131 .AND. NINT(j1(ij)) == 76) then
!         WRITE(6,*) "letkf_tools.f90"
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "i1(ij) = ", i1(ij)
!         WRITE(6,*) "j1(ij) = ", j1(ij)
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!         WRITE(6,*) "should be > 20 in this region"
!       endif
!       if (debug_sfckmt .AND. ABS(gues2d(ij,m,iv2d_sst)-gues3d(ij,1,m,iv3d_t)) > TINY(1.0)) then
!         WRITE(6,*) "letkf_tools.f90:: SST does not equal SFC T" 
!         WRITE(6,*) "ij, m = ", ij, m
!         WRITE(6,*) "gues2d(ij,m,iv2d_sst) = ", gues2d(ij,m,iv2d_sst)
!         WRITE(6,*) "gues3d(ij,1,m,iv3d_t) = ", gues3d(ij,1,m,iv3d_t)
!!        stop 4
!       endif
!     enddo
!     endif

      !-------------------------------------------------------------------------
      ! Loop through all prognostic variables (e.g. temp, salt, u, v, etc.)
      !-------------------------------------------------------------------------
      DO n=1,nv3d
        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        ELSE
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

          !STEVE: some critical checks to make sure letkf_core equations are valid:
          !STEVE: check to make sure input inflation parameter is valid
!         if ( isnan(parm) ) then
!           WRITE(6,*) "parm = work3d(ij,ilev,n) = ", parm
!           WRITE(6,*) "ij = ", ij
!           WRITE(6,*) "ilev = ", ilev
!           WRITE(6,*) "n = ", n
!         endif
          !STEVE: this shouldn't really happen, but fix it if it does...
          if ( parm == 0 ) then
            WRITE(6,*) "parm = work3d(ij,ilev,n) = ", parm
            WRITE(6,*) "ij = ", ij
            WRITE(6,*) "n = ", n
            WRITE(6,*) "letkf_tools.f90:: pre-letkf_core, parm changed to ABS(cov_infl_mul)"
            parm = ABS(cov_infl_mul)
          endif
          !STEVE: check rdiag for > 0
          if (MINVAL(rdiag(1:nobsl)) .le. 0) then
            WRITE(6,*) "letkf.f90:: after obs_local, before letkf_core, for 3D, MINVAL(rdiag(1:nobsl)) ≤ 0"
            WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
          endif

          !STEVE: debug
          if ( debug_hdxf_0 .AND. MINVAL(hdxf(1:nobsl,1:nbv)) == 0 ) then
            WRITE(6,*) "letkf_tools.f90:: (3D) ij = ", ij
            WRITE(6,*) "letkf_tools.f90:: inputs to letkf_core:"
            WRITE(6,*) "nobstotal = ", nobstotal
            WRITE(6,*) "nobsl = ", nobsl
            WRITE(6,*) "hdxf(1:nobsl,1:nbv) = ", hdxf(1:nobsl,:)
            WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
            WRITE(6,*) "rloc(1:nobsl) = ", rloc(1:nobsl)
            WRITE(6,*) "dep(1:nobsl) = ", dep(1:nobsl) 
            WRITE(6,*) "parm = ", parm
            WRITE(6,*) "trans(:,:,n) = ", trans(:,:,n)
          endif
          !STEVE: end

          !-------------------------------------------------------------------------
          ! Call LETKF MAIN subroutine
          !-------------------------------------------------------------------------
          CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))   !STEVE: need to change for RIP

          ! (if doing adaptive inflation)
          work3d(ij,ilev,n) = parm

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

          DO k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) + gues3d(ij,ilev,k,n) * trans(k,m,n)

            !STEVE: debug - check for bad values
!           if ( anal3d(ij,ilev,m,n) < -10 ) then
!             WRITE(6,*) "Problem in letkf_das after letkf_core. k = ", k
!             WRITE(6,*) "ij, ilev, m, n = ", ij,ilev,m,n
!             WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d(ij,ilev,m,n)
!             WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d(ij,ilev,m,n)
!             STOP 6 
!           endif

            IF(DO_NO_VERT_LOC .and. ilev .eq. 1) THEN
              !STEVE: match up ij to ij at other vertical levels
              DO klev=2,nlev
                if (kmt1(ij) < klev) then
                  anal3d(ij,klev,:,:) = 0.0
                  work3d(ij,klev,:) = 0.0
                else
                  anal3d(ij,klev,m,n) = anal3d(ij,klev,m,n) + gues3d(ij,klev,k,n) * trans(k,m,n)
                  work3d(ij,klev,:) = work3d(ij,ilev,:)
                endif
              ENDDO
            ENDIF
          END DO
          !STEVE: debug
!         if ( i1(ij) .eq. 456 .and. j1(ij) .eq. 319 .and. ilev .eq. 5) then
!           WRITE(6,*) "------------------------------------------------------------"
!           WRITE(6,*) "ij,ilev,m,n = ", ij,ilev,m,n
!           WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d(ij,ilev,m,n) + mean3d(ij,ilev,n)
!           WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d(ij,ilev,m,n)
!           WRITE(6,*) "A-B                 = ", anal3d(ij,ilev,m,n) - gues3d(ij,ilev,m,n) - mean3d(ij,ilev,n)
!           WRITE(6,*) "------------------------------------------------------------"
!         endif
        END DO

      END DO ! n=1,nv3d

      !-------------------------------------------------------------------------
      ! Go through the 2d variables
      !-------------------------------------------------------------------------
      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          IF(var_local_n2n(nv3d+n) <= nv3d) THEN
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n))
          ELSE IF(var_local_n2n(nv3d+n) < nv3d+n) THEN
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)
          ELSE
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
          END IF

          !STEVE: process 2D SFC variables here:
          DO m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            DO k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m,nv3d+n)
            END DO
          END DO
        END DO

      END IF !(ilev == 1)

    END DO !ij
  END DO !ilev

  DEALLOCATE(hdxf,rdiag,rloc,dep)
! DEALLOCATE(obs_useidx) !STEVE: for debugging...

  !-------------------------------------------------------------------------
  ! Write out the adaptive inflation
  !-------------------------------------------------------------------------
  adaptive_inflation : IF(cov_infl_mul < 0.0d0) THEN
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',infloutfile
      CALL write_bingrd4(infloutfile,work3dg,work2dg)
      !STEVE: check
!     do n=1,nv3d
!       do k=1,nlev
!         do j=1,nlat
!           do i=1,nlon
!             if (isnan4(work3dg(i,j,k,n))) then
!               WRITE(6,*) "writing work3dg(i,j,k,n) = ", work3dg(i,j,k,n)
!               WRITE(6,*) "i,j,k,n = ", i,j,k,n
!               stop 2
!             endif
!           enddo
!         enddo
!       enddo
!     enddo
      !STEVE: end
 
    END IF
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  ENDIF adaptive_inflation

  !-------------------------------------------------------------------------
  ! Compute and apply the additive inflation
  !-------------------------------------------------------------------------
  additive_inflation : IF(sp_infl_add > 0.0d0) THEN
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
      DO m=1,nbv
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
      END DO
    END DO

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,nbv
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        END DO
      END DO
    END DO
  ENDIF additive_inflation

  ! STEVE: I'd like to output the mean field in the gues spot
! gues3d(:,:,1,:) = mean3d(:,:,:)
! gues2d(:,1,:)   = mean2d(:,:)
  DEALLOCATE(mean3d,mean2d)

  RETURN
END SUBROUTINE das_letkf

SUBROUTINE adapt_obserr(hdxa,oer3d,oer2d)
!================================================================================
!STEVE: Use the following subroutines for doing adaptive observation error:
!================================================================================
  REAL(r_size), INTENT(IN) :: hdxa(nobs,nbv)
  REAL(r_size), INTENT(INOUT) :: oer3d(nij1,nlev,nv3d)
  REAL(r_size), INTENT(INOUT) :: oer2d(nij1,nv2d)
  REAL(r_size),ALLOCATABLE :: hdxb(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:), dep_ij(:)
  REAL(r_size),ALLOCATABLE :: rdiag_b(:), rdiag_a(:)
  REAL(r_size),ALLOCATABLE :: rloc_b(:), rloc_a(:)
  REAL(r_size),ALLOCATABLE :: dep_b(:), dep_a(:)
  INTEGER :: nobsl, nobsl_b, nobsl_a
  REAL(r_size), DIMENSION(nij1) :: oer, old_oer, new_oer
  INTEGER :: ij, ilev, n, j, nn, i, l
  INTEGER :: prntmod = HUGE(1)
  CHARACTER(7) :: analfile='analNNN'
  REAL(r_size) :: var_local_oer(nv3d+nv2d,nid_obs+nid_sfcflxobs) = 0.0d0 !(DO_SFCFLUXES)
!STEVE: need to use this version below to make this work. Commented out to
!support testing the assimilation of surface fluxes
! REAL(r_size),PARAMETER :: var_local_oer(nv3d+nv2d,nid_obs) = RESHAPE( &
!!           U      V      T      S    SSH    SST    SSS
!   & (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  & ! U
!   &    0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  & ! V
!   &    0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  & ! T
!   &    0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,  & ! S
!   &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,  & ! SSH
!   &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,  & ! SST
!   &    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /)& ! SSS
!   & ,(/nv3d+nv2d,nid_obs/))
  LOGICAL :: dodebug = .true.

  do i=1,nv3d+nv2d
    var_local_oer(i,i) = 1.0d0
  enddo
  if (dodebug) prntmod = NINT(nij1/3.0) !100000 !NINT(nij1/3.0)

  ! Compute per grid point (in parallel)
  WRITE(6,*) "letkf_tools::adapt_obserr: Allocating arrays..."
  WRITE(6,*) "hdxb..."
  ALLOCATE(hdxb(1:nobs,1:nbv))
  WRITE(6,*) "rdiag_b..."
  ALLOCATE(rdiag_b(1:nobs))
  WRITE(6,*) "rloc_b..."
  ALLOCATE(rloc_b(1:nobs))
  WRITE(6,*) "dep_b..."
  ALLOCATE(dep_b(1:nobs))
  WRITE(6,*) "dep_a..."
  ALLOCATE(dep_a(1:nobs)) !,rdiag_a(1:nobs),rloc_a(1:nobs),dep_a(1:nobs))
  WRITE(6,*) "obs_useidx..."
  ALLOCATE(obs_useidx(1:nobs))
  WRITE(6,*) "adapt_obserr:: ... done."

  WRITE(6,*) "------------------------------------------------------------------"
  WRITE(6,*) "letkf_tools.f90::adapt_obserr: calling obs_local and desroziers..."

  !STEVE: Just in case we didn't call das_letkf to get this initialized:
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local_oer(i,:)-var_local_oer(n,:))) < TINY(var_local_oer)) EXIT
    END DO
  END DO
  WRITE(6,*) "var_local_n2n = ", var_local_n2n

  DO ilev=1,nlev ! Cycle through vertical levels:
    WRITE(6,*) "ilev = ", ilev
    DO ij=1,nij1 ! Cycle through grid points divided up for this processor
      !(OCEAN) STEVE: it's on land, so just assign undef values and CYCLE
      !STEVE: NEED to define kmt1 as ocean depth
      if (kmt1(ij) < ilev) then
        oer3d(ij,ilev,:) = 0.0
        if (ilev == 1) oer2d(ij,:) = 0.0
        CYCLE
      endif

      if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "ij = ", ij

      DO n=1,nv3d ! Cycle through model variables: ssh, sst, sss
        if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "var_local_n2n(n), n  = ", var_local_n2n(n), n
        !STEVE: this is inefficient, and thus temporary... (don't want to call
        !this twice, and I would rather do this processing in letkf_tools)
        
        if (cnt_obs(n) < 1) CYCLE 
        CALL obs_local(ij,ilev,var_local_oer(n,:),hdxb,rdiag_b,rloc_b,dep_b,nobsl_b,nobs)
        if (nobsl_b < 1) CYCLE
          
        do nn = 1 ,nobsl_b
          dep_a(nn) = obsdat(obs_useidx(nn)) - hdxa(obs_useidx(nn),1)
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "----------"
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "nn = ", nn
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obs_useidx(nn) = ", obs_useidx(nn)
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "ilev,ij,n,var_local_n2n(n),nn = ", ilev,ij,n,var_local_n2n(n),nn
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obselm(obs_useidx(nn)) = ", obselm(obs_useidx(nn))
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obsdat(obs_useidx(nn)) = ", obsdat(obs_useidx(nn))
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "hdxa(obs_useidx(nn),1) = ", hdxa(obs_useidx(nn),1)
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "dep_a(nn) = ", dep_a(nn)
          !STEVE: or,
          !dep_a(nn) = obsdep_a(obs_useidx(nn))
        enddo

        CALL desroziers(nobsl_b, dep_a(1:nobsl_b), dep_b(1:nobsl_b), rloc_b(1:nobsl_b), &
                                                   rdiag_b(1:nobsl_b), oer3d(ij,ilev,n))
        !STEVE: debug hijack:
        !oer3d(ij,ilev,n) = nobsl_b
        if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "letkf_tools.f90:: post-desroziers, oer3d(ij,ilev,n) = ", oer3d(ij,ilev,n)
      ENDDO

      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          if (cnt_obs(nv3d+n) < 1) CYCLE 
          CALL obs_local(ij,ilev,var_local_oer(nv3d+n,:),hdxb,rdiag_b,rloc_b,dep_b,nobsl_b,nobs)
          if (nobsl_b < 1) CYCLE

          do nn = 1 ,nobsl_b
            dep_a(nn) = obsdat(obs_useidx(nn)) - hdxa(obs_useidx(nn),1)
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "----------"
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "ILEV=1 SFC, nv3d+n=",nv3d+n
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "nn = ", nn
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obs_useidx(nn) = ", obs_useidx(nn)
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "ilev,ij,n,var_local_n2n(n),nn = ", ilev,ij,n,var_local_n2n(n),nn
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obselm(obs_useidx(nn)) = ", obselm(obs_useidx(nn))
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "obsdat(obs_useidx(nn)) = ", obsdat(obs_useidx(nn))
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "hdxa(obs_useidx(nn),1) = ", hdxa(obs_useidx(nn),1)
            if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "dep_a(nn) = ", dep_a(nn)
          enddo

          CALL desroziers(nobsl_b, dep_a(1:nobsl_b), dep_b(1:nobsl_b), rloc_b(1:nobsl_b), &
                                                     rdiag_b(1:nobsl_b), oer2d(ij,n))
          if (MOD(ij,prntmod) .eq. 0) WRITE(6,*) "letkf_tools.f90:: post-desroziers, oer2d(ij,n) = ", oer2d(ij,n)

        ENDDO
      ENDIF

    ENDDO
  ENDDO

  WRITE(6,*) "letkf_tools::adapt_obserr: Deallocating..."
  WRITE(6,*) "hdxb..."
  DEALLOCATE(hdxb)
  WRITE(6,*) "rdiag_b..."
  DEALLOCATE(rdiag_b)
  WRITE(6,*) "rloc_b..."
  DEALLOCATE(rloc_b)
  WRITE(6,*) "dep_b..."
  DEALLOCATE(dep_b)
  WRITE(6,*) "dep_a..."
  DEALLOCATE(dep_a)
  WRITE(6,*) "obs_useidx..."
  DEALLOCATE(obs_useidx)
  WRITE(6,*) "... done."

END SUBROUTINE adapt_obserr

SUBROUTINE desroziers(nobsl,dep_a,dep_b,rloc,rdiag,oer)
!================================================================================
! The Desroziers-type (Desroziers, 2005) statistics for the O-F, O-A, and A-B
!================================================================================
INTEGER, INTENT(IN) :: nobsl
REAL(r_size), DIMENSION(nobsl), INTENT(IN) :: dep_a, dep_b, rloc, rdiag
REAL(r_size), INTENT(INOUT) :: oer
REAL(r_size), DIMENSION(nobsl) :: dep
REAL(r_size) :: old_oer, new_oer, rlsum
INTEGER :: j
REAL(r_size) :: gain = 0.01 ! (default)
LOGICAL, SAVE :: dodebug = .true.

  if ( nobsl < 1 ) then
    RETURN 
  endif
  ! Watch out for under-sampled regions
! if (nobsl < 100) RETURN

  do j = 1,nobsl
    dep(j) = dep_a(j) * dep_b(j)
  enddo

  ! Average old error around this grid point:
  !STEVE: at some point, it would be interesting to
  !       use some technique to compare to this value
  !       as well. Maybe (but it's not changing with time)
! old_oer = oer 
! old_oer = SUM(rdiag(1:nobsl)) / nobsl
  rlsum = SUM(rloc(1:nobsl))
  old_oer = SUM(rdiag(1:nobsl)*(rloc(1:nobsl))) / rlsum

  ! New error at this grid point:
! new_oer = SUM(dep(1:nobsl)) / nobsl
  new_oer = SUM(  dep(1:nobsl)*(rloc(1:nobsl))) / rlsum

  ! Scale the gain proportional to the number of samples
  ! (i.e. # of local obs / # of global obs)
  gain = 0.01

  if (dodebug) then
    WRITE(6,*) "nobsl = ", nobsl
    WRITE(6,*) "nobs  = ", nobs
    WRITE(6,*) "dep_a = ", dep_a
    WRITE(6,*) "dep_b = ", dep_b
    WRITE(6,*) "rloc  = ", rloc
    WRITE(6,*) "rdiag = ", rdiag
    WRITE(6,*) "oer   = ", oer
    WRITE(6,*) "old_oer = ", old_oer
    WRITE(6,*) "new_oer = ", new_oer
    WRITE(6,*) "gain  = ", gain
  endif

  ! Smoothed error adjustment: (STEVE: change to kalman filter update by
  ! using the ensemble information in hdxb & hdxa.)
  oer = SQRT(new_oer * gain + (1 - gain) * old_oer)

  if (dodebug) WRITE(6,*) "oer (final) = ", oer
  dodebug = .false.

END SUBROUTINE desroziers

SUBROUTINE create_oer_init(infile,oer3dg,oer2dg)
!===============================================================================
! Initialize the observation error adaptive estimation
!===============================================================================
CHARACTER(*), INTENT(IN) :: infile
REAL(r_sngl), INTENT(OUT) :: oer3dg(nlon,nlat,nlev,nv3d)
REAL(r_sngl), INTENT(OUT) :: oer2dg(nlon,nlat,nv2d)
INTEGER :: k

! ERROR PROFILES: Taken from a high-resolution SODA analysis
REAL, DIMENSION(nlev), PARAMETER :: t_eprof = (/0.56, 0.61, 0.66, 0.70, 0.74, 0.78, 0.80, &
0.82, 0.82, 0.84, 0.83, 0.82, 0.81, 0.80, 0.78, 0.76, 0.74, 0.72, 0.70, 0.68, &
0.67, 0.66, 0.65, 0.63, 0.60, 0.57, 0.54, 0.49, 0.43, 0.36, 0.31, 0.27, 0.27, &
0.27, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27/)

REAL, DIMENSION(nlev), PARAMETER :: s_eprof = (/1.05, 0.98, 0.90, 0.84, 0.79, 0.76, 0.74, &
0.73, 0.72, 0.70, 0.69, 0.67, 0.65, 0.64, 0.61, 0.57, 0.54, 0.52, 0.50, 0.48, &
0.47, 0.46, 0.44, 0.42, 0.39, 0.36, 0.34, 0.31, 0.27, 0.22, 0.19, 0.16, 0.16, &
0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16/)

! Create initial observation error profile for each ob type
oer3dg = 0.0
do k=1,nlev
  where(kmt0 .ge. k)
    oer3dg(:,:,k,iv3d_u) = 1.0 !m/s
    oer3dg(:,:,k,iv3d_v) = 1.0 !m/s
    oer3dg(:,:,k,iv3d_t) = t_eprof(k) 
    oer3dg(:,:,k,iv3d_s) = s_eprof(k)
  end where
enddo

oer2dg = 0.0
where(kmt0 .gt. 0)
  oer2dg(:,:,iv2d_sst) = 2.0 !deg C
  oer2dg(:,:,iv2d_sss) = 1.0 !psu
  oer2dg(:,:,iv2d_ssh) = 0.1 !meters
end where

END SUBROUTINE create_oer_init

END MODULE letkf_tools
