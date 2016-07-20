MODULE letkf_tools
!===============================================================================
! MODULE: letkf_tools
! 
! USES:
!  use common
!  use common_mpi
!  use common_oceanmodel
!  use common_mpi_oceanmodel
!  use common_letkf
!  use letkf_obs !contains debug_hdxf_0
!  use letkf_local
!  use params_letkf, ONLY: nbv, cov_infl_mul, sp_infl_add
!  use vars_letkf
!
! PUBLIC TYPES:
!                 das_letkf
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


  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

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

  USE common
  USE common_mpi
  USE common_oceanmodel
  USE common_mpi_oceanmodel
  USE common_letkf
  USE letkf_obs !contains debug_hdxf_0
  USE letkf_local !STEVE: separating localization functions
  USE params_letkf, ONLY: nbv, cov_infl_mul, sp_infl_add
  USE vars_letkf,   ONLY: var_local, var_local_n2n
  USE params_obs,   ONLY: nobs
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: iv2d_mld
  USE params_model, ONLY: iv3d_t !STEVE: for debugging
  USE params_letkf, ONLY: DO_MLD, DO_NO_VERT_LOC, DO_MLD_MAXSPRD
  USE vars_model,   ONLY: lev

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
  INTEGER :: ij,n,m,i,j,k,nobsl,ierr
  INTEGER :: ilev
  !STEVE: for debugging
  LOGICAL :: debug_sfckmt = .false.
  LOGICAL :: dodebug = .true.
  LOGICAL :: debug_zeroinc = .false.
  INTEGER :: nn
  REAL(r_size) :: maxdep_val
  INTEGER :: maxdep_nn
  REAL(r_size) :: mindep_val
  INTEGER :: mindep_nn
  !STEVE: for DO_NO_VERT_LOC
  INTEGER :: klev, mlev !mlev is the level for the mixed layer depth
  !STEVE: for DO_MLD
  REAL(r_size) :: max_sprd, mld_sprd
  REAL(r_size) :: mld=0

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
  do n=2,nv3d+nv2d
    do i=1,n
      var_local_n2n(n) = i
      if(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) exit
    enddo
  enddo
  WRITE(6,*) "var_local_n2n = ", var_local_n2n

  !-----------------------------------------------------------------------------
  ! Forecast perturbations
  !-----------------------------------------------------------------------------
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)

  do n=1,nv3d
    do m=1,nbv
      do k=1,nlev
        do i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        enddo
      enddo
    enddo
  enddo
  do n=1,nv2d
    do m=1,nbv
      do i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      enddo
    enddo
  enddo

  !STEVE: Check to make sure user supplied a set of non-identical ensemble members:
  if (dodebug) then
    WRITE(6,*) "letkf_tools.f90:: After computing perturbations..."
    WRITE(6,*) "min val for level gues3d(:,1,:,iv3d_t) = ", MINVAL(gues3d(:,1,:,iv3d_t))
    WRITE(6,*) "max val for level gues3d(:,1,:,iv3d_t) = ", MAXVAL(gues3d(:,1,:,iv3d_t))
  endif
  if (MAXVAL(gues3d(:,1,:,iv3d_t)) == MINVAL(gues3d(:,1,:,iv3d_t))) then
    WRITE(6,*) "letkf_tools.f90::das_letkf:: It appears that all the ensemble members are identical. EXITING..."
    STOP(833)
  endif

  !-----------------------------------------------------------------------------
  ! multiplicative inflation
  !-----------------------------------------------------------------------------
  ALLOCATE( work3d(nij1,nlev,nv3d) )
  ALLOCATE( work2d(nij1,nv2d) )
  if (cov_infl_mul > 0.0d0) then ! fixed multiplicative inflation parameter
    work3d = cov_infl_mul
    work2d = cov_infl_mul
  endif
  if (cov_infl_mul < 0.0d0) then ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    INQUIRE(FILE=inflinfile,EXIST=ex)
    if (ex) then
      if (myrank == 0) then
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflinfile
!       CALL read_grd4(inflinfile,work3dg,work2dg)
        CALL read_grd(inflinfile,work3dg,work2dg)
      endif
      CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
    else
      WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflinfile
      work3d = -1.0d0 * cov_infl_mul
      work2d = -1.0d0 * cov_infl_mul
    endif
  endif

  !-----------------------------------------------------------------------------
  ! MAIN ASSIMILATION LOOP
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
! ALLOCATE(obs_useidx(1:nobs)) !STEVE: for debugging...
! obs_useidx=0
  WRITE(6,*) "... done."

  do ij=1,nij1 !STEVE: go through every possible coordinate of the grid in list form...
               !NOTE: I switched the loops for ij and ilev (below) based on an indication
               !      by T. Sluka that this improved performance due to caching issues (3/22/16)
    if (dodebug) WRITE(6,*) "ij = ", ij

    !(OCEAN) The gridpoint is on land, so just assign undef values and CYCLE
    if (kmt1(ij) < 1) then
      anal3d(ij,:,:,:) = 0.0
      work3d(ij,:,:) = 0.0
      anal2d(ij,:,:) = 0.0
      work2d(ij,:) = 0.0
      CYCLE
    endif

    ! Initialize the mixed layer depth at level 1 (STEVE: could initialize at background ensemble mean)
    mlev=1
    if (DO_MLD) then
      if (DO_MLD_MAXSPRD) then
        max_sprd=0
        do k=1,nlev
!         mld_sprd = SQRT(SUM(gues3d(ij,k,:,iv3d_t)**2)/(nbv-1))
          mld_sprd = SUM(gues3d(ij,k,:,iv3d_t)**2) !STEVE: since we're just finding the max, sqrt and division are not really necessary.
          if (mld_sprd > max_sprd) then
            ! find the max spread in the column, use this as the 'mixed layer depth' for SST assimilation
            max_sprd = mld_sprd
            mlev = k
          endif
        enddo
      else
        ! Initialize the mixed layer depth at the mean MLD of the background ensemble
        ! Update the model-derived mixed layer depth
        mld = mean2d(ij,iv2d_mld)
        do k=1,nlev-1
          if (lev(k+1) > mld) then
            mlev = k
            exit
          endif
        enddo
        !NOTE: this will be updated after the analysis at the first level of this column
      endif
    endif

    ilev=0
    do while (ilev < nlev) !STEVE: cycle through each level.
      ilev=ilev+1
!     if (dodebug) WRITE(6,'(A,I3)') 'ilev = ',ilev

      !(OCEAN) STEVE: it's on land, so just assign undef values and CYCLE
      !NOTE: need to define kmt1 as ocean depth during startup/initialization
      if (kmt1(ij) < ilev) then
        anal3d(ij,ilev:nlev,:,:) = 0.0
        work3d(ij,ilev:nlev,:) = 0.0
        if (ilev == 1) then
          anal2d(ij,:,:) = 0.0
          work2d(ij,:) = 0.0 !STEVE: added
        endif
        !STEVE: debug
        !WRITE(6,*) "CYCLE: kmt1(ij) = ", kmt1(ij)
        EXIT !STEVE: skip all deeper levels
      endif
      if (dodebug) WRITE(6,*) "ilev = ", ilev
      if (dodebug) WRITE(6,*) "mlev = ", mlev

      !-------------------------------------------------------------------------
      ! Loop through all prognostic variables (e.g. temp, salt, u, v, etc.)
      !-------------------------------------------------------------------------
      if (dodebug) WRITE(6,*) "Doing 3D variables..."
      do n=1,nv3d
        if (var_local_n2n(n) < n) then
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))  !(inflation)
        else
          !STEVE: need to change localization for obs below mlev - added to input arguments (3/25/2016)
          CALL obs_local(ij,ilev,mlev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)
          if (dodebug) WRITE(6,*) "letkf_tools.f90::post-obs_local(3d):: Assimilating ", nobsl, " observations."

          parm = work3d(ij,ilev,n)

          !STEVE: debug
          if ( debug_hdxf_0 .and. MINVAL(hdxf(1:nobsl,1:nbv)) == 0 ) then
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
            STOP(630)
          endif
          !STEVE: end

          !STEVE: debug
          if ( MINVAL(rdiag(1:nobsl)) .le. 0.0 ) then
            WRITE(6,*) "common_letkf.f90:: ERROR: rdiag <=0 (i.e. there is an obserr <= 0)"
            WRITE(6,*) "nbv = ", nbv
            WRITE(6,*) "MINVAL(rdiag) = ",MINVAL(rdiag)
            STOP(1)
          endif

          !-------------------------------------------------------------------------
          ! Call LETKF MAIN subroutine
          !-------------------------------------------------------------------------
!         if (dodebug) WRITE(6,*) "Calling letkf_core..."
          CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))

          if ( debug_zeroinc .and. nobsl > 0 ) then 
            WRITE(6,*) "letkf_tools.f90::das_letkf:: post-letkf_core"
            WRITE(6,*) "n = ", n
            WRITE(6,*) "trans(:,:,n) = ", trans(:,:,n)
            WRITE(6,*) "letkf_tools.f90::das_letkf:: EXITING on purpose..."
            STOP(318)
          endif

          ! (if doing adaptive inflation) 
          work3d(ij,ilev,n) = parm

        endif

        !STEVE: Use the trans matrix computed in letkf_core to form the analysis ensemble
!       if (dodebug) WRITE(6,*) "Setting anal3d..."
        do m=1,nbv
          anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
          do k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) + gues3d(ij,ilev,k,n) * trans(k,m,n)

            !STEVE: debug - check for bad values
            if ( anal3d(ij,ilev,m,n) < -10 ) then
              WRITE(6,*) "Problem in letkf_das after letkf_core. k = ", k
              WRITE(6,*) "ij, ilev, m, n = ", ij,ilev,m,n
              WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d(ij,ilev,m,n)
              WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d(ij,ilev,m,n)
              !ij = 895
              WRITE(6,*) "lon = ", lon1(ij)
              WRITE(6,*) "lat = ", lat1(ij)
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
              STOP(6)
            endif
          enddo
        enddo

      enddo ! n=1,nv3d

      !-------------------------------------------------------------------------
      ! Go through the 2d variables
      !-------------------------------------------------------------------------
      if (ilev == 1) then !update 2d variable at ilev=1
        if (dodebug) WRITE(6,*) "Doing 2D variables..."
        do n=1,nv2d
          if (var_local_n2n(nv3d+n) <= nv3d) then
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n))      !(inflation)
          elseif (var_local_n2n(nv3d+n) < nv3d+n) then
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d) !(inflation)
          else
            CALL obs_local(ij,ilev,mlev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)
            if (dodebug) WRITE(6,*) "letkf_tools.f90::post-obs_local(2d):: Assimilating ", nobsl, " observations."

            parm = work2d(ij,n)
            !STEVE: check rdiag for > 0
            if (MINVAL(rdiag(1:nobsl)) .le. 0) then
              WRITE(6,*) "letkf.f90:: after obs_local, before letkf_core, for 2D, MINVAL(rdiag(1:nobsl)) <= 0"
              WRITE(6,*) "rdiag(1:nobsl) = ", rdiag(1:nobsl)
              STOP(43)
            endif

            CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,nv3d+n))

            work2d(ij,n) = parm
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

      if (ilev==1 .and. DO_MLD .and. .not. DO_MLD_MAXSPRD) then
        !If analyzing the mixed layer depth, update the layer associated with the mld here (was initialized at mlev=1 for this column)
        do k=1,nlev-1
          mld = anal2d(ij,k,iv2d_mld)
          if (lev(k+1) > mld) then
            mlev = k
            exit
          endif
        enddo
      endif

      if (DO_NO_VERT_LOC .or. DO_MLD) then
        if (dodebug) WRITE(6,*) "===================================================================="
        if (dodebug) WRITE(6,*) "Doing DO_NO_VERT_LOC projection to lower levels..."
        if (dodebug) WRITE(6,*) "===================================================================="
        !STEVE: Assimilating SSTs: I want to assimilate ONLY profiles below the mixed layer,
        !       and BOTH profiles and SSTs in the mixed layer.
        ! If we are below the mixed layer depth, remove all SST data, do one more level, then
        ! cycle and apply identical weights to the rest of the depths.
        ! This should allow SSTs to influence ONLY the mixed layer, while profiles influence all levels.

        ! After assimilating the mld in v2d, update the mld based on the analysis ensemble mean
        !STEVE: match up ij to ij at other vertical levels
        !STEVE: reset analysis to mean for all levels

        if (DO_NO_VERT_LOC .and. .not. DO_MLD) mlev=nlev
        if (DO_MLD .and. ilev > mlev) mlev=nlev

        do n=1,nv3d
!         if (dodebug) WRITE(6,*) "n = ", n
          do klev=ilev+1,mlev
!           if (dodebug) WRITE(6,*) "set mean: klev = ", klev
            anal3d(ij,klev,:,n) = mean3d(ij,klev,n)
          enddo

          do m=1,nbv
!           if (dodebug) WRITE(6,*) "row member m = ", m
            !-------------------------------------------
            ! Update the whole layer with these weights:
            !-------------------------------------------
            do k=1,nbv
!             if (dodebug) WRITE(6,*) "col member: k = ", k
              do klev=ilev+1,mlev
!              if (dodebug) WRITE(6,*) "set members: k, klev = ", k,klev
                if (klev <= kmt1(ij)) then !STEVE: apply the weights from this level to all deeper levels equally
!                 if (dodebug) WRITE(6,*) "applying weights..."
                  anal3d(ij,klev,m,n) = anal3d(ij,klev,m,n) + gues3d(ij,klev,k,n) * trans(k,m,n)
                  work3d(ij,klev,:) = work3d(ij,ilev,:)   !adaptive inflation
                else !STEVE: this depth is beyond the ocean floor
!                 if (dodebug) WRITE(6,*) "hit land."
                  anal3d(ij,klev:nlev,m,n) = 0.0
                  work3d(ij,klev:nlev,n) = 0.0
                  EXIT
                endif
              enddo !klev
            enddo !k

          enddo !m
        enddo !n
        ilev=mlev
        if (dodebug) WRITE(6,*) "final mlev=ilev = ", ilev
        if (dodebug) WRITE(6,*) "===================================================================="
      endif !DO_MLD

    enddo !ilev
  enddo !ij

  DEALLOCATE(hdxf,rdiag,rloc,dep)
! DEALLOCATE(obs_useidx) !STEVE: for debugging...

  !-------------------------------------------------------------------------
  ! Write out the adaptive inflation
  !-------------------------------------------------------------------------
  adaptive_inflation : if (cov_infl_mul < 0.0d0) then
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    if (myrank == 0) then
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',infloutfile
!     CALL write_grd4(infloutfile,work3dg,work2dg)
      CALL write_grd(infloutfile,work3dg,work2dg)
    endif
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  endif adaptive_inflation

  !-------------------------------------------------------------------------
  ! Compute and apply the additive inflation
  !-------------------------------------------------------------------------
  additive_inflation : if (sp_infl_add > 0.0d0) then
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    do n=1,nv3d
      do m=1,nbv
        do k=1,nlev
          do i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          enddo
        enddo
      enddo
    enddo
    do n=1,nv2d
      do m=1,nbv
        do i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        enddo
      enddo
    enddo

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(6,'(A)') '========================================='
    do n=1,nv3d
      do m=1,nbv
        do ilev=1,nlev
          do ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          enddo
        enddo
      enddo
    enddo
    do n=1,nv2d
      do m=1,nbv
        do ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        enddo
      enddo
    enddo
  endif additive_inflation

  ! STEVE: I'd like to output the mean field in the gues spot
! gues3d(:,:,1,:) = mean3d(:,:,:)
! gues2d(:,1,:)   = mean2d(:,:)
  DEALLOCATE(mean3d,mean2d)

END SUBROUTINE das_letkf

END MODULE letkf_tools
