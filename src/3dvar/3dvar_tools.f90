MODULE 3dvar_tools
!===============================================================================
! MODULE: 3dvar_tools
! 
! USES:
!  use common
!  use common_mpi
!  use common_oceanmodel
!  use common_mpi_oceanmodel
!  use common_3dvar
!  use 3dvar_obs !contains debug_hdxf_0
!  use 3dvar_local
!  use params_3dvar, ONLY: nbv, cov_infl_mul, sp_infl_add
!  use vars_3dvar
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
!   The 3dvar_core is called for each grid point. In this routine,
!   the grid points have already been distributed across processes.
!   Note that the 3dvar_core computes the transformation matrix
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
  PUBLIC ::  das_3dvar

  INTEGER,SAVE :: nobstotal

CONTAINS

!-------------------------------------------------------------------------------
! Data Assimilation
!-------------------------------------------------------------------------------

SUBROUTINE das_3dvar(gues3d,gues2d,anal3d,anal2d)
!===============================================================================
! The main control structure for the data assimilation algorithm.
! Calls the 3dvar_core algorithm for a global minimization
!===============================================================================

  USE common
  USE common_mpi
  USE common_oceanmodel
  USE common_mpi_oceanmodel
  USE common_3dvar
  USE 3dvar_obs !contains debug_hdxf_0
  USE 3dvar_local !STEVE: separating localization functions
  USE params_3dvar, ONLY: nbv, cov_infl_mul, sp_infl_add
  USE vars_3dvar,   ONLY: var_local, var_local_n2n
  USE params_obs,   ONLY: nobs
  USE params_model, ONLY: nlon, nlat, nlev, nv3d, nv2d
  USE params_model, ONLY: iv2d_mld
  USE params_model, ONLY: iv3d_t !STEVE: for debugging
  USE params_3dvar, ONLY: DO_MLD, DO_NO_VERT_LOC, DO_MLD_MAXSPRD
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
  LOGICAL :: doverbose = .false.
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
  ! Forecast perturbations
  !-----------------------------------------------------------------------------
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))

  !-----------------------------------------------------------------------------
  ! MAIN ASSIMILATION LOOP
  !-----------------------------------------------------------------------------
  WRITE(6,*) "Allocating hdxf, rdiag, rloc, and dep..."
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
! ALLOCATE(obs_useidx(1:nobs)) !STEVE: for debugging...
! obs_useidx=0
  WRITE(6,*) "... done."

  CALL 3dvar_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))

  DEALLOCATE(hdxf,rdiag,rloc,dep)


END SUBROUTINE das_3dvar

END MODULE 3dvar_tools
