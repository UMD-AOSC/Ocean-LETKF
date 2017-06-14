module common_debug_oceanmodel
! This module is for debugging, but also for standard checks
! that should be maintained to catch any problems immediately and exit
! rather than continuing with bad values.
!
! Author: Steve Penny, June 2017 (visiting scientist, ECMWF)

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC :: debug_post_obslocal, debug_post_obslocal2d, debug_post_letkfcore, debug_ens_diversity, debug_post_anal3d

PRIVATE

contains

subroutine debug_post_obslocal(calling_routine,ij,debug_hdxf_0,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans)
!Inputs:
!  ij           :: index of vectorized model grid
!  debug_hdxf_0 :: flag to decide if this check should be performed
!  nobsl        :: total number of local observations
!  hdxf         :: model equivalent for each ensemble member
!  rdiag        :: diagonal of R matrix
!  rloc         :: weighting for lcoalization of R matrix diagonal elements
!  dep          :: innovation departures
!  parm         :: inflation parameter
!  trans        :: transform matrix for a given variable

  CHARACTER(*) :: calling_routine
  INTEGER,      INTENT(IN) :: ij
  LOGICAL,      INTENT(IN) :: debug_hdxf_0
  INTEGER,      INTENT(IN) :: nobstotal,nobsl
  REAL(r_size), DIMENSION(:,:), INTENT(IN) :: hdxf
  REAL(r_size), DIMENSION(nobsl), INTENT(IN) :: rdiag,rloc,dep
  REAL(r_size), INTENT(IN) :: parm
  REAL(r_size), DIMENSION(:,:), INTENT(IN) :: trans

  WRITE(6,*) "From: ", calling_routine

  ! Make sure there is not a zero value in hdxf, this was a sign of problems in the input fields
  if ( debug_hdxf_0 .and. MINVAL(hdxf) == 0 ) then
    WRITE(6,*) "(3D) ij = ", ij
    WRITE(6,*) "inputs to letkf_core:"
    WRITE(6,*) "nobstotal = ", nobstotal
    WRITE(6,*) "nobsl = ", nobsl
    WRITE(6,*) "hdxf(1:nobsl,1:nbv) = ", hdxf
    WRITE(6,*) "rdiag(1:nobsl) = ", rdiag
    WRITE(6,*) "rloc(1:nobsl) = ", rloc
    WRITE(6,*) "dep(1:nobsl) = ", dep
    WRITE(6,*) "parm = ", parm
    WRITE(6,*) "trans(:,:,n) = ", trans(:,:)
    STOP("STOP::debug_post_obslocal::1")
  endif

  ! Make sure there is no observation error that was input as less than 0
  if ( MINVAL(rdiag(1:nobsl)) .le. 0.0 ) then
    WRITE(6,*) "ERROR: rdiag <=0 (i.e. there is an obserr <= 0)"
    WRITE(6,*) "MINVAL(rdiag) = ",MINVAL(rdiag)
    STOP("STOP::debug_post_obslocal::2")
  endif

end subroutine debug_post_obslocal

subroutine debug_post_obslocal2d(calling_routine,nobsl,rdiag)
!Inputs:
!  nobsl        :: total number of local observations
!  rdiag        :: diagonal of R matrix

  CHARACTER(*) :: calling_routine
  INTEGER,      INTENT(IN) :: nobsl
  REAL(r_size), DIMENSION(nobsl), INTENT(IN) :: rdiag

  WRITE(6,*) "From: ", calling_routine

  ! Make sure there is no observation error that was input as less than 0
  if ( MINVAL(rdiag(1:nobsl)) .le. 0.0 ) then
    WRITE(6,*) "ERROR: rdiag <=0 (i.e. there is an obserr <= 0)"
    WRITE(6,*) "MINVAL(rdiag) = ",MINVAL(rdiag)
    STOP("STOP::debug_post_obslocal2d")
  endif

end subroutine debug_post_obslocal2d

subroutine debug_post_letkfcore(calling_routine,debug_zeroinc,nobsl,n,trans)
!Inputs:
!  debug_zeroinc :: flag to decide if this check should be performed
!  nobsl      :: total number of local observations
!  n          :: variable type index for trans array: trans(:,:,n)
!  trans      :: transform matrix element

  CHARACTER(*) :: calling_routine
  LOGICAL,      INTENT(IN) :: debug_zeroinc
  INTEGER,      INTENT(IN) :: nobsl, n
  REAL(r_size), DIMENSION(:,:), INTENT(IN) :: trans

  WRITE(6,*) "From: ", calling_routine

  if ( debug_zeroinc .and. nobsl > 0 ) then
    WRITE(6,*) "letkf_tools.f90::das_letkf:: post-letkf_core"
    WRITE(6,*) "n = ", n
    WRITE(6,*) "trans(:,:,n) = ", trans
    WRITE(6,*) "letkf_tools.f90::das_letkf:: EXITING on purpose..."
    STOP("STOP::debug_post_letkfcore")
  endif

end subroutine debug_post_letkfcore

subroutine debug_ens_diversity(calling_routine,gues3d,dodebug)
! INPUTS:
!  gues3d  :: the 1st guess array
!  dodebug :: flag to identify print level
 
  CHARACTER(*) :: calling_routine
  REAL(r_size), DIMENSION(:,:,:,:), INTENT(IN) :: gues3d  
  LOGICAL, INTENT(IN) :: dodebug

  WRITE(6,*) "From: ", calling_routine

  if (dodebug) then
    WRITE(6,*) "from letkf_tools.f90:: After computing perturbations..."
    WRITE(6,*) "min val for level gues3d(:,1,:,1) = ", MINVAL(gues3d(:,1,:,1))
    WRITE(6,*) "max val for level gues3d(:,1,:,1) = ", MAXVAL(gues3d(:,1,:,1))
  endif
  if (MAXVAL(gues3d(:,1,:,1)) == MINVAL(gues3d(:,1,:,1))) then
    WRITE(6,*) "from letkf_tools.f90::das_letkf:: It appears that all the ensemble members are identical. EXITING..."
    STOP("STOP::debug_ens_diversity")
  endif

end subroutine debug_ens_diversity

subroutine debug_post_anal3d(calling_routine,anal3d,gues3d,lon1,lat1,k,ij,ilev,m,n,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans)
! INPUTS: 
!  anal3d     :: the analysis value
!  gues3d     :: the 1st guess value
!  lon1       :: the longitude of the model grid point
!  lat1       :: the latitude of the model grid point
!  ij         :: the 2d concatenated model grid array index
!  ilev       :: the vertical level in the model grid
!  m          :: index of ensemble members
!  n          :: 2nd index of ensemble members
!--Inputs to letkf_core:--
!  nobstotal  :: total number of observations
!  nobsl      :: total number of local observations
!  hdxf       :: innovations
!  rdiag      :: diagonal of R matrix
!  rloc       :: localizing scaling factor for diagonal of R matrix
!  dep        :: innovations
!  parm       :: inflation parameter
!  trans      :: transform matrix element

!USE params_letkf, ONLY: nbv

  CHARACTER(*) :: calling_routine
  REAL(r_size), INTENT(IN) :: anal3d,gues3d,lon1,lat1
  INTEGER,      INTENT(IN) :: k,ij,ilev,m,n,nobstotal,nobsl
  REAL(r_size), DIMENSION(:,:), INTENT(IN) :: hdxf
  REAL(r_size), DIMENSION(nobsl), INTENT(IN) :: rdiag,rloc,dep
  REAL(r_size), INTENT(IN) :: parm, trans

  REAL(r_size) :: bad_value_low = -10

  return !STEVE: cice5-specific checks can be put below if desired

  !STEVE: debug - check for bad values
  if ( anal3d < bad_value_low ) then
    WRITE(6,*) "From: ", calling_routine
    WRITE(6,*) "Problem in letkf_das after letkf_core. k = ", k
    WRITE(6,*) "ij, ilev, m, n = ", ij,ilev,m,n
    WRITE(6,*) "anal3d(ij,ilev,m,n) = ", anal3d
    WRITE(6,*) "gues3d(ij,ilev,m,n) = ", gues3d
    WRITE(6,*) "lon = ", lon1
    WRITE(6,*) "lat = ", lat1
    WRITE(6,*) "(3D) ij = ", ij
    WRITE(6,*) "--inputs to letkf_core:"
    WRITE(6,*) "nobstotal = ", nobstotal
    WRITE(6,*) "nobsl = ", nobsl
    WRITE(6,*) "hdxf(1:nobsl,1:nbv) = ", hdxf
    WRITE(6,*) "rdiag(1:nobsl) = ", rdiag
    WRITE(6,*) "rloc(1:nobsl) = ", rloc
    WRITE(6,*) "dep(1:nobsl) = ", dep
    WRITE(6,*) "parm = ", parm
    WRITE(6,*) "trans(k,m,n) = ", trans
    STOP("STOP::debug_post_anal3d")
  endif

end subroutine debug_post_anal3d


end module common_debug_oceanmodel
