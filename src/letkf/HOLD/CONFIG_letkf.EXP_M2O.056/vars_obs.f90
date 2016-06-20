MODULE vars_obs

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshour(:)
REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
INTEGER     ,ALLOCATABLE,SAVE :: obsqc(:)
INTEGER     ,ALLOCATABLE,SAVE :: obsqc0(:,:)
INTEGER     ,ALLOCATABLE,SAVE :: obs_useidx(:) !STEVE: general version of nobs_use()
INTEGER     ,ALLOCATABLE,SAVE :: nobsgrd(:,:) !(nlon,nlat)
!STEVE: for (DRIFTERS)
REAL(r_size),ALLOCATABLE,SAVE :: obsid(:)
REAL(r_size),ALLOCATABLE,SAVE :: obstime(:)


END MODULE vars_obs
