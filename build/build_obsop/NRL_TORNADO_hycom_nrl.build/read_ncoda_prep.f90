MODULE read_ncoda_prep
!===============================================================================
! This program reads the NCODA prep files:
! coda.MVO_prp.*
! coda.SSH_prp.*
! coda.MOV_ENS_obs.*
! coda.SSH_ENS_obs.*
!
! These containt the observations (prp) and ensemble member innovations (ENS_obs)
!
! This routine reads both files and recovers the observation and model
! equivalent for a specified member.
!
! The obsop.f90 program uses this module to read these in directly.
!
! Observation errors are provided in the 'prp' files.
!
! Author:
! 12/22/17:  Stephen G. Penny, University of Maryland, College Park
!                              Visiting Scientist, Naval Research Laboratory - Stennis Space Center
!
!===============================================================================

USE common,                     ONLY: r_sngl, r_size, slen
USE params_obs,                 ONLY: id_t_obs, id_s_obs, id_u_obs, id_v_obs, id_eta_obs
USE compute_profile_error,      ONLY: cmpTz

IMPLICIT NONE

PUBLIC :: read_ncoda_prp, ncoda_prep_data

PRIVATE

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n
REAL(r_size) :: se0, seF

TYPE ncoda_prep_data
  REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: hxb       ! Model equivalent ! (NEMO) added
  REAL(r_size) :: lev       ! grid z level
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  CHARACTER(9) :: plat      ! Platform
  CHARACTER(3) :: ptyp      ! Profile type
  CHARACTER(3) :: sid       ! Source id
  CHARACTER(1) :: qkey      ! Quality key
  INTEGER :: qc             ! Quality control flag  ! (NEMO) added
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
  INTEGER :: id     ! id number used in observation files to identify the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE ncoda_prep_data

LOGICAL :: DO_READ_OBE      ! Read the observation error estimate from the file
LOGICAL :: DO_COMPUTE_OBE   ! Compute the observation error estimate using the NCEP approach w/ vertical gradient

INTEGER :: undef = 99999

!! Write letkf file
!do i=1,nobs
!!STEVE: the following are required for miyoshi's letkf observation input format:
!!1 = obelm
!!2 = lon
!!3 = lat
!!4 = lev
!!5 = value
!!6 = oberr
! wk(1) = obs_data(i)%typ
! wk(2) = obs_data(i)%x_grd(1)
! wk(3) = obs_data(i)%x_grd(2)
! wk(4) = obs_data(i)%x_grd(3)
! wk(5) = obs_data(i)%value
! wk(6) = obs_data(i)%oerr
! WRITE(fid) wk
!enddo

CONTAINS

SUBROUTINE read_ncoda_prp(infile,infile2,obs_data,nobs,mem,obid_in)
!===============================================================================
! Read the ncoda prep observation data and model equivalent
!===============================================================================

USE params_letkf, ONLY: nbv

CHARACTER(*), INTENT(IN) :: infile, infile2
TYPE(ncoda_prep_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs
INTEGER, INTENT(IN) :: mem ! Specify ensemble member to process
INTEGER, INTENT(IN), OPTIONAL :: obid_in

! Other variables:
INTEGER :: i,j,k,n
REAL(r_sngl), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour
CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: ptyp, sid
CHARACTER(9), ALLOCATABLE, DIMENSION(:) :: plat
CHARACTER(1), ALLOCATABLE, DIMENSION(:) :: qkey
INTEGER, ALLOCATABLE, DIMENSION(:) :: qc
REAL(r_sngl), ALLOCATABLE, DIMENSION(:) :: ovals, mvals, stde
REAL(r_sngl), ALLOCATABLE, DIMENSION(:) :: depth ! profile depths  (dimensioned: nobs, depth)
REAL(r_sngl) :: val, hxb, err
INTEGER :: cnt, nlv
LOGICAL :: dodebug=.true.
REAL(r_sngl) :: missing_value=-999  ! Missing value flag in the data
REAL(r_sngl) :: max_value=999       ! Maximum value condition
REAL(r_sngl) :: max_depth = 9999    ! Maximum depth condition

INTEGER :: fid = 21
LOGICAL :: ex
INTEGER :: obid

REAL(r_sngl) :: obs_anm
REAL(r_sngl), ALLOCATABLE :: obs_xi (:)
REAL(r_sngl), ALLOCATABLE :: obs_yj (:)
REAL(r_sngl), ALLOCATABLE :: obs_zk (:)
INTEGER,      ALLOCATABLE :: otype (:)
INTEGER,      ALLOCATABLE :: otime (:)
INTEGER,      ALLOCATABLE :: odate (:)
REAL(r_sngl), ALLOCATABLE :: ens_anm(:)
REAL(r_sngl), ALLOCATABLE :: dum_anm(:)
INTEGER,      ALLOCATABLE :: etyp (:)
INTEGER,      ALLOCATABLE :: obs_var (:)

!var is the variable 1-6 according to:
!     data      var_lbl / 'seatmp', 'salint', 'geoptl',
!    *                    'uucurr', 'vvcurr', 'lyrprs' /

! Check optional argument
if (present(obid_in)) then
  obid = obid_in
else
  obid = id_t_obs ! (Default, reads both id_t_obs and id_s_obs)
endif

!-------------------------------------------------------------------------------
! Open obs prep file
!-------------------------------------------------------------------------------
inquire (file=trim(infile), exist=ex)
if (ex) then
  open (unit=fid,file=trim(infile),status='old',form='unformatted')
! rewind(fid)
else
  print *, "infile = ", trim(infile) 
  STOP('ERROR: prp file is missing. EXITING...')
endif

!-------------------------------------------------------------------------------
! Read the number of obs
!-------------------------------------------------------------------------------
read (fid) cnt
nobs=cnt !STEVE: ?
print *, 'read innovation file 1: ', trim(infile)
print *, "nobs, cnt, cnt-nobs = ", nobs, cnt, cnt-nobs

!-------------------------------------------------------------------------------
! Allocate arrays for obs data
!-------------------------------------------------------------------------------
ALLOCATE (ovals(nobs))
ALLOCATE (stde(nobs))
ALLOCATE (xlon(nobs))
ALLOCATE (ylat(nobs))
ALLOCATE (depth(nobs))
ALLOCATE (hour(nobs))
ALLOCATE (etyp(nobs))
ALLOCATE (qc(nobs))

ALLOCATE (obs_xi(nobs))
ALLOCATE (obs_yj(nobs))
ALLOCATE (obs_zk(nobs))
ALLOCATE (obs_var(nobs))
ALLOCATE (otype(nobs))
ALLOCATE (otime(nobs))
ALLOCATE (odate(nobs))

!-------------------------------------------------------------------------------
! Read the obs data
!-------------------------------------------------------------------------------
WRITE(6,*) "Reading obs data, cnt = ", cnt
if (cnt .gt. 0) then
  do i = 1, nobs
!   print *, "start reading obs i = ", i
    read (fid) odate(i), otime(i), obs_var(i), otype(i), &
     &     ovals(i), obs_anm, stde(i), depth(i), ylat(i), &
     &     xlon(i), obs_xi(i),  obs_yj(i),  obs_zk(i)

    hour(i) = otime(i)
    
    ! Specify the observation type:
    qc(i) = 1
    if (obid == id_t_obs) then
      select case (obs_var(i))
        case(1)
          etyp(i) = id_t_obs
        case(2)
          etyp(i) = id_s_obs
        case default
!         print *, "Observation variable type not recognized: ", obs_var(i)
!         print *, "REMOVING..."
          qc(i) = 0
      end select
    else
      ! Use this for ssh or sea ice ncoda-prep data:
      etyp(i) = obid
    endif
  enddo

  print *, "final nobs, xlon,ylat,depth = ", nobs, xlon(nobs),ylat(nobs),depth(nobs)
  print *, "obs_var,ovals,stde,obs_anm = ", obs_var(nobs),ovals(nobs),stde(nobs),obs_anm

else
  WRITE (6,*) "cnt==0, EXITING..."
  STOP(13)
endif 

!-------------------------------------------------------------------------------
! Close the obs file
!-------------------------------------------------------------------------------
close(fid)

!-------------------------------------------------------------------------------
! Open ens obs prep file
!-------------------------------------------------------------------------------
inquire (file=trim(infile2), exist=ex)
if (ex) then
  open (unit=fid,file=trim(infile2),status='old',form='unformatted')
! rewind(fid)
else
  print *, "infile2 = ", trim(infile2) 
  STOP('ERROR: ENS_obs file is missing. EXITING...')
endif
print *, 'read innovation file 2: ', trim(infile2)

!-------------------------------------------------------------------------------
! Allocate arrays for innovation data
!-------------------------------------------------------------------------------
ALLOCATE (ens_anm(nobs))
ALLOCATE (dum_anm(4)) !nbv))

!-------------------------------------------------------------------------------
! Read the ensemble innovation file
!-------------------------------------------------------------------------------
if (cnt .gt. 0) then
  do i = 1, nobs
!   print *, "start reading ens_anm i = ", i
    read (fid) dum_anm
    ens_anm(i) = dum_anm(mem)
!   print *, "after reading ens_anm i = ", i
  enddo
  print *, "ens_anm(1) = ",  ens_anm(1)
  print *, "ens_anm(nobs-1) = ", ens_anm(nobs-1)
endif

!-------------------------------------------------------------------------------
! Close the ensemble innovation file
!-------------------------------------------------------------------------------
close(fid)

! CONVERT data to ncoda_prep_data format:
print *, "Finished reading NCODA prep files, formatting data..."

! Loop through all profiles and create a new observation for each unique data point
!nobs = cnt * nlv
print *, "ALLOCATING obs_data with nobs = ", nobs
ALLOCATE(obs_data(nobs))
n = 0
do i=1,cnt
  if (qc(i)==0) CYCLE  ! Skip any obs deemed unusable
  if (dodebug) print *, "i = ", i

  ! Assign observation value:
  val = ovals(i)
  err = stde(i)
! if (etyp(i)==id_eta_obs) then
!   hxb = ens_anm(i)
! else
  hxb = val - ens_anm(i)
! endif

  if (abs(val) < abs(max_value) .and. abs(val) < abs(missing_value)-1 .and. depth(i) < max_depth) then
    n = n+1
    if (dodebug) print *, "n,lon,lat,depth,hour,val,err,hxb,qc = ", n,xlon(i),ylat(i),depth(i),hour(i),val,err,hxb,qc(i)

    obs_data(n)%typ = etyp(i)
    obs_data(n)%x_grd(1) = xlon(i)
    obs_data(n)%x_grd(2) = ylat(i)
    obs_data(n)%x_grd(3) = depth(i)
    obs_data(n)%hour = hour(i)
    obs_data(n)%value = val
    obs_data(n)%hxb   = hxb
    obs_data(n)%oerr = err
    obs_data(n)%rid = obs_var(i) !STEVE: ?
    obs_data(n)%lid = obs_zk(i)  ! level id
  endif
enddo

nobs = n
if (dodebug) print *, "nobs = ", nobs

! Explicitly deallocate temporary arrays
if (ALLOCATED(depth)) DEALLOCATE(depth)
if (ALLOCATED(xlon))  DEALLOCATE(xlon)
if (ALLOCATED(ylat))  DEALLOCATE(ylat)
if (ALLOCATED(hour))  DEALLOCATE(hour)
if (ALLOCATED(qc))    DEALLOCATE(qc)
if (ALLOCATED(ovals))  DEALLOCATE(ovals)
if (ALLOCATED(ens_anm))  DEALLOCATE(ens_anm)
if (ALLOCATED(stde))  DEALLOCATE(stde)

if (dodebug) print *, "Temporary arrays deallocated."
if (dodebug) print *, "Returning from read_ncoda_prp..."

END SUBROUTINE read_ncoda_prp

END MODULE read_ncoda_prep
