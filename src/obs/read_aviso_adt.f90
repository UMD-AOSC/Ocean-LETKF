MODULE read_aviso_adt
!===============================================================================
! This program reads netcdf to convert delayed-time absolute dynamic topography 
! (ADT) from netcdf to a format readable by letkf. The obsop_adt.f90 program 
! uses this module to read these in directly.
!
! Observation errors should be computed independently.
!
! The global mean volume bias can be subtracted by removing the mean from the
! dataset after reading in via the observation operator.
!===============================================================================

USE common,                     ONLY: r_sngl, r_size, slen
USE params_model,               ONLY: nlev
USE params_obs,                 ONLY: id_eta_obs
USE compute_profile_error,      ONLY: cmpTz

IMPLICIT NONE

PUBLIC :: read_aviso_adt_nc, aviso_adt_data

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n
REAL(r_size) :: se0, seF

TYPE aviso_adt_data
  REAL(r_size) :: x_grd(2)  ! longitude, latitude
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  CHARACTER(9) :: plat      ! Platform
  CHARACTER(3) :: sid       ! Source id
  CHARACTER(1) :: qkey      ! Quality key
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: id     ! id number used in observation files to identify the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: cycle  ! satellite cycle
  INTEGER :: track  ! satellite track
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE aviso_adt_data

TYPE(aviso_adt_data), ALLOCATABLE, DIMENSION(:) :: obs_data

CONTAINS

SUBROUTINE read_aviso_adt_nc(infile,days_since,obs_data,nobs)
!===============================================================================
! Read the aviso adt netcdf file
!===============================================================================

USE netcdf

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: days_since
TYPE(aviso_adt_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs

! Other variables:
REAL :: err
INTEGER :: i,j,k,n
INTEGER :: ncid,istat,varid,dimid
CHARACTER(NF90_MAX_NAME) :: dimname
INTEGER :: grid_z_dim, len1, len2
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: xlon, ylat, hour, day
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: vals, stde
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: times
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: acycle,track
REAL(r_size) :: val, scale_factor
INTEGER :: cnt
LOGICAL :: dodebug=.false.
REAL(r_size) :: missing_value=-99.0
REAL(r_sngl) :: mvc=999
INTEGER :: start_idx, end_idx, subcnt

!-------------------------------------------------------------------------------
! Open netcdf file
!-------------------------------------------------------------------------------
istat = NF90_OPEN(infile,NF90_NOWRITE,ncid)
if (istat /= NF90_NOERR) then
  WRITE(6,'(A)') 'netCDF OPEN ERROR on ', infile
  STOP
endif

!-------------------------------------------------------------------------------
! Read the number of records
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'time',dimid)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID count failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid,dimname,cnt)
if (istat /= NF90_NOERR) then
  print *, "count NF90_INQUIRE_DIMENSION failed"
  STOP
endif
print *, "**************************"
print *, "Record Count is = ", cnt
print *, "**************************"

!STEVE: There SHOULD be a way to read in only the records on the date of interest...

!-------------------------------------------------------------------------------
! Read the times
!-------------------------------------------------------------------------------
ALLOCATE(times(cnt))
istat = NF90_INQ_VARID(ncid,'time',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID time failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,times)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR depth failed"
  STOP
ELSE
  print *, "times(1:10) = ", times(1:10)
!  STOP 0
endif

! Scan the times and pick out the start and end records of the date of interest:
!start_idx=0
!do i=1,cnt
!  if (FLOOR(times(i)) >= days_since) then
!    start_idx = i 
!    EXIT
!  endif
!enddo
!do j=start_idx,cnt-1
!  if (FLOOR(times(j+1)) > days_since) then
!    EXIT
!  endif
!enddo
!if (FLOOR(times(cnt)) == days_since) then
!  end_idx = cnt
!else
!  end_idx = j 
!endif
!subcnt = end_idx - start_idx
start_idx = 1
end_idx = cnt
subcnt = cnt
i = 1
j = cnt
print *, "---------------------------------------------------------------------"
print *, "start_idx  = ", start_idx
print *, "end_idx    = ", end_idx
print *, "days_since (INPUT,but not currently used) = ", days_since
print *, "subcnt     = ", subcnt
print *, "times(i)   = ", times(i)
print *, "times(j)   = ", times(j)
print *, "---------------------------------------------------------------------"

!-------------------------------------------------------------------------------
! Read the longitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(xlon(cnt))
istat = NF90_INQ_VARID(ncid,'longitude',varid)  
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID longitude failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,xlon,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR longitude failed"
  STOP
endif

if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR scale_factor failed"
  STOP
else
  print *, "xlon(1:10) = ", xlon(1:10)
endif

istat = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
print *, "scale_factor = ", scale_factor
xlon = xlon * scale_factor
print *, "xlon(1:10) = ", xlon(1:10)

!-------------------------------------------------------------------------------
! Read the latitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(ylat(cnt))
istat = NF90_INQ_VARID(ncid,'latitude',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID latitude failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,ylat,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR latitude failed"
  STOP
endif

if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR scale_factor failed"
  STOP
else
  print *, "ylat(1:10) = ", ylat(1:10)
endif

istat = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
print *, "scale_factor = ", scale_factor
ylat = ylat * scale_factor
print *, "ylat(1:10) = ", ylat(1:10)

!-------------------------------------------------------------------------------
! Read the time as 'day'
!-------------------------------------------------------------------------------
ALLOCATE(day(cnt))
ALLOCATE(hour(cnt))
istat = NF90_INQ_VARID(ncid,'time',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID time failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,day,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR day failed"
  STOP
endif
hour = (day - FLOOR(day))*24.0d0

!-------------------------------------------------------------------------------
! Read the cycle
!-------------------------------------------------------------------------------
ALLOCATE(acycle(cnt))
istat = NF90_INQ_VARID(ncid,'cycle',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID cycle failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,acycle,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR cycle failed"
  STOP
else
  print *, "cycle(1:10) = ", acycle(1:10)
endif

!-------------------------------------------------------------------------------
! Read the track
!-------------------------------------------------------------------------------
ALLOCATE(track(cnt))
istat = NF90_INQ_VARID(ncid,'track',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID track failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,track,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR track failed"
  STOP
else
  print *, "track(1:10) = ", track(1:10)
endif

!-------------------------------------------------------------------------------
! Read the observed profile data (temperature or salinity)
!-------------------------------------------------------------------------------
ALLOCATE(vals(cnt))
istat = NF90_INQ_VARID(ncid,'ADT',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID ADT failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,vals,(/ start_idx /),(/ subcnt /))
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR ADT failed"
  STOP
else
  print *, "ADT(1:10) = ", vals(1:10)
endif

istat = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
vals = vals * scale_factor
print *, "ADT(1:10) = ", vals(1:10)
istat = nf90_get_att(ncid, varid, '_FillValue', missing_value)
missing_value = missing_value*scale_factor

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
print *, "Finished reading netCDF file, closing..."
istat = NF90_CLOSE(ncid)
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! STEVE: need to compute elsewhere, but estimate a standard error:
!-------------------------------------------------------------------------------
ALLOCATE(stde(cnt))
stde=0.0 !STEVE: placeholder for now. Perhaps better to do this in the calling function (e.g. obsop_tprof.f90)

! CONVERT netcdf data to aviso_adt_data format:
print *, "Finished reading netCDF file, formatting data..."

nobs = end_idx - start_idx + 1

print *, "ALLOCATING obs_data with nobs = ", nobs
ALLOCATE(obs_data(nobs))
n = 0
do i=1,nobs
  if (dodebug) print *, "i = ", i
! if (dodebug) STOP(1)

  val = vals(i)
  err = stde(i)
  if (val < missing_value) then
    n = n+1
    if (dodebug) print *, "n,lon,lat,hour,val,err = ", n,xlon(i),ylat(i),hour(i),val,err
    obs_data(n)%typ = id_eta_obs
    obs_data(n)%x_grd(1) = xlon(i)
    obs_data(n)%x_grd(2) = ylat(i)
    obs_data(n)%hour = hour(i)
    obs_data(n)%value = val
    obs_data(n)%oerr = err
    obs_data(n)%cycle = acycle(i)
    obs_data(n)%track = track(i)
  endif
enddo
nobs = n
if (dodebug) print *, "nobs = ", nobs

END SUBROUTINE read_aviso_adt_nc

END MODULE read_aviso_adt
