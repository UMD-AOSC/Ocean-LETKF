MODULE read_avhrr_pathfinder
!===============================================================================
! This program reads netcdf to convert avhrr pathfinder data from
! http://data.nodc.noaa.gov/pathfinder/Version5.2/
! to a format readable by letkf. The obsop_sst.f90 program uses this module to read these
! in directly.
!
! Observation errors should be computed externally or read from an observation error file
! (e.g. when doing adaptive obs error)
!
! I assume nighttime data is preferable for climate applications. Bias correction may
! be needed on the satellite fields to match avhrr with in situ floats.
!
!===============================================================================

USE common,       ONLY: r_sngl, r_size, slen
USE params_obs,   ONLY: id_sst_obs, id_sic_obs

IMPLICIT NONE

PUBLIC :: read_avhrr_pathfinder_nc, avhrr_pathfinder_data

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n

TYPE avhrr_pathfinder_data
  REAL(r_size) :: x_grd(2)  ! longitude, latitude
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  INTEGER :: qkey      ! Quality key
  INTEGER :: typ    ! type of observation (elem)
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE avhrr_pathfinder_data

TYPE(avhrr_pathfinder_data), ALLOCATABLE, DIMENSION(:) :: obs_data

!! Write letkf file
!do i=1,nobs
!!STEVE: the following are required for miyoshi's letkf observation input:
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

SUBROUTINE read_avhrr_pathfinder_nc(infile,typ,min_quality_level,obs_data,nobs)
!===============================================================================
! Read the argo profile data
!===============================================================================

USE netcdf

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
INTEGER, INTENT(IN) :: min_quality_level
TYPE(avhrr_pathfinder_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs

! Other variables:
REAL(r_sngl) :: err
INTEGER :: i,j,k,n
INTEGER :: cnt
INTEGER :: ncid,istat,varid,dimid1,dimid2
CHARACTER(NF90_MAX_NAME) :: dimname
INTEGER :: time
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: alon, alat
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: sea_surface_temperature
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: sea_ice_fraction
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: stde
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: quality_level
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sst_dtime
REAL(r_size) :: val
INTEGER :: nlons,nlats
REAL(r_size) :: hour
LOGICAL :: dodebug=.true.

!-------------------------------------------------------------------------------
! Open netcdf file
!-------------------------------------------------------------------------------
istat = NF90_OPEN(infile,NF90_NOWRITE,ncid)
if (istat /= NF90_NOERR) then
  WRITE(6,'(A)') 'netCDF OPEN ERROR on ', infile
  STOP
endif

!-------------------------------------------------------------------------------
! Read the dimension of the lon/lat grid
!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'lon',dimid1)
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_DIMID lon failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid1,dimname,nlons)
if (istat /= NF90_NOERR) then
  print *, "lon NF90_INQUIRE_DIMENSION failed"
  STOP
endif
print *, "**************************"
print *, "lon dimension is = ", nlons
print *, "**************************"

!-------------------------------------------------------------------------------
istat = NF90_INQ_DIMID(ncid,'lat',dimid2)
if (istat /= NF90_NOERR) then 
  print *, "NF90_INQ_DIMID lat failed"
  STOP
endif
istat = NF90_INQUIRE_DIMENSION(ncid,dimid2,dimname,nlats)
if (istat /= NF90_NOERR) then
  print *, "lat NF90_INQUIRE_DIMENSION failed"
  STOP
endif
print *, "**************************"
print *, "lat dimension is = ", nlats
print *, "**************************"

!-------------------------------------------------------------------------------
! Read the longitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(alon(nlons))
istat = NF90_INQ_VARID(ncid,'lon',varid)  
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID lon failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,alon)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR lon failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the latitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(alat(nlats))
istat = NF90_INQ_VARID(ncid,'lat',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID lat failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,alat)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR lat failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the base time and time offset
!-------------------------------------------------------------------------------
istat = NF90_INQ_VARID(ncid,'time',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID time failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,time)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR time failed"
  STOP
endif

ALLOCATE(sst_dtime(nlons,nlats))
istat = NF90_INQ_VARID(ncid,'sst_dtime',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID sst_dtime failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,sst_dtime)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR sst_dtime failed"
  STOP
endif

!-------------------------------------------------------------------------------
! Read the Quality Key
! "Note, the native Pathfinder processing system returns quality levels ranging 
! from 0 to 7 (7 is best quality; -1 represents missing data) and has been 
! converted to the extent possible into the six levels required by the GDS2 
! (ranging from 0 to 5, where 5 is best). Below is the conversion table: 
! \n GDS2 required quality_level 5  =  native Pathfinder quality level 7 == best_quality 
! \n GDS2 required quality_level 4  =  native Pathfinder quality level 4-6 == acceptable_quality 
! \n GDS2 required quality_level 3  =  native Pathfinder quality level 2-3 == low_quality 
! \n GDS2 required quality_level 2  =  native Pathfinder quality level 1 == worst_quality 
! \n GDS2 required quality_level 1  =  native Pathfinder quality level 0 = bad_data 
! \n GDS2 required quality_level 0  =  native Pathfinder quality level -1 = missing_data 
! \n The original Pathfinder quality level is recorded in the optional variable pathfinder_quality_level."
!-------------------------------------------------------------------------------
ALLOCATE(quality_level(nlons,nlats))
istat = NF90_INQ_VARID(ncid,'quality_level',varid)   
if (istat /= NF90_NOERR) then
  print *, "NF90_INQ_VARID quality_level failed"
  STOP
endif
istat = NF90_GET_VAR(ncid,varid,quality_level)
if (istat /= NF90_NOERR) then
  print *, "NF90_GET_VAR quality_level failed"
  STOP
else
  print *, "quality_level(1,1) = ", quality_level(1,1)
endif

!-------------------------------------------------------------------------------
! Read the observed SST data
!-------------------------------------------------------------------------------
if (typ .eq. id_sst_obs) then
  ALLOCATE(sea_surface_temperature(nlons,nlats))
  istat = NF90_INQ_VARID(ncid,'sea_surface_temperature',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID sea_surface_temperature failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,sea_surface_temperature)
  if (dodebug) print *, "sea_surface_temperature(1,1) = ", sea_surface_temperature(1,1)

  !STEVE: read these attributes from file:
  ! sea_surface_temperature:add_offset = 273.15 ;
  ! sea_surface_temperature:scale_factor = 0.01 ;
  ! sea_surface_temperature:_FillValue = -32768s ;
  sea_surface_temperature = sea_surface_temperature + 273.15
  sea_surface_temperature = sea_surface_temperature * 0.01
  
elseif (typ .eq. id_sic_obs) then ! Sea ice concentration (fraction)
  ALLOCATE(sea_ice_fraction(nlons,nlats))
  istat = NF90_INQ_VARID(ncid,'sea_ice_fraction',varid)   
  if (istat /= NF90_NOERR) then
    print *, "NF90_INQ_VARID sea_ice_fraction failed"
    STOP
  endif
  istat = NF90_GET_VAR(ncid,varid,sea_ice_fraction)
  if (dodebug) print *, "sea_ice_fraction(1,1) = ", sea_ice_fraction(1,1)
endif
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
print *, "Finished reading netCDF file, closing..."
istat = NF90_CLOSE(ncid)
if (istat /= NF90_NOERR) STOP

!-------------------------------------------------------------------------------
! STEVE: need to compute the standard error by some logical means in the calling function
!-------------------------------------------------------------------------------
ALLOCATE(stde(nlons,nlats))
stde=0.0 !STEVE: placeholder for now. Perhaps better to do this in the calling function (e.g. obsop_sst.f90)

! CONVERT netcdf data to argo_data format:
print *, "Finished reading netCDF file, formatting data..."

cnt = nlons*nlats
print *, "ALLOCATING obs_data with nobs = ", cnt
ALLOCATE(obs_data(cnt))
n = 0
do j=1,nlats
  do i=1,nlons
    if (typ .eq. id_sst_obs) then
      val = sea_surface_temperature(i,j)
      err = stde(i,j)
    elseif (typ .eq. id_sic_obs) then
      val = sea_ice_fraction(i,j)
      err = stde(i,j)
    endif
    if (quality_level(i,j) >= min_quality_level) then
      n = n+1
      hour = ( sst_dtime(i,j) ) / 3600.0d0
      print *, "n,lon,lat,hour,val,err,qkey = ", n,alon(i),alat(j),hour,val,err,quality_level(i,j)
      obs_data(n)%typ = typ
      obs_data(n)%x_grd(1) = alon(i)
      obs_data(n)%x_grd(2) = alat(j)
      obs_data(n)%hour = hour
      obs_data(n)%value = val
      obs_data(n)%oerr = err
      obs_data(n)%qkey = quality_level(i,j)
    endif
  enddo
enddo
nobs = n
if (dodebug) print *, "nobs = ", nobs

END SUBROUTINE read_avhrr_pathfinder_nc

END MODULE read_avhrr_pathfinder
