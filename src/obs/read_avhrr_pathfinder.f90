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
  USE m_ncio,     ONLY: nc_get_fid, nc_close_fid, nc_rddim, nc_rdvar2d, nc_rdvar1d, &
                        nc_rdatt, nc_fndvar
IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
INTEGER, INTENT(IN) :: min_quality_level
TYPE(avhrr_pathfinder_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs
REAL(r_size) :: scale_error = 2.0 !STEVE: scale the sst_dtime as an estimate of obs error. Make this an (optional) input.

! Other variables:
REAL(r_sngl) :: err
INTEGER :: i,j,k,n
INTEGER :: cnt
INTEGER :: ncid,istat,varid
INTEGER :: time(1)
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: alon, alat
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: sea_surface_temperature
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: sea_ice_fraction
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: dt_analysis
REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: stde
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: quality_level
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sst_dtime
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: valid
REAL(r_size) :: val
REAL(r_size) :: rFillValue
INTEGER :: iFillValue
INTEGER :: nlons,nlats
REAL(r_size) :: hour
REAL(r_size) :: add_offset, scale_factor
INTEGER :: obsdate(8), refdate(8)
REAL :: tdelta(5)
LOGICAL :: dodebug=.true.

!-------------------------------------------------------------------------------
! Open netcdf file
!-------------------------------------------------------------------------------
CALL nc_get_fid(trim(infile), ncid)

!-------------------------------------------------------------------------------
! Read the dimension of the lon/lat grid
!-------------------------------------------------------------------------------
CALL nc_rddim(ncid, "lon", nlons)
print *, "**************************"
print *, "lon dimension is = ", nlons
print *, "**************************"

!-------------------------------------------------------------------------------
CALL nc_rddim(ncid, "lat", nlats)
print *, "**************************"
print *, "lat dimension is = ", nlats
print *, "**************************"

!-------------------------------------------------------------------------------
! Read the longitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(alon(nlons))
CALL nc_rdvar1d(ncid, "lon", alon)
print*, "[msg] read_avhrr_pathfinder_nc::lon: min, max=", &
        minval(alon), maxval(alon)

!-------------------------------------------------------------------------------
! Read the latitude coordinates
!-------------------------------------------------------------------------------
ALLOCATE(alat(nlats))
CALL nc_rdvar1d(ncid, "lat", alat)
print*, "[msg] read_avhrr_pathfinder_nc::lat: min, max=", &
        minval(alat), maxval(alat)

!-------------------------------------------------------------------------------
! allocate vars to filter out grids with missing values for any vars
!-------------------------------------------------------------------------------
ALLOCATE(valid(nlons,nlats))
valid = .true.

!-------------------------------------------------------------------------------
! Read the base time and time offset
!-------------------------------------------------------------------------------
CALL nc_rdvar1d(ncid, "time", time)
refdate = [1981, 1, 1, 0, 0, 0, 0, 0] ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
tdelta  = [0.0, 0.0, 0.0, REAL(time(1)), 0.0] ! DAY/HR/MIN/SEC/MSEC
CALL w3movdat(tdelta, refdate, obsdate)
print*, "[msg] read_avhrr_pathfinder_nc::obs time=", obsdate(1:3), obsdate(5:7)

ALLOCATE(sst_dtime(nlons,nlats))
CALL nc_rdvar2d(ncid, "sst_dtime", sst_dtime)
CALL nc_rdatt(  ncid, "sst_dtime", "_FillValue", iFillValue)

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
CALL nc_rdvar2d(ncid, "quality_level", quality_level)
CALL nc_rdatt(  ncid, "quality_level", "_FillValue", iFillValue)
where (quality_level == iFillValue)
    valid = .false.
end where
print*, "[msg] read_avhrr_pathfinder_nc::quality_level: min, max=", &
        minval(quality_level,mask=valid), maxval(quality_level,mask=valid)

!-------------------------------------------------------------------------------
! Read the change from yesterday, as a proxy for error estimate
!-------------------------------------------------------------------------------
ALLOCATE(dt_analysis(nlons,nlats))
CALL nc_rdvar2d(ncid, "dt_analysis", dt_analysis)
CALL nc_rdatt(  ncid, "dt_analysis", "scale_factor", scale_factor)
CALL nc_rdatt(  ncid, "dt_analysis", "_FillValue", iFillValue)
dt_analysis = dt_analysis * scale_factor
where (dt_analysis == iFillValue)
    valid = .false.
end where
print*, "[msg] read_avhrr_pathfinder::dt_analysis: min, max=", &
        minval(dt_analysis,mask=valid), maxval(dt_analysis,mask=valid)

!-------------------------------------------------------------------------------
! Read the observed SST data
!-------------------------------------------------------------------------------
if (typ .eq. id_sst_obs) then
  ALLOCATE(sea_surface_temperature(nlons,nlats))
  CALL nc_rdvar2d(ncid, "sea_surface_temperature", sea_surface_temperature)
  CALL nc_rdatt(  ncid, "sea_surface_temperature", "scale_factor", scale_factor)
  CALL nc_rdatt(  ncid, "sea_surface_temperature", "add_offset",   add_offset)
  CALL nc_rdatt(  ncid, "sea_surface_temperature", "_FillValue",   iFillValue)
  sea_surface_temperature = sea_surface_temperature * scale_factor + add_offset
  sea_surface_temperature = cvt_temp_K2C(sea_surface_temperature) ! From K to degC, since model state uses degC
  where (NINT(sea_surface_temperature)==iFillValue)
       valid = .false.
  end where
  print*, "[msg] read_avhrr_pathfinder::sea_surface_temperature: min, max=", &
           minval(sea_surface_temperature,mask=valid), &
           maxval(sea_surface_temperature,mask=valid)
  
elseif (typ .eq. id_sic_obs) then ! Sea ice concentration (fraction)
  ALLOCATE(sea_ice_fraction(nlons,nlats))
  CALL nc_rdvar2d(ncid, "sea_ice_fraction", sea_ice_fraction)
  CALL nc_rdatt(  ncid, "sea_ice_fraction", "scale_factor", scale_factor)
  CALL nc_rdatt(  ncid, "sea_ice_fraction", "add_offset",   add_offset)
  CALL nc_rdatt(  ncid, "sea_ice_fraction", "_FillValue",   iFillValue)
  sea_ice_fraction = sea_ice_fraction * scale_factor + add_offset
  where (NINT(sea_ice_fraction)==iFillValue)
        valid = .false.
  end where
  print*, "[msg] read_avhrr_pathfinder::sea_ice_fraction: min, max=", &
           minval(sea_ice_fraction,mask=valid), &
           maxval(sea_ice_fraction,mask=valid)
endif

!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
print *, "Finished reading netCDF file, closing..."
CALL nc_close_fid(ncid)

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
      err = 1.0 + scale_error*ABS(dt_analysis(i,j)) !stde(i,j)
    elseif (typ .eq. id_sic_obs) then
      val = sea_ice_fraction(i,j)
      err = stde(i,j)
    endif
    if (quality_level(i,j) >= min_quality_level) then
      n = n+1
      hour = ( sst_dtime(i,j) ) / 3600.0d0
      !if (dodebug) print *, "n,lon,lat,hour,val,err,qkey = ", n,alon(i),alat(j),hour,val,err,quality_level(i,j)
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
CALL inspect_obs_data(obs_data(1:nobs))
if (dodebug) print *, "nobs = ", nobs

END SUBROUTINE read_avhrr_pathfinder_nc

SUBROUTINE inspect_obs_data(obs_data)
  IMPLICIT NONE
  TYPE(avhrr_pathfinder_data),INTENT(IN) :: obs_data(:)

  print*, "[msg] read_avhrr_pathfinder_nc::info"
  print*, "                nobs=", size(obs_data)
  print*, "  x_grd(1): min, max=", minval(obs_data(:)%x_grd(1)), maxval(obs_data(:)%x_grd(1))
  print*, "  x_grd(2): min, max=", minval(obs_data(:)%x_grd(2)), maxval(obs_data(:)%x_grd(2))
  print*, "      hour: min, max=", minval(obs_data(:)%hour), maxval(obs_data(:)%hour)
  print*, "     value: min, max=", minval(obs_data(:)%value), maxval(obs_data(:)%value)
  print*, "      oerr: min, max=", minval(obs_data(:)%oerr), maxval(obs_data(:)%oerr)
  print*, "      qkey: min, max=", minval(obs_data(:)%qkey), maxval(obs_data(:)%qkey)

END SUBROUTINE


ELEMENTAL FUNCTION cvt_temp_K2C(tf) RESULT (tc)
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: tf   ! temperature in Kelvin
  REAL(r_size)            :: tc   ! temperature in degree Celsius

  REAL(r_size),PARAMETER :: ADD_OFFSET_K2C = -273.15_r_size

  tc = tf + ADD_OFFSET_K2C

END FUNCTION

END MODULE read_avhrr_pathfinder
