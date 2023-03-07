MODULE read_viirs
  USE common,       ONLY: r_sngl, r_dble, r_size, slen
  USE params_obs,   ONLY: id_sst_obs
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sst_viirs_data
  PUBLIC :: read_viirs_nc
  PRIVATE :: inspect_obs_data

  TYPE sst_viirs_data
    REAL(r_size) :: x_grd(2)  ! longitude, latitude
    REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
    REAL(r_size) :: oerr      ! observation standard error
    REAL(r_size) :: hour      ! Hour of observation
    INTEGER :: qkey           ! Quality key
    INTEGER :: typ            ! type of observation (elem)
    LOGICAL :: kept           ! tells letkf whether this obs is kept for assimilation
  END TYPE sst_viirs_data

CONTAINS

SUBROUTINE read_viirs_nc(obsinfile, min_quality_level, obs_data, nobs)
  USE m_ncio,     ONLY: nc_get_fid, nc_close_fid, nc_rddim, nc_rdvar2d, nc_rdvar1d, &
                        nc_rdatt, nc_fndvar
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: obsinfile               
  INTEGER,     INTENT(IN) :: min_quality_level
  TYPE(sst_viirs_data),ALLOCATABLE,INTENT(INOUT) :: obs_data(:)
  INTEGER,     INTENT(OUT) :: nobs

  INTEGER :: fid, ni, nj, i, j, n, qlev

  REAL(r_size),ALLOCATABLE :: alon2d(:,:), alat2d(:,:)
  REAL(r_size),ALLOCATABLE :: sea_surface_temperature(:,:), stde(:,:), sst_dtime(:,:)
  INTEGER(4),  ALLOCATABLE :: quality_level(:,:)
  INTEGER(1),  ALLOCATABLE :: i1buf2d(:,:)
  INTEGER(2),  ALLOCATABLE :: i2buf2d(:,:)
  LOGICAL,     ALLOCATABLE :: valid(:,:) ! .true. if no missing info at this grid by checking different FillValues in nc file
  INTEGER :: FillValue ! missing value indicated by vars in nc file
  REAL(r_size) :: offset, scale_factor
  LOGICAL :: dodebug = .true.
  INTEGER :: obsdate(8), sdate(8)
  REAL :: tdelta(5)
  INTEGER :: time(1)

!-------------------------------------------------------------------------------
! Read dimension & geo info:
!-------------------------------------------------------------------------------
  CALL nc_get_fid(trim(obsinfile),fid)
  CALL nc_rddim(fid,"ni", ni)
  print *, "**************************"
  print *, "ni dimension is = ", ni
  print *, "**************************"

  CALL nc_rddim(fid,"nj", nj)
  print *, "**************************"
  print *, "nj dimension is = ", nj
  print *, "**************************"

  ALLOCATE(alon2d(ni,nj),alat2d(ni,nj))
  ALLOCATE(valid(ni,nj))
  valid = .true.

  CALL nc_rdvar2d(fid,"lon",alon2d)
  WRITE(6,*) "[msg] read_viirs_nc::alon2d: min, max=", &
             minval(alon2d),maxval(alon2d)
  CALL nc_rdvar2d(fid,"lat",alat2d)
  WRITE(6,*) "[msg] read_viirs_nc::alat2d: min, max=", &
             minval(alat2d),maxval(alat2d)

!-------------------------------------------------------------------------------
! Read the observed SST data & other info from obs file (obsinfile)
! [NOTE]:
! The SST data from file uses K as its units, while MOM6 uses degC
! We convert SST from K to degC here.
!-------------------------------------------------------------------------------

  ALLOCATE(i1buf2d(ni,nj))
  ALLOCATE(i2buf2d(ni,nj))


  CALL nc_rdvar2d(fid, "sea_surface_temperature", i2buf2d)
  CALL nc_rdatt(fid, "sea_surface_temperature","scale_factor", scale_factor)
  CALL nc_rdatt(fid, "sea_surface_temperature","add_offset", offset)
  CALL nc_rdatt(fid, "sea_surface_temperature","_FillValue", FillValue)
  if (dodebug) WRITE(6,*) "_FillValue, scale_factor, add_offset=", &
               FillValue, scale_factor, offset
  where (i2buf2d == FillValue)
      valid = .false.
  end where
  ALLOCATE(sea_surface_temperature(ni,nj))
  sea_surface_temperature = REAL(i2buf2d,r_size)*scale_factor + offset
  WRITE(6,*) "[msg] read_viirs_nc::sst: min, max=", &
             minval(sea_surface_temperature), maxval(sea_surface_temperature), &
             minval(sea_surface_temperature, mask=valid),&
             maxval(sea_surface_temperature, mask=valid)

  where (valid)
      sea_surface_temperature = cvt_temp_K2C(sea_surface_temperature)
  end where
  WRITE(6,*) "[msg] read_viirs_nc::sst: Converting from K to degC"
  WRITE(6,*) "[msg] read_viirs_nc::sst: min, max=", &
             minval(sea_surface_temperature), maxval(sea_surface_temperature), &
             minval(sea_surface_temperature, mask=valid),&
             maxval(sea_surface_temperature, mask=valid)

! SST error
  CALL nc_rdvar2d(fid, "sses_standard_deviation", i1buf2d)
  CALL nc_rdatt(fid, "sses_standard_deviation", "scale_factor", scale_factor)
  CALL nc_rdatt(fid, "sses_standard_deviation", "add_offset", offset)
  CALL nc_rdatt(fid, "sses_standard_deviation", "_FillValue", FillValue)
  if (dodebug) WRITE(6,*) "_FillValue, scale_factor, add_offset=", &
               FillValue, scale_factor, offset
  where (i1buf2d == FillValue)
     valid = .false.
  end where
  ALLOCATE(stde(ni,nj))
  stde = REAL(i1buf2d,r_size)*scale_factor + offset
  WRITE(6,*) "[msg] read_viirs_nc::stde: min, max=", &
             minval(stde), maxval(stde), &
             minval(stde, mask=valid),maxval(stde, mask=valid)

!-------------------------------------------------------------------------------
! Read the base time and time offset
!-------------------------------------------------------------------------------
  CALL nc_rdvar1d(fid, "time", time)
  WRITE(6,*) "[msg] read_viirs_nc::time=",time(1)
  sdate = [1981, 1, 1, 0, 0, 0, 0, 0] ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
  tdelta = [0.0,0.0,0.0,REAL(time(1)),0.0] ! DAY/HR/MIN/SEC/MSEC
  CALL w3movdat(tdelta, sdate, obsdate)
  WRITE(6,*) "[msg] read_viirs_nc::obstime=", obsdate(1:3),obsdate(5:7)

! sst_dtime
  CALL nc_rdvar2d(fid, "sst_dtime", i2buf2d)
  CALL nc_rdatt(fid, "sst_dtime", "scale_factor", scale_factor)
  CALL nc_rdatt(fid, "sst_dtime", "add_offset", offset)
  CALL nc_rdatt(fid, "sst_dtime", "_FillValue", FillValue)
  if (dodebug) WRITE(6,*) "_FillValue, scale_factor, add_offset=", &
               FillValue, scale_factor, offset
  where (i2buf2d == FillValue)
     valid = .false.
  end where
  ALLOCATE(sst_dtime(ni,nj))
  sst_dtime = REAL(i2buf2d,r_size)*scale_factor + offset
  WRITE(6,*) "[msg] read_viirs_nc::sst_dtime: min, max=", &
             minval(sst_dtime), maxval(sst_dtime), &
             minval(sst_dtime, mask=valid),maxval(sst_dtime, mask=valid)

!-------------------------------------------------------------------------------
! Read the Quality Key (5 is best quality; 0 missing data)
! 5    == best_quality
! 4-2  == not used
! 1    == bad data
! 0    == missing data
!-------------------------------------------------------------------------------
  ALLOCATE(quality_level(ni,nj))
  CALL nc_rdvar2d(fid, "quality_level", quality_level)
  CALL nc_rdatt(fid, "quality_level", "_FillValue", FillValue)
  if (dodebug) then
     WRITE(6,*) "_FillValue=", FillValue
     ! check num of obs for each quality level
     do i = -1, 5
        i1buf2d = 0
        if (i==-1) then
           qlev = FillValue 
        else
           qlev = i
        end if
        where (quality_level == qlev)
            i1buf2d = 1
        end where
        n = sum(int(i1buf2d,4))
        WRITE(6,*) "[msg] read_viirs_nc::quality_level, # of obs, %=", i, n, n*100.0/(ni*nj)
     end do
  end if
  where (quality_level == FillValue)
      valid = .false.
  end where
  !WRITE(6,*) "[msg] read_viirs_nc::quality_level: min, max=", &
  !           minval(quality_level, mask=valid), maxval(quality_level, mask=valid)
 
!-------------------------------------------------------------------------------
! Close the netcdf file
!-------------------------------------------------------------------------------
  CALL nc_close_fid(fid)

!-------------------------------------------------------------------------------
! CONVERT netcdf data to argo_data format:
!-------------------------------------------------------------------------------
  WRITE(6,*) "Finished reading netCDF file, formatting data..."

! determine num of valid obs
  WRITE(6,*) "[msg] read_viirs_nc::use obs with min_quality_level=", min_quality_level
  nobs = 0
  do j = 1, nj; do i = 1, ni
    if (quality_level(i,j) >= min_quality_level.and.valid(i,j)) then
       nobs = nobs + 1
    end if
  end do; end do
  WRITE(6,*) "[msg] read_viirs_nc::nobs_retained, %=", nobs, nobs*100.0/(ni*nj)
     
! fill into data struct
  ALLOCATE(obs_data(nobs))
  n = 0
  do j = 1, nj
     do i = 1, ni
        if (quality_level(i,j) >= min_quality_level .and. valid(i,j)) then
           n = n + 1
           obs_data(n)%typ      = id_sst_obs
           obs_data(n)%x_grd(1) = alon2d(i,j)
           obs_data(n)%x_grd(2) = alat2d(i,j)
           obs_data(n)%hour     = sst_dtime(i,j)/3600.0
           obs_data(n)%value    = sea_surface_temperature(i,j)
           obs_data(n)%oerr     = stde(i,j)
           obs_data(n)%qkey     = quality_level(i,j)
        end if
     end do
  end do
  CALL inspect_obs_data(obs_data)
  if (n/=nobs) then
     WRITE(6,*) "[err] read_viirs_nc::n/=nobs: n, nobs=", n, nobs
     STOP (26)
  end if

END SUBROUTINE read_viirs_nc

SUBROUTINE inspect_obs_data(obs_data)
  IMPLICIT NONE
  TYPE(sst_viirs_data),INTENT(IN) :: obs_data(:)

  WRITE(6,*) "[msg] read_viirs_nc::info"
  WRITE(6,*) "                nobs=", size(obs_data)
  WRITE(6,*) "  x_grd(1): min, max=", minval(obs_data(:)%x_grd(1)), maxval(obs_data(:)%x_grd(1))
  WRITE(6,*) "  x_grd(2): min, max=", minval(obs_data(:)%x_grd(2)), maxval(obs_data(:)%x_grd(2))
  WRITE(6,*) "      hour: min, max=", minval(obs_data(:)%hour), maxval(obs_data(:)%hour)
  WRITE(6,*) "     value: min, max=", minval(obs_data(:)%value), maxval(obs_data(:)%value)
  WRITE(6,*) "      oerr: min, max=", minval(obs_data(:)%oerr), maxval(obs_data(:)%oerr)
  WRITE(6,*) "      qkey: min, max=", minval(obs_data(:)%qkey), maxval(obs_data(:)%qkey)

END SUBROUTINE


ELEMENTAL FUNCTION cvt_temp_K2C(tf) RESULT (tc)
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: tf   ! temperature in Kelvin
  REAL(r_size)            :: tc   ! temperature in degree Celsius

  REAL(r_size),PARAMETER :: ADD_OFFSET_K2C = -273.15_r_size

  tc = tf + ADD_OFFSET_K2C

END FUNCTION

END MODULE read_viirs
