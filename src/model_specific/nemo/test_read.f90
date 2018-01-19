PROGRAM test_read
  USE netcdf
  IMPLICIT NONE
  character(128) :: gridfile = "mesh_mask.nc"
  character(128) :: grid_nlon_name = "x"
  character(128) :: grid_nlat_name = "y"
  character(128) :: grid_nlev_name = "z"
  character(128) :: grid_lon2d_name = "nav_lon"
  character(128) :: grid_lat2d_name = "nav_lat"
  character(128) :: grid_lev_name = "nav_lev"
  integer :: ncid, dimid, varid
  integer :: nlon, nlat, nlev
  real(kind=4), dimension(:,:), allocatable :: lon2d, lat2d
  real(kind=4), dimension(:), allocatable :: lev

  ! Open file:
  call check( NF90_OPEN(gridfile,NF90_NOWRITE,ncid) )

  ! Read dimension sizes:
  WRITE(6,*) "Reading x dimension (lon)..."
  call check( NF90_INQ_DIMID(ncid,grid_nlon_name,dimid) )   ! Longitude for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlon) )
  WRITE(6,*) "nlon = ", nlon

  WRITE(6,*) "Reading y dimension (lat)..."
  call check( NF90_INQ_DIMID(ncid,grid_nlat_name,dimid) )   ! Longitude for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlat) )
  WRITE(6,*) "nlat = ", nlat

  WRITE(6,*) "Reading z dimension (lev)..."
  call check( NF90_INQ_DIMID(ncid,grid_nlev_name,dimid) )   ! Longitude for T-cell
  call check( NF90_INQUIRE_DIMENSION(ncid,dimid,len=nlev) )
  WRITE(6,*) "nlev = ", nlev

  ! Allocate memory:
  allocate(lon2d(nlon,nlat), lat2d(nlon,nlat), lev(nlev))

  ! Read values:
  WRITE(6,*) "Reading x coordinate (lon2d)..."
  call check( NF90_INQ_VARID(ncid,grid_lon2d_name,varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lon2d) )
  WRITE(6,*) "lon2d(1,1) = ", lon2d(1,1)
  WRITE(6,*) "lon2d(nlon,nlat) = ", lon2d(nlon,nlat)

  WRITE(6,*) "Reading y coordinate (lat2d)..."
  call check( NF90_INQ_VARID(ncid,grid_lat2d_name,varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lat2d) )
  WRITE(6,*) "lat2d(1,1) = ", lat2d(1,1)
  WRITE(6,*) "lat2d(nlon,nlat) = ", lat2d(nlon,nlat)

  WRITE(6,*) "Reading depth (lev)..."
  call check( NF90_INQ_VARID(ncid,grid_lev_name,varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,lev) )
  WRITE(6,*) "lev(1) = ", lev(1)
  WRITE(6,*) "lev(nlev) = ", lev(nlev)

  ! Close file:
  call check( NF90_CLOSE(ncid) )

CONTAINS

subroutine check(status)
  USE netcdf

  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check


END PROGRAM test_read
