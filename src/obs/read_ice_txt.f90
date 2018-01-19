MODULE read_ice_txt
!===============================================================================
! This program reads a sea-ice text file to convert to
! a format readable by letkf. The obsop_icefrac.f90 program uses this module to read these
! in directly.
!
! Observation errors should be computed externally
! on the temperature gradients, or read from an observation error file
! (e.g. when doing adaptive obs error)
!
!===============================================================================

USE common,                     ONLY: r_sngl, r_size, slen
USE params_obs,                 ONLY: id_cn_obs

IMPLICIT NONE

PUBLIC :: read_icetxt, ice_data

PUBLIC

INTEGER :: nobs, nobs0
INTEGER :: i,j,k,n
REAL(r_size) :: se0, seF

TYPE ice_data
  REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
  REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
  REAL(r_size) :: lev       ! grid z level
  REAL(r_size) :: oerr      ! observation standard error
  REAL(r_size) :: hour      ! Hour of observation
  CHARACTER(9) :: plat      ! Platform
  CHARACTER(3) :: ptyp      ! Profile type
  CHARACTER(3) :: sid       ! Source id
  CHARACTER(1) :: qkey      ! Quality key
  INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
  INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
  INTEGER :: id     ! id number used in observation files to identify the observation
  INTEGER :: rid    ! id of the record, in order that it is read in
  INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
  LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
END TYPE ice_data

TYPE(ice_data), ALLOCATABLE, DIMENSION(:) :: obs_data


CONTAINS


SUBROUTINE read_icetxt(infile,typ,obs_data,nobs,err)  !(infile,obs_data,nobs,typ,err)
!===============================================================================
! Read the sea-ice concentration analysis data from a prepared txt file
!===============================================================================
IMPLICIT NONE
TYPE(ice_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
INTEGER, INTENT(OUT) :: nobs
CHARACTER(*), INTENT(IN) :: infile
INTEGER, INTENT(IN) :: typ
REAL(r_size), INTENT(IN) :: err

!INTEGER, PARAMETER :: nx = 720
!INTEGER, PARAMETER :: ny = 360
INTEGER :: nx, ny
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: lon !nx, linear -0.25 0.5
REAL(r_size), ALLOCATABLE, DIMENSION(:) :: lat !ny, linear -89.75 0.5

INTEGER :: nu=53
REAL(r_size) :: val !, err
INTEGER :: n,i,j,k,icnt,oidx
CHARACTER(32) dummy1,dummy2
LOGICAL :: dodebug=.true.

print *, "Reading from file: ", infile

!Read the input file
print *, "Opening file..."
open(nu,FILE=trim(infile),ACCESS='SEQUENTIAL',FORM='FORMATTED')

print *, "Reading dimensions header..."
read (nu,*) nx,ny
print *, "dim = ", nx, ny
nobs=nx*ny
ALLOCATE(obs_data(nobs))
ALLOCATE(lon(nx),lat(ny))

! Setup lons and lats:
do i=1,nx
  lon(i) = -0.25 + (i-1)*0.5
enddo
do j=1,ny
  lat(j) = -89.75 + (j-1)*0.5
enddo

oidx=0
do j=1,ny
  do i=1,nx
!   read(nu,*) cnrev(i,j)
!   read(nu,*) cn(i,ny-j+1)
    read(nu,*) val

    oidx=oidx+1
    if (dodebug) print *, "oidx, val, err = ", oidx, val, err

    obs_data(oidx)%typ = typ
    obs_data(oidx)%x_grd(1) = lon(i)
    obs_data(oidx)%x_grd(2) = lat(ny-j+1)
    obs_data(oidx)%x_grd(3) = 0
    obs_data(oidx)%hour = 0.0d0 !hour_a*1.0d0 + (minute_a/60.0d0)
    obs_data(oidx)%value = val
    obs_data(oidx)%oerr = err
    obs_data(oidx)%rid = oidx
    obs_data(oidx)%lid = 0
  enddo
enddo

! Output number of observations:
nobs = oidx
print *, "nobs = ", nobs


END SUBROUTINE read_icetxt

END MODULE read_ice_txt
