PROGRAM test_io_subgrid

USE common
USE common_oceanmodel
USE input_nml_oceanmodel, ONLY: read_input_namelist

IMPLICIT NONE

REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:,:) :: v3d
REAL(r_sngl), ALLOCATABLE, DIMENSION(:,:,:)   :: v2d
INTEGER :: i,j,k
CHARACTER(slen) :: infile, outfile

!-----------------------------------------------------------------------------
! Read in namelist parameters
!-----------------------------------------------------------------------------
WRITE(6,*) "test_io_subgrid:: calling read_input_namelist..."
CALL read_input_namelist

!-----------------------------------------------------------------------------
! Set model grid and variables
!-----------------------------------------------------------------------------
WRITE(6,*) "test_io_subgrid:: calling set_common_oceanmodel..."
CALL set_common_oceanmodel

!-----------------------------------------------------------------------------
! Read a sample forecast restart file
!-----------------------------------------------------------------------------
infile = 'gues'
WRITE(6,*) "test_io_subgrid:: calling read_restart..."
CALL read_restart(trim(infile),v3d,v2d,1)

i=12
j=12
k=1
print *, "v3d(i,j,k,1) = ", v3d(i,j,k,1)

!-----------------------------------------------------------------------------
! Write a sample analysis patch in an template restart file
!-----------------------------------------------------------------------------
v3d = v3d + 1
v2d = v2d + 1
outfile='anal'
WRITE(6,*) "test_io_subgrid:: calling write_restart..."
CALL write_restart(trim(outfile),v3d,v2d)


END PROGRAM test_io_subgrid
