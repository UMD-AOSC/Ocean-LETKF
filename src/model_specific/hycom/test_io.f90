PROGRAM test_io

USE common
USE common_hycom

REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
REAL(r_size) :: v2d(nlon,nlat,nv2d)
CHARACTER(32) :: file

CALL set_common_hycom

file = 'gues'
CALL read_diag(file,v3d,v2d,1)
file = 'gues_me.grd'
CALL write_bingrd(file,v3d,v2d)
v3d = v3d + 1
v2d = v2d + 1
file = 'anal'
CALL write_restart(file,v3d,v2d)
file = 'anal_me.grd'
CALL write_bingrd(file,v3d,v2d)


END PROGRAM test_io
