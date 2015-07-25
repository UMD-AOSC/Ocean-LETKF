program test

  !STEVE:
  CHARACTER(30) :: outfile2
  !STEVE: for outputting restart data:
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: outdata
  INTEGER, DIMENSION(3) :: outshape
  REAL, PARAMETER :: Tmin=-4
  REAL, PARAMETER :: Tmax=50

  REAL, DIMENSION(100,200,10) :: template
 
  outfile2 = 'RESTART/ocean_temp_salt.res.nc'
  outshape = SHAPE(template)
  ALLOCATE(outdata(outshape(1),outshape(2),outshape(3)))  

  print *, "outshape = ", outshape(1), outshape(2), outshape(3)

  outdata = -5

  if (.true.) then
    where(outdata < Tmin) outdata = Tmin
    where(outdata > Tmax) outdata = Tmax
  endif

  print *, "MAXVAL(outdata) = ", MAXVAL(outdata)
  print *, "MINVAL(outdata) = ", MINVAL(outdata)

end program test
