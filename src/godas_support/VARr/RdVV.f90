  program RdVV
!
  integer :: imx, jmx, kass
  integer :: year, month, day
  real, allocatable, dimension(:,:,:) :: buf
!
  open(11, form='unformatted', status='old', access='sequential')
  read(11, err=10) year, month, day
  read(11, err=11) imx, jmx, kass
  allocate (buf(imx,jmx,kass))
  read(11, err=12) buf
!
  write(6,'(i4,1x,i2,1x,i2)') year, month, day
  write(6,'(3i4)') imx, jmx, kass
!
  i=imx/2
  j=jmx/2
  k=kass/2
  write(6,'(3i4,1pe15.3)') i,j,k, buf(i,j,k)
!
  i=imx/4
  j=jmx/4
  k=kass/4
  write(6,'(3i4,1pe15.3)') i,j,k, buf(i,j,k)
!
  i=3*i
  j=3*j
  k=3*k
  write(6,'(3i4,1pe15.3)') i,j,k, buf(i,j,k)

  call exit

10 write(6,'(a)') 'Error reading date'
  call exit

11 write(6,'(a)') 'Error reading dimensions'
  call exit

12 write(6,'(a)') 'Error reading field'
  call exit

  end program RdVV
