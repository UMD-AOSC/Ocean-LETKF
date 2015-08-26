PROGRAM prof2letkf

 !USE common_mom4 ! Use check subroutine
 USE netcdf
 USE common, ONLY: r_size, r_sngl
 USE params_obs

 IMPLICIT NONE

!SUBROUTINE read_grid(infile,v3d,v2d)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid,dimid
  
  !INTEGER,PARAMETER :: nlon=192
  !INTEGER,PARAMETER :: nlat=189
  !INTEGER,PARAMETER :: nlev=31
  INTEGER :: nlon, nlat, nlev
  REAL(r_size), ALLOCATABLE :: lon(:)
  REAL(r_size), ALLOCATABLE :: lat(:)
  REAL(r_size), ALLOCATABLE :: lev(:)
  REAL(r_size), ALLOCATABLE :: salt(:,:,:)
  REAL(r_size), ALLOCATABLE :: temp(:,:,:)
  REAL :: time, err
  REAL(r_sngl) :: wk(6)

  CHARACTER (len=*), PARAMETER :: tsfile="ocean_temp_salt.res.nc"
  !CHARACTER (len=*), PARAMETER :: uvfile="ocean_velocity.res.nc"
  !CHARACTER (len=*), PARAMETER :: sffile="ocean_sbc.res.nc"
  CHARACTER (len=*), PARAMETER :: outfile="obsin.dat"

  call check(NF90_OPEN(tsfile,NF90_NOWRITE,ncid))
  print *, 'finish open the file: ocean_temp_salt.res.nc.'
  ! Read Dimensions
  call check(NF90_INQ_DIMID(ncid,'xaxis_1',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlon))
  call check(NF90_INQ_DIMID(ncid,'yaxis_1',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlat))
  call check(NF90_INQ_DIMID(ncid,'zaxis_1',varid))
  call check(NF90_INQUIRE_DIMENSION(ncid,varid,len=nlev))

  ALLOCATE(lon(nlon))
  ALLOCATE(lat(nlat))
  ALLOCATE(lev(nlev))
  ALLOCATE(salt(nlon,nlat,nlev))
  ALLOCATE(temp(nlon,nlat,nlev))

  ! Read Data
  call check(NF90_INQ_VARID(ncid,'xaxis_1',varid))
  call check(NF90_GET_VAR(ncid,varid,lon))
  print *, 'finish setting up xaxis_1.'
  call check(NF90_INQ_VARID(ncid,'yaxis_1',varid))
  call check(NF90_GET_VAR(ncid,varid,lat))
  print *, 'finish setting up yaxis_1.'
  call check(NF90_INQ_VARID(ncid,'zaxis_1',varid))
  call check(NF90_GET_VAR(ncid,varid,lev))
  print *, 'finish setting up zaxis_1.'
  call check(NF90_INQ_VARID(ncid,'Time',varid))
  call check(NF90_GET_VAR(ncid,varid,time))
  print *, 'finish setting up Time.'

  call check(NF90_INQ_VARID(ncid,'salt',varid))
  call check(NF90_GET_VAR(ncid,varid,salt))
  print *, 'finish setting up salt.'
  call check(NF90_INQ_VARID(ncid,'temp',varid))
  call check(NF90_GET_VAR(ncid,varid,temp))
  print *, 'finish setting up temp.' 
  call check(NF90_CLOSE(ncid))
 

  OPEN(91,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
        call RANDOM_NUMBER(err)
        wk(1)=id_s_obs
        wk(2)=lon(i)
        wk(3)=lat(j)
        wk(4)=lev(k)
        wk(5)=salt(i,j,k)
        wk(6)=err
        WRITE(91) wk
        wk(1)=id_t_obs
        wk(5)=temp(i,j,k) 
        WRITE(91) wk
        !print *, id_s_obs,lon(i),lat(j),lev(k),salt(i,j,k),err   
      END DO
    END DO
  END DO

  CLOSE(91)
!  RETURN
!END SUBROUTINE read_grid
CONTAINS

SUBROUTINE check(status)
!===============================================================================
! Check the error status of the netcdf command
!===============================================================================
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
END SUBROUTINE check
  
END PROGRAM prof2letkf
