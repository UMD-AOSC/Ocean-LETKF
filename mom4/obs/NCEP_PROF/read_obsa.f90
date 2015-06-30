PROGRAM read_g4p1_mom
! STEVE: 
! Read the godas input observation/covariance files

INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_dble=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, PARAMETER :: ni=720, nj=410, nk=40
REAL, DIMENSION(ni,nj,nk) :: buf
CHARACTER(64) :: filename
REAL, DIMENSION(ni) :: lons
REAL, DIMENSION(nj) :: lats
REAL, DIMENSION(nk) :: levs
REAL, DIMENSION(ni,nj,nk) :: v3d, e3d

!!STEVE: read tvv.mom file
!filename='tvv.mom'
!CALL read_vv(filename,buf)
!!STEVE: convert to tvv.grd grads file
!filename='tvv.grd'
!CALL write_vv(filename,buf)

!!STEVE: and for salinity,
!filename='svv.mom'
!CALL read_vv(filename,buf)
!filename='svv.grd'
!CALL write_vv(filename,buf)

!STEVE: read temperature profile observations
!filename='tmpa.mom'
!CALL read_obsa_stn(filename) !STEVE: trying to use grads station format, not very successful

CALL init_grid_spec
filename='tmpa.mom'
CALL read_obsa(filename,v3d,e3d)
filename='tmpa.grd'
CALL write_obsgrd(filename,v3d,e3d)

filename='sala.mom'
CALL read_obsa(filename,v3d,e3d)
filename='sala.grd'
CALL write_obsgrd(filename,v3d,e3d)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_vv(filename,buf)
CHARACTER(*),INTENT(IN) :: filename
REAL(r_sngl), DIMENSION(ni,nj,nk), INTENT(OUT) :: buf
INTEGER :: nu = 54
REAL :: tbvrf = 1.0
INTEGER :: ig, jg, kg
INTEGER :: year, month, day

!call mpp_open(nu,file_name,action=MPP_RDONLY,form=MPP_IEEE32, &
!                       access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)

print *, "Reading from file: ", filename

open(nu,FILE=trim(filename),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

print *, "Reading year, month, day = "
read (nu) year, month, day
print *, year, month, day
print *, "Reading ig, jg, kg = "
read (nu) ig, jg, kg
print *, ig,jg,kg
print *, "Reading buf (sample) = " 
buf=0.0
read (nu) buf(:,:,1:kg)
print *, buf(1:10,1:10,1)
close(nu)

return

END SUBROUTINE read_vv

SUBROUTINE write_vv(filename,buf)
CHARACTER(*),INTENT(IN) :: filename
REAL, DIMENSION(ni,nj,nk), INTENT(IN) :: buf
REAL(r_sngl) :: buf4(ni,nj)
INTEGER :: iunit,iolen
INTEGER :: k,n,irec

iunit=55
INQUIRE(IOLENGTH=iolen) iolen
OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=ni*nj*iolen)

print *, "Writing to file: ", filename

irec=1
DO k=1,nk
  buf4 = 0.0
  buf4 = REAL(buf(:,:,k),r_sngl)
  print *, "Writing level :: ", k
  WRITE(iunit,REC=irec) buf4
  irec = irec + 1
END DO

CLOSE(iunit)

END SUBROUTINE write_vv

SUBROUTINE read_obsa_stn(filename)
IMPLICIT NONE
CHARACTER(*),INTENT(IN) :: filename
INTEGER :: nsgobs = 1
CHARACTER(64) :: obfile
INTEGER :: nu, nwu
INTEGER :: pe = 1
INTEGER :: plti, pltf
INTEGER :: year, month, day, hour, minute, kd
REAL, DIMENSION(nk) :: val, err
INTEGER :: n, k, icnt
!
INTEGER :: date
!
CHARACTER*8 :: stnid
INTEGER*4 :: nlev, nflag
REAL*4 :: lat, lon, time
REAL*4 :: outval, outerr

!call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
!                       access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)

print *, "Reading from file: ", filename

!Read the input file
open(nu,FILE=trim(filename),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!Write out the data in a format that can be read in grads
open(21, file='Tz.stn', form='unformatted', status='unknown',access='stream')

! - stnr has to be integer*8
! - all other integer values have to be *4 (single precision)
! - all real values have to be *4 (single precision)

read (nu) icnt
do n=1,icnt
    val = 0.0
    err = 0.0
    read (nu) plti, pltf
    read (nu) year, month, day, hour, minute, kd
    read (nu) lon, lat
!   print *, "lat, lon = ", lat, lon
    read (nu) (val(k), err(k), k=1,kd)
!   print *, "val(1), err(1) = ", val(1), err(1)

    time = 0.0 
!   print *, "time = ", time

!   WRITE THE START HEADER
    nlev=1  ! the number of level-dependent groups
    nflag=1 ! 1 := no surface variables following header
    write(stnid,'(i8)') n
    print *, "stnid,lat,lon,time,nlev,nflag = ", stnid,lat,lon,time,nlev,nflag
    write(21) stnid,lat,lon,time,nlev,nflag
!   WRITE THE VARIABLES
!   val(1) = 999.999
    outval = 1.0
    write(21) outval !,err(1)
!   write(21) val,err
enddo

!WRITE THE TERMINATE TAIL  
nlev=0
write(21) stnid,lat,lon,time,nlev,nflag

close(nu)
close(21)

END SUBROUTINE read_obsa_stn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_obsa(filename, v3d, e3d)
IMPLICIT NONE
CHARACTER(*), INTENT(IN) :: filename
REAL, DIMENSION(ni,nj,nk), INTENT(OUT) :: v3d, e3d
INTEGER :: nu=53
REAL, DIMENSION(nk) :: val, err
INTEGER :: plti, pltf
INTEGER :: year, month, day, hour, minute, kd
INTEGER :: n,i,j,k,icnt
INTEGER, DIMENSION(ni,nj,nk) :: cnt3d
REAL :: lon, lat

!call
!mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
!                       access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)

print *, "Reading from file: ", filename

!Read the input file
print *, "Opening file..."
open(nu,FILE=trim(filename),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

print *, "Reading icnt..."
read (nu) icnt
print *, "icnt = ", icnt
do n=1,icnt
    val = 0.0
    err = 0.0
!   STEVE: the old version needs these:
!   print *, "Reading plti, pltf..."
!   read (nu) plti, pltf
    print *, "plti, pltf = ", plti, pltf
    read (nu) year, month, day, hour, minute, kd
    print *, "Reading year, month, day, hour, minute, kd ..."
    print *, "year, month, day, hour, minute, kd = ", year, month, day, hour, minute, kd
    print *, "Reading lon, lat"
    read (nu) lon, lat
    print *, "lon, lat = ", lon, lat
    print *, "Reading val, err..."
    read (nu) (val(k), err(k), k=1,kd)
    print *, "val = ", val
    print *, "err = ", err
!   print *, "val(1), err(1) = ", val(1), err(1)
    
    !find the grid cell this belongs to
    do i=1,ni-1
      if (lon < lons(i+1) ) then
        print *, "i, lons(i), lon, lons(i+1) = ", i, lons(i), lon, lons(i+1) 
        exit
      endif
    enddo

    do j=1,nj-1
      if (lat < lats(j+1) ) then
        print *, "j, lats(j), lat, lats(j+1) = ", j, lats(j), lat, lats(j+1)
        exit
      endif
    enddo

    !add the value to compute mean
    print *, "n = ", n
    do k=1,kd
      if (val(k) < 99) then !STEVE: if not NaN or empty value
        v3d(i,j,k) = v3d(i,j,k) + val(k)
        e3d(i,j,k) = e3d(i,j,k) + 1/SQRT(err(k))
        cnt3d(i,j,k) = cnt3d(i,j,k) + 1
        print *, "v3d(i,j,k) = ", v3d(i,j,k)
        print *, "cnt3d(i,j,k) = ", cnt3d(i,j,k)
      endif
    enddo
enddo

! Compute the mean field
where (cnt3d > 1) v3d = v3d / cnt3d
where (cnt3d > 1) e3d = e3d / cnt3d

END SUBROUTINE read_obsa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_obsgrd(filename, v3d, e3d)
IMPLICIT NONE
! Write out obs data to a grads binary file
CHARACTER(*),INTENT(IN) :: filename
REAL, DIMENSION(ni,nj,nk), INTENT(IN) :: v3d, e3d
REAL(r_sngl) :: buf4(ni,nj)
INTEGER :: iunit,iolen
INTEGER :: k,n,irec

iunit=55
INQUIRE(IOLENGTH=iolen) iolen
OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=ni*nj*iolen)

print *, "Writing to file: ", filename

irec=1
DO k=1,nk
  buf4 = 0.0
  buf4 = REAL(v3d(:,:,k),r_sngl)
  print *, "Writing level :: ", k
  WRITE(iunit,REC=irec) buf4
  irec = irec + 1
END DO

DO k=1,nk
  buf4 = 0.0
  buf4 = REAL(e3d(:,:,k),r_sngl)
  print *, "Writing level :: ", k
  WRITE(iunit,REC=irec) buf4
  irec = irec + 1
END DO

CLOSE(iunit)

END SUBROUTINE write_obsgrd

SUBROUTINE init_grid_spec
  USE netcdf
  INTEGER :: ncid,istat,varid,dimid
  CHARACTER(NF90_MAX_NAME) :: dimname

  WRITE(6,'(A)') '  >> accessing file: grid_spec.nc'
  call check( NF90_OPEN('grid_spec.nc',NF90_NOWRITE,ncid) )
  call check( NF90_INQ_VARID(ncid,'grid_x_T',varid) )   ! Longitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lons) )
  WRITE(6,*) "lon(1) = ", lons(1)
  WRITE(6,*) "lon(nlon) = ", lons(ni)
  call check( NF90_INQ_VARID(ncid,'grid_y_T',varid) )   ! Latitude for T-cell
  call check( NF90_GET_VAR(ncid,varid,lats) )
  WRITE(6,*) "lat(1) = ", lats(1)
  WRITE(6,*) "lat(nlat) = ", lats(nj)
  call check( NF90_INQ_VARID(ncid,'zt',varid) )      ! depth of T-cell
  call check( NF90_GET_VAR(ncid,varid,levs) )
  WRITE(6,*) "lev(1) = ", levs(1)
  WRITE(6,*) "lev(nlev) = ", levs(nk)


END SUBROUTINE init_grid_spec

subroutine check(status)
  USE netcdf
  IMPLICIT NONE
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

END PROGRAM read_g4p1_mom
