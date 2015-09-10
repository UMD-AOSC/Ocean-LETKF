PROGRAM create_drifters

USE netcdf
! The following parameters can be called from "common" module...
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: slen=512

REAL(r_size), ALLOCATABLE, DIMENSION(:,:), SAVE :: v4d 
REAL(r_size), ALLOCATABLE, DIMENSION(:), SAVE :: drifter_ids 
INTEGER, SAVE :: num_drifters 
INTEGER, SAVE :: num_times
INTEGER, SAVE :: nv4d=3

!STEVE: add drifters here...
  ! Create new drifters netcdf drfile ('drifters_inp.nc')
  ! Output drifter positions in format that can be read by mom4p1
! CHARACTER(32) :: drfile = "drifters_test.nc"
  CHARACTER(32) :: infile = "drifters_inpc.txt"
  CHARACTER(32) :: drfile = "drifters_inp.nc"
  INTEGER :: nd_dimid, np_dimid, dimids(2)
  INTEGER :: pos_varid, ids_varid
  REAL, DIMENSION(:,:), ALLOCATABLE :: positions
  INTEGER, DIMENSION(:), ALLOCATABLE :: ids
  ! For netcdf:
  INTEGER :: ncid

  call read_dimension(infile,num_drifters,num_times)
  nd = nv4d !LUYU: most likely this fix. And it means we only read the position of each drifter.
  np = num_drifters

  ALLOCATE(v4d(np,nd))

  call read_drifters(infile,v4d)
  ALLOCATE(positions(nd,np),ids(np))

  !WRITE(6,*) "#### BEFORE####"
  !WRITE(6,*) v4d_all(1:np,num_times,1:nd)
  positions = RESHAPE(TRANSPOSE(v4d),(/nd,np/)) 
  !WRITE(6,*) "#### AFTER####"
  !WRITE(6,*) positions
  ids = drifter_ids

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  call check( nf90_create(drfile, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "nd", nd, nd_dimid) )
  call check( nf90_def_dim(ncid, "np", np, np_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ nd_dimid, np_dimid /)

  ! Define the variable.
  call check( nf90_def_var(ncid, "positions", NF90_DOUBLE, dimids, pos_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, pos_varid, "names", "lon lat depth") )
  call check( nf90_put_att(ncid, pos_varid, "units", "deg_E deg_N meters") )

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).
  call check( nf90_def_var(ncid, "ids", NF90_INT, np_dimid, ids_varid) )

  ! Add global attributes with NF90_GLOBAL
  call check( nf90_put_att(ncid, NF90_GLOBAL, "velocity_names", "u v w") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_names", "lon lat depth temp salt") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "field_units", "deg_E deg_N meters Celsius PSU") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "time_units", "seconds") )
  call check( nf90_put_att(ncid, NF90_GLOBAL, "title", "LETKF analyzed positions for drifters, for input into MOM4p1") )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the data to the file.
  call check( nf90_put_var(ncid, pos_varid, positions) )
  call check( nf90_put_var(ncid, ids_varid, ids) )

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

CONTAINS

subroutine check(status)
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check

SUBROUTINE read_dimension(infile,num_drifters,num_times)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile
  INTEGER, INTENT(OUT) :: num_drifters,num_times
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: fid = 33
  CHARACTER(slen) :: drfile
 
  print *, 'Hello from read_dimension.'
  drfile = infile ! (DRIFTERS)

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  print *, 'num_drifters=', num_drifters, 'num_times=', num_times
  close(fid)

  ALLOCATE(drifter_ids(num_drifters))

  RETURN
END SUBROUTINE read_dimension

SUBROUTINE read_drifters(infile,v4d_all)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: infile 
  REAL(r_size),INTENT(OUT) :: v4d_all(num_drifters,nv4d)
  REAL :: dlon, dlat, ddepth, dtemp, dsalt, dtime
  INTEGER :: ditime, dids
  CHARACTER(16) :: dummy_char
  CHARACTER(8*12) :: header_line
  INTEGER :: di
  INTEGER :: fid = 33, dodebug = 0
  CHARACTER(24) :: drfile

  drfile = infile
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Open the XYZ drifters positions file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *, 'Hello from read_drifters.'

  ! Open the DRIFTERS file postprocessed from mom4p1 netcdf output files
  open(fid,FILE=drfile,ACCESS='sequential')
  read(fid,'(A16,I8,A16,I8)')  dummy_char, num_drifters, dummy_char, num_times
  read(fid,*) header_line

  ! Read all positions (and possibly temp and salt)  of each drifter:
  DO di=1,num_drifters
      read(fid,'(I12,6F12.4,I12)') dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      IF (dodebug .eq. 1) THEN
        print *, dids, dlon, dlat, ddepth, dtemp, dsalt, dtime, ditime
      END IF
      drifter_ids(di) = dids
      v4d_all(di,1) = dlon
      v4d_all(di,2) = dlat
      v4d_all(di,3) = ddepth
      IF (nv4d .ge. 5) THEN
        ! If we have temperature and salinity observations at each position,
        ! we can assimilate this data too
        v4d_all(di,4) = dtemp
        v4d_all(di,5) = dsalt
      END IF
  END DO
  close(fid)

END SUBROUTINE read_drifters

END PROGRAM create_drifters
