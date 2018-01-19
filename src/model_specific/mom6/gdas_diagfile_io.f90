MODULE gdas_diagfile_io
!
! Include the I/O for the GDAS netcdf-formatted diag files here
!
IMPLICIT NONE

PUBLIC :: read_obs_nc, write_obs_nc

PRIVATE

CONTAINS

!=============================================================================== ===============================================
! NETCDF Obs I/O following NCEP/NASA GDAS netcdf diagfile format
!=============================================================================== ===============================================

!STEVE: adding for GDAS-type netcdf format diag file
SUBROUTINE read_obs_nc(infile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
  USE netcdf
  CHARACTER(*),INTENT(IN) :: infile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: ohx(nn)
  REAL(r_size),INTENT(OUT) :: obhr(nn)
  INTEGER,INTENT(OUT) :: oqc(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: i,j,k,n
  INTEGER :: ncid,istat,varid
  INTEGER :: iunit,iolen,irec
  REAL(r_size), DIMENSION(nn) :: vardata

  !-----------------------------------------------------------------------------
  ! Open netcdf file
  !-----------------------------------------------------------------------------
  WRITE(6,*) "read_obs_nc:: opening file: ",infile
  call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "read_obs_nc :: just opened file ", infile

  !-----------------------------------------------------------------------------
  ! Read observation info
  !-----------------------------------------------------------------------------

  varname='Observation_Class'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  elem = vardata

! varname='Observation_Type'
! CALL read_obs_nc_var(varname,vardata,ncid,prec_in)

! varname='Observation_Subtype'
! CALL read_obs_nc_var(varname,vardata,ncid,prec_in)

  !-----------------------------------------------------------------------------
  ! Read position fields
  !-----------------------------------------------------------------------------

  varname='Longitude'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  rlon = vardata

  varname='Latitude'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  rlat = vardata

  varname='Depth'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  rlev = vardata

  !-----------------------------------------------------------------------------
  ! Read observation data
  !-----------------------------------------------------------------------------

  varname='Observation'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  odat = vardata

  varname='Input_Observation_Error'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  oerr = vardata

  varname='QC_Flag'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  oqc = vardata

  varname='Obs_Time'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  obhr = vardata

  !-----------------------------------------------------------------------------
  ! Read 'model equivalent' of observation (for a single ensemble member)
  !-----------------------------------------------------------------------------

  varname='Hxb'
  CALL read_obs_nc_var(varname,vardata,ncid,prec_in)
  ohx = vardata
  
  !-----------------------------------------------------------------------------
  ! Close the netcdf file
  !-----------------------------------------------------------------------------
  call check( NF90_CLOSE(ncid) )

END SUBROUTINE read_obs_nc

SUBROUTINE read_obs_nc_var(varname,vardata,ncid,prec_in)
  CHARACTER(*), INTENT(IN) :: varname
  REAL(r_size), DIMENSION(:), INTENT(OUT) :: vardata
  INTEGER, INTENT(IN) :: ncid
  INTEGER, INTENT(IN), OPTIONAL :: prec_in = 1
  
  select case(prec)
    case(1)
      ALLOCATE(buf4(nlon,nlat,nlev))
    case(2)
      ALLOCATE(buf8(nlon,nlat,nlev))
  end select

  call check( NF90_INQ_VARID(ncid,trim(varname),varid) )

  if (dodebug) WRITE(6,*) "read_diag:: just got data for variable :: ",trim(varname)
  select case (prec)
    case(1)
      buf4=0.0
      call check( NF90_GET_VAR(ncid,varid,buf4) )
      vardata = buf4
    case(2)
      buf8=0.0d0
      call check( NF90_GET_VAR(ncid,varid,buf8) )
      vardata = buf8
  end select

  if (dodebug) WRITE(6,*) "read_diag:: finished processing data for variable :: ",trim(varname)

END SUBROUTINE read_obs_nc_var


SUBROUTINE write_obs_nc(infile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,qcflag_in)
  CHARACTER(*),INTENT(IN) :: infile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: odat(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: ohx(nn)
  REAL(r_size),INTENT(IN) :: obhr(nn)
  LOGICAL, INTENT(IN), OPTIONAL :: qcflag_in
  LOGICAL :: qcflag
  INTEGER,INTENT(IN) :: oqc(nn)
  REAL(r_sngl) :: wk(9)
  INTEGER :: n,iunit
  LOGICAL, PARAMETER :: dodebug=.false.
  INTEGER :: cmode

  if (PRESENT(qcflag_in)) then
    qcflag = qcflag_in
  else
    qcflag = .false.
  endif

  !-----------------------------------------------------------------------------
  ! Create and Open netcdf file
  !-----------------------------------------------------------------------------
  cmode = NF90_HDF5
  WRITE(6,*) "write_obs_nc:: creating file: ",infile
  CALL NF90_CREATE(infile, cmode, ncid)
  WRITE(6,*) "write_obs_nc:: opening file: ",infile
  call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
  WRITE(6,*) "write_obs_nc :: just opened file ", infile

! if (qcflag .and. oqc(n)==0) CYCLE  !STEVE: find a way to eliminate entries that were qc'd out
!                                    !       This is especially useful for regional models

  !-----------------------------------------------------------------------------
  ! Write observation info
  !-----------------------------------------------------------------------------
  ! ID for observation type
  varname = 'Observation_Class'
  vardata = elem
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  !-----------------------------------------------------------------------------
  ! Write position fields
  !-----------------------------------------------------------------------------

  ! Ob lon
  varname = 'Longitude'
  vardata = rlon
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  ! Ob lat
  varname = 'Latitude'
  vardata = rlat
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  ! Ob lev
  varname = 'Depth'
  vardata = rlev
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  !-----------------------------------------------------------------------------
  ! Write observation data
  !-----------------------------------------------------------------------------

  ! Observed data quantity
  varname = 'Observation'
  vardata = odat
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  ! Estimated observation error
  varname = 'Input_Observation_Error'
  vardata = oerr
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)
   
  ! Model forecast transformed to observation space: H(xb)
  varname = 'Hxb'
  vardata = oerr
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)
   
  ! Quality control ID (1==keep, 0==discard) for use in assimilation
  varname = 'QC_Flag'
  vardata = oerr
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  varname = 'Obs_Time'
  vardata = obhr
  CALL write_obs_nc_var(varname,vardata,ncid,prec_in)

  !-----------------------------------------------------------------------------
  ! Close the netcdf file
  !-----------------------------------------------------------------------------
  call check( NF90_CLOSE(ncid) )

    
END SUBROUTINE write_obs_nc

SUBROUTINE write_obs_nc_var(varname,vardata,ncid,prec_in)
  CHARACTER(*), INTENT(IN) :: varname
  REAL(r_size), DIMENSION(:), INTENT(IN) :: vardata
  INTEGER, INTENT(IN) :: ncid
  INTEGER, INTENT(IN), OPTIONAL :: prec_in = 1
  LOGICAL :: dodebug = .true.
  
  if (dodebug) WRITE(6,*) "write_obs_nc_var:: Get id ..."
  call check( NF90_INQ_VARID(ncid,rsrt_temp_name,varid) )

  if (dodebug) WRITE(6,*) "write_obs_nc_var:: Write data ..."
  call check( NF90_PUT_VAR(ncid,varid,vardata) )

END SUBROUTINE write_obs_nc_var

END MODULE gdas_diagfile_io

