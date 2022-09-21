MODULE read_sss
  USE common,       ONLY: r_sngl, r_dble, r_size, slen
  USE params_obs,   ONLY: id_sss_obs
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sss_data
  PUBLIC :: read_jpl_smap_sss_h5
  PRIVATE :: inspect_obs_data

  TYPE sss_data
    REAL(r_size) :: x_grd(2)  ! longitude, latitude
    REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
    REAL(r_size) :: oerr      ! observation standard error
    REAL(r_size) :: hour      ! Hour of observation
    INTEGER :: qkey           ! Quality key
    INTEGER :: typ            ! type of observation (elem)
    LOGICAL :: kept           ! tells letkf whether this obs is kept for assimilation
  END TYPE sss_data

CONTAINS

SUBROUTINE read_jpl_smap_sss_h5(obsinfile,min_quality_level,obs_data,nobs)
  USE m_h5io,    ONLY: h5_get_fid, h5_close_fid, h5_rdvarshp, h5_rdvar2d, &
                       h5_rdvar1d, h5_rdatt, &
                       HID_T, HSIZE_T
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: obsinfile
  INTEGER,     INTENT(IN) :: min_quality_level
  TYPE(sss_data),ALLOCATABLE,INTENT(INOUT) :: obs_data(:)
  INTEGER,       INTENT(OUT)   :: nobs
  
  INTEGER(HID_T) :: fid
  INTEGER :: istat

!-------------------------------------------------------------------------------
! Open the hdf5 file 
!-------------------------------------------------------------------------------
  write(6,*) "[msg] read_jpl_smap_sss_h5 :: obsinfile=",trim(obsinfile)

  call h5_get_fid(trim(obsinfile), fid)
  print*, "fid=", fid
  call h5_close_fid(fid)

  nobs = 1

END SUBROUTINE 


SUBROUTINE inspect_obs_data(obs_data)
  IMPLICIT NONE
  TYPE(sss_data),INTENT(IN) :: obs_data(:)

  WRITE(6,*) "[msg] read_sss_nc::info"
  WRITE(6,*) "                nobs=", size(obs_data)
  WRITE(6,*) "  x_grd(1): min, max=", minval(obs_data(:)%x_grd(1)), maxval(obs_data(:)%x_grd(1))
  WRITE(6,*) "  x_grd(2): min, max=", minval(obs_data(:)%x_grd(2)), maxval(obs_data(:)%x_grd(2))
  WRITE(6,*) "      hour: min, max=", minval(obs_data(:)%hour), maxval(obs_data(:)%hour)
  WRITE(6,*) "     value: min, max=", minval(obs_data(:)%value), maxval(obs_data(:)%value)
  WRITE(6,*) "      oerr: min, max=", minval(obs_data(:)%oerr), maxval(obs_data(:)%oerr)
  WRITE(6,*) "      qkey: min, max=", minval(obs_data(:)%qkey), maxval(obs_data(:)%qkey)

END SUBROUTINE


END MODULE read_sss
