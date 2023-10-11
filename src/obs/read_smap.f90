!===============================================================================
! This module includes subs to read in JPL SMAP L2/L3 salinity data
!
! Level-2 data are from:
! - Data access:
!   https://opendap.jpl.nasa.gov/opendap/hyrax/SalinityDensity/smap/L2/JPL/V5.0/contents.html
! - Data documentation webpage:
!   https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L2B_SSS_CAP_V5
!
! Level-3 data are from:
! - Data access:
!   https://opendap.jpl.nasa.gov/opendap/SalinityDensity/smap/L3/JPL/V5.0/8day_running/
! - Data documentation webpage:
!   https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5
!-------------------------------------------------------------------------------

MODULE read_smap
  USE common,       ONLY: r_sngl, r_dble, r_size, slen
  USE params_obs,   ONLY: id_sss_obs
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sss_data
  PUBLIC :: read_jpl_smap_l2_sss_h5
  PUBLIC :: read_jpl_smap_l3_sss_nc
  PRIVATE :: inspect_obs_data

  TYPE sss_data
    REAL(r_size) :: x_grd(2)  ! longitude, latitude
    REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
    REAL(r_size) :: oerr      ! observation standard error
    REAL(r_size) :: hour      ! Hour of observation
    INTEGER :: typ            ! type of observation (elem)
    LOGICAL :: kept           ! tells letkf whether this obs is kept for assimilation
  END TYPE sss_data

CONTAINS

!-------------------------------------------------------------------------------
! Read JPL SMAP Level-3 8-day-runing-meaning 
!
! - Data access:
!   https://opendap.jpl.nasa.gov/opendap/SalinityDensity/smap/L3/JPL/V5.0/8day_running/
! - Data documentation webpage:
!   https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5
!-------------------------------------------------------------------------------
SUBROUTINE read_jpl_smap_l3_sss_nc(obsinfile, obs_data, nobs)
  USE m_ncio,   ONLY: nc_get_fid, nc_close_fid, nc_rddim, nc_rdatt, &
                      nc_rdvar1d, nc_rdvar2d
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: obsinfile
  TYPE(sss_data),ALLOCATABLE,INTENT(INOUT) :: obs_data(:)
  INTEGER,     INTENT(OUT) :: nobs

  REAL(r_size),PARAMETER :: MAX_ALLOWED_LAND_FRACTION = 0.001_r_size ! 0.0-1.0
  REAL(r_size),PARAMETER :: MAX_ALLOWED_ICE_FRACTION  = 0.001_r_size ! 0.0-1.0
  INTEGER :: fid
  INTEGER :: nlon, nlat, nt
  REAL(r_size),ALLOCATABLE :: alon1d(:), alat1d(:) 
  REAL(r_size),ALLOCATABLE :: sea_surface_salinity(:,:), stde(:,:) ! (nlon,nlat)
  REAL(r_size),ALLOCATABLE :: land_fraction(:,:), ice_fraction(:,:) ! (nlon,nlat)
  REAL(r_size),ALLOCATABLE :: sss_time(:) ! (nt)
  LOGICAL,     ALLOCATABLE :: valid(:,:) ! (nlon,nlat)
  REAL(r_size) :: rFillValue
  LOGICAL :: dodebug = .true.
  INTEGER :: obsdate(8), sdate(8)
  REAL :: tdelta(5)
  INTEGER :: n, i ,j 


!-------------------------------------------------------------------------------
! Open the nc file 
!-------------------------------------------------------------------------------
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::obsinfile=", trim(obsinfile)
  CALL nc_get_fid(trim(obsinfile), fid)
  
!-------------------------------------------------------------------------------
! Read dimension & geo vars
!-------------------------------------------------------------------------------
  CALL nc_rddim(fid, "longitude", nlon)
  CALL nc_rddim(fid, "latitude", nlat)
  CALL nc_rddim(fid, "time",     nt)
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::nlon,nlat,nt=", nlon, nlat, nt

  ALLOCATE(alon1d(nlon), alat1d(nlat))
  CALL nc_rdvar1d(fid, "longitude", alon1d)
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::lon: min, max=", &
             minval(alon1d), maxval(alon1d)
  CALL nc_rdvar1d(fid, "latitude",  alat1d)
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::lat: min, max=", &
             minval(alat1d), maxval(alat1d)

!-------------------------------------------------------------------------------
! Read obs time info
!-------------------------------------------------------------------------------
  ALLOCATE(sss_time(nt))
  CALL nc_rdvar1d(fid, "time", sss_time)
  sdate = [2015, 1, 1, 0, 0, 0, 0, 0] ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
  tdelta = [0.0,0.0,0.0,REAL(sss_time(1)),0.0] ! DAY/HR/MIN/SEC/MSEC
  call w3movdat(tdelta, sdate, obsdate)
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::obsdate=", &
             obsdate(1:3),obsdate(5:7)

!-------------------------------------------------------------------------------
! Read sea surface salinity and its uncertainty
!-------------------------------------------------------------------------------
  ALLOCATE(valid(nlon,nlat))
  valid = .true.

  ALLOCATE(sea_surface_salinity(nlon,nlat))
  CALL nc_rdvar2d(fid, "smap_sss", sea_surface_salinity)
  CALL nc_rdatt(fid,   "smap_sss", "_FillValue", rFillValue)
  if (dodebug) WRITE(6,*) "rFillValue=", rFillValue
  where (NINT(sea_surface_salinity)==NINT(rFillValue)) 
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::smap_sss: min, max=", &
             minval(sea_surface_salinity, mask=valid), &
             maxval(sea_surface_salinity, mask=valid)

  ALLOCATE(stde(nlon,nlat))
  CALL nc_rdvar2d(fid, "smap_sss_uncertainty", stde)
  CALL nc_rdatt(fid,   "smap_sss_uncertainty", "_FillValue", rFillValue)
  if (dodebug) WRITE(6,*) "rFillValue=", rFillValue
  where (NINT(stde)==NINT(rFillValue)) 
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::smap_sss_uncertainty: min, max=", &
             minval(stde, mask=valid), &
             maxval(stde, mask=valid)

!-------------------------------------------------------------------------------
! Read land & ice fraction for additional QC
!-------------------------------------------------------------------------------
  ALLOCATE(land_fraction(nlon,nlat))
  CALL nc_rdvar2d(fid, "land_fraction", land_fraction)
  CALL nc_rdatt(fid,   "land_fraction", "_FillValue", rFillValue)
  if (dodebug) WRITE(6,*) "rFillValue=", rFillValue
  where (NINT(land_fraction)==NINT(rFillValue)) 
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::land_fraction: min, max=", &
             minval(land_fraction, mask=valid), &
             maxval(land_fraction, mask=valid)

  ALLOCATE(ice_fraction(nlon,nlat))
  CALL nc_rdvar2d(fid, "ice_fraction", ice_fraction)
  CALL nc_rdatt(fid,   "ice_fraction", "_FillValue", rFillValue)
  if (dodebug) WRITE(6,*) "rFillValue=", rFillValue
  where (NINT(ice_fraction)==NINT(rFillValue)) 
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc::ice_fraction: min, max=", &
             minval(ice_fraction, mask=valid), &
             maxval(ice_fraction, mask=valid)

!-------------------------------------------------------------------------------
! Close the nc file 
!-------------------------------------------------------------------------------
  CALL nc_close_fid(fid)

!-------------------------------------------------------------------------------
! Perform QC based on land & ice fraction 
! [FIXME]: Need to determine the best fraction 
!-------------------------------------------------------------------------------
  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc:: remove obs with land_fraction >", &
             MAX_ALLOWED_LAND_FRACTION*100, "%"
  where (land_fraction >= MAX_ALLOWED_LAND_FRACTION)
    valid = .false.
  end where

  WRITE(6,*) "[msg] read_jpl_smap_l3_sss_nc:: remove obs with ice_fraction >", &
             MAX_ALLOWED_ICE_FRACTION*100, "%"
  where (ice_fraction >= MAX_ALLOWED_ICE_FRACTION)
    valid = .false.
  end where

!-------------------------------------------------------------------------------
! Convert NC data to the data format required by obsop
!-------------------------------------------------------------------------------
  WRITE(6,*) "Finished reading NC file, formatting data..."
! determine num of valid obs
  nobs = 0
  do j = 1, nlat; do i = 1, nlon
    if (valid(i,j)) then
       nobs = nobs + 1
    end if
  end do; end do
  WRITE(6,*) "[msg] read_jpl_smap_l3_smap_h5::nobs_retained, %=", nobs, nobs*100.0/(nlon*nlat)

! fill into data struct
  ALLOCATE(obs_data(nobs))
  n = 0
  do j = 1, nlat
     do i = 1, nlon
        if (valid(i,j)) then
           n = n + 1
           obs_data(n)%typ      = id_sss_obs
           obs_data(n)%x_grd(1) = alon1d(i)
           obs_data(n)%x_grd(2) = alat1d(j)
           obs_data(n)%hour     = 0.0  ! [FIXME]: need to determine later what hours to use here
           obs_data(n)%value    = sea_surface_salinity(i,j)
           obs_data(n)%oerr     = stde(i,j)
        end if
     end do
  end do
  CALL inspect_obs_data(obs_data,subname="read_jpl_smap_l3_sss_h5")
  if (n/=nobs) then
     WRITE(6,*) "[err] read_jpl_smap_l2_sss_h5::n/=nobs: n, nobs=", n, nobs
     STOP (27)
  end if
END SUBROUTINE 

!-------------------------------------------------------------------------------
! Read JPL SMAP Level-2 SSS swath data
!
! - Data access:
!   https://opendap.jpl.nasa.gov/opendap/hyrax/SalinityDensity/smap/L2/JPL/V5.0/contents.html
! - Data documentation webpage:
!   https://podaac.jpl.nasa.gov/dataset/SMAP_JPL_L2B_SSS_CAP_V5
!-------------------------------------------------------------------------------
SUBROUTINE read_jpl_smap_l2_sss_h5(obsinfile, obs_data, nobs, Syyyymmddhh, delta_seconds)
  USE m_h5io,    ONLY: h5_get_fid, h5_close_fid, h5_rdvarshp, h5_rdvar2d, &
                       h5_rdvar1d, h5_rdatt, &
                       HID_T, HSIZE_T, i4, r4
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: obsinfile
  TYPE(sss_data),ALLOCATABLE,INTENT(INOUT) :: obs_data(:)
  INTEGER,       INTENT(OUT)   :: nobs
  CHARACTER(10),INTENT(IN), OPTIONAL :: Syyyymmddhh
                                        !1234567890
  INTEGER,      INTENT(IN), OPTIONAL :: delta_seconds  ! qc'ed if abs(obstime-Syyyymmddhh)>delta_seconds
  
  INTEGER(HID_T) :: fid
  INTEGER(HSIZE_T),allocatable :: varsizes(:)  ! (2)
  REAL(r_size),ALLOCATABLE :: alon2d(:,:), alat2d(:,:)              ! (nlines,npixels)
  REAL(r_size),ALLOCATABLE :: sea_surface_salinity(:,:), stde(:,:)  ! (nlines,npixels)
  REAL(r_size),ALLOCATABLE :: sss_time(:)                          ! (nlines) 
  INTEGER,ALLOCATABLE :: sss_time_in_seconds_since19780101(:)      ! (nlines)
  INTEGER(i4), ALLOCATABLE :: quality_flag(:,:)                     ! (nlines,npixels), need to ensure bits>16. Use 32bits here.
  LOGICAL,     ALLOCATABLE :: valid(:,:) ! .true. if no missing info at this grid by checking different FillValues in ncfile
  REAL(r4)     :: r4FillValue ! need to ensure bits=32
  INTEGER(i4)  :: i4FillValue ! need to ensure bits>16. Use 32bits here
  INTEGER :: sss_sdate(8), sss_edate(8)
  INTEGER :: sdate(8), sdate_in_seconds_since19780101, sdate_in_min_since19780101
  INTEGER :: qc_sdate(8), qc_edate(8)
  INTEGER :: qcdate(5), qcdate_in_seconds_since19780101, qcdate_in_min_since19780101
  REAL :: tdelta(5)
  INTEGER :: nlines, npixels, i, j, n, istat
  LOGICAL :: dodebug = .true.
  INTEGER,PARAMETER :: QUAL_FLAG_SSS_USABLE_GOOD = 0
  INTEGER,PARAMETER :: QUAL_FLAG_SSS_USABLE_BAD  = 1
  INTEGER,PARAMETER :: QUAL_FLAG_SSS_HAS_LAND    = 1
  INTEGER,PARAMETER :: QUAL_FLAG_SSS_HAS_ICE     = 1

!-------------------------------------------------------------------------------
! Open the hdf5 file 
!-------------------------------------------------------------------------------
  WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::read obs from obsinfile=",trim(obsinfile)
  CALL h5_get_fid(trim(obsinfile), fid)

!-------------------------------------------------------------------------------
! get dim info
!-------------------------------------------------------------------------------
  CALL h5_rdvarshp(fid, "/lon", varsizes)
  nlines = varsizes(1); npixels = varsizes(2)
  WRITE(6,*) "nlines, npixels=", nlines, npixels

  ALLOCATE(alon2d(nlines,npixels), alat2d(nlines,npixels))
  ALLOCATE(valid(nlines,npixels))
  valid = .true.
!-------------------------------------------------------------------------------
! Read row time
!-------------------------------------------------------------------------------
   sdate = [2015, 1, 1, 0, 0, 0, 0, 0] ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
   CALL W3FS21(sdate(1:5), sdate_in_min_since19780101)
   sdate_in_seconds_since19780101 = 60 * sdate_in_min_since19780101

   ALLOCATE(sss_time(nlines))
   ALLOCATE(sss_time_in_seconds_since19780101(nlines))

   CALL h5_rdvar1d(fid, "/row_time", sss_time) ! elapse seconds since sdate
   tdelta = [0.0,0.0,0.0,REAL(sss_time(1)),0.0] ! DAY/HR/MIN/SEC/MSEC
   CALL w3movdat(tdelta, sdate, sss_sdate)
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::start date=", sss_sdate(1:3),sss_sdate(5:7)
   tdelta = [0.0,0.0,0.0,REAL(sss_time(nlines)),0.0] ! DAY/HR/MIN/SEC/MSEC
   CALL w3movdat(tdelta, sdate, sss_edate)
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::end date=", sss_edate(1:3),sss_edate(5:7)

   sss_time_in_seconds_since19780101(:) = sdate_in_seconds_since19780101 + &
                                          CEILING(sss_time(:))
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::sss_dtime_since19780101: min, max=", &
             minval(sss_time_in_seconds_since19780101), &
             maxval(sss_time_in_seconds_since19780101)

  ! time filtering
  if (PRESENT(Syyyymmddhh) .and. PRESENT(delta_seconds)) then
     ! convert qc center time
     qcdate = 0
     READ(Syyyymmddhh,"(I4,I2,I2,I2)") qcdate(1),qcdate(2),qcdate(3),qcdate(4)
     WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::Syyyymmddhh=", Syyyymmddhh, &
                "qcdate(:)=",qcdate(:)
     call W3FS21(qcdate,qcdate_in_min_since19780101)
     qcdate_in_seconds_since19780101 = 60 * qcdate_in_min_since19780101
     WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::qcdate_in_seconds_since19780101, delta_seconds=", &
                qcdate_in_seconds_since19780101, delta_seconds
     
     sdate =  [1978, 1, 1, 0, 0, 0, 0, 0] ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
     tdelta = [0.0, &
                REAL((qcdate_in_seconds_since19780101-delta_seconds)/3600), &
                0.0, &
                REAL(mod(qcdate_in_seconds_since19780101-delta_seconds,3600)), &
                0.0]
     CALL w3movdat(tdelta, sdate, qc_sdate)
     tdelta = [0.0,&
               REAL((qcdate_in_seconds_since19780101+delta_seconds)/3600), &
               0.0,&
               REAL(mod(qcdate_in_seconds_since19780101+delta_seconds,3600)),&
               0.0]
     CALL w3movdat(tdelta, sdate, qc_edate)
     WRITE(6,*) "[DEBUG]: qc_sdate=", qc_sdate(1:3), qc_sdate(5:7)
     WRITE(6,*) "[DEBUG]: qc_edate=", qc_edate(1:3), qc_edate(5:7)

     ! calculate time difference
     do i = 1, nlines
        if (ABS(sss_time_in_seconds_since19780101(i) - &
               qcdate_in_seconds_since19780101) > delta_seconds) then
           valid(i,:) = .false.
        end if
     end do
  end if


!-------------------------------------------------------------------------------
! Read lat & lon info
!-------------------------------------------------------------------------------
  CALL h5_rdvar2d(fid, "/lon", alon2d)
  CALL h5_rdatt(fid,   "/lon", "_FillValue", r4FillValue)
  if (dodebug) WRITE(6,*) "_FillValue=", NINT(r4FillValue)
  where (NINT(alon2d)==NINT(r4FillValue))
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::lon: min, max=", &
             minval(alon2d,mask=valid), maxval(alon2d,mask=valid)

  CALL h5_rdvar2d(fid, "/lat", alat2d)
  CALL h5_rdatt(fid,   "/lat", "_FillValue", r4FillValue)
  if (dodebug) WRITE(6,*) "_FillValue=", NINT(r4FillValue)
  where (NINT(alat2d)==NINT(r4FillValue))
    valid = .false.
  end where
  WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::lat: min, max=", &
             minval(alat2d,mask=valid), maxval(alat2d,mask=valid)


!-------------------------------------------------------------------------------
! Read sea surface salinity & its uncertainty
!-------------------------------------------------------------------------------
   ALLOCATE(sea_surface_salinity(nlines,npixels), stde(nlines,npixels))
   CALL h5_rdvar2d(fid, "/smap_sss", sea_surface_salinity)
   CALL h5_rdatt(fid,   "/smap_sss", "_FillValue", r4FillValue)
   where (NINT(sea_surface_salinity)==NINT(r4FillValue))
      valid = .false.
   end where
   CALL h5_rdatt(fid, "/smap_sss", "valid_max", r4FillValue)
   WRITE(6,*) "valid_max_yo=", r4FillValue
   where(sea_surface_salinity > r4FillValue)
      valid = .false.
   end where
   CALL h5_rdatt(fid, "/smap_sss", "valid_min", r4FillValue)
   WRITE(6,*) "valid_min_yo=", r4FillValue
   where( sea_surface_salinity < r4FillValue)
      valid = .false.
   endwhere
 
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::smap_sss: min, max=", &
              minval(sea_surface_salinity, mask=valid), &
              maxval(sea_surface_salinity, mask=valid)

   CALL h5_rdvar2d(fid, "/smap_sss_uncertainty", stde)
   CALL h5_rdatt(fid,   "/smap_sss_uncertainty", "_FillValue", r4FillValue)
   where (NINT(stde)==NINT(r4FillValue)) ! 
      valid = .false.
   end where
   CALL h5_rdatt(fid, "/smap_sss_uncertainty", "valid_max", r4FillValue)
   WRITE(6,*) "valid_max_error=", r4FillValue
   where(stde > r4FillValue)
      valid = .false.
   end where
   CALL h5_rdatt(fid, "/smap_sss_uncertainty", "valid_min", r4FillValue)
   WRITE(6,*) "valid_min_error=", r4FillValue
   where(stde < r4FillValue)
      valid = .false.
   end where
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::smap_sss_uncertainty: min, max=", &
              minval(stde, mask=valid), maxval(stde, mask=valid)
 
!-------------------------------------------------------------------------------
! Read quality control flag
! Based on
! https://podaac-tools.jpl.nasa.gov/drive/files/allData/smap/docs/%20JPL-CAP_V5/SMAP-SSS_JPL_V5.0_Documentation.pdf
!
! Bit | Definition                        | Bit Significance Text
! 0   | QUAL_FLAG_SSS_USABLE              | 0 - Overall SSS quality good
!                                         |  1 - Overall SSS quality bad
! 1   | QUAL_FLAG_FOUR_LOOKS              |  0 - Data from all four looks available
!                                         | 1 - Data from all four looks not available
! 2   | QUAL_FLAG_POINTING                | 0 - Nominal incidence angles (within 0.2◦ of 40◦)
!                                         | 1 - Non-nominal incidence angles
! 4   | QUAL_FLAG_LARGE_GALAXY_CORRECTION | 0 - All galaxy corrections < 5 K
!                                         | 1 - At least one TB flavor had galaxy correction > 5 K
! 5   | QUAL_FLAG_ROUGHNESS_CORRECTION    | 0 - Ancillary wind speed < 20 m/s
!                                         | 1 - Ancillary wind speed > 20 m/s
! 6   | QUAL_FLAG_SST_TOO_COLD            | 0 - SST > 5 deg C
!                                         | 1 - SST < 5 deg C
! 7   | QUAL_FLAG_LAND                    | 0 - No land detected in SWC
!                                         | 1 - Land detected in SWC
! 8   | QUAL_FLAG_ICE                     | 0 - No ice detected in SWC
!                                         | 1 - ice detected in SWC
! 9   | QUAL_FLAG_HIGH_SPEED_USABLE       | 0 - Overall high speed quality good
!                                         | 1 - Overall high speed quality bad
!
! Bits 3 and 10-15 are reserved for possible future use.
!
! Note Fortran does not support unsigned integer, so we will use 32-bit signed 
! integer (i.e., INTEGER(4)) to contains values from 16-bit unsigned integer in 
! the file. This ensures the bit value from the data does not change through fortran
! reading. 
! For example, If we use INTEGER(2) to read 65535_uint16, then the outputs of IBITS
! will be 0111111111111111, which is incorrect. 
! If using INTEGER(4) to read 65535_uint16, then the ouputs with IBITS will be
! 1111111111111111, which is correct.
! 
! Use IBITS(var,pos,len) to get the bit of the integer. Bit is counted from the 
! the right side. Alignments of bit position are shown below.
!
! for 16-bits (2-bytes)             |15|...|3|2|1|0|
! for 32-bits (4-bytes) |31|30|..|16|15|...|3|2|1|0|
!-------------------------------------------------------------------------------
   ALLOCATE(quality_flag(nlines,npixels))
   CALL h5_rdvar2d(fid, "/quality_flag", quality_flag)
   CALL h5_rdatt(fid,   "/quality_flag", "_FillValue", i4FillValue)
   if (dodebug) WRITE(6,*) "i4FillValue", i4FillValue
   where (quality_flag == i4FillValue .or. &
          IBITS(quality_flag,0,1)==QUAL_FLAG_SSS_USABLE_BAD .or. &   ! overall bad
          IBITS(quality_flag,7,1)==QUAL_FLAG_SSS_HAS_LAND .or. &   ! has land
          IBITS(quality_flag,8,1)==QUAL_FLAG_SSS_HAS_ICE ) !      ! has ice
     valid = .false.
   end where
   WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::quality_flag: min, max=", &
              minval(quality_flag, mask=valid), maxval(quality_flag, mask=valid)

!-------------------------------------------------------------------------------
! Close the hdf5 file 
!-------------------------------------------------------------------------------
  CALL h5_close_fid(fid)

!-------------------------------------------------------------------------------
! Convert H5 data to the data format required by obsop
!-------------------------------------------------------------------------------
  WRITE(6,*) "Finished reading HDF5 file, formatting data..."
! determine num of valid obs
  nobs = 0
  do j = 1, npixels; do i = 1, nlines
    if (valid(i,j)) then
       nobs = nobs + 1
    end if
  end do; end do
  WRITE(6,*) "[msg] read_jpl_smap_l2_sss_h5::nobs_retained, %=", nobs, nobs*100.0/(npixels*nlines)

! fill into data struct
  ALLOCATE(obs_data(nobs))
  n = 0
  do j = 1, npixels
     do i = 1, nlines
        if (valid(i,j)) then
           n = n + 1
           obs_data(n)%typ      = id_sss_obs
           obs_data(n)%x_grd(1) = alon2d(i,j)
           obs_data(n)%x_grd(2) = alat2d(i,j)
           obs_data(n)%hour     = sss_time_in_seconds_since19780101(i)/3600.
           obs_data(n)%value    = sea_surface_salinity(i,j)
           obs_data(n)%oerr     = stde(i,j)
        end if
     end do
  end do
  CALL inspect_obs_data(obs_data,subname="read_jpl_smap_l2_sss_h5")
  if (n/=nobs) then
     WRITE(6,*) "[err] read_jpl_smap_l2_sss_h5::n/=nobs: n, nobs=", n, nobs
     STOP (27)
  end if
END SUBROUTINE 


SUBROUTINE inspect_obs_data(obs_data, subname)
  IMPLICIT NONE
  TYPE(sss_data),INTENT(IN) :: obs_data(:)
  CHARACTER(*),INTENT(IN) :: subname

  WRITE(6,*) "[msg] "//trim(subname)//"::info"
  WRITE(6,*) "                nobs=", size(obs_data)
  WRITE(6,*) "  x_grd(1): min, max=", minval(obs_data(:)%x_grd(1)), maxval(obs_data(:)%x_grd(1))
  WRITE(6,*) "  x_grd(2): min, max=", minval(obs_data(:)%x_grd(2)), maxval(obs_data(:)%x_grd(2))
  WRITE(6,*) "      hour: min, max=", minval(obs_data(:)%hour), maxval(obs_data(:)%hour)
  WRITE(6,*) "     value: min, max=", minval(obs_data(:)%value), maxval(obs_data(:)%value)
  WRITE(6,*) "      oerr: min, max=", minval(obs_data(:)%oerr), maxval(obs_data(:)%oerr)

END SUBROUTINE


END MODULE read_smap
