!------------------------------------------------------------
! CFS common observations
! This contains the definitions for observations used by the
! CFSv2-LETKF. For strongly coupled data assimilation, observation
! ID's are required to be unique between the ocean and atmosphere
!
! Author: Travis Sluka
!------------------------------------------------------------
module common_obs

  USE common

  IMPLICIT NONE

  PUBLIC

  !!------------------------------------------------------------
  !! Observations
  !!------------------------------------------------------------

  !!------------------------------------------------------------
  !! unique ID's for observations
  !!------------------------------
  INTEGER, PARAMTER :: obsid_num = 16

  !! atmosphere obs
  INTEGER, PARAMTER :: obsid_atm_min  = 1000
  INTEGER, PARAMTER :: obsid_atm_max  = 1999
  INTEGER, PARAMTER :: obsid_atm_num  = 8
  INTEGER, PARAMTER :: obsid_atm_offset = 0 
  
  INTEGER, PARAMTER :: obsid_atm_ps   = 1100
  INTEGER, PARAMTER :: obsid_atm_rain = 1110
  INTEGER, PARAMTER :: obsid_atm_t    = 1210
  INTEGER, PARAMTER :: obsid_atm_tv   = 1211
  INTEGER, PARAMTER :: obsid_atm_q    = 1220
  INTEGER, PARAMTER :: obsid_atm_rh   = 1221
  INTEGER, PARAMTER :: obsid_atm_u    = 1250
  INTEGER, PARAMTER :: obsid_atm_v    = 1251

  !! ocean obs
  INTEGER, PARAMTER :: obsid_ocn_min  = 2000
  INTEGER, PARAMTER :: obsid_ocn_max  = 2999
  INTEGER, PARAMTER :: obsid_ocn_num  = 8
  INTEGER, PARAMTER :: obsid_ocn_offset = 8

  INTEGER, PARAMTER :: obsid_ocn_ssh  = 2100
  INTEGER, PARAMTER :: obsid_ocn_eta  = 2101
  INTEGER, PARAMTER :: obsid_ocn_sst  = 2110
  INTEGER, PARAMTER :: obsid_ocn_sss  = 2120
  INTEGER, PARAMTER :: obsid_ocn_t    = 2210
  INTEGER, PARAMTER :: obsid_ocn_s    = 2220
  INTEGER, PARAMTER :: obsid_ocn_u    = 2250
  INTEGER, PARAMTER :: obsid_ocn_v    = 2251

  !! arrays holding all observation id's and names, for easy iteration
  !! in loops that want to print stats for obs
  INTEGER :: obsid_ids(obsid_num) = (/&
       obsid_atm_ps, obsid_atm_rain, obsid_atm_t, obsid_atm_tv, &
       obsid_atm_q, obsid_atm_rh, obsid_atm_u, obsid_atm_v, &
       obsid_ocn_ssh, obsid_ocn_eta, obsid_ocn_sst, obsid_ocn_sss, &
       obsid_ocn_t, obsid_ocn_s, obsid_ocn_u, obsid_ocn_v/)
  CHARACTER (len=10) :: obsid_names(obsid_num) = (/&
       "ATM_PS  ", "ATM_RAIN", "ATM_T   ", "ATM_TV  ", &
       "ATM_Q   ", "ATM_RH  ", "ATM_U   ", "ATM_V   ", &
       "OCN_SSH ", "OCN_ETA ", "OCN_SST ", "OCN_SSS ",&
       "OCN_T   ", "OCN_S   ", "OCN_U   ", "OCN_V   "/)


  !! ------------------------------------------------------------
  !! Structure to hold LETKF formated observations,
  !! see the wiki for further documentation. This is the
  !! common format produced by the observation operators
  !! ------------------------------------------------------------
  TYPE :: Observation
     INTEGER      :: id     !! one of the obsid_* values from above
     REAL(r_size) :: lon    !! longitude (degrees)
     REAL(r_size) :: lat    !! latitude (degrees)
     REAL(r_size) :: lev    !! height/depth (mb, meters)
     REAL(r_size) :: odat   !! observation value
     REAL(r_size) :: oerr   !! estimate observation error
     INTEGER      :: platform !! platform ID (currently not used)
     
     REAL(r_size) :: time   !! time realtive to analysis
     REAL(r_size) :: ohx    !! model value in observation space
     INTEGER      :: qc     !! quality control (1 = valid, 0 = invalid)
  END TYPE Observation


  !!------------------------------------------------------------
  !!------------------------------------------------------------
  !! Platform types, used mainly for QC, atmosphere only now
  !!------------------------------------------------------------
  !!------------------------------------------------------------
  INTEGER, PARAMTER :: platform_num = 21
  CHARACTER (6), PARAMTER :: platform_name(platform_num)= (/&
       'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR'/)

  
CONTAINS

  

  !! ------------------------------------------------------------
  !! given an INTEGER observation ID, returns its location in the
  !! above obsid_ids and obdsid_name arrays (and presumably
  !! from any user specified arrays involving obsids as well)
  !! returns -1 if not found
  !! ------------------------------------------------------------
  FUNCTION obsid_id2idx(id)
    INTEGER, intent(in) ::id
    INTEGER ::  obsid_id2idx
    INTEGER n

    obsid_id2idx = -1
    do n=1,obsid_num
       if (obsid_ids(n) .eq. id) then
          obsid_id2idx = n
          return
       endif
    end do
  END FUNCTION obsid_id2idx



  !! ------------------------------------------------------------
  !! ------------------------------------------------------------
  FUNCTION obs_getnum(file, extended) result(nobs)
    CHARACTER(*), intent(in) :: file
    logical, intent(in), optional :: extended
    logical :: ext
    INTEGER :: nobs, unit, ios, recl, n
    REAL(r_sngl),allocatable :: wk(:)
    INTEGER :: obs_count(obsid_num),inv_count
    
    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    endif
    allocate(wk(recl))

    !! open the file and count the number of observations
    nobs = 0
    obs_count(:) = 0
    inv_count = 0
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do
       read(unit, iostat=ios) wk
       if (ios /= 0) exit  !! end of the file found, exit loop
       !! another record found, increment count
       if ( nint(wk(1)) .gt. 0) then          
          nobs = nobs + 1
          n = obsid_id2idx( nint(wk(1)) )
          if (n < 0 ) then
             inv_count = inv_count + 1
          else
             obs_count(n) = obs_count(n) +1
          endif
          
       else
          inv_count = inv_count + 1
       endif
       
    end do

    !!done, cleanup
    deallocate(wk)
    close(unit)

    !! print statistics
    write (6,*) ""    
    write (6,*) "--------------------------"
    write (6,*) "Observation count check"
    write (6,*) "File: ",file
    write (6,*) "--------------------------"
    do n=1,obsid_num
       if (obs_count(n) > 0) then
          write (6,*) obsid_names(n), obs_count(n)
       endif
    end do
    write (6,*) "--------------------------"    
    write (6,*) "Total Valid: ",nobs
    write (6,*) "Invalid:     ",inv_count
    write (6,*) ""
    
  END FUNCTION obs_getnum


  
  !! ------------------------------------------------------------
  !! reads in an LETKF formatted observation file.
  !!  This file contains a series of observations each with
  !!  7 values, or 10 values (if an extended format)
  !!  The short format is what is read into the obsop programs.
  !!  The extended format is read into/out of the LETKF program
  !!  Memory for "obs" will be allocated with this subroutine.
  !!  When the user is finished with obs, they should deallocate it.
  !! ------------------------------------------------------------
  SUBROUTINE obs_read(file, obs, extended)
    CHARACTER(*), intent(in) :: file
    TYPE(Observation), allocatable, intent(out) :: obs(:)
    logical, intent(in), optional :: extended
    REAL(r_sngl),allocatable :: wk(:)
    logical :: ex, ext
    INTEGER :: count, unit, ios, n, recl

    !! make sure the file exists
    inquire(file=file, exist=ex)
    if ( .not. ex) then
       write (6,*) file, ' does not exist -- skipped'
    endif

    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    endif
    allocate(wk(recl))

    !! count the number of records, and allocate space
    count = obs_getnum(file, ext)
    write (6,*) 'Reading in ', count,' observations'
    allocate( obs(count) )

    !! read in the observations
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do n=1,count
       read(unit) wk
       obs(n)%id        = nint(wk(1))
       obs(n)%lon       = wk(2)
       obs(n)%lat       = wk(3)
       obs(n)%lev       = wk(4)
       obs(n)%odat      = wk(5)
       obs(n)%oerr      = wk(6)
       obs(n)%platform  = wk(7)
       if (ext) then
          obs(n)%time   = wk(8)
          obs(n)%ohx    = wk(9)
          obs(n)%qc     = wk(10)
       else
          obs(n)%time   = 0
          obs(n)%ohx    = 0
          obs(n)%qc     = 0
       endif
       if (obsid_id2idx(obs(n)%id) < 0) then
          write (6,*) "ERROR: Invalid observation ID", obs(n)%id
          stop 1
       endif
    end do

    write (*,*) size(obs)

    !! all done, close up
    deallocate(wk)
    close(unit)
  END SUBROUTINE obs_read


  !! ------------------------------------------------------------
  !! Writes an array of observations into a file with the LETKF format
  !!  See "read_obs" for further comments.
  !! ------------------------------------------------------------
  SUBROUTINE obs_write(file, obs, extended)
    CHARACTER(*), intent(in) :: file
    TYPE(Observation), intent(in) :: obs(:)
    logical, intent(in), optional :: extended
    REAL(r_sngl),allocatable :: wk(:)
    logical :: ex, ext
    INTEGER :: recl, n, unit

    !! see if file already exists
    inquire(file=file, exist=ex)
    if (ex) then
       write(6,*) "WARNING: ",file," is being overwritten"
    endif

    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    endif
    allocate(wk(recl))

    !! open file and start writing records to it
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do n=1,size(obs)
       wk(1)  = obs(n)%id
       wk(2)  = obs(n)%lon
       wk(3)  = obs(n)%lat
       wk(4)  = obs(n)%lev
       wk(5)  = obs(n)%odat
       wk(6)  = obs(n)%oerr
       wk(7)  = obs(n)%platform
       if (ext) then
          wk(8)  = obs(n)%time
          wk(9)  = obs(n)%ohx
          wk(10) = obs(n)%qc
       endif
       write(unit) wk
    end do
    close(unit)
    
  END SUBROUTINE obs_write

  
  !!------------------------------------------------------------
  !! converts from an array of Observations types to several
  !! arrays of each variables. mainly for dealing with legacy code.
  !!------------------------------------------------------------
  !! TODO, try to get rid of references to this in the code
  SUBROUTINE obs_new2old(obs, elem,rlon,rlat,rlev,odat,oerr,otyp)
    TYPE(Observation),intent(in) :: obs(:)
    REAL(r_size), intent(inout) :: elem(:)
    REAL(r_size), intent(inout) :: rlon(:)
    REAL(r_size), intent(inout) :: rlat(:)
    REAL(r_size), intent(inout) :: rlev(:)
    REAL(r_size), intent(inout) :: odat(:)
    REAL(r_size), intent(inout) :: oerr(:)
    REAL(r_size), intent(inout) :: otyp(:)

    INTEGER n

    do n=1,size(obs)
       elem(n) = obs(n)%id
       rlon(n) = obs(n)%lon
       rlat(n) = obs(n)%lat
       rlev(n) = obs(n)%lev
       odat(n) = obs(n)%odat
       oerr(n) = obs(n)%oerr
       otyp(n) = obs(n)%platform

       !! convert pressure for atmosphere to Pa
       select case(obs(n)%id)
       case(obsid_atm_u, obsid_atm_v, obsid_atm_t,&
            obsid_atm_tv,obsid_atm_q)
          rlev(n) = rlev(n) * 100.0
       case ( obsid_atm_rh)
          rlev(n) = rlev(n) * 100.0
          odat(n) = odat(n) * 0.01
          oerr(n) = oerr(n) * 0.01
       case(obsid_atm_ps)
          odat(n) = odat(n) * 100.0
          oerr(n) = oerr(n) * 100.0
       end select       
    end do   
  END SUBROUTINE obs_new2old
  
end module common_obs
