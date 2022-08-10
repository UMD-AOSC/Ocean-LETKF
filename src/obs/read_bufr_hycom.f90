MODULE read_bufr
! Module contains read subroutine for TESAC/BUFR data for the ocean:
!
! Available subroutines for bufr here:
! http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/toc/
!
! Authors:
! Steve Penny
! Jack Woollen

  USE common,                     ONLY: r_sngl, r_size, slen, com_interp_spline
  USE params_model,               ONLY: nlev
  USE params_obs,                 ONLY: id_t_obs, id_s_obs
  USE compute_profile_error,      ONLY: cmpTz

  IMPLICIT NONE

  REAL(r_size) :: k2c = -273.15

  PUBLIC :: read_bufr_dumpjb, bufr_data, maxobs !, obs_data

  INTEGER :: nobs, nobs0
  INTEGER :: i,j,k,n
  REAL(r_size) :: se0, seF

  TYPE bufr_data
    REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
    REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
    REAL(r_size) :: lev       ! grid z level
    REAL(r_size) :: oerr      ! observation standard error
    REAL(r_size) :: hour      ! Hour of observation
    CHARACTER(9) :: plat      ! Platform
    CHARACTER(3) :: ptyp      ! Profile type
    CHARACTER(3) :: sid       ! Source id
    CHARACTER(1) :: qkey      ! Quality key
    INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
    INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
    INTEGER :: id     ! id number used in observation files to identify the observation
    INTEGER :: rid    ! id of the record, in order that it is read in
    INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
    LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
  END TYPE bufr_data

! TYPE(bufr_data), ALLOCATABLE, DIMENSION(:), SAVE :: obs_data
  INTEGER, PARAMETER :: maxobs = 499999 !Maximum number of observations in one day

  CONTAINS

  SUBROUTINE read_bufr_dumpjb(lunit,typ,obs_data,nobs)                                                                                        

    IMPLICIT NONE  !STEVE: this isn't working with bufr subroutines

    INTEGER, INTENT(IN) :: lunit
    INTEGER, INTENT(IN) :: typ
    TYPE(bufr_data), INTENT(INOUT), DIMENSION(:) :: obs_data
    INTEGER, INTENT(OUT) :: nobs
    REAL(8) :: date(5),sloc(3),misc(4)
    INTEGER, PARAMETER :: nvars = 5
    INTEGER, PARAMETER :: nlevs = 1000
    REAL(8),DIMENSION(nvars,nlevs) :: btocn
    REAL(8), PARAMETER :: bmiss=10e10
    INTEGER :: rdate(5)
    CHARACTER(8) :: subset,stype
    REAL(r_size) :: dbss,stmp,saln,droc,spoc
    CHARACTER(9) :: plat
    REAL(r_size), DIMENSION(nlevs) :: vals, stde
    REAL(r_size) :: missing_value=-999 !10.0E10 !STEVE: specify from bufr
    REAL(r_size) :: val,err
    INTEGER :: idate,iret,nlvl,l
    INTEGER :: ireadmg,ireadsb
    LOGICAL :: dodebug=.false.
    INTEGER :: n


    !STEVE: for lack of a better approach, for now I'm just guessing at the
    !possible number of observations in one given day:
!   nobs = maxobs
!   print *, "ALLOCATING obs_data with nobs = ", nobs
!   ALLOCATE(obs_data(nobs))
    stde=0.0 !STEVE: placeholder for now. Perhaps better to do this in the
             !       calling function (e.g. obsop_tprof.f90)
    if (typ .eq. id_t_obs) then
      se0 = 1.0
      seF = 1.5
    elseif (typ .eq. id_s_obs) then
      se0 = 0.05
      seF = 0.15
    endif
 
    ! open a bufr file to read
    ! http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/toc/intro/#openbf
    !
    ! This subroutine identifies to the BUFRLIB software a BUFR file that is
    ! connected to logical unit LUBFR.
    !
    ! In the case of an existing BUFR file, [the definition of the DX BUFR tables] 
    ! may be embedded within the first few BUFR messages of the file itself, 
    ! in which case the user can denote this fact to the subroutine by setting 
    ! LUNDX to the same value as LUBFR.
    !
    ! Note that LUBFR and LUNDX are logical unit numbers;
    ! therefore, the user within his or her application program must have
    ! already associated these logical unit numbers with actual filenames on the
    ! local system, typically via a FORTRAN "OPEN" statement.
    !
    ! CALL OPENBF  ( LUBFR, CIO, LUNDX )
    CALL openbf(lunit,'IN',lunit)

    ! read a message
    ! http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/toc/intro/#readmg
    !
    ! Subroutine READMG reads the next BUFR message from the given BUFR file
    ! pointed to by LUBFR. The associated function IREADMG does the same thing,
    ! but returns IRET as its function value which can then, e.g. be directly
    ! utilized as the target variable in an iterative program loop. 
    !
    ! IRET = IREADMG  ( LUBFR, CSUBSET, IDATE )
    nobs = 0
    do while(ireadmg(lunit,subset,idate)==0)

      ! filter by subset type
      stype=''
      ! JILI modified code number
      !if(subset=='NC031001') stype='BATHY'
      !if(subset=='NC031002') stype='TESAC'
      !if(subset=='NC031003') stype='TRKOB'
      if(subset=='TESAC') stype='TESAC'
      if(subset=='BATHY') stype='BATHY'
      if(subset=='TRKOB') stype='TRKOB'
      if(stype=='') CALL bort('unknown message type '//subset)
      if (dodebug) print"(80('-'))"
      if (dodebug) print*

      ! read a subset
      ! http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/toc/intro/#readsb
      !
      ! the next step is to do the following in order to read a subset from
      ! that internal message:
      !
      ! CALL READSB  ( LUBFR, IRET )
      do while(ireadsb(lunit)==0)

        ! read the date time
        ! http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/toc/intro/#ufb
        !
        !-----------------------------------------------------------------------
        ! CALL UFBINT  ( LUBFR, R8ARR, MXMN, MXLV, NLV, CMNSTR )
        !-----------------------------------------------------------------------
        CALL ufbint(lunit,date,5,1,iret,'YEAR MNTH DAYS HOUR MINU')
        rdate=date 


        !-----------------------------------------------------------------------
        ! read the id and location - read clon/clat if clonh/clath not available
        !-----------------------------------------------------------------------
        CALL ufbint(lunit,sloc,3,1,iret,'RPID CLONH CLATH')
        if(sloc(2)>=bmiss) CALL ufbint(lunit,sloc,3,1,iret,'RPID CLON CLAT ')


        !-----------------------------------------------------------------------
        ! read the measurement indicators
        !-----------------------------------------------------------------------
        !| IDGT     | 002032 | INDICATOR FOR DIGITIZATION
        !| IWTEMP   | 022067 | INSTRUMENT TYPE FOR WATER TEMPERATURE PROFILE MEAS
        !| WTEMPR   | 022068 | WATER TEMPERATURE PROFILE RECORDER TYPES
        !| MSDM     | 002033 | METHOD OF SALINITY/DEPTH MEASUREMENT
        CALL ufbint(lunit,misc,4,1,iret,'IDGT IWTEMP WTEMPR MSDM')
!       WRITE (plat,*) misc(2)

        !-----------------------------------------------------------------------
        ! print the ob header
        !-----------------------------------------------------------------------
        if (dodebug) print'(a8,1x,i4,4i2.2,2x,a8,2f10.4,4(1x,f6.0))',stype,rdate,sloc,misc

        !-----------------------------------------------------------------------
        ! read/print sub-surface levels of data
        !-----------------------------------------------------------------------
        !| DBSS     | 007062 | DEPTH BELOW SEA/WATER SURFACE                            |
        !| STMP     | 022193 | SEA TEMPERATURE AT SPECIFIED DEPTH                       |
        !| SALN     | 022062 | SALINITY                                                 |
        !| DROC     | 022004 | DIRECTION OF CURRENT                                     |
        !| SPOC     | 022031 | SPEED OF CURRENT                                         |
        CALL ufbint(lunit,btocn,5,1000,nlvl,'DBSS STMP SALN DROC SPOC')
        if (dodebug) then
          print*
          print'(4x,5(1x,4x,a4,1x))','dbss','stmp','saln','droc','spoc'
          print* 
        endif

        !STEVE: first, cubic spline interpolate the observations to the model levels
        !INTEGER,INTENT(IN) :: ndim         ! number of grid points (i.e. number of observation levels)
        !REAL(r_size),INTENT(IN) :: x(ndim) ! coordinate
        !REAL(r_size),INTENT(IN) :: y(ndim) ! variable
        !INTEGER,INTENT(IN) :: n            ! number of targets
        !REAL(r_size),INTENT(IN) :: x5(n)   ! target coordinates (i.e. model levels)
        !REAL(r_size),INTENT(OUT) :: y5(n)  ! target values (i.e. observations interpolated to model levels)
!       CALL com_interp_spline(ndim,x,y,n,x5,y5)

        ! Identify the type of observation (e.g. temperature, salinity, etc.)
        if (typ==id_t_obs) then
          vals = btocn(2,:) + k2c ! (convert from Kelvin to Celcius)
        elseif (typ==id_s_obs) then
          vals = btocn(3,:)
        endif

        ! Estimate the observation error based on the profile
        if (nlvl>1) then
          !STEVE: I would prefer to call this in the calling function, but it's
          !       easier here while the data is still organized into profiles
          CALL cmpTz(stde(1:nlvl),se0,seF,vals(1:nlvl),btocn(1,1:nlvl),nlvl,missing_value)
        else
          stde(1) = se0 + seF
        endif

       if (sloc(2) .lt. 0.0) sloc(2)=modulo(sloc(2),360.0)

        do l=1,nlvl
          dbss = btocn(1,l)
          stmp = btocn(2,l)
          saln = btocn(3,l)
          droc = btocn(4,l)
          spoc = btocn(5,l)
          err  = stde(l)

          if (typ==id_t_obs) then
            val = stmp + k2c
          elseif (typ==id_s_obs) then
            val = saln
          endif

          if (-4 < val .and. val < 50) then
            nobs = nobs+1
            n = nobs
            if (nobs > maxobs) then
              print *, "read_bufr_dumpjb:: ERROR! Must increase size of maxobs = ", maxobs                 
              STOP (99)
            endif
          else
            if (dodebug) then
              print *, "This is not a valid observation value :: val = ", val
              print *, "Cylcing..."
            endif
            CYCLE
          endif

          if (dodebug) print'(i4,5(1x,f8.2,1x))',l,dbss,stmp,saln,droc,spoc
          obs_data(n)%typ = typ
          obs_data(n)%x_grd(1) = sloc(2)
          obs_data(n)%x_grd(2) = sloc(3)
          obs_data(n)%x_grd(3) = dbss
          obs_data(n)%hour = date(4)+date(5)/60.0d0  ! real-valued hour of day
          obs_data(n)%value = val
          obs_data(n)%oerr = err
          obs_data(n)%rid = i         ! record id
          obs_data(n)%lid = l         ! level id
          obs_data(n)%plat = plat     !platform
          obs_data(n)%ptyp = plat     !platform type
          obs_data(n)%sid  = ""       !
          obs_data(n)%qkey = ""       !quality key

          if (dodebug) print *, "obs_data(n) = ", obs_data(n)
        enddo



        if (dodebug) then
          print*
          print"(80('-'))"
          !pause
          print"(80('-'))"
        endif
      enddo

    enddo

  END SUBROUTINE read_bufr_dumpjb

END MODULE read_bufr
