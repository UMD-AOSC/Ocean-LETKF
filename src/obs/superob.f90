PROGRAM superob
  USE common
  USE params_model
  USE vars_model
  USE common_oceanmodel
  USE params_obs,                ONLY: nobs, id_sst_obs
  USE params_obs,                ONLY: DO_REMOVE_65N
  USE vars_obs
  USE common_obs_oceanmodel
#ifdef DYNAMIC
  USE input_nml_oceanmodel,      ONLY: read_input_namelist
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Command line inputs:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  REAL(r_size) :: obserr_scaling=1.0d0 !STEVE: use this to scale the input observations
  REAL(r_size) :: obs_randselect=1.0d0 !STEVE: use this to select random input observations

  !-----------------------------------------------------------------------------
  ! Obs data arrays
  !-----------------------------------------------------------------------------
  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  REAL(r_size), ALLOCATABLE :: obhr(:)
  INTEGER     , ALLOCATABLE :: oqc(:)

  !-----------------------------------------------------------------------------
  ! Miscellaneous
  !-----------------------------------------------------------------------------
  REAL(r_size) :: dk,tg,qg
  REAL(r_size) :: ri,rj,rk
  INTEGER :: n, i,j,k
  REAL(r_size), DIMENSION(1) :: rand
  
  LOGICAL :: remap_obs_coords = .true.


  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
  LOGICAL :: dodebug1 = .false.
  LOGICAL :: print1st = .true.

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  INTEGER :: min_quality_level=5  ! CDA
  INTEGER :: typ = id_sst_obs
  LOGICAL :: DO_SUPEROBS = .true. 
  REAL(r_size), DIMENSION(:,:), ALLOCATABLE :: superobs, delta, M2 ! for online computation of the mean and variance
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: supercnt
  INTEGER :: idx
  REAL(r_size) :: min_oerr = 0.2 !(K)

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
#ifdef DYNAMIC
  CALL read_input_namelist
#endif
  CALL set_common_oceanmodel
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)


  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  nobs = 0
  CALL get_nobs(trim(obsinfile),8,nobs)

  if (nobs<0) then ! [NOBS <= 0]

     print*, "[warning] superob :: nobs = 0. will output an empty file" 
     call create_empty_obsfile(trim(obsoutfile))

  else ! [NOBS > 0 ]

      ALLOCATE( elem(nobs) )
      ALLOCATE( rlon(nobs) )
      ALLOCATE( rlat(nobs) )
      ALLOCATE( rlev(nobs) )
      ALLOCATE( odat(nobs) )
      ALLOCATE( oerr(nobs) )
      ALLOCATE( ohx(nobs) )
      ALLOCATE( oqc(nobs) )
      ALLOCATE( obhr(nobs) )

      CALL read_obs3(trim(obsinfile), nobs, elem, rlon, rlat, rlev, odat, oerr, obhr)
      oqc(:) = 0    ! now set all obs as bad. 
      ohx(:) = 0.d0

      if (print1st) then  
        i=1
        print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
        i=nobs
        print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
        print *, "odat(1:40) = ", odat(1:40)
        print *, "oerr(1:40) = ", oerr(1:40)
      endif

      !-----------------------------------------------------------------------------
      ! Update the coordinate to match the model grid
      !-----------------------------------------------------------------------------
      if (remap_obs_coords) then
        CALL center_obs_coords(rlon,oerr,nobs)
      endif

      !-----------------------------------------------------------------------------
      ! Bin the obs and estimate the obs error based on bin standard deviations
      !-----------------------------------------------------------------------------
      if (DO_SUPEROBS) then
        ALLOCATE(superobs(nlon,nlat),delta(nlon,nlat),M2(nlon,nlat),supercnt(nlon,nlat))
        superobs=0.0; delta = 0.0; M2 = 0.0; supercnt = 0

        print *, "Computing superobs..."

        do n=1,nobs ! for each ob,
    !     if (dodebug1) print *, "n = ", n

          if (DO_REMOVE_65N .and. rlat(n) > 65) CYCLE

          !CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
          !print*, "ijk: ri, rj=", ri, rj
          CALL phys2ij_nearest(rlon(n),rlat(n),ri,rj) ! OCEAN

          if (CEILING(ri) < 1 .OR. nlon < CEILING(ri)) CYCLE
          if (CEILING(rj) < 1 .OR. nlat < CEILING(rj)) CYCLE

          i = NINT(ri); j = NINT(rj)
          supercnt(i,j) = supercnt(i,j) + 1
          delta(i,j)    = odat(n) - superobs(i,j)
          superobs(i,j) = superobs(i,j) + delta(i,j)/supercnt(i,j)
          M2(i,j)       = M2(i,j) + delta(i,j)*(odat(n) - superobs(i,j))
    !     if (dodebug1) print *, "supercnt(",i,",",j,") = ", supercnt(i,j)
        enddo
        !"superobs" contains the mean
        !M2 contains the variance:
        WHERE(supercnt > 1) M2 = M2 / (supercnt - 1)

        idx=0
        !do j=1,nlat-1
        !  do i=1,nlon-1
        do j=1,nlat
          do i=1,nlon
            if (supercnt(i,j) > 1 .and. wet(i,j)>0.5) then ! only retain ocn pts defined by MOM
              idx = idx+1
              if (dodebug1) print *, "idx = ", idx
              odat(idx) = superobs(i,j)
              oerr(idx) = obserr_scaling*(min_oerr + SQRT(M2(i,j)))
              rlon(idx) = lon(i)
              rlat(idx) = lat(j)
              rlev(idx) = 0.d0
              elem(idx) = id_sst_obs
              oqc(idx) = 1   ! set good qc as default
              obhr(idx) = 0.d0

              if (dodebug1) print *, "odat(idx) = ", odat(idx)
              if (dodebug1) print *, "oerr(idx) = ", oerr(idx)
              if (dodebug1) print *, "rlon(idx) = ", rlon(idx)
              if (dodebug1) print *, "rlat(idx) = ", rlat(idx)
              if (dodebug1) print *, "ocnt(idx) = ", supercnt(i,j)
            endif
          enddo
        enddo
        print *, "DO_SUPEROBS:: superobs reducing from ",nobs," to ",idx, " observations." 
        print *, "with min obs error = ", MINVAL(oerr(1:idx))
        print *, "with max obs error = ", MAXVAL(oerr(1:idx))
        nobs = idx
      endif

      call write_obs3(trim(obsoutfile),nobs,elem(1:nobs), &
                                            rlon(1:nobs), &
                                            rlat(1:nobs), &
                                            rlev(1:nobs), &
                                            odat(1:nobs), &
                                            oerr(1:nobs), &
                                            obhr(1:nobs), &
                                             oqc(1:nobs), &
                                           qcflag_in=.true.)

      DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr )
  
  end if ! [IF NOBS > 0]

CONTAINS

SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In obsop_sst_podaac.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case('-superob')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_SUPEROBS
    case('-scale')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) obserr_scaling
    case('-minerr')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) min_oerr
    case('-minqc')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) min_quality_level
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug1
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM superob
