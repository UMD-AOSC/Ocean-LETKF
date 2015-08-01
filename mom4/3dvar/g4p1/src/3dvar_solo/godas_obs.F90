module godas_obs_mod
!
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module manages observations for GODAS.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module manages observations for GODAS.
! Initialization should be called after GODAS initialization.
!</DESCRIPTION>
!
use fms_mod,                    only: FATAL, WARNING, NOTE, stdout, stdlog
use mpp_mod,                    only: mpp_error, mpp_pe
use mpp_io_mod,                 only: mpp_open, mpp_close
use mpp_io_mod,                 only: MPP_WRONLY, MPP_RDONLY, MPP_IEEE32, MPP_SEQUENTIAL
use mpp_io_mod,                 only: MPP_SINGLE, MPP_MULTI
use time_manager_mod,           only: time_type, set_date, increment_time
use time_manager_mod,           only: get_date, get_time
use ocean_domains_mod,          only: get_local_indices, get_global_indices
use ocean_util_mod,             only: write_timestamp
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,            only: ocean_time_type

use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
use godas_data_mod,             only: kass, nsgobs, nsgsobs
use godas_data_mod,             only: temp_code, xbt_code, tao_code, argo_code
use godas_data_mod,             only: salt_code, stao_code, sargo_code
use godas_data_mod,             only: sst_code, sss_code
use godas_data_mod,             only: altm_code, tp_code, j1_code
use godas_data_mod,             only: num_obsz, num_obs0, num_obsa, jemx
use godas_data_mod,             only: gds_freq, tovrf, sovrf, tov0f, sov0f
use godas_data_mod,             only: obs_trk_cnt
use godas_data_mod,             only: restore_sfc
use godas_data_mod,             only: no_lat_mx
!
implicit none

private

logical :: godas_obsz_module_initialized = .false.
logical :: godas_obs0_module_initialized = .false.
logical :: godas_obsa_module_initialized = .false.

character(len=256) :: version = '$Id: godas_obs.F90,v 1.0 2007/03/20 10:08:00 gtn Exp $'
character(len=256) :: tagname = 'Tag $Name: gds2p0d $'
character(len=48), parameter          :: mod_name = 'godas_obs_mod'

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

type data_type
   character(len=3) :: gridname
   character(len=128) :: fieldname_code ! used in user's code (e.g. temp, salt)
   character(len=128) :: fieldname_file ! fieldname used in the data file (not used)
   character(len=128) :: file_name      ! name of data file
   logical :: ongrid                    ! false, not relevant, here for compatibility
   real :: factor                       ! For unit conversion, default=1
end type data_type

integer, parameter :: max_table=10

type(data_type), dimension(max_table) :: data_table

type(time_type)  :: o_time
integer          :: o_sec
real             :: rovrf

! for ascii output
integer :: unit=6

!STEVE: for debugging:
LOGICAL :: dodebug=.false.
LOGICAL :: file_exists

public  godas_obsz_init
public  godas_obs0_init
public  godas_obsa_init
public  godas_obs_track
contains


!#######################################################################
! <FUNCTION NAME="godas_obsz_init">
!
! <DESCRIPTION>
! Initialization code for profile observations, returning a pointer to
! the obs_Z array.
! </DESCRIPTION>
!
function godas_obsz_init (Grid, Domain, Time, num_obs_z, debug) &
                    result (obs_Z)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  integer, intent(out)                        :: num_obs_z

  logical, intent(in), optional               :: debug

  ! return value
  type(ocean_obsz_type), dimension(:), pointer :: obs_Z

  integer               :: n, nf, nsg, ntable, num_obs_files
  integer               :: i, j, k, icnt, nu, nwu, nru, nobs
  integer               :: year, month, day, hour, minute, kd, icode
  character(len=8)      :: plti
  character(len=2)      :: pltf
  integer               :: kda
  real(kind=4)          :: xo, yo
  real(kind=4), allocatable, dimension(:) :: val, err
  real(kind=4)          :: xs, xe, ys, ye
  real(kind=4), allocatable, dimension(:) :: xlc, ylc
  real(kind=4)          :: dci, dcj, di, dip, dj, djp
  integer               :: pe, ierr
  character(len=256)    :: record
  type(data_type)       :: default_table, data_entry

  character(len=32)     :: obfile
  integer               :: ios

  character(len=48),  parameter :: sub_name = 'godas_obsz_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  if (godas_obsz_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' GODAS ObsZ already initialized')
  endif

  pe = mpp_pe()

  ! set local array indices
  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Dom, isg, ieg, jsg, jeg)
  nk=Grd%nk

  j=(jsd+jed)/2
  xs=Grd%xt(isd,j)
  if (xs .gt. Grd%xt(isc,j)) xs = xs - 360.0
  xe=Grd%xt(ied,j)
  if (xe .lt. Grd%xt(iec,j)) xe = xe + 360.0

  allocate (xlc(isd:ied))
  j=(jsd+jed)/2
  xlc(isd) = xs
  do i=isc,iec
    xlc(i)=Grd%xt(i,j)
  enddo
  xlc(ied) = xe
  allocate (ylc(jsd:jed))
  i=(isd+ied)/2
  do j=jsd,jed
    ylc(j)=Grd%yt(i,j)
  enddo
  xs=xlc(isd)
  xe=xlc(ied)
  ys=ylc(jsd)
  ye=ylc(jed)
  if (jed .ge. jemx) then
!   ye = ylc(jemx-4)
    if (dodebug .and. jemx-4 < jsd) then
      print *, "godas_obs.F90::"
      print *, "xs = ", xs
      print *, "xe = ", xe
      print *, "isd:ied = ", isd,ied
      print *, "xlc(isd:ied) = ", xlc(isd:ied)
      print *, "ys = ", ys
      print *, "ye = ", ye
      print *, "jemx = ", jemx
      print *, "jsd:jed = ", jsd,jed
      print *, "ylc(jsd:jed) = ", ylc(jsd:jed)
    endif
!   ye = ylc(jemx-4) !STEVE: replaced with below
    if (jemx-4 < jsd) then
      ye = ylc(jsd)
    elseif (jemx-4 > jed) then
      ye = ylc(jed)
    else
      ye = ylc(jemx-4)
    endif
  endif
! xs=xlc(isc)
! xe=xlc(iec)
! ys=ylc(jsc)
! ye=ylc(jec)
! if (jsc .eq. jsg) then
!   ys = ylc(jsc+1)
! endif
! if (jec .eq. jeg) then
!   ye = ylc(jec-3)
! endif

  allocate (val(nk))
  allocate (err(nk))

  nullify(obs_Z)

  write( stdlog(),'(/a/)') trim(version)

! initialize data table
  default_table%gridname = 'none'
  default_table%fieldname_code = 'none'
  default_table%fieldname_file = 'none'
  default_table%file_name = 'none'
  default_table%ongrid = .FALSE.
  default_table%factor = 1.
  do n = 1,max_table
    data_table(n) = default_table
  enddo

! read observations table
  if (dodebug) print *, "Opening data_table..."
  call mpp_open(nu, 'data_table', action=MPP_RDONLY)
  ntable = 0
  num_obs_files = 0
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'OBS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no OBS entry in data_table')
  endif
  num_obs_files = 0
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') num_obs_files = num_obs_files + 1
    if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') num_obs_files = num_obs_files + 1
  enddo

  nobs = 0
  if (num_obs_files .gt. 0) then
!
!  could not get readable files written with mpp_open
!
!   call mpp_open(nwu,'tmpZ',action=MPP_WRONLY,form=MPP_IEEE32,access=MPP_SEQUENTIAL, &
!                       threading=MPP_MULTI,fileset=MPP_MULTI)
!
    write(obfile,'(a,i4.4)') 'tmpZ.',pe
    nwu=1000+pe
    open(unit=nwu,file=trim(obfile),status='NEW',access='SEQUENTIAL',form='UNFORMATTED', &
         action='WRITE',iostat=ios)
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') then
        INQUIRE(FILE=trim(data_table(nf)%file_name), EXIST=file_exists) !STEVE: checking for file before opening
        if (file_exists) then
          if (dodebug) print *, "godas_obs.F90:: Opening ", trim(data_table(nf)%file_name)
          call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        else
          print *, "File does not exist: ", trim(data_table(nf)%file_name)
          print *, "Exiting..."
          STOP(111)
        endif
        icode = temp_code
        do nsg=1,nsgobs
          read (nu) icnt
          do n=1,icnt
            val = 0.0
            err = 0.0
! STEVE: temporarily not reading this data as a test of new obs prep scripts
! STEVE: NOW READING this data (plti,pltf) again:
            read (nu) plti, pltf
!           if (dodebug) print *, "godas_obs.F90::read temp:: nsg,n=",nsg,n
            read (nu) year, month, day, hour, minute, kd
!           if (dodebug) print *, "year, month, day, hour, minute, kd = ", year, month, day, hour, minute, kd
            read (nu) xo, yo
!           if (dodebug) print *, "xo, yo = ", xo, yo
            read (nu) (val(k), err(k), k=1,kd)
!           if (dodebug) print *, "val,err = ", val(1:kd), err(1:kd)

            !STEVE: I'm assuming this partitions the domain into different sets for each processor:
            if (dodebug) print *, "in local? ::", (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye)
            if (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye) then
              dci = get_dec_index(xo,xlc,isd,ied)
              dcj = get_dec_index(yo,ylc,jsd,jed)
              i = int(dci)
              j = int(dcj)
              kda = min(Grd%kmt(i,j), Grd%kmt(i,j+1), Grd%kmt(i+1,j), Grd%kmt(i+1,j+1), kd)
              if (dodebug) print *, "kda = ", kda
              !STEVE: write data to new tmpZ.pe file, specifically for godas MPI processing
              if (kda .gt. 0) then
                nobs = nobs + 1
                if (dodebug) print *, "Writing temp to obfile = ", obfile, ":: nobs = ", nobs
                write (nwu) plti, pltf
                write (nwu) icode, year, month, day, hour, minute, kda
                write (nwu) xo, yo, dci, dcj
                write (nwu) (val(k), err(k), k=1,kda)
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      else if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') then
        if (dodebug) print *, "godas_obs.F90:: Opening ", trim(data_table(nf)%file_name)
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = salt_code
        do nsg=1,nsgobs
          read (nu) icnt
          do n=1,icnt
            val = 0.0
            err = 0.0
! STEVE: as above with temperature...
            read (nu) plti, pltf
            read (nu) year, month, day, hour, minute, kd
            read (nu) xo, yo
            read (nu) (val(k), err(k), k=1,kd)
            if (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye) then
              dci = get_dec_index(xo,xlc,isd,ied)
              dcj = get_dec_index(yo,ylc,jsd,jed)
              i = int(dci)
              j = int(dcj)
              kda = min(Grd%kmt(i,j), Grd%kmt(i,j+1), Grd%kmt(i+1,j), Grd%kmt(i+1,j+1), kd)
              if (kda .gt. 0) then
                nobs = nobs + 1
                if (dodebug) print *, "Writing salt to obfile = ", obfile, ":: nobs = ", nobs
                write (nwu) plti, pltf
                write (nwu) icode, year, month, day, hour, minute, kda
                write (nwu) xo, yo, dci, dcj
                write (nwu) (val(k), err(k), k=1,kda)
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      endif
    enddo
    close(nwu)
    num_obsz = nobs
    num_obs_z = nobs

! allocate obs_Z
    allocate( obs_Z (num_obsz) )

    do n=1,num_obsz
#ifndef STATIC_MEMORY
      allocate( obs_Z(n)%val(kass) )
      allocate( obs_Z(n)%err(kass) )
      allocate( obs_Z(n)%aerr(kass) )
      allocate( obs_Z(n)%inv(kass) )
      allocate( obs_Z(n)%inc(kass) )
      allocate( obs_Z(n)%bke(kass) )
#endif
      obs_Z(n)%val(:) = 0.0
      obs_Z(n)%err(:) = 0.0
      obs_Z(n)%aerr(:) = 0.0
      obs_Z(n)%inv(:) = 0.0
      obs_Z(n)%inc(:) = 0.0
      obs_Z(n)%bke(:) = 0.0
    enddo

! fill obs_Z data array
!
!  could not read files with mpp_open.  see above.
!
!   call mpp_open(nwu,'tmpZ',action=MPP_RDONLY,form=MPP_IEEE32,access=MPP_SEQUENTIAL, &
!                       threading=MPP_MULTI,fileset=MPP_MULTI)
!
    nru=1000+pe
    if (dodebug) print *, "godas_obs.F90:: Opening ", trim(obfile)
    open(unit=nru,file=trim(obfile),status='OLD',access='SEQUENTIAL',form='UNFORMATTED', &
         position='REWIND',action='READ',iostat=ios)
    do n=1,num_obsz
!STEVE: temporary
      if (dodebug) print *, "ln370-402dbg:: n = ", n, '/', num_obsz
      read (nru) plti, pltf
      if (dodebug) print *, "plti, pltf = ", plti, pltf
      read (nru) icode, year, month, day, hour, minute, kd
      if (dodebug) print *, "icode, year, month, day, hour, minute, kd", icode, year, month, day, hour, minute, kd
      read (nru) xo, yo, dci, dcj
      if (dodebug) print *, "xo, yo, dci, dcj = ", xo, yo, dci, dcj
      read (nru) (val(k), err(k), k=1,kd)
      if (dodebug) print *, "val = ", val(1:kd)
      if (dodebug) print *, "err = ", err(1:kd)
      if (icode .eq. temp_code) then
        obs_Z(n)%name = 'temp'
        obs_Z(n)%units = 'deg-C'
        obs_Z(n)%code = temp_code
        rovrf = 1.0 / tovrf
      else if (icode .eq. salt_code) then
        obs_Z(n)%name = 'salt'
        obs_Z(n)%units = 'psu'
        obs_Z(n)%code = salt_code
        rovrf = 1.0 / sovrf
      else
        call mpp_error(FATAL,'observation code does not match known values')
      endif
      obs_Z(n)%platid = plti
      obs_Z(n)%platform = pltf
      o_sec = (hour*60 + minute)*60
      o_time = set_date(year,month,day,0,0,0)
      obs_Z(n)%obs_time = increment_time(o_time,o_sec,0)
      obs_Z(n)%kd  = kd
      obs_Z(n)%olng  = xo
      obs_Z(n)%olat  = yo
      di = dci - int(dci)
      dip = 1.0 - di
      dj = dcj - int(dcj)
      djp = 1.0 - dj
      obs_Z(n)%io  = int(dci)
      obs_Z(n)%jo  = int(dcj)
      obs_Z(n)%a00 = dip*djp
      obs_Z(n)%a01 = dip*dj
      obs_Z(n)%a11 = di*dj
      obs_Z(n)%a10 = di*djp
      do k=1,kd
        obs_Z(n)%val(k) = val(k)
        obs_Z(n)%err(k) = err(k) * rovrf
      enddo
      obs_Z(n)%stat = 0
    enddo

    close(nru)

  else
    num_obsz = 0
    num_obs_z = 0
  endif
 
  deallocate (val)
  deallocate (err)
  deallocate (xlc)
  deallocate (ylc)

  godas_obsz_module_initialized = .true.

end function godas_obsz_init
! </FUNCTION> NAME="godas_obsz_init">


!#######################################################################
! <FUNCTION NAME="godas_obs0_init">
!
! <DESCRIPTION>
! Initialization code for surface observations, returning a pointer to
! the obs_0 array.
! </DESCRIPTION>
!
function godas_obs0_init (Grid, Domain, Time, num_obs_0, debug) &
                    result (obs_0)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  integer, intent(out)                        :: num_obs_0

  logical, intent(in), optional               :: debug

  ! return value
  type(ocean_obs0_type), dimension(:), pointer :: obs_0

  integer               :: n, nf, nsg, ntable, num_obs_files
  integer               :: i, j, k, icnt, nu, nwu, nru, nobs
  integer               :: year, month, day, hour, minute, kd, icode
  character(len=8)      :: plti
  character(len=2)      :: pltf
  integer               :: kda
  real(kind=4)          :: val, err, xo, yo
  real(kind=4)          :: xs, xe, ys, ye
  real(kind=4), allocatable, dimension(:) :: xlc, ylc
  real(kind=4)          :: dci, dcj, di, dip, dj, djp
  integer               :: pe, ierr
  character(len=256)    :: record
  type(data_type)       :: default_table, data_entry

  character(len=32)     :: obfile
  integer               :: ios

  character(len=48),  parameter :: sub_name = 'godas_obs0_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  if (godas_obs0_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' GODAS Obs0 already initialized')
  endif

  pe = mpp_pe()

  ! set local array indices
  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Dom, isg, ieg, jsg, jeg)
  nk=Grd%nk

  j=(jsd+jed)/2
  xs=Grd%xt(isd,j)
  if (xs .gt. Grd%xt(isc,j)) xs = xs - 360.0
  xe=Grd%xt(ied,j)
  if (xe .lt. Grd%xt(iec,j)) xe = xe + 360.0

  allocate (xlc(isd:ied))
  j=(jsd+jed)/2
  xlc(isd) = xs
  do i=isc,iec
    xlc(i)=Grd%xt(i,j)
  enddo
  xlc(ied) = xe
  allocate (ylc(jsd:jed))
  i=(isd+ied)/2
  do j=jsd,jed
    ylc(j)=Grd%yt(i,j)
  enddo
  xs=xlc(isd)
  xe=xlc(ied)
  ys=ylc(jsd)
  ye=ylc(jed)
  if (jed .ge. jemx) then
!   ye = ylc(jemx-4)
    if (dodebug .and. jemx-4 < jsd) then
      print *, "godas_obs.F90::"
      print *, "xs = ", xs
      print *, "xe = ", xe
      print *, "isd:ied = ", isd,ied
      print *, "xlc(isd:ied) = ", xlc(isd:ied)
      print *, "ys = ", ys
      print *, "ye = ", ye
      print *, "jemx = ", jemx
      print *, "jsd:jed = ", jsd,jed
      print *, "ylc(jsd:jed) = ", ylc(jsd:jed)
    endif
!   ye = ylc(jemx-4) !STEVE: replaced with below
    if (jemx-4 < jsd) then
      ye = ylc(jsd)
    elseif (jemx-4 > jed) then
      ye = ylc(jed)
    else
      ye = ylc(jemx-4)
    endif

  endif
! xs=xlc(isc)
! xe=xlc(iec)
! ys=ylc(jsc)
! ye=ylc(jec)
! if (jsc .eq. jsg) then
!   ys = ylc(jsc+1)
! endif
! if (jec .eq. jeg) then
!   ye = ylc(jec-3)
! endif

  nullify(obs_0)

  write( stdlog(),'(/a/)') trim(version)

! initialize data table
  default_table%gridname = 'none'
  default_table%fieldname_code = 'none'
  default_table%fieldname_file = 'none'
  default_table%file_name = 'none'
  default_table%ongrid = .FALSE.
  default_table%factor = 1.
  do n = 1,max_table
    data_table(n) = default_table
  enddo

! read observations table
  call mpp_open(nu, 'data_table', action=MPP_RDONLY)
  ntable = 0
  num_obs_files = 0
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'OBS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no OBS entry in data_table')
  endif
  num_obs_files = 0
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:3) .eq. 'sst') num_obs_files = num_obs_files + 1
    if (data_table(nf)%fieldname_code(1:3) .eq. 'sss') num_obs_files = num_obs_files + 1
  enddo

  nobs = 0
  if (num_obs_files .gt. 0) then
    write(obfile,'(a,i4.4)') 'tmp0.',pe
    nwu=1000+pe
    open(unit=nwu,file=trim(obfile),status='NEW',access='SEQUENTIAL',form='UNFORMATTED', &
         action='WRITE',iostat=ios)
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:3) .eq. 'sst') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = sst_code
        do nsg=1,nsgsobs
          read (nu) icnt
          do n=1,icnt
!STEVE: temporary
!           read (nu) plti, pltf
            read (nu) year, month, day, hour, minute
            read (nu) xo, yo, val, err
            if (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye) then
              dci = get_dec_index(xo,xlc,isd,ied)
              dcj = get_dec_index(yo,ylc,jsd,jed)
              i = int(dci)
              j = int(dcj)
              kda = min(Grd%kmt(i,j), Grd%kmt(i,j+1), Grd%kmt(i+1,j), Grd%kmt(i+1,j+1), kd)
              if (kda .gt. 0) then
                nobs = nobs + 1
!               write (nwu) plti, pltf
                write (nwu) icode, year, month, day, hour, minute
                write (nwu) xo, yo, dci, dcj, val, err
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      else if (data_table(nf)%fieldname_code(1:3) .eq. 'sss') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = sss_code
        do nsg=1,nsgsobs
          read (nu) icnt
          do n=1,icnt
!           read (nu) plti, pltf
            read (nu) year, month, day, hour, minute
            read (nu) xo, yo, val, err
            if (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye) then
              dci = get_dec_index(xo,xlc,isd,ied)
              dcj = get_dec_index(yo,ylc,jsd,jed)
              i = int(dci)
              j = int(dcj)
              kda = min(Grd%kmt(i,j), Grd%kmt(i,j+1), Grd%kmt(i+1,j), Grd%kmt(i+1,j+1), kd)
              if (kda .gt. 0) then
                nobs = nobs + 1
!               write (nwu) plti, pltf
                write (nwu) icode, year, month, day, hour, minute
                write (nwu) xo, yo, dci, dcj, val, err
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      endif
    enddo
    close(nwu)
    num_obs0 = nobs
    num_obs_0 = nobs

! allocate obs_0
    allocate( obs_0 (num_obs0) )

! fill obs_0

    nru=1000+pe
    open(unit=nru,file=trim(obfile),status='OLD',access='SEQUENTIAL',form='UNFORMATTED', &
         position='REWIND',action='READ',iostat=ios)
    do n=1,num_obs0
!     read (nru) plti, pltf
      read (nru) icode, year, month, day, hour, minute
      read (nru) xo, yo, dci, dcj, val, err
      if (icode .eq. sst_code) then
        obs_0(n)%name = 'sst'
        obs_0(n)%units = 'deg-C'
        obs_0(n)%code = sst_code
        rovrf = 1.0 / tov0f
      else if (icode .eq. sss_code) then
        obs_0(n)%name = 'sss'
        obs_0(n)%units = 'psu'
        obs_0(n)%code = sss_code
        rovrf = 1.0 / sov0f
      else
        call mpp_error(FATAL,'observation code does not match known values')
      endif
      obs_0(n)%platid = plti
      obs_0(n)%platform = pltf
      o_sec = (hour*60 + minute)*60
      o_time = set_date(year,month,day,0,0,0)
      obs_0(n)%obs_time = increment_time(o_time,o_sec,0)
      obs_0(n)%kd  = 1
      obs_0(n)%olng  = xo
      obs_0(n)%olat  = yo
      di = dci - int(dci)
      dip = 1.0 - di
      dj = dcj - int(dcj)
      djp = 1.0 - dj
      obs_0(n)%io  = int(dci)
      obs_0(n)%jo  = int(dcj)
      obs_0(n)%a00 = dip*djp
      obs_0(n)%a01 = dip*dj
      obs_0(n)%a11 = di*dj
      obs_0(n)%a10 = di*djp
      obs_0(n)%val    = val
      obs_0(n)%err    = err * rovrf
      obs_0(n)%inv    = 0.0
      obs_0(n)%inc    = 0.0
      obs_0(n)%bke    = 0.0
      obs_0(n)%stat   = 0
    enddo

    close(nru)

  else
    num_obs0 = 0
    num_obs_0 = 0
  endif

  deallocate (xlc)
  deallocate (ylc)

  godas_obs0_module_initialized = .true.

end function godas_obs0_init
! </FUNCTION> NAME="godas_obs0_init">

!#######################################################################
! <FUNCTION NAME="godas_obsa_init">
!
! <DESCRIPTION>
! Initialization code for altimetry observations, returning a pointer to
! the obs_A array.
! </DESCRIPTION>
!
function godas_obsa_init (Grid, Domain, Time, num_obs_a, debug) &
                    result (obs_A)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  integer, intent(out)                        :: num_obs_a

  logical, intent(in), optional               :: debug

  ! return value
  type(ocean_obs0_type), dimension(:), pointer :: obs_A

  integer               :: n, nf, nsg, ntable, num_obs_files
  integer               :: i, j, k, icnt, nu, nwu, nru, nobs
  integer               :: year, month, day, hour, minute, kd, icode
  character(len=8)      :: plti
  character(len=2)      :: pltf
  integer               :: kda
  real(kind=4)          :: val, err, xo, yo
  real(kind=4)          :: xs, xe, ys, ye
  real(kind=4), allocatable, dimension(:) :: xlc, ylc
  real(kind=4)          :: dci, dcj, di, dip, dj, djp
  integer               :: pe, ierr
  character(len=256)    :: record
  type(data_type)       :: default_table, data_entry

  character(len=32)     :: obfile
  integer               :: ios

  character(len=48),  parameter :: sub_name = 'godas_obsa_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  if (godas_obsa_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' GODAS ObsA already initialized')
  endif

  pe = mpp_pe()

  ! set local array indices
  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Dom, isg, ieg, jsg, jeg)
  nk=Grd%nk

  j=(jsd+jed)/2
  xs=Grd%xt(isd,j)
  if (xs .gt. Grd%xt(isc,j)) xs = xs - 360.0
  xe=Grd%xt(ied,j)
  if (xe .lt. Grd%xt(iec,j)) xe = xe + 360.0

  allocate (xlc(isd:ied))
  j=(jsd+jed)/2
  xlc(isd) = xs
  do i=isc,iec
    xlc(i)=Grd%xt(i,j)
  enddo
  xlc(ied) = xe
  allocate (ylc(jsd:jed))
  i=(isd+ied)/2
  do j=jsd,jed
    ylc(j)=Grd%yt(i,j)
  enddo
  xs=xlc(isd)
  xe=xlc(ied)
  ys=ylc(jsd)
  ye=ylc(jed)
  if (jed .ge. jemx) then
!   ye = ylc(jemx-4)
    if (dodebug .and. jemx-4 < jsd) then
      print *, "godas_obs.F90::"
      print *, "xs = ", xs
      print *, "xe = ", xe
      print *, "isd:ied = ", isd,ied
      print *, "xlc(isd:ied) = ", xlc(isd:ied)
      print *, "ys = ", ys
      print *, "ye = ", ye
      print *, "jemx = ", jemx
      print *, "jsd:jed = ", jsd,jed
      print *, "ylc(jsd:jed) = ", ylc(jsd:jed)
    endif
!   ye = ylc(jemx-4) !STEVE: replaced with below
    if (jemx-4 < jsd) then
      ye = ylc(jsd)
    elseif (jemx-4 > jed) then
      ye = ylc(jed)
    else
      ye = ylc(jemx-4)
    endif
  endif
! xs=xlc(isc)
! xe=xlc(iec)
! ys=ylc(jsc)
! ye=ylc(jec)
! if (jsc .eq. jsg) then
!   ys = ylc(jsc+1)
! endif
! if (jec .eq. jeg) then
!   ye = ylc(jec-3)
! endif

  nullify(obs_A)

  write( stdlog(),'(/a/)') trim(version)

! initialize data table
  default_table%gridname = 'none'
  default_table%fieldname_code = 'none'
  default_table%fieldname_file = 'none'
  default_table%file_name = 'none'
  default_table%ongrid = .FALSE.
  default_table%factor = 1.
  do n = 1,max_table
    data_table(n) = default_table
  enddo

! read observations table
  call mpp_open(nu, 'data_table', action=MPP_RDONLY)
  ntable = 0
  num_obs_files = 0
  do while (.true.)
    read(nu,'(a)',end=9) record
    if (record(1:1) == '#') cycle
    if (record(1:10) == '          ') cycle
    read(record,*,err=7) data_entry
    if (data_entry%gridname(1:3) .eq. 'OBS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no OBS entry in data_table')
  endif
  num_obs_files = 0
  do nf=1,ntable
    if (data_table(nf)%fieldname_code(1:4) .eq. 'altm') num_obs_files = num_obs_files + 1
  enddo

  nobs = 0
  if (num_obs_files .gt. 0) then
    write(obfile,'(a,i4.4)') 'altm.',pe
    nwu=1000+pe
    open(unit=nwu,file=trim(obfile),status='NEW',access='SEQUENTIAL',form='UNFORMATTED', &
         action='WRITE',iostat=ios)
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:4) .eq. 'altm') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = altm_code
        do nsg=1,nsgobs
          read (nu) icnt
          do n=1,icnt
!STEVE: temporary:
!           read (nu) plti, pltf
            read (nu) year, month, day, hour, minute
            read (nu) xo, yo, val, err
            if (yo .lt. no_lat_mx .and. xo .ge. xs .and. xo .lt. xe .and. yo .ge. ys .and. yo .lt. ye) then
              dci = get_dec_index(xo,xlc,isd,ied)
              dcj = get_dec_index(yo,ylc,jsd,jed)
              i = int(dci)
              j = int(dcj)
              kda = min(Grd%kmt(i,j), Grd%kmt(i,j+1), Grd%kmt(i+1,j), Grd%kmt(i+1,j+1), kd)
              if (kda .gt. 0) then
                nobs = nobs + 1
!               write (nwu) plti, pltf
                write (nwu) icode, year, month, day, hour, minute
                write (nwu) xo, yo, dci, dcj, val, err
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      endif
    enddo
    close(nwu)
    num_obsa = nobs
    num_obs_a = nobs

! allocate obs_A
    allocate( obs_A (num_obsa) )

! fill obs_A

    nru=1000+pe
    open(unit=nru,file=trim(obfile),status='OLD',access='SEQUENTIAL',form='UNFORMATTED', &
         position='REWIND',action='READ',iostat=ios)
    do n=1,num_obsa
!     read (nru) plti, pltf
      read (nru) icode, year, month, day, hour, minute
      read (nru) xo, yo, dci, dcj, val, err
      if (icode .eq. altm_code) then
        obs_A(n)%name = 'altimetry'
        obs_A(n)%units = 'meter'
        obs_A(n)%code = altm_code
      else
        call mpp_error(FATAL,'observation code does not match known values')
      endif
      obs_A(n)%platid = plti
      obs_A(n)%platform = pltf
      o_sec = (hour*60 + minute)*60
      o_time = set_date(year,month,day,0,0,0)
      obs_A(n)%obs_time = increment_time(o_time,o_sec,0)
      obs_A(n)%kd  = 1
      obs_A(n)%olng  = xo
      obs_A(n)%olat  = yo
      di = dci - int(dci)
      dip = 1.0 - di
      dj = dcj - int(dcj)
      djp = 1.0 - dj
      obs_A(n)%io  = int(dci)
      obs_A(n)%jo  = int(dcj)
      obs_A(n)%a00 = dip*djp
      obs_A(n)%a01 = dip*dj
      obs_A(n)%a11 = di*dj
      obs_A(n)%a10 = di*djp
      obs_A(n)%val    = val
      obs_A(n)%err    = err
      obs_A(n)%inv    = 0.0
      obs_A(n)%inc    = 0.0
      obs_A(n)%bke    = 0.0
      obs_A(n)%stat   = 0
    enddo

    close(nru)

  else
    num_obsa = 0
    num_obs_a = 0
  endif

  deallocate (xlc)
  deallocate (ylc)

  godas_obsa_module_initialized = .true.

end function godas_obsa_init
! </FUNCTION> NAME="godas_obsa_init">

!#######################################################################
! <SUBROUTINE NAME="godas_obs_track">
!
! <DESCRIPTION>
! Write innovations etc.
! </DESCRIPTION>
!
subroutine godas_obs_track(Time, obs_Z, obs_0, obs_A)

  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_obsz_type), intent(in)                 :: obs_Z(:)
  type(ocean_obs0_type), intent(in)                 :: obs_0(:)
  type(ocean_obs0_type), intent(in)                 :: obs_A(:)

  integer            :: n, k, nws, ios, pe
  integer            :: myr, mmo, mdy, mhr, mmi, msc
  integer            :: year, month, day, hour, minute, second
  integer            :: msec, mdys, gsec, gdys
  real               :: olng, inv, ose, oase, bse
  character(len=20)  :: ofile

  pe = mpp_pe()
  call get_date(Time%model_time, myr, mmo, mdy, mhr, mmi, msc)
  call get_time(Time%Time_step,msec,mdys)
  call get_time(gds_freq,gsec,gdys)

  obs_trk_cnt = obs_trk_cnt + 1

  if (num_obsz .gt. 0) then

    write(ofile,'(a,i4.4,a,i3.3)') 'TstatZ.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obsz
      if (obs_Z(n)%code .eq. temp_code .and. obs_Z(n)%stat .eq. 3) then
        write(nws,'(a,1x,a)') trim(obs_Z(n)%platid), trim(obs_Z(n)%platform)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_Z(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_Z(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3,a,i4)') 'Lng: ',olng,'  Lat: ',obs_Z(n)%olat,'  Lvls: ',obs_Z(n)%kd
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        do k=1,obs_Z(n)%kd
          bse = sqrt(obs_Z(n)%bke(k))
          ose = sqrt(obs_Z(n)%err(k))
          oase = sqrt(obs_Z(n)%aerr(k))
          inv = -obs_Z(n)%inv(k)
          if (obs_Z(n)%val(k) .gt. 999.9) then
            inv = 0.0
            ose = 0.0
            oase = 0.0
          else
             ose = 1.0/ose
             oase = 1.0/oase
          endif
          write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(k),obs_Z(n)%val(k),ose,oase,inv,obs_Z(n)%inc(k),bse
        enddo
      endif
    enddo
    close(nws)
    
    write(ofile,'(a,i4.4,a,i3.3)') 'SstatZ.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obsz
      if (obs_Z(n)%code .eq. salt_code .and. obs_Z(n)%stat .eq. 3) then
        write(nws,'(a,1x,a)') trim(obs_Z(n)%platid), trim(obs_Z(n)%platform)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_Z(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_Z(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3,a,i4)') 'Lng: ',olng,'  Lat: ',obs_Z(n)%olat,'  Lvls: ',obs_Z(n)%kd
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        do k=1,obs_Z(n)%kd
          bse = sqrt(obs_Z(n)%bke(k))
          ose = sqrt(obs_Z(n)%err(k))
          oase = sqrt(obs_Z(n)%aerr(k))
          inv = -obs_Z(n)%inv(k)
          if (obs_Z(n)%val(k) .gt. 999.9) then
            inv = 0.0
            ose = 0.0
            oase = 0.0
          else
             ose = 1.0/ose
             oase = 1.0/oase
          endif
          write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(k),obs_Z(n)%val(k),ose,oase,inv,obs_Z(n)%inc(k),bse
        enddo
      endif
    enddo
    close(nws)

  endif

  if ( .not.restore_sfc .and. num_obs0 .gt. 0) then

    write(ofile,'(a,i4.4,a,i3.3)') 'Tstat0.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obs0
      if (obs_0(n)%code .eq. sst_code .and. obs_0(n)%stat .eq. 3) then
        write(nws,'(a,1x,a)') trim(obs_0(n)%platid), trim(obs_0(n)%platform)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_0(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_0(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3)') 'Lng: ',olng,'  Lat: ',obs_0(n)%olat
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        bse = sqrt(obs_0(n)%bke)
        ose = sqrt(obs_0(n)%err)
        oase = sqrt(obs_0(n)%aerr)
        inv = -obs_0(n)%inv
        if (obs_0(n)%val .gt. 999.9) then
          inv = 0.0
          ose = 0.0
          oase = 0.0
        else
          ose = 1.0/ose
          oase = 1.0/oase
        endif
        write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(1),obs_0(n)%val,ose,oase,inv,obs_0(n)%inc,bse
      endif
    enddo
    close(nws)

    write(ofile,'(a,i4.4,a,i3.3)') 'Sstat0.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obs0
      if (obs_0(n)%code .eq. sss_code .and. obs_0(n)%stat .eq. 3) then
        write(nws,'(a,1x,a)') trim(obs_0(n)%platid), trim(obs_0(n)%platform)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_0(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_0(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3)') 'Lng: ',olng,'  Lat: ',obs_0(n)%olat
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        bse = sqrt(obs_0(n)%bke)
        ose = sqrt(obs_0(n)%err)
        oase = sqrt(obs_0(n)%aerr)
        inv = -obs_0(n)%inv
        if (obs_0(n)%val .gt. 999.9) then
          inv = 0.0
          ose = 0.0
          oase = 0.0
        else
           ose = 1.0/ose
           oase = 1.0/oase
        endif
        write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(1),obs_0(n)%val,ose,oase,inv,obs_0(n)%inc,bse
      endif
    enddo
    close(nws)

  endif

end subroutine godas_obs_track
! </SUBROUTINE> NAME="godas_obs_track">


!#######################################################################
! <FUNCTION NAME="get_dec_index">
!
! <DESCRIPTION>
! Finds the decimal index of p within the array pa
! </DESCRIPTION>
!
function get_dec_index (p,pa,is,ie)  result (pi)
!
  real(kind=4), intent(in)       :: p
  real(kind=4), intent(in)       :: pa(is:ie)
  integer, intent(in)            :: is, ie
! result
  real(kind=4)                   :: pi
!
  integer     :: i
!
  pi = -1.0
  if (p .le. pa(is)) then
    pi = is
  else if (p .ge. pa(ie)) then
    pi = ie
  else
    do i=is+1,ie
      if (p .ge. pa(i-1) .and. p .le. pa(i)) then
        pi = real(i-1) + (p - pa(i-1))/(pa(i) - pa(i-1))
        exit
      endif
    enddo
  endif

end function get_dec_index
! </FUNCTION> NAME="get_dec_index">

end module godas_obs_mod
