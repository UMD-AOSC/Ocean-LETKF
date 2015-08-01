module godas_obs_bias_mod
!
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!</CONTACT>
!
!<OVERVIEW>
! This module manages climate observations for bias control in GODAS.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module manages climate observations for bias control in GODAS.
! Initialization should be called after GODAS and GODAS BIAS initialization.
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
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,            only: ocean_time_type

use godas_types_mod,            only: ocean_obsz_type
use godas_data_mod,             only: apply_bias_correction
use godas_data_mod,             only: kass, nsgobs, nsgsobs, jemx
use godas_data_mod,             only: temp_code, xbt_code, tao_code, argo_code
use godas_data_mod,             only: salt_code, stao_code, sargo_code
use godas_data_mod,             only: sst_code, sss_code
use godas_data_mod,             only: altm_code, tp_code, j1_code
use godas_data_mod,             only: num_obsbiasz
use godas_data_mod,             only: gds_freq, tovrf, sovrf, tov0f, sov0f
use godas_data_mod,             only: obs_trk_cnt
use godas_data_mod,             only: no_lat_mx
!
implicit none

private

logical :: godas_obs_bias_module_initialized = .false.

character(len=256) :: version = '$Id: godas_obs_bias.F90,v 1.0 2007/03/20 10:08:00 gtn Exp $'
character(len=256) :: tagname = 'Tag $Name: gds2p0d $'
character(len=48), parameter          :: mod_name = 'godas_obs_bias_mod'

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

public  godas_obsz_bias_init
public  godas_obs_bias_track
contains


!#######################################################################
! <FUNCTION NAME="godas_obsz_bias_init">
!
! <DESCRIPTION>
! Initialization code for profile observations for bias correction, returning a 
! pointer to the obs_bias_Z array.
! </DESCRIPTION>
!
function godas_obsz_bias_init (Grid, Domain, Time, num_obs_bias_z, debug) &
                    result (obs_bias_Z)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  integer, intent(out)                        :: num_obs_bias_z

  logical, intent(in), optional               :: debug

  ! return value
  type(ocean_obsz_type), dimension(:), pointer :: obs_bias_Z

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

  character(len=48),  parameter :: sub_name = 'godas_obsz_bias_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  if (.not. apply_bias_correction) then
    return
  endif

  if (godas_obs_bias_module_initialized) then
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
    ye = ylc(jemx-4)
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

  nullify(obs_bias_Z)

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
    if (data_entry%gridname(1:3) .eq. 'BAS') then
      ntable=ntable+1
      data_table(ntable) = data_entry
    endif
  enddo
7 call mpp_error(FATAL,'error reading data_table')
9 continue
  call mpp_close(nu)
  if (ntable .eq. 0) then
    call mpp_error(FATAL,'no BAS entry in data_table')
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
!   call mpp_open(nwu,'tbiasZ',action=MPP_WRONLY,form=MPP_IEEE32,access=MPP_SEQUENTIAL, &
!                       threading=MPP_MULTI,fileset=MPP_MULTI)
!
    write(obfile,'(a,i4.4)') 'tbiasZ.',pe
    nwu=1000+pe
    open(unit=nwu,file=trim(obfile),status='NEW',access='SEQUENTIAL',form='UNFORMATTED', &
         action='WRITE',iostat=ios)
    do nf=1,ntable
      if (data_table(nf)%fieldname_code(1:4) .eq. 'temp') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = temp_code
        do nsg=1,nsgobs
          read (nu) icnt
          do n=1,icnt
            val = 0.0
            err = 0.0
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
                write (nwu) icode, year, month, day, hour, minute, kda
                write (nwu) xo, yo, dci, dcj
                write (nwu) (val(k), err(k), k=1,kda)
              endif
            endif
          enddo
        enddo
        call mpp_close(nu)
      else if (data_table(nf)%fieldname_code(1:4) .eq. 'salt') then
        call mpp_open(nu,trim(data_table(nf)%file_name),action=MPP_RDONLY,form=MPP_IEEE32, &
                        access=MPP_SEQUENTIAL,threading=MPP_MULTI,fileset=MPP_SINGLE)
        icode = salt_code
        do nsg=1,nsgobs
          read (nu) icnt
          do n=1,icnt
            val = 0.0
            err = 0.0
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
    num_obsbiasz = nobs
    num_obs_bias_z = nobs

! allocate obs_bias_Z
    allocate( obs_bias_Z (num_obsbiasz) )

    do n=1,num_obsbiasz
#ifndef STATIC_MEMORY
      allocate( obs_bias_Z(n)%val(kass) )
      allocate( obs_bias_Z(n)%err(kass) )
      allocate( obs_bias_Z(n)%aerr(kass) )
      allocate( obs_bias_Z(n)%inv(kass) )
      allocate( obs_bias_Z(n)%inc(kass) )
      allocate( obs_bias_Z(n)%bke(kass) )
#endif
      obs_bias_Z(n)%val(:) = 0.0
      obs_bias_Z(n)%err(:) = 0.0
      obs_bias_Z(n)%aerr(:) = 0.0
      obs_bias_Z(n)%inv(:) = 0.0
      obs_bias_Z(n)%inc(:) = 0.0
      obs_bias_Z(n)%bke(:) = 0.0
    enddo

! fill obs_bias_Z data array
!
!  could not read files with mpp_open.  see above.
!
!   call mpp_open(nwu,'tbiasZ',action=MPP_RDONLY,form=MPP_IEEE32,access=MPP_SEQUENTIAL, &
!                       threading=MPP_MULTI,fileset=MPP_MULTI)
!
    nru=1000+pe
    open(unit=nru,file=trim(obfile),status='OLD',access='SEQUENTIAL',form='UNFORMATTED', &
         position='REWIND',action='READ',iostat=ios)
    do n=1,num_obsbiasz
      read (nru) plti, pltf
      read (nru) icode, year, month, day, hour, minute, kd
      read (nru) xo, yo, dci, dcj
      read (nru) (val(k), err(k), k=1,kd)
      if (icode .eq. temp_code) then
        obs_bias_Z(n)%name = 'temp'
        obs_bias_Z(n)%units = 'deg-C'
        obs_bias_Z(n)%code = temp_code
        rovrf = 1.0 / tovrf
      else if (icode .eq. salt_code) then
        obs_bias_Z(n)%name = 'salt'
        obs_bias_Z(n)%units = 'psu'
        obs_bias_Z(n)%code = salt_code
        rovrf = 1.0 / sovrf
      else
        call mpp_error(FATAL,'observation code does not match known values')
      endif
      obs_bias_Z(n)%platid = plti
      obs_bias_Z(n)%platform = pltf
      o_sec = (hour*60 + minute)*60
      o_time = set_date(year,month,day,0,0,0)
      obs_bias_Z(n)%obs_time = increment_time(o_time,o_sec,0)
      obs_bias_Z(n)%kd  = kd
      obs_bias_Z(n)%olng  = xo
      obs_bias_Z(n)%olat  = yo
      di = dci - int(dci)
      dip = 1.0 - di
      dj = dcj - int(dcj)
      djp = 1.0 - dj
      obs_bias_Z(n)%io  = int(dci)
      obs_bias_Z(n)%jo  = int(dcj)
      obs_bias_Z(n)%a00 = dip*djp
      obs_bias_Z(n)%a01 = dip*dj
      obs_bias_Z(n)%a11 = di*dj
      obs_bias_Z(n)%a10 = di*djp
      do k=1,kd
        obs_bias_Z(n)%val(k) = val(k)
        obs_bias_Z(n)%err(k) = err(k) * rovrf
      enddo
      obs_bias_Z(n)%stat = 0
    enddo

    close(nru)

  else
    num_obsbiasz = 0
    num_obs_bias_z = 0
  endif
 
  deallocate (val)
  deallocate (err)
  deallocate (xlc)
  deallocate (ylc)

  godas_obs_bias_module_initialized = .true.

end function godas_obsz_bias_init
! </FUNCTION> NAME="godas_obsz_bias_init">


!#######################################################################
!#######################################################################
! <SUBROUTINE NAME="godas_obs_bias_track">
!
! <DESCRIPTION>
! Write innovations etc.
! </DESCRIPTION>
!
subroutine godas_obs_bias_track(Time, obs_bias_Z)

  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_obsz_type), intent(in)                 :: obs_bias_Z(:)

  integer            :: n, k, nws, ios, pe
  integer            :: myr, mmo, mdy, mhr, mmi, msc
  integer            :: year, month, day, hour, minute, second
  integer            :: msec, mdys, gsec, gdys
  real               :: olng, inv, ose, oase, bse
  character(len=20)  :: ofile

  if (.not. apply_bias_correction) then
    return
  endif

  pe = mpp_pe()
  call get_date(Time%model_time, myr, mmo, mdy, mhr, mmi, msc)
  call get_time(Time%Time_step,msec,mdys)
  call get_time(gds_freq,gsec,gdys)

  obs_trk_cnt = obs_trk_cnt + 1

  if (num_obsbiasz .gt. 0) then

    write(ofile,'(a,i4.4,a,i3.3)') 'TBstatZ.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obsbiasz
      if (obs_bias_Z(n)%code .eq. temp_code .and. obs_bias_Z(n)%stat .eq. 3) then
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_bias_Z(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_bias_Z(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3,a,i4)') 'Lng: ',olng,'  Lat: ',obs_bias_Z(n)%olat,'  Lvls: ',obs_bias_Z(n)%kd
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        do k=1,obs_bias_Z(n)%kd
          bse = sqrt(obs_bias_Z(n)%bke(k))
          ose = sqrt(obs_bias_Z(n)%err(k))
          oase = sqrt(obs_bias_Z(n)%aerr(k))
          inv = -obs_bias_Z(n)%inv(k)
          if (obs_bias_Z(n)%val(k) .gt. 999.9) then
            inv = 0.0
            ose = 0.0
            oase = 0.0
          else
             ose = 1.0/ose
             oase = 1.0/oase
          endif
          write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(k),obs_bias_Z(n)%val(k),ose,oase,inv,obs_bias_Z(n)%inc(k),bse
        enddo
      endif
    enddo
    close(nws)
    
    write(ofile,'(a,i4.4,a,i3.3)') 'SBstatZ.',pe,'.',obs_trk_cnt
    nws=1000+pe
    open(unit=nws,file=trim(ofile),status='NEW',access='SEQUENTIAL',form='FORMATTED', &
         action='WRITE',iostat=ios)
    do n=1,num_obsbiasz
      if (obs_bias_Z(n)%code .eq. salt_code .and. obs_bias_Z(n)%stat .eq. 3) then
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2),a,i5,a,i6)')  &
          'Mdl Date: ',myr,'/',mmo,'/',mdy,'  Time: ',mhr,':',mmi,':',msc,'  TS:',msec,'  AI:',gsec
        call get_date(obs_bias_Z(n)%obs_time, year, month, day, hour, minute, second)
        write(nws,'(a,i4,2(a,i2),a,i2,2(a,i2))')  &
          'Obs Date: ',year,'/',month,'/',day,'  Time: ',hour,':',minute,':',second
        olng = obs_bias_Z(n)%olng
        if (olng .le. -180.0) olng = olng + 360.0
        if (olng .gt. 180.0)  olng = olng - 360.0
        write(nws,'(a,f9.3,a,f9.3,a,i4)') 'Lng: ',olng,'  Lat: ',obs_bias_Z(n)%olat,'  Lvls: ',obs_bias_Z(n)%kd
        write(nws,'(7x,a,6x,a,5(9x,a))') 'Z','Obs','oEr','aEr','Inv','Inc','bEr'
        do k=1,obs_bias_Z(n)%kd
          bse = sqrt(obs_bias_Z(n)%bke(k))
          ose = sqrt(obs_bias_Z(n)%err(k))
          oase = sqrt(obs_bias_Z(n)%aerr(k))
          inv = -obs_bias_Z(n)%inv(k)
          if (obs_bias_Z(n)%val(k) .gt. 999.9) then
            inv = 0.0
            ose = 0.0
            oase = 0.0
          else
             ose = 1.0/ose
             oase = 1.0/oase
          endif
          write(nws,'(f9.2,f11.5,1p5e12.3)') Grd%zt(k),obs_bias_Z(n)%val(k),ose,oase,inv,obs_bias_Z(n)%inc(k),bse
        enddo
      endif
    enddo
    close(nws)

  endif

end subroutine godas_obs_bias_track
! </SUBROUTINE> NAME="godas_obs_bias_track">


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

end module godas_obs_bias_mod
