MODULE mpp_sgp



!***********************************************************************
!
!     public interface from mpp_util.h
!
!***********************************************************************
  ! <INTERFACE NAME="mpp_error">
  !  <OVERVIEW>
  !    Error handler.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    It is strongly recommended that all error exits pass through
  !    <TT>mpp_error</TT> to assure the program fails cleanly. An individual
  !    PE encountering a <TT>STOP</TT> statement, for instance, can cause the
  !    program to hang. The use of the <TT>STOP</TT> statement is strongly
  !    discouraged.
  !    
  !    Calling mpp_error with no arguments produces an immediate error
  !    exit, i.e:
  !    <PRE>
  !    call mpp_error
  !    call mpp_error(FATAL)
  !    </PRE>
  !    are equivalent.
  !    
  !    The argument order
  !    <PRE>
  !    call mpp_error( routine, errormsg, errortype )
  !    </PRE>
  !    is also provided to support legacy code. In this version of the
  !    call, none of the arguments may be omitted.
  !    
  !    The behaviour of <TT>mpp_error</TT> for a <TT>WARNING</TT> can be
  !    controlled with an additional call <TT>mpp_set_warn_level</TT>.
  !    <PRE>
  !    call mpp_set_warn_level(ERROR)
  !    </PRE>
  !    causes <TT>mpp_error</TT> to treat <TT>WARNING</TT>
  !    exactly like <TT>FATAL</TT>.
  !    <PRE>
  !    call mpp_set_warn_level(WARNING)
  !    </PRE>
  !    resets to the default behaviour described above.
  !
  !    <TT>mpp_error</TT> also has an internal error state which
  !    maintains knowledge of whether a warning has been issued. This can be
  !    used at startup in a subroutine that checks if the model has been
  !    properly configured. You can generate a series of warnings using
  !    <TT>mpp_error</TT>, and then check at the end if any warnings has been
  !    issued using the function <TT>mpp_error_state()</TT>. If the value of
  !    this is <TT>WARNING</TT>, at least one warning has been issued, and
  !    the user can take appropriate action:
  !    
  !    <PRE>
  !    if( ... )call mpp_error( WARNING, '...' )
  !    if( ... )call mpp_error( WARNING, '...' )
  !    if( ... )call mpp_error( WARNING, '...' )
  !    ...
  !    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )
  !    </PRE>
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !    call mpp_error( errortype, routine, errormsg )
  !  </TEMPLATE>
  !  <IN NAME="errortype">
  !    One of <TT>NOTE</TT>, <TT>WARNING</TT> or <TT>FATAL</TT> 
  !    (these definitions are acquired by use association).
  !    <TT>NOTE</TT> writes <TT>errormsg</TT> to <TT>STDOUT</TT>. 
  !    <TT>WARNING</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>.
  !    <TT>FATAL</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>,
  !    and induces a clean error exit with a call stack traceback.
  !  </IN>
  ! </INTERFACE>
  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
     module procedure mpp_error_is
     module procedure mpp_error_rs
     module procedure mpp_error_ia
     module procedure mpp_error_ra
     module procedure mpp_error_ia_ia
     module procedure mpp_error_ia_ra
     module procedure mpp_error_ra_ia
     module procedure mpp_error_ra_ra
     module procedure mpp_error_ia_is
     module procedure mpp_error_ia_rs
     module procedure mpp_error_ra_is
     module procedure mpp_error_ra_rs
     module procedure mpp_error_is_ia
     module procedure mpp_error_is_ra
     module procedure mpp_error_rs_ia
     module procedure mpp_error_rs_ra
     module procedure mpp_error_is_is
     module procedure mpp_error_is_rs
     module procedure mpp_error_rs_is
     module procedure mpp_error_rs_rs
  end interface


CONTAINS

!#####################################################################
subroutine mpp_set_warn_level(flag)
  integer, intent(in) :: flag

  if( flag.EQ.WARNING )then
    warnings_are_fatal = .FALSE.
  else if( flag.EQ.FATAL )then
    warnings_are_fatal = .TRUE.
  else
    call mpp_error( FATAL, 'MPP_SET_WARN_LEVEL: warning flag must be set to WARNING or FATAL.' )
  end if
end subroutine mpp_set_warn_level

!#####################################################################
function mpp_error_state()
  integer :: mpp_error_state
  mpp_error_state = error_state
end function mpp_error_state

!#####################################################################
!overloads to mpp_error_basic
!support for error_mesg routine in FMS
subroutine mpp_error_mesg( routine, errormsg, errortype )
  character(len=*), intent(in) :: routine, errormsg
  integer,          intent(in) :: errortype

  call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
end subroutine mpp_error_mesg

!#####################################################################
subroutine mpp_error_noargs()
  call mpp_error(FATAL)
end subroutine mpp_error_noargs

!#####################################################################
subroutine mpp_error_Is(errortype, errormsg1, value, errormsg2)
  integer,          intent(in) :: errortype
  INTEGER,          intent(in) :: value
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/value/), errormsg2)
end subroutine mpp_error_Is
!#####################################################################
subroutine mpp_error_Rs(errortype, errormsg1, value, errormsg2)
  integer,          intent(in) :: errortype
  REAL,             intent(in) :: value
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/value/), errormsg2)
end subroutine mpp_error_Rs
!#####################################################################
subroutine mpp_error_Ia(errortype, errormsg1, array, errormsg2)
  integer,               intent(in) :: errortype
  INTEGER, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ia

!#####################################################################
subroutine mpp_error_Ra(errortype, errormsg1, array, errormsg2)
  integer,            intent(in) :: errortype
  REAL, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),   intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ra
!#####################################################################
function iarray_to_char(iarray) result(string)
integer, intent(in) :: iarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(iarray)
   write(chtmp,'(i16)') iarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function iarray_to_char
!#####################################################################
function rarray_to_char(rarray) result(string)
real, intent(in) :: rarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(rarray)
   write(chtmp,'(G16.9)') rarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function rarray_to_char



END MODULE mpp_sgp
