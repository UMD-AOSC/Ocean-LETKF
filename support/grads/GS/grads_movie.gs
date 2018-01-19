*'set mproj off'
'set display color white'
*'set vpage 1.0 10.0 1.0 7.0'
*'set parea 0.1 11.0 0.4 7.3'
'clear'
'set grads off'
'set grid off'
* =======================================================================
* SET DIMENSIONS FOR PLOT
* - default is to plot entire domain
* =======================================================================
*'set lev 0 1.75'
*'set lon 0 25'
* =======================================================================

* =======================================================================
* this next command sets the color table
* =======================================================================
'run rgbset.gs'

say 'Create gif images as well (1=yes ; 0=no)'
pull ans
frame = 1
'q file'
rec=sublin(result,5)
_endtime=subwrd(rec,12)
runscript = 1
* start at time 1
dis_t = 1

* movie loop
while(runscript)

'set t ' dis_t
'q dims'
rec=sublin(result,5)
_analysis=subwrd(rec,6)

say 'Time is ' _analysis

'clear'
*'set yflip on'
'set grads off'
'set strsiz 0.16'
'set string 1 c 2'
'set xlopts 1 2 0.16'
'set ylopts 1 2 0.16'
* =======================================================================
* SHADED PLOT MADE HERE
* - this example makes a plot of field "thp"
* - revise xaxis, yaxis, clevs to fit your needs
* - set ccols selects color from the color table
* =======================================================================
'set gxout shaded'
**'set xaxis 0 25 5'
**'set yaxis 0 1.75 .25'
*'set clevs -5.0 -4.5 -4.0 -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5'
*'set ccols   49 48 48 47 46 45 44 43 42 41 0 0 61 62 63 65 67 69'

** TEMPERATURE **
'scl 0 .5 .05'
*'scl 0 1 .1'
*'scl 0 .1 .01'
'd t.3'

** SALINITY **
*'scl 0 0.05 0.005'
*'d s.3'

** U **
*'scl 0 .1 .01'
*'d u.3'

** V **
*'scl 0 .1 .01'
*'d v.3'

** UV **
*'scl 0 .1 .01'
*'d u.3;v.3;sqrt(u.3*u.3+v.3*v.3)'

** SSH **
*'scl 0 0.04 0.004'
*'d ssh.3'

** SST **
*'scl 0 0.5 0.05'
*'d sst.3'

** SSS **
*'scl 0 0.05 0.005'
*'d sss.3'

*'set strsiz .2'
*'set string 1 l 4'
** hardwired 60 sec plotting interval
**STEVE: for horizontal plot:
*'draw string 9.4 7.4 ' _analysis
**STEVE: for vertical plot:
'draw string 9.4 8.2 ' _analysis
'cbarm'
*'cbarn 1.0 0 6.0'

* =======================================================================
* CONTOUR PLOT OVERLAY MADE HERE
* change cint, ccolor, cthick to fit your needs
* set black command attempts to suppress zero contour
* =======================================================================
*'set gxout contour'
**'set cthick 5'
**'set cint 2.0'
**'set ccolor 1'
**'set black 0 0'
**'set clab masked'
*'set clab off'
*'scl 0 .5 .05'
*'d t.3'

* =======================================================================
* FINISH
* =======================================================================
if(ans)

if( frame < 10 )
 'printim movie000'frame'.gif gif '
else 
 if ( frame < 100 )
  'printim movie00'frame'.gif gif '
 else
  if ( frame < 1000 )
    'printim movie0'frame'.gif gif '
  else
    'printim movie'frame'.gif gif '
  endif
 endif
endif

frame=frame+1
endif

**STEVE: wait here for user to ask for next timestep...
*pull dummy

if ( dis_t=_endtime )
 runscript=0
endif 

dis_t = dis_t + 1

endwhile

