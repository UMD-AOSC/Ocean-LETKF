*
*  Script to plot a colorbar
*
*  The script will assume a colorbar is wanted even if there is 
*  not room -- it will plot on the side or the bottom if there is
*  room in either place, otherwise it will plot along the bottom and
*  overlay labels there if any.  This can be dealt with via 
*  the 'set parea' command. 
*
*	- option to draw extreme colors as triangles
*	- the colors are boxed in white
*	- input arguments in during a run execution:
* 
*	cbarm sf vert tri
*
*	sf   - scale the whole bar 1.0 = original 0.5 half the size, etc.
*	vert - 0 FORCES a horizontal bar = 1 a vertical bar
*       tri  - draw extreme colors as triangles (default is rectangles)
*  

function colorbar (args)

* parse arguments
sf = subwrd(args,1)
vert = subwrd(args,2)
tri = subwrd(args,3)

* set defaults
if (sf=''); sf=1.0; endif
if (tri=''); tri=0; endif

* Check shading information
'query shades'
shdinfo = result
if (subwrd(shdinfo,1)='None') 
  say 'Cannot plot color bar: No shading information'
  return
endif

* Get plot size info
'query gxinfo'
rec2 = sublin(result,2)
rec3 = sublin(result,3)
rec4 = sublin(result,4)
* page size
pagex = subwrd(rec2,4)
pagey = subwrd(rec2,6)
* width of plot area
xmin = subwrd(rec3,4)
xmax = subwrd(rec3,6)
xsiz = xmax-xmin
* x midpoint of plot
gxmid = xsiz/2 + xmin   
* width and midpoint of area between right edge of plot and end of page
xd = pagex - xmax
xmid = xd/2 + xmax  
* height of plot area
ymin = subwrd(rec4,4)
ymax = subwrd(rec4,6)
ysiz = ymax-ymin
* y midpoint of plot
gymid = ysiz/2 + ymin
* midpoint of area between right edge of plot and end of page
ymid = ymin/2


ylolim=0.6*sf
xdlim1=1.0*sf
xdlim2=1.5*sf  
*barsf=0.8*sf
barsf=sf
yoffset=0.2*sf
stroff=0.05*sf
strxsiz=0.09*sf
strysiz=0.11*sf
bdrthk=1  ;* border thickness

* Decide if horizontal or vertical color bar
if (ymin<ylolim & xd<xdlim1) 
  say "Not enough room in plot for a colorbar"
  return
endif

* logic for setting the bar orientation with user overides
if (ymin<ylolim | xd>xdlim1)       ;* there's room on the right side
  vchk = 1                         ;* default vertical colorbar
  if(vert = 0) ; vchk = 0 ; endif  ;* user override
else
  vchk = 0                         ;* default horizontal colorbar
  if(vert = 1) ; vchk = 1 ; endif  ;* user override
endif

* set up constants
'set strsiz 'strxsiz' 'strysiz
cnum = subwrd(shdinfo,5)
* vertical bar
if (vchk = 1 )        
  xwid = 0.25*sf
  ywid = 0.5*sf
  xl = xmid-xwid/2
  xr = xl + xwid
  if (ywid*cnum > ysiz*barsf) 
    ywid = ysiz*barsf/cnum
  endif
  yb = gymid - ywid*cnum/2
  'set string 1 l 1'
  vert = 1
else
* horizontal bar
  ywid = 0.4*sf
  xwid = 0.8*sf
  yt = ymid + yoffset
  yb = ymid
  if (xwid*cnum > xsiz*barsf)
    xwid = xsiz*barsf/cnum
  endif
  xl = gxmid - xwid*cnum/2
  'set string 1 tc 1'
  vert = 0
endif
* save the initial values of yb and xl for later use
yb0=yb
xl0=xl

* Start plotting
num = 0
maxlen=0
while (num< cnum) 
  rec = sublin(shdinfo,num+2)
  col = subwrd(rec,1)
  hi  = subwrd(rec,3)
  len = math_strlen(hi)
  if (len>maxlen); maxlen=len; endif
  if (vert) 
    yt = yb + ywid
  else 
    xr = xl + xwid
  endif
* save the corner locations for drawing a border
  if (num=0); blcorner=''xl' 'yb; endif
  if (num=cnum-1); trcorner=''xr' 'yt; endif
* save the edge locations for drawing a border with triangles 
  if (tri=1)
    if (num=1)      ;* bottom/left edge points
      if (vert=1)
        ledgex1=xl; ledgey1=yb
        redgex1=xr; redgey1=yb
      else
        tedgex1=xl; tedgey1=yt
        bedgex1=xl; bedgey1=yb
      endif
    endif
    if (num=cnum-1) ;* top/right edge points
      if (vert=1)
        ledgex2=xl; ledgey2=yb
        redgex2=xr; redgey2=yb
      else
        tedgex2=xl; tedgey2=yt
        bedgex2=xl; bedgey2=yb
      endif
    endif
  endif

* draw the filled color polygon (triangle or rectangle)
  if (num=0)
   if (tri=1)
    if (vert=1)
*     bottom triangle    
      xm = (xl+xr)*0.5
      'set line 'col
      'draw polyf 'xl' 'yt' 'xm' 'yb' 'xr' 'yt' 'xl' 'yt
      'set line 1 1 'bdrthk
      'draw line 'xl' 'yt' 'xm' 'yb
      'draw line 'xm' 'yb' 'xr' 'yt
    else
*     left triangle
      ym = (yb+yt)*0.5
      'set line 'col
      'draw polyf 'xl' 'ym' 'xr' 'yb' 'xr' 'yt' 'xl' 'ym
      'set line 1 1 'bdrthk
      'draw line 'xl' 'ym' 'xr' 'yb
      'draw line 'xr' 'yt' 'xl' 'ym
    endif
   else
    'set line 'col
    'draw recf 'xl' 'yb' 'xr' 'yt
   endif
  endif
 
  if (num>0 & num<cnum-1)
    'set line 'col
    'draw recf 'xl' 'yb' 'xr' 'yt
  endif

  if (num=cnum-1)
    if (tri=1)
      if (vert=1)
*       top triangle    
        xm = (xl+xr)*0.5
        'set line 'col
        'draw polyf 'xl' 'yb' 'xm' 'yt' 'xr' 'yb' 'xl' 'yb
        'set line 1 1 'bdrthk
        'draw line 'xl' 'yb' 'xm' 'yt
        'draw line 'xm' 'yt' 'xr' 'yb
      else
*       right triangle
        ym = (yb+yt)*0.5
        'set line 'col
        'draw polyf 'xr' 'ym' 'xl' 'yb' 'xl' 'yt' 'xr' 'ym
        'set line 1 1 5'
        'draw line 'xr' 'ym' 'xl' 'yb
        'draw line 'xl' 'yt' 'xr' 'ym
      endif
    else
      'set line 'col
      'draw recf 'xl' 'yb' 'xr' 'yt
    endif
  endif

* Reset variables for next loop execution
  if (vert) 
    yb = yt
  else
    xl = xr
  endif
  num = num + 1
endwhile


* draw the border
'set line 1 1 'bdrthk
if (tri=1)
  if (vert=1)
    'draw line 'ledgex1' 'ledgey1' 'ledgex2' 'ledgey2
    'draw line 'redgex1' 'redgey1' 'redgex2' 'redgey2
  else
    'draw line 'tedgex1' 'tedgey1' 'tedgex2' 'tedgey2
    'draw line 'bedgex1' 'bedgey1' 'bedgex2' 'bedgey2
  endif
else
  'draw rec 'blcorner' 'trcorner
endif

* now the tic marks and number strings

* reset xl and yb 
xl = subwrd(blcorner,1)
yb = subwrd(blcorner,2)

dist = 0
ticwid = 0.5
num = 0
while (num<cnum-1) 
  rec = sublin(shdinfo,num+2)
  col = subwrd(rec,1)
  hi = subwrd(rec,3)
  if (vert) 
    yt = yb + ywid
    dist = dist + ywid
  else 
    xr = xl + xwid
    dist = dist + xwid
  endif

* draw tic marks and numbers centered on boundaries of each segment of the color key
* make sure there will be room for the tic label, if not don't draw a tic or label

* length of tic marks is 20% of width of colorbar
  xdelt=(xr-xl)*0.2
  ydelt=(yt-yb)*0.2
 
  if (vert)
*   distance between vertical tics is 4 times the height of the string size
    if (dist > 4*strysiz) 
      'draw line 'xl' 'yt' 'xl+xdelt' 'yt   ;* left tic
      'draw line 'xr-xdelt' 'yt' 'xr' 'yt   ;* right tic
      xp=xr+stroff
      'draw string 'xp' 'yt' 'hi
      dist=0
    endif
  else
*   distance between horizontal tics is 1.5 times the width of the widest label
    if (dist > 1.5*maxlen*strxsiz) 
      'draw line 'xr' 'yt' 'xr' 'yt-ydelt   ;* top tic
      'draw line 'xr' 'yb' 'xr' 'yb+ydelt   ;* bottom tic
      yp=yb-(yt-yb)/2
      'draw string 'xr' 'yp' 'hi
      dist=0
    endif
  endif

* Reset variables for next loop execution
  if (vert) 
    yb = yt
  else
    xl = xr
  endif
  num = num + 1
endwhile


return
