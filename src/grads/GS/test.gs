function main(args)
it=128
itf=267
var=subwrd(args,1)
flip=subwrd(args,2)
if (flip=''); flip = off; endif
while (it < itf)
  it=it+1
  say 'On iteration 'it
* say 'Using var = 'var
  'set t 'it
  'set gxout shaded'
  'clear'
  'set yflip 'flip
  if (var = t)
    'scl 0 .5 .05'
  endif
  if (var = s)
    'scl 0 .08 .008'
  endif
* 'd ('var'.1 - 'var'.5)*('var'.1 - 'var'.5)'
  'd 'var'.3'
  'cbar'
  'set gxout contour'
  'set ccolor 0'
  'd 'var'.1'
* 'set ccolor 1'
* 'd 'var'.3'
* 'set ccolor 3'
* 'd 'var'.2-'var'.1'
  if (var = t)
    'set ccolor 7'
    'd s'
  endif
  if (var = s)
    'set ccolor 6'
    'd t'
  endif
  pull dummy
endwhile
