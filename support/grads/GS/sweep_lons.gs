function sweep_lons(arg)
'reset'
'set lev 0 1000'
'set gxout shaded'
'set grads off'
'set t 460'

input1=subwrd(arg,1)
vin=input1
input2=subwrd(arg,2)
input3=subwrd(arg,3)
input4=subwrd(arg,4)
min=input2
max=input3
inc=input4

say 'plotting: "'vin'"'
count = 1 
while (count <= 720) 
  'set x 'count
  'c'
  'set yflip on'
  'scl 'min' 'max' 'inc
* 'd 'vin
  'set gxout shaded'
  'd 'vin'.5'
  'set gxout contour'
  'd 'vin'.1'
  'cbar'
  say count
* 'q dim'
* pull pause
  count = count + 1 
endwhile

return
