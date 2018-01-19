function scl(arg)


lb = subwrd(arg,1)
ub = subwrd(arg,2)
inc = subwrd(arg,3)

clstring = ''

i=lb
while (i <= ub)
   clstring = clstring%' '%i%' '
   i=i+inc
endwhile

'set clevs '%clstring
