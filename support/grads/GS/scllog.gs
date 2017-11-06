function scllog(arg)


lb = subwrd(arg,1)
ub = subwrd(arg,2)
inc = subwrd(arg,3)

clstring = ''

base=2
idx=0
i=lb
while (i <= ub)
   val=math_pow(base,idx)
   val=val*i
   say val
   clstring = clstring%' '%val%' '
   i=i+inc
   idx=idx+1
endwhile

'set clevs '%clstring
