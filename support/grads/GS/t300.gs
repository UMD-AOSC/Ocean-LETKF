'set time 1jan1997 1jan1999'

'set lev 20'
'set lat -5 5'

* Pacific:
'set lon -230 -80'
* Indian:
*'set lon 35 120'

'eq20 = ave(t,lat=-5,lat=5)'
*'eq20 = ave(t,lat=-1,lat=1)'

'set lat 0'
'd eq20'

