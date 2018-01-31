#!/usr/bin/env python

#import ncodalib
from ncodalib import ncodaField2D, ncodaField3D
from coamps_grid import COAMPSGrid
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
from scipy.interpolate import griddata

#NCODA flatfile output - restart directory
rdir='/u/scrh/smedstad/GLBocn0.08/work/restart'
wdir='/home/spenny/Research/osprey/TESTCASE/expt_75.1m'

figname="test_letkf_ainc.png"

#-------------------------------------------------------------------------------
# set DTG, tau, parameter name, field type, depth index to plot
#-------------------------------------------------------------------------------
mem=1
vartype='t'
k_index=1
dtg='2014020112'

# Get figure title and filename:
title='vartype=%s mem=%s dtg=%s k=%d'%(vartype,mem,dtg,k_index)
fn=wdir+'/'+'anal'+'%3.3d'%(mem)+'.'+vartype+'.dat'

# Outfile name:
figname='test_letkf_anal.'+vartype+'.'+'mem%3.3d'%(mem)+'.'+'lvl%2.2d'%(k_index)+'.png'


nest=1
doreduce=0   # doreduce = 0, keep the input dimensions
doreduce=400 # doreduce > 0, x/y skip factors for a grid approx doreduce square
doreduce=-8  # doreduce < 0, x/y skip factors set to abs(doreduce)

# set colorbar limits
clim=(-2.,36.)
mycmap='jet'

# Initialize NCODA field type
if vartype == 'ssh':
  field=ncodaField2D()
else:
  field=ncodaField3D()

#-------------------------------------------------------------------------------
# load the map grid info
#-------------------------------------------------------------------------------
dtg1='2014020212'
glbl=True
fn1=rdir+'/datahd_pre_000000_000000_1o2000x0001_'+dtg1+'_00000000_infofld'
mygrid=COAMPSGrid('grid')
mygrid.datahd(fn1)

if mygrid.nproj < 0:
  mygrid.restartdir=rdir
  mygrid.dtg=dtg
mygrid.glbl=glbl
(grdlon,grdlat,f,hx,hy,xpos,ypos) = mygrid.grid(nest)

field.grid(mygrid,nest)

#-------------------------------------------------------------------------------
# Read data from file
#-------------------------------------------------------------------------------
field.read(fn)

#-------------------------------------------------------------------------------
# Process grid and data to reduce display output
#-------------------------------------------------------------------------------
xstep=1; ystep=1
if doreduce < 0:
  xstep = abs(doreduce)
  ystep = abs(doreduce)
else:
  if doreduce > 0:
    xstep = np.int(np.floor(field.m/doreduce))
    ystep = np.int(np.floor(field.n/doreduce))

if mygrid.glbl:
  bbox=(0.,360.,-85.,90)
else:
  bbox=mygrid.boundbox(nest)
print(bbox)

lonskip = grdlon[0::xstep,0::ystep]
latskip = grdlat[0::xstep,0::ystep]

if vartype == 'ssh':
  data=field.data[0::xstep,0::ystep]
else:
  data=field.data[0::xstep,0::ystep,k_index]

print(data.min())
print(data.max())
data[data < -900]=np.nan

#-------------------------------------------------------------------------------
# open the figure
#-------------------------------------------------------------------------------

fig = plt.figure(num=1,figsize=(8,5),dpi=120,facecolor='w',edgecolor='k')

#-------------------------------------------------------------------------------
# Interpolate the global grid to reduce the data density for the plot
#-------------------------------------------------------------------------------
if mygrid.glbl:
  ma = bm.Basemap(projection='eck4',lon_0=0.)

# reduced grid in map space
  lonskip[np.where(lonskip < -180)] = lonskip[np.where(lonskip < -180)] + 360
  lonskip[np.where(lonskip >  180)] = lonskip[np.where(lonskip >  180)] - 360
  (xi,yi) = ma(lonskip,latskip)

# reduced lat/lon grid in map space
  rm,rn = data.shape 
  lonr = -180. + np.arange(0,rm,1) * 360. / (rm-1)
  latr = -90. + np.arange(0,rn,1) * 180. / (rn-1)
  Xr,Yr = np.meshgrid(lonr,latr)
  (xx,yy) = ma(Xr,Yr)

# interpolate to xx,yy
  points = np.vstack([xi.flatten(),yi.flatten()]).transpose()
  data = griddata(points, data.flatten(), (xx,yy), method='nearest')

else:
# no interpolation, just plot the reduced grid

  if bbox[1] > 360.:
    lonskip[np.where(lonskip < -180)] = lonskip[np.where(lonskip < -180)] + 360
    lonskip[np.where(lonskip >  180)] = lonskip[np.where(lonskip >  180)] - 360
    bbox[0]=lonskip.min()
    bbox[1]=lonskip.max()

  ma = bm.Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3], \
                  llcrnrlon=bbox[0],urcrnrlon=bbox[1],lat_ts=0,resolution='i')

  (xx,yy) = ma(lonskip,latskip)

#-------------------------------------------------------------------------------
#  Setup plot specifics and save file
#-------------------------------------------------------------------------------
ax = fig.add_subplot(111)
img = ma.pcolor(xx,yy,data,cmap=mycmap,vmin=clim[0],vmax=clim[1])
ma.drawcoastlines(linewidth=.2)
ma.fillcontinents(color='white')
ma.drawmapboundary(linewidth=.2)
ma.drawparallels(np.arange(-90.,90.,30.))
ma.drawmeridians(np.arange(-180.,180.,30.))
l, b, w, h = ax.get_position().bounds
cax = plt.axes([l, b-.05, w, .04])
plt.colorbar(img,cax=cax,orientation='horizontal')
ax.set_title(title)
plt.draw()
plt.savefig(figname)
plt.close(fig)

