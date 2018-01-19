#!/usr/bin/env python

import ncodalib
from ncodalib import ncodaField3D
from ncodalib import ncodaField2D
from coamps_grid import COAMPSGrid
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm

# set restart dir path and dtg here

rdir='/u/rpsea/relo/RPSEA3X/scr/gom_rec/analysis/restart'
dtg='2015042000'

# set parameter name and field type and depth index here
param='seatmp'
sfctype='pre'
k_index=0
ftype='analinc'
# ftype='fcsterr'
# sfctype='sfc'

# set colorbar limits here
clim=(-2.,2.)

# that should be all you need to change

# get grid info

fn=rdir+'/datahd_pre_000000_000000_1o2000x0001_'+dtg+'_00000000_infofld'
mygrid=COAMPSGrid('gom')
mygrid.datahd(fn)

if sfctype == 'pre':
  field=ncodaField3D()
else:
  field=ncodaField2D()
field.grid(mygrid,1)
(grdlon,grdlat,f,hx,hy,xpos,ypos) = mygrid.grid(0)

# get file

title='param=%s ftype=%s dtg=%s k=%d'%(param,ftype,dtg,k_index)

if sfctype == 'pre':
  fn=rdir+'/'+param+ \
  '_pre_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(mygrid.lev1,mygrid.lev2,field.m,field.n) \
  +dtg+'_00000000_'+ftype
else:
  fn=rdir+'/'+param+ \
  '_sfc_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(0,0,field.m,field.n) \
  +dtg+'_00000000_'+ftype

field.read(fn)

# reduce field size for plotting

xstep=1; ystep=1
if field.m > 600:
  xstep = np.int(np.floor(field.m/300))
if field.n > 600:
  ystep = np.int(np.floor(field.n/300))

bbox=mygrid.boundbox()

lonskip = grdlon[0::xstep,0::ystep]
latskip = grdlat[0::xstep,0::ystep]

if sfctype == 'pre':
  data=field.data[0::xstep,0::ystep,k_index]
else:
  data=field.data[0::xstep,0::ystep]

data[data < -900]=np.nan

# open the figure

fig = plt.figure(num=1,figsize=(8,5),dpi=120,facecolor='w',edgecolor='k')

# mercator projection using the grid bounding box

ma = bm.Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3], \
                llcrnrlon=bbox[0],urcrnrlon=bbox[1],lat_ts=0,resolution='i')
(xx,yy)=ma(lonskip,latskip)

ax = fig.add_subplot(111)
# color map min,max set here
img = ma.pcolor(xx,yy,data,vmin=clim[0],vmax=clim[1])
ma.drawcoastlines(linewidth=.2)
ma.fillcontinents(color='white')
ma.drawmapboundary(linewidth=.2)
ma.drawparallels(np.arange(10.,30.,10.))
ma.drawmeridians(np.arange(260.,300.,10.))
l, b, w, h = ax.get_position().bounds
cax = plt.axes([l, b-.05, w, .04])
plt.colorbar(img,cax=cax,orientation='horizontal')
ax.set_title(title)
figname="ncoda.test.png"
plt.draw()
plt.savefig(figname)
plt.close(fig)



