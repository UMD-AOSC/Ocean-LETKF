#!/usr/bin/env python

#import ncodalib
from ncodalib import ncodaField2D, ncodaField3D
from coamps_grid import COAMPSGrid
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#NCODA flatfile output - restart directory
rdir='/u/scrh/smedstad/GLBocn0.08/work/restart/'

# load the grid info
dtg='2014020212'
glbl=True
fn=rdir+'/datahd_pre_000000_000000_1o2000x0001_'+dtg+'_00000000_infofld'
mygrid=COAMPSGrid('grid')
mygrid.datahd(fn)

# set DTG, tau, parameter name, field type, depth index to plot
param='seatmp'
sfctype='pre'
k_index=0
ftype='fcstfld'
dtg='2014020112'
tau=24
nest=1

# set colorbar limits
clim=(-2.,36.)

if sfctype == 'pre':
  field=ncodaField3D()
else:
  field=ncodaField2D()
field.grid(mygrid,nest)

if mygrid.nproj < 0:
  mygrid.restartdir=rdir
  mygrid.dtg=dtg

mygrid.glbl=glbl

(grdlon,grdlat,f,hx,hy,xpos,ypos) = mygrid.grid(nest)

# get the file
title='param=%s ftype=%s dtg=%s k=%d'%(param,ftype,dtg,k_index)

if sfctype == 'pre':
  fn=rdir+'/'+param+ \
  '_pre_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(mygrid.lev1,mygrid.lev2,field.m,field.n) \
  +dtg+'_%4.4d0000_'%(tau)+ftype
else:
  fn=rdir+'/'+param+ \
  '_sfc_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(0,0,field.m,field.n) \
  +dtg+'_%4.4d0000_'%(tau)+ftype

field.read(fn)

# make the reduced plotting grid
[grdi,grdj]=np.meshgrid(np.arange(field.m),np.arange(field.n),indexing='ij')

xstep=1; ystep=1
if field.m > 600:
  xstep = np.int(np.floor(field.m/300))
if field.n > 600:
  ystep = np.int(np.floor(field.n/300))

bbox=(1.,np.real(field.m),1.,np.real(field.n))

xskip = grdi[0::xstep,0::ystep]
yskip = grdj[0::xstep,0::ystep]

# select the level for plotting, mask -999 as nan
if sfctype == 'pre':
  data=field.data[0::xstep,0::ystep,k_index]
else:
  data=field.data[0::xstep,0::ystep]

data[data < -900]=np.nan
data=np.ma.masked_invalid(data)

# lon/lat for contouring
lons=grdlon[0::xstep,0::ystep]
lons[lons>360]=lons[lons>360]-360
lons[lons<0  ]=lons[lons<0  ]+360
lats=grdlat[0::xstep,0::ystep]

# make the figure
fig = plt.figure
img=plt.pcolor(xskip,yskip,data,cmap='jet',vmin=clim[0],vmax=clim[1])
plt.contour(xskip,yskip,lons,np.arange(0.,370.,30.),colors='w',linewidths=1)
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.contour(xskip,yskip,lats,np.arange(-90.,90.,30.),colors='w',linewidths=1)
plt.colorbar(img,orientation='horizontal')
plt.title(title)
figname="ncoda.test.png"
plt.draw()
plt.savefig(figname)
plt.close()

# # make the figure
# fig = plt.figure(num=1,figsize=(8,5),dpi=120,facecolor='w',edgecolor='k')
# ax = fig.add_subplot(111)
# img=plt.pcolor(xskip,yskip,data,cmap='jet')
# plt.contour(xskip,yskip,lons,np.arange(0.,370.,30.),colors='w',linewidths=1)
# plt.rcParams['contour.negative_linestyle'] = 'solid'
# plt.contour(xskip,yskip,lats,np.arange(-90.,90.,30.),colors='w',linewidths=1)
# l, b, w, h = ax.get_position().bounds
# cax = plt.axes([l, b-.1, w, .04])
# plt.colorbar(img,cax=cax,orientation='horizontal')
# ax.set_title(title)
# figname="ncoda.test.png"
# plt.draw()
# plt.savefig(figname)
# plt.close(fig)
#
