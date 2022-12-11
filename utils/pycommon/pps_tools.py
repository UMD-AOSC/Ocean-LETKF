# copied from Cheng's SPEEDY-LETKF

import numpy as np

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.util import add_cyclic_point

import matplotlib.ticker as mticker
import matplotlib.colorbar as mcolorbar



def set_cartopy_gridlines(ax        = None, \
                          crs       = ccrs.PlateCarree(), \
                          xtick     = np.arange(-150,150+60,60), \
                          #xtick     = np.arange(-180,180+60,60), \
                          ytick     = np.arange(-90,90+30,30), \
                          showxlabel = True, \
                          showylabel = True, \
                          xlabelrotation = 0, \
                          ylabelrotation = 0, \
                          linewidth = 1, \
                          linecolor = (0.5,0.5,0.5), \
                          linestyle = '--',\
                          fontsize  = 10, \
                          fontcolor ='k', \
                          fontweight='normal'):

    gl = ax.gridlines(crs=crs, draw_labels=True,linewidth=linewidth,color=linecolor,linestyle=linestyle)
    gl.top_labels = False
    gl.right_labels = False

    if linewidth > 0:
        gl.xlines = True
        gl.ylines = True
    else:
        gl.xlines = False
        gl.ylines = False

    if showxlabel:
        gl.x_inline = False
    else:
        gl.x_inline = True

    if showylabel:
        gl.y_inline = False
    else:
        gl.y_inline = True

    gl.xlocator = mticker.FixedLocator(xtick)
    gl.ylocator = mticker.FixedLocator(ytick)
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LongitudeFormatter(dateline_direction_label=True,degree_symbol='')
    gl.yformatter = LatitudeFormatter(degree_symbol='')
    gl.xlabel_style={'size':fontsize,'color':fontcolor,'weight':fontweight,'rotation':xlabelrotation}
    gl.ylabel_style={'size':fontsize,'color':fontcolor,'weight':fontweight,'rotation':ylabelrotation}

    return gl


def set_cartopy_colorbar(ax, surf, fig, \
                         location='right',\
                         pad=0.02, \
                         shrink=1, \
                         widthratio=0.1, \
                         fontsize=None, \
                         title=None):
    cax, kw = mcolorbar.make_axes(ax,\
                 location=location, \
                 pad=pad, \
                 shrink=shrink, \
                 aspect=1.0/widthratio)
    cbar = fig.colorbar(surf,cax=cax, **kw)
    if fontsize is not None:
        cbar.ax.tick_params(labelsize=fontsize)
    if title is not None:
        cbar.set_label(title)

    return cax, cbar
