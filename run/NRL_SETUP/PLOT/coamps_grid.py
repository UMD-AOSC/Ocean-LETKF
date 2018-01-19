#!/usr/bin/env python

import numpy as np
import os

class Grid:
    """Base class of all grid types.  This contains the basic parameters that
    are common to all grid types.
    """

    def __init__(self, name):
        self.name = name

class COAMPSGrid(Grid):
    """Class for COAMPS grids.  It inherits the Grid
    class, and includes parameters specific to this type.
    """

    def __init__(self, name):
        Grid.__init__(self, name)
        self.type = 'COAMPS Grid'
        self.mx_grids = 7
        self.nproj = 2
        self.reflat = 42.5
        self.reflon = 16.5
        self.iref = np.ones(self.mx_grids,dtype=int)
        self.jref = np.ones(self.mx_grids,dtype=int)
        self.ii = np.ones(self.mx_grids,dtype=int)
        self.jj = np.ones(self.mx_grids,dtype=int)
        self.stdlt1 = 60.
        self.stdlt2 = 30.
        self.stdlon = 240.
        self.delx = 81000. * np.ones(self.mx_grids,dtype=float)
        self.dely = 81000. * np.ones(self.mx_grids,dtype=float)
        self.m = 61 * np.ones(self.mx_grids,dtype=int)
        self.n = 61 * np.ones(self.mx_grids,dtype=int)
        self.nnest = 1
        self.npgrid = np.array([1] + list(range(1,self.mx_grids)))
        self.rnest = 3.
        self.lnmove = [ False ] * self.mx_grids
        self.lncom = False
        self.kka = 30
        self.kko = 1
        self.kkom = 1
        self.kkosm = 0
        self.restartdir = ''
        self.dtg = ''
        self.glbl = False

    def setdefault(self):
        self.mx_grids = 7
        self.nproj = 0
        self.reflat = 0.
        self.reflon = 0.
        self.iref = np.zeros(self.mx_grids,dtype=int)
        self.jref = np.zeros(self.mx_grids,dtype=int)
        self.ii = np.zeros(self.mx_grids,dtype=int)
        self.jj = np.zeros(self.mx_grids,dtype=int)
        self.stdlt1 = 0.
        self.stdlt2 = 0.
        self.stdlon = 0.
        self.delx = np.zeros(self.mx_grids,dtype=float)
        self.dely = np.zeros(self.mx_grids,dtype=float)
        self.m = np.zeros(self.mx_grids,dtype=int)
        self.n = np.zeros(self.mx_grids,dtype=int)
        self.nnest = 1
        self.npgrid = np.ones(self.mx_grids,dtype=int)
        self.rnest = 3.
        self.restartdir = ''
        self.dtg = ''
        self.glbl = False

    def datahd(self,datahdfile):
      import os
      from struct import unpack

      if not os.path.isfile(datahdfile):
        print('missing or bad file ' + datahdfile)

#     fid=file(datahdfile,'rb');
      fid=open(datahdfile,'rb');
      datao=np.float32(unpack('>2000f',fid.read(2000*4)))
      fid.close()

      self.kko     = np.int(np.floor(datao[ 0]))
      self.kkom    = np.int(np.floor(datao[ 0]))
      self.kkosm   = np.int(np.floor(datao[ 1]))
      self.nproj   = np.int(np.floor(datao[ 2]))
      self.stdlt1  =          datao[ 3]
      self.stdlt2  =          datao[ 4]
      self.stdlon  =          datao[ 5]
      self.reflat  =          datao[ 6]
      self.reflon  =          datao[ 7]
      self.nnest   = np.int(np.floor(datao[10]))
      self.lev1    =          datao[20]
      self.lev2    =          datao[21]

      for i in np.arange(self.nnest):
        k = 30 + i * 30 - 1;
        self.m[i]      = datao[k+0];
        self.n[i]      = datao[k+1];
        self.ii[i]     = datao[k+2];
        self.jj[i]     = datao[k+3];
        self.iref[i]   = datao[k+4];
        self.jref[i]   = datao[k+5];
        self.npgrid[i] = datao[k+6];
        self.delx[i]   = datao[k+7];
        self.dely[i]   = datao[k+8];

    def define(self,nproj,reflat,reflon,iref,jref,ii,jj,
               stdlt1,stdlt2,stdlon,delx,dely,m,n,nnest,npgrid,rnest):
        self.nproj = nproj
        self.reflat = reflat
        self.reflon = reflon
        iref=np.array((iref,),dtype=int)
        self.iref[:iref.size]=iref
        self.iref[iref.size:]=0.
        jref=np.array((jref,),dtype=int)
        self.jref[:jref.size]=jref
        self.jref[jref.size:]=0.
        ii=np.array((ii,),dtype=int)
        self.ii[:ii.size]=ii
        self.ii[ii.size:]=0.
        jj=np.array((jj,),dtype=int)
        self.jj[:jj.size]=jj
        self.jj[jj.size:]=0.
        self.stdlt1 = stdlt1
        self.stdlt2 = stdlt2
        self.stdlon = stdlon
        delx=np.array((delx,),dtype=float)
        self.delx[:delx.size]=delx
        self.delx[delx.size:]=0.
        dely=np.array((dely,),dtype=float)
        self.dely[:dely.size]=dely
        self.dely[dely.size:]=0.
        m=np.array((m,),dtype=int)
        self.m[:m.size]=m
        self.m[m.size:]=0.
        n=np.array((n,),dtype=int)
        self.n[:n.size]=n
        self.n[n.size:]=0.
        self.nnest = nnest
        npgrid=np.array((npgrid,),dtype=int)
        self.npgrid[:npgrid.size]=npgrid
        self.npgrid[npgrid.size:]=0.
        self.rnest = rnest

    def definenml(self, nml):

        delx     = nml['delx']  
        dely     = nml['dely']  
        ii       = nml['ii']    
        iref     = nml['iref']  
        jj       = nml['jj']    
        jref     = nml['jref']  
        m        = nml['m']     
        n        = nml['n']     
        npgrid   = nml['npgrid']

        self.kka      = nml['kka']   
        self.kko      = nml['kko']   
        self.kkom     = nml['kkom']  
        self.kkosm    = nml['kkosm']
        self.lncom    = nml['lncom']
        self.lnmove   = nml['lnmove']
        self.nnest    = nml['nnest']
        self.nproj    = nml['nproj']
        self.stdlt1   = nml['phnt1']
        self.stdlt2   = nml['phnt2']
        self.reflat   = nml['rlat']  
        self.reflon   = nml['rlon']  
        self.rnest    = nml['rnest']
        self.nproj    = nml['nproj']
        self.stdlon   = nml['alnnt']

        iref=np.array((iref,),dtype=int)
        self.iref[:iref.size]=iref
        self.iref[iref.size:]=0.
        jref=np.array((jref,),dtype=int)
        self.jref[:jref.size]=jref
        self.jref[jref.size:]=0.
        ii=np.array((ii,),dtype=int)
        self.ii[:ii.size]=ii
        self.ii[ii.size:]=0.
        jj=np.array((jj,),dtype=int)
        self.jj[:jj.size]=jj
        self.jj[jj.size:]=0.
        delx=np.array((delx,),dtype=float)
        self.delx[:delx.size]=delx
        self.delx[delx.size:]=0.
        dely=np.array((dely,),dtype=float)
        self.dely[:dely.size]=dely
        self.dely[dely.size:]=0.
        m=np.array((m,),dtype=int)
        self.m[:m.size]=m
        self.m[m.size:]=0.
        n=np.array((n,),dtype=int)
        self.n[:n.size]=n
        self.n[n.size:]=0.
        npgrid=np.array((npgrid,),dtype=int)
        self.npgrid[:npgrid.size]=npgrid
        self.npgrid[npgrid.size:]=0.

    def grid(self,nest):

#     Name        Type      Usage            Description
#   ---------   --------   -------   -------------------------------
#    delx       real        input    grid spacing in x-direction
#    dely       real        input    grid spacing in y-direction
#    f          real        output   coriolis force at point
#    grdlat     real        output   latitude of points
#    grdlon     real        output   longitude of points
#    hx         real        output   map factor in x-direction
#    hy         real        output   map factor in y-direction
#    igrid      integer     input    type of grid projection:
#                                    1, mercator projection
#                                    2, lambert conformal projection
#                                    3, polar stereographic projection
#                                    4, cartesian coordinates
#                                    5, spherical projection
#    iref       integer     input    i-coordinate of reference point
#    jref       integer     input    j-coordinate of reference point
#    m          integer     input    number of points in x-direction
#    n          integer     input    number of points in y-direction
#    reflat     real        input    latitude at reference point 
#    reflon     real        input    longitude at reference point
#    stdlt1     real        input    standard latitude of grid
#    stdlt2     real        input    second standard latitude of grid 
#                                    (required for lambert conformal)
#    stdlon     real        input    standard longitude of grid
#                                    (points to the north)
#    xpos       real        output   x-position in meters from line 
#                                    extending from pole through 
#                                    standard longitue (reflat)
#    ypos       real        output   y-position in meters from the pole

# grid parameters for this nest
      igrid=self.nproj
      reflat=self.reflat
      reflon=self.reflon
      iref=(self.iref[nest-1]).astype(float)
      jref=(self.jref[nest-1]).astype(float)
      stdlt1=self.stdlt1
      stdlt2=self.stdlt2
      stdlon=self.stdlon
      delx=self.delx[nest-1]
      dely=self.dely[nest-1]
      m=(self.m[nest-1]).astype(int)
      n=(self.n[nest-1]).astype(int)

#     ..local constants
      flat = 0.
      pi2 = np.pi / 2.
      pi4 = np.pi / 4.
      d2r = np.pi / 180.
      r2d = 180. / np.pi
      radius = 6371229.
      omega4 = 4. * np.pi / 86400.
      onedeg = radius * 2. * np.pi / 360.

#     ..set computation points
      [grdi,grdj]=np.meshgrid(np.arange(m),np.arange(n),indexing='ij')

#     ..compute distances on grid
      distx = (grdi - iref) * delx
      disty = (grdj - jref) * dely

      if igrid == 1:

#        ..mercator projection
        gcon = 0.
        deg = np.abs(stdlt1) * d2r
        con1 = np.cos(deg)
        con2 = radius * con1
        deg = reflat * 0.5 * d2r
        rih = con2 * np.log (np.tan (pi4 + deg))
        xpos = np.ones((m,n))
        ypos = np.ones((m,n))
        rr = rih + (grdj - jref) * dely
        grdlat = (2. * np.arctan (np.exp (rr / con2)) - pi2) * r2d
        grdlon = reflon + (grdi - iref) * r2d * delx / con2
        grdlon[grdlon>360]=grdlon[grdlon>360]-360
        grdlon[grdlon<0  ]=grdlon[grdlon<0  ]+360
        deg = grdlat * d2r
        f = omega4 * np.sin (deg)
        hx = np.cos (deg) / con1
        hy = hx

      elif (igrid == 2 or igrid == 3):

#        ..lambert conformal or polar stereographic
        if igrid == 2:
          if stdlt1 == stdlt2:
            gcon = np.sin (np.abs (stdlt1) * d2r)
          else:
            gcon = (
                     (np.log (np.sin ((90. - np.abs (stdlt1)) * d2r))
                     -np.log (np.sin ((90. - np.abs (stdlt2)) * d2r))) /
                     (np.log (np.tan ((90. - np.abs (stdlt1)) * .5 * d2r))
                     -np.log (np.tan ((90. - np.abs (stdlt2)) * .5 * d2r)))
                   )
        else:
          gcon = 1.
        ogcon = 1. / gcon
        ihem = np.round (np.abs (stdlt1) / stdlt1)
        deg = (90. - np.abs (stdlt1)) * d2r
        cn1 = np.sin (deg)
        cn2 = radius * cn1 * ogcon
        deg = deg * .5
        cn3 = np.tan (deg)
        deg = (90. - np.abs (reflat)) * .5 * d2r
        cn4 = np.tan (deg)
        rih = cn2 * (cn4 / cn3)**gcon
        deg = (reflon - stdlon) * d2r * gcon
        xih =  rih * np.sin (deg)
        yih = -rih * np.cos (deg) * ihem
        x = xih + distx
        y = yih + disty
        xpos = x * ihem
        ypos = y * ihem
        rr = np.sqrt (x*x + y*y)
        grdlat = r2d * (pi2 - 2. * np.arctan (cn3 * (rr / cn2)**ogcon)) * ihem
        yy = -y * ihem
        angle = np.arctan2 (x,yy) * r2d
        angle[(yy==0) & (x == 0)]=-90
        angle[(yy==0) & (x  > 0)]= 90
        grdlon = stdlon + angle * ogcon
        deg = grdlat * d2r
        f = omega4 * np.sin (deg)
        grdlon[grdlon>360]=grdlon[grdlon>360]-360
        grdlon[grdlon<0  ]=grdlon[grdlon<0  ]+360
        deg = (90. - grdlat * ihem) * d2r
        cn5 = np.sin (deg)
        deg = deg * .5
        cn6 = np.tan (deg)
        if igrid == 2:
           hx = cn5 / cn1 * (cn6 / cn3)**(-gcon)
        else:
           hx = ( (1. + np.sin (np.abs (grdlat) * d2r)) /
                  (1. + np.sin (np.abs (stdlt1) * d2r)) )
        hy = hx

      elif igrid == 4:

#        ..analytic grid
        gcon = 2.
        cn2 = delx / onedeg
        xpos = np.nan*np.ones(m,n)
        ypos = np.nan*np.ones(m,n)
        grdlat = (grdj - jref) * dely / onedeg + reflat
        grdlon = (grdi - iref) * delx / onedeg + reflon
        grdlon[grdlon>360]=grdlon[grdlon>360]-360
        grdlon[grdlon<0  ]=grdlon[grdlon<0  ]+360
        hx = np.ones(m,n)
        hy = np.ones(m,n)
        f = sin (flat * d2r) * omega4

      elif igrid == 5:

#        ..spherical grid
        xpos = np.nan*np.ones(m,n)
        ypos = np.nan*np.ones(m,n)
        grdlat = (grdj - jref) * dely / onedeg + reflat
        grdlon = (grdi - iref) * delx / onedeg + reflon
        grdlon[grdlon>360]=grdlon[grdlon>360]-360
        grdlon[grdlon<0  ]=grdlon[grdlon<0  ]+360
        f = np.sin (grdlat * d2r) * omega4
        hx = np.cos (grdlat * d2r)
        hx[hx<0]=0
        hy = np.ones(m,n)

      elif igrid == 6:

#        ..gnomic
        cn1=np.sin(stdlt1*d2r)
        cn2=np.cos(stdlt1*d2r)
        cn3=np.sin(stdlt2*d2r)
        cn4=np.cos(stdlt2*d2r)
#       ..find x,y of reference lat,lon in xih,yih
        cn5=np.sin((reflon-stdlon)*d2r)
        cn6=np.cos((reflon-stdlon)*d2r)
        cn7=cn1*np.sin(reflat*d2r)+cn2*np.cos(reflat*d2r)*cn6
        cn7=1/cn7
        xx=radius*cn7*cos(reflat*d2r)*cn5
        yy=radius*cn7*(cn2*np.sin(reflat*d2r)-cn1*np.cos(reflat*d2r)*cn6)
        xih= xx*cn4+yy*cn3
        yih=-xx*cn3+yy*cn4
#
        x = xih + distx
        y = yih + disty
        xpos = x
        ypos = y
        xx = cn4*x-cn3*y
        yy = cn3*x+cn4*y
        rr = np.sqrt(xx*xx + yy*yy)
        cn5 = np.arctan(rr/radius)
        cn6 = np.cos(cn5)
        cn5 = np.sin(cn5)
        if rr == 0:
           grdlon = stdlon
           grdlat = stdlt1
        else:
           grdlon = stdlon + r2d * np.arctan(xx*cn5/(rr*cn2*cn6-yy*cn1*cn5))
           grdlat = r2d * np.arcsin(cn6*cn1 + yy*cn5*cn2/rr)
        end
        deg = grdlat * d2r
        f = omega4 * np.sin (deg)
        grdlon[grdlon>360]=grdlon[grdlon>360]-360
        grdlon[grdlon<0  ]=grdlon[grdlon<0  ]+360

        hx = np.ones(m,n)
        hy = hx

      else:
        rdir = self.restartdir
        dtg = self.dtg

        grdlon = grdi * np.nan
        grdlat = grdlon
        f = grdlon
        hx = grdlon
        hy = grdlon
        xpos = grdlon
        ypos = grdlon

        if os.path.isdir(rdir):
          fn=rdir+'/'+'grdlon'+ \
            '_sfc_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(0,0,m,n) \
            +dtg+'_00000000_'+'datafld'
          f = open(fn, "rb")
          data = np.fromfile(f, dtype=">f4", count=-1)
          grdlon = np.reshape(data, (m, n), order='F')
          f.close()
          fn=rdir+'/'+'grdlat'+ \
            '_sfc_%6.6d_%6.6d_1o%4.4dx%4.4d_'%(0,0,m,n) \
            +dtg+'_00000000_'+'datafld'
          f = open(fn, "rb")
          data = np.fromfile(f, dtype=">f4", count=-1)
          grdlat = np.reshape(data, (m, n), order='F')
          f.close()
        else:
          print(' no dir %s\n',rdir)

      return (grdlon,grdlat,f,hx,hy,xpos,ypos)

    def grid_outline(self):

      import matplotlib.pyplot as plt
      from mpl_toolkits.basemap import Basemap

      igrid=self.nproj
      reflat=self.reflat
      reflon=self.reflon
      stdlt1=self.stdlt1
      stdlt2=self.stdlt2
      stdlon=self.stdlon

      x = dict()

      for nest in range(0,self.nnest):

        iref=(self.iref[nest]).astype(float)
        jref=(self.jref[nest]).astype(float)
        delx=self.delx[nest]
        dely=self.dely[nest]
        m=self.m[nest]
        n=self.n[nest]
        gr=self.rnest

        if nest > 0:
          npg=self.npgrid[nest]-1
          self.delx[nest]=self.delx[npg]/gr
          self.dely[nest]=self.dely[npg]/gr
          self.iref[nest]=np.int((self.iref[npg]-self.ii[nest])*gr+1)
          self.jref[nest]=np.int((self.jref[npg]-self.jj[nest])*gr+1)

        (grdlon,grdlat,f,hx,hy,xpos,ypos) = self.grid(nest)

        if nest == 0:
          n0=[grdlon.min(),grdlon.max(),grdlat.min(),grdlat.max()]
          n0=self.boundbox()
          dx=n0[1]-n0[0]
          dy=n0[3]-n0[2]
          n1=np.array([-dx,dx,-dy,dy])
          n1=n0+.2*n1

        xb=np.concatenate((grdlon[:,0],grdlon[m-1,1:],grdlon[m-2::-1,n-1],\
                           grdlon[0,n-2::-1]))
        yb=np.concatenate((grdlat[:,0],grdlat[m-1,1:],grdlat[m-2::-1,n-1],\
                           grdlat[0,n-2::-1]))

        xb[xb<n1[0]]=xb[xb<n1[0]]+360

        x0 = np.zeros(xb.size, dtype=[('x','f4'),('y','f4')])
        x0['x']=xb
        x0['y']=yb

        x[nest]=x0

      return x

    def boundbox(self,nest):

      from matplotlib.path import Path

      igrid=self.nproj
      reflat=self.reflat
      reflon=self.reflon
      stdlt1=self.stdlt1
      stdlt2=self.stdlt2
      stdlon=self.stdlon

      iref=(self.iref[nest-1]).astype(float)
      jref=(self.jref[nest-1]).astype(float)
      delx=self.delx[nest-1]
      dely=self.dely[nest-1]
      m=self.m[nest-1]
      n=self.n[nest-1]

      (grdlon,grdlat,f,hx,hy,xpos,ypos) = self.grid(nest)

    # lat range of bounding box for this domain
      y1=grdlat.min()
      y2=grdlat.max()
      x1=grdlon.min()
      x2=grdlon.max()

    # grid has pole
      if np.abs(grdlat).max() > np.max([np.abs(y1), np.abs(y2)]):
        x1=-180.
        x2=180.
        return np.array([-180., 180., y1, y2])

    # grab the four grid corners
      xcor=np.zeros([4,1])
      ycor=np.zeros([4,1])
      xcor[0]=grdlon[0,0]
      xcor[1]=grdlon[m-1,0]
      xcor[2]=grdlon[m-1,n-1]
      xcor[3]=grdlon[0,n-1]
      ycor[0]=grdlat[0,0]
      ycor[1]=grdlat[m-1,0]
      ycor[2]=grdlat[m-1,n-1]
      ycor[3]=grdlat[0,n-1]

    # cyclic global if two corners are the same....
      for j in range(0,4):
        for i in range(j+1,4):
          if xcor[i] == xcor[j] and ycor[i] == ycor[j]:
            x1=-180.
            x2=180.
            return np.array([x1, x2, y1, y2])

    # check if grid has 0E cut with 0-360 lons -
    # flip corners to -180:180
      xcor[xcor > 180.] = xcor[xcor > 180.] - 360.

      cor= np.squeeze([[xcor[0],ycor[0]],
                     [xcor[1],ycor[1]],
                     [xcor[2],ycor[2]],
                     [xcor[3],ycor[3]],
                     [xcor[0],ycor[0]]])

      polygon = Path(cor,[Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY])

    # check for points on a line from min to max lat along 0E in the
    # polygon of the four corners.  for kicks, do this at increasing (1, 0.1,
    # 0.01, 0.001 degree) resolution to save time, perhaps
      a=10.
      inp=-1
      for l in range(0,4):
        a=.1*a
        k=np.int(np.floor((y2-y1)/a)+1)
        for j in range(0,k):
          if polygon.contains_point([0., y1+j*a]):
            inp=1
            break
        if inp == 1:
          break

    # lon range of bounding box for this domain crossing 0E
      maybe_0 = False
      if inp == 1:
        x1=np.min(xcor[xcor < 0])+360.
        x2=np.max(xcor[xcor > 0])+360.
        if x2-x1 > 180.:
          maybe_0 = True
        else:
          print('box crosses 0E cut: %d %d' % (x1,x2))
          return np.array([x1, x2, y1, y2])

    # check if grid has 180E cut with -180:180 lons -
    # flip corners to 0:360
      xcor[xcor < 0.] = xcor[xcor < 0.] +360.
      xcor[xcor > 360.] = xcor[xcor > 360.] -360.

      cor= np.squeeze([[xcor[0],ycor[0]],
                     [xcor[1],ycor[1]],
                     [xcor[2],ycor[2]],
                     [xcor[3],ycor[3]],
                     [xcor[0],ycor[0]]])

      polygon = Path(cor,[Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY])

    # next, check for points on a line from min to max lat along 0E in the
    # polygon of the four corners.  for kicks, do this at increasing (1, 0.1,
    # 0.01, 0.001 degree) resolution to save time, perhaps
      a=10.
      inp=-1
      for l in range(0,4):
        a=.1*a
        k=np.int(np.floor((y2-y1)/a)+1)
        for j in range(0,k):
          if polygon.contains_point([0., y1+j*a]):
            inp=1
            break
        if inp == 1:
          break

    # yes, cut 180E
      if inp == 1:
    # lon range of bounding box for this domain crossing 180E
        x1=xcor.min()
        x2=xcor.max()
        if np.abs(x2-x1) > 180:
          if maybe_0:
            xcor[xcor > 180.] = xcor[xcor > 180.]-360.
            x1=np.min(xcor[xcor < 0.])
            x2=np.max(xcor[xcor > 0.])
            print('box crosses 0E cut: %d %d' % (x1,x2))
        else:
            print('box crosses 180E cut: %d %d' % (x1,x2))

      return np.array([x1, x2, y1, y2])

      def gridnl_default(self):

        import f90nml
        nml=f90nml.Namelist()

        mx_grids = 7

        nml['alnnt'] =  240.
        nml['delx']  =  81000.
        nml['dely']  =  81000.
        nml['ii']    =  [1] * mx_grids
        nml['iref']  =  [1] * mx_grids
        nml['jj']    =  [1] * mx_grids
        nml['jref']  =  [1] * mx_grids
        nml['kka']   =  30
        nml['kko']   =  1
        nml['kkom']  =  1
        nml['kkosm'] =  0
        nml['lncom'] =  False
        nml['lnmove']=  [ False ] * mx_grids
        nml['m']     =  61
        nml['n']     =  61
        nml['nnest'] =  1
        nml['npgrid']=  [1] + list(range(1,mx_grids))
        nml['nproj'] =  2
        nml['phnt1'] =  60.
        nml['phnt2'] =  30.
        nml['rlat']  =  42.5
        nml['rlon']  =  16.5
        nml['rnest'] =  3.

        return nml

    def get_gridnl(self):

        import f90nml
        gridnl=f90nml.Namelist()
        nml=f90nml.Namelist()

        gridnl['delx']   = self.delx.tolist()
        gridnl['dely']   = self.dely.tolist()
        gridnl['ii']     = self.ii.tolist()
        gridnl['iref']   = self.iref.tolist()
        gridnl['jj']     = self.jj.tolist()
        gridnl['jref']   = self.jref.tolist()
        gridnl['m']      = self.m.tolist()
        gridnl['n']      = self.n.tolist()
        gridnl['nproj']  = self.nproj
        gridnl['alnnt']  = self.stdlon
        gridnl['npgrid'] = self.npgrid.tolist()
        gridnl['kka']    = self.kka
        gridnl['kko']    = self.kko
        gridnl['kkom']   = self.kkom
        gridnl['kkosm']  = self.kkosm
        gridnl['lncom']  = self.lncom
        gridnl['lnmove'] = list(self.lnmove)
        gridnl['nnest']  = self.nnest
        gridnl['phnt1']  = self.stdlt1
        gridnl['phnt2']  = self.stdlt2
        gridnl['rlat']   = self.reflat
        gridnl['rlon']   = self.reflon
        gridnl['rnest']  = self.rnest

        nml['gridnl']=gridnl
        return nml

def test0():

  import matplotlib.pyplot as plt
  import f90nml
  import json

  agrid = COAMPSGrid(" ")
  agrid.define(1,24.878,270.704,38,36,(38,14,21),(36,15,21),22.5,-99.0,89.296,\
               54000.0,54000.0,(75,145,322),(71,124,250),3,(1,1,2,3),3.)

  ogrid = COAMPSGrid(" ")
  ogrid.define(1,24.572,270.704,291,251,291,251,22.5,-99.0,270.704,\
               3000.0,3000.0,670,500,1,1,3.)

  wgrid = COAMPSGrid(" ")
  wgrid.define(1,24.572,270.704,291,251,291,251,22.5,30.0,270.509,\
               6000.0,6000.0,350,250,1,1,3.)
 
  agridnml = agrid.get_gridnl()
  ogridnml = ogrid.get_gridnl()
  wgridnml = wgrid.get_gridnl()

  try:
    agridnml.write('atmos.gridnl.test',force=True)
    ogridnml.write('ocean.gridnl.test',force=True)
    wgridnml.write('wave.gridnl.test',force=True)
  except:
    print(json.dumps(agridnml,indent=2))
    print(json.dumps(ogridnml,indent=2))
    print(json.dumps(wgridnml,indent=2))

  print("atmos grid bounding box: ",agrid.boundbox())
  print("ocean grid bounding box: ",ogrid.boundbox())
  print("wave  grid bounding box: ",wgrid.boundbox())

if __name__ == "__main__":
  test0()

