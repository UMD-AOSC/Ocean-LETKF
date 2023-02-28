#
# the grid file control structure (CTL)
#
import numpy as np
import numpy.matlib
from scipy import interpolate

class grdctl:
    def __init__(self,nx,ny,nz,nv3d,nv2d,v3dnames,v2dnames,v3dncnames=None,v2dncnames=None):
        self.nx       = nx
        self.ny       = ny
        self.nz       = nz

        self.nv3d     = nv3d
        self.nv2d     = nv2d
        self.v3dnames = v3dnames.copy()
        self.v2dnames = v2dnames.copy()

        self.v3dncnames = []
        if v3dncnames is not None:
            self.v3dncnames = v3dncnames.copy()
        self.v2dncnames = []
        if v2dncnames is not None:
            self.v2dncnames = v2dncnames.copy()

        self.lats     = np.zeros(ny)
        self.lons     = np.zeros(nx)
        self.lats2d   = np.zeros([ny,nx])
        self.lons2d   = np.zeros([ny,nx])
        self.wts2d    = np.zeros([ny,nx])
        self.levs     = np.zeros(nz)

    def __str__(self):
        print("nx=",self.nx)
        print("ny=",self.ny)
        print("nz=",self.nz)
        print("nv3d=",self.nv3d)
        print("nv2d=",self.nv2d)
        print("v3dnames=",self.v3dnames)
        print("v2dnames=",self.v2dnames)
        print("v3dncnames=",self.v3dncnames)
        print("v2dncnames=",self.v2dncnames)
        print("lats=",self.lats)
        print("lons=",self.lons)
        print("levs=",self.levs)
        return 'type::grdctl'

    def genGauLatsLons(self,N2S=False):
        self.lons[:]=np.linspace(0,360,self.nx,endpoint=False)
        slats, wts =np.polynomial.legendre.leggauss(self.ny)
        self.lats[:] = 180.0/np.pi*np.arcsin(slats)

        latsm1 = np.insert(0.5*(self.lats[1:]+self.lats[:-1]),0,-90.0)
        latsm2 = np.append(0.5*(self.lats[1:]+self.lats[:-1]),90.0)
        wts2 = np.abs(np.sin(latsm2*np.pi/180.0)-np.sin(latsm1*np.pi/180.0))
        buf = np.matlib.repmat(wts2,self.nx,1)
        self.wts2d = buf.transpose().copy()

        if N2S:
            self.lats[:] = self.lats[::-1].copy()
            self.wts2d[:,:] = self.wts2d[:,::-1].copy()
        self.lons2d,self.lats2d = np.meshgrid(self.lons,self.lats)

    def setLevs(self, levs):
        self.levs = np.array(levs).copy()


    def remapTo3d(self, v3dFrom, ctlTo):
        v3dTo={}
        for var in ctlTo.v3dnames:
            buf = np.zeros([ctlTo.nz,ctlTo.ny,ctlTo.nx])
            v3dTo.update({var: buf.copy()})
            for lev in ctlTo.levs:
                #print("var, lev=",var,lev)
                kn = ctlTo.levs[:].tolist().index(lev)
                #print("kn=",kn,ctlTo.levs[kn])
                ko = self.levs[:].tolist().index(lev)
                #print("ko=",ko,self.levs[ko])
                o2n = interpolate.interp2d(self.lons,self.lats,v3dFrom[var][ko,:,:],kind='cubic')
                v3dTo[var][kn,:,:] = o2n(ctlTo.lons,ctlTo.lats).copy()

        return v3dTo


    def sprd(self,field2d,simpleAvg=False,verbose=True):
        """
        calculate the average spread of analaysis for a 2-D field
        global  -90 ->  90
        tropics -20 ->  20
        nh       20 ->  90
        sh      -90 -> -20
        """

        if simpleAvg:
            sprd_gl = np.sqrt( np.mean(field2d**2) )

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array(field2d**2, mask=m)
            sprd_nh = np.sqrt( md.mean() )

            m = self.lats2d > -20 # removed points
            md = np.ma.array(field2d**2, mask=m)
            sprd_sh = np.sqrt( md.mean() )

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array(field2d**2, mask=m)
            sprd_tr = np.sqrt( md.mean() )
        else:
            md = self.wts2d*field2d**2
            mw = self.wts2d
            sprd_gl = np.sqrt( md.sum()/mw.sum() )

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array(self.wts2d*field2d**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            sprd_nh = np.sqrt( md.sum()/mw.sum() )

            m = self.lats2d > -20 # removed points
            md = np.ma.array(self.wts2d*field2d**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            sprd_sh = np.sqrt( md.sum()/mw.sum() )

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array(self.wts2d*field2d**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            sprd_tr = np.sqrt( md.sum()/mw.sum() )
           
        if verbose: 
            print("sprd: ",sprd_gl, sprd_tr, sprd_nh, sprd_sh)
        return sprd_gl, sprd_tr, sprd_nh, sprd_sh


    def aave(self,field2d,simpleAvg=False,verbose=True):
        """
        calculate the area average of a 2-D field
        global  -90 ->  90
        tropics -20 ->  20
        nh       20 ->  90
        sh      -90 -> -20
        """

        if simpleAvg:
            aave_gl = field2d.mean()

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array(field2d, mask=m)
            aave_nh = md.mean()

            m = self.lats2d > -20 # removed points
            md = np.ma.array(field2d, mask=m)
            aave_sh = md.mean()

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array(field2d, mask=m)
            aave_tr = md.mean()
        else:
            md = self.wts2d*field2d
            mw = self.wts2d
            aave_gl = md.sum()/mw.sum() 

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array(self.wts2d*field2d, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            aave_nh = md.sum()/mw.sum()

            m = self.lats2d > -20 # removed points
            md = np.ma.array(self.wts2d*field2d, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            aave_sh = md.sum()/mw.sum()

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array(self.wts2d*field2d, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            aave_tr = md.sum()/mw.sum()
                
        if verbose: 
            print("aave: ",aave_gl, aave_tr, aave_nh, aave_sh)
        return aave_gl, aave_tr, aave_nh, aave_sh


    def rmse(self,ref2d,ana2d,simpleAvg=False,verbose=True):
        """
        calculate the rmse difference of analaysis from the reference for a 2-D field
        global  -90 ->  90
        tropics -20 ->  20
        nh       20 ->  90
        sh      -90 -> -20
        """

        if simpleAvg:
            rmse_gl = np.sqrt( np.mean((ref2d-ana2d)**2) )

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array((ref2d-ana2d)**2, mask=m)
            rmse_nh = np.sqrt(md.mean())

            m = self.lats2d > -20 # removed points
            md = np.ma.array((ref2d-ana2d)**2, mask=m)
            rmse_sh = np.sqrt(md.mean())

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array((ref2d-ana2d)**2, mask=m)
            rmse_tr = np.sqrt(md.mean())
        else:
            md = self.wts2d*(ref2d-ana2d)**2
            mw = self.wts2d
            rmse_gl = np.sqrt( md.sum()/mw.sum() )

            m = self.lats2d < 20 # removed pointA
            md = np.ma.array(self.wts2d*(ref2d-ana2d)**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            rmse_nh = np.sqrt(md.sum()/mw.sum())

            m = self.lats2d > -20 # removed points
            md = np.ma.array(self.wts2d*(ref2d-ana2d)**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            rmse_sh = np.sqrt(md.sum()/mw.sum())

            m = np.logical_or(self.lats2d >= 20, self.lats2d <= -20) # removed points
            md = np.ma.array(self.wts2d*(ref2d-ana2d)**2, mask=m)
            mw = np.ma.array(self.wts2d, mask=m)
            rmse_tr = np.sqrt(md.sum()/mw.sum())
                
        if verbose: 
            print("rmse: ",rmse_gl, rmse_tr, rmse_nh, rmse_sh)
        return rmse_gl, rmse_tr, rmse_nh, rmse_sh


        
if __name__ == "__main__":
    
    gfsctl = grdctl(nx=192,ny=94,nz=64,nv3d=8,nv2d=2,v3dnames=[],v2dnames=[])
    print(gfsctl)
    gfsctl.genGauLatLon()
    print(gfsctl.lats)
    print(gfsctl.lons)

