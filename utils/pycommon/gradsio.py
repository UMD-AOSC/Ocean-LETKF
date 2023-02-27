import numpy as np
from grdctl import grdctl


class gradsctl (grdctl):

    def readVars3d2d(self, filename, vtype='>f4'):
        #print("nx,ny,nz=",(nx,ny,nz))
        #print("nv3d,nv2d=",(nv3d,nv2d))
        nrec=self.nx*self.ny*(self.nv3d*self.nz+self.nv2d)
        nlevall=self.nv3d*self.nz+self.nv2d
        f=open(filename,'rb')
        buf=np.fromfile(f,dtype=vtype,count=nrec)
        f.close()

        buf2=np.reshape(buf,(nlevall,self.ny,self.nx))
        #print(type(buf2))
        #print(buf2.shape)
    
        v3d={}
        for ivar in range(0,self.nv3d):
            js=self.nz*ivar
            je=js+self.nz
            #print("ivar3d: js,je=",v3dnames[ivar],js,je)
            v3d.update({self.v3dnames[ivar]: buf2[js:je,:,:].copy() })

        v2d={}
        for ivar in range(0,self.nv2d):
            js=self.nv3d*self.nz+ivar
            je=js+1
            #print("ivar2d: js,je=",v2dnames[ivar],js,je)
            v2d.update({self.v2dnames[ivar]: np.squeeze(buf2[js:je,:,:]).copy() })
            #print("var: min,max=",v2dnames[ivar],v2d[v2dnames[ivar]].min(), v2d[v2dnames[ivar]].max())

        return v3d, v2d

    def readVars3d2d_nm(self,filename,vtype='>f4'):
        """ read 3d,2dvar from Grads stream files with no name attachement to the output
        """
        nrec=self.nx * self.ny *(self.nv3d*self.nz+self.nv2d)
        nlevall=self.nv3d*self.nz+self.nv2d
        f=open(filename,'rb')
        buf=np.fromfile(f,dtype=vtype,count=nrec)
        f.close()

        buf2=np.reshape(buf,(nlevall,self.ny,self.nx))
        #print(type(buf2))
        #print(buf2.shape)
        
        v3d=np.zeros((self.nv3d,self.nz,self.ny,self.nx))
        v2d=np.zeros((self.nv2d,self.ny,self.nx))
        for ivar in range(0,v3d.shape[0]):
            js=self.nz*ivar
            je=js+self.nz
            #print("ivar3d: js,je=",ivar,js,je)
            v3d[ivar,0:self.nz,:,:] = buf2[js:je,:,:].copy()

        nz3d=self.nv3d*self.nz
        for ivar in range(0,v2d.shape[0]):
            js=nz3d+ivar
            je=js+1
            #print("ivar2d: js,je=",ivar,js,je)
            v2d[ivar,:,:] = buf2[js:je,:,:].copy()

        return v3d, v2d


def read_grads2(filename,nv3d,nv2d,v3dnames,v2dnames,nx,ny,nz,vtype='>f4'):
    #print("nx,ny,nz=",(nx,ny,nz))
    #print("nv3d,nv2d=",(nv3d,nv2d))
    nrec=nx*ny*(nv3d*nz+nv2d)
    nlevall=nv3d*nz+nv2d
    f=open(filename,'rb')
    buf=np.fromfile(f,dtype=vtype,count=nrec)
    f.close()

    buf2=np.reshape(buf,(nlevall,ny,nx))
    #print(type(buf2))
    #print(buf2.shape)
    
    v3d={}
    for ivar in range(0,nv3d):
        js=nz*ivar
        je=js+nz
        #print("ivar3d: js,je=",v3dnames[ivar],js,je)
        v3d.update({v3dnames[ivar]: buf2[js:je,:,:].copy() })

    v2d={}
    for ivar in range(0,nv2d):
        js=nv3d*nz+ivar
        je=js+1
        #print("ivar2d: js,je=",v2dnames[ivar],js,je)
        v2d.update({v2dnames[ivar]: np.squeeze(buf2[js:je,:,:]).copy() })
        #print("var: min,max=",v2dnames[ivar],v2d[v2dnames[ivar]].min(), v2d[v2dnames[ivar]].max())

    return v3d, v2d


def read_grads(filename,nv3d,nv2d,nx,ny,nz,vtype='>f4'):
    #print("nx,ny,nz=",(nx,ny,nz))
    #print("nv3d,nv2d=",(nv3d,nv2d))
    nrec=nx*ny*(nv3d*nz+nv2d)
    nlevall=nv3d*nz+nv2d
    f=open(filename,'rb')
    buf=np.fromfile(f,dtype=vtype,count=nrec)
    f.close()

    buf2=np.reshape(buf,(nlevall,ny,nx))
    #print(type(buf2))
    #print(buf2.shape)
    
    v3d=np.zeros((nv3d,nz,ny,nx))
    v2d=np.zeros((nv2d,ny,nx))
    for ivar in range(0,v3d.shape[0]):
        js=nz*ivar
        je=js+nz
        #print("ivar3d: js,je=",ivar,js,je)
        v3d[ivar,0:nz,:,:] = buf2[js:je,:,:].copy()

    nz3d=nv3d*nz
    for ivar in range(0,v2d.shape[0]):
        js=nz3d+ivar
        je=js+1
        #print("ivar2d: js,je=",ivar,js,je)
        v2d[ivar,:,:] = buf2[js:je,:,:].copy()

    return v3d, v2d



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    import os

    #define CTL
    gfs = gradsctl( nx=192, ny=94, nz=12, nv3d=8, nv2d=10, \
                    v3dnames=['u','v','t','q','rh','qc','gph','o3'], \
                    v2dnames=['ps','slp','orog','slmsk','tsea','u10m','v10m','t2m','q2m','tprcp'])
    gfs.genGauLatsLons()
    gfs.setLevs(levs=[1000,925,850,700,500,300,200,100,50,20,10,5])

    # load data
    v3d,v2d = gfs.readVars3d2d(fnin)

    for var in v3d:
        print("name: min, max=",var, v3d[var].min(), v3d[var].max())

    for var in v2d:
        print("name: min, max=",var, v2d[var].min(), v2d[var].max())

    # plot
    pngdir='checkall'
    os.mkdir(pngdir)

    for var in v2d:
        print("plot v2d:",var)
        fig = plt.figure()
        m = Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=0,urcrnrlat=90,urcrnrlon=360)
        m.drawcoastlines(linewidth=1)
        h = m.contourf(gfs.lon2d,gfs.lat2d,v2d[var][:,:],cmap=plt.cm.jet)
        cb = m.colorbar(h)
        #h.set_clim(vmin=,vmax=)
        plt.title(var)
        #plt.show()
        fnpng = pngdir+'/'+var+'.png'
        fig.savefig(fnpng)
        plt.close(fig)


    for var in v3d:
        print("plot v3d:",var)
        os.mkdir(pngdir+'/'+var)
        for lv in range(0,len(gfs.levs)):   
            fig = plt.figure()
            m = Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=0,urcrnrlat=90,urcrnrlon=360)
            m.drawcoastlines(linewidth=1)
            h = m.contourf(gfs.lon2d,gfs.lat2d,v3d[var][lv,:,:],cmap=plt.cm.jet)
            cb = m.colorbar(h)
            plt.title(var+' at '+str(levs[lv]) )
            fnpng=pngdir+'/'+var+'/'+var+'_'+str(levs[lv])+'.png'
            fig.savefig(fnpng)
            plt.close(fig)


