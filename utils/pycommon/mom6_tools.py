import yaml
import os, sys
import datetime as dt
from netCDF4 import Dataset
import numpy as np

def parse_time_units(tstr):

    if tstr.lower().find("seconds") >= 0:
       dtime_units = dt.timedelta(seconds=1)
    elif tstr.lower().find("days") >= 0:
       dtime_units = 24*3600*dt.timedelta(seconds=1)
    elif tstr.lower().find("minutes") >= 0:
       dtime_units = 60*dt.timedelta(seconds=1)
    elif tstr.lower().find("hours") >= 0:
       dtime_units = 3600*dt.timedelta(seconds=1)
    else:
       raise RuntimeError("cannot determine the units in {}".format(tstr))

    ib = tstr.lower().find("since")
    if ib > 0:
        ib += 6 # "since 1920-01-02 12:00:00"
        try:
            slen = 19
            bdate = dt.datetime.strptime(tstr[ib:ib+slen],"%Y-%m-%d %H:%M:%S")
        except:
            try:
                slen = 17
                bdate = dt.datetime.strptime(tstr[ib:ib+slen],"%Y-%m-%d %H:%M")
            except:
                try:
                    slen = 15
                    bdate = dt.datetime.strptime(tstr[ib:ib+slen],"%Y-%m-%d %H")
                except:
                    try:
                        slen = 13
                        bdate = dt.datetime.strptime(tstr[ib:ib+slen],"%Y-%m-%d")
                    except:
                        raise RuntimeError("cannot identify base time {}".format(tstr))
    return dtime_units, bdate


class RestoreTemplate:
    def __init__(self):
        self.tpl = None

    def read(self,fnin):
        with open(fnin, "r") as f:
             self.tpl = yaml.safe_load(f)

    def create_restore(self, fnout, input_dims={"time":1,"i":10,"q":5,"IQ":11,"JQ":6}, input_vars=None):
        if os.path.exists(os.path.abspath(fnout)):
            raise RuntimeError("fnout ({}) already exists.".format(os.path.abspath(fnout)) )
            sys.exit(1)

        fnout = os.path.abspath(fnout)
        f = Dataset(fnout,mode='w',format='NETCDF4_CLASSIC')

        tpl_dims = self.tpl['dimensions']
        for dname in tpl_dims.keys():
            if tpl_dims[dname]["unlimited"]:
               dlen = 0
            else:
               if tpl_dims[dname]['len_from_file']:
                  dlen = tpl_dims[dname]['len_from_file'] 
               else:
                  dlen = input_dims[dname]
            f.createDimension(dname, dlen)

        tpl_vars = self.tpl['variables']
        dict_ncvars = {}
        for vname in tpl_vars.keys():
            vdim  = tpl_vars[vname]['dims']
            vtype = tpl_vars[vname]['type']
            vfill = tpl_vars[vname]['_FillValue']
            dict_ncvars[vname] = f.createVariable(vname, vtype, vdim, fill_value=vfill)

            vatts = tpl_vars[vname]['atts']
            for att in vatts.keys():
                setattr(dict_ncvars[vname], att, vatts[att])

        if input_vars is not None:
           for v in input_vars.keys():
               if v not in dict_ncvars.keys():
                   raise RuntimeError("variable {} in input_vars not found in the yaml template".format(v))
                   sys.exit(2)

               dict_ncvars[v][:] = input_vars[v].astype(dict_ncvars[v].dtype)

        f.close()


    def get_restore_static(self, ocean_hgrid_file, ocean_static_file):
        f = Dataset(ocean_hgrid_file)
        x2d = f.variables["x"][:,:]
        y2d = f.variables["y"][:,:]
        f.close()
        lon_crnr = x2d[::2,::2]
        lat_crnr = y2d[::2,::2]

        f = Dataset(ocean_static_file)
        lon = f.variables["geolon"][:,:]
        lat = f.variables["geolat"][:,:]
        wet = f.variables["wet"][:,:]
        area = f.variables["area_t"][:,:]
        f.close()

        yh, xh = lon.shape

        i  = np.arange(0.5, xh+0.5, 1.0)
        j  = np.arange(0.5, yh+0.5, 1.0)
        IQ = np.arange(0.0, xh+1.0, 1.0)
        JQ = np.arange(0.0, yh+1.0, 1.0)

        input_dims = {"time":0, "i": len(i), "j": len(j), \
                      "IQ": len(IQ), "JQ": len(JQ) }

        input_vars = {"i":i, "j":j, "IQ":IQ, "JQ":JQ, \
                      "lon": lon, "lat": lat, \
                      "lon_crnr": lon_crnr, "lat_crnr": lat_crnr, \
                      "wet": wet, "area": area}

        return input_dims, input_vars


if __name__ == '__main__':
    tpl = RestoreTemplate()
    tpl.read("static_restore_file.yaml")
    print(tpl.tpl)
    input_dims, input_vars = tpl.get_restore_static(\
                              ocean_static_file=os.path.abspath("./ocean_static.nc"), \
                              ocean_hgrid_file=os.path.abspath("./ocean_hgrid.nc"))
    tpl.create_restore("test_static.nc",input_dims, input_vars)
        
