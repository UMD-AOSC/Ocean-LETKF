#!/usr/bin/env python3

import argparse, os
import numpy as np
import xesmf as xe
from netCDF4 import Dataset
from config_tools import Config

class ConfigRegrid(Config):
    def __init__(self):
        # static config applied for each processing
        self.ocean_static_path = "ocean_static.nc"
        self.mom_lat_var = "geolat"
        self.mom_lon_var = "geolon"
        self.sst_grid_path = "OSTIA.nc"
        self.sst_lat_var = "lat"
        self.sst_lon_var = "lon"
        self.regridder_file_path = "xesmf_wts.nc"


def parseCommandLine():
    parser = argparse.ArgumentParser(description=("Generate regrid weight file by xesmf"))
    parser.add_argument("config_path", help=())
    args = parser.parse_args()
    print(args)
    return args

def main(args):

    #
    # read configurations
    #
    cfg = ConfigRegrid()
    cfg.read_cfg(args.config_path,namelist="regrid")

    # MOM latlon grids as output grids
    fn_grd_out = os.path.abspath(cfg.ocean_static_path)
    f = Dataset(fn_grd_out)
    lat2d_grd_out = f.variables[cfg.mom_lat_var][:].astype(np.float32)
    lon2d_grd_out = f.variables[cfg.mom_lon_var][:].astype(np.float32)
    f.close()

    print("lat2d_grd_out: shape, min, max=",lat2d_grd_out.shape, np.min(lat2d_grd_out), np.max(lat2d_grd_out))
    print("lon2d_grd_out: shape, min, max=",lon2d_grd_out.shape, np.min(lon2d_grd_out), np.max(lon2d_grd_out))

    # L4 latlon grids as input grids
    fn_grd_in = os.path.abspath(cfg.sst_grid_path)
    f = Dataset(fn_grd_in)
    lat1d_grd_in = f.variables[cfg.sst_lat_var][:]
    lon1d_grd_in = f.variables[cfg.sst_lon_var][:]
    f.close()

    print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
    print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))

    # generate the regridder 
    grd_in = {"lon": lon1d_grd_in, "lat": lat1d_grd_in}  
    grd_out = {"lon": lon2d_grd_out, "lat": lat2d_grd_out}
    interp_method = "bilinear"
    periodic = True
    regridder = xe.Regridder(grd_in, grd_out, interp_method, periodic=periodic)
  
    fn_regridder_out = os.path.abspath(cfg.regridder_file_path)
    print("write regrid info into the file: ", fn_regridder_out)
    regridder.filename = fn_regridder_out
    regridder.to_netcdf()

    print(regridder)

    return None

if __name__ == '__main__':
    args = parseCommandLine()
    main(args)

